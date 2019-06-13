# Markus Reinert, May 2019
#
# The fluid2d Experiment Management System (EMS).
#
# This module provides a way to handle big sets of fluid2d-experiments.
# It creates a database in which information about all the experiments
# is stored.  This allows to easily show, compare, search for, filter
# and save comments to already executed experiments.  Furthermore, it
# provides, after set up for a specific experiment, a simpler interface
# to modify the parameteres of interest.  To avoid overwriting files,
# it automatically assigns a unique identifier to every new experiment.
#
# Its usage is explained in the breaking waves experiment.

import os
import datetime
import sqlite3 as dbsys


# Columns with the following names are automatically added to every table of the database.
# Therefore they must not be used as parameter names.
RESERVED_COLUMN_NAMES = [
    'id',
    'datetime',
    'duration',
    'size_total',
    'size_mp4',
    'size_his',
    'size_diag',
    'size_flux',
    'comment',
]


class InputMismatch(Exception):
    """Incompatibility of given experiment file and database."""
    pass


class ParamError(Exception):
    """Non-conformance of user-set attributes of param."""
    pass


class ExpFileError(Exception):
    """Malformed experiment file."""
    pass


class EMS:
    """Experiment Management System"""

    def __init__(self, param: "param.Param", experiment_file: str=""):
        """Constructor of the Experiment Management System.

        Set up a connection to the database.  If an experiment file is
        given, a new entry in the database is created for this file.
        Otherwise an existing entry corresponding to `param.expname` is
        used.  The parameters of this entry are loaded.
        They can be accessed via the method `get_parameters`."""

        # Get the path to the database and create a connection
        datadir = param.datadir
        if datadir.startswith("~"):
            datadir = os.path.expanduser(param.datadir)
        dbpath = os.path.join(datadir, "experiments.db")
        print("-"*50)
        print(" Opening database {}.".format(dbpath))
        self.connection = dbsys.connect(dbpath)
        cursor = self.connection.cursor()

        if experiment_file:
            # If an experiment file is given, parse it and add a new entry to the database.
            self.exp_class, self.description, param_list = parse_experiment_file(experiment_file)
            if not self.exp_class:
                raise ExpFileError("no name given in file {}.".format(experiment_file))
            # Add static parameter fields
            param_list = (
                [
                    ["INTEGER", "id", -1],
                    ["TEXT", "datetime",
                     datetime.datetime.now().isoformat(timespec="microseconds")],
                ]
                + param_list
                + [
                    ["REAL", "duration", -1],
                    ["REAL", "size_total", -1],
                    ["REAL", "size_mp4", -1],
                    ["REAL", "size_his", -1],
                    ["REAL", "size_diag", -1],
                    ["REAL", "size_flux", -1],
                    ["TEXT", "comment", self.description],
                ]
                # Everything added here should be in RESERVED_COLUMN_NAMES
            )
            # Check whether table exists already
            if table_exists(cursor, self.exp_class):
                # Check if table has the same columns
                cursor.execute('PRAGMA table_info("{}")'.format(self.exp_class))
                for column in cursor.fetchall():
                    col_index = column[0]
                    col_name = column[1]
                    col_type = column[2]
                    if (col_type != param_list[col_index][0] or
                        col_name != param_list[col_index][1]):
                        raise InputMismatch(
                            "column {} of the database ({} {}) does not match "
                            "the corresponding parameter ({} {}) in the file {}."
                            .format(col_index + 1, col_type, col_name,
                                    *param_list[col_index][:2], experiment_file)
                        )
                # Get the highest index
                cursor.execute('SELECT id from "{}" ORDER BY id DESC'.format(self.exp_class))
                highest_entry = cursor.fetchone()
                if highest_entry:
                    self.id_ = highest_entry[0] + 1
                else:
                    self.id_ = 1
            else:
                # Create a new table
                print(' Creating new table "{}".'.format(self.exp_class))
                column_string = ", ".join(['"{}" {}'.format(n, t) for t, n, v in param_list])
                sql_command = 'CREATE TABLE "{}" ({})'.format(self.exp_class, column_string)
                cursor.execute(sql_command)
                # First entry has index 1 (one)
                self.id_ = 1
            # Set id in the parameter list
            param_list[0][2] = self.id_
            # Add a new entry to the table
            print(' Adding new entry #{} to table "{}".'.format(self.id_, self.exp_class))
            value_list = [v for t, n, v in param_list]
            sql_command = (
                'INSERT INTO "{}" VALUES ('.format(self.exp_class)
                + ', '.join(['?'] * len(value_list))
                + ')'
            )
            cursor.execute(sql_command, value_list)
            # Save the database
            self.connection.commit()
            # Set the name of the experiment
            param.expname = "{}_{:03}".format(self.exp_class, self.id_)
        else:
            # Get name and id of the experiment
            expname_parts = param.expname.split('_')
            if len(expname_parts) == 1:
                raise ParamError(
                    'param.expname is not a valid database entry: "{}"'.format(param.expname)
                )
            self.exp_class = '_'.join(expname_parts[:-1])
            try:
                self.id_ = int(expname_parts[-1])
            except ValueError:
                raise ParamError(
                    'param.expname is not a valid database entry: "{}"'.format(param.expname)
                )
            print(' Reading from entry #{} of table "{}".'.format(self.id_, self.exp_class))

        # Remember directory of the experiment
        self.output_dir = os.path.join(datadir, param.expname)

        # Get columns of the table and their value for the current experiment
        self.params = dict()
        cursor.execute('PRAGMA table_info("{}")'.format(self.exp_class))
        columns = cursor.fetchall()
        cursor.execute('SELECT * FROM "{}" WHERE id = ?'.format(self.exp_class), (self.id_,))
        values = cursor.fetchone()
        for column in columns:
            col_index = column[0]
            col_name = column[1]
            value = values[col_index]
            if col_name in RESERVED_COLUMN_NAMES:
                continue
            else:
                if value == "True":
                    value = True
                elif value == "False":
                    value = False
                self.params[col_name] = value

    def __del__(self):
        """Destructor of the Experiment Management System.

        Save the database and close the connection to it."""
        # Save and close database
        print(" Closing database.")
        print("-"*50)
        self.connection.commit()
        self.connection.close()

    def get_parameters(self):
        """Get user-set experiment parameters and values as dictionary."""
        return self.params

    def finalize(self, fluid2d):
        """Save integration time and size of output files in database.

        This method must be called when the simulation is finished,
        that means, after the line with `f2d.loop()`.
        It also sets the field `datetime` to the current time."""

        # Write duration of the run into the database
        self.connection.execute(
            'UPDATE "{}" SET duration = ? WHERE id = ?'.format(self.exp_class),
            (round(fluid2d.t, 2), self.id_,)
        )
        # Update date and time in the database
        self.connection.execute(
            'UPDATE "{}" SET datetime = ? WHERE id = ?'.format(self.exp_class),
            (datetime.datetime.now().isoformat(timespec="microseconds"), self.id_)
        )
        # Write size of output into the database
        # Divide size by 1000*1000 = 1e6 to get value in MB
        try:
            # Get list of files in the output directory
            output_files = [os.path.join(self.output_dir, f) for f in os.listdir(self.output_dir)]
            total_size = sum(os.path.getsize(path) for path in output_files)
        except FileNotFoundError as e:
            print(" Error getting total size of output:", e)
        else:
            self.connection.execute(
                'UPDATE "{}" SET size_total = ? WHERE id = ?'.format(self.exp_class),
                (round(total_size / 1e6, 3), self.id_,)
            )
        # History file
        try:
            his_size = os.path.getsize(fluid2d.output.hisfile)
        except FileNotFoundError as e:
            print(" Error getting size of his-file:", e)
        else:
            self.connection.execute(
                'UPDATE "{}" SET size_his = ? WHERE id = ?'.format(self.exp_class),
                (round(his_size / 1e6, 3), self.id_,)
            )
        # Diagnostics file
        try:
            diag_size = os.path.getsize(fluid2d.output.diagfile)
        except FileNotFoundError:
            print(" Error getting size of diag-file:", e)
        else:
            self.connection.execute(
                'UPDATE "{}" SET size_diag = ? WHERE id = ?'.format(self.exp_class),
                (round(diag_size / 1e6, 3), self.id_,)
            )
        # MP4 file
        mp4_size = -1.0
        if fluid2d.plot_interactive and hasattr(fluid2d.plotting, 'mp4file'):
            try:
                mp4_size = os.path.getsize(fluid2d.plotting.mp4file)
            except FileNotFoundError:
                print(" Error getting size of mp4-file:", e)
        else:
            mp4_size = 0.0
        if mp4_size >= 0:
            self.connection.execute(
                'UPDATE "{}" SET size_mp4 = ? WHERE id = ?'.format(self.exp_class),
                (round(mp4_size / 1e6, 3), self.id_,)
            )
        # Flux file
        flux_size = -1.0
        if fluid2d.diag_fluxes:
            try:
                flux_size = os.path.getsize(fluid2d.output.flxfile)
            except FileNotFoundError:
                print(" Error getting size of flux-file:", e)
        else:
            flux_size = 0.0
        if flux_size >= 0:
            self.connection.execute(
                'UPDATE "{}" SET size_flux = ? WHERE id = ?'.format(self.exp_class),
                (round(flux_size / 1e6, 3), self.id_,)
            )
        # Save the database
        self.connection.commit()


def parse_experiment_file(path: str):
    """Parse the given experiment file.

    An experiment file
     - must provide a name,
     - can provide a description,
     - can provide parameters with values.

    The name is written in a line starting with "Name:".  It must be a
    valid string to be used as a filename;  in particular, it must not
    contain the slash (/) or the quotation mark (").

    The description begins in or after a line starting with
    "Description:" and goes until the beginning of the parameters
    or the end of the file.

    The parameters follow after a line starting with "Parameters:".
    Every parameter is written in its own line.  This line contains the
    datatype, the name and the value of the parameter seperated by one
    or several whitespaces.
    The datatype must be one of "int", "float", "bool" or "str".
    The name must not contain whitespace characters or quotation marks
    and must not be in the list of reserved column names.
    The value must be a valid value for the given datatype.  If the
    value is omitted, it defaults to zero for numbers, True for booleans
    and the empty string for strings.  The values "True" and "False" can
    only be used for parameters of boolean type, not as strings.
    It is planned to include the possibility of providing several values
    in order to run multiple experiments from one experiment file.

    Lines in the experiment file starting with the #-symbol are ignored.
    It is not possible to write in-line comments."""
    # TODO: allow parameter names and string-values with whitespace like "long name".
    # TODO: allow to give several values to run multiple experiments.
    with open(path) as f:
        experiment_lines = f.readlines()
    name = ""
    description = ""
    param_list = []
    reading_params = False
    reading_description = False
    for line in experiment_lines:
        if line.startswith("#"):
            # Skip comments
            continue
        elif not line.strip():
            # Skip empty lines except in the description
            if reading_description:
                description += line
            else:
                continue
        elif line.lower().startswith("name:"):
            if not name:
                name = line[5:]
                if '"' in name:
                    raise ExpFileError('name must not contain the symbol ".')
            else:
                raise ExpFileError(
                    "name defined more than once in file {}.".format(path)
                )
        elif line.lower().startswith("description:"):
            reading_description = True
            reading_params = False
            description += line[12:]
        elif line.lower().startswith("parameters:"):
            reading_description = False
            reading_params = True
        elif reading_description:
            description += line
        elif reading_params:
            words = line.split()
            if len(words) < 2:
                raise ExpFileError(
                    'type or name missing for parameter "{}" in file {}.'.format(words[0], path)
                )
            param_type = words[0].lower()
            param_name = words[1]
            param_values = words[2:]
            # Check values
            if len(param_values) == 0:
                param_value = None
            elif len(param_values) == 1:
                param_value = param_values[0]
            else:
                raise NotImplementedError(
                    'multiple values given for parameter "{}" in file {}.'
                    .format(param_name, path)
                )
            # Check type
            if param_type == "int":
                sql_type = "INTEGER"
                if param_value is None:
                    param_value = 0
                else:
                    param_value = int(param_value)
            elif param_type == "float":
                sql_type = "REAL"
                if param_value is None:
                    param_value = 0.0
                else:
                    param_value = float(param_value)
            elif param_type == "bool":
                if param_value is None:
                    param_value = "True"
                elif param_value.lower() == "true":
                    param_value = "True"
                elif param_value.lower() == "false":
                    param_value = "False"
                else:
                    raise ExpFileError(
                        'boolean parameter "{}" is neither "True" nor "False" '
                        'but "{}" in file {}.'.format(param_name, param_value, path)
                    )
                sql_type = "TEXT"
            elif param_type == "str":
                sql_type = "TEXT"
                if param_value is None:
                    param_value = ""
                elif param_value == "True" or param_value == "False":
                    raise ExpFileError(
                        'the words "True" and "False" cannot be used as value for '
                        'parameter "{}" of type string.  Instead, use boolean type '
                        'parameter in file {}.'.format(param_name, path)
                    )
            else:
                raise ExpFileError(
                    'unknown parameter type "{}" in file {}.'.format(param_type, path)
                )
            # Check name
            if param_name in RESERVED_COLUMN_NAMES:
                raise ExpFileError(
                    'reserved name used for parameter "{}" in file {}.'
                    .format(param_name, path)
                )
            if '"' in param_name:
                raise ExpFileError(
                    'name of parameter "{}" must not contain the symbol ".'
                    .format(param_name)
                )
            param_list.append([sql_type, param_name, param_value])
        else:
            raise ExpFileError(
                "unexpected line in file {}:\n{}".format(path, repr(line))
            )
    return name.strip(), description.strip(), param_list


def table_exists(cursor: dbsys.Cursor, name: str):
    """Check if table with given name exists in the connected database."""
    # Get all tables
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    # This gives a list like that: [('table1',), ('table2',)]
    for table in cursor.fetchall():
        if table[0] == name:
            return True
    return False
