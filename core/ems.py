"""EMS: The fluid2d Experiment Management System

The Experiment Management System provides a convenient way to handle big
sets of fluid2d-experiments.  It creates a database in which information
about all the performed experiments is stored.  This allows to easily
show, compare, search, filter, organise, analyse them, and a lot more.
Furthermore, the EMS offers a simpler way to modify the parameters of
interest in an experiment, including the automation of running an
experiment several times with different parameter values.  To prevent
overwriting files, the EMS automatically assigns a unique identifier for
every new run of an experiment.  However, if fluid2d runs on multiple
cores, an ID has to be specified initially.

Its usage is explained in detail in the experiment on breaking waves.

With the EMShell, a command line interface is provided to access the
experiment-database created by the EMS.

Author: Markus Reinert, May/June 2019
"""

import os
import datetime
import sqlite3 as dbsys
from collections import namedtuple
from itertools import product


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
    """Incompatibility between experiment file and corresponding table in the database."""
    pass


class ExpFileError(Exception):
    """Malformed experiment file."""
    pass


class DatabaseError(Exception):
    """Malformed table in the experiment database."""
    pass


class NotInitializedError(Exception):
    """Call to `EMS.initialize` required before this action."""
    pass


class EMS:
    """Experiment Management System"""

    def __init__(self, experiment_file: str):
        """Load the parameters from the experiment file."""
        self.connection = None
        self.cursor = None
        self.exp_class, self.exp_id, self.description, self.param_name_list, param_values_list = parse_experiment_file(experiment_file)
        self.param_combinations = product(*param_values_list)
        self.setup_next_parameters(increase_id=False)

    def __del__(self):
        """Close the connection to the experiment database if any."""
        if self.connection:
            print(" Closing database.")
            print("-"*50)
            self.connection.close()

    def get_expname(self):
        if self.exp_id is None:
            raise NotInitializedError("ID not set.")
        return "{}_{:03}".format(self.exp_class, self.exp_id)

    def setup_next_parameters(self, increase_id=True):
        try:
            self.parameters = {
                name: val for name, val in
                zip(self.param_name_list, next(self.param_combinations))
            }
        except StopIteration:
            self.parameters = None
        if increase_id:
            self.exp_id += 1

    def initialize(self, data_directory: str):
        if not self.connection:
            self.connect(data_directory)
        # All of the names in the next two lists should be in RESERVED_COLUMN_NAMES
        DBColumn = namedtuple("DBColumn", ["sql_type", "name", "value"])
        static_columns_start = [
            DBColumn("INTEGER", "id", -1),
            DBColumn("TEXT", "datetime", datetime.datetime.now().isoformat()),
        ]
        static_columns_end = [
            DBColumn("REAL", "duration", -1),
            DBColumn("REAL", "size_total", -1),
            DBColumn("REAL", "size_mp4", -1),
            DBColumn("REAL", "size_his", -1),
            DBColumn("REAL", "size_diag", -1),
            DBColumn("REAL", "size_flux", -1),
            DBColumn("TEXT", "comment", self.description),
        ]
        new_columns = (
            static_columns_start
            + [DBColumn(sql_type(val), name, val) for name, val in self.parameters.items()]
            + static_columns_end
        )
        # If the table exists, check its columns, otherwise create it
        if self.table_exists(self.exp_class):
            self.cursor.execute('PRAGMA table_info("{}")'.format(self.exp_class))
            columns = self.cursor.fetchall()
            # Check static columns
            for st_col in static_columns_start:
                column = columns.pop(0)
                col_index = column[0]
                col_name = column[1]
                col_type = column[2]
                if col_name != st_col.name or col_type != st_col.sql_type:
                    raise DatabaseError(
                        "expected column {} in the database to be {} {} but is {} {}."
                        .format(col_index, st_col.sql_type, st_col.name, col_type, col_name)
                    )
            for st_col in reversed(static_columns_end):
                column = columns.pop()
                col_index = column[0]
                col_name = column[1]
                col_type = column[2]
                if col_name != st_col.name or col_type != st_col.sql_type:
                    raise DatabaseError(
                        "expected column {} in the database to be {} {} but is {} {}."
                        .format(col_index, st_col.sql_type, st_col.name, col_type, col_name)
                    )
            # Check user-defined columns
            # TODO: be more flexible here:
            #  - allow to skip columns which are no longer needed
            #  - allow to add new columns if needed
            #  - allow to convert types
            for name, val in self.parameters.items():
                column = columns.pop(0)
                col_index = column[0]
                col_name = column[1]
                col_type = column[2]
                type_ = sql_type(val)
                if col_name != name or col_type != type_:
                    print(repr(val), sql_type(val))
                    raise InputMismatch(
                        "parameter {} of type {} does not fit into column {} "
                        "with name {} and type {}."
                        .format(name, type_, col_index, col_name, col_type)
                    )
        else:
            print(' Creating new table "{}".'.format(self.exp_class))
            column_string = ", ".join(
                ['"{}" {}'.format(col.name, col.sql_type) for col in new_columns]
            )
            self.cursor.execute(
                'CREATE TABLE "{}" ({})'.format(self.exp_class, column_string)
            )
        # Set the experiment ID if it was not defined in the experiment file
        if self.exp_id is None:
            new_entry = True
            # Get the highest index of the table or start with ID 1 if table is empty
            self.cursor.execute(
                'SELECT id from "{}" ORDER BY id DESC'.format(self.exp_class)
            )
            highest_entry = self.cursor.fetchone()
            self.exp_id = highest_entry[0] + 1 if highest_entry else 1
        else:
            # If no entry with this ID exists, create a new one
            self.cursor.execute(
                'SELECT id from "{}" WHERE id = ?'.format(self.exp_class),
                [self.exp_id],
            )
            new_entry = self.cursor.fetchone() is None
        new_columns[0] = DBColumn("INTEGER", "id", self.exp_id)
        if new_entry:
            print(' Adding new entry #{} to table "{}".'.format(self.exp_id, self.exp_class))
            # Use the question mark as a placeholder to profit from string formatting of sqlite
            value_string = ', '.join(['?'] * len(new_columns))
            self.cursor.execute(
                'INSERT INTO "{}" VALUES ({})'.format(self.exp_class, value_string),
                [sql_value(col.value) for col in new_columns],
            )
        else:
            print(' Overwriting entry #{} of table "{}".'.format(self.exp_id, self.exp_class))
            # Use the question mark as a placeholder to profit from string formatting of sqlite
            column_name_string = ", ".join('"{}" = ?'.format(col.name) for col in new_columns)
            self.cursor.execute(
                'UPDATE "{}" SET {} WHERE id = ?'.format(self.exp_class, column_name_string),
                [sql_value(col.value) for col in new_columns] + [self.exp_id],
            )
        # Save the database
        self.connection.commit()

    def finalize(self, fluid2d):
        """Save information about the completed run in the database.

        This method must be called when the simulation is finished,
        that means, after the line `f2d.loop()`.
        It writes the integration time and the sizes of the created
        output files into the database.  If a blow-up was detected, this
        is stored in the comment-field.  Furthermore, the datetime-field
        of the database entry is set to the current time."""
        if not self.connection:
            return
        # Divide every size by 1000*1000 = 1e6 to get the value in MB
        # Total size
        output_dir = os.path.dirname(fluid2d.output.hisfile)
        try:
            output_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir)]
            total_size = sum(os.path.getsize(path) for path in output_files) / 1e6
        except FileNotFoundError as e:
            print(" Error getting total size of output:", e)
            total_size = -1
        # History file
        try:
            his_size = os.path.getsize(fluid2d.output.hisfile) / 1e6
        except FileNotFoundError as e:
            print(" Error getting size of his-file:", e)
            his_size = -1
        # Diagnostics file
        try:
            diag_size = os.path.getsize(fluid2d.output.diagfile) / 1e6
        except FileNotFoundError:
            print(" Error getting size of diag-file:", e)
            diag_size = -1
        # MP4 file
        if fluid2d.plot_interactive and hasattr(fluid2d.plotting, 'mp4file'):
            try:
                mp4_size = os.path.getsize(fluid2d.plotting.mp4file) / 1e6
            except FileNotFoundError:
                print(" Error getting size of mp4-file:", e)
                mp4_size = -1
        else:
            mp4_size = 0
        # Flux file
        if fluid2d.diag_fluxes:
            try:
                flux_size = os.path.getsize(fluid2d.output.flxfile) / 1e6
            except FileNotFoundError:
                print(" Error getting size of flux-file:", e)
                flux_size = -1
        else:
            flux_size = 0
        # Check for blow-up
        if hasattr(fluid2d, "blow_up") and fluid2d.blow_up:
            comment = "Blow-up! " + self.description
        else:
            comment = self.description
        # Update and save the database
        self.cursor.execute(
            """UPDATE "{}" SET
                duration = ?,
                datetime = ?,
                size_total = ?,
                size_mp4 = ?,
                size_his = ?,
                size_diag = ?,
                size_flux = ?,
                comment = ?
            WHERE id = ?""".format(self.exp_class),
            (
                round(fluid2d.t, 2),
                datetime.datetime.now().isoformat(),
                round(total_size, 3),
                round(mp4_size, 3),
                round(his_size, 3),
                round(diag_size, 3),
                round(flux_size, 3),
                comment,
                self.exp_id,
            ),
        )
        self.connection.commit()

    def connect(self, data_dir: str):
        print("-"*50)
        if data_dir.startswith("~"):
            data_dir = os.path.expanduser(data_dir)
        if not os.path.isdir(data_dir):
            print(" Creating directory {}.".format(data_dir))
            os.makedirs(data_dir)
        dbpath = os.path.join(data_dir, "experiments.db")
        print(" Opening database {}.".format(dbpath))
        self.connection = dbsys.connect(dbpath)
        self.cursor = self.connection.cursor()

    def table_exists(self, name: str):
        """Check if table with given name exists in the connected database."""
        # Get all tables
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        # This gives a list like that: [('table1',), ('table2',)]
        for table in self.cursor.fetchall():
            if table[0] == name:
                return True
        return False


def parse_experiment_file(path: str):
    """Parse the given experiment file.

    An experiment file
     - must provide a name,
     - can provide a description,
     - can provide parameters with values.

    The name is written in a line starting with "Name:".  It must be a
    valid string to be used as a filename;  in particular, it must not
    contain a slash (/) and it must not contain a quotation mark (").
    It is advised to use underscores instead of spaces in the name.

    The description begins in or after a line starting with
    "Description:" and goes until the beginning of the parameters.
    If no parameters are defined below the description, it goes until
    the end of the file.

    The parameters follow after a line starting with "Parameters:".
    Every parameter is written in its own line.  This line contains a
    name and one value or several values, each separated by one or
    several whitespaces.  The name must not contain any whitespace
    characters (otherwise it is not recognised as such) and must not
    contain a quotation mark.  The name must not be in the list of
    reserved column names.  If the value is True, False, an integer or a
    floation point number, it is interpreted as the corresponding Python
    type, otherwise it is considered a string.  Quotation marks in the
    value are taken literally and are not interpreted.

    Everything after a #-symbol is considered a comment and ignored.

    TODO:
     - allow to use quotation marks for multi-word strings as values
    """
    with open(path) as f:
        experiment_lines = f.readlines()
    name = ""
    id_ = None
    description_lines = []
    param_name_list = []
    param_values_list = []
    reading_params = False
    reading_description = False
    for line in experiment_lines:
        # Remove comments and whitespace at the beginning and the end
        line = line.split("#")[0].strip()
        if not line:
            # Skip empty lines except in the description
            if reading_description:
                description_lines.append(line)
            else:
                continue
        elif line.lower().startswith("name:"):
            if not name:
                name = line[5:].strip()
                if '"' in name:
                    raise ExpFileError("Name must not contain a quotation mark.")
            else:
                raise ExpFileError(
                    "Name defined more than once."
                )
        elif line.lower().startswith("id:"):
            if id_ is None:
                value = line[3:].strip()
                if value != "":
                    try:
                        id_ = int(value)
                    except ValueError:
                        raise ExpFileError("ID is not an integer.")
            else:
                raise ExpFileError("ID defined more than once.")
        elif line.lower().startswith("description:"):
            reading_description = True
            reading_params = False
            description_lines.append(line[12:].lstrip())
        elif line.lower().startswith("parameters:"):
            reading_description = False
            reading_params = True
        elif reading_description:
            description_lines.append(line)
        elif reading_params:
            # Parse the parameter name
            param_name = line.split()[0]
            if param_name in RESERVED_COLUMN_NAMES:
                raise ExpFileError(
                    "reserved name {} must not be used as a parameter.".format(param_name)
                )
            if '"' in param_name:
                raise ExpFileError("parameter name must not contain a quotation mark.")
            # Parse the value(s) of the parameter
            param_val_text = line[len(param_name):].lstrip()
            if not param_val_text:
                raise ExpFileError("no value given for parameter {}.".format(param_name))
            # TODO: allow multi-word strings containing quotation marks
            param_values = [cast_string(val) for val in param_val_text.split()]
            param_name_list.append(param_name)
            param_values_list.append(param_values)
        else:
            raise ExpFileError("unexpected line:\n{!r}".format(line))
    if not name:
        raise ExpFileError("Name must be defined in every experiment file.")
    return name, id_, "\n".join(description_lines).strip(), param_name_list, param_values_list


def sql_type(value):
    if isinstance(value, bool):
        return "TEXT"
    elif isinstance(value, int):
        return "INTEGER"
    elif isinstance(value, float):
        return "REAL"
    return "TEXT"


def sql_value(value):
    if isinstance(value, bool):
        return str(value)
    return value


def cast_string(value: str):
    """Cast into the most specialised Python type possible.

    This method can recognize "True", "False", integers and floats,
    everything else is treated as a string and returned unchanged."""
    if value == "True":
        return True
    if value == "False":
        return False
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value
