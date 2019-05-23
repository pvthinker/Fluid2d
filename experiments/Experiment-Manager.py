# Markus Reinert, May 2019
#
# Command-line interface for the fluid2d Experiment Management System (EMS).
#
# This programme can be started directly with Python if fluid2d is activated.
# Otherwise, the path to the folder with the experiment-files (that is the
# value of param.datadir) must be specified as a command-line argument.

import os
import sys
import cmd
import sqlite3 as dbsys
import subprocess


# Command to open mp4-files
MP4_PLAYER = "mplayer"
# Command to open NetCDF (his or diag) files
NETCDF_VIEWER = "ncview"

# Settings for the output of a table
COMMENT_MAX_LENGTH = 21
FLOAT_FORMAT = "{:.4f}"
LIMITER = "  "

# Hide the following information in the table
# They can be activated with "enable" during runtime.
# More information can be hidden with "disable".
hidden_information = {
    'size_diag',
    'size_flux',
    'datetime_seconds',
}

class EMDBConnection:
    """Experiment Management Database Connection"""

    def __init__(self, dbpath: str):
        # Create a connection to the given database
        print("-"*50)
        print("Opening database {}.".format(dbpath))
        self.connection = dbsys.connect(dbpath)
        cursor = self.connection.cursor()

        # Get all tables of the database
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        # This gives a list like that: [('table1',), ('table2',)]
        self.tables = [EMDBTable(cursor, t[0]) for t in cursor.fetchall()]

    def __del__(self):
        print("Closing database.")
        print("-"*50)
        self.connection.close()

    def get_table_overview(self):
        if self.tables:
            text = "Experiments in database:"
            for table in self.tables:
                text += (
                    "\n - {}: {} experiments, {} columns"
                    .format(table.name, table.get_length(), len(table.columns))
                )
        else:
            text = "No experiments in database."
        return text

    def show_full_tables(self):
        print("-"*50)
        for table in self.tables:
            print(table)
            print("-"*50)

    def show_table(self, name):
        for table in self.tables:
            if table.name == name:
                print("-"*50)
                print(table)
                print("-"*50)
                return True
        return False

    def show_filtered_table(self, table_name, statement):
        for table in self.tables:
            if table.name == table_name or table_name == "":
                print("-"*50)
                table.print_selection(statement)
        print("-"*50)

    def show_sorted_table(self, table_name, statement):
        for table in self.tables:
            if table.name == table_name or table_name == "":
                print("-"*50)
                table.print_sorted(statement)
        print("-"*50)

    def is_valid_column(self, name):
        for table in self.tables:
            for n, t in table.columns:
                if name == n:
                    return True
        return False


class EMDBTable:
    """Experiment Management Database Table"""

    def __init__(self, cursor: dbsys.Cursor, name: str):
        self.name = str(name)
        self.c = cursor
        # Get columns
        self.c.execute('PRAGMA table_info("{}")'.format(self.name))
        self.columns = [(column[1], column[2]) for column in self.c.fetchall()]
        # column[0]: index from 0, column[1]: name, column[2]: type

    def __str__(self):
        # Get entries
        self.c.execute('SELECT * from "{}"'.format(self.name))
        return string_format_table(self.name, self.columns, self.c.fetchall())

    def get_length(self):
        self.c.execute('SELECT Count(*) FROM "{}"'.format(self.name))
        return self.c.fetchone()[0]

    def entry_exists(self, id_):
        self.c.execute('SELECT id FROM "{}" WHERE id = ?'.format(self.name), (id_,))
        if self.c.fetchone() is not None:
            return True
        return False

    def print_selection(self, statement):
        try:
            self.c.execute('SELECT * FROM "{}" WHERE {}'
                           .format(self.name, statement))
        except dbsys.OperationalError as e:
            print('SQL error for experiment "{}":'.format(self.name), e)
        else:
            print(string_format_table(self.name, self.columns, self.c.fetchall()))

    def print_sorted(self, statement):
        try:
            self.c.execute('SELECT * FROM "{}" ORDER BY {}'
                           .format(self.name, statement))
        except dbsys.OperationalError as e:
            print('SQL error for experiment "{}":'.format(self.name), e)
        else:
            print(string_format_table(self.name, self.columns, self.c.fetchall()))


# https://docs.python.org/3/library/cmd.html
class EMShell(cmd.Cmd):
    """Experiment Management (System) Shell"""
    
    intro = (
        "-" * 50
        + "\nType help or ? to list available commands."
        + "\nType exit or Ctrl+D or Ctrl+C to exit."
        + "\nPress Tab-key for auto-completion of commands, table names or column names."
        + "\n"
    )
    prompt = "(EMS) "

    def __init__(self, experiments_dir: str):
        super().__init__()
        print("-"*50)
        print("Fluid2d Experiment Management System (EMS)")
        self.exp_dir = experiments_dir
        self.con = EMDBConnection(os.path.join(self.exp_dir, "experiments.db"))
        self.intro += "\n" + self.con.get_table_overview() + "\n"
        self.selected_table = ""
        self.silent_mode = True

    ### Functionality to MODIFY how the programme acts
    def do_verbose(self, params):
        if params == "" or params.lower() == "on":
            self.silent_mode = False
        elif params.lower() == "off":
            self.silent_mode = True
        else:
            print('Unknown parameter.  Please use "verbose on" or "verbose off".')

    def complete_verbose(self, text, line, begidx, endidx):
        if text == "" or text == "o":
            return ["on", "off"]
        if text == "of":
            return ["off",]

    def help_verbose(self):
        print(
            "Toggle between verbose- and silent-mode.\n"
            'Use "verbose on" or "verbose" to see the output of external programmes '
            'started from this shell.\n'
            'Use "verbose off" to hide all output of external programmes (default).\n'
            'No error message is displayed in silent-mode when the opening of a file fails.\n'
            'In any case, external programmes are started in the background, so the shell can '
            'still be used, even if the input prompt is polluted.'
        )

    def do_enable(self, params):
        global hidden_information
        if params == "all":
            hidden_information.clear()
        else:
            for param in params.split():
                if param == "size":
                    hidden_information.difference_update(
                        {'size_diag', 'size_flux', 'size_his', 'size_mp4', 'size_total'}
                    )
                else:
                    hidden_information.discard(param)

    def complete_enable(self, text, line, begidx, endidx):
        parameters = hidden_information.copy()
        if len(parameters) > 0:
            parameters.add('all')
        if not parameters.isdisjoint({'size_diag', 'size_flux', 'size_his',
                                      'size_mp4', 'size_total'}):
            parameters.add('size')
        completions = []
        for p in parameters:
            if p.startswith(text):
                completions.append(p)
        return completions

    def help_enable(self):
        print(
            'Make hidden information in the experiment table visible.\n'
            'See the help of "disable" for further explanations.'
        )

    def do_disable(self, params):
        global hidden_information
        for param in params.split():
            if param == "size":
                # Short notation for the union of two sets
                hidden_information |= {'size_diag', 'size_flux', 'size_his',
                                       'size_mp4', 'size_total'}
            elif (param == "datetime_year"
                  or param == "datetime_seconds"
                  or self.con.is_valid_column(param)
            ):
                hidden_information.add(param)
            else:
                print("Unknown argument:", param)

    def complete_disable(self, text, line, begidx, endidx):
        parameters = [
            'datetime_year',
            'datetime_seconds',
        ]
        # If not all size columns are hidden, add shorthand for all sizes
        if not hidden_information.issuperset({'size_diag', 'size_flux', 'size_his',
                                              'size_mp4', 'size_total'}):
            parameters += ['size']
        # Add columns of selected table to the list of auto-completions
        for table in self.con.tables:
            if self.selected_table == "" or self.selected_table == table.name:
                parameters += [n for n,t in table.columns]
        completions = []
        for p in parameters:
            if p.startswith(text) and p not in hidden_information:
                completions.append(p)
        return completions

    def help_disable(self):
        print(
            'Hide unneeded information in the experiment table.\n'
            'Every parameter/column can be hidden, for example:\n'
            ' - disable comment\n'
            'Furthermore, the year and number of seconds can be hidden from the datetime '
            'column by using "disable datetime_year" and "disable datetime_seconds".\n'
            'Multiple parameters can be specified at once, for example:\n'
            ' - disable size_his size_diag size_flux\n'
            'To hide all file sizes, use "disable size".\n'
            'To show hidden parameters again, use "enable".\n'
            'In addition to the behaviour described here, it also supports the command '
            '"enable all".'
        )

    ### Functionality to SHOW the content of the database
    def do_list(self, params):
        """List all experiment classes (tables) in the database."""
        print(self.con.get_table_overview())

    def do_show(self, table_name):
        """Show the content of a table.

        If no table name is specified or selected, the content of every table is shown."""
        if table_name:
            # Try to open the specified table.
            if not self.con.show_table(table_name):
                # If it fails, print a message.
                print('Unknown experiment: "{}"'.format(name))
        else:
            if self.selected_table:
                self.con.show_table(self.selected_table)
            else:
                # No table name given and no table selected
                self.con.show_full_tables()

    def complete_show(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def help_show(self):
        print(
            "Show all information in the database about a class of experiments.\n"
            "If no experiment is specified and no experiment is selected, "
            "information about all the experiments is shown."
        )

    def do_filter(self, params):
        self.con.show_filtered_table(self.selected_table, params)

    def complete_filter(self, text, line, begidx, endidx):
        return self.column_name_completion(self.selected_table, text)

    def help_filter(self):
        print(
            'Filter experiments by the given value of a given parameter.\n'
            'Any valid SQLite WHERE-statement can be used as a filter, for example:\n'
            ' - filter intensity <= 0.2\n'
            ' - filter slope = 0.5\n'
            ' - filter diffusion = "True"\n'
            ' - filter perturbation = "gauss" AND duration > 20\n'
            'It is necessary to put the value for string or boolean arguments in quotation '
            'marks like this: "True", "False", "Some text".\n'
            'Quotation marks are also necessary if the name of the parameter contains '
            'whitespace.\n'
            'To sort the experiments, use "sort".\n'
            'To filter and sort the experiments, use SQLite syntax.\n'
            'Examples:\n'
            ' - filter intensity > 0.1 ORDER BY intensity\n'
            ' - filter perturbation != "gauss" ORDER BY size_total DESC'
        )

    def do_sort(self, params):
        self.con.show_sorted_table(self.selected_table, params)

    def complete_sort(self, text, line, begidx, endidx):
        return self.column_name_completion(self.selected_table, text)

    def help_sort(self):
        print(
            'Sort the experiments by the value of a given parameter.\n'
            'To invert the order, add "desc" (descending) to the command.\n'
            'Example usages:\n'
            ' - sort intensity: show experiments with lowest intensity on top of the table.\n'
            ' - sort size_total desc: show experiments with biggest file size on top.\n'
            'Quotation marks are necessary if the name of the parameter contains whitespace.\n'
            'To filter and sort the experiments, see the help of "filter".'
        )

    ### Functionality to OPEN experiment files
    def do_open_mp4(self, params):
        """Open the mp4-file for an experiment specified by its name and ID."""
        expname = self.parse_params_to_experiment(params)
        if not expname:
            return
        dir_ = os.path.join(self.exp_dir, expname)
        for f in os.listdir(dir_):
            if f.endswith(".mp4") and f.startswith(expname):
                break
        else:
            print("No mp4-file found in folder:", dir_)    
        path = os.path.join(self.exp_dir, expname, f)
        self.open_file(MP4_PLAYER, path)

    def complete_open_mp4(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def help_open_mp4(self):
        print("Open the mp4-file for an experiment specified by its name and ID.\n"
              "It is not possible to interact with the video player via input in the shell, "
              "for example with mplayer.\n"
              "User input via the graphical interface is not affected by this."
        )

    def do_open_his(self, params):
        """Open the his-file for an experiment specified by its name and ID."""
        expname = self.parse_params_to_experiment(params)
        if not expname:
            return
        path = os.path.join(self.exp_dir, expname, expname + "_his.nc")
        if not os.path.isfile(path):
            print("File does not exist:", path)
            return
        self.open_file(NETCDF_VIEWER, path)

    def complete_open_his(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def do_open_diag(self, params):
        """Open the diag-file for an experiment specified by its name and ID."""
        expname = self.parse_params_to_experiment(params)
        if not expname:
            return
        path = os.path.join(self.exp_dir, expname, expname + "_diag.nc")
        if not os.path.isfile(path):
            print("File does not exist:", path)
            return
        self.open_file(NETCDF_VIEWER, path)

    def complete_open_diag(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    ### Functionality to SELECT a specific table
    def do_select(self, params):
        if params == "":
            self.prompt = "(EMS) "
            self.selected_table = params
        elif self.check_table_exists(params):
            self.prompt = "({}) ".format(params)
            self.selected_table = params
        else:
            print('Unknown experiment: "{}"'.format(params))

    def complete_select(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def help_select(self):
        print(
            'Select an experiment by specifing its name.\n'
            'When an experiment is selected, every operation is automatically '
            'executed for this experiment, except if specified differently.\n'
            'To unselect again, use the "select"-command with no argument or press Ctrl+D.'
        )

    ### Functionality to QUIT the program
    def do_exit(self, params):
        return True

    def do_EOF(self, params):
        print("")
        if self.selected_table:
            self.prompt = "(EMS) "
            self.selected_table = ""
        else:
            return True

    ### Behaviour for empty input
    def emptyline(self):
        pass

    ### Helper functions
    def table_name_completion(self, text):
        completions = []
        for table in self.con.tables:
            if table.name.startswith(text):
                completions.append(table.name)
        return completions

    def column_name_completion(self, table_name, text):
        completions = []
        for table in self.con.tables:
            if table.name == table_name or table_name == "":
                for column_name, column_type in table.columns:
                    if column_name.startswith(text):
                        completions.append(column_name)
        return completions

    def check_table_exists(self, name):
        return name in [table.name for table in self.con.tables]

    def parse_params_to_experiment(self, params):
        params = params.split(" ")
        # Check for correct parameter input
        if len(params) < 2:
            if self.selected_table:
                if len(params) == 1:
                    name = self.selected_table
                else:
                    print("Exactly 1 ID must be specified.")
                    return
            else:
                print("Name and ID of experiment must be specified.")
                return
        else:
            # Extract name from input
            name = " ".join(params[:-1]).strip()
        # Check name
        for table in self.con.tables:
            if table.name == name:
                break
        else:
            print('Unknown experiment: "{}"'.format(name))
            return
        # Extract and check ID
        try:
            id_ = int(params[-1])
        except ValueError:
            print('Last argument "{}" is not a valid ID.'.format(params[-1]))
            return
        if not table.entry_exists(id_):
            print("No entry exists in table {} with ID {}.".format(name, id_))
            return
        # Return name of experiment folder
        return "{}_{:03}".format(name, id_)

    def open_file(self, command, path):
        if self.silent_mode:
            print("Opening file {} with {} in silent-mode.".format(path, command))
            subprocess.Popen(
                [command, path],
                # Disable standard input via the shell, for example with mplayer.
                stdin=subprocess.DEVNULL,
                # Throw away output and error messages.
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        else:
            print("Opening file {} with {} in verbose-mode.".format(path, command))
            subprocess.Popen(
                [command, path],
                # Disable standard input via the shell, for example with mplayer.
                stdin=subprocess.DEVNULL,
            )


def string_format_table(table_name, columns, rows):
    # Convert rows to list of lists (since fetchall() returns a list of tuples)
    rows = [list(row) for row in rows]
    columns = columns.copy()
    # Remove ignored columns
    indices_to_remove = []
    for i, (n, t) in enumerate(columns):
        if n in hidden_information:
            indices_to_remove.append(i)
    for i in sorted(indices_to_remove, reverse=True):
        columns.pop(i)
        for row in rows:
            row.pop(i)
    # To calculate the total size of the files associated with the table
    table_size = 0
    # Get necessary length for each column and process table
    lengths = [len(n) for n, t in columns]
    for row in rows:
        for i, val in enumerate(row):
            # Check name of column
            if columns[i][0] == "comment":
                # Cut comments which are too long and end them with an ellipsis
                if len(val) > COMMENT_MAX_LENGTH:
                    val = val[:COMMENT_MAX_LENGTH-1] + "…"
                    row[i] = val
            elif columns[i][0] == "datetime":
                # Cut unnecessary parts of the date and time
                if "datetime_year" in hidden_information:
                    val = val[5:]
                if "datetime_seconds" in hidden_information:
                    val = val[:-10]
                row[i] = val.replace('T', ',')
            elif columns[i][0] == "size_total":
                if val > 0:
                    table_size += val

            # Check type of column and adopt
            if columns[i][1] == "TEXT" or columns[i][1] == "INTEGER":
                lengths[i] = max(lengths[i], len(str(val)))
            elif columns[i][0] in ["size_diag", "size_flux", "size_his",
                                   "size_mp4", "size_total"]:
                lengths[i] = max(lengths[i], len("{:.3f}".format(val)))
            elif columns[i][1] == "REAL":
                lengths[i] = max(lengths[i], len(FLOAT_FORMAT.format(val)))
            else:
                # This is an unexpected situation, which probably means that
                # sqlite3 does not work as it was, when this script was written.
                raise Exception(
                    "unknown type {} of column {} in table {}."
                    .format(*columns[i], table_name)
                )
    # Create top line of the text to be returned
    text = "Experiment: " + table_name
    if "size_total" not in hidden_information:
        text += " ({:.3f} MB)".format(table_size)
    if len(indices_to_remove) == 1:
        text += " (1 parameter hidden)"
    elif len(indices_to_remove) > 1:
        text += " ({} parameters hidden)".format(len(indices_to_remove))
    text += "\n"
    # Add column name
    text += LIMITER.join([
        ("{:^" + str(l) + "}").format(n)
        for (n, t), l in zip(columns, lengths)
    ]) + "\n"
    format_strings = []
    for (n, t), l in zip(columns, lengths):
        # Numbers right justified,
        # comments left justified,
        # text centered
        if n in ["size_diag", "size_flux", "size_his", "size_mp4", "size_total"]:
            format_strings.append("{:>" + str(l) + ".3f}")
        elif t == "REAL":
            format_strings.append(FLOAT_FORMAT.replace(":", ":>"+str(l)))
        elif t == "INTEGER":
            format_strings.append("{:>" + str(l) + "}")
        elif n == "comment":
            format_strings.append("{:" + str(l) + "}")
        else:
            format_strings.append("{:^" + str(l) + "}")
    for row in rows:
        text_cols = [f_str.format(val) for f_str, val in zip(format_strings, row)]
        text += LIMITER.join(text_cols) + "\n"
    return text.strip()


if len(sys.argv) == 1:
    try:
        from param import Param
    except ModuleNotFoundError:
        raise Exception(
            "When fluid2d is not available, this programme has to be started "
            "with the experiments-folder as argument."
        )
    param = Param(None)  # it is not necessary to specify a defaultfile for Param
    datadir = param.datadir
    if datadir.startswith("~"):
        datadir = os.path.expanduser(param.datadir)
    del param
else:
    datadir = sys.argv[1]
ems_cli = EMShell(datadir)
ems_cli.cmdloop()
