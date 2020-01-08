"""EMShell: Command-line interface for the EMS of fluid2d

EMS is the Experiment Management System of fluid2d, a powerful and
convenient way to handle big sets of experiments, cf. core/ems.py.
The EMShell is a command line interface to access, inspect, modify
and analyse the database and its entries.

To start this programme, activate fluid2d, then run this script with
Python (version 3.6 or newer).  Alternatively, without the need to
activate fluid2d, specify the path to the experiment folder as a
command-line argument when launching this script with Python.

This code uses f-strings, which were introduced in Python 3.6:
  https://docs.python.org/3/whatsnew/3.6.html#whatsnew36-pep498
Unfortunately, they create a SyntaxError in older Python versions.

It is possible to access the experiment database from multiple
processes at the same time, so one can add new entries to the database
with the EMS of fluid2d while looking at the database and performing
data analysis in the EMShell;  the EMShell shows automatically the
updated version of the database.  However, this works less well if the
database is accessed by processes on different computers.  In this
case it can be necessary to restart the EMShell to see changes in the
database.  It is always necessary to restart the EMShell when a new
class of experiments was added to the database.

Author: Markus Reinert, May/June 2019
"""

import os
import re
import sys
import cmd
import time
import numpy as np
import shutil
import sqlite3 as dbsys
import datetime
import readline
import itertools
import subprocess
import matplotlib.pyplot as plt
# Local imports
import EMShellExtensions as EMExt


# Command to open mp4-files
MP4_PLAYER = "mplayer"
# Command to open NetCDF (his or diag) files
NETCDF_VIEWER = "ncview"

### Extensions
# Make new shell extensions available by adding them to this dictionary.
# The key is the name under which the function is called in the shell,
# which must not contain any whitespace.
# The value is the function, which takes as only argument the path to the
# his-file of an experiment and can return the calculated value.
extra_tools = {
    "wavenumber": EMExt.get_strongest_wavenumber_y,
    "wavelength": lambda hisname: 1/EMExt.get_strongest_wavenumber_y(hisname),
}

### Settings for the output of a table
# Maximal length of comments (possible values: "AUTO", "FULL" or a positive integer)
# This setting is ignored if display_all is True.
COMMENT_MAX_LENGTH = "AUTO"
# Format-string for real numbers
FLOAT_FORMAT = "{:.4f}"
# Format-string for file sizes (MB) or disk space (GiB)
# Two decimals are used instead of three, to avoid confusion between the dot as
# a decimal separator and a delimiter after 1000.
SIZE_FORMAT = "{:.2f}"
# Symbol to separate two columns of the table
LIMITER = "  "
# Symbol of the linebreak and its replacement in comments
# This is ignored if display_all is True or COMMENT_MAX_LENGTH is "FULL".
LINEBREAK_REPLACE = ("\n", "|")
# Show date and time in ISO-format or easy-to-read-format
# The easy-to-read-format does not include seconds, independent whether seconds
# are hidden or not.  This setting is ignored if display_all is True.
ISO_DATETIME = False

### Settings to print tables in colour
# Set COLOURS to False or None to disable colours (this does not disable highlighting).
# Otherwise, COLOURS must be list, of which each element describes the colours used in one row.
# The colours of every row are described by a list of colour codes,
# cf. https://en.wikipedia.org/wiki/ANSI_escape_code#Colors.
# An arbitrary number of colours can be specified.
# Apart from background colours, also text colours and font styles can be specified.
# Example to fill rows alternately with white and default colour:
# COLOURS = (
#     ("\033[47m",),
#     ("\033[49m",),
# )
# Example to fill columns alternately with white and default colour:
# COLOURS = (
#     ("\033[47m", "\033[49m",),
# )
# Example for a check-pattern:
# COLOURS = (
#     ("\033[47m", "\033[49m",),
#     ("\033[49m", "\033[47m",),
# )
# Example for alternating background colours in rows and alternating font colours in columns:
# COLOURS = (
#    ("\033[31;47m", "\033[39;47m",),
#    ("\033[31;49m", "\033[39;49m",),
# )
COLOURS = (
    ("\033[107m",),
    ("\033[49m",),
)
# Colour-code to highlight columns (used by "compare -h")
COLOURS_HIGHLIGHT = "\033[103m"
# Colour-code to reset to default colour and default font
COLOURS_END = "\033[39;49m"


### Language settings (currently only for dates in easy-to-read-format)
# Definition of languages
class English:
    JUST_NOW = "just now"
    AGO_MINUTES = "{} min ago"
    AGO_HOURS = "{}:{:02} h ago"
    YESTERDAY = "yesterday"
    FUTURE = "in the future"

class French:
    JUST_NOW = "maintenant"
    AGO_MINUTES = "il y a {} min"
    AGO_HOURS = "il y a {}h{:02}"
    YESTERDAY = "hier"
    FUTURE = "à l'avenir"

class German:
    JUST_NOW = "gerade eben"
    AGO_MINUTES = "vor {} Min."
    AGO_HOURS = "vor {}:{:02} Std."
    YESTERDAY = "gestern"
    FUTURE = "in der Zukunft"

# Selection of a language
LANG = English


### Global variables modifiable during runtime
# Hide the following information in the table
# They can be activated with the command "enable" during runtime.
# More information can be hidden with the command "disable".
hidden_information = {
    'size_diag',
    'size_flux',
    'datetime_seconds',
}

# Temporarily show all information, including full comments
# This is set to True by the argument "-v" to commands which print tables.
display_all = False


class EMDBConnection:
    """Experiment Management Database Connection"""

    def __init__(self, dbpath: str):
        self.connection = None
        if not os.path.isfile(dbpath):
            raise FileNotFoundError(f"Database file {dbpath} does not exist.")
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
        if self.connection:
            print("Closing database.")
            print("-"*50)
            self.connection.close()

    def save_database(self):
        self.connection.commit()

    def get_table_overview(self):
        if self.tables:
            text = "Experiments in database:"
            for table in self.tables:
                n_entries = len(table)
                n_columns = len(table.columns)
                text += f"\n - {table.name}: "
                text += f"{n_entries} experiment{'s' if n_entries != 1 else ''}, "
                text += f"{n_columns} columns"
        else:
            text = "No experiments in database."
        return text

    def get_data(self, table_name, column_name, condition=""):
        for table in self.tables:
            if table_name == table.name:
                return table.get_data(column_name, condition)
        print(f'No table with name {table_name} found.')
        return []

    def get_column_names(self, table_name):
        for table in self.tables:
            if table.name == table_name:
                return [n for n,t in table.columns]

    def get_table_length(self, table_name):
        for table in self.tables:
            if table.name == table_name:
                return len(table)

    def get_latest_entry(self, table_name):
        for table in self.tables:
            if table.name == table_name:
                return table.get_latest_entry()

    def set_comment(self, table_name, id_, new_comment):
        for table in self.tables:
            if table.name == table_name:
                return table.set_comment(id_, new_comment)
        print(f'Unknown experiment: "{table_name}"')
        return False

    def show_all_tables(self):
        for table in self.tables:
            print("-"*50)
            print(table)

    def show_table(self, name):
        print("-"*50)
        for table in self.tables:
            if table.name == name:
                print(table)
                return
        print('Unknown experiment: "{}"'.format(name))

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

    def show_comparison(self, table_name, ids, highlight=False):
        for table in self.tables:
            if table.name == table_name:
                print("-"*50)
                table.print_comparison(ids, highlight)
                print("-"*50)
                return

    def is_valid_column(self, column_name, table_name=""):
        """Check whether a column with the given name exists.

        If no table_name is specified, check whether any table has such
        a column, otherwise look only in the table with the given name."""
        for table in self.tables:
            if table.name == table_name or table_name == "":
                for n, t in table.columns:
                    if column_name == n:
                        return True
        return False

    def table_exists(self, table_name):
        for table in self.tables:
            if table.name == table_name:
                return True
        return False

    def entry_exists(self, table_name, id_):
        for table in self.tables:
            if table.name == table_name:
                return table.entry_exists(id_)

    def delete_entry(self, table_name, id_):
        for table in self.tables:
            if table.name == table_name:
                table.delete_entry(id_)
                return

    def delete_table(self, table_name):
        self.connection.execute(f'DROP TABLE "{table_name}"')


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

    def __len__(self):
        self.c.execute('SELECT Count(*) FROM "{}"'.format(self.name))
        return self.c.fetchone()[0]

    def get_data(self, column_name, condition=""):
        try:
            self.c.execute(
                f'SELECT {column_name} FROM "{self.name}" WHERE {condition}'
                if condition else
                f'SELECT {column_name} FROM "{self.name}"'
            )
        except dbsys.OperationalError as e:
            print(f'SQL error for experiment "{self.name}":', e)
            return []
        else:
            return [e[0] for e in self.c.fetchall()]

    def get_latest_entry(self):
        self.c.execute('SELECT id from "{}" ORDER BY datetime DESC'.format(self.name))
        result = self.c.fetchone()
        self.c.fetchall()  # otherwise the database stays locked
        return result[0] if result else None

    def entry_exists(self, id_):
        self.c.execute('SELECT id FROM "{}" WHERE id = ?'.format(self.name), (id_,))
        if self.c.fetchone() is not None:
            return True
        return False

    def delete_entry(self, id_):
        self.c.execute('DELETE FROM "{}" WHERE id = ?'.format(self.name), (id_,))

    def print_selection(self, statement):
        try:
            self.c.execute('SELECT * FROM "{}" WHERE {}'.format(self.name, statement))
        except dbsys.OperationalError as e:
            print('SQL error for experiment "{}":'.format(self.name), e)
        else:
            print(string_format_table(self.name, self.columns, self.c.fetchall()))

    def print_sorted(self, statement):
        try:
            self.c.execute('SELECT * FROM "{}" ORDER BY {}'.format(self.name, statement))
        except dbsys.OperationalError as e:
            print('SQL error for experiment "{}":'.format(self.name), e)
        else:
            print(string_format_table(self.name, self.columns, self.c.fetchall()))

    def print_comparison(self, ids, highlight=False):
        try:
            self.c.execute(
                'SELECT * FROM "{}" WHERE id IN {}'.format(self.name, tuple(ids))
                if ids else
                'SELECT * FROM "{}"'.format(self.name)
            )
        except dbsys.OperationalError as e:
            print('SQL error for experiment "{}":'.format(self.name), e)
            return
        full_entries = self.c.fetchall()
        different_columns = []
        for i, row in enumerate(full_entries):
            if i == 0:
                continue
            for c, value in enumerate(row):
                if c not in different_columns and value != full_entries[i-1][c]:
                    different_columns.append(c)
        if not highlight:
            # Print only the columns of the table with differences
            print(string_format_table(
                self.name,
                [col for c, col in enumerate(self.columns) if c in different_columns],
                [
                    [val for c, val in enumerate(row) if c in different_columns]
                    for row in full_entries
                ]
            ))
        else:
            # Print the full table and highlight the columns with differences
            print(string_format_table(
                self.name, self.columns, full_entries,
                [col[0] for c, col in enumerate(self.columns) if c in different_columns]
            ))

    def set_comment(self, id_, new_comment):
        try:
            self.c.execute(
                'UPDATE "{}" SET comment = ? WHERE id = ?'.format(self.name),
                (new_comment, id_,)
            )
        except dbsys.OperationalError as e:
            print('SQL error for experiment "{}":'.format(self.name), e)
            return False
        else:
            return True


# https://docs.python.org/3/library/cmd.html
class EMShell(cmd.Cmd):
    """Experiment Management (System) Shell"""

    intro = (
        "-" * 50
        + "\nType help or ? to list available commands."
        + "\nType exit or Ctrl+D or Ctrl+C to exit."
        + "\nPress Tab-key for auto-completion of commands, table names, etc."
        + "\n"
    )
    prompt = "(EMS) "
    ruler = "-"

    def __init__(self, experiments_dir: str):
        self.initialized = False
        super().__init__()
        print("-"*50)
        print("Fluid2d Experiment Management System (EMS)")
        self.exp_dir = experiments_dir
        self.con = EMDBConnection(os.path.join(self.exp_dir, "experiments.db"))
        self.intro += "\n" + self.con.get_table_overview() + "\n"
        self.selected_table = ""
        # Settings for saving the command history
        self.command_history_file = os.path.join(self.exp_dir, ".emshell_history")
        readline.set_history_length(1000)
        # Load previously saved command history
        if os.path.exists(self.command_history_file):
            readline.read_history_file(self.command_history_file)
        self.initialized = True

    def __del__(self):
        if self.initialized:
            print("Saving command history.")
            readline.write_history_file(self.command_history_file)

    ### Functionality to MODIFY how the programme acts
    def do_enable(self, params):
        global hidden_information
        if not params:
            print("Please specify information to enable.")
        elif params == "all":
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
        return [p for p in parameters if p.startswith(text)]

    def help_enable(self):
        print(
            'Make hidden information in the experiment table visible.\n'
            'See the help of "disable" for further explanations.'
        )

    def do_disable(self, params):
        global hidden_information
        if not params:
            print("Please specify information to disable.")
            return
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
        return [p for p in parameters if p.startswith(text) and p not in hidden_information]

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

    def do_show(self, params):
        """Show the content of a table."""
        global display_all
        params = set(params.split())
        if "-v" in params:
            display_all = True
            params.remove('-v')
        else:
            display_all = False
        if len(params) == 0:
            # No table name given
            if self.selected_table:
                self.con.show_table(self.selected_table)
            else:
                self.con.show_all_tables()
        else:
            for table_name in sorted(params):
                # Try to open the specified table.
                self.con.show_table(table_name)
        print("-"*50)

    def complete_show(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def help_show(self):
        print("""> show [name(s) of experiment class(es)] [-v]
    Show all entries in the database for one or more classes of experiments.
    Specify as parameter the name of the experiment class to show.  Multiple
    names may be specified.  If no name is specified, then the experiments of
    the currently selected class are shown (cf. "select").  If no name name is
    specified and no experiment selected, then all entries are shown.
    Use the commands "enable" and "disable" to specify which information
    (i.e., which column) is displayed.  To display all information instead, use
    "show" with the parameter "-v".

    See also: filter, sort.""")

    def do_filter(self, params):
        global display_all
        if params.endswith(" -v"):
            params = params[:-3]
            display_all = True
        else:
            display_all = False
        if not params:
            print('No condition to filter given.  Type "help filter" for further information.')
        else:
            self.con.show_filtered_table(self.selected_table, params)

    def complete_filter(self, text, line, begidx, endidx):
        return self.column_name_completion(self.selected_table, text)

    def help_filter(self):
        print("""> filter <condition> [-v]
    Filter experiments of the currently selected class according to the given
    condition.  Any valid SQLite WHERE-statement can be used as a filter.
    Examples:
     - filter intensity <= 0.2
     - filter slope = 0.5
     - filter diffusion = "True"
     - filter datetime >= "2019-03-21"
     - filter perturbation != "gauss" AND duration > 20
    It is necessary to put the value for string, datetime or boolean argument
    in quotation marks as shown.
    To sort the filtered experiments, use SQLite syntax.
    Examples:
     - filter intensity > 0.1 ORDER BY intensity
     - filter perturbation != "gauss" ORDER BY size_total DESC
    The see all information about the filtered experiments, add "-v" at the end
    of the command.

    See also: sort, show.""")

    def do_sort(self, params):
        global display_all
        if params.endswith(" -v"):
            params = params[:-3]
            display_all = True
        else:
            display_all = False
        if not params:
            print('No parameter to sort given.  Type "help sort" for further information.')
        else:
            self.con.show_sorted_table(self.selected_table, params)

    def complete_sort(self, text, line, begidx, endidx):
        return self.column_name_completion(self.selected_table, text)

    def help_sort(self):
        print("""> sort <parameter> [desc] [-v]
    Sort the experiments by the value of a given parameter.
    To invert the order, add "desc" (descending) to the command.
    Examples:
     - sort intensity: show experiments with lowest intensity on top of the table.
     - sort size_total desc: show experiments with biggest file size on top.
    The see all information about the sorted experiments, add "-v" at the end
    of the command.

    See also: filter, show.""")

    def do_compare(self, params):
        # Check the run condition
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to make a comparison.')
            return
        # Parse the arguments
        global display_all
        display_all = False
        highlight = False
        other_params = []
        for p in params.split():
            if p == "-v":
                display_all = True
            elif p == "-h":
                highlight = True
            else:
                other_params.append(p)
        if other_params:
            ids = self.parse_multiple_ids(other_params)
            if ids is None:
                return
            if len(ids) < 2:
                print("Please specify at least 2 different IDs.")
                return
        elif self.con.get_table_length(self.selected_table) < 2:
            print("Selected experiment contains less than 2 entries.")
            return
        else:
            # No IDs means no filtering, i.e., all entries are compared
            ids = []
        self.con.show_comparison(self.selected_table, ids, highlight)

    def help_compare(self):
        print("""> compare [IDs] [-h] [-v]
    Show the difference between two or more entries of the selected experiment.
    This prints a table which includes only the parameters in which the
    specified entries differ.  Alternatively, add "-h" to show all parameters
    and highlight the differences in colour.  Disabled columns are not shown by
    default.  Add the argument "-v" to show disabled columns, the full date-time
    and the full comment.  Instead of an ID, "last" can be used to compare with
    the latest entry.  If no ID is specified, all entries of the selected
    experiment class are compared.""")

    ### Functionality to OPEN experiment files
    def do_open_mp4(self, params):
        """Open the mp4-file for an experiment specified by its name and ID."""
        if params.endswith(" -v"):
            verbose = True
            params = params[:-3].rstrip()
        else:
            verbose = False
        expname_id = self.parse_params_to_experiment(params)
        if not expname_id:
            return
        expname = "{}_{:03}".format(*expname_id)
        dir_ = os.path.join(self.exp_dir, expname)
        try:
            files = os.listdir(dir_)
        except FileNotFoundError:
            print("Folder does not exist:", dir_)
            return
        for f in files:
            if f.endswith(".mp4") and f.startswith(expname):
                break
        else:
            print("No mp4-file found in folder:", dir_)
            return
        path = os.path.join(self.exp_dir, expname, f)
        self.open_file(MP4_PLAYER, path, verbose)

    def complete_open_mp4(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def help_open_mp4(self):
        print(f"""> open_mp4 [experiment] <ID> [-v]
    Open the mp4-file of an experiment with {MP4_PLAYER}.""")
        self.print_param_parser_help()
        print("""
    Add "-v" to the command to see the output of the external programme.

    The programme to open mp4-files with can be configured in the Python script
    of the Experiment-Manager with the constant "MP4_PLAYER".""")

    def do_open_his(self, params):
        """Open the his-file for an experiment specified by its name and ID."""
        if params.endswith(" -v"):
            verbose = True
            params = params[:-3].rstrip()
        else:
            verbose = False
        expname_id = self.parse_params_to_experiment(params)
        if not expname_id:
            return
        expname = "{}_{:03}".format(*expname_id)
        path = os.path.join(self.exp_dir, expname, expname + "_his.nc")
        if not os.path.isfile(path):
            print("File does not exist:", path)
            return
        self.open_file(NETCDF_VIEWER, path, verbose)

    def complete_open_his(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def help_open_his(self):
        print(f"""> open_his [experiment] <ID> [-v]
    Open the NetCDF history-file of an experiment with {NETCDF_VIEWER}.""")
        self.print_param_parser_help()
        print("""
    Add "-v" to the command to see the output of the external programme.

    The programme to open NetCDF-files with can be configured in the Python
    script of the Experiment-Manager with the constant "NETCDF_VIEWER".""")

    def do_open_diag(self, params):
        """Open the diag-file for an experiment specified by its name and ID."""
        if params.endswith(" -v"):
            verbose = True
            params = params[:-3].rstrip()
        else:
            verbose = False
        expname_id = self.parse_params_to_experiment(params)
        if not expname_id:
            return
        expname = "{}_{:03}".format(*expname_id)
        path = os.path.join(self.exp_dir, expname, expname + "_diag.nc")
        if not os.path.isfile(path):
            print("File does not exist:", path)
            return
        self.open_file(NETCDF_VIEWER, path, verbose)

    def complete_open_diag(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def help_open_diag(self):
        print(f"""> open_diag [experiment] <ID> [-v]
    Open the NetCDF diagnostics-file of an experiment with {NETCDF_VIEWER}.""")
        self.print_param_parser_help()
        print("""
    Add "-v" to the command to see the output of the external programme.

    The programme to open NetCDF-files with can be configured in the Python
    script of the Experiment-Manager with the constant "NETCDF_VIEWER".""")

    ### Functionality for Data Analysis
    def do_calculate(self, params):
        if " " not in params:
            print("At least two arguments are needed.")
            return
        # Parse and check toolname
        i = params.find(" ")
        toolname = params[:i]
        if not extra_tools:
            print('There are no tools for calculations loaded.  '
                  'Add some tools in the code and restart the programme.')
            return
        if toolname not in extra_tools:
            print(f'Unknown tool: "{toolname}".  The available tools are:')
            for t in extra_tools:
                print(" -", t)
            return
        # Parse ID and optionally table name
        experiment_name_id = self.parse_params_to_experiment(params[i+1:])
        if experiment_name_id is None:
            return
        table_name, id_ = experiment_name_id
        expname = "{}_{:03}".format(table_name, id_)
        path = os.path.join(self.exp_dir, expname, expname + '_his.nc')
        try:
            val = extra_tools[toolname](path)
        except Exception as e:
            print(f'Tool "{toolname}" did not succeed on experiment {id_}.  '
                  'The error message is:')
            print(e)
            return
        if val is not None:
            print(f'-> {toolname}({table_name}:{id_}) = {val}')
        else:
            print(f'-> {toolname}({table_name}:{id_}) succeeded without return value.')

    def complete_calculate(self, text, line, begidx, endidx):
        return [p for p in extra_tools if p.startswith(text)]

    def help_calculate(self):
        print("""> calculate <tool> [experiment] <ID>
    Call an EMShell-Extension on an experiment and print the result.""")
        self.print_param_parser_help()
        print("""
    To make a tool available, load it into the Python script of the
    Experiment-Manager and add it to the dictionary "extra_tools".
    Its name must not contain any whitespace.""")

    def do_plot(self, params):
        # Check the run condition
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to plot its data.')
            return
        # Parse the arguments
        parameters = self.parse_plot_parameters(params, 2)
        if not parameters:
            return
        variables = parameters["variables"]
        # Get the data belonging to the parameters
        datas = self.get_multiple_data(self.selected_table, variables, parameters["condition"])
        if not datas:
            return
        # Print and plot the data
        for i, c in enumerate(["x", "y"]):
            print(f'-> {c} ({variables[i]}):', datas[i])
        plot_label = variables[1]
        if parameters["condition"]:
            plot_label += ' (' + parameters["condition"] + ')'
        plt.title(self.selected_table)
        plt.xlabel(variables[0])
        try:
            plt.plot(
                datas[0], datas[1], parameters["format_string"], label=plot_label,
            )
        except Exception as e:
            print("Plot did not succeed.  Error message:")
            print(e)
        else:
            plt.legend()
            if parameters["xmin"] is not None or parameters["xmax"] is not None:
                plt.xlim(parameters["xmin"], parameters["xmax"])
            if parameters["ymin"] is not None or parameters["ymax"] is not None:
                plt.ylim(parameters["ymin"], parameters["ymax"])
            if parameters["grid"]:
                plt.grid(True)
            if parameters["save_as"]:
                plt.savefig(get_unique_save_filename("plot_{}." + parameters["save_as"]))
            plt.show()

    def complete_plot(self, text, line, begidx, endidx):
        return self.plot_attribute_completion(text, [
            'f=', 'grid',
            'png', 'pdf', 'svg',
            'xmin=', 'xmax=', 'ymin=', 'ymax=',
        ])

    def help_plot(self):
        print("""> plot <variable> <variable> [condition] [-grid] [-f=<format>] [-{x|y}{min|max}=<number>] [-{png|pdf|svg}]
    Make a diagram showing the relation between the two variables for the
    selected experiment.""")
        self.print_general_plot_help("2d")
        print("""
    Example:
     - plot size_his duration id <= 10 -f=g.- -grid
       This command creates a plot in green with dots and lines on a grid
       comparing the integration time of the first ten experiments in the
       selected class of experiments with the size of their history-file.""")

    def do_scatter(self, params):
        # Check the run condition
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to plot its data.')
            return
        # Parse the arguments
        parameters = self.parse_plot_parameters(params, 3)
        if not parameters:
            return
        variables = parameters["variables"]
        # Get the data belonging to the parameters
        datas = self.get_multiple_data(self.selected_table, variables, parameters["condition"])
        if not datas:
            return
        # Print and plot the data
        for i, c in enumerate(["x", "y", "z"]):
            print('-> {} ({}):'.format(c, variables[i]), datas[i])
        plot_title = self.selected_table + ": " + variables[2]
        if parameters["condition"]:
            plot_title += ' (' + parameters["condition"] + ')'
        plt.title(plot_title)
        plt.xlabel(variables[0])
        plt.ylabel(variables[1])
        try:
            plt.scatter(
                datas[0], datas[1], c=datas[2],
                cmap=parameters["cmap"],
                vmin=parameters["zmin"], vmax=parameters["zmax"],
            )
        except Exception as e:
            print("Scatter did not succeed.  Error message:")
            print(e)
        else:
            plt.colorbar()
            if parameters["xmin"] is not None or parameters["xmax"] is not None:
                plt.xlim(parameters["xmin"], parameters["xmax"])
            if parameters["ymin"] is not None or parameters["ymax"] is not None:
                plt.ylim(parameters["ymin"], parameters["ymax"])
            if parameters["grid"]:
                plt.grid(True)
            if parameters["save_as"]:
                plt.savefig(get_unique_save_filename("scatter_{}." + parameters["save_as"]))
            plt.show()

    def complete_scatter(self, text, line, begidx, endidx):
        return self.plot_attribute_completion(text, [
            'cmap=', 'grid',
            'png', 'pdf', 'svg',
            'xmin=', 'xmax=', 'ymin=', 'ymax=', 'zmin=', 'zmax=',
        ])

    def help_scatter(self):
        print("""> scatter <variable> <variable> <variable> [condition] [-grid] [-cmap=<name>] [-{x|y|z}{min|max}=<number>] [-{png|pdf|svg}]
    Make a scatter plot showing the values of the third variable in colour in
    relation to the other two variables for the selected experiment.""")
        self.print_general_plot_help("3d")
        print("""
    Example:
     - scatter intensity sigma duration diffusion="False" OR Kdiff=0 -zmin=0
       This command creates a scatter plot with the default colourmap of
       matplotlib on a colour-axis starting from 0 which shows the integration
       time in relation to the intensity and the value of sigma for experiments
       of the current class without diffusion.""")

    def do_pcolor(self, params):
        # Check the run condition
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to plot its data.')
            return
        # Parse the arguments
        parameters = self.parse_plot_parameters(params, 3)
        if not parameters:
            return
        variables = parameters["variables"]
        # Get the data belonging to the parameters
        datas = self.get_multiple_data(self.selected_table, variables, parameters["condition"])
        if not datas:
            return
        # Arrange the data in a grid
        xvalues = sorted(set(datas[0]))
        yvalues = sorted(set(datas[1]))
        print(f'There are {len(xvalues)} unique x-values and {len(yvalues)} unique y-values.')
        data_grid = np.full((len(yvalues), len(xvalues)), np.nan)
        for i, zval in enumerate(datas[2]):
            data_grid[yvalues.index(datas[1][i]), xvalues.index(datas[0][i])] = zval
        # If no shading is applied, extend the axes to have each rectangle
        # centred at the corresponding (x,y)-value
        if parameters["shading"] == 'flat':
            xdiff2 = np.diff(xvalues)/2
            xaxis = [xvalues[0] - xdiff2[0], *(xvalues[:-1] + xdiff2), xvalues[-1] + xdiff2[-1]]
            ydiff2 = np.diff(yvalues)/2
            yaxis = [yvalues[0] - ydiff2[0], *(yvalues[:-1] + ydiff2), yvalues[-1] + ydiff2[-1]]
        else:
            xaxis = xvalues
            yaxis = yvalues
        # Print and plot the data
        print(f"-> x ({variables[0]}):", xvalues)
        print(f"-> y ({variables[1]}):", yvalues)
        print(f"-> z ({variables[2]}):")
        print(data_grid)
        plot_title = self.selected_table + ": " + variables[2]
        if parameters["condition"]:
            plot_title += ' (' + parameters["condition"] + ')'
        plt.title(plot_title)
        plt.xlabel(variables[0])
        plt.ylabel(variables[1])
        try:
            plt.pcolormesh(
                xaxis, yaxis, data_grid,
                cmap=parameters["cmap"], shading=parameters["shading"],
                vmin=parameters["zmin"], vmax=parameters["zmax"],
            )
        except Exception as e:
            print("Pcolormesh did not succeed.  Error message:")
            print(e)
        else:
            plt.colorbar()
            if parameters["xmin"] is not None or parameters["xmax"] is not None:
                plt.xlim(parameters["xmin"], parameters["xmax"])
            if parameters["ymin"] is not None or parameters["ymax"] is not None:
                plt.ylim(parameters["ymin"], parameters["ymax"])
            if parameters["grid"]:
                plt.grid(True)
            if parameters["save_as"]:
                plt.savefig(get_unique_save_filename("pcolor_{}." + parameters["save_as"]))
            plt.show()

    def complete_pcolor(self, text, line, begidx, endidx):
        return self.plot_attribute_completion(text, [
            'cmap=', 'shading', 'grid',
            'png', 'pdf', 'svg',
            'xmin=', 'xmax=', 'ymin=', 'ymax=', 'zmin=', 'zmax=',
        ])

    def help_pcolor(self):
        print("""> pcolor <variable> <variable> <variable> [condition] [-grid] [-shading] [-cmap=<name>] [-{x|y|z}{min|max}=<number>] [-{png|pdf|svg}]
    Make a pseudo-colour plot showing the values of the third variable in colour
    in relation to the two other variables for the selected experiment.""")
        self.print_general_plot_help("3d", shading=True)
        print("""
    Example:
     - pcolor intensity sigma duration diffusion="False" OR Kdiff=0 -zmin=0
       This command creates a pseudo-colour plot with the default colourmap of
       matplotlib on a colour-axis starting from 0 which shows the integration
       time in relation to the intensity and the value of sigma for experiments
       of the current class without diffusion.""")

    def do_save_figure(self, params):
        if not params:
            params = "png"
        if params in ["png", "pdf", "svg"]:
            plt.savefig(get_unique_save_filename("figure_{}." + params))
        else:
            plt.savefig(params)

    def complete_save_figure(self, text, line, begidx, endidx):
        if "." not in text:
            return [e for e in ("png", "pdf", "svg") if e.startswith(text)]
        stub = text[:text.rfind(".")]
        return [n for n in (stub + ".png", stub + ".pdf", stub + ".svg") if n.startswith(text)]

    def help_save_figure(self):
        print("""> save_figure [filename or -type]
    Save the currently opened figure in a file.  One can either specify the full
    filename or one of the filetypes png, pdf or svg alone.  If only the type is
    specified, then the name is automatically set to "figure_#.type", where
    "#" is an integer such that the filename is unique and "type" is the given
    filetype.  If no type and no name is given, the figure is saved as png.""")

    def do_new_figure(self, params):
        """Open a window to draw the next plot in a new figure."""
        plt.figure()

    ### Functionality to MODIFY the entries
    def do_new_comment(self, params):
        # Check and parse parameters
        experiment_name_id = self.parse_params_to_experiment(params)
        if experiment_name_id is None:
            return
        table_name, id_ = experiment_name_id
        # Print the current entry fully
        global display_all
        display_all = True
        self.con.show_filtered_table(table_name, f"id = {id_}")
        # Ask for user input
        print("Write a new comment for this entry (Ctrl+D to finish, Ctrl+C to cancel):")
        comment_lines = []
        while True:
            try:
                line = input()
            except EOFError:
                break
            except KeyboardInterrupt:
                print("Cancelling.")
                return
            comment_lines.append(line)
        comment = "\n".join(comment_lines).strip()
        # Update the comment
        if self.con.set_comment(table_name, id_, comment):
            self.con.save_database()
            print("New comment was saved.")
        else:
            print("An error occured.")

    def help_new_comment(self):
        print("""> new_comment [experiment] <ID>
    Ask the user to enter a new comment for an experiment.""")
        self.print_param_parser_help()

    ### Functionality to CLEAN up
    def do_remove(self, params):
        # Check conditions and parameters
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to remove entries from it.')
            return
        if not params:
            print('No IDs given.  Specify at least one ID to remove its entry.')
            return
        # Parse parameters
        ids = self.parse_multiple_ids(params.split())
        if ids is None:
            print("No data removed.")
            return

        # Print full information of selected entries
        global display_all
        display_all = True
        if len(ids) == 1:
            statement = "id = {}".format(*ids)
        else:
            statement = "id IN {}".format(tuple(ids))
        print('WARNING: the following entries will be DELETED:')
        self.con.show_filtered_table(self.selected_table, statement)

        # Print full information of related folders
        folders = []
        print('WARNING: the following folders and files will be DELETED:')
        for id_ in ids:
            expname = "{}_{:03}".format(self.selected_table, id_)
            folder = os.path.join(self.exp_dir, expname)
            try:
                files = os.listdir(folder)
            except FileNotFoundError:
                print(' x Folder does not exist:', folder)
            else:
                folders.append(folder)
                print(" -", folder)
                for f in sorted(files):
                    print("   -", f)

        # Final check
        print('Do you really want to permanently delete these files, folders and entries?')
        print('This cannot be undone.')
        answer = input('Continue [yes/no] ? ')
        if answer == 'yes':
            # Remove entries
            for id_ in ids:
                self.con.delete_entry(self.selected_table, id_)
                print('Deleted entry', id_, 'from experiment "{}".'.format(self.selected_table))
            self.con.save_database()
            # Remove files
            for folder in folders:
                try:
                    shutil.rmtree(folder)
                except OSError as e:
                    print('Error deleting folder {}:'.format(folder), e)
                else:
                    print('Deleted folder {}.'.format(folder))
        elif answer == "no":
            # Do nothing.
            pass
        else:
            print('Answer was not "yes".  No data removed.')

    def help_remove(self):
        print("""> remove <IDs>
    Delete the entry with the given ID from the currently selected class of
    experiments and all files associated with it.  Multiple IDs can be specified
    to remove several entries and their folders at once.  Instead of an ID, the
    argument "last" can be used to choose the latest entry.
    Before the data is deleted, the user is asked to confirm, which must be
    answered with "yes".
    The remove an empty class of experiments, use "remove_selected_class".""")

    def do_remove_selected_class(self, params):
        global display_all
        display_all = True

        # Check conditions and parameters
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to remove it.')
            return
        if params:
            print('This command takes no attributes.  Cancelling.')
            return

        if self.con.get_table_length(self.selected_table):
            print('Selected experiment class contains experiments.  Remove these '
                  'experiments first before removing the class.  Cancelling.')
            return

        print('Do you really want to permanently delete the experiment class '
              f'"{self.selected_table}"?')
        print('This cannot be undone.')
        print('The EMShell will exit after the deletion.')
        answer = input('Continue [yes/no] ? ')
        if answer == 'yes':
            self.con.delete_table(self.selected_table)
            self.con.save_database()
            print(f'Deleted experiment class "{self.selected_table}".')
            print(f'EMShell must be restarted to update the list of tables.')
            return True
        else:
            print('Answer was not "yes".  No data removed.')

    def help_remove_selected_class(self):
        print("""> remove_selected_class
    Delete the currently selected class of experiments from the database.
    This works only if no entries are associated with this experiment class.
    Before the class is removed, the user is asked to confirm, which must be
    answered with "yes".  To remove entries from the selected class, use the
    command "remove".""")

    def do_print_disk_info(self, params):
        # Explanation: https://stackoverflow.com/a/12327880/3661532
        statvfs = os.statvfs(self.exp_dir)
        available_space = statvfs.f_frsize * statvfs.f_bavail
        print(
            "Available disk space in the experiment directory:",
            SIZE_FORMAT.format(available_space / 1024**3),
            "GiB"
        )
        print("Experiment directory:", self.exp_dir)

    def help_print_disk_info(self):
        print("""> print_disk_info
    Display the available disk space in the experiment directory and its path.
    The available disk space is printed in Gibibytes (GiB), this means:
        1 GiB = 1024^3 bytes.
    This is used instead of the Gigabyte (GB), which is defined as:
        1 GB = 1000^3 bytes,
    because the GiB underestimates the free disk space, which is considered
    safer than overestimating it.  In contrast, the size used by experiments is
    displayed in Megabytes (MB), where
        1 MB = 1000^2 bytes,
    since the used space should rather be overestimated.
    More information about the unit: https://en.wikipedia.org/wiki/Gibibyte ."""
        )

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
        if self.selected_table:
            print("select")
            self.prompt = "(EMS) "
            self.selected_table = ""
        else:
            print("exit")
            return True

    ### Behaviour for empty input
    def emptyline(self):
        pass

    ### Helper functions
    def table_name_completion(self, text):
        return [table.name for table in self.con.tables if table.name.startswith(text)]

    def plot_attribute_completion(self, text, parameters):
        if not self.selected_table:
            print("\nError: select an experiment first!")
            return []
        parameters += self.con.get_column_names(self.selected_table)
        parameters.extend(extra_tools.keys())
        return [p for p in parameters if p.startswith(text)]

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

    def parse_params_to_experiment(self, params: str) -> [str, int]:
        """Parse and check input of the form "[experiment] <ID>".

        The argument "experiment" is optional, the ID is necessary.
        Instead of an ID, "last" can be used to refer to the latest entry.
        If no experiment is given, the selected experiment is taken.
        Return None and print a message if input is not valid, otherwise
        return the experiment name and the ID."""
        if not params:
            print("No ID given.")
            return None
        params = params.split(" ")
        # Parse the ID
        specifier = params.pop()
        if specifier == "last":
            id_ = "last"
        else:
            try:
                id_ = int(specifier)
            except ValueError:
                print(f'Last argument is not a valid ID: "{specifier}".')
                return None
        # Parse and check the table name
        table_name = " ".join(params).strip()
        if table_name:
            if not self.con.table_exists(table_name):
                print(f'No experiment class with the name "{table_name}" exists.')
                return None
        else:
            if not self.selected_table:
                print("No experiment selected and no experiment name given.")
                return None
            table_name = self.selected_table
        # Check the ID
        if id_ == "last":
            id_ = self.con.get_latest_entry(table_name)
            if id_ is None:
                print(f'No entry exists for the experiment class "{table_name}".')
                return None
        elif not self.con.entry_exists(table_name, id_):
            print(f'No entry with ID {id_} exists for the experiment class "{table_name}".')
            return None
        return table_name, id_

    @staticmethod
    def parse_plot_parameters(params, n_necessary):
        parameters = {
            "variables": [],
            "condition": "",
            "save_as": None,
            "grid": False,
            "xmin": None,
            "xmax": None,
            "ymin": None,
            "ymax": None,
            "zmin": None,
            "zmax": None,
            "cmap": None,
            "shading": "flat",
            "format_string": "",
        }
        param_list = params.split(" ")
        # Necessary parameters are at the beginning
        if len(param_list) < n_necessary:
            print(f'Not enough parameters given, {n_necessary} are necessary.')
            return None
        for i in range(n_necessary):
            parameters["variables"].append(param_list.pop(0))
        # The last parameters can be used to modify the plot
        for param in reversed(param_list):
            if not param:
                # Ignore empty strings
                pass
            elif param in ["-png", "-pdf", "-svg"]:
                parameters["save_as"] = param[1:]
            elif param == "-grid":
                parameters["grid"] = True
            elif param == "-shading":
                parameters["shading"] = "gouraud"
            elif param.startswith("-f="):
                parameters["format_string"] = param[3:]
            elif param.startswith("-cmap="):
                parameters["cmap"] = param[6:]
            elif xyzminmax.match(param):
                param_name = param[1:5]
                param_value = param[6:]
                try:
                    parameters[param_name] = float(param_value)
                except ValueError:
                    print(f'Value for {param_name} cannot be converted to '
                          f'float: "{param_value}".')
            else:
                # End of extra parameters, beginning of the SQL statement
                break
            param_list.pop()
        # Every other parameter is treated as an SQL condition to filter the data
        parameters["condition"] = " ".join(param_list).strip()
        return parameters

    def parse_multiple_ids(self, param_list):
        """Transform the given list of strings into a list of IDs.

        This method also checks that the IDs are valid entries of the
        selected database and returns None if an invalid ID is given.
        Otherwise, a sorted list of unique values is returned.
        It parses "last" as the latest entry."""
        ids = []
        for id_ in param_list:
            if id_ == "last":
                id_ = self.con.get_latest_entry(self.selected_table)
                if id_ is None:
                    print(f'No entry exists for the experiment class "{self.selected_table}".')
                    return
                elif id_ in ids:
                    continue
            else:
                try:
                    id_ = int(id_)
                except ValueError:
                    print('Parameter', id_, 'is not a valid ID.')
                    return
                if id_ in ids:
                    continue
                if not self.con.entry_exists(self.selected_table, id_):
                    print('No entry with ID', id_, 'exists for the selected experiment.')
                    return
            ids.append(id_)
        return sorted(ids)

    def open_file(self, command, path, verbose=False):
        if verbose:
            print("Opening file {} with {} in verbose-mode.".format(path, command))
            subprocess.Popen(
                [command, path],
                # Disable standard input via the shell, for example with mplayer.
                stdin=subprocess.DEVNULL,
            )
        else:
            print("Opening file {} with {} silently.".format(path, command))
            subprocess.Popen(
                [command, path],
                # Disable standard input via the shell, for example with mplayer.
                stdin=subprocess.DEVNULL,
                # Throw away output and error messages.
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    def get_multiple_data(self, table_name, variables, condition):
        datas = []
        for var in variables:
            data = self.get_data(table_name, var, condition)
            if not data:
                print(f'Get data for parameter "{var}" ', end="")
                if condition:
                    print(f'under condition "{condition}" ', end="")
                print('failed.')
                return None
            datas.append(data)
        return datas

    def get_data(self, table_name, parameter, condition):
        if self.con.is_valid_column(parameter, table_name):
            return self.con.get_data(table_name, parameter, condition)
        elif parameter in extra_tools:
            data = []
            ids = self.con.get_data(table_name, "id", condition)
            print(f'Calculating data with tool "{parameter}" for {len(ids)} experiments.')
            tstart = time.time()
            for id_ in ids:
                expname = "{}_{:03}".format(table_name, id_)
                path = os.path.join(self.exp_dir, expname, expname + '_his.nc')
                try:
                    val = extra_tools[parameter](path)
                except Exception as e:
                    print(f'Tool "{parameter}" did not succeed on experiment {id_}.  '
                          'The error message is:')
                    print(e)
                    print('Using value 0 (zero) as fallback.')
                    val = 0
                data.append(val)
            print('Calculation finished in {:.3f} seconds.'.format(time.time() - tstart))
            return data
        else:
            print(f'Parameter "{parameter}" is neither a parameter of the '
                  'selected experiment nor a loaded extension.')
            return []

    @staticmethod
    def print_general_plot_help(plot_type, shading=False):
        print("""
    The variables can be read from the database or calculated using shell-
    extensions.  Type "help calculate" for further information on extensions.

    A condition can be used to filter or sort the data.  Type "help filter" for
    further information on filtering and sorting.

    To draw a grid on the plot, add "-grid" to the command."""
        )
        if shading:
            print("""
    To make a smooth plot instead of rectangles of solid colour, add "-shading"
    to the command.  This enables the Gouraud-shading of matplotlib."""
            )
        if plot_type == "2d":
            print("""
    To specify the type of plot, use "-f=<format>", where <format> is a format
    string for matplotlib like "-f=o" to plot with dots instead of a line, or
    "-f=rx" to plot with red crosses.  For further information on format strings
    in matplotlib, see the Notes section on the following website:
    https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.plot.html

    To specify the range in x- or y-direction, use the following attributes:
        -xmin=<v>, -xmax=<v>, -ymin=<v>, -ymax=<v>
    with a number as the argument <v>."""
            )
        elif plot_type == "3d":
            print("""
    To specify the colour map of the plot, use "-cmap=<name>", where <name> is
    the name of a colour map for matplotlib like "-cmap=bwr" to plot in blue-
    white-red or "-cmap=jet" for a colourful rainbow.  For further examples,
    consider the website https://matplotlib.org/users/colormaps.html .

    To specify the range in x-, y- or z-direction, use the following attributes:
        -xmin=<v>, -xmax=<v>, -ymin=<v>, -ymax=<v>, -zmin=<v>, -zmax=<v>
    with a number as the argument <v>.  The z-axis refers to the colourbar."""
            )
        print("""
    Add either "-png" or "-pdf" or "-svg" to the command to save the figure in
    a file of the corresponding filetype.  Alternatively, use the command 
    "save_figure" after creating the plot.

    The order of the qualifiers starting with a dash ("-") does not play a role,
    but they must be at the end of the command, i.e., after the variables and
    after the filter condition (if any).  If a qualifier is given more than
    once, only the first occurence is taken into account.

    While the window with the plot is open, the shell can be used as usual.
    If another plot command is executed with the figure still open, the new plot
    is superimposed onto the previous one.  To draw the next plot command in a
    new figure instead, while keeping the previous window open, use the command
    "new_figure"."""
        )

    @staticmethod
    def print_param_parser_help():
        print("""
    If only an ID and no experiment name is given, take the currently selected
    class of experiments.  Instead of an ID, the value "last" can be used to
    choose the latest entry."""
        )


def get_unique_save_filename(template):
    id_ = 0
    filename = template.format(id_)
    while os.path.exists(filename):
        id_ += 1
        filename = template.format(id_)
    print(f'Saving as "{filename}".')
    return filename


def string_format_table(table_name, columns, rows, highlight_columns=None):
    # Get current date and time to make datetime easy to read
    if not ISO_DATETIME:
        dt_now = datetime.datetime.now()
    # Convert rows to list of lists (because fetchall() returns a list of tuples
    # and to work on a copy of it)
    rows = [list(row) for row in rows]
    # Remove ignored columns
    columns = columns.copy()
    indices_to_remove = []
    if not display_all:
        for i, (n, t) in enumerate(columns):
            if n in hidden_information:
                indices_to_remove.append(i)
        for i in sorted(indices_to_remove, reverse=True):
            columns.pop(i)
            for row in rows:
                row.pop(i)
    # Get width necessary for each column, get total size of the files
    # associated with the table and process the entries of the table.
    lengths = [len(n) for n, t in columns]
    table_size = 0
    for row in rows:
        for i, val in enumerate(row):
            # Check name of column
            if columns[i][0] == "comment":
                if display_all or COMMENT_MAX_LENGTH == "FULL":
                    # No formatting needed
                    pass
                else:
                    # Replace linebreaks
                    val = val.replace(*LINEBREAK_REPLACE)
                    if type(COMMENT_MAX_LENGTH) is int and COMMENT_MAX_LENGTH > 0:
                        # Cut comments which are too long and end them with an ellipsis
                        if len(val) > COMMENT_MAX_LENGTH:
                            val = val[:COMMENT_MAX_LENGTH-1] + "…"
                    elif COMMENT_MAX_LENGTH == "AUTO":
                        # Cut comments later
                        pass
                    else:
                        raise ValueError(
                            'COMMENT_MAX_LENGTH has to be "AUTO" or "FULL" '
                            'or a positive integer, not "{}".'.format(COMMENT_MAX_LENGTH)
                        )
                    row[i] = val
            elif columns[i][0] == "datetime":
                if display_all:
                    # No formatting needed
                    pass
                elif ISO_DATETIME:
                    # Cut unnecessary parts of the date and time
                    if "datetime_year" in hidden_information:
                        val = val[5:]
                    if "datetime_seconds" in hidden_information:
                        val = val.split(".")[0][:-3]
                    val = val.replace('T', ',')
                else:
                    # Create datetime-object from ISO-format
                    dt_obj = datetime.datetime.strptime(val, "%Y-%m-%dT%H:%M:%S.%f")
                    val = make_nice_time_string(dt_obj, dt_now)
                row[i] = val
            elif columns[i][0] == "size_total":
                if val > 0:
                    table_size += val
            # Check type of column and adopt
            if columns[i][1] == "TEXT" or columns[i][1] == "INTEGER":
                lengths[i] = max(lengths[i], len(str(val)))
            elif columns[i][0] in ["size_diag", "size_flux", "size_his",
                                   "size_mp4", "size_total"]:
                lengths[i] = max(lengths[i], len(SIZE_FORMAT.format(val)))
            elif columns[i][1] == "REAL":
                lengths[i] = max(lengths[i], len(FLOAT_FORMAT.format(val)))
            else:
                # This is an unexpected situation, which probably means that
                # sqlite3 does not work as it was, when this script was written.
                raise Exception(
                    "unknown type {} of column {} in table {}."
                    .format(*columns[i], table_name)
                )
    if (
            COMMENT_MAX_LENGTH == "AUTO"
            and not display_all
            and "comment" == columns[-1][0]
    ):
        line_length = shutil.get_terminal_size().columns
        total_length = sum(lengths[:-1]) + len(LIMITER) * len(lengths[:-1])
        comment_length = line_length - total_length % line_length
        lengths[-1] = comment_length
        for row in rows:
            comment = row[-1]
            # Cut comments which are too long and end them with an ellipsis
            if len(comment) > comment_length:
                comment = comment[:comment_length-1] + "…"
                row[-1] = comment
    # Create top line of the text to be returned
    text = "Experiment: " + table_name
    if "size_total" not in hidden_information or display_all:
        text += " (" + SIZE_FORMAT.format(table_size) + " MB)"
    if len(indices_to_remove) == 1:
        text += " (1 parameter hidden)"
    elif len(indices_to_remove) > 1:
        text += " ({} parameters hidden)".format(len(indices_to_remove))
    text += "\n"
    column_names_formatted = []
    format_strings = []
    for (n, t), l in zip(columns, lengths):
        # Column names centred, except comment, since it is the last column.
        cname_form = f"{n:^{l}}"
        # Numbers right justified,
        # comments left justified,
        # other texts centred.
        if n in ["size_diag", "size_flux", "size_his", "size_mp4", "size_total"]:
            f_string = SIZE_FORMAT.replace(":", ":>"+str(l))
        elif t == "REAL":
            f_string = FLOAT_FORMAT.replace(":", ":>"+str(l))
        elif t == "INTEGER":
            f_string = "{:>" + str(l) + "}"
        elif n == "comment":
            f_string = "{:" + str(l) + "}"
            cname_form = n
        else:
            f_string = "{:^" + str(l) + "}"
        if highlight_columns is not None and n in highlight_columns:
            f_string = COLOURS_HIGHLIGHT + f_string + COLOURS_END
            cname_form = COLOURS_HIGHLIGHT + cname_form + COLOURS_END
        format_strings.append(f_string)
        column_names_formatted.append(cname_form)
    text += LIMITER.join(column_names_formatted) + "\n"
    if COLOURS:
        row_colours = itertools.cycle(COLOURS)
    for row in rows:
        if COLOURS:
            col_colours = itertools.cycle(next(row_colours))
        text_cols = [next(col_colours) + f_str.format(val) if COLOURS
                     else f_str.format(val) for f_str, val in zip(format_strings, row)]
        text += LIMITER.join(text_cols)
        text += "\n"
    if text.endswith("\n"):
        text = text[:-1]
    if COLOURS:
        text = text + COLOURS_END
    return text


def make_nice_time_string(datetime_object, datetime_reference):
    """Create an easy to read representation of the datetime object."""
    if datetime_object > datetime_reference:
        return LANG.FUTURE
    elif datetime_object.date() == datetime_reference.date():
        # same day
        dt_diff = datetime_reference - datetime_object
        dt_diff_minutes = dt_diff.seconds / 60
        if dt_diff_minutes < 1:
            return LANG.JUST_NOW
        elif dt_diff_minutes < 60:
            return LANG.AGO_MINUTES.format(int(dt_diff_minutes))
        else:
            return LANG.AGO_HOURS.format(int(dt_diff_minutes / 60), int(dt_diff_minutes % 60))
    elif datetime_object.date() == (datetime_reference - datetime.timedelta(days=1)).date():
        # yesterday
        return datetime_object.strftime(LANG.YESTERDAY + ", %H:%M")
    elif datetime_object.date() > (datetime_reference - datetime.timedelta(days=7)).date():
        # this week, i.e., within the last six days
        return datetime_object.strftime("%a, %H:%M")
    else:
        format_string = "%b-%d, %H:%M"
        if "datetime_year" not in hidden_information:
            format_string = "%Y-" + format_string
        return datetime_object.strftime(format_string)


# Get the directory of the experiments
if len(sys.argv) == 1:
    try:
        from param import Param
    except ModuleNotFoundError:
        raise Exception(
            "Please activate fluid2d or specify the experiment-folder as "
            "argument when starting this programme."
        )
    param = Param(None)  # it is not necessary to specify a defaultfile for Param
    datadir = param.datadir
    del param
    if datadir.startswith("~"):
        datadir = os.path.expanduser(datadir)
elif len(sys.argv) == 2:
    datadir = sys.argv[1]
else:
    raise Exception("More than one argument given.")

# Use fancy colours during 6 days of Carnival
try:
    from dateutil.easter import easter
except ModuleNotFoundError:
    # No fun at carnival possible
    pass
else:
    date_today = datetime.date.today()
    carnival_start = easter(date_today.year)-datetime.timedelta(days=46+6)  # Weiberfastnacht
    carnival_end = easter(date_today.year)-datetime.timedelta(days=46)  # Aschermittwoch
    if carnival_start <= date_today < carnival_end or COLOURS == "HAPPY":
        print("It's carnival!  Let's hope your terminal supports colours!")
        # This looks like a rainbow on many terminals, e.g. xterm.
        COLOURS = (
            ("\033[41m",),
            ("\033[101m",),
            ("\033[43m",),
            ("\033[103m",),
            ("\033[102m",),
            ("\033[106m",),
            ("\033[104m",),
            ("\033[105m",),
        )

# Activate interactive plotting
plt.ion()

# Compile regular expression to match parameters
# This matches -xmin=, -xmax=, -ymin=, -ymax=, -zmin=, -zmax=
xyzminmax = re.compile("-[xyz]m(in|ax)=")

# Start the shell
ems_cli = EMShell(datadir)
ems_cli.cmdloop()
