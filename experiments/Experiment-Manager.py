# Markus Reinert, May 2019
#
# Command-line interface for the fluid2d Experiment Management System (EMS)
#
# This programme can be started directly with Python if fluid2d is activated.
# Otherwise, the path to the folder with the experiment-files (that is the
# value of param.datadir) must be specified as a command-line argument.

import os
import cmd
import time
import numpy as np
import shutil
import sqlite3 as dbsys
import datetime
import itertools
import subprocess
import matplotlib.pyplot as plt

import EMShell_Extensions as EMExt

try:
    from param import Param
except ModuleNotFoundError:
    raise Exception(
        "Please activate fluid2d to use this programme!"
    )


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
    "wavenumber": EMExt.get_strongest_wavenumber,
    "wavelength": lambda hisname: 1/EMExt.get_strongest_wavenumber(hisname),
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
# Set COLOURS to False or None to disable colours.
# Otherwise, COLOURS must be list, of which each element describes the colours used in one row.
# The colours of every row are described by a list of colour codes,
# cf. https://en.wikipedia.org/wiki/ANSI_escape_code#Colors.
# An arbitrary number of colours can be specified.
# Apart from background colours, also text colours and font styles can be specified.
# Example to colour rows alternately with white and default colour:
# COLOURS = (
#     ("\033[47m",),
#     ("\033[49m",),
# )
# Example to colour columns alternately with white and default colour:
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

    def save_database(self):
        self.connection.commit()

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

    def is_valid_column(self, column_name, table_name=""):
        for table in self.tables:
            if not table_name or table_name == table.name:
                for n, t in table.columns:
                    if column_name == n:
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

    def get_data(self, column_name, condition=""):
        try:
            self.c.execute(f'SELECT {column_name} FROM "{self.name}" WHERE {condition}' if condition else
                           f'SELECT {column_name} FROM "{self.name}"')
        except dbsys.OperationalError as e:
            print(f'SQL error for experiment "{self.name}":', e)
            return []
        else:
            return [e[0] for e in self.c.fetchall()]

    def entry_exists(self, id_):
        self.c.execute('SELECT id FROM "{}" WHERE id = ?'.format(self.name), (id_,))
        if self.c.fetchone() is not None:
            return True
        return False

    def delete_entry(self, id_):
        self.c.execute('DELETE FROM "{}" WHERE id = ?'.format(self.name), (id_,))

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
        + "\nPress Tab-key for auto-completion of commands, table names, etc."
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
        return []

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
     - filter datetime >= "2019-03-2
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

    ### Functionality to OPEN experiment files
    def do_open_mp4(self, params):
        """Open the mp4-file for an experiment specified by its name and ID."""
        expname = self.parse_params_to_experiment(params)
        if not expname:
            return
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

    ### Functionality for Data Analysis
    def do_calculate(self, params):
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to do calculations on it.')
            return
        try:
            toolname, id_ = params.split()
        except ValueError:
            print('Invalid syntax.  The correct syntax is: calculate <tool> <id>.  '
                  'Make sure the name of the tool contains no space.')
            return
        if not extra_tools:
            print('There are no tools for calculations loaded.  '
                  'Add some tools in the code and restart the programme.')
            return
        if toolname not in extra_tools:
            print(f'Unknown tool: "{toolname}".  The available tools are:')
            for t in extra_tools:
                print(" -", t)
            return
        try:
            id_ = int(id_)
        except ValueError:
            print(f'Second parameter is not a valid id: "{id_}".')
            return
        if not self.con.entry_exists(self.selected_table, id_):
            print(f'No entry with the id {id_} exists for the selected experiment.')
            return
        expname = "{}_{:03}".format(self.selected_table, id_)
        path = os.path.join(self.exp_dir, expname, expname + '_his.nc')
        try:
            val = extra_tools[toolname](path)
        except Exception as e:
            print(f'Tool "{toolname}" did not succeed on experiment {id_}.  '
                  'The error message is:')
            print(e)
            return
        if val:
            print(f'-> {toolname}({id_}) = {val}')
        else:
            print(f'-> {toolname}({id_}) succeeded without return value.')

    def complete_calculate(self, text, line, begidx, endidx):
        return [p for p in extra_tools if p.startswith(text)]

    def help_calculate(self):
        print("""> calculate <tool> <id>
    Call an EMShell-Extension on an experiment and print the result.
    This works always on the entries of the currently selected experiment.
    To make a tool available, load it into the Python script of the
    Experiment-Manager and add it to the dictionary "extra_tools".
    Its name must not contain any whitespace.""")

    def do_plot(self, params):
        # Check the run condition
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to plot its data.')
            return
        # Parse the arguments
        param_list = params.split(" ")
        # Two parameters are necessary
        if len(param_list) < 2:
            print('No two parameters given.  Type "help plot" for further information.')
            return
        param1 = param_list[0]
        param_list.pop(0)
        param2 = param_list[0]
        param_list.pop(0)
        # The last parameters can be used to modify the plot
        save_as = None
        draw_grid = False
        format_string = ""
        for param in reversed(param_list):
            if param == "-png":
                save_as = "PNG"
                param_list.pop()
            elif param == "-pdf":
                save_as = "PDF"
                param_list.pop()
            elif param == "-grid":
                draw_grid = True
                param_list.pop()
            elif param.startswith("-f="):
                format_string = param[3:]
                param_list.pop()
            else:
                # End of extra parameters, beginning of the SQL statement
                break
        # Every other parameter is an SQL condition to filter the data
        condition = " ".join(param_list).strip()
        # Get the data belonging to the parameters
        data1 = self.get_data(self.selected_table, param1, condition)
        if not data1:
            print(f'Get data for first parameter "{param1}" ', end="")
            if condition:
                print(f'under condition "{condition}" ', end="")
            print('failed.')
            return
        data2 = self.get_data(self.selected_table, param2, condition)
        if not data2:
            print(f'Get data for second parameter "{param2}" ', end="")
            if condition:
                print(f'under condition "{condition}" ', end="")
            print('failed.')
            return

        print("-> x:", data1)
        print("-> y:", data2)
        plt.xlabel(param1)
        plt.ylabel(param2)
        if draw_grid:
            plt.grid()
        try:
            plt.plot(
                data1, data2, format_string,
                label=f'{self.selected_table} ({condition if condition else "all data"})',
            )
        except Exception as e:
            print("Plot did not succeed.  Error message:")
            print(e)
        else:
            plt.legend()
            if save_as == "PDF":
                figid = 0
                filename = f"figure_{figid}.pdf"
                while os.path.exists(filename):
                    figid += 1
                    filename = f"figure_{figid}.pdf"
                print(f'Saving plot as "{filename}".')
                plt.savefig(filename)
            elif save_as == "PNG":
                figid = 0
                filename = f"figure_{figid}.png"
                while os.path.exists(filename):
                    figid += 1
                    filename = f"figure_{figid}.png"
                print(f'Saving plot as "{filename}".')
                plt.savefig(filename)
            plt.show()

    def complete_plot(self, text, line, begidx, endidx):
        if not self.selected_table:
            return []
        parameters = ["f=", "grid", "png", "pdf"]
        parameters += self.con.get_column_names(self.selected_table)
        parameters.extend(extra_tools.keys())
        return [p for p in parameters if p.startswith(text)]

    def help_plot(self):
        print("""> plot <parameter_1> <parameter_2> [condition] [-f=<format>] [-grid] [-{png|pdf}]
    Make a diagram showing the relation between the two parameters for the
    selected experiment.

    The parameters can be read from the database or calculated using shell-
    extensions.  Type "help calculate" for further information on extensions.

    A condition can be used to filter or sort the data.  Type "help filter" for
    further information on filtering and sorting.

    To specify the type of plot, use "-f=<format>", where <format> is a format
    string for matplotlib like "-f=o" to plot with dots instead of a line, or
    "-f=rx" to plot with red crosses.  For further information on format strings
    in matplotlib, see the Notes section on the following website:
    https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.plot.html
    If multiple format strings are specified, only the first one is considered.

    To draw a grid on the plot, add "-grid" to the command.

    Add either "-png" or "-pdf" to the command to save the figure in a file of
    the corresponding filetype.  If both are specified, only the first argument
    is taken into account.

    The order of the qualifiers starting with a dash ("-") does not play a role,
    but they must be at the end of the command.

    While the window with the plot is open, the shell can be used as usual.
    If another plot command is executed with the figure still open, the new plot
    is superimposed onto the previous one.  Matplotlib automatically selects a
    new colour if no colour is specified explicitely.  If the -grid argument is
    specified a second time in another plot command, the grid is switched off.

    Example:
     - plot size_his duration id <= 10 -f=g.- -grid
       This command creates a plot in green with dots and lines on a grid
       comparing the integration time of the first ten experiments in the
       selected class of experiments with the size of their history-file.""")

    def do_scatter(self, params):
        # Check the run condition
        if not self.selected_table:
            print('No experiment selected.  '
                  'Select an experiment to make a scatter plot of its data.')
            return
        # Parse the arguments
        param_list = params.split(" ")
        # Two parameters are necessary
        if len(param_list) < 3:
            print('No three parameters given.  Type "help scatter" for further information.')
            return
        parameters = []
        for i in range(3):
            parameters.append(param_list.pop(0))
        # The last parameters can be used to modify the plot
        save_as = None
        draw_grid = False
        cmap = None
        xmin = None
        xmax = None
        ymin = None
        ymax = None
        zmin = None
        zmax = None
        for param in reversed(param_list):
            if param == "-png":
                save_as = "PNG"
                param_list.pop()
            elif param == "-pdf":
                save_as = "PDF"
                param_list.pop()
            elif param == "-grid":
                draw_grid = True
                param_list.pop()
            elif param.startswith("-cmap="):
                cmap = param[6:]
                param_list.pop()
            elif param.startswith("-xmin="):
                try:
                    xmin = float(param[6:])
                except ValueError:
                    print(f'Value for xmin cannot be converted to float: "{xmin}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-xmax="):
                try:
                    xmax = float(param[6:])
                except ValueError:
                    print(f'Value for xmax cannot be converted to float: "{xmax}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-ymin="):
                try:
                    ymin = float(param[6:])
                except ValueError:
                    print(f'Value for ymin cannot be converted to float: "{ymin}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-ymax="):
                try:
                    ymax = float(param[6:])
                except ValueError:
                    print(f'Value for ymax cannot be converted to float: "{ymax}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-zmin="):
                try:
                    zmin = float(param[6:])
                except ValueError:
                    print(f'Value for zmin cannot be converted to float: "{zmin}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-zmax="):
                try:
                    zmax = float(param[6:])
                except ValueError:
                    print(f'Value for zmax cannot be converted to float: "{zmax}" '
                          '-- ignoring it.')
                param_list.pop()
            else:
                # End of extra parameters, beginning of the SQL statement
                break
        # Every other parameter is an SQL condition to filter the data
        condition = " ".join(param_list).strip()
        # Get the data belonging to the parameters
        datas = []
        for parameter in parameters:
            data = self.get_data(self.selected_table, parameter, condition)
            if not data:
                print(f'Get data for parameter "{parameter}" ', end="")
                if condition:
                    print(f'under condition "{condition}" ', end="")
                print('failed.')
                return
            datas.append(data)
        print("-> x:", datas[0])
        print("-> y:", datas[1])
        print("-> z:", datas[2])
        plt.xlabel(parameters[0])
        plt.ylabel(parameters[1])
        plt.title(parameters[2])
        if xmin is not None:
            plt.xlim(left=xmin)
        if xmax is not None:
            plt.xlim(right=xmax)
        if ymin is not None:
            plt.ylim(bottom=ymin)
        if ymax is not None:
            plt.ylim(top=ymax)
        if draw_grid:
            plt.grid()
        try:
            plt.scatter(
                datas[0], datas[1], c=datas[2], cmap=cmap, vmin=zmin, vmax=zmax,
            )
        except Exception as e:
            print("Plot did not succeed.  Error message:")
            print(e)
        else:
            plt.colorbar()
            if save_as == "PDF":
                figid = 0
                filename = f"figure_{figid}.pdf"
                while os.path.exists(filename):
                    figid += 1
                    filename = f"figure_{figid}.pdf"
                print(f'Saving plot as "{filename}".')
                plt.savefig(filename)
            elif save_as == "PNG":
                figid = 0
                filename = f"figure_{figid}.png"
                while os.path.exists(filename):
                    figid += 1
                    filename = f"figure_{figid}.png"
                print(f'Saving plot as "{filename}".')
                plt.savefig(filename)
            plt.show()

    def complete_scatter(self, text, line, begidx, endidx):
        if not self.selected_table:
            return []
        parameters = ["cmap=", "grid", "png", "pdf"]
        parameters += self.con.get_column_names(self.selected_table)
        parameters.extend(extra_tools.keys())
        return [p for p in parameters if p.startswith(text)]

    def help_scatter(self):
        print("""> scatter <parameter_1> <parameter_2> <parameter_3> [condition] [-cmap=<name>] [-grid] [-{x|y|y}{min|max}=<number>] [-{png|pdf}]
    Make a scatter plot showing the values of parameter_3 in colour in relation
    to parameter_1 and parameter_2 for the selected experiment.

    The parameters can be read from the database or calculated using shell-
    extensions.  Type "help calculate" for further information on extensions.

    A condition can be used to filter or sort the data.  Type "help filter" for
    further information on filtering and sorting.

    To specify the colour map of plot, use "-cmap=<name>", where <name> is the
    name of a colour_map for matplotlib like "-cmap=bwr" to plot in blue-white-
    red or "-cmap=jet" for a colourful rainbow.  For further examples, consider
    https://matplotlib.org/users/colormaps.html .
    If multiple colour maps are specified, only the first one is considered.

    To draw a grid on the plot, add "-grid" to the command.

    To specify the range in either x-, y- or z-direction, use the attributes
    -xmin=, -xmax=, -ymin=, -ymax=, -zmin=, -zmax= with a number as argument.
    The z-axis refers to the colourbar.

    Add either "-png" or "-pdf" to the command to save the figure in a file of
    the corresponding filetype.  If both are specified, only the first argument
    is taken into account.

    The order of the qualifiers starting with a dash ("-") does not play a role,
    but they must be at the end of the command.

    While the window with the plot is open, the shell can be used as usual.
    If another plot command is executed with the figure still open, the new plot
    is superimposed onto the previous one.  Matplotlib automatically selects a
    new colour if no colour is specified explicitely.

    Example:
     - scatter intensity sigma duration diffusion="False" OR Kdiff=0 -zmin=0
       This command creates a scatter plot with the default colourmap of
       matplotlib on a colour-axis starting from 0 which shows the integration
       time in relation to the intensity and the value of sigma for experiments
       of the current class without diffusion.""")

    def do_pcolor(self, params):
        # Check the run condition
        if not self.selected_table:
            print('No experiment selected.  '
                  'Select an experiment to make a pseudocolour plot of its data.')
            return
        # Parse the arguments
        param_list = params.split(" ")
        # Two parameters are necessary
        if len(param_list) < 3:
            print('No three parameters given.  Type "help pcolor" for further information.')
            return
        parameters = []
        for i in range(3):
            parameters.append(param_list.pop(0))
        # The last parameters can be used to modify the plot
        save_as = None
        draw_grid = False
        shading = 'flat'
        cmap = None
        xmin = None
        xmax = None
        ymin = None
        ymax = None
        zmin = None
        zmax = None
        for param in reversed(param_list):
            if param == "-png":
                save_as = "PNG"
                param_list.pop()
            elif param == "-pdf":
                save_as = "PDF"
                param_list.pop()
            elif param == "-grid":
                draw_grid = True
                param_list.pop()
            elif param == "-shading":
                shading = 'gouraud'
                param_list.pop()
            elif param.startswith("-cmap="):
                cmap = param[6:]
                param_list.pop()
            elif param.startswith("-xmin="):
                try:
                    xmin = float(param[6:])
                except ValueError:
                    print(f'Value for xmin cannot be converted to float: "{xmin}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-xmax="):
                try:
                    xmax = float(param[6:])
                except ValueError:
                    print(f'Value for xmax cannot be converted to float: "{xmax}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-ymin="):
                try:
                    ymin = float(param[6:])
                except ValueError:
                    print(f'Value for ymin cannot be converted to float: "{ymin}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-ymax="):
                try:
                    ymax = float(param[6:])
                except ValueError:
                    print(f'Value for ymax cannot be converted to float: "{ymax}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-zmin="):
                try:
                    zmin = float(param[6:])
                except ValueError:
                    print(f'Value for zmin cannot be converted to float: "{zmin}" '
                          '-- ignoring it.')
                param_list.pop()
            elif param.startswith("-zmax="):
                try:
                    zmax = float(param[6:])
                except ValueError:
                    print(f'Value for zmax cannot be converted to float: "{zmax}" '
                          '-- ignoring it.')
                param_list.pop()
            else:
                # End of extra parameters, beginning of the SQL statement
                break
        # Every other parameter is an SQL condition to filter the data
        condition = " ".join(param_list).strip()
        # Get the data belonging to the parameters
        datas = []
        for parameter in parameters:
            data = self.get_data(self.selected_table, parameter, condition)
            if not data:
                print(f'Get data for parameter "{parameter}" ', end="")
                if condition:
                    print(f'under condition "{condition}" ', end="")
                print('failed.')
                return
            datas.append(data)
        # Arrange the data in a grid
        xvalues = sorted(set(datas[0]))
        yvalues = sorted(set(datas[1]))
        print(f'There are {len(xvalues)} unique x-values and {len(yvalues)} unique y-values.')
        data_grid = np.full((len(yvalues), len(xvalues)), np.nan)
        for i, zval in enumerate(datas[2]):
            data_grid[yvalues.index(datas[1][i]), xvalues.index(datas[0][i])] = zval
        print(f"-> x ({parameters[0]}):", xvalues)
        print(f"-> y ({parameters[1]}):", yvalues)
        print(f"-> z ({parameters[2]}):")
        print(data_grid)
        if shading == 'flat':
            xdiff2 = np.diff(xvalues)/2
            xaxis = [xvalues[0] - xdiff2[0], *(xvalues[:-1] + xdiff2), xvalues[-1] + xdiff2[-1]]
            ydiff2 = np.diff(yvalues)/2
            yaxis = [yvalues[0] - ydiff2[0], *(yvalues[:-1] + ydiff2), yvalues[-1] + ydiff2[-1]]
        else:
            xaxis = xvalues
            yaxis = yvalues
        plt.xlabel(parameters[0])
        plt.ylabel(parameters[1])
        plt.title(f'{parameters[2]} ({condition if condition else "all data"})'  )
        if xmin is not None:
            plt.xlim(left=xmin)
        if xmax is not None:
            plt.xlim(right=xmax)
        if ymin is not None:
            plt.ylim(bottom=ymin)
        if ymax is not None:
            plt.ylim(top=ymax)
        if draw_grid:
            plt.grid()
        try:
            plt.pcolormesh(
                xaxis, yaxis, data_grid, cmap=cmap, vmin=zmin, vmax=zmax,
                shading=shading,
            )
        except Exception as e:
            print("Plot did not succeed.  Error message:")
            print(e)
        else:
            plt.colorbar()
            if save_as == "PDF":
                figid = 0
                filename = f"figure_{figid}.pdf"
                while os.path.exists(filename):
                    figid += 1
                    filename = f"figure_{figid}.pdf"
                print(f'Saving plot as "{filename}".')
                plt.savefig(filename)
            elif save_as == "PNG":
                figid = 0
                filename = f"figure_{figid}.png"
                while os.path.exists(filename):
                    figid += 1
                    filename = f"figure_{figid}.png"
                print(f'Saving plot as "{filename}".')
                plt.savefig(filename)
            plt.show()

    def complete_pcolor(self, text, line, begidx, endidx):
        if not self.selected_table:
            return []
        parameters = ["cmap=", "grid", "png", "pdf", "shading"]
        parameters += self.con.get_column_names(self.selected_table)
        parameters.extend(extra_tools.keys())
        return [p for p in parameters if p.startswith(text)]

    def help_scatter(self):
        print("""> scatter <parameter_1> <parameter_2> <parameter_3> [condition] [-cmap=<name>] [-grid] [-{x|y|y}{min|max}=<number>] [-{png|pdf}]
    Make a scatter plot showing the values of parameter_3 in colour in relation
    to parameter_1 and parameter_2 for the selected experiment.

    The parameters can be read from the database or calculated using shell-
    extensions.  Type "help calculate" for further information on extensions.

    A condition can be used to filter or sort the data.  Type "help filter" for
    further information on filtering and sorting.

    To specify the colour map of plot, use "-cmap=<name>", where <name> is the
    name of a colour_map for matplotlib like "-cmap=bwr" to plot in blue-white-
    red or "-cmap=jet" for a colourful rainbow.  For further examples, consider
    https://matplotlib.org/users/colormaps.html .
    If multiple colour maps are specified, only the first one is considered.

    To draw a grid on the plot, add "-grid" to the command.

    To specify the range in either x-, y- or z-direction, use the attributes
    -xmin=, -xmax=, -ymin=, -ymax=, -zmin=, -zmax= with a number as argument.
    The z-axis refers to the colourbar.

    Add either "-png" or "-pdf" to the command to save the figure in a file of
    the corresponding filetype.  If both are specified, only the first argument
    is taken into account.

    The order of the qualifiers starting with a dash ("-") does not play a role,
    but they must be at the end of the command.

    While the window with the plot is open, the shell can be used as usual.
    If another plot command is executed with the figure still open, the new plot
    is superimposed onto the previous one.  Matplotlib automatically selects a
    new colour if no colour is specified explicitely.

    Example:
     - scatter intensity sigma duration diffusion="False" OR Kdiff=0 -zmin=0
       This command creates a scatter plot with the default colourmap of
       matplotlib on a colour-axis starting from 0 which shows the integration
       time in relation to the intensity and the value of sigma for experiments
       of the current class without diffusion.""")

    ### Functionality to CLEAN up
    def do_remove(self, params):
        global display_all
        display_all = True

        # Check conditions and parameters
        if not self.selected_table:
            print('No experiment selected.  Select an experiment to remove entries from it.')
            return
        if not params:
            print('No ids given.  Specify ids of the entries to remove.')
            return

        # Parse parameters
        ids = set()
        for p in params.split():
            try:
                id_ = int(p)
            except ValueError:
                print('Parameter is not a valid id: ' + p + '.  No data removed.')
                return
            if not self.con.entry_exists(self.selected_table, id_):
                print('No entry with id', id_, 'exists in the selected table.  No data removed.')
                return
            ids.add(id_)

        # Print full information of selected entries
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
        else:
            print('Answer was not "yes".  No data removed.')

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


def string_format_table(table_name, columns, rows):
    # Get current date and time to make datetime easy to read
    if not ISO_DATETIME:
        dt_now = datetime.datetime.now()
    # Convert rows to list of lists (because fetchall() returns a list of tuples)
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
                    if type(COMMENT_MAX_LENGTH) == int and COMMENT_MAX_LENGTH > 0:
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
                        val = val[:-10]
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
            and "comment" not in hidden_information
            and not display_all
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
    # Add column name centred, except the last one
    text += LIMITER.join([
        ("{:^" + str(l) + "}").format(n)
        for (n, t), l in zip(columns[:-1], lengths[:-1])
    ] + [columns[-1][0]]) + "\n"
    format_strings = []
    for (n, t), l in zip(columns, lengths):
        # Numbers right justified,
        # text centred,
        # comments unformatted
        if n in ["size_diag", "size_flux", "size_his", "size_mp4", "size_total"]:
            format_strings.append(SIZE_FORMAT.replace(":", ":>"+str(l)))
        elif t == "REAL":
            format_strings.append(FLOAT_FORMAT.replace(":", ":>"+str(l)))
        elif t == "INTEGER":
            format_strings.append("{:>" + str(l) + "}")
        elif n == "comment":
            format_strings.append("{:" + str(l) + "}")
        else:
            format_strings.append("{:^" + str(l) + "}")
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
param = Param(None)  # it is not necessary to specify a defaultfile for Param
datadir = param.datadir
del param
if datadir.startswith("~"):
    datadir = os.path.expanduser(datadir)

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

# Start the shell
ems_cli = EMShell(datadir)
ems_cli.cmdloop()
