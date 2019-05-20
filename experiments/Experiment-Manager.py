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


MP4_PLAYER = "mplayer"


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

    def show_table_overview(self, title=True):
        if title:
            print("Experiments in database:")
        for table in self.tables:
            print(
                " - {}: {} experiments, {} columns"
                .format(table.name, table.get_length(), len(table.columns))
            )

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
        COMMENT_MAX_LENGTH = 15
        LIMITER = "  "
        CUT_DATE = False
        # If CUT_DATE is True, CUT_MONTH is ignored
        CUT_MONTH = True
        # If CUT_MONTH is True, CUT_YEAR is ignored
        CUT_YEAR = False
        CUT_SECONDS = True
        # Get entries
        self.c.execute('SELECT * from "{}"'.format(self.name))
        rows = [list(row) for row in self.c.fetchall()]
        lengths = [len(c[0]) if c[0] != "comment" else COMMENT_MAX_LENGTH
                   for c in self.columns]
        table_size = 0
        for row in rows:
            for i, val in enumerate(row):
                # Check name of column
                if self.columns[i][0] == "comment":
                    # Cut comments that are too long and end them with an ellipsis
                    if len(val) > COMMENT_MAX_LENGTH:
                        row[i] = val[:COMMENT_MAX_LENGTH-1] + "…"
                    # No need to determine length of value
                    continue
                elif self.columns[i][0] == "datetime":
                    if CUT_DATE:
                        val = val.split("T")[1]
                    elif CUT_MONTH:
                        val = val[8:]
                    elif CUT_YEAR:
                        val = val[5:]
                    if CUT_SECONDS:
                        val = val[:-10]
                    val = val.replace("T", ",")
                    row[i] = val
                elif self.columns[i][0] == "size_total":
                    table_size += val

                # Check type of column
                if self.columns[i][1] == "TEXT" or self.columns[i][1] == "INTEGER":
                    lengths[i] = max(len(str(val)), lengths[i])
                elif self.columns[i][1] == "REAL":
                    lengths[i] = max(len("{:.3f}".format(val)), lengths[i])
                else:
                    raise Exception(
                        "unknown type {} of column {} in table {}."
                        .format(*self.columns[i], self.name)
                    )
        text = "Experiment: " + self.name + " ({:.3f} MB)\n".format(table_size)
        text_cols = []
        for i, (n, t) in enumerate(self.columns):
            text_cols.append(("{:^" + str(lengths[i]) + "}").format(n))
        text += LIMITER.join(text_cols) + "\n"
        format_strings = []
        for i, l in enumerate(lengths):
            if self.columns[i][1] == "REAL":
                format_strings.append("{:^" + str(l) + ".3f}")
            elif self.columns[i][0] == "comment":
                format_strings.append("{:" + str(l) + "}")
            else:
                format_strings.append("{:^" + str(l) + "}")
        for row in rows:
            text_cols = [f_str.format(val) for f_str, val in zip(format_strings, row)]
            text += LIMITER.join(text_cols) + "\n"
        return text.strip()

    def get_length(self):
        self.c.execute('SELECT Count(*) FROM "{}"'.format(self.name))
        return self.c.fetchone()[0]

    def entry_exists(self, id_):
        self.c.execute('SELECT id FROM "{}" WHERE id = ?'.format(self.name), (id_,))
        if self.c.fetchone() is not None:
            return True
        return False


# https://docs.python.org/3/library/cmd.html
class EMShell(cmd.Cmd):
    """Experiment Management (System) Shell"""
    
    intro = (
        "-" * 50
        + "\nFluid2d Experiment Management System (EMS)\n"
        + "-" * 50
        + "\nType help or ? to list available commands."
        + "\nType exit or Ctrl+D or Ctrl+C to exit."
        + "\nPress Tab-key for auto-completion of commands or table names.\n"
    )
    prompt = "(EMS) "

    def __init__(self, experiments_dir: str):
        super().__init__()
        self.exp_dir = experiments_dir
        self.con = EMDBConnection(os.path.join(self.exp_dir, "experiments.db"))
        self.selected_table = ""

    ### Functionality to SHOW the content of the database
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
        os.system('{} "{}" &'.format(MP4_PLAYER, path))

    def complete_open_mp4(self, text, line, begidx, endidx):
        return self.table_name_completion(text)

    def do_open_his(self, params):
        """Open the his-file for an experiment specified by its name and ID."""
        expname = self.parse_params_to_experiment(params)
        if not expname:
            return
        path = os.path.join(self.exp_dir, expname, expname + "_his.nc")
        if not os.path.isfile(path):
            print("File does not exist:", path)
            return
        os.system('ncview "{}" &'.format(path))

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
        os.system('ncview "{}" &'.format(path))

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
        self.con.show_table_overview()

    ### Helper functions
    def table_name_completion(self, text):
        completions = []
        for table in self.con.tables:
            if table.name.startswith(text):
                completions.append(table.name)
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
