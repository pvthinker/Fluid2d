import json
import os.path as path
import sys
import getopt

# Local import
import ems


class Param(object):
    """class to set up the default parameters value of the model

    the default parameters are stored in a json file along with
    their definition

    launch the code with the -h to get the help on all the model
    parameters, e.g.

    python vortex.py -h

    """

    def __init__(self, defaultfile=None, ems_file=""):
        """Load default parameters and optionally experiment parameters.

        The parameter `defaultfile` is no longer used and exists only
        for backwards compatibility.  The default file is always
        the file `defaults.json` located in the core-folder of fluid2d.

        The parameter `ems_file` takes optionally the name of an
        experiment file.  If given, the Experiment Management System
        (EMS) is activated and the experiment file is parsed.
        """

        import grid
        d = path.dirname(grid.__file__)
        jasonfile = d+'/defaults.json'

        with open(jasonfile) as f:
            namelist = json.load(f)

        self.set_parameters(namelist)

        opts, args = getopt.getopt(str(sys.argv[1:]), 'h:v', [''])
        if '-h' in args:
            self.manall()
            sys.exit()

        if '-v' in args:
            self.print_param = True
        else:
            self.print_param = False

        if ems_file:
            self.ems = ems.EMS(ems_file)
        else:
            self.ems = None

    def set_parameters(self, namelist):
        avail = {}
        doc = {}
        for d in namelist.keys():
            dd = namelist[d]
            for name in dd.keys():
                val = dd[name]['default']
                # print(name, val)
                setattr(self, name, val)
                if 'avail' in dd[name]:
                    avail[name] = dd[name]['avail']
                if 'doc' in dd[name]:
                    doc[name] = dd[name]['doc']
        self.avail = avail
        self.doc = doc

    def man(self, name):
        if name in self.doc:
            helpstr = self.doc[name]
            if name in self.avail:
                availstr = ', '.join([str(l) for l in self.avail[name]])
                helpstr += ' / available values = ['+availstr+']'
        else:
            helpstr = 'no manual for this parameter'

        name = '\033[0;32;40m' + name + '\033[0m'
        print('  - "%s" : %s\n' % (name, helpstr))

    def manall(self):
        ps = self.listall()
        for p in ps:
            self.man(p)

    def checkall(self):
        if self.ems:
            # Only create a new database entry once, not by every core
            if self.myrank == 0:
                self.ems.initialize(self.datadir)
            self.expname = self.ems.get_expname()
        for p, avail in self.avail.items():
            if getattr(self, p) in avail:
                # the parameter 'p' is well set
                pass
            else:
                msg = 'parameter "%s" should in ' % p
                msg += str(avail)
                raise ValueError(msg)

    def listall(self):
        """ return the list of all the parameters"""
        ps = [d for d in self.__dict__ if not(d in ['avail', 'doc'])]
        return ps

    def printvalues(self):
        """ print the value of each parameter"""
        for d in self.__dict__.keys():
            if not(d in ['avail', 'doc']):
                print('%20s :' % d, getattr(self, d))

    def copy(self, obj, list_param):
        """ copy attributes listed in list_param to obj

        On output it returns missing attributes
        """
        missing = []
        for k in list_param:
            if hasattr(self, k):
                setattr(obj, k, getattr(self, k))
            else:
                missing.append(k)
        return missing

    def get_experiment_parameters(self):
        """Return the experiment parameters dictionary loaded by the EMS.

        The EMS must be activated in the constructor of `Param` to use
        this method.  It returns the dictionary of experiment parameters
        and exits.  It is advised to use, when possible, the method
        `loop_experiment_parameters` instead.
        """
        return self.ems.parameters

    def loop_experiment_parameters(self):
        """Iterate over the experiment parameters loaded by the EMS.

        In every iteration, this method returns a new dictionary of
        experiment parameters containing a combination of the values
        specified in the experiment file.  This experiment file for the
        EMS must be specified in the constructor of `Param`.  If only
        one value is given for every parameter in the experiment file,
        the method `get_experiment_parameters` can be used instead.
        The ID of the experiment is increased in every iteration.
        """
        while self.ems.parameters:
            yield self.ems.parameters
            self.ems.setup_next_parameters()

    def finalize(self, fluid2d):
        """Invoke the finalize method of the EMS if activated."""
        if self.ems:
            self.ems.finalize(fluid2d)


if __name__ == "__main__":
    param = Param('default.xml')
    print('liste of parameters')
    print(param.listall())

    # to have the documentation on one particular parameter
    param.man('beta')

    # to get the documentation on all the parameters
    param.manall()

    # to check that all parameters that should a value taken from a list
    # have an acceptable value
    param.checkall()
