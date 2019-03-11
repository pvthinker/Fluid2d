import numpy as np
from param import Param


class Var(object):
    def __init__(self, param):
        """Define a model state variable

        for which all variables are packed together
        in a larger single array
        the var also has a time axis available """

        self.list_param = ['sizevar', 'varname_list']
        param.copy(self, self.list_param)

        self.nvar = len(self.varname_list)

        if type(self.sizevar) != list:
            print('sizevar has to be a list')
            print('Abort!!!')
            exit(0)

        self.sizestate = [self.nvar]+self.sizevar

        self.state = np.zeros(self.sizestate)

    def get(self, name):
        """ extract variable 'name' from the stack """
        k = self.varname_list.index(name)
        return self.state[k]


if __name__ == "__main__":

    param = Param()
    param.varname_list = ['vort', 'psi', 'u']
    param.sizevar = [4, 2]

    var = Var(param)

    print(np.shape(var.state))

    vor = var.get('vort')
    vor[:, 0] = 1.
    print(np.shape(vor))

    print(var.state)
