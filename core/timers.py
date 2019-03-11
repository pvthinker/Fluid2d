from time import time as clock


class Timers(object):
    def __init__(self, param):
        self.t0 = {}
        self.elapse = {}
        self.ncalls = {}

    def tic(self, name):
        if not(name in self.t0.keys()):
            self.elapse[name] = 0.
            self.ncalls[name] = 0.
        self.t0[name] = clock()

    def toc(self, name):
        dt = clock()-self.t0[name]
        self.elapse[name] += dt
        self.ncalls[name] += 1

    def _print(self):
        keys = self.elapse.keys()
        for key in sorted(keys):
            print('%10s : %6.2f s / %6i calls / %6.2e' %
                  (key, self.elapse[key], self.ncalls[key],
                   self.elapse[key]/self.ncalls[key]))
