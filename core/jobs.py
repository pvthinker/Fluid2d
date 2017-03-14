from restart import Restart
import os
from importlib import import_module

class Jobs(object):
    def __init__(self,expname):
        
        self.logfile = '%s_log_jobs'%expname


        if os.path.isfile(self.logfile):
            # the log exists, we resume from the last job
            f = open(self.logfile,'r')
            lines = f.readlines()
            ele = lines[-1].split(' ')
            self.kjob = ele[0]
            f.close()
        else:
            # let's create the log and start the first job
            self.kjob = 0
            f = open(self.logfile,'w')
            f.writelines('%i'%self.kjob)
            f.close()
        
        #self.e = import_module( expname ) 
        #self.e.init()
        script = "%s.py"%expname
        #execfile(script)

        import subprocess
        
        subprocess.call(script)
        self.f2d = f2d

        if self.e.param.myrank==0:
            if self.kjob==0:
                print('start new experiment')
            if self.kjob>0:
                print('resume experiment from restart %i',self.kjob)

    def launch(self):

        f2d = self.f2d

        self.restart=Restart(f2d.param,f2d.grid)

        if self.kjob >0:
            # overwrite initial conditions
            t,dt,kt,tnextdiag,tnexthis = f2d.restart.read(f2d.model.var)
            f2d.t = t
            f2d.dt = dt
            f2d.kt = kt
            f2d.output.tnextdiag = tnextdiag
            f2d.output.tnexthis  = tnexthis
               
            f2d.tend += f2d.param.tend
            
            print('t = %f / kt = %i / nextdiag = %f / nexthist =%f'%(t,kt,tnextdiag,tnexthis))


        f2d.loop()

        self.restart.write(f2d.t,f2d.dt,f2d.kt,f2d.output.tnextdiag,f2d.output.tnexthis,self.model.var)
