import matplotlib
font = {'size'   : 16}

matplotlib.use('TkAgg')
matplotlib.rc('font', **font)

from param import Param
import matplotlib.pyplot as plt
from numpy import nan
import subprocess
import os


class Plotting(Param):
    def __init__(self,param,grid,var,diag):
        self.list_param=['nh','cax','plot_var','colorscheme','Lx','Ly','expname','generate_mp4']
        param.copy(self,self.list_param)

        nh = self.nh

        self.msk = grid.msk[nh:-nh,nh:-nh]

        self.z2d = var.get(self.plot_var)[nh:-nh,nh:-nh]
        

    def create_fig(self,t):
        fig  = plt.figure(figsize=(10*self.Lx/self.Ly,8),facecolor='white')
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.cla()
        #ax1.hold(True) # <- deprecated in recent versions

        self.time_str = 'time = %-6.2f'
        ax1.set_title( self.time_str%t )
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')

        self.ax1 = ax1

        plt.ion()

        self.im=ax1.imshow( self.z2d, 
                       vmin=self.cax[0],vmax=self.cax[1],
                       cmap=plt.get_cmap('jet'),origin='lower',
                       interpolation='nearest')

        cb = plt.colorbar(self.im,ax=self.ax1)

        fig.show()
        fig.canvas.draw()
        self.fig = fig

        if self.generate_mp4:
            canvas_width, canvas_height = fig.canvas.get_width_height()
            # Open an ffmpeg process
            outf = self.expname+'.mp4'
            cmdstring = ('avconv', 
                '-y', '-r', '30', # overwrite, 30fps
                '-s', '%dx%d' % (canvas_width, canvas_height), # size of image string
                '-pix_fmt', 'argb', # format
                '-f', 'rawvideo',  '-i', '-', # tell ffmpeg to expect raw video from the pipe
                '-vcodec', 'libx264', outf) # output encoding

            devnull = open(os.devnull, 'wb')
            self.process = subprocess.Popen(cmdstring, stdin=subprocess.PIPE,stdout=devnull,stderr=devnull)


    def update_fig(self,t,dt,kt):
        
        self.im.set_array( self.z2d )
        self.ax1.set_title( self.time_str%t )
        self.set_cax(self.z2d)

        self.im.set_clim(vmin=self.cax[0],vmax=self.cax[1])
        self.fig.canvas.draw()
        if self.generate_mp4:
            string = self.fig.canvas.tostring_argb()
            self.process.stdin.write(string)

    def finalize(self):
        # Finish up
        if self.generate_mp4:
            self.process.communicate()


    def set_cax(self,z):
        if self.colorscheme=='minmax':
            self.cax = [min(z.ravel()),max(z.ravel())]
        if self.colorscheme=='symmetric':
            mm = max(abs(z.ravel()))
            self.cax = [-mm,+mm]
        if self.colorscheme=='imposed':
                #self.cax is already defined
            pass
