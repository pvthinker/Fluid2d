import matplotlib
font = {'size': 16}

matplotlib.use('TkAgg')
matplotlib.rc('font', **font)

import os
import subprocess
from numpy import nan
import matplotlib.pyplot as plt
from param import Param
from sys import platform as _platform


class Plotting(Param):
    def __init__(self, param, grid, var, diag):
        self.list_param = ['nh', 'cax', 'plot_var', 'colorscheme',
                           'Lx', 'Ly', 'expname', 'nx', 'ny',
                           'imshow_interpolation',
                           'plot_psi', 'cmap', 'plot_pvback',
                           'plot_ua', 'nh', 'expdir',
                           'generate_mp4', 'myrank', 'modelname']

        param.copy(self, self.list_param)

        nh = self.nh
        self.msk = grid.msk[nh:-nh, nh:-nh]
        self.qE = var.get('qE')[nh:-nh, nh:-nh]
        self.buoy = var.get('buoyancy')[nh:-nh, nh:-nh]
        self.Vjet = var.get('V')[nh:-nh, nh:-nh]
        self.vor = var.get('vorticity')[nh:-nh, nh:-nh]

    def create_fig(self, t):
        fig, axes = plt.subplots(
            nrows=2, ncols=2, figsize=(12, 8), facecolor='white')
        for i in range(2):
            for j in range(2):
                # axes[j, i].hold(True)
                axes[j, i].set_ylabel('Y')
                axes[j, i].set_xlabel('X')

        self.time_str = 'time = %-6.2f'
        self.txt = plt.annotate(self.time_str % t, xy=(
            0.5, .9), xycoords='figure fraction', fontsize=16)
#        axes[0,0].set_title( self.time_str%t  )

        self.axes = axes

        plt.ion()

        fs = 16
        axes[0, 0].set_title('V', fontsize=fs)
        axes[0, 1].set_title('b', fontsize=fs)
        axes[1, 0].set_title('Ertel PV', fontsize=fs)
        axes[1, 1].set_title(r'$\omega$', fontsize=fs)

        self.im0 = axes[0, 0].imshow(self.Vjet,
                                     vmin=self.cax[0], vmax=self.cax[1],
                                     cmap=plt.get_cmap('jet'), origin='lower',
                                     interpolation='nearest')

        self.im1 = axes[0, 1].imshow(self.buoy,
                                     vmin=self.cax[0], vmax=self.cax[1],
                                     cmap=plt.get_cmap('jet'), origin='lower',
                                     interpolation='nearest')

        self.im2 = axes[1, 0].imshow(self.qE,
                                     vmin=self.cax[0], vmax=self.cax[1],
                                     cmap=plt.get_cmap('jet'), origin='lower',
                                     interpolation='nearest')

        self.im3 = axes[1, 1].imshow(self.vor,
                                     vmin=self.cax[0], vmax=self.cax[1],
                                     cmap=plt.get_cmap('jet'), origin='lower',
                                     interpolation='nearest')

    #cb = plt.colorbar(self.im0)

        fig.show()
        fig.canvas.draw()
        self.fig = fig

        if self.generate_mp4 and (self.myrank == 0):
            canvas_width, canvas_height = fig.canvas.get_width_height()
            # Open an ffmpeg process
            outf = '%s/%s_%s.mp4' % (self.expdir, self.expname, self.plot_var)
            self.mp4file = outf
            videoencoder = None
            for v in ['avconv', 'ffmpeg']:
                if subprocess.call(['which', v], stdout=subprocess.PIPE) == 0:
                    videoencoder = v

            if videoencoder is None:
                print('\n')
                print('Neither avconv or ffmpeg was found')
                print('Install one of them or set param.generate_mp4 = False')
                exit(0)

            cmdstring = (videoencoder,
                         '-y', '-r', '30',  # overwrite, 30fps
                         # size of image string
                         '-s', '%dx%d' % (canvas_width, canvas_height),
                         '-pix_fmt', 'argb',  # format
                         '-f', 'rawvideo',
                         # tell ffmpeg to expect raw video from the pipe
                         '-i', '-',
                         '-vcodec', 'libx264', outf)  # output encoding

            devnull = open(os.devnull, 'wb')
            self.process = subprocess.Popen(cmdstring,
                                            stdin=subprocess.PIPE,
                                            stdout=devnull,
                                            stderr=devnull)

    def update_fig(self, t, dt, kt):

        self.txt.set_text(self.time_str % t)

        self.im0.set_array(self.Vjet)
        self.set_cax(self.Vjet)
        self.im0.set_clim(vmin=self.cax[0], vmax=self.cax[1])

        self.im1.set_array(self.buoy)
        self.set_cax(self.buoy)
        self.im1.set_clim(vmin=self.cax[0], vmax=self.cax[1])

        self.im2.set_array(self.qE)
        self.set_cax(self.qE)
        self.im2.set_clim(vmin=self.cax[0], vmax=self.cax[1])

        self.im3.set_array(self.vor)
        self.set_cax(self.vor)
        self.im3.set_clim(vmin=self.cax[0], vmax=self.cax[1])

        # self.fig.canvas.draw() does not work with MacOSX, instead
        # use plt.pause() BUT plt.pause() prevents to use the other windows
        # during the animation
        if _platform in ["linux", "linux2"]:
            self.fig.canvas.draw()
            plt.pause(1e-4)
        else:
            plt.pause(1e-4)

        if self.generate_mp4:
            string = self.fig.canvas.tostring_argb()
            self.process.stdin.write(string)

    def finalize(self):
        # Finish up
        if self.generate_mp4:
            self.process.communicate()

    def set_cax(self, z):
        if self.colorscheme == 'minmax':
            self.cax = [min(z.ravel()), max(z.ravel())]
        if self.colorscheme == 'symmetric':
            mm = max(abs(z.ravel()))
            self.cax = [-mm, +mm]
        if self.colorscheme == 'imposed':
                # self.cax is already defined
            pass
