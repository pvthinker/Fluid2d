import subprocess
import os
import numpy as np
import matplotlib
from sys import platform as _platform
font = {'size': 16}

matplotlib.use('TkAgg')
matplotlib.rc('font', **font)

# the backend 'TkAgg' has to be set before pyplot is imported
import matplotlib.pyplot as plt


class Plotting(object):
    def __init__(self, param, grid, var, diag):
        self.list_param = ['nh', 'cax', 'plot_var', 'colorscheme',
                           'Lx', 'Ly', 'expname', 'nx', 'ny',
                           'imshow_interpolation',
                           'plot_psi', 'cmap', 'plot_pvback',
                           'plot_ua', 'nh', 'expdir',
                           'generate_mp4', 'myrank', 'modelname',
                           'adaptable_dt', 'cfl', 'timestepping', 'order']
        param.copy(self, self.list_param)

        nh = self.nh
        nx = param.nx
        ny = param.ny

        dx2 = grid.dx/2
        self.xmin, self.xmax = grid.xr[nh, nh]-dx2, grid.xr[0, -nh]+dx2
        self.ymin, self.ymax = grid.yr[nh, nh]-dx2, grid.yr[-nh, 0]+dx2

        cmap_list = plt.cm.cmap_d.keys()
        if self.cmap in cmap_list:
            self.cmap_obj = plt.get_cmap(self.cmap)
        else:
            raise ValueError('cmap should be in ', cmap_list)

        # solid cells (nan) look gray
        self.cmap_obj.set_bad('#AAAAAA', 1.)

        self.msk = grid.msk[nh:-nh, nh:-nh]*1.
        self.msk[self.msk == 0] = np.nan
        self.z2d = var.get(self.plot_var)[nh:-nh, nh:-nh]

        if self.plot_psi:
            self.xp = grid.xr[nh:-nh, nh:-nh]+0.5*grid.dx
            self.yp = grid.yr[nh:-nh, nh:-nh]+0.5*grid.dx
            self.psi = var.get('psi')[nh:-nh, nh:-nh]

        if self.plot_pvback:
            self.xr = grid.xr[nh:-nh, nh:-nh]
            self.yr = grid.yr[nh:-nh, nh:-nh]

        if self.plot_ua:
            skip = 4
            self.U = var.get('ua')[nh:-nh:skip, nh:-nh:skip]
            self.V = var.get('va')[nh:-nh:skip, nh:-nh:skip]
            self.xQ = grid.xr[nh:-nh:skip, nh:-nh:skip]
            self.yQ = grid.yr[nh:-nh:skip, nh:-nh:skip]

    def create_fig(self, t):
        """create the figure for interactive plotting, we use imshow(), the
        axes are given in user scale (Lx, Ly) thanks to the 'extent' kwarg

        """

        if self.nx >= self.ny:
            width, height = 10*self.ny/self.ny, 8
        else:
            width, height = 11*self.nx/self.ny, 8

        fig = plt.figure(figsize=(width, height), facecolor='white')
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.cla()

        self.time_str = self.plot_var+' / time = %-6.2f\n'
        if self.adaptable_dt:
            self.time_str += 'cfl=%g / %s / order=%i' % (
                self.cfl, self.timestepping, self.order)
        else:
            self.time_str += 'dt=%g / %s / order=%i' % (
                self.dt, self.timestepping, self.order)

        ax1.set_title(self.time_str % t)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')

        self.ax1 = ax1


        self.im = ax1.imshow(self.z2d*self.msk,
                             vmin=self.cax[0], vmax=self.cax[1],
                             cmap=self.cmap_obj,
                             origin='lower',
                             extent=[self.xmin, self.xmax,
                                     self.ymin, self.ymax],
                             interpolation=self.imshow_interpolation)

        cb = plt.colorbar(self.im, ax=self.ax1)

        if self.plot_psi:
            self.psicolor = '#556b2f'  # psi contours are green
            maxpsi = np.max(np.abs(self.psi.ravel()))
            self.psilevs = np.linspace(-1., 1., 21)*maxpsi*1.1
            self.psilevs[self.psilevs == 0.] = maxpsi*1e-3

            extent = [self.xmin, self.xmax, self.ymin, self.ymax]
            self.cont = self.ax1.contour(self.xp, self.yp, self.psi,
                                         self.psilevs,
                                         colors=self.psicolor,
                                         extent=extent)

        if self.plot_pvback:
            nh = self.nh
            maxpvback = np.max(self.pvback.ravel())
            minpvback = np.min(self.pvback.ravel())
            self.pvblevs = np.linspace(minpvback, maxpvback, 6)
            self.pvblevs[self.pvblevs == 0.] = (maxpvback-minpvback)*1e-3
            ax1.contour(self.xr, self.yr, self.pvback[nh:-nh, nh: -nh],
                        self.pvblevs,
                        colors='k', linestyles=':', alpha=0.5,
                        extent=[self.xmin, self.xmax, self.ymin, self.ymax])

        self.set_cax(self.z2d)
        fig.show()
        if _platform in ["linux", "linux2"]:
            fig.canvas.draw()
        else:
            plt.pause(1e-4)
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
        """update the figure during the interactive plotting to speed up the
        animation, graphical objects are updated, not recreated.

        """
        self.im.set_array(self.z2d*self.msk)
        self.ax1.set_title(self.time_str % t)
        self.set_cax(self.z2d)

        self.im.set_clim(vmin=self.cax[0], vmax=self.cax[1])

        if self.plot_ua:
            self.Q.set_UVC(self.U, self.V)

        if self.plot_psi:
            # remove the previous contours
            for c in self.cont.collections:
                c.remove()
            # draw new ones
            self.cont = self.ax1.contour(self.xp, self.yp, self.psi,
                                         self.psilevs,
                                         colors=self.psicolor)

        # self.fig.canvas.draw() does not work with MacOSX, instead
        # use plt.pause() BUT plt.pause() prevents to use the other windows
        # during the animation
        if _platform in ["linux", "linux2"]:
            self.fig.canvas.draw()
            plt.pause(1e-4)            
        else:
            plt.pause(1e-4)


        if self.generate_mp4 and self.myrank == 0:
            string = self.fig.canvas.tostring_argb()
            self.process.stdin.write(string)

    def finalize(self):
        """ do nothing but close the mp4 thread if any"""
        if self.generate_mp4:
            self.process.communicate()

    def set_cax(self, z):
        """ set self.cax, the color range"""
        if self.colorscheme == 'minmax':
            self.cax = [min(z.ravel()), max(z.ravel())]

        if self.colorscheme == 'symmetric':
            mm = max(abs(z.ravel()))
            self.cax = [-mm, +mm]

        if self.colorscheme == 'imposed':
                #  self.cax is already defined, imposed by the user
            pass
