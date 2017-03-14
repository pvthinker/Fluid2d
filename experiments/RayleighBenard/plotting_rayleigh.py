from param import Param
import matplotlib.pyplot as plt

class Plotting(Param):
    def __init__(self,param,grid,var,diag):
        self.list_param=['nh','cax','plot_var','colorscheme']
        param.copy(self,self.list_param)

        nh = self.nh
        self.msk = grid.msk[nh:-nh,nh:-nh]
        self.z2d = var.get(self.plot_var)[nh:-nh,nh:-nh]

    def create_fig(self):
        fig  = plt.figure(figsize=(16,6))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.cla()
        ax1.hold(True)

        self.time_str = 'time = %-6.2f'
        ax1.set_title( '' )
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')

        self.ax1 = ax1

        plt.ion()

        self.im=ax1.imshow( self.z2d, 
                       vmin=self.cax[0],vmax=self.cax[1],
                       cmap=plt.get_cmap('jet'),origin='lower',
                       interpolation='nearest')

        cb = plt.colorbar(self.im)

        fig.show()
        fig.canvas.draw()
        self.fig = fig

    def update_fig(self,t,dt,kt):
        
        self.im.set_array( self.z2d )
        self.ax1.set_title( self.time_str%t )
        self.set_cax(self.z2d)

        self.im.set_clim(vmin=self.cax[0],vmax=self.cax[1])
        self.fig.canvas.draw()

    def set_cax(self,z):
        if self.colorscheme=='minmax':
            self.cax = [min(z.ravel()),max(z.ravel())]
        if self.colorscheme=='symmetric':
            mm = max(abs(z.ravel()))
            self.cax = [-mm,+mm]
        if self.colorscheme=='imposed':
                #self.cax is already defined
            pass
