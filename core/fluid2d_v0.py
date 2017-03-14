from param import Param
from diagnostics import Diag
from plotting import Plotting
from output import Output
from numpy import nan,sqrt,abs

class Fluid2d(Param):
    def __init__(self,param,grid):

        self.list_param=['modelname','tend','fixed_dt','dt','cfl',
                         'plot_var','cax','colorscheme',
                         'plot_interactive','fixed_dt','dtmax',
                         'freq_save','freq_plot']

        param.copy(self,self.list_param)

        self.list_grid=['dx','nh','msk']
        grid.copy(self,self.list_grid)

        if param.modelname=='euler':
            from euler import Euler
            self.model = Euler(param,grid)

        if param.modelname=='advection':
            from advection import Advection
            self.model = Advection(param,grid)

        if param.modelname=='boussinesq':
            from boussinesq import Boussinesq
            self.model = Boussinesq(param,grid)

        if param.modelname=='quasigeostrophic':
            from quasigeostrophic import QG
            self.model = QG(param,grid)


        self.diag = Diag(param,grid)
        self.plotting = Plotting(param)
        
        # here's a shortcut to the model state
        self.state = self.model.var.state
        
        self.t  = 0.
        self.kt = 0
        
        self.output=Output(param,grid,self.diag)


    def update_fig(self,*args):
        self.diag.integrals(self.model.var)
#        print(self.diag.ke)
        self.set_dt()
        self.model.step(self.t,self.dt)
        self.kt +=1
        self.t = self.t+self.dt
#        print(self.t,self.dt)
        if(self.kt%self.freq_plot)==0:#self.plot_freq)==0:
            self.z2d[self.plotting_msk==0]=0.
            #            self.cax=self.get_cax(self.z2d)
            self.set_cax(self.z2d)
            self.im.set_array( self.z2d )
            self.im.set_clim(vmin=self.cax[0],vmax=self.cax[1])


        #self.output.do(self.model.var.state,self.diag,self.t,self.kt)

        return self.im,

    def update_fig2(self,*args):
        self.diag.integrals(self.model.var)
#        print(self.diag.ke)
        self.set_dt()
        self.model.step(self.t,self.dt)
        self.kt +=1
        self.t = self.t+self.dt
#        print(self.t,self.dt)
        if(self.kt%self.freq_plot)==0:#self.plot_freq)==0:
            self.z2d[self.plotting_msk==0]=0.
            self.set_cax(self.z2d)
            self.im.set_array( self.z2d )
            self.im.set_clim(vmin=self.cax[0],vmax=self.cax[1])


        self.output.do(self.model.var.state,self.diag,self.t,self.kt)
        self.ti.label='%f4.2'%self.t
        #print(self.ti.label)
        
        return self.im, self.ti,

    def loop(self):
        nh = self.nh
        self.plotting_msk=self.msk[nh:-nh,nh:-nh]
        self.z2d  = self.model.var.get(self.plot_var)[nh:-nh,nh:-nh]

        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
        if not(self.plot_interactive):
            

            fig  = plt.figure()
            self.ti = plt.title('')
            self.ti.animated=True
            self.im = plt.imshow( self.z2d, 
                                  vmin=self.cax[0],vmax=self.cax[1],
                                  cmap=plt.get_cmap('jet'),origin='lower',
                                  interpolation='nearest')

            plt.colorbar()
            ani = animation.FuncAnimation(fig,self.update_fig,interval=0.0000001,blit=True)
            plt.show()
            self.plotting.set_im(self.z2d,self.plotting_msk,self.cax,ani=True,update=self.update_fig)

        else:
            self.diag.integrals(self.model.var)
            time_str = 'time = %-6.2f'

            fig  = plt.figure(figsize=(10,8))
            ax1 = fig.add_subplot(1, 1, 1)
            ax1.cla()
            ax1.hold(True)

            ax1.set_title( time_str%self.t )
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')

            # ax2 = fig.add_subplot(1, 2, 2)
            # ax2.cla()
            # ax2.hold(True)
            # ax2.set_title( 'integral quantities')
            # ax2.set_xlabel('t')
            # ax2.set_ylabel('Ke')

            def on_press(event):
                print('stop')
                if self.ion:
                    plt.ioff()
                else:
                    plt.ion()
                plt.pause(1)
            self.ion = True

            #plt.ion()
            #plt.show()#False)
            im=ax1.imshow( self.z2d, 
                             vmin=self.cax[0],vmax=self.cax[1],
                             cmap=plt.get_cmap('jet'),origin='lower',
                             interpolation='nearest')

            time=[0]

            z = self.diag.ke
            y= [z]

            #line1, = ax2.plot(time, y, '.-', alpha=0.8, color="gray", markerfacecolor="red")
            cb = plt.colorbar(im)

            fig.show()
            fig.canvas.draw()
            background1 = fig.canvas.copy_from_bbox(ax1.bbox)

            cid = fig.canvas.mpl_connect('button_press_event',on_press)
            fig.canvas.key_press_event=on_press


            #background2 = fig.canvas.copy_from_bbox(ax2.bbox)
            self.output.do(self.model.var.state,self.diag,self.t,self.kt)

            while (self.t<self.tend):        
                self.diag.integrals(self.model.var)
                self.set_dt()
                self.model.step(self.t,self.dt)
                self.t  += self.dt
                self.kt += 1
                self.output.do(self.model.var.state,self.diag,self.t,self.kt)
                
                #print(self.t)
                if self.kt%self.freq_plot ==0:
                    time.append(self.t)
                    y.append(self.diag.ke)

                    if len(time)>1000:
                        del time[0]
                        del y[0]
                        
                    #line1.set_xdata(time)
                    #line1.set_ydata(y)

                    im.set_array( self.z2d )
                    ti = ax1.title
                    ax1.set_title( time_str%self.t )
                    self.set_cax(self.z2d)

                    im.set_clim(vmin=self.cax[0],vmax=self.cax[1])
                    if False:
                        fig.canvas.restore_region(background1)    # restore background
                        ax1.draw_artist(im)                   # redraw just the points
                        #ax1.draw_artist(ti)                   # redraw just the points
                        fig.canvas.blit(ax1.bbox)                # fill in the axes rectangle

                        ax2.axis([min(time),max(time),min(y),max(y)])
                        fig.canvas.restore_region(background2)    # restore background
                        ax2.draw_artist(line1)
                        fig.canvas.blit(ax2.bbox)                # fill in the axes rectangle
                    else:                    
                        fig.canvas.draw()

                    


    def set_dt(self):

        if ((self.fixed_dt == 0) & (self.diag.maxspeed != 0)):
            dt = self.cfl * self.dx / self.diag.maxspeed
            # filter in time
            self.filter_coef=1.
            self.dt = (1.-self.filter_coef)*self.dt + self.filter_coef*dt
            if self.dt>self.dtmax:
#                print('dt=%g'%self.dt)
                self.dt=self.dtmax
#            else:
#                print(self.dt)
        else:
            pass
        # transfer this information to the advection scheme
        self.model.ope.cst[2]=self.diag.maxspeed
#        print('dt=%g / cfl=%g'%(self.dt,self.cfl))


    def set_cax(self,z):
        if self.colorscheme=='minmax':
            self.cax = [min(z.ravel()),max(z.ravel())]
        if self.colorscheme=='symmetric':
            mm = max(abs(z.ravel()))
            self.cax = [-mm,+mm]
        if self.colorscheme=='imposed':
                #self.cax is already defined
            pass

