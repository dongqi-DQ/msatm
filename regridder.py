import numpy as np
import matplotlib.pyplot as plt

def regrid(dat,x,y,nx,ny,make_plot=False):
    """
    dat - data array to regrid
    x - x values
    y - y values
    nx - new x vals
    ny - new y vals
    """
    xx = np.tile(  np.expand_dims(x,axis=1), (1,len(y)))
    yy = np.tile(  y, (len(x),1))

    nxx = np.tile(  np.expand_dims(nx,axis=1), (1,len(ny)))
    nyy = np.tile(  ny, (len(nx),1))

    interped = griddata((xx.flatten(),yy.flatten()),dat.flatten(),(nxx.flatten(),nyy.flatten()),method='linear').reshape((len(nx),len(ny)))

    if make_plot:
        clevs = np.linspace(293,313,21)
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax1.pcolor(yy,xx,dat,cmap='jet',vmin=150,vmax=310)
        ax1.set_title('original')
        ax1.invert_yaxis()

        ax2 = fig.add_subplot(122)
        ax2.pcolor(nyy,nxx,interped,cmap='jet',vmin=150,vmax=310)
        ax2.set_title('regridded')
        ax2.invert_yaxis()
        fig.show()

    return interped

def regrid_1d(dat,old,new,ax,make_plot=False):
    interpfunc = interp1d(old,dat,axis=ax)
    interped = interpfunc(new)
    
    return interped 


