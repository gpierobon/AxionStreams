import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Line3D
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from galpy import potential,df
from galpy.orbit import Orbit
from astropy import units as u

plt.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Export units 
kpc = u.kpc
kms = u.km/u.s
deg = u.deg
Gyr = u.Gyr

def plot_settings(ax,xlim,ylim,zlim):
    ax.set_xlim3d([-xlim,xlim])
    ax.set_ylim3d([-ylim,ylim])
    ax.set_zlim3d([-zlim,zlim])
    ax.set_xlabel('Galactic $X$ [kpc]',fontsize=20,labelpad=30)
    ax.set_ylabel('Galactic $Y$ [kpc]',fontsize=20,labelpad=30)
    ax.set_zlabel('Galactic $Z$ [kpc]',fontsize=20,labelpad=30)
    ax.tick_params(which='major',direction='in',width=3,length=10,right=True,top=True,pad=10,labelsize=20)
    ax.tick_params(which='minor',direction='in',width=1,length=7,right=True,top=True)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0
    ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0
    ax.grid(linestyle='--', linewidth=0.5)


def MyGrid(ax):
    ax.set_axis_off()
    ax.set_facecolor('black')
    xline = np.array([(0,0,0),(-100,0,0)])
    yline = np.array([(0,0,0),(0,100,0)])
    zline = np.array([(0,0,0),(0,0,100)])
    xline = Line3D(xline[:,0], xline[:,1],xline[:,2],color='w',alpha=0.15)
    yline = Line3D(yline[:,0], yline[:,1],yline[:,2],color='w',alpha=0.15)
    zline = Line3D(zline[:,0], zline[:,1],zline[:,2],color='w',alpha=0.15)
    ax.add_line(xline)
    ax.add_line(yline)
    ax.add_line(zline)
    theta = np.linspace(0, 2 * np.pi, 201)
    for i in range(0,100,10):
        x = i*kpc*np.cos(theta)
        y = i*kpc*np.sin(theta)
        ax.plot(x,y,0,c='w',alpha=0.5)
    

def plot_Sun(ax,label=False):
    Sun = np.array([8.122,0.0,0.005])
    o_sun1 = Orbit(vxvv=[Sun[0]*kpc,0.0*kms,232.0*kms,0.0*kpc,0.0*kms,0.0*deg]).flip()
    ts = np.linspace(0.0,250.0*u.Myr,200)
    o_sun1.integrate(ts,potential.MWPotential2014)
    Orb0 = np.column_stack((o_sun1.x(ts),o_sun1.y(ts),o_sun1.z(ts)))
    ax.plot(Orb0[:,0],Orb0[:,1],Orb0[:,2],'k-')
    if label == True:
        plt.gcf().text(0.15, 0.72,r'{\bf Sun}',fontsize=20,color='k')

def plot_single_orbit(orbit,azim=0,save=False,sun=True,disk=True,points=False,lim=10,size_x=10,size_y=10):
    fig = plt.figure(figsize=(size_x,size_y))
    ax = plt.axes(projection='3d')
    plot_settings(ax,lim,lim,lim)
    col = 'teal'
    ax.view_init(azim, 130)
    ax.plot(orbit[:,0],orbit[:,1],orbit[:,2],'-',alpha=1,c=col)
    if points == True:
        ax.scatter(orbit[:,0],orbit[:,1],orbit[:,2],'-',alpha=0.5,c=col)
    ax.plot(orbit[0,0],orbit[0,1],orbit[0,2],'-',marker='o',markersize=5,c=col)
    ax.plot(orbit[-1,0],orbit[-1,1],orbit[-1,2],'-',marker='X',markersize=5,c=col)

    if sun == True:
        plot_Sun(ax)

    if disk == True:
        theta = np.linspace(0, 2 * np.pi, 201)
        x = 30*kpc*np.cos(theta)
        y = 30*kpc*np.sin(theta)
        ax.plot(x,y,0,c='gray',ls='--',alpha=0.75)

    if save == True:
        fig.savefig('3D/single_orbit.png',bbox_inches='tight')
    else:
        plt.show() 

def single_plot_3d(mygrid=False,elev=10,rot=130,lim=30,size_x=10,size_y=10):
    fig = plt.figure(figsize=(size_x,size_y))
    ax = plt.axes(projection='3d')
    plot_settings(ax,lim,lim,lim)
    ax.view_init(elev,rot)
    if mygrid == True:
        ax.grid(False)
        MyGrid(ax)
    return fig,ax

def plot_stream(stream,mygrid=False,elev=10,rot=130,lim=30,size_x=10,size_y=10):
    fig,ax = single_plot_3d(mygrid=mygrid,elev=elev,rot=rot,lim=lim,size_x=size_x,size_y=size_y)
    ax.plot(stream.x,stream.y,stream.z)
    plt.show()

def plot_frame(streams,frame,mygrid=False,elev=10,rot=130,lim=30,size_x=10,size_y=10):
    fig,ax = single_plot_3d(mygrid=mygrid,elev=elev,rot=rot,lim=lim,size_x=size_x,size_y=size_y)
    
    x = np.array([streams[i].x[frame] for i in range(len(streams))])
    y = np.array([streams[i].y[frame] for i in range(len(streams))])
    z = np.array([streams[i].z[frame] for i in range(len(streams))])
    ax.scatter(x,y,z,c='k')
    
    return fig

def plot_frame_with_stream(streams,frame,mygrid=False,elev=10,rot=130,lim=30,sun=False,size_x=10,size_y=10):
    fig,ax = single_plot_3d(mygrid=mygrid,elev=elev,rot=rot,lim=lim,size_x=size_x,size_y=size_y)
    
    xr = [streams[i].x[:frame] for i in range(len(streams))]
    yr = [streams[i].y[:frame] for i in range(len(streams))]
    zr = [streams[i].z[:frame] for i in range(len(streams))]
    
    x = np.array([streams[i].x[frame] for i in range(len(streams))])
    y = np.array([streams[i].y[frame] for i in range(len(streams))])
    z = np.array([streams[i].z[frame] for i in range(len(streams))])
    
    col = 'k'
    alpha = 0.2

    if mygrid == True:
        col = 'teal'
        alpha = 1
    
    if sun == True:
        plot_Sun(ax)
    
    ax.scatter(x,y,z,c=col)
    for i in range(len(streams)):
        ax.plot(xr[i],yr[i],zr[i],alpha=alpha,c=col,lw=0.5)
    return fig



def plot_single_orbit_time(orbit,i,azim=0,save=False,sun=True,disk=True,points=False,lim=10,size_x=10,size_y=10):
    fig = plt.figure(figsize=(size_x,size_y))
    ax = plt.axes(projection='3d')
    plot_settings(ax,lim,lim,lim)
    col = 'teal'

def plot_single_orbit_time(orbit,i,azim=0,save=False,sun=True,disk=True,points=False,lim=10,size_x=10,size_y=10):
    fig = plt.figure(figsize=(size_x,size_y))
    ax = plt.axes(projection='3d')
    plot_settings(ax,lim,lim,lim)
    col = 'teal'
    ax.view_init(azim, 130)
    ax.plot(orbit[:i,0],orbit[:i,1],orbit[:i,2],'-',alpha=1,c=col)
    if points == True:
        ax.scatter(orbit[:,0],orbit[:,1],orbit[:,2],'-',alpha=0.5,c=col)
    ax.plot(orbit[0,0],orbit[0,1],orbit[0,2],'-',marker='o',markersize=5,c=col)
    #ax.plot(orbit[-1,0],orbit[-1,1],orbit[-1,2],'-',marker='X',markersize=5,c=col)
    if sun == True:
        plot_Sun(ax)

    if disk == True:
        theta = np.linspace(0, 2 * np.pi, 201)
        x = 30*kpc*np.cos(theta)
        y = 30*kpc*np.sin(theta)
        ax.plot(x,y,0,c='gray',ls='--',alpha=0.75)
    
    if save == True:
        fig.savefig('3D/fig_%.3d'%i,bbox_inches='tight')
    else:
        plt.show() 

def single_plot(xlab='',ylab='',\
                 lw=1.5,lfs=25,tfs=18,size_x=13,size_y=8,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)

    fig = plt.figure(figsize=(size_x,size_y))
    ax = fig.add_subplot(111)
    ax.set_xlabel(xlab,fontsize=lfs)
    ax.set_ylabel(ylab,fontsize=lfs)
    ax.tick_params(which='major',direction='in',width=0.8,length=8,right=True,top=True,pad=7)
    ax.tick_params(which='minor',direction='in',width=0.5,length=5,right=True,top=True)
    if Grid: ax.grid()
    return fig,ax
