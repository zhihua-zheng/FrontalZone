#!/usr/bin/env python3

import sys
import warnings
import argparse
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from viztool import FormatScalarFormatter
from matplotlib.colors import LinearSegmentedColormap

###patch start###
from mpl_toolkits.mplot3d.axis3d import Axis
if not hasattr(Axis, "_get_coord_info_old"):
    def _get_coord_info_new(self, renderer):
        mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(renderer)
        mins += deltas / 4
        maxs -= deltas / 4
        return mins, maxs, centers, deltas, tc, highs
    Axis._get_coord_info_old = Axis._get_coord_info  
    Axis._get_coord_info = _get_coord_info_new
###patch end###

def main():
    # process input arguments
    parser = argparse.ArgumentParser(description="""
            Animate 3D buoyancy field from Oceananigans simualtion.""")
    parser.add_argument('-f', '--file', action='store', dest='fname',
            metavar='FILENAME', help='Input simulation data')
    #parser.add_argument('-f2', '--file2', action='store', dest='fname2',
    #        metavar='FILENAME', help='Input GOTM data 2')
    #parser.add_argument('-v', '--variable', action='store', dest='vname',
    #        metavar='VARNAME', nargs='+', help='Variable name')
    parser.add_argument('-o', '--output', action='store', dest='fname_out',
            metavar='FIGNAME', help='Output figure name')
    #parser.add_argument('-ds', '--date_start', action='store', dest='date_start',
    #        metavar='STARTDATE',
    #        help='Starting date of input data, in the format of YYYYMMDD')
    #parser.add_argument('-de', '--date_end', action='store', dest='date_end',
    #        metavar='ENDDATE',
    #        help='Ending date of input data, in the format of YYYYMMDD')
    parser.add_argument('--version', action='version', version='%(prog)s: 1.0')
    # parsing arguments and save to args
    args = parser.parse_args()
    
    # check input
    if not args.fname or not args.fname_out:# or not args.vname
        print('Oceananigans netCDF data and output figure name are required. Stop.\n')#, variable name,
        parser.print_help()
        sys.exit(1)
    
    # specify file path
    if sys.platform == 'linux' or sys.platform == 'linux2':
        data_dir = '/glade/work/zhihuaz/Data/FrontalZone/'
        figs_dir = '/glade/u/home/zhihuaz/Projects/TRACE-SEAS/FrontalZone/Figures/'
    elif sys.platform == 'darwin':
        data_dir = '/Users/zhihua/Documents/Work/Research/Projects/TRACE-SEAS/FrontalZone/Data/'
        figs_dir = '/Users/zhihua/Documents/Work/Research/Projects/TRACE-SEAS/FrontalZone/Figures/'
    else:
        print('OS not supported.')

    # customized colormap
    colorlist = ['xkcd:navy', 'xkcd:cerulean', 'xkcd:pale grey', 'xkcd:melon', 'xkcd:deep red']
    nodes = [0.0, 0.75, 0.82, 0.89, 1.0]
    cmap = LinearSegmentedColormap.from_list('buoyancy', list(zip(nodes, colorlist)))

    # read data
    fpath = data_dir+args.fname
    top = xr.open_dataset(fpath, group='slice/top').load()#.chunk('auto')
    top.close()
    z_top_slice = top.zC[0]
    south = xr.open_dataset(fpath, group='slice/south').load().sel(zC=slice(None,z_top_slice))#.chunk('auto')
    south.close()
    y_south_slice = south.yC[0]
    east = xr.open_dataset(fpath, group='slice/east').load().sel(zC=slice(None,z_top_slice))#.chunk('auto')
    east.close()
    x_east_slice = east.xC[0]

    # total buoyancy
    top['bt']   = (-top.attrs['M²']   * top.xC   + top.b).squeeze().transpose('yC','xC',...)
    south['bt'] = (-south.attrs['M²'] * south.xC + south.b).squeeze().transpose('xC','zC',...)
    east['bt']  = (-east.attrs['M²']  * east.xC  + east.b).squeeze().transpose('yC','zC',...)

    # time interval and initial timestamp in unit of hours
    dhr = int(((top.time[2] - top.time[1])/np.timedelta64(1,'h')).data)
    hr0 = int((top.time[0]/np.timedelta64(1,'h')).data)

    X, Y, Z = np.meshgrid(top.xC, top.yC, east.zC)

    # limits and contour values
    bmin = min([top.bt.min(), east.bt.min(), south.bt.min()])
    bmax = max([top.bt.max(), east.bt.max(), south.bt.max()])
    bmin = np.ceil(bmin*1e6)/1e6
    bmax = np.floor(bmax*1e6)/1e6
    blines = np.concatenate([np.arange(0,1.4,0.2), np.arange(1.3,2,0.05)])*1e-4
    xmin, xmax = -500, 500
    ymin, ymax = 0, 1000
    zmin, zmax = -140, z_top_slice.data

    Ckw = {'vmin': bmin,
           'vmax': bmax,
           'levels': np.linspace(bmin, bmax, 256),
           'cmap': cmap
          }
    Lkw = {'linewidths': 0.8, 
       'colors': 'xkcd:almost black'
      }
    edges_kw = {'color': 'xkcd:charcoal', 
                'linewidth': 1.7, 
                'zorder': 2
               }
    
    fig = plt.figure(figsize=(6,4.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection='3d', computed_zorder=False)
    itime = 0
    xroll = 0
    
    Ct = ax.contourf(X[:, :, -1], Y[:, :, -1], top.bt.isel(time=itime).roll(xC=xroll), zdir='z', offset=z_top_slice, **Ckw)
    Cs = ax.contourf(X[0, :, :], south.bt.isel(time=itime).roll(xC=xroll), Z[0, :, :], zdir='y', offset=y_south_slice, **Ckw)
    # xroll doesn't work for east slice!!
    Ce = ax.contourf(east.bt.isel(time=itime), Y[:, -1, :], Z[:, -1, :], zdir='x', offset=x_east_slice, **Ckw)

    Lt = ax.contour(X[:, :, -1], Y[:, :, -1], top.bt.isel(time=itime).roll(xC=xroll), blines, zdir='z', offset=z_top_slice, **Lkw)
    Ls = ax.contour(X[0, :, :], south.bt.isel(time=itime).roll(xC=xroll), Z[0, :, :], blines, zdir='y', offset=y_south_slice, **Lkw)
    Le = ax.contour(east.bt.isel(time=itime), Y[:, -1, :], Z[:, -1, :], blines, zdir='x', offset=x_east_slice, **Lkw)

    ax.plot([xmax, xmax], [ymin, ymax], zmin, **edges_kw)
    ax.plot([xmin, xmax], [ymin, ymin], zmin, **edges_kw)
    ax.plot([xmax, xmax], [ymin, ymax], zmax, **edges_kw)
    ax.plot([xmin, xmax], [ymin, ymin], zmax, **edges_kw)
    ax.plot([xmin, xmin], [ymin, ymax], zmax, **edges_kw)
    ax.plot([xmin, xmax], [ymax, ymax], zmax, **edges_kw)
    ax.plot([xmax, xmax], [ymax, ymax], [zmin, zmax], **edges_kw)
    ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)
    ax.plot([xmin, xmin], [ymin, ymin], [zmin, zmax], **edges_kw)
    ax.set(xlabel='X [m]',
           ylabel='Y [m]',
           zlabel='Z [m]',
           zticks=[-10, -60, -120],
           xticks=[-400, -200, 0, 200, 400],
           yticks=[100, 300, 500, 700, 900],
           xlim=[xmin, xmax], 
           ylim=[ymin, ymax], 
           zlim=[zmin, zmax])
    ax.tick_params(axis='both', labelsize=8)
    ax.view_init(20, -70, 0)
    ax.set_box_aspect((1,1,0.58), zoom=1.17)
    ax.set_title(f't = {(itime*dhr + hr0)//24} days, {(itime*dhr + hr0)%24} hours', y=0.995);
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    ax.grid(False)

    cbar = fig.colorbar(Ct, ax=ax, fraction=0.025, pad=0.1, ticks=np.arange(0,2,0.2)*1e-4, 
                        format=FormatScalarFormatter('%.1f'))
    cbar.set_label(r'Buoyancy [m s$^{-2}$]', labelpad=-42, fontsize=10)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.formatter.set_useMathText(True)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.ax.tick_params(labelsize=8) 

    def update(frame):
        nonlocal Ct, Cs, Ce, Lt, Ls, Le
        # for each frame, get new data and clear old data stored on each artist
        top_field   = top.bt.isel(time=frame).roll(xC=xroll)
        south_field = south.bt.isel(time=frame).roll(xC=xroll)
        east_field  = east.bt.isel(time=frame)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for obj in (Ct, Cs, Ce, Lt, Ls, Le):
                for coll in obj.collections:
                    coll.remove()
        # update the plot
        Ct = ax.contourf(X[:, :, -1], Y[:, :, -1], top_field,   zdir='z', offset=z_top_slice,   **Ckw)
        Cs = ax.contourf(X[0, :, :], south_field,  Z[0, :, :],  zdir='y', offset=y_south_slice, **Ckw)
        Ce = ax.contourf(east_field, Y[:, -1, :],  Z[:, -1, :], zdir='x', offset=x_east_slice,  **Ckw)
        Lt = ax.contour(X[:, :, -1], Y[:, :, -1], top_field,   blines, zdir='z', offset=z_top_slice,   **Lkw)
        Ls = ax.contour(X[0, :, :],  south_field, Z[0, :, :],  blines, zdir='y', offset=y_south_slice, **Lkw)
        Le = ax.contour(east_field,  Y[:, -1, :], Z[:, -1, :], blines, zdir='x', offset=x_east_slice,  **Lkw)
        ax.set_title(f't = {(frame*dhr + hr0)//24} days, {(frame*dhr + hr0)%24} hours', y=0.995)
        return Ct, Cs, Ce, Lt, Ls, Le
    
    ani = animation.FuncAnimation(fig=fig, func=update, frames=range(top.dims['time']), 
                                  interval=400, repeat_delay=800)#, blit=True
    # save animation
    ani.save(figs_dir+args.fname_out, writer='ffmpeg', fps=3, dpi=200)


if __name__ == "__main__":
    main()