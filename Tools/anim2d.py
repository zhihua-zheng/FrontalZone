#!/usr/bin/env python3

import warnings
import argparse
import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from viztool import FormatScalarFormatter

def main():
    # process input arguments
    parser = argparse.ArgumentParser(description="""
        Animate top slice and cross front transect from Oceananigans simualtion. Accept multiple variables.""")
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
    
    data_dir = '/glade/work/zhihuaz/Data/FrontalZone/'
    figs_dir = '/glade/u/home/zhihuaz/Projects/TRACE-SEAS/FrontalZone/Figures/'
    
    # read data
    ds = xr.open_dataset(data_dir+args.fname).load()
    ds.close()
    
    # add background fields
    H = abs(ds.zF[0])
    ds['Bt'] = (-ds.M2 * ds.xC + ds.B).transpose('time','zC','xC')
    ds['Vt'] = (-ds.M2/ds.f * (ds.zC + H) + ds.V).transpose('time','zC','xC')
    
    # time interval in unit of hours
    ds = int(((ds.time[1] - ds.time[0])/np.timedelta64(1,'h')).data)
    
    fig = plt.figure(figsize=(4.8,6), constrained_layout=True)
    gs = gridspec.GridSpec(8,1, figure=fig, left=0.1, right=0.9, bottom=0.1, top=0.9,
                           wspace=0.05, hspace=0.05)

    ax1 = fig.add_subplot(gs[:6, 0])
    kw_w1 = {'norm': mcolors.CenteredNorm(),
           'levels': np.linspace(-3e-3,3e-3,31),
           'extend': 'both',
           'cmap': 'RdBu_r'
          }
    kw_b1 = {'levels': np.arange(12,17,0.5)*1e-5,
           'colors': 'k',
           'linewidths': 1
          }
    Cw1 = ax1.contourf(ds.xC, ds.yC, ds.w.isel(time=itime), **kw_w1)
    Cb1 = ax1.contour(ds.xC, ds.yC, ds.bt.isel(time=itime), **kw_b1)
    cbar1 = fig.colorbar(Cw1, ax=ax1, fraction=0.027, format=FormatScalarFormatter('%.0f'))

    ax2 = fig.add_subplot(gs[6:, 0], sharex=ax1)
    kw_w2 = {'norm': mcolors.CenteredNorm(),
           'levels': np.linspace(-6e-4,6e-4,31),
           'extend': 'both',
           'cmap': 'RdBu_r'
          }
    kw_b2 = {'levels': np.arange(-1,17,0.5)*1e-5, #np.arange(-1,2.36,0.07)*1e-5, 
           'colors': 'k', 
           'linewidths': 1
          }
    Cw2 = ax2.contourf(ds.xC, ds.zF, ds.W.isel(time=itime), **kw_w2)
    Cb2 = ax2.contour(ds.xC, ds.zC, ds.Bt.isel(time=itime), **kw_b2)
    cbar2 = fig.colorbar(Cw2, ax=ax2, fraction=0.5, format=FormatScalarFormatter('%.0f'))
    
    ani = animation.FuncAnimation(fig=fig, func=update, init_func=init, frames=np.arange(1, ds.dims['time']),
                                  interval=400, repeat_delay=800, blit=True)
    
    # save animation
    ani.save(figs_dir+args.fname_out, writer='ffmpeg', fps=3, dpi=200)

    def init():
        cbar1.set_label(r'w|$_{z=-10 \:m}$ [m $s^{-1}$]', labelpad=-50)
        cbar1.formatter.set_powerlimits((0, 0))
        cbar1.formatter.set_useMathText(True)
        cbar1.ax.yaxis.set_offset_position('left')
        cbar1.ax.set_yticks([-3e-3, -2e-3, -1e-3, 0, 1e-3, 2e-3, 3e-3])
        ax1.set_title(f't = {(itime*dhr)//24} days, {(itime*dhr)%24} hours')
        ax1.set(aspect=1)
        ax1.set_ylabel('Y [m]')

        cbar2.set_label(r'$\langle w\rangle_{y}$ [m $s^{-1}$]', labelpad=-50)
        cbar2.formatter.set_powerlimits((0, 0))
        cbar2.formatter.set_useMathText(True)
        cbar2.ax.yaxis.set_offset_position('left')
        cbar2.ax.set_yticks([-6e-4, -4e-4, -2e-4, 0, 2e-4, 4e-4, 6e-4])
        ax2.set_ylim(-100,0)
        ax2.set_ylabel('Z [m]')
        ax2.set_xlabel('X [m]')
        return Cw1, Cb1, Cw2, Cb2

    def update(frame):
        global Cw1, Cb1, Cw2, Cb2
        # for each frame, get new data and clear old data stored on each artist
        wtop = ds.w.isel(time=frame)
        btop = ds.bt.isel(time=frame)
        w_ymean = ds.W.isel(time=frame)
        b_ymean = ds.Bt.isel(time=frame)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for obj in (Cw1, Cb1, Cw2, Cb2):
                for coll in obj.collections:
                    coll.remove()
        # update the plot
        Cw1 = ax1.contourf(ds.xC, ds.yC, wtop, **kw_w1)
        Cb1 = ax1.contour(ds.xC, ds.yC, btop, **kw_b1)
        Cw2 = ax2.contourf(ds.xC, ds.zF, w_ymean, **kw_w2)
        Cb2 = ax2.contour(ds.xC, ds.zC, b_ymean, **kw_b2)
        ax1.set_title(f't = {(frame*dhr)//24} days, {(frame*dhr)%24} hours');
        return Cw1, Cb1, Cw2, Cb2

if __name__ == "__main__":
    main()