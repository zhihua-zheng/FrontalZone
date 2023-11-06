#!/usr/bin/env python3

import sys
import argparse
import xarray as xr
import numpy as np
from dask.diagnostics import ProgressBar

def main():
    # process input arguments
    parser = argparse.ArgumentParser(description="""
            Merge several outputs from Oceananigans simualtion into one data file.""")
    parser.add_argument('-c', '--case', action='store', dest='cname',
            metavar='CASENAME', help='Simulation case name')
    # parser.add_argument('-o', '--output', action='store', dest='fname_out',
    #         metavar='FIGNAME', help='Output figure name')
    parser.add_argument('--version', action='version', version='%(prog)s: 1.0')
    # parsing arguments and save to args
    args = parser.parse_args()

    # check input
    if not args.cname:
        print('Oceananigans simulation case name are required. Stop.\n')
        parser.print_help()
        sys.exit(1)
    
    # specify file path
    if sys.platform == 'linux' or sys.platform == 'linux2':
        data_dir = '/glade/work/zhihuaz/Data/FrontalZone/'
    elif sys.platform == 'darwin':
        data_dir = '/Users/zhihua/Documents/Work/Research/Projects/TRACE-SEAS/FrontalZone/Data/'
    else:
        print('OS not supported.')

    # read data
    with xr.open_dataset(data_dir+args.cname+'_top_slice.nc').chunk('auto').isel(zC=0,zF=0) as ds_top:
        with xr.open_dataset(data_dir+args.cname+'_averages.nc').chunk('auto') as ds_average:
             ds = xr.merge([ds_average, ds_top])
    # ds_east = xr.open_dataset(data_dir+args.cname+'_east_slice.nc').load()
    # ds_east.close()
    # ds_south = xr.open_dataset(data_dir+args.cname+'_south_slice.nc').load()
    # ds_south.close()

    # # make the restart time as day 0
    # if args.cname == 'spinup':
    #     ds = ds.assign_coords(time=(ds.time - np.timedelta64(24*3+12,'h')))
    
    # add background fields
    H = abs(ds.zF[0])
    ds['bt'] = (-ds.attrs['M²'] * ds.xC + ds.b).transpose('time','yC','xC')
    ds['vt'] = (-ds.attrs['M²']/ds.f * (ds.zC.isel(zC=-6).data + H) + ds.v).transpose('time','yF','xC')
    ds['Bt'] = (-ds.attrs['M²'] * ds.xC + ds.B).transpose('time','zC','xC')
    ds['Vt'] = (-ds.attrs['M²']/ds.f * (ds.zC + H) + ds.V).transpose('time','zC','xC')
    
    delayed_obj = ds.to_netcdf(data_dir+args.cname+'.nc', compute=False)
    with ProgressBar():
        results = delayed_obj.compute()


if __name__ == "__main__":
    main()