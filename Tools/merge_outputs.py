#!/usr/bin/env python3

import sys
import argparse
import xarray as xr
import numpy as np
# from datatree import DataTree
from pathlib import Path
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
    with xr.open_dataset(data_dir+args.cname+'_averages.nc').chunk('auto') as dsa:
        fpath_extra = data_dir+args.cname+'_averages_extra.nc'
        if Path(fpath_extra).is_file():
            with xr.open_dataset(fpath_extra).chunk('auto') as ds_extra:
                 dsa = xr.merge([ds_extra, dsa])

    # with xr.open_dataset(data_dir+args.cname+'_top_slice.nc').chunk('auto') as ds_top:
        # ds_top = ds_top.rename_vars(dict([(i, i+'_top') for i in list(ds_top.keys())]))
        # with xr.open_dataset(data_dir+args.cname+'_south_slice.nc').chunk('auto') as ds_south:
            # ds_south = ds_south.rename_vars(dict([(i, i+'_south') for i in list(ds_south.keys())]))
            # with xr.open_dataset(data_dir+args.cname+'_east_slice.nc').chunk('auto') as ds_east:
                # ds_east = ds_east.rename_vars(dict([(i, i+'_east') for i in list(ds_east.keys())]))
                # ds = xr.merge([ds_top, ds_south, ds_east, dsm])
                # dt = DataTree.from_dict({'slice/top': ds_top, 'slice/south': ds_south, 'slice/east': ds_east, 'average': dsm})
    dst = xr.open_dataset(data_dir+args.cname+'_top_slice.nc').chunk('auto')
    dst.close()
    dss = xr.open_dataset(data_dir+args.cname+'_south_slice.nc').chunk('auto')
    dss.close()
    dse = xr.open_dataset(data_dir+args.cname+'_east_slice.nc').chunk('auto')
    dse.close()
    # # make the restart time as day 0
    # if args.cname == 'spinup':
    #     ds = ds.assign_coords(time=(ds.time - np.timedelta64(24*3+12,'h')))
    
    # add background fields
    # H = abs(ds_east.zF[0])
    # z_top_slice = ds_top.zC[0]
    # ds['bt'] = (-ds.attrs['M²'] * ds.xC + ds.b).transpose('time','yC','xC')
    # ds['vt'] = (-ds.attrs['M²'] / ds.f * (z_top_slice + H) + ds.v).transpose('time','yF','xC')
    # ds['Bt'] = (-ds.attrs['M²'] * ds.xC + ds.B).transpose('time','zC','xC')
    # ds['Vt'] = (-ds.attrs['M²'] / ds.f * (ds.zC + H) + ds.V).transpose('time','zC','xC')

    fpath = data_dir+args.cname+'.nc'
    delayed_obj_a = dsa.to_netcdf(fpath, group='average',     compute=False)
    delayed_obj_t = dst.to_netcdf(fpath, group='slice/top',   mode='a', compute=False)
    delayed_obj_s = dss.to_netcdf(fpath, group='slice/south', mode='a', compute=False)
    delayed_obj_e = dse.to_netcdf(fpath, group='slice/east',  mode='a', compute=False)
    with ProgressBar():
        delayed_obj_a.compute()
        delayed_obj_t.compute()
        delayed_obj_s.compute()
        delayed_obj_e.compute()


if __name__ == "__main__":
    main()