#!/usr/bin/env python3

import sys
import warnings
import argparse
import xarray as xr
import numpy as np
# from datatree import DataTree
from pathlib import Path
from xgcm import Grid
from dask.diagnostics import ProgressBar
from spectrum import Gaussian_filter_2d

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
    with xr.open_dataset(data_dir+args.cname+'_averages.nc').chunk('auto') as doa:
        fpath_extra = data_dir+args.cname+'_averages_extra.nc'
        if Path(fpath_extra).is_file():
            with xr.open_dataset(fpath_extra).chunk('auto') as doa_extra:
                 doa = xr.merge([doa_extra, doa])

    # with xr.open_dataset(data_dir+args.cname+'_top_slice.nc').chunk('auto') as ds_top:
        # ds_top = ds_top.rename_vars(dict([(i, i+'_top') for i in list(ds_top.keys())]))
        # with xr.open_dataset(data_dir+args.cname+'_south_slice.nc').chunk('auto') as ds_south:
            # ds_south = ds_south.rename_vars(dict([(i, i+'_south') for i in list(ds_south.keys())]))
            # with xr.open_dataset(data_dir+args.cname+'_east_slice.nc').chunk('auto') as ds_east:
                # ds_east = ds_east.rename_vars(dict([(i, i+'_east') for i in list(ds_east.keys())]))
                # ds = xr.merge([ds_top, ds_south, ds_east, dsm])
                # dt = DataTree.from_dict({'slice/top': ds_top, 'slice/south': ds_south, 'slice/east': ds_east, 'average': dsm})

    dso = xr.open_dataset(data_dir+args.cname+'_full.nc').drop_vars('wfrc').chunk({'time': 4})
    dso.close()
    # # make the restart time as day 0
    # if args.cname == 'spinup':
    #     ds = ds.assign_coords(time=(ds.time - np.timedelta64(24*3+12,'h')))

    periodic_coords = { dim : dict(left=f'{dim}F', center=f'{dim}C') for dim in 'xyz' }
    bounded_coords = { dim : dict(outer=f'{dim}F', center=f'{dim}C') for dim in 'xyz' }
    coords = { dim : periodic_coords[dim] if tpl=='P' else bounded_coords[dim] for dim, tpl in zip('xyz', 'PPN') }
    grid = Grid(dso, coords=coords, periodic=['x','y'])
    dx = dso.xF.diff('xF').data[0]
    dy = dso.yF.diff('yF').data[0]
    Lz = dso.zF[-1] - dso.zF[0]
    H  = dso.zC[-1] - dso.zC[0]

    # bulk PV in original field
    doa['timeTf'] = doa.time/np.timedelta64(int(np.around(2*np.pi/doa.f)), 's')
    doa['bhm']    = dso.b.mean(['xC','yC']).interp(zF=doa.zC).drop_vars('zF')
    doa['PVvm_z'] = (doa.PVfz.isel(zC=-1) - doa.PVfz.isel(zC=0))/H
    doa['PVvm_x'] = (doa.RVx.isel(xF=0) + doa.attrs['M²']/doa.f)*(-doa.attrs['M²'])
    doa['PVvm_f'] = doa.f*(doa.bhm.isel(zC=-1) - doa.bhm.isel(zC=0))/H
    doa['PVvm']   = doa.PVvm_z + doa.PVvm_x + doa.PVvm_f

    # filter original field
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dso['ub'] = grid.interp(dso.u, axis='z') * grid.interp(dso.b, axis=['x','z'])
        dso['vb'] = grid.interp(dso.v, axis='z') * grid.interp(dso.b, axis=['y','z'])
        dso['wb'] = dso.w * dso.b
        dso['uu'] = grid.interp(dso.u, axis='z')**2
        dso['vu'] = grid.interp(dso.v, axis='z') * grid.interp(dso.u, axis=['x','y','z'])
        dso['wu'] = dso.w * grid.interp(dso.u, axis='x')
        dso['uv'] = grid.interp(dso.u, axis='z') * grid.interp(dso.v, axis=['x','y','z'])
        dso['vv'] = grid.interp(dso.v, axis='z')**2
        dso['wv'] = dso.w * grid.interp(dso.v, axis='y')
        # dso['uw'] = grid.interp(dso.u, axis='z') * grid.interp(dso.w, axis=['x','z'])
        # dso['vw'] = grid.interp(dso.v, axis='z') * grid.interp(dso.w, axis=['y','z'])
        # dso['ww'] = dso.w**2
    dsf = dso.map(Gaussian_filter_2d, args=(60, 'xy'), keep_attrs=True)

    dfa = xr.full_like(doa, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dz = grid.diff(dsf.zF, axis='z')
        v_cfc = grid.interp(dsf.v, axis='z')
        u_fcc = grid.interp(dsf.u, axis='z')
        v_fcf = grid.interp(dsf.v, axis=['x','y'])
        w_ffc = grid.interp(dsf.w, axis=['x','y','z'])
        bt    = grid.interp(dsf.b, axis='z') - dsf.attrs['M²'] * dsf.xC
        omegaz = grid.diff(v_cfc, axis='x')/dx - grid.diff(u_fcc, axis='y')/dy
        omegaz = grid.interp(omegaz, axis=['x','y'])
        omegax = grid.diff(w_ffc, axis='y')/dy - grid.diff(v_fcf, axis='z')/dz

    # bulk PV in filtered field
    dfa['PVfz'] = (omegaz*bt).mean(['xC','yC'])
    dfa['RVx']  = ((omegax*dz).sum('zC')/Lz).mean('yC')
    dfa['bhm']  = dsf.b.mean(['xC','yC']).interp(zF=dfa.zC).drop_vars('zF')
    dfa['PVvm_z'] = (dfa.PVfz.isel(zC=-1) - dfa.PVfz.isel(zC=0))/H
    dfa['PVvm_x'] = (dfa.RVx.isel(xF=0) + dfa.attrs['M²']/dfa.f)*(-dfa.attrs['M²'])
    dfa['PVvm_f'] = dfa.f*(dfa.bhm.isel(zC=-1) - dfa.bhm.isel(zC=0))/H
    dfa['PVvm']   = dfa.PVvm_z + dfa.PVvm_x + dfa.PVvm_f

    # PV flux from filtered field
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        bfdia = dsf.bdia - grid.diff(dsf.ub - grid.interp(dsf.u, axis='z') * grid.interp(dsf.b, axis=['x','z']), axis='x')/dx \
                         - grid.diff(dsf.vb - grid.interp(dsf.v, axis='z') * grid.interp(dsf.b, axis=['y','z']), axis='y')/dy \
                         - grid.diff(dsf.wb - dsf.w * dsf.b, axis='z')/dz
        uffrc = dsf.ufrc - grid.diff(dsf.uu - grid.interp(dsf.u, axis='z')**2, axis='x')/dx \
                         - grid.diff(dsf.vu - grid.interp(dsf.v, axis='z') * grid.interp(dsf.u, axis=['x','y','z']), axis='y')/dy \
                         - grid.diff(dsf.wu - dsf.w * grid.interp(dsf.u, axis='x'), axis='z')/dz
        vffrc = dsf.vfrc - grid.diff(dsf.uv - grid.interp(dsf.u, axis='z') * grid.interp(dsf.v, axis=['x','y','z']), axis='x')/dx \
                         - grid.diff(dsf.vv - grid.interp(dsf.v, axis='z')**2, axis='y')/dy \
                         - grid.diff(dsf.wv - dsf.w * grid.interp(dsf.v, axis='y'), axis='z')/dz
        dbfdy = grid.interp(grid.diff(grid.interp(dsf.b, axis='z'), axis='y')/dy, axis='y')
        dbfdx = grid.interp(grid.diff(grid.interp(dsf.b, axis='z'), axis='x')/dx, axis='x')

    dfs = xr.Dataset()
    dfs['JzD'] = - (omegaz + dsf.f) * bfdia
    dfs['JzF'] = - (uffrc*dbfdy - vffrc*dbfdx)
    dfs = dfs.isel(zC=-1)

    # add background fields
    # H = abs(ds_east.zF[0])
    # z_top_slice = ds_top.zC[0]
    # ds['bt'] = (-ds.attrs['M²'] * ds.xC + ds.b).transpose('time','yC','xC')
    # ds['vt'] = (-ds.attrs['M²'] / ds.f * (z_top_slice + H) + ds.v).transpose('time','yF','xC')
    # ds['Bt'] = (-ds.attrs['M²'] * ds.xC + ds.B).transpose('time','zC','xC')
    # ds['Vt'] = (-ds.attrs['M²'] / ds.f * (ds.zC + H) + ds.V).transpose('time','zC','xC')

    fpath = data_dir+args.cname+'_PVbulk.nc'
    delayed_oa = doa.to_netcdf(fpath, group='original',           compute=False)
    # delayed_so = dso.to_netcdf(fpath, group='original',         mode='a', compute=False)
    delayed_fa = dfa.to_netcdf(fpath, group='filtered', mode='a', compute=False)
    fpath_pv = data_dir+args.cname+'_PVflux.nc'
    delayed_fs = dfs.to_netcdf(fpath_pv, compute=False)
    with ProgressBar():
        delayed_oa.compute()
        # delayed_so.compute()
        delayed_fa.compute()
        delayed_fs.compute()

    # Path(data_dir+args.cname+'_averages.nc').unlink()
    # if Path(fpath_extra).is_file():
    #     Path(fpath_extra).unlink()
    # Path(data_dir+args.cname+'_top_slice.nc').unlink()
    # Path(data_dir+args.cname+'_south_slice.nc').unlink()
    # Path(data_dir+args.cname+'_east_slice.nc').unlink()


if __name__ == "__main__":
    main()