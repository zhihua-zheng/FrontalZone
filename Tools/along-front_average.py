#!/usr/bin/env python3

import sys
import warnings
import argparse
import xarray as xr
import numpy as np
from xgcm import Grid
from dask.diagnostics import ProgressBar


def main():
    # process input arguments
    parser = argparse.ArgumentParser(description="""
            Average coarse grained fields in aloing-front direction.""")
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
        data_dir = '/glade/derecho/scratch/zhihuaz/FrontalZone/Output/'
    elif sys.platform == 'darwin':
        data_dir = '/Users/zhihua/Documents/Work/Research/Projects/TRACE-SEAS/FrontalZone/Data/'
    else:
        print('OS not supported.')

    # read data
    # dsa = xr.open_dataset(data_dir+args.cname+'_averages.nc').drop_vars(['us','vs','dusdz','dvsdz','Vbak','Bbak']).chunk('auto')
    # dsa.close()
    dsc = xr.open_dataset(data_dir+args.cname+'_coarse.nc')
    dsc.close()

    dsc_afm = dsc.mean('yC')
    varlist = ['u','w','b','c']
    dsc_sbmeso_prime = (dsc - dsc_afm)[varlist]
    dsc_sbmeso_covar = xr.Dataset()
    dsc_sbmeso_covar['wbs'] = dsc_sbmeso_prime.w * dsc_sbmeso_prime.b
    dsc_sbmeso_covar['ubs'] = dsc_sbmeso_prime.u * dsc_sbmeso_prime.b
    dsc_sbmeso_covar['wcs'] = dsc_sbmeso_prime.w * dsc_sbmeso_prime.c
    dsc_sbmeso_covar['ucs'] = dsc_sbmeso_prime.u * dsc_sbmeso_prime.c
    dsc_sbmeso_covar = dsc_sbmeso_covar.mean('yC')

    periodic_coords = {dim : dict(left=f'{dim}F', center=f'{dim}C') for dim in 'xyz'}
    bounded_coords = {dim : dict(outer=f'{dim}F', center=f'{dim}C') for dim in 'xyz'}
    coords = {dim : periodic_coords[dim] if tpl=='P' else bounded_coords[dim] for dim, tpl in zip('xyz', 'PPN')}
    grid = Grid(dsc, coords=coords, periodic=['x', 'y'])
    dxC = dsc.xC.diff('xC').data[0]
    dzF = dsc.zF.diff('zF').data

    # diag = xr.Dataset()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ba_fc = grid.interp(dsc_afm.b, axis='x')
        ca_fc = grid.interp(dsc_afm.c, axis='x')
        ba_cf = grid.interp(dsc_afm.b, axis='z', boundary='extend') # note boundary conditions
        ca_cf = grid.interp(dsc_afm.c, axis='z', boundary='extend')
        dsc_afm['dbadx'] = grid.diff(ba_fc, axis='x') / dxC
        dsc_afm['dcadx'] = grid.diff(ca_fc, axis='x') / dxC
        dsc_afm['dbadz'] = grid.diff(ba_cf, axis='z') / dzF[None,:,None]
        dsc_afm['dcadz'] = grid.diff(ca_cf, axis='z') / dzF[None,:,None]

    # diag['psiEb'] = -dsc_sbmeso_covar.wbs / (diag.dbadx - 3e-8)
    # # psiEc = -dsc_sbmeso_covar.wcs / dcadx
    # diag['Kb'] = -dsc_afm.wbt / diag.dbadz
    # diag['Kc'] = -dsc_afm.wct / diag.dcadz

    dsafm = xr.merge([dsc_afm, dsc_sbmeso_covar])
    fpath = data_dir+args.cname+'_afm.nc'
    dsafm.to_netcdf(fpath)


if __name__ == "__main__":
    main()
