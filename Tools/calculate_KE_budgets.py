#!/usr/bin/env python3

import os
import sys
import time
import warnings
import argparse
import numpy as np
import xarray as xr
from xgcm import Grid


def main():
    # process input arguments
    parser = argparse.ArgumentParser(description="""
            Calculate kinetic energy budgets.""")
    parser.add_argument('-c', '--case', action='store', dest='cname',
            metavar='CASENAME', help='simulation case name')
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

    t0 = time.time()

    # specify file path
    if sys.platform == 'linux' or sys.platform == 'linux2':
        data_dir = '/glade/derecho/scratch/zhihuaz/FrontalZone/Output/'
    elif sys.platform == 'darwin':
        data_dir = '/Users/zhihua/Documents/Work/Research/Projects/TRACE-SEAS/FrontalZone/Data/'
    else:
        print('OS not supported.')

    # read data
    dsa = xr.open_dataset(data_dir+args.cname+'_averages.nc')#.drop_vars(['uym','vym','wym','bym','cym'])
    dsa.close()

    # construct coordinates
    periodic_coords = {dim : dict(left=f'{dim}F', center=f'{dim}C') for dim in 'z'}
    bounded_coords = {dim : dict(outer=f'{dim}F', center=f'{dim}C') for dim in 'z'}
    coords = {dim : periodic_coords[dim] if tpl=='P' else bounded_coords[dim] for dim, tpl in zip('z', 'N')}
    grid = Grid(dsa, coords=coords)
    dzF  = dsa.zF.diff('zF').data

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        #dzC = grid.diff(dsa.zC, axis='z', boundary='extend')
        #dzC = dzC.where(dzC.zF != dzC.zF[-1]).fillna(dzF[-1])
        #dzC = dzC.where(dzC.zF != dzC.zF[0]).fillna(dzF[0]).data
        u_f = grid.interp(dsa.u,    axis='z', boundary='extend')
        v_f = grid.interp(dsa.v,    axis='z', boundary='extend')
        Vgf = grid.interp(dsa.Vbak, axis='z', boundary='extend')
        #dvdz_f = grid.diff(dsa.v, axis='z', boundary='extend') / dzC 
        #Tau23 = -dsa['νₑ']*dvdz_f
        #dTau23dz = grid.diff(Tau23, axis='z') / dzF
        dsa['dudz'] = grid.diff(u_f, axis='z') / dzF
        dsa['dvdz'] = grid.diff(v_f, axis='z') / dzF
        dsa['dVdz'] = grid.diff(Vgf, axis='z') / dzF
    dsa['GSP']  = -dsa.wvt * dsa.dVdz
    dsa['ASPy'] = -dsa.wvt * dsa.dvdz
    dsa['ASPx'] = -dsa.wut * dsa.dudz
    dsa['SSPy'] = -dsa.wvt * dsa.dvsdz
    dsa['SSPx'] = -dsa.wut * dsa.dusdz
    dsa['ASP']  =  dsa.ASPx + dsa.ASPy
    dsa['SSP']  =  dsa.SSPx + dsa.SSPy
    dsa['CKE']  =  dsa.v * dsa.Vbak
    dsa['MKE']  = (dsa.u**2 + dsa.v**2) / 2
    dsa['TKE']  = (dsa.uut + dsa.vvt + dsa.wwt) / 2
    dsa['TKE_vis'] = dsa.TKE_sgs - dsa.TKE_eps

    MKE_stress_top = (dsa.u.isel(zC=-1) * dsa.Qu + dsa.v.isel(zC=-1) * dsa.Qv) / dzF[-1]
    dsa['MKE_stress_top'] = xr.zeros_like(dsa.MKE_sgs)
    dsa['MKE_stress_top'] = dsa.MKE_stress_top.where(dsa.zC != dsa.zC[-1], MKE_stress_top) 
    dsa['MKE_stress'] = dsa.MKE_sgs + dsa.MKE_stress_top 

    CKE_stress_top = dsa.Vbak.isel(zC=-1) * dsa.Qv / dzF[-1]
    dsa['CKE_stress_top'] = xr.zeros_like(dsa.CKE_sgs)
    dsa['CKE_stress_top'] = dsa.CKE_stress_top.where(dsa.zC != dsa.zC[-1], CKE_stress_top)
    dsa['CKE_stress'] = dsa.CKE_sgs + dsa.CKE_stress_top
    dsa['CKE_GPW'] = -dsa.u * dsa.Vbak * dsa.f

    #dsa['wusgs'] = dsa.wusgs.interp(zF=dsa.zC).drop_vars('zF')
    #dsa['wvsgs'] = dsa.wvsgs.interp(zF=dsa.zC).drop_vars('zF')
    #dsa['wbsgs'] = dsa.wbsgs.interp(zF=dsa.zC).drop_vars('zF')
    #dsa['wcsgs'] = dsa.wcsgs.interp(zF=dsa.zC).drop_vars('zF')

    fpath = data_dir+args.cname+'_KE_budgets.nc'
    dsa.to_netcdf(fpath)
    t1 = time.time()
    print(f'Computation finished in {((t1-t0)/60):.1f} minutes')


if __name__ == "__main__":
    main()
