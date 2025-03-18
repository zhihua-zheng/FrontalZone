#!/usr/bin/env python3

import os
import sys
import time
import warnings
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from xgcm import Grid
from scipy import interpolate
from dask.distributed import Client
from dask_jobqueue import PBSCluster


def bld_from_Rib(Rib, d, Ribc=0.3):
    if np.nanmin(Rib) > Ribc:
        bld = d[-1]
    elif np.nanmax(Rib) < Ribc:
        bld = np.nan
    else:
        f = interpolate.interp1d(Rib, d, assume_sorted=False)
        bld = f(Ribc)
    return bld


def get_vt2(Ribc, d, dbdz, ustar, B0, h):
    wstar3 = B0*d
    if B0 == 0:
        zeta = np.zeros_like(d)
    else:
        LObukhov = -ustar**3/0.4/B0
        zeta = d/LObukhov if B0 < 0 else np.maximum(0.1*h/LObukhov, d/LObukhov)
    if B0 <= 0:
        phis = 1 + 5*zeta
    else:
        phis = (1 - 16*zeta)**(-1/2)
        idx_fc = zeta < -1
        phis[idx_fc] = (-28.86 - 98.96*zeta[idx_fc])**(-1/3)
    ws = 0.4*ustar/phis
    # ws[d>h] = 1e-16
    k_osbl = np.argmin(np.abs(d - h))
    dbdz_e = dbdz[k_osbl] if k_osbl==(len(d) - 1) else np.maximum(dbdz[k_osbl], dbdz[k_osbl+1])
    Ne  = np.sqrt(np.abs(dbdz_e)) * np.sign(dbdz_e)
    Cv  = 2.1 - 200*np.maximum(0, np.minimum(Ne, 0.002))
    rL  = 1 #(1+0.49*La_sl**(-2))
    if B0 <= 0: # Large et al. 1994
        vt2 = Cv*d*Ne*ws*np.sqrt(0.2/98.96/0.1) / (0.4**2) / Ribc
    else: # Li & Fox-Kemper 2017
        vt2 = Cv*d*Ne*np.sqrt((0.15*wstar3 + 0.17*ustar**3*rL) / ws) / Ribc
    return np.maximum(0, vt2)


def get_Rib_bld(ds, Ribc=0.3):
    bld_array = xr.full_like(ds.timeTf, fill_value=np.nan).assign_attrs(units='m', long_name='Boundary layer depth')
    Rib_array = xr.full_like(ds.b,      fill_value=np.nan).assign_attrs(units='',  long_name='Bulk Richardson number')
    d  = np.abs(ds.zC).data
    dz = np.diff(ds.zF)
    Vbak = ds.Vbak.isel(time=0).data
    for i in range(ds.sizes['time']):
        b = ds.b.isel(time=i).data
        u = ds.u.isel(time=i).data
        v = ds.v.isel(time=i).data
        ustar = np.sqrt(ds.ustar2.isel(time=i).data)
        B0    = ds.Qb.isel(time=i).data
        dbdz  = ds.dbdz.isel(time=i).data
        if i == 0:
            duv2 = ((u[-1] - u)**2 + (v[-1] - v + Vbak[-1] - Vbak)**2)
            duv2[duv2==0]  = np.nan
            Rib_array[i,:] = d*(b[-1] - b) / duv2
        else:
            in_sl = d <= previous_bld/10 if previous_bld/10 > d[-1] else d == d[-1]

            bsl   = np.average(b[in_sl], weights=dz[in_sl])
            usl   = np.average(u[in_sl], weights=dz[in_sl])
            vsl   = np.average(v[in_sl] + Vbak[in_sl], weights=dz[in_sl])
            duv2  = ((usl - u)**2 + (vsl - v - Vbak)**2)
            vt2   = get_vt2(Ribc, d, dbdz, ustar, B0, previous_bld)
            Duv2  = duv2 + vt2
            Duv2[Duv2==0]  = np.nan
            Rib_array[i,:] = d*(bsl - b) / Duv2
        previous_bld = bld_from_Rib(Rib_array[i,:], d)
        bld_array[i] = previous_bld
    return Rib_array, bld_array


def get_mld_ufunc(b, z, criteria=5.4e-6):
    if isinstance(criteria, str) and criteria.startswith('max_over'):
        threshold = float(criteria.split('_')[-1])
        idx_mld = arg_local_max_last(b, threshold=threshold)#np.argmax(b)
        mld = -z[idx_mld]
    else:
        if criteria == '5-percent':
            rtau = 0.05
            delb = b - rtau
        elif criteria == 'critical':
            delb = b - 1e-8
        elif criteria == 0:
            delb = b - 0
        else:
            delb = b - (b[-1] - criteria)

        if np.all(delb >= 0) or (delb[-1] <= 0) or np.all(np.isnan(b)):
            mld = np.nan
        else:
            last_idx = np.max(np.where(delb <= 0))
            crossing = range(last_idx, (last_idx + 2))
            f = interpolate.interp1d(delb[crossing], z[crossing], assume_sorted=True)
            mld = -f(0)
    return mld


def get_mld(ds, cvar='b', dims=['zC']):
    a = ds[cvar]
    if cvar == 'b':
        criteria = 0.03/ds.attrs['ρ₀']*9.81 #ds.attrs['N₀²'] * ds.attrs['hᵢ']
    elif cvar == 'dbdz':
        N2_base = ds.attrs['N₁²']
        criteria = f'max_over_{N2_base}'
    elif cvar == 'wuvn':
        criteria = '5-percent'
    elif cvar == 'TKE_eps' or cvar == 'eps':
        criteria = 'critical'
    elif cvar == 'wb':
        criteria = 0
    return xr.apply_ufunc(get_mld_ufunc, a, ds.zC,
                          input_core_dims=[dims, dims],
                          output_core_dims=[[]],
                          output_dtypes=[float],
                          kwargs=dict(criteria=criteria),
                          dask='parallelized',
                          vectorize=True)


def arg_local_max_last(y, threshold):
    """
    https://gist.github.com/ben741/d8c70b608d96d9f7ed231086b237ba6b
    """
    idx_peaks = np.where((y[1:-1] > y[0:-2]) * (y[1:-1] > y[2:]) * (y[1:-1] > threshold))[0] + 1
    return idx_peaks[-1]


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
    dsf = xr.open_dataset(data_dir+args.cname+'_full.nc').drop_vars(['u', 'v', 'w', 'q'])\
            .chunk({'xC':200, 'yC':200, 'time':1})
    dsf.close()
    dsf = dsf.where(((dsf.time / np.timedelta64(int(dsf.out_interval_slice), 's')) % 1) == 0, drop=True)
    dsf = dsf.drop_duplicates(dim='time', keep='last')

    dsa = xr.open_dataset(data_dir+args.cname+'_averages.nc')
    dsa.close()
    dsa = dsa.where(((dsa.time / np.timedelta64(int(dsa.out_interval_mean), 's')) % 1) == 0, drop=True).transpose('time',...)
    dsa = dsa.drop_duplicates(dim='time', keep='last')

    # construct coordinates
    periodic_coords = {dim : dict(left=f'{dim}F', center=f'{dim}C') for dim in 'z'}
    bounded_coords = {dim : dict(outer=f'{dim}F', center=f'{dim}C') for dim in 'z'}
    coords = {dim : periodic_coords[dim] if tpl=='P' else bounded_coords[dim] for dim, tpl in zip('z', 'N')}
    grid = Grid(dsa, coords=coords)
    dzF  = dsa.zF.diff('zF').data

    dsa['wusgs'] = dsa.wusgs.where(dsa.zF != dsa.zF[-1], dsa.Qu)
    dsa['wvsgs'] = dsa.wvsgs.where(dsa.zF != dsa.zF[-1], dsa.Qv)
    dsa['wbsgs'] = dsa.wbsgs.where(dsa.zF != dsa.zF[-1], dsa.Qb)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wusgs_c = grid.interp(dsa.wusgs, axis='z')
        wvsgs_c = grid.interp(dsa.wvsgs, axis='z')
        wbsgs_c = grid.interp(dsa.wbsgs, axis='z')
        #wcsgs_c = grid.interp(dsa.wcsgs, axis='z')
        wut_f = grid.interp(dsa.wut, axis='z', boundary='fill', fill_value=0)
        wvt_f = grid.interp(dsa.wvt, axis='z', boundary='fill', fill_value=0)
        wbt_f = grid.interp(dsa.wbt, axis='z', boundary='fill', fill_value=0)
        dsa['dwudz'] = grid.diff((wut_f + dsa.wusgs), axis='z') / dzF
        dsa['dwvdz'] = grid.diff((wvt_f + dsa.wvsgs), axis='z') / dzF
        dsa['dwbdz'] = grid.diff((wbt_f + dsa.wbsgs), axis='z') / dzF
        b3f = grid.interp(dsf.b,    axis='z', boundary='extend').transpose(...,'zF')
        b_f = grid.interp(dsa.b,    axis='z', boundary='extend')
        u_f = grid.interp(dsa.u,    axis='z', boundary='extend')
        v_f = grid.interp(dsa.v,    axis='z', boundary='extend')
        Vgf = grid.interp(dsa.Vbak, axis='z', boundary='extend')
        dsf['dbdz'] = grid.diff(b3f, axis='z') / dzF
        dsa['dbdz'] = grid.diff(b_f, axis='z') / dzF
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

    dsa['wu'] = dsa.wut + wusgs_c
    dsa['wv'] = dsa.wvt + wvsgs_c
    dsa['wb'] = dsa.wbt + wbsgs_c
    #dsa['wc'] = dsa.wct + wcsgs_c
    dsa['ustar2'] = np.sqrt(dsa.Qu**2 + dsa.Qv**2)
    dsa['wuv']    = np.sqrt(dsa.wu**2 + dsa.wv**2)
    dsa['wuvn']   = dsa.wuv / dsa.ustar2
    dsa['EBF']    = dsa.Qv * dsa.attrs['M²'] / dsa.f
    dsa['GSP_ful'] = -dsa.wv * dsa.dVdz
    dsa['ASP_ful'] = -dsa.wv * dsa.dvdz -dsa.wu * dsa.dudz

    # fix Ertel PV at the surface, where the model thinks dVdz = 0
    qsurf = dsa.q.isel(zF=-1)
    qbak  = - dsa.attrs['M²']**2 / dsa.f
    dsa['q'] = dsa.q.where(dsa.zF != dsa.zF[-1], qsurf - qbak)

    # bulk Richardson number, boundary layer depth, and mixed layer depth
    dsa['dbdz']   = dsa.dbdz.where(dsa.zC != dsa.zC[0], dsa.attrs['N₁²'])
    dsa['timeTf'] = dsa.time/np.timedelta64(int(np.around(2*np.pi/dsa.f)), 's')
    dsa['Rib'], dsa['hRib'] = get_Rib_bld(dsa)
    dsa['heps'] = get_mld(dsa, cvar='TKE_eps')
    dsa['bld']  = get_mld(dsa, cvar='dbdz')#.rolling(time=15, center=True, min_periods=1).median()
    dsa['mld']  = get_mld(dsa, cvar='b')#.rolling(time=15, center=True, min_periods=1).median()
    dsa['htau'] = get_mld(dsa, cvar='wuvn')
    dsa['htau'] = dsa.htau.where(dsa.timeTf > 0.05)

    USER = os.getenv('USER')
    TMPDIR = f'/glade/derecho/scratch/{USER}/temp'
    job_script_prologue = [f'export TMPDIR={TMPDIR}', 'mkdir -p $TMPDIR']
    cluster_kw = dict(job_name=args.cname+'_hm',
                      cores=1,
                      memory='4GiB',
                      processes=1,
                      local_directory=f'{TMPDIR}/pbs.$PBS_JOBID/dask/spill',
                      log_directory=f'{TMPDIR}/pbs.$PBS_JOBID/dask/worker_logs',
                      job_extra_directives=['-j oe', '-A UMCP0036'],
                      job_script_prologue=job_script_prologue,
                      resource_spec='select=1:ncpus=1:mem=4GB',
                      queue='casper',
                      walltime='30:00',
                      interface='ext')
    with PBSCluster(**cluster_kw) as cluster, Client(cluster) as client:
        cluster.scale(32)
        #dsa['bld'] = get_mld(dsf, cvar='dbdz').mean(['xC','yC']).compute()
        #dsa['mld'] = get_mld(dsf, cvar='b').mean(['xC','yC']).compute()
        dsft = dsf[['eps', 'dbdz']].where(dsf.eps >= 1e-10)
        meps = dsft.eps.mean(['xC','yC']).compute()
        mNsq = dsft.dbdz.mean(['xC','yC']).compute()

    dsa['LOz'] = 2*np.pi*np.sqrt(meps / mNsq**(3/2)) 
    dsa['sigmaC'] = dsa.zC / dsa.bld

    # integrate GSP and dissipation within the boundary layer
    dsa['iGSP'] = (dsa.GSP.where(dsa.sigmaC >= -1)*dzF).sum('zC')
    dsa['iASP'] = (dsa.ASP.where(dsa.sigmaC >= -1)*dzF).sum('zC')
    dsa['iEPS'] = (dsa.TKE_eps.where(dsa.sigmaC >= -1)*dzF).sum('zC')
    dsa['iGSP_ful'] = (dsa.GSP_ful.where(dsa.sigmaC >= -1)*dzF).sum('zC')
    dsa['iASP_ful'] = (dsa.ASP_ful.where(dsa.sigmaC >= -1)*dzF).sum('zC')

    # average statistics in the last inertial period
    dsa = dsa.assign_coords(sC=('sC', np.arange(-1.3, 0, 0.01)))
    time_interval = (dsa.timeTf >= 2) & (dsa.timeTf < 3)
    dsali = dsa.where(time_interval, drop=True)
    apuf_kwargs = dict(input_core_dims=[['sC'], ['zC'], ['zC']],
                       output_core_dims=[['sC']],
                       output_dtypes=[float],
                       vectorize=True,
                       dask='parallelized')
    dsa['u_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.u, **apuf_kwargs).mean('time')
    dsa['v_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.v, **apuf_kwargs).mean('time')
    dsa['b_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.b, **apuf_kwargs).mean('time')
    dsa['dwudz_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.dwudz, **apuf_kwargs).mean('time')
    dsa['dwvdz_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.dwvdz, **apuf_kwargs).mean('time')
    dsa['dwbdz_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.dwbdz, **apuf_kwargs).mean('time')
    dsa['GSP_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.GSP, **apuf_kwargs).mean('time')
    dsa['ASP_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.ASP, **apuf_kwargs).mean('time')
    dsa['wbt_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.wbt, **apuf_kwargs).mean('time')
    dsa['wut_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.wut, **apuf_kwargs).mean('time')
    dsa['wvt_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.wvt, **apuf_kwargs).mean('time')
    dsa['wb_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.wb, **apuf_kwargs).mean('time')
    dsa['wu_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.wu, **apuf_kwargs).mean('time')
    dsa['wv_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.wv, **apuf_kwargs).mean('time')
    dsa['GSP_ful_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.GSP_ful, **apuf_kwargs).mean('time')
    dsa['ASP_ful_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.ASP_ful, **apuf_kwargs).mean('time')
    dsa['TKE_tur_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.TKE_tur, **apuf_kwargs).mean('time')
    dsa['TKE_prs_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.TKE_prs, **apuf_kwargs).mean('time')
    dsa['TKE_vis_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.TKE_vis, **apuf_kwargs).mean('time')
    dsa['TKE_eps_sc'] = xr.apply_ufunc(np.interp, dsa.sC, dsali.sigmaC, dsali.TKE_eps, **apuf_kwargs).mean('time')

    fpath = data_dir+args.cname+'_KE_budgets.nc'
    dsa.to_netcdf(fpath)

    t1 = time.time()
    print(f'Computation finished in {((t1-t0)/60):.1f} minutes')


if __name__ == "__main__":
    main()
