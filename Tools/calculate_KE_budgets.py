#!/usr/bin/env python3

import os
import sys
import time
import warnings
import argparse
import numpy as np
import xarray as xr
from xgcm import Grid
from scipy import interpolate


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
    Vbak = ds.Vbak.data
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
    f = interpolate.interp1d(b, z, assume_sorted=False)
    if criteria == '5-percent':
        rtau = 0.05
        delb = b - rtau
        if np.all(delb >= 0) or np.all(np.isnan(b)) or (delb[-1] <= 0):
            mld = np.nan
        else:
            last_idx = np.max(np.where(delb <= 0))
            crossing = range(last_idx, (last_idx + 2))
            f = interpolate.interp1d(b[crossing], z[crossing], assume_sorted=True)
            mld = -f(rtau)
    elif criteria == 'critical':
        if max(b) < 1e-8:
            mld = -z[-1]
        else:
            mld = -f(1e-8)
    elif criteria == 'max':
        idx_mld = np.argmax(b)
        mld = -z[idx_mld]
    else:
        mld = -f(b[-1] - criteria)
    return mld


def get_mld(ds, cvar='b', dims=['zC']):
    a = ds[cvar]
    if cvar == 'b':
        criteria = 0.001/ds.attrs['ρ₀']*9.81 #ds.attrs['N₀²'] * ds.attrs['hᵢ']
    elif cvar == 'dbdz':
        criteria = 'max'
    elif cvar == 'wuvn':
        criteria = '5-percent'
    elif cvar == 'TKE_eps':
        criteria = 'critical'
    return xr.apply_ufunc(get_mld_ufunc, a, ds.zC,
                          input_core_dims=[dims, dims],
                          output_core_dims=[[]],
                          output_dtypes=[float],
                          kwargs=dict(criteria = criteria),
                          vectorize=True)


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
    dsa = xr.open_dataset(data_dir+args.cname+'_averages.nc').isel(time=slice(1,None))#.drop_vars(['uym','vym','wym','bym','cym'])
    dsa.close()

    # fix initial mean buoyancy profile (only problematic for time averaged outputs)
    if dsa.b.isel(time=0).mean() == 0:
        b_ini = dsa.attrs['N₁²']*(dsa.zC + dsa.Lz) + \
               (dsa.attrs['N₀²'] - dsa.attrs['N₁²'])*np.maximum(dsa.zC + dsa.attrs['hᵢ'], 0)
        dsa['b'] = dsa.b.where(dsa.time != dsa.time[0], b_ini)

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
        b_f = grid.interp(dsa.b,    axis='z', boundary='extend')
        u_f = grid.interp(dsa.u,    axis='z', boundary='extend')
        v_f = grid.interp(dsa.v,    axis='z', boundary='extend')
        Vgf = grid.interp(dsa.Vbak, axis='z', boundary='extend')
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

    # fix Ertel PV at the surface, where the model thinks dVdz = 0
    qsurf = dsa.q.isel(zF=-1)
    qbak  = - dsa.attrs['M²']**2 / dsa.f
    dsa['q'] = dsa.q.where(dsa.zF != dsa.zF[-1], qsurf - qbak)

    # bulk Richardson number, boundary layer depth, and mixed layer depth
    dsa['dbdz']   = dsa.dbdz.where(dsa.zC != dsa.zC[0], dsa.attrs['N₁²'])
    dsa['timeTf'] = dsa.time/np.timedelta64(int(np.around(2*np.pi/dsa.f)), 's')
    dsa['Rib'], dsa['hRib'] = get_Rib_bld(dsa)
    #dsa['mld']    = get_mld(dsa)
    dsa['bld']    = get_mld(dsa, cvar='TKE_eps')
    dsa['hNsq']   = get_mld(dsa, cvar='dbdz')
    dsa['htau']   = get_mld(dsa, cvar='wuvn')
    dsa['htau']   = dsa.htau.where(dsa.timeTf > 0.05)
    dsa['sigmaC'] = dsa.zC / dsa.bld

    # integrate GSP and dissipation within the boundary layer
    dsa['iGSP'] = (dsa.GSP.where(dsa.sigmaC >= -1)*dzF).sum('zC')
    dsa['iEPS'] = (dsa.TKE_eps.where(dsa.sigmaC >= -1)*dzF).sum('zC')

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
