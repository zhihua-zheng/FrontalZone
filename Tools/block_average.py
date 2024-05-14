#!/usr/bin/env python3

import os
import sys
import argparse
import xarray as xr
from dask.distributed import Client
from dask_jobqueue import PBSCluster
from xgcm.autogenerate import generate_grid_ds


# for xr.apply_ufunc
def block_mean_ufunc(u, v, w, b, c, block_width=64):
    # print(f'u: {u.shape} | v: {v.shape} | w: {w.shape} | b: {b.shape} | c: {c.shape}')
    def _decompose(a, block_width):
        # Be careful with row-major!
        ny, nx = u.shape
        nyc, nxc = ny//block_width, nx//block_width
        facet = a.reshape(nyc, block_width, nxc, block_width)
        tilde = facet.mean(axis=(1, 3), keepdims=1)
        prime = facet - tilde
        return tilde.squeeze(), prime

    u_tilde, u_prime = _decompose(u, block_width)
    v_tilde, v_prime = _decompose(v, block_width)
    w_tilde, w_prime = _decompose(w, block_width)
    b_tilde, b_prime = _decompose(b, block_width)
    c_tilde, c_prime = _decompose(c, block_width)
    wbt = (w_prime*b_prime).mean(axis=(1, 3))
    wct = (w_prime*c_prime).mean(axis=(1, 3))
    wut = (w_prime*u_prime).mean(axis=(1, 3))
    wvt = (w_prime*v_prime).mean(axis=(1, 3))
    uvt = (u_prime*v_prime).mean(axis=(1, 3))
    ubt = (u_prime*b_prime).mean(axis=(1, 3))
    uct = (u_prime*c_prime).mean(axis=(1, 3))
    uut = (u_prime**2).mean(axis=(1, 3))
    vvt = (v_prime**2).mean(axis=(1, 3))
    wwt = (w_prime**2).mean(axis=(1, 3))
    bbt = (b_prime**2).mean(axis=(1, 3))
    cct = (c_prime**2).mean(axis=(1, 3))
    return u_tilde, v_tilde, w_tilde, b_tilde, c_tilde, \
        wbt, wct, wut, wvt, uvt, ubt, uct, uut, vvt, wwt, bbt, cct


def block_mean(ds, dims, block_width):
    nxc, nyc = ds[dims[0]].size//block_width, ds[dims[1]].size//block_width
    return xr.apply_ufunc(block_mean_ufunc, ds.u, ds.v, ds.w, ds.b, ds.c,
                          input_core_dims=[dims, dims, dims, dims, dims],
                          output_core_dims=[dims, dims, dims, dims, dims, dims,
                                            dims, dims, dims, dims, dims, dims,
                                            dims, dims, dims, dims, dims],
                          output_dtypes=[float, float, float, float, float,
                                         float, float, float, float, float,
                                         float, float, float, float, float,
                                         float, float],
                          exclude_dims=set(dims), vectorize=True,
                          dask='parallelized', kwargs={'block_width': block_width},
                          dask_gufunc_kwargs=dict(output_sizes={dims[0]: nxc, dims[1]: nyc}))


def main():
    # process input arguments
    parser = argparse.ArgumentParser(description="""
            Apply block averaging to 3D Oceananigans fields.""")
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
    dsf = xr.open_dataset(data_dir+args.cname+'_full.nc').drop_vars(['eps','PV','xF','yF']).chunk({'time':3, 'zC':16})
    dsf.close()

    # with dask.config.set(**{'array.slicing.split_large_chunks': False}):
    #     dsf_coarsened = dsf.coarsen(xC=64, yC=64) # correspond to 125 m, separation between 2D and 3D turbulence
    #     dsf_tilde = dsf_coarsened.mean()#.persist()
    #     dsf_facet = dsf_coarsened.construct(xC=('xC', 'xC_fine'), yC=('yC', 'yC_fine'))

    # dsf_prime = dsf_facet - dsf_tilde
    # dsf_variance = (dsf_prime**2).mean(['xC_fine', 'yC_fine']).rename_vars({'u':'uut', 'v':'vvt', 'w':'wwt', 'b':'bbt', 'c':'cct'})
    # dsf_covariance = xr.Dataset()
    # dsf_covariance['wbt'] = dsf_prime.w * dsf_prime.b
    # dsf_covariance['wct'] = dsf_prime.w * dsf_prime.c
    # dsf_covariance['wut'] = dsf_prime.w * dsf_prime.u
    # dsf_covariance['wvt'] = dsf_prime.w * dsf_prime.v
    # dsf_covariance['uvt'] = dsf_prime.u * dsf_prime.v
    # dsf_covariance['ubt'] = dsf_prime.u * dsf_prime.b
    # dsf_covariance['uct'] = dsf_prime.u * dsf_prime.c
    # dsf_covariance = dsf_covariance.mean(['xC_fine', 'yC_fine'])
    # dsc = xr.merge([dsf_tilde, dsf_variance, dsf_covariance]).drop_vars(['xF','yF'])

    USER = os.getenv('USER')
    TMPDIR = f'/glade/derecho/scratch/{USER}/temp'
    job_script_prologue = ['export TMPDIR=/glade/scratch/derecho/$USER/temp', 'mkdir -p $TMPDIR']
    cluster_kw = dict(job_name=args.cname+'_hmean',
                      cores=1,
                      memory='4GiB',
                      processes=1,
                      local_directory=f'{TMPDIR}/pbs.$PBS_JOBID/dask/spill',
                      log_directory=f'{TMPDIR}/pbs.$PBS_JOBID/dask/worker_logs',
                      job_extra_directives=['-j oe'],
                      job_script_prologue=job_script_prologue,
                      resource_spec='select=1:ncpus=1:mem=4GB',
                      queue='casper',
                      walltime='30:00',
                      interface='ext')
    # dask.config.config.get('distributed').get('dashboard').update({'link':'{JUPYTERHUB_SERVICE_PREFIX}/proxy/{port}/status'})
    with PBSCluster(**cluster_kw) as cluster, Client(cluster) as client:
        print(cluster.job_script())
        cluster.scale(8)
        print(cluster.dashboard_link.replace(':8787', ':1212/proxy/8787'))
        u, v, w, b, c, wbt, wct, wut, wvt, uvt, ubt, uct, uut, vvt, wwt, bbt, cct = block_mean(dsf, ['xC', 'yC'], 64)
        u.name = 'u'
        v.name = 'v'
        w.name = 'w'
        b.name = 'b'
        c.name = 'c'
        wbt.name = 'wbt'
        wct.name = 'wct'
        wut.name = 'wut'
        wvt.name = 'wvt'
        uvt.name = 'uvt'
        ubt.name = 'ubt'
        uct.name = 'uct'
        uut.name = 'uut'
        vvt.name = 'vvt'
        wwt.name = 'wwt'
        bbt.name = 'bbt'
        cct.name = 'cct'

        dsc = xr.merge([u, v, w, b, c, wbt, wct, wut, wvt, uvt, ubt, uct,
                        uut, vvt, wwt, bbt, cct])
        dsc['zF'] = dsf.zF
        dsc['xC'] = dsf.xC.coarsen(xC=64).mean()
        dsc['yC'] = dsf.yC.coarsen(yC=64).mean()
        dsc = generate_grid_ds(dsc, {'x':'xC', 'y':'yC'}).rename({'xC_left':'xF', 'yC_left':'yF'})
        fpath = data_dir+args.cname+'_coarse.nc'
        delayed_coarse = dsc.to_netcdf(fpath, compute=False)
        delayed_coarse.compute()
        client.close()


if __name__ == "__main__":
    main()
