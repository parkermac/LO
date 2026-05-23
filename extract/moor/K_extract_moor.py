"""
S3-backed mooring extractor. Mirrors extract_moor.py but reads ROMS history/
average files from kopah S3 instead of local disk, replacing the ncks/ncrcat
subprocess pipeline with xarray + h5netcdf + fsspec + ThreadPoolExecutor.

Requires AWS credentials in the environment:
    AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY

Test on mac in ipython:
run K_extract_moor -gtx cas7_t0_x4b -test True

The same job would be run with flags as:
python K_extract_moor.py -gtx cas7_t2_x11b -ro 0 -0 2025.07.04 -1 2025.07.05 -lt hourly -sn test -lon ' -125' -lat 47 -get_all True -bucket liveocean-pmacc > test.log &
run K_extract_moor.py -gtx cas7_t2_x11b -ro 0 -0 2025.07.04 -1 2025.07.05 -lt hourly -sn test -lon ' -125' -lat 47 -get_all True -bucket liveocean-pmacc
run K_extract_moor.py -gtx cas7_t1_x11ab -ro 0 -0 2025.07.04 -1 2025.07.05 -lt average -sn test -lon ' -125' -lat 47 -get_all True -bucket liveocean-pmacc

NOTE: the quotes and space are required to feed it a negative longitude.
"""

import sys
import os
import argparse
from concurrent.futures import ThreadPoolExecutor
from time import time

import numpy as np
import pandas as pd
import xarray as xr
import s3fs

from lo_tools import Lfun, zrfun, zfun

if __name__ == '__main__':

    # ── command line arguments ────────────────────────────────────────────────────
    parser = argparse.ArgumentParser()
    # which run to use
    parser.add_argument('-gtx', '--gtagex', type=str)
    parser.add_argument('-ro', '--roms_out_num', type=int)
    # select time period and frequency
    parser.add_argument('-0', '--ds0', type=str)
    parser.add_argument('-1', '--ds1', type=str)
    parser.add_argument('-lt', '--list_type', type=str)
    # select mooring location and name
    parser.add_argument('-lon', type=float)
    parser.add_argument('-lat', type=float)
    parser.add_argument('-sn', type=str)
    # select categories of variables to extract
    parser.add_argument('-get_tsa',      type=zfun.boolean_string, default=False)
    parser.add_argument('-get_vel',      type=zfun.boolean_string, default=False)
    parser.add_argument('-get_bio',      type=zfun.boolean_string, default=False)
    parser.add_argument('-get_surfbot',  type=zfun.boolean_string, default=False)
    parser.add_argument('-get_pressure', type=zfun.boolean_string, default=False)
    parser.add_argument('-get_all',      type=zfun.boolean_string, default=False)
    # S3 options
    parser.add_argument('-bucket', type=str, default='liveocean')
    # number of parallel threads for S3 reads
    parser.add_argument('-Nproc', type=int, default=10)
    # optional testing
    parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)
    args = parser.parse_args()

    argsd = args.__dict__
    for a in ['gtagex']:
        if argsd[a] is None:
            print('*** Missing required argument: ' + a)
            sys.exit()

    gridname, tag, ex_name = args.gtagex.split('_')
    Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
    for a in argsd.keys():
        if a not in Ldir.keys():
            Ldir[a] = argsd[a]

    if Ldir['testing']:
        Ldir['roms_out_num'] = 0
        Ldir['ds0'] = '2017.07.04'
        Ldir['ds1'] = '2017.07.06'
        Ldir['list_type'] = 'hourly'
        Ldir['lon'] = -125
        Ldir['lat'] = 47
        Ldir['sn'] = 'test'
        Ldir['get_all'] = True

    if Ldir['get_all']:
        Ldir['get_tsa']     = True
        Ldir['get_vel']     = True
        Ldir['get_bio']     = True
        Ldir['get_surfbot'] = True

    tt00 = time()

    # ── output location ───────────────────────────────────────────────────────────
    out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'moor'
    Lfun.make_dir(out_dir)
    moor_fn = out_dir / (Ldir['sn'] + '_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '.nc')
    moor_fn.unlink(missing_ok=True)
    print(moor_fn)

    # ── S3 setup ──────────────────────────────────────────────────────────────────
    ENDPOINT = 'https://s3.kopah.uw.edu'
    fs = s3fs.S3FileSystem(
        key=os.environ['AWS_ACCESS_KEY_ID'],
        secret=os.environ['AWS_SECRET_ACCESS_KEY'],
        client_kwargs={'endpoint_url': ENDPOINT},
    )

    BUCKET = Ldir['bucket']
    GTAGEX = Ldir['gtagex']

    # ── build S3 file list ────────────────────────────────────────────────────────
    dates = pd.date_range(
        Ldir['ds0'].replace('.', '-'),
        Ldir['ds1'].replace('.', '-'),
        freq='D',
    )

    if Ldir['list_type'] == 'hourly':
        # ocean_his_0001 repeats the last snapshot of the previous day; skip it.
        fn_list = [
            f's3://{BUCKET}/LO_roms/{GTAGEX}/f{d.strftime("%Y.%m.%d")}/ocean_his_{n:04d}.nc'
            for d in dates
            for n in range(2, 26)
        ]
    elif Ldir['list_type'] in ('daily', 'average', 'lowpass'):
        fn_list = [
            f's3://{BUCKET}/LO_roms/{GTAGEX}/f{d.strftime("%Y.%m.%d")}/ocean_avg_0001.nc'
            for d in dates
        ]
    else:
        print('*** Unsupported list_type: ' + Ldir['list_type'])
        sys.exit()

    print('Files to open: %d' % len(fn_list))

    # ── read grid info and sigma params from first S3 file ────────────────────────
    # Replaces zrfun.get_basic_info(), which expects a local path.
    # block_size=2**22 (4 MB) fetches HDF5 metadata in one range request instead
    # of the ~dozens of small round trips h5py needs by default over S3.
    with fs.open(fn_list[0], 'rb', cache_type='blockcache', block_size=2**22) as f:
        ds0 = xr.open_dataset(f, engine='h5netcdf')

        G = {
            'lon_rho':  ds0['lon_rho'].values,
            'lat_rho':  ds0['lat_rho'].values,
            'mask_rho': ds0['mask_rho'].values,
            'mask_u':   ds0['mask_u'].values,
            'mask_v':   ds0['mask_v'].values,
        }
        S = {
            's_rho':       ds0['s_rho'].values,
            's_w':         ds0['s_w'].values,
            'N':           ds0.sizes['s_rho'],
            'Cs_r':        ds0['Cs_r'].values,
            'Cs_w':        ds0['Cs_w'].values,
            'hc':          float(ds0['hc'].values),
            'Vtransform':  int(ds0['Vtransform'].values),
            'Vstretching': int(ds0['Vstretching'].values),
        }

        if 'NH4' in ds0.data_vars:
            bio_list = ['NO3','NH4','phytoplankton','zooplankton',
                        'SdetritusN','LdetritusN','SdetritusC','LdetritusC',
                        'oxygen','alkalinity','TIC','rho']
        else:
            bio_list = ['NO3','phytoplankton','zooplankton','detritus',
                        'Ldetritus','oxygen','alkalinity','TIC','rho']

        do_wetdry = 'wetdry_mask_rho' in ds0.data_vars
        all_vars  = list(ds0.data_vars)
        ds0.close()

    # ── build variable list ───────────────────────────────────────────────────────
    vn_list = ['h', 'zeta']
    if do_wetdry:
        vn_list += ['wetdry_mask_rho']
    if Ldir['get_tsa']:
        vn_list += ['salt', 'temp', 'AKs', 'AKv']
    if Ldir['get_vel']:
        vn_list += ['u', 'v', 'w', 'ubar', 'vbar']
        if do_wetdry:
            vn_list += ['wetdry_mask_u', 'wetdry_mask_v']
    if Ldir['get_bio']:
        vn_list += bio_list
    if Ldir['get_surfbot']:
        vn_list += ['Pair','Uwind','Vwind','shflux','ssflux','latent','sensible',
                    'lwrad','swrad','sustr','svstr','bustr','bvstr']
    if Ldir['get_pressure']:
        vn_list += ['salt','temp','u','v','Pair','Uwind','Vwind']

    # deduplicate, preserving order; drop variables absent from the model output
    seen = set()
    vn_list = [v for v in vn_list if v in all_vars and not (v in seen or seen.add(v))]

    # ── find nearest wet grid cell (identical logic to extract_moor.py) ───────────
    lon = Ldir['lon']
    lat = Ldir['lat']
    Lon = G['lon_rho'][0, :]
    Lat = G['lat_rho'][:, 0]

    if (lon < Lon[0]) or (lon > Lon[-1]):
        print('ERROR: lon out of bounds ' + moor_fn.name); sys.exit()
    if (lat < Lat[0]) or (lat > Lat[-1]):
        print('ERROR: lat out of bounds ' + moor_fn.name); sys.exit()

    ilon = zfun.find_nearest_ind(Lon, lon)
    ilat = zfun.find_nearest_ind(Lat, lat)

    def find_good(ilat, ilon, mask):
        if mask == 'rho':
            jig_list = [[0,0],[1,0],[-1,0],[0,1],[0,-1]]
        elif mask == 'u':
            jig_list = [[0,0],[0,-1]]
        elif mask == 'v':
            jig_list = [[0,0],[-1,0]]
        for jig in jig_list:
            Ilat = ilat + jig[0]
            Ilon = ilon + jig[1]
            if G['mask_' + mask][Ilat, Ilon] == 1:
                print('   - %s: (%d, %d) => (%d, %d), jig = [%d, %d]' %
                      (mask, ilat, ilon, Ilat, Ilon, jig[0], jig[1]))
                return Ilat, Ilon
        print('ERROR: no good nearby point found on mask for ' + mask)
        sys.exit()

    ilat_rho, ilon_rho = find_good(ilat,     ilon,     'rho')
    ilat_u,   ilon_u   = find_good(ilat_rho, ilon_rho, 'u')
    ilat_v,   ilon_v   = find_good(ilat_rho, ilon_rho, 'v')

    # ── per-file fetch: open, slice to mooring point, load into memory ────────────
    # ThreadPoolExecutor (not Dask) is used because the work is I/O-bound: each
    # thread independently opens one S3 file, extracts the single grid point, and
    # returns a tiny in-memory dataset. block_size=2**22 (4 MB) front-loads the
    # HDF5 metadata in one HTTP range request instead of dozens of small ones.
    def fetch_one(fn):
        try:
            with fs.open(fn, 'rb', cache_type='blockcache', block_size=2**22) as f:
                with xr.open_dataset(f, engine='h5netcdf') as ds:
                    ds = ds[[v for v in vn_list if v in ds.data_vars]]
                    isel_kw = {}
                    if 'eta_rho' in ds.dims: isel_kw['eta_rho'] = ilat_rho
                    if 'xi_rho'  in ds.dims: isel_kw['xi_rho']  = ilon_rho
                    if 'eta_u'   in ds.dims: isel_kw['eta_u']   = ilat_u
                    if 'xi_u'    in ds.dims: isel_kw['xi_u']    = ilon_u
                    if 'eta_v'   in ds.dims: isel_kw['eta_v']   = ilat_v
                    if 'xi_v'    in ds.dims: isel_kw['xi_v']    = ilon_v
                    return ds.isel(**isel_kw).load()
        except Exception as e:
            print('WARNING: skipping %s: %s' % (fn, e))
            return None

    tt_open = time()
    with ThreadPoolExecutor(max_workers=Ldir['Nproc']) as pool:
        raw = list(pool.map(fetch_one, fn_list))
    n_skip = sum(1 for x in raw if x is None)
    if n_skip:
        print('WARNING: %d/%d files skipped (missing or corrupt)' % (n_skip, len(fn_list)))
    slices = [x for x in raw if x is not None]
    if not slices:
        print('ERROR: no files successfully read'); sys.exit()
    print(' - fetch all files: %0.2f sec' % (time() - tt_open))

    # ── concatenate and write ─────────────────────────────────────────────────────
    # Hierarchical concat: first join within each day-sized chunk, then join
    # chunks. For 8000+ files this is much faster than one flat xr.concat because
    # xarray's coordinate alignment cost is O(N) in the number of datasets.
    tt_write = time()
    CHUNK = 100
    if len(slices) > CHUNK:
        chunks = [xr.concat(slices[i:i+CHUNK], dim='ocean_time')
                  for i in range(0, len(slices), CHUNK)]
        ds = xr.concat(chunks, dim='ocean_time')
    else:
        ds = xr.concat(slices, dim='ocean_time')
    ds = ds.squeeze()
    ds.to_netcdf(moor_fn)
    ds.close()
    print(' - concat + write: %0.2f sec' % (time() - tt_write))

    # ── add z coordinates (mirrors extract_moor.py exactly) ──────────────────────
    ds   = xr.load_dataset(moor_fn)
    zeta = ds.zeta.values
    NT   = len(zeta)
    hh   = ds.h.values * np.ones(NT)
    z_rho, z_w = zrfun.get_z(hh, zeta, S)
    ds['z_rho'] = (('ocean_time', 's_rho'), np.transpose(z_rho.data))
    ds['z_w']   = (('ocean_time', 's_w'),   np.transpose(z_w.data))
    ds.z_rho.attrs['units']     = 'm'
    ds.z_w.attrs['units']       = 'm'
    ds.z_rho.attrs['long name'] = 'vertical position on s_rho grid, positive up'
    ds.z_w.attrs['long name']   = 'vertical position on s_w grid, positive up'
    if 'salt' in ds.data_vars:
        ds.salt.attrs['units'] = 'g kg-1'
    ds.ocean_time.attrs['long_name'] = 'Time [UTC]'
    ds.attrs['format'] = 'netCDF-4'
    ds.to_netcdf(moor_fn)
    ds.close()

    # ── copy attributes from source file ─────────────────────────────────────────
    ds = xr.load_dataset(moor_fn)  # load into memory so we can overwrite the file
    with fs.open(fn_list[0], 'rb', cache_type='blockcache', block_size=2**22) as f:
        ds_src = xr.open_dataset(f, engine='h5netcdf')
        for vn in list(ds_src.data_vars):
            if vn not in ds.data_vars:
                continue
            for attr in ('long_name', 'units'):
                try:
                    ds[vn].attrs[attr] = ds_src[vn].attrs[attr]
                except KeyError:
                    pass
        ds_src.close()
    ds.to_netcdf(moor_fn)
    ds.close()

    print('- total elapsed time was %0.2f sec' % (time() - tt00))
    print('Path to file:\n%s' % str(moor_fn))
