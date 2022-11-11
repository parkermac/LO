"""
This creates a layers file for a single time.

To test on mac:

run make_layers -test True
"""

from pathlib import Path
import argparse
import xarray as xr
from time import time
from PyCO2SYS import CO2SYS
import seawater as sw
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages from CO2YS
import numpy as np
import sys

from lo_tools import plotting_functions as pfun

# defaults for testing
from lo_tools import Lfun
Ldir = Lfun.Lstart()
in_fn = (Ldir['roms_out'] /
    'cas6_v0_live' / 'f2019.07.04' / 'ocean_his_0002.nc')
out_fn = (Ldir['parent'] / 'LO_output' / 'post' /
    'cas6_v0_live' / 'f2019.07.04' / 'layers0' / 'tempfiles' / 'layers_000000.nc')
# in_fn = (Ldir['roms_out'] /
#     'cas6_v00_uu0mb' / 'f2021.07.04' / 'ocean_his_0002.nc')
# out_fn = (Ldir['parent'] / 'LO_output' / 'post' /
#     'cas6_v00_uu0mb' / 'f2021.07.04' / 'layers0' / 'tempfiles' / 'layers_000000.nc')
# make sure output directory exists
Lfun.make_dir(out_fn.parent)

parser = argparse.ArgumentParser()
parser.add_argument('-in_fn', default=in_fn, type=str)
parser.add_argument('-out_fn', default=out_fn, type=str)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)

args = parser.parse_args()
in_fn = Path(args.in_fn)
out_fn = Path(args.out_fn)
testing = args.testing

out_fn.unlink(missing_ok=True)

in_ds = xr.open_dataset(in_fn, decode_times=False)
out_ds = xr.Dataset()

# add time
out_ds['ocean_time'] = (('ocean_time'), in_ds.ocean_time.values)
out_ds['ocean_time'].attrs['units'] = Ldir['roms_time_units']
# add coordinates and other static 2-D fields
vn_list = [ 'lon_rho', 'lat_rho', 'mask_rho', 'h']
for vn in vn_list:
    out_ds[vn] = (('eta_rho', 'xi_rho'), in_ds[vn].values)
    out_ds[vn].attrs['long_name'] = in_ds[vn].attrs['long_name']
    try:
        out_ds[vn].attrs['units'] = in_ds[vn].attrs['units']
    except KeyError:
        pass

# Get ready to initialize the data layers.
if testing:
    depth_list = ['10']
else:
    depth_list = ['surface', '10', '20', '30', '50', 'bottom']

tag_dict = dict()
for depth in depth_list:
    if depth in ['surface', 'bottom']:
        tag_dict[depth] = ' at ' + depth
    else:
        tag_dict[depth] = ' at ' + depth + ' m depth'

if testing:
    # vn_in_list = ['temp', 'salt']
    # vn_out_list = ['temp', 'salt']
    vn_in_list = ['temp', 'salt', 'alkalinity', 'TIC']
    vn_out_list = ['temp', 'salt' , 'PH', 'ARAG']
else:
    # vn_in_list = ['temp', 'salt', 'phytoplankton', 'NO3', 'oxygen', 'rho', 'alkalinity', 'TIC']
    vn_in_list = ['temp', 'salt', 'phytoplankton', 'NO3', 'oxygen', 'alkalinity', 'TIC']
    vn_out_list = ['temp', 'salt', 'phytoplankton', 'NO3', 'oxygen', 'PH', 'ARAG']

do_carbon = True
# this flag is introduced to handle different testgin choices
if 'TIC' not in vn_in_list:
    do_carbon = False

# create zfull to use with the pfun.get_laym() function
zfull = pfun.get_zfull(in_ds, in_fn, 'rho')
in_mask_rho = in_ds.mask_rho.values # 1 = water, 0 = land

def get_layer(vn, depth, in_ds, zfull, in_mask_rho):
    if depth == 'surface':
        L = in_ds[vn][0,-1,:,:].values
    elif depth == 'bottom':
        L = in_ds[vn][0,0,:,:].values
    else:
        L = pfun.get_laym(in_ds, zfull, in_mask_rho, vn, -float(depth))
    return L

def get_Ld(depth, in_ds, in_mask_rho):
    # makes a field of layer depth (m)
    if depth == 'surface':
        Ld = 0 * np.ones(in_mask_rho.shape)
    elif depth == 'bottom':
        Ld = in_ds.h.values
    else:
        Ld = float(depth) * np.ones(in_mask_rho.shape)
    Ld[in_mask_rho == 0] = np.nan
    return Ld

# Fill the layers
for depth in depth_list:
    if testing:
        print(' - Depth = ' + depth)
    suffix = '_' + depth
    v_dict = dict()
    tt0 = time()
    for vn in vn_in_list:
        v_dict[vn] = get_layer(vn, depth, in_ds, zfull, in_mask_rho)
    if testing:
        print('   -- fill v_dict took %0.2f sec' % (time()-tt0))
        
    if do_carbon:
        # ------------- the CO2SYS steps -------------------------
        tt0 = time()
        # ===============================================================
        # OLD
        # Ld = get_Ld(depth, in_ds, in_mask_rho)
        # lat = in_ds.lat_rho.values
        # # create pressure
        # Lpres = sw.pres(Ld, lat)
        # # get in situ temperature from potential temperature
        # Ltemp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, Lpres)
        # # convert from umol/L to umol/kg using in situ dentity
        # Lalkalinity = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
        # Lalkalinity[Lalkalinity < 100] = np.nan
        # LTIC = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
        # LTIC[LTIC < 100] = np.nan
        # CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
        #     Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # ===============================================================
        # NEW: using gsw to create in-situ density
        import gsw
        lon = in_ds.lon_rho.values
        lat = in_ds.lat_rho.values
        Ld = get_Ld(depth, in_ds, in_mask_rho)
        Lpres = gsw.p_from_z(-Ld, lat) # pressure [dbar]
        SA = gsw.SA_from_SP(v_dict['salt'], Lpres, lon, lat)
        CT = gsw.CT_from_pt(SA, v_dict['temp'])
        rho = gsw.rho(SA, CT, Lpres) # in situ density
        Ltemp = gsw.t_from_CT(SA, CT, Lpres) # in situ temperature
        # convert from umol/L to umol/kg using in situ dentity
        Lalkalinity = 1000 * v_dict['alkalinity'] / rho
        Lalkalinity[Lalkalinity < 100] = np.nan
        LTIC = 1000 * v_dict['TIC'] / rho
        LTIC[LTIC < 100] = np.nan
        CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
            Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # ===============================================================
        # NOTE that the "in" and "out" versions of the returned variables will be
        # identical because we pass the same pressure and temperature for
        # input and output (in-situ in both cases)
        PH = CO2dict['pHout']
        v_dict['PH'] = PH.reshape((v_dict['salt'].shape))
        ARAG = CO2dict['OmegaARout']
        v_dict['ARAG'] = ARAG.reshape((v_dict['salt'].shape))
        if testing:
            print('   -- carbon took %0.2f sec' % (time()-tt0))
        # --------------------------------------------------------
    
    # Write data to the output file.
    for vn in vn_out_list:
        out_ds[vn + suffix] = (('ocean_time', 'eta_rho', 'xi_rho'), v_dict[vn].reshape((1,) + in_mask_rho.shape))
        if vn == 'PH':
            out_ds[vn + suffix].attrs['long_name'] = 'pH' + tag_dict[depth]
        elif vn == 'ARAG':
            out_ds[vn + suffix].attrs['long_name'] = 'Aragonite Saturation State' + tag_dict[depth]
        else:
            out_ds[vn + suffix].attrs['long_name'] = in_ds[vn].attrs['long_name'] + tag_dict[depth]
        try:
            out_ds[vn + suffix].attrs['units'] = in_ds[vn].attrs['units']
        except KeyError:
            pass
    sys.stdout.flush()


in_ds.close()

enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}
Enc_dict = {vn:enc_dict for vn in out_ds.data_vars if 'ocean_time' in out_ds[vn].dims}
out_ds.to_netcdf(out_fn, unlimited_dims='ocean_time', encoding=Enc_dict)
# NOTE: we use unlimited_dims so that ncrcat knows what to do later
out_ds.close()