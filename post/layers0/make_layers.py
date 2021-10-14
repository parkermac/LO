"""
This creates a layers file for a subset of the fn_list.  The goal is to run
several instances of this in parallel to speed up the process.

To test on mac:

run make_layers -testing True -verbose True
"""

# imports
import sys, os
sys.path.append(os.path.abspath('../../alpha/'))
import zfun
sys.path.append(os.path.abspath('../../plotting/'))
import pfun

import layer_fun

import argparse
import netCDF4 as nc
from time import time
from PyCO2SYS import CO2SYS
import seawater as sw
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
import numpy as np

# defaults for standalone testing
in_dir = '/Users/pm8/Documents/LiveOcean_roms/output/cas6_v3_lo8b/f2019.07.04/'

parser = argparse.ArgumentParser()
parser.add_argument('-in_dir', type=str, default=in_dir)
parser.add_argument('-nn', type=int, default=0)
parser.add_argument('-istart', type=int, default=0)
parser.add_argument('-iend', type=int, default=0)
parser.add_argument('-testing', type=zfun.boolean_string, default=True)
parser.add_argument('-verbose', type=zfun.boolean_string, default=True)

args = parser.parse_args()
in_dir = args.in_dir
nn = args.nn
istart = args.istart
iend = args.iend
testing = args.testing
verbose = args.verbose

# ======== Create an output file for SCOOT, NANOOS and etc. =============================

fn_list = layer_fun.get_fn_list(in_dir, testing)
    
# trim fn_list
fn_list = fn_list[istart:iend+1]

# Initialize
out_dir = in_dir
out_fn = out_dir + 'temp_layer' + str(nn) + '.nc'
if verbose:
    print(' - Writing to: ' + out_fn)
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
out_ds = nc.Dataset(out_fn, 'w')

# Copy dimensions
dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
in_ds = nc.Dataset(fn_list[0])
for dname, the_dim in in_ds.dimensions.items():
    if dname in dlist:
        out_ds.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# Copy grid variables
for vn in vn_list2:
    varin = in_ds[vn]
    vv = out_ds.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    vv[:] = in_ds[vn][:]
    
# Create and copy time variable
vn = 'ocean_time'
varin = in_ds[vn]
vv = out_ds.createVariable(vn, varin.dtype, varin.dimensions)
vv.long_name = varin.long_name
vv.units = varin.units
in_ds.close()
counter = 0
for fn in fn_list:
    in_ds = nc.Dataset(fn)
    vv[counter] = in_ds[vn][:]
    counter += 1
    in_ds.close()

# Get ready to initialize the data layers.
if testing:
    depth_list = ['10']
else:
    depth_list = ['surface', '10', '20', '30', '50', 'bottom']
    
tag_dict = {}
for depth in depth_list:
    if depth in ['surface', 'bottom']:
        tag_dict[depth] = ' at ' + depth
    else:
        tag_dict[depth] = ' at ' + depth + ' m depth'

if testing:
    vn_in_list = ['temp', 'salt' , 'rho', 'alkalinity', 'TIC']
    vn_out_list = ['temp', 'salt' , 'PH', 'ARAG']
else:
    vn_in_list = ['temp', 'salt', 'phytoplankton', 'NO3', 'oxygen', 'rho', 'alkalinity', 'TIC']
    vn_out_list = ['temp', 'salt', 'phytoplankton', 'NO3', 'oxygen', 'PH', 'ARAG']
    
# Create the data layer objects
in_ds = nc.Dataset(fn_list[0])
for vn in vn_out_list:
    try:
        varin = in_ds[vn]
    except IndexError:
        pass
    for depth in depth_list:
        suffix = '_' + depth
        vv = out_ds.createVariable(vn + suffix, float, ('ocean_time', 'eta_rho', 'xi_rho'))
        if vn == 'PH':
            vv.long_name = 'pH' + tag_dict[depth]
        elif vn == 'ARAG':
            vv.long_name = 'Aragonite Saturation State' + tag_dict[depth]
        else:
            vv.long_name = varin.long_name + tag_dict[depth]
        try:
            vv.units = varin.units
        except AttributeError:
            pass
        vv.time = 'ocean_time'
in_ds.close()

# create zfull to use with the pfun.get_laym() function
fn0 = fn_list[0]
in_ds = nc.Dataset(fn0)
zfull = pfun.get_zfull(in_ds, fn0, 'rho')
in_mask_rho = in_ds['mask_rho'][:] # 1 = water, 0 = land
in_ds.close()

def get_layer(vn, depth, in_ds, zfull, in_mask_rho):
    if depth == 'surface':
        L = zfun.fillit(in_ds[vn][0,-1,:,:])
    elif depth == 'bottom':
        L = zfun.fillit(in_ds[vn][0,0,:,:])
    else:
        L = pfun.get_laym(in_ds, zfull, in_mask_rho, vn, -float(depth))
        L = zfun.fillit(L)
    return L
    
def get_Ld(depth, in_ds, in_mask_rho):
    # makes a field of layer depth (m)
    if depth == 'surface':
        Ld = 0 * np.ones_like(in_mask_rho)
    elif depth == 'bottom':
        Ld = in_ds['h'][:]
    else:
        Ld = float(depth) * np.ones_like(in_mask_rho)
    Ld[in_mask_rho == 0] = np.nan
    return Ld

# Fill the layers
tt = 0
for in_fn in fn_list:
    if verbose:
        print('Working on: ' + in_fn.split('/')[-2] + ' ' + in_fn.split('/')[-1])
    in_ds = nc.Dataset(fn)
    for depth in depth_list:
        if verbose:
            print(' - Depth = ' + depth)
        suffix = '_' + depth
        v_dict = {}
        tt0 = time()
        for vn in vn_in_list:
            v_dict[vn] = get_layer(vn, depth, in_ds, zfull, in_mask_rho)
        if verbose:
            print('   -- fill v_dict took %0.2f sec' % (time()-tt0))
        # ------------- the CO2SYS steps -------------------------
        tt0 = time()
        Ld = get_Ld(depth, in_ds, in_mask_rho)
        lat = in_ds['lat_rho'][:]
        # create pressure
        Lpres = sw.pres(Ld, lat)
        # get in situ temperature from potential temperature
        Ltemp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, Lpres)
        # convert from umol/L to umol/kg using in situ dentity
        Lalkalinity = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
        Lalkalinity[Lalkalinity < 100] = np.nan
        LTIC = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
        LTIC[LTIC < 100] = np.nan
        Lpres = zfun.fillit(Lpres)
        Ltemp = zfun.fillit(Ltemp)
        CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
            Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # NOTE that the "in" and "out" versions of the returned variables will be
        # identical because we pass the same pressure and temperature for
        # input and output (in-situ in both cases)
        PH = CO2dict['pHout']
        v_dict['PH'] = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
        ARAG = CO2dict['OmegaARout']
        v_dict['ARAG'] = zfun.fillit(ARAG.reshape((v_dict['salt'].shape)))
        if verbose:
            print('   -- carbon took %0.2f sec' % (time()-tt0))
        # --------------------------------------------------------
        # Write data to the output file.
        for vn in vn_out_list:
            out_ds[vn + suffix][tt,:,:] = v_dict[vn]
        sys.stdout.flush()
    in_ds.close()
    tt += 1

# Close the output file
out_ds.close()