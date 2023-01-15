"""
This is refactored from ccast_functions.get_cast() to be a stand-alone
program run by command line arguments, so that we can parallelize extract_casts.py
using subprocess.

"""

from lo_tools import Lfun, zfun, zrfun
import subprocess
import xarray as xr
import numpy as np
from datetime import timedelta
import argparse
import sys
from pathlib import Path

parser = argparse.ArgumentParser()

# select mooring location and name
parser.add_argument('-out_fn', type=str)   # full path to output file
parser.add_argument('-fn', type=str)   # full path to history file
parser.add_argument('-lon', type=str) # longitude
parser.add_argument('-lat', type=str) # latitude
parser.add_argument('-npzd', type=str) # new, old, or none
args = parser.parse_args()

out_fn = Path(args.out_fn)
fn = Path(args.fn)
lon = float(args.lon)
lat = float(args.lat)
npzd = args.npzd
    
# This function does the cast extraction and saves it to a NetCDF file.
G, S, T = zrfun.get_basic_info(fn)
Lon = G['lon_rho'][0,:]
Lat = G['lat_rho'][:,0]

# error checking
if (lon < Lon[0]) or (lon > Lon[-1]):
    print('ERROR: lon out of bounds ' + out_fn.name)
    sys.exit()
if (lat < Lat[0]) or (lat > Lat[-1]):
    print('ERROR: lat out of bounds ' + out_fn.name)
    sys.exit()

ix = zfun.find_nearest_ind(Lon, lon)
iy = zfun.find_nearest_ind(Lat, lat)

# error checking
if G['mask_rho'][iy,ix] == 0:
    print('ERROR: point on land mask ' + out_fn.name)
    sys.exit()
    
v_list = 'AKs,salt,temp,h'
if npzd == 'new':
    v_list += ',phytoplankton,chlorophyll,zooplankton,SdetritusN,LdetritusN,oxygen,alkalinity,TIC,NO3,NH4'
elif npzd == 'old':
    v_list += ',phytoplankton,zooplankton,detritus,Ldetritus,oxygen,alkalinity,TIC,NO3'
elif npzd == 'none':
    pass
else:
    print('Error: unrecognized npzd')
    sys.exit()

# Run ncks to do the extraction, overwriting any existing file
cmd_list = ['ncks', '-d', 'xi_rho,'+str(ix), '-d', 'eta_rho,'+str(iy),
    '-v', v_list, '-O', str(fn), str(out_fn)]
# Note: We get AKs so that the s_w dimension is retained
proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# and check on the results
stdout, stderr = proc.communicate()
if len(stdout) > 0:
    print('\n' + ' sdtout '.center(60,'-'))
    print(stdout.decode())
if len(stderr) > 0:
    print('\n' + ' stderr '.center(60,'-'))
    print(stderr.decode())
    
# Add z-coordinates to the file using xarray
foo = xr.load_dataset(out_fn)
foo = foo.squeeze()
z_rho, z_w = zrfun.get_z(foo['h'].values, np.array([0.]), S)
foo['z_rho'] = (('s_rho'), z_rho)
foo['z_w'] = (('s_w'), z_w)
foo.s_rho.attrs['long_name'] = 'vertical position on s_rho grid, positive up, zero at surface'
foo.s_rho.attrs['units'] = 'm'
foo.s_w.attrs['long_name'] = 'vertical position on s_w grid, positive up, zero at surface'
foo.s_w.attrs['units'] = 'm'
foo.salt.attrs['units'] = 'g kg-1'
foo.to_netcdf(out_fn)
foo.close()