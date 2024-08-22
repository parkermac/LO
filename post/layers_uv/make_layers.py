"""
This creates a layers file for a single time.

It is designed for the needs of some APL glider researchers who
want velocities in the Strait of Juan de Fuca.

One difference from layers1 is that it
packs the fields as 3-D (t,z,y,x) instead of individual files for
individual layers. Time is a singleton dimension, but unlimited in
case we wanted to concatenate later.

To test on mac:

run make_layers -test True
"""

from pathlib import Path
import argparse
import xarray as xr
from time import time
import numpy as np
import sys

from lo_tools import plotting_functions as pfun
from lo_tools import zfun

# defaults for testing
from lo_tools import Lfun
Ldir = Lfun.Lstart()
in_fn = (Ldir['roms_out'] /
    'cas7_t0_x4b' / 'f2017.07.04' / 'ocean_his_0002.nc')
in_num_str = in_fn.name[-7:-3]
out_fn = (Ldir['parent'] / 'LO_output' / 'post' /
    'cas7_t0_x4b' / 'f2017.07.04' / 'layers_uv' / ('layers_' + in_num_str + '.nc'))

# command line arguments, passed from post_main.py
parser = argparse.ArgumentParser()
parser.add_argument('-in_fn', default=in_fn, type=str)
parser.add_argument('-out_fn', default=out_fn, type=str)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)

args = parser.parse_args()
in_fn = Path(args.in_fn)
out_fn = Path(args.out_fn)
testing = args.testing

# make sure output directory exists
out_dir = out_fn.parent
Lfun.make_dir(out_dir)
out_fn.unlink(missing_ok=True)

"""
The algorithm strategy is to first use ncks to create a version of the
history file that oonly has a specific region and only the needed variables.
This makes things much faster when we go to extract layers, and it makes
the output files much smaller.
"""
# specify region to extract and get its indices
#
# just the Strait of Juan de Fuca
aa = [-125.3, -122.3, 47.9, 48.8]
#
in_ds = xr.open_dataset(in_fn, decode_times=False)
Lon = in_ds.lon_rho[0,:].values
Lat = in_ds.lat_rho[:,0].values
def check_bounds(lon, lat):
    # error checking
    if (lon < Lon[0]) or (lon > Lon[-1]):
        print('ERROR: lon out of bounds ')
        sys.exit()
    if (lat < Lat[0]) or (lat > Lat[-1]):
        print('ERROR: lat out of bounds ')
        sys.exit()
    # get indices
    ilon = zfun.find_nearest_ind(Lon, lon)
    ilat = zfun.find_nearest_ind(Lat, lat)
    return ilon, ilat
lon0, lon1, lat0, lat1 = aa
ilon0, ilat0 = check_bounds(lon0, lat0)
ilon1, ilat1 = check_bounds(lon1, lat1)
in_ds.close()

Ncenter = 30
verbose = False
def messages(stdout, stderr, mtitle, verbose):
    # utility function for displaying subprocess info
    if verbose:
        print((' ' + mtitle + ' ').center(Ncenter,'='))
        if len(stdout) > 0:
            print(' sdtout '.center(Ncenter,'-'))
            print(stdout.decode())
    if len(stderr) > 0:
        print((' ' + mtitle + ' ').center(Ncenter,'='))
        # always print errors
        print(' stderr '.center(Ncenter,'-'))
        print(stderr.decode())
    sys.stdout.flush()

# Here we do the ncks extraction to get a smaller file to work with. Note
# that we take care to have all the grid limits just like in a history file,
# meaning that the rho-grid extends to the corners and the u, v, and psi
# grids are smaller.
out_fn_small = out_dir / out_fn.name.replace('layers','small')
# use ncks to extract a sub-domain
from subprocess import Popen as Po
from subprocess import PIPE as Pi
vn_list = 'h,f,pm,pn,mask_rho,mask_u,mask_v,mask_psi,salt,temp,zeta,u,v,w,Vtransform'
cmd_list1 = ['ncks',
    '-v', vn_list,
    '-d', 'xi_rho,'+str(ilon0)+','+str(ilon1), '-d', 'eta_rho,'+str(ilat0)+','+str(ilat1),
    '-d', 'xi_u,'+str(ilon0)+','+str(ilon1-1), '-d', 'eta_u,'+str(ilat0)+','+str(ilat1),
    '-d', 'xi_v,'+str(ilon0)+','+str(ilon1), '-d', 'eta_v,'+str(ilat0)+','+str(ilat1-1),
    '-d', 'xi_psi,'+str(ilon0)+','+str(ilon1-1), '-d', 'eta_psi,'+str(ilat0)+','+str(ilat1-1),
    '--mk_rec_dim', 'ocean_time']
cmd_list1 += ['-O', str(in_fn), str(out_fn_small)]
proc = Po(cmd_list1, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
messages(stdout, stderr, 'ncks', verbose)

# rename "in_fn" so that now we will do the rest of the work using our small file.
in_fn = out_fn_small
in_ds = xr.open_dataset(in_fn, decode_times=False)
# NOTE: now "in_ds" refers to our small temporary file

# initiate the output Dataset
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

# Specify the z-levels to use
if testing:
    z_vec = np.array([-4000,-50,-10,0]) 
else:
    z_vec = np.linspace(-300,0,16)

# add z as a coordinate and as a data_var
out_ds['z'] = (('z'), z_vec)
out_ds['z'].attrs['units'] = 'm'
out_ds['z'].attrs['long_name'] = 'vertical coordinate relative to flat sea surface, positive up'

# set list of variables we will interpolate to layers
vn_in_list = ['temp', 'salt', 'u', 'v']
vn_out_list = ['temp', 'salt' , 'ur', 'vr']
# we introduce the variable names ur and vr for the u and v velocities
# interpolated to the rho grid.

# create zfull to use with the pfun.get_laym() function
zfull_rho = pfun.get_zfull(in_ds, in_fn, 'rho')
zfull_u = pfun.get_zfull(in_ds, in_fn, 'u')
zfull_v = pfun.get_zfull(in_ds, in_fn, 'v')
# also get masks
in_mask_rho = in_ds.mask_rho.values # 1 = water, 0 = land
in_mask_u = in_ds.mask_u.values # 1 = water, 0 = land
in_mask_v = in_ds.mask_v.values # 1 = water, 0 = land

def get_layer(vn, z, in_ds, zfull, in_mask):
    # convenience function that wraps pfun.get_laym() and makes sure
    # that z = 0 is interpreted as the topmost (-1) bin.
    if vn in ['salt','temp']:
        if z == 0:
            L = in_ds[vn][0,-1,:,:].values
        else:
            L = pfun.get_laym(in_ds, zfull_rho, in_mask_rho, vn, z)
    elif vn == 'u':
        if z == 0:
            L = in_ds[vn][0,-1,:,:].values
        else:
            L = pfun.get_laym(in_ds, zfull_u, in_mask_u, vn, z)
    elif vn == 'v':
        if z == 0:
            L = in_ds[vn][0,-1,:,:].values
        else:
            L = pfun.get_laym(in_ds, zfull_v, in_mask_v, vn, z)
    return L

# prepare a dict of nan arrays for output
out_dict = dict()
NZ = len(z_vec)
NY, NX = in_ds.lon_rho.shape
out_mat = np.nan * np.ones((1,NZ,NY,NX))
for vn in vn_out_list:
    out_dict[vn] = out_mat.copy()

# Fill the output arrays by interpolating to each layer z. We also interpolate
# u and v to be on the rho-grid. This leaves a band of nan on the outer edge
# of the ur and vr fields, but we assume the box is bigger than what the user
# needs.
for vn in vn_in_list:
    tt0 = time()
    ii = 0
    for z in z_vec:
        if vn in ['salt','temp']:
            out_dict[vn][0,ii,:,:] = get_layer(vn, z, in_ds, zfull_rho, in_mask_rho)
        elif vn == 'u':
            u = get_layer(vn, z, in_ds, zfull_u, in_mask_u)
            out_dict['ur'][0,ii,1:-1,1:-1] = u[1:-1,:-1] + (u[1:-1,1:]-u[1:-1,:-1])/2
        elif vn == 'v':
            v = get_layer(vn, z, in_ds, zfull_v, in_mask_v)
            out_dict['vr'][0,ii,1:-1,1:-1] = v[:-1,1:-1] + (v[1:,1:-1]-v[:-1,1:-1])/2
        ii += 1
    if testing:
        print('   -- vn = %s: fill out_dict took %0.2f sec' % (vn,time()-tt0))
        sys.stdout.flush()
    
# Write data to the output file.
for vn in vn_out_list:
    out_ds[vn] = (('ocean_time', 'z', 'eta_rho', 'xi_rho'), out_dict[vn])
    if vn == 'ur':
        out_ds[vn].attrs['long_name'] = 'u-velocity, positive to east'
        out_ds[vn].attrs['units'] = 'm s-1'
    elif vn == 'vr':
        out_ds[vn].attrs['long_name'] = 'v-velocity, positive to north'
        out_ds[vn].attrs['units'] = 'm s-1'
    else:
        out_ds[vn].attrs['long_name'] = in_ds[vn].attrs['long_name']
        try:
            out_ds[vn].attrs['units'] = in_ds[vn].attrs['units']
        except KeyError:
            pass
out_ds['salt'].attrs['units'] = 'PSU'


# compress and write to output NetCDF file
enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}
Enc_dict = {vn:enc_dict for vn in out_ds.data_vars if 'ocean_time' in out_ds[vn].dims}
out_ds.to_netcdf(out_fn, unlimited_dims='ocean_time', encoding=Enc_dict)
# NOTE: we use unlimited_dims so that ncrcat knows what to do later, if we decide
# to concatenate to a larger file.

# close open Datasets
in_ds.close()
out_ds.close()

if testing == False:
    # clean up temp file
    in_fn.unlink(missing_ok=True)

if testing:
    # a visual check that things are not way off

    out_ds = xr.open_dataset(out_fn)
    # plotting
    import matplotlib.pyplot as plt
    plt.close('all')
    pfun.start_plot(figsize=(20,10))
    lon = out_ds.lon_rho.values
    lat = out_ds.lat_rho.values
    plon, plat = pfun.get_plon_plat(lon,lat)

    lim_dict = {'salt':(25,34),'temp':(8,18),'ur':(-1,1),'vr':(-1,1)}
    cmap_dict = {'salt':'Spectral_r','temp':'jet','ur':'RdBu_r','vr':'RdBu_r'}

    for vn in vn_out_list:
        fig = plt.figure()
        vmin, vmax = lim_dict[vn]
        NT, NZ, NY, NX = out_ds[vn].shape
        zz = out_ds.z.values
        ii = 0
        for z in zz:
            ax = fig.add_subplot(2,2,ii+1)
            cs = ax.pcolormesh(plon,plat,out_ds[vn][0,ii,:,:].values,
                vmin=vmin,vmax=vmax,cmap=cmap_dict[vn])
            fig.colorbar(cs,ax=ax)
            ax.set_title('%s, z = %d m' % (vn,z))
            pfun.dar(ax)
            ii += 1

    plt.show()
    pfun.end_plot()

