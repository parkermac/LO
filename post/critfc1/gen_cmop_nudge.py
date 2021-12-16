"""
 Generate SELFE nudging file for CRITFC CMOP forecast from LiveOcean data
 
 Charles Seaton (cseaton@crtifc.org), Oct 30, 2020
"""
import netCDF4 as nc
import numpy as np
import os, sys
from dateutil import parser
from scipy import interp
from scipy.spatial import cKDTree
from scipy.io import FortranFile
from scipy import interpolate
import argparse
from datetime import datetime, timedelta
from pathlib import Path

from lo_tools import Lfun, zrfun
from time import time

import critfc_functions as crfun

# get the grid  and startdate
argparser = argparse.ArgumentParser()
argparser.add_argument('hgrid', help='Name of hgrid.gr3 file')
argparser.add_argument('vgrid', help='Name of vgrid.in file')
argparser.add_argument('depthfile', help='Name of ocean_depths.nc file')
argparser.add_argument('basedir', help='Base directory for LiveOcean files')
argparser.add_argument('outdir', help='Output directory for CMOP nudging files')
argparser.add_argument('startdate', help='YYYY-MM-DD start date of nudging file')
argparser.add_argument('--ndays', type=int, help='number of days of forecast, defaults to 3',default=3)
argparser.add_argument('-test', '--testing', type=Lfun.boolean_string, default=False)
args = argparser.parse_args()

hgrid = args.hgrid
vgrid = args.vgrid
depthfile = args.depthfile
startdate = args.startdate
basedir = args.basedir
outdir = args.outdir
ndays = args.ndays
testing = args.testing

try:
    grid = crfun.readHGrid(hgrid, True)
    x = grid.nodes[:,1]
    y = grid.nodes[:,2]
    d = grid.nodes[:,3]
    nnodes = len(x)
except:
    print('''Unable to open file %s''' % (hgrid,))
    sys.exit(1)
    
try:
    vCoords = crfun.vc.fromVGridFile(vgrid)
except Exception as e:
    print(e)
    
try:
    # construct vertical grid
    print('Vert grid')
    Z, kbp2, iwet = vCoords.computeVerticalCoordinates(d * 0.0,d)
except Exception as e:
    print(e)
    print('''Unable to open vertical grid file %s''' % (vgrid,))
    sys.exit(1)

nvrt = vCoords.nvrt

print(hgrid, startdate, nnodes)

#print('ocean depths')
pncid = nc.Dataset(depthfile, 'r');
xrho = pncid.variables['lon_rho'][:]
yrho = pncid.variables['lat_rho'][:]
xr=xrho[:]
yr=yrho[:]

# PM Edit
#depths = pncid.variables['depths'][:]
# make the depths file
G, S, T = zrfun.get_basic_info(depthfile)
depths = -zrfun.get_z(G['h'], 0*G['h'], S, only_rho=True)
depn = depths.shape[0]

# ROMS depths interpolated to SELFE grid
dep=np.zeros((depn,nnodes))
for i in range(0,depn):
    fz=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],depths[i,:,:].T,kx=1,ky=1)
    dep[i,:]=fz.ev(x,y)

(year, month, day) = startdate.split('-')
year = int(year)
month = int(month)
day = int(day)
t1 = parser.parse(startdate+'T00:00:00-00')
t2 = parser.parse(startdate+'T00:00:00-08')
dt= t2 -t1
startf = dt.total_seconds() / 3600.0
nfiles = ndays * 24 + 1

basedir0 = Path(basedir)
#basedir =  basedir0 + '/f%4d.%02d.%02d' % (year,month,day)
#basedir += '''/f%4d.%02d.%02d''' % (year,month,day)
#print(basedir)

# PM edits to make history file list for forecast in 1-day folders
dt0 = datetime(year,month,day)
date_list = []
for dd in [0,1,2,3]:
    this_dt = dt0 + timedelta(days=dd)
    date_list.append(this_dt.strftime(Lfun.ds_fmt))
fn_list = []
for dl in date_list:
    hourmax = 24
    f_string = 'f' + dl
    if dl == date_list[0]:
        hourmin = 0
    else:
        hourmin = 1
        # skip hour zero on subsequent days because it is a repeat
    for nhis in range(hourmin+1, hourmax+2):
        nhiss = ('0000' + str(nhis))[-4:]
        fn = basedir0 / f_string / ('ocean_his_' + nhiss + '.nc')
        fn_list.append(fn)
fn_list = fn_list[8:81]

# initialize output files - NetCDF version
# salt_fn = '%s/cmop_salt_nu.%4d-%02d-%02d.nc' % (outdir,year,month,day)
# temp_fn = '%s/cmop_temp_nu.%4d-%02d-%02d.nc' % (outdir,year,month,day)
salt_fn = '%s/critfc_salt.nc' % (outdir)
temp_fn = '%s/critfc_temp.nc' % (outdir)
print('=== OUTPUT FILE NAMES:=========================')
print(salt_fn)
print(temp_fn)
print('===============================================')
try:
    os.remove(salt_fn)
except OSError:
    pass
try:
    os.remove(temp_fn)
except OSError:
    pass
sf = nc.Dataset(salt_fn, 'w')
tf = nc.Dataset(temp_fn, 'w')
# create dimensions
NZ, NX = Z.shape
sf.createDimension('t', None)
sf.createDimension('z', NZ)
sf.createDimension('x', NX)
tf.createDimension('t', None)
tf.createDimension('z', NZ)
tf.createDimension('x', NX)
# create variables
sf.createVariable('time', np.float32, ('t'))
sf.createVariable('data', np.float32, ('t', 'z', 'x'))
tf.createVariable('time', np.float32, ('t'))
tf.createVariable('data', np.float32, ('t', 'z', 'x'))

# original version - slow and buggy
# outname = '%s/cmop_%%s_nu.%4d-%02d-%02d.in' % (outdir,year,month,day)
# sf = FortranFile(outname % 'salt', 'w')
# tf = FortranFile(outname % 'temp', 'w')
    
temp_prev=np.zeros((depn,nnodes))
salt_prev=np.zeros((depn,nnodes))
print('processing')

if testing:
    Nfiles = 3
else:
    Nfiles = nfiles
    
# if testing:
#     print(len(fn_list))
#     # this little section is to develop the code so that we map onto the new format
#     # where the forecast is in individual day folders.
#     for i in range(nfiles):
#         print(str(fn_list[i]))
#     sys.exit()
    
for i in range(Nfiles):
# there is no file #1 use #2 twice
    tt0 = time()
    
    ncfile = str(fn_list[i])
    print(ncfile)
    
    tfile = startf + i + 1 #add 1 because there is an offset between file # and time step
    # ncfile = '''%s/ocean_his_%04d.nc''' % (basedir, tfile,)
    
    otime = i * 3600
    #for i in range(len(x)):
    temp_new=np.zeros((depn,nnodes))
    salt_new=np.zeros((depn,nnodes))
    
    print('Step %d: salt_new.max = %0.2f' % (1, np.nanmax(salt_new)))
    
    try:
        psncid = nc.Dataset(ncfile, 'r');
        ocean_time = psncid.variables['ocean_time'][0]
        salt = psncid.variables['salt'][0]
        temp = psncid.variables['temp'][0]
        
        # This ensures that we feed huge values to RectBivariateSpline for masked places
        # which are then trimmed by hand in the "salt_new>35" step.
        Salt = salt.data
        Salt[salt.mask] = 1e37
        Temp = temp.data
        Temp[temp.mask] = 1e37
        # With the original LiveOcean output (NetCDF3) I think RectBivariateSpline would
        # misinterpret the mask as 1e37, but with NetCDF4 (from the compression step, but
        # presumably also maybe with the cas6_v0_u0kb output) it interprets the masked
        # values as nans, and the interpolation fails.
        
        # In general, this is not how I would code this step!  A kDTree would be cleaner.
        
        psncid.close()
        
        for j in range(0,depn):
            ft=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],Temp[j,:,:].T,kx=1,ky=1)
            temp_new[j,:]=ft.ev(x,y)
            fs=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],Salt[j,:,:].T,kx=1,ky=1)
            salt_new[j,:]=fs.ev(x,y)
        
        print('Step %d: salt_new.max = %0.2f' % (2, np.nanmax(salt_new)))
        
        temp_new[np.where(salt_new>35)]=float('nan')
        salt_new[np.where(salt_new>35)]=float('nan')
        
        print('Step %d: salt_new.max = %0.2f' % (3, np.nanmax(salt_new)))
        
        good=~np.isnan(salt_new[0,:])
        bad=np.isnan(salt_new[0,:])
        xygood=np.array((x[good],y[good])).T
        xybad=np.array((x[bad],y[bad])).T    
        ind=np.arange(nnodes)
        map_tr=ind[good]
        #kdt=cKDTree(zip(x[good],y[good]))
        XY = np.array((x[good],y[good])).T
        kdt=cKDTree(XY)
        dists,neighs=kdt.query(xybad)
        min_index=map_tr[neighs]
        
        ix_min=np.unravel_index(min_index,x.shape)
        salt_new[:,bad]=salt_new[:,ix_min][:,0,:]
        temp_new[:,bad]=temp_new[:,ix_min][:,0,:]
        
        print('Step %d: salt_new.max = %0.2f' % (4, np.nanmax(salt_new)))
        
        salt_prev = salt_new
        temp_prev = temp_new
        
        
    except Exception as e:
        print(e)
        print('Exception in loading new, using previous time step')
        salt_new = salt_prev
        temp_new = temp_prev
                
    print(i+1, ' of ', nfiles, tfile, otime, ocean_time, ncfile)
    
    # write time to output files
    sf['time'][i] = otime
    tf['time'][i] = otime
    # sf.write_record(np.array([otime], dtype=np.float32))
    # tf.write_record(np.array([otime], dtype=np.float32))
    
    # ** the slowest step? **
    s1arr = np.empty((NZ,NX), dtype=np.float32)
    t1arr = np.empty((NZ,NX), dtype=np.float32)
    for j in range(len(x)):
        s1=np.zeros(nvrt)
        t1=np.zeros(nvrt)
        ss=np.interp(Z[:,j],dep[:,j],salt_new[:,j])
        tt=np.interp(Z[:,j],dep[:,j],temp_new[:,j])
        s1[Z[:,j].mask==False]=ss[Z[:,j].mask==False]
        t1[Z[:,j].mask==False]=tt[Z[:,j].mask==False]
        
        s1arr[:,j] = s1
        t1arr[:,j] = t1
        # sf.write_record(s1.astype(np.float32))
        # tf.write_record(t1.astype(np.float32))
            
    sf['data'][i,:,:] = s1arr
    tf['data'][i,:,:] = t1arr

    print(' >> time to process this hour %0.1f seconds' % (time()-tt0))
    sys.stdout.flush()
    
sf.close()
tf.close()
    
