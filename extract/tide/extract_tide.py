"""
Code to extract tide height time series from multiple locations,
from a user-specified model run.

Test on my mac:
run extract_tide -test True

Run for real on apogee:
python extract_tide.py -gtx cas7_t0_x4b -ro 2 -0 2017.01.01 -1 2017.12.31 > test.log &

Performance:
20 sec for three days, both mac and apogee
2000 sec for a year on apogee

"""

# imports
import sys
from lo_tools import Lfun, zrfun, zfun
import argparse
from time import time
import numpy as np
import xarray as xr
import pandas as pd

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', default='cas7_t0_x4b', type=str)   # e.g. cas7_t0_x4b
parser.add_argument('-ro', '--roms_out_num', default=0, type=int) # 1 = Ldir['roms_out1'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', default='2017.07.04', type=str)
parser.add_argument('-1', '--ds1', default='2017.07.06', type=str)
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=zfun.boolean_string)    
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided (other)
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]

year_str = Ldir['ds0'][:4]
info_year_str = year_str
# testing
if Ldir['testing']:
    info_year_str = '2022'

# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

# do the extractions

# set input and output locations
in_dir = Ldir['LOo'] / 'obs' / 'tide' # where to find station location info
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tide'
Lfun.make_dir(out_dir)

# To keep things clear I will stick to the ordering convention:
# (lat,lon) (j,i) (row,col)

# combine station info into a DataFrame, with station number as the index
sn_df = pd.DataFrame(columns=['name','lat','lon','jgood','igood'])

sn_df_noaa = pd.read_pickle(in_dir / ('sn_df_noaa_' + info_year_str + '.p'))
sn_df_dfo = pd.read_pickle(in_dir / ('sn_df_dfo_' + info_year_str + '.p'))
for sn in sn_df_noaa.index:
    sn_df.loc[sn,['name','lat','lon']] = sn_df_noaa.loc[sn,['name','lat','lon']]
for sn in sn_df_dfo.index:
    sn_df.loc[sn,['name','lat','lon']] = sn_df_dfo.loc[sn,['name','lat','lon']]

# get indices for extraction
in_dir0 = Ldir['roms_out'] / Ldir['gtagex']
G, S, T = zrfun.get_basic_info(in_dir0 / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')

def get_ji_good(G, lat00, lon00):
    # Function to find the indices of the nearest unmasked point.
    
    # G is the grid info dict
    # lat00, lon00 is the mooring/cast location

    # useful things
    lat = G['lat_rho']
    lon = G['lon_rho']
    yvec = G['lat_rho'][:,0]
    xvec = G['lon_rho'][0,:]
    mask = G['mask_rho']

    # Out of bounds error checking
    if (lat00 < yvec[0]) or (lat00 > yvec[-1]):
        print(' ERROR: lat00 out of bounds ')
        sys.exit()
    if (lon00 < xvec[0]) or (lon00 > xvec[-1]):
        print(' ERROR: lon00 out of bounds ')
        sys.exit()

    # initial guess
    j0 = zfun.find_nearest_ind(yvec, lat00)
    i0 = zfun.find_nearest_ind(xvec, lon00)

    # starting point
    lat0 = lat[j0,i0]
    lon0 = lon[j0,i0]

    pad = 5 # how far to look (points)

    # indices of box to search over
    jmax = len(yvec)-1
    imax = len(xvec)-1
    J = np.arange(j0-pad,j0+pad)
    I = np.arange(i0-pad, i0+pad)

    # account for out-of-range points (not needed?)
    if I[0] < 0:
        I = I - I[0]
    if I[-1] > imax:
        I = I - (I[-1] - imax)
    if J[0] < 0:
        J = J - J[0]
    if J[-1] > jmax:
        J = J - (J[-1] - jmax)

    # array of indices to search
    ii, jj = np.meshgrid(I, J)

    # sub arrays, distances
    llat = lat[jj,ii]
    llon = lon[jj,ii]
    xxx, yyy = zfun.ll2xy(llon, llat, lon0, lat0)
    ddd = np.sqrt(xxx**2 + yyy**2) # distance from original point
    mmask = mask[jj,ii] 
    mm = mmask==1 # Boolean array of good points
    dddm = ddd[mm] # vector of good distances

    # indices of best point
    jgood = jj[mm][dddm==dddm.min()][0]
    igood = ii[mm][dddm==dddm.min()][0]
    
    return jgood, igood

for sn in sn_df.index:
    # Find the closest good grid indices and pack everything
    # in a single DataFrame.
    lat00 = sn_df.loc[sn,'lat']
    lon00 = sn_df.loc[sn,'lon']
    name = sn_df.loc[sn,'name']
    # get indices
    jgood, igood = get_ji_good(G, lat00, lon00)
    sn_df.loc[sn,'jgood'] = jgood
    sn_df.loc[sn,'igood'] = igood
    print('%s (j,i)=(%d,%d)' % (name, jgood, igood))

if Ldir['testing']:
    # a plot to check on things
    import matplotlib.pyplot as plt
    from lo_tools import plotting_functions as pfun
    plt.close('all')
    pfun.start_plot(figsize=(12,8))
    plon, plat = pfun.get_plon_plat(G['lon_rho'],G['lat_rho'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    zz = -G['h']
    zz[G['mask_rho']==0] = np.nan
    cs = ax.pcolormesh(plon,plat,zz)
    pfun.dar(ax)
    pfun.add_coast(ax)
    aa = [plon[0,0], plon[0,-1], plat[0,0], plat[-1,0]]
    ax.axis(aa)
    for sn in sn_df.index:
        # Find the closest good grid indices and pack everything
        # in a single DataFrame.
        lat00 = sn_df.loc[sn,'lat']
        lon00 = sn_df.loc[sn,'lon']
        name = sn_df.loc[sn,'name']
        jgood = sn_df.loc[sn,'jgood']
        igood = sn_df.loc[sn,'igood']
        ax.plot(lon00,lat00,'or')
        ax.plot(G['lon_rho'][jgood,igood],G['lat_rho'][jgood,igood],'ob')
        ax.text(lon00,lat00,name)
    plt.show()
    # RESULT: looks good for NOAA stations

# generate list of files to open
fn_list = Lfun.get_fn_list('hourly', Ldir, Ldir['ds0'], Ldir['ds1'])

# initialize a DataFrame to hold results
ssh_df = pd.DataFrame(columns=sn_df.index)
# do the extraction
tt0 = time()
jj = sn_df.jgood.to_numpy(dtype=int)
ii = sn_df.igood.to_numpy(dtype=int)
for fn in fn_list:
    ds = xr.open_dataset(fn)
    ssh = ds.zeta[0,:,:].values
    # extract using fancy indexing
    sn_ssh = ssh[jj,ii]
    tt = pd.Timestamp(ds.ocean_time.values[0])
    ssh_df.loc[tt,:] = sn_ssh
    ds.close()
print('\nTime to do extractions = %0.1f sec' % (time()-tt0))
# performance: 21 sec for three days on my mac

# save the output in a single file
if Ldir['testing'] == False:
    ssh_df.to_pickle(out_dir / ('ssh_df_' + year_str + '.p'))
# and save station info
sn_df.to_pickle(out_dir / ('sn_df_' + year_str + '.p'))
sn_df.to_csv(out_dir / ('sn_df_' + year_str + '.csv'))
