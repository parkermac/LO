"""
Code to test TRAPS forcing in river.nc file

From ipython:
run test_traps_forcing.py -gridname cas7 -frc trapsN00 -dstr0 2019.07.10 -dstr1 2019.07.10
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
import pandas as pd
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
from time import time
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

#################################################################################
#                               Argument parsing                                #
#################################################################################

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# info to find a rivers.nc file
parser.add_argument('-gridname', default='cas7', type=str)
parser.add_argument('-frc', default='trapsV00', type=str)
parser.add_argument('-dstr0',default='2017.01.01', type=str)
parser.add_argument('-dstr1',default='2017.01.31', type=str)
args = parser.parse_args()
gridname = args.gridname
dstr0 = args.dstr0
dstr1 = args.dstr1
frc = args.frc

#################################################################################
#                                   Get data                                    #
#################################################################################

# input and output locations
Ldir = Lfun.Lstart(gridname=gridname)
grid_fn = Ldir['grid'] / 'grid.nc'
riv_fn = Ldir['LOo'] / 'forcing' / gridname / ('f' + dstr0) / frc / 'rivers.nc'
out_dir = Ldir['LOo'] / 'testing' / 'traps' / (gridname + '_' + frc)
Lfun.make_dir(out_dir)
    
# get grid data
ds = xr.open_dataset(grid_fn)
h = ds.h.values
m = ds.mask_rho.values
# depth for plotting
h[m==0] = np.nan
# coordinates for plotting
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
aa = pfun.get_aa(ds)
# coordinates for convenience in plotting
lor = ds.lon_rho[0,:].values
lar = ds.lat_rho[:,0].values
lou = ds.lon_u[0,:].values
lau = ds.lat_u[:,0].values
lov = ds.lon_v[0,:].values
lav = ds.lat_v[:,0].values
ds.close

# get river info
rds = xr.open_dataset(riv_fn)

#################################################################################
#         Check river direction (LO/extract/tef2/create_river_info.py)          #
#################################################################################

# Pull out river info into a pandas DataFrame
df = pd.DataFrame(columns=['name','ii','jj','dir','sgn', 'iu','ju', 'iv','jv', 'irho','jrho'])
df.loc[:,'name'] = list(rds.river_name.values)
df.loc[:,'ii'] = list(rds.river_Xposition.values)
df.loc[:,'jj'] = list(rds.river_Eposition.values)
df.loc[:,'dir'] = list(rds.river_direction.values) # 0 = u-grid, 1 = v-grid, 2 = rho-grid

# also need to know if the source is positive or negative in order to assign it to the
# rho grid point it flows into
q = rds.river_transport.values
df.loc[:,'sgn'] = np.sign(q.mean(axis=0))

# translate u-grid sources to dataset indices
df.loc[df.dir==0,'iu'] = df.loc[df.dir==0,'ii'] - 1
df.loc[df.dir==0,'ju'] = df.loc[df.dir==0,'jj']

# translate v-grid sources to dataset indices
df.loc[df.dir==1,'iv'] = df.loc[df.dir==1,'ii']
df.loc[df.dir==1,'jv'] = df.loc[df.dir==1,'jj'] - 1

# use the sign of the source to figure out which rho grid point to associate with the source
# i.e. which rho point does it flow into

df.loc[(df.dir==0) & (df.sgn==1), 'irho'] = df.loc[(df.dir==0) & (df.sgn==1), 'iu'] + 1
df.loc[(df.dir==0) & (df.sgn==1), 'jrho'] = df.loc[(df.dir==0) & (df.sgn==1), 'ju'] 

df.loc[(df.dir==0) & (df.sgn==-1), 'irho'] = df.loc[(df.dir==0) & (df.sgn==-1), 'iu']
df.loc[(df.dir==0) & (df.sgn==-1), 'jrho'] = df.loc[(df.dir==0) & (df.sgn==-1), 'ju']

df.loc[(df.dir==1) & (df.sgn==1), 'irho'] = df.loc[(df.dir==1) & (df.sgn==1), 'iv']
df.loc[(df.dir==1) & (df.sgn==1), 'jrho'] = df.loc[(df.dir==1) & (df.sgn==1), 'jv'] + 1

df.loc[(df.dir==1) & (df.sgn==-1), 'irho'] = df.loc[(df.dir==1) & (df.sgn==-1), 'iv']
df.loc[(df.dir==1) & (df.sgn==-1), 'jrho'] = df.loc[(df.dir==1) & (df.sgn==-1), 'jv']

df.loc[(df.dir==2), 'irho'] = df.loc[(df.dir==2),'ii']
df.loc[(df.dir==2), 'jrho'] = df.loc[(df.dir==2),'jj']


# check if there are any bad rivers
# (i.e. pointing in the wrong direction or located on land)
ir = df.irho.to_numpy(dtype=int)
jr = df.jrho.to_numpy(dtype=int)
source_mask = m[jr,ir]
ngood = (source_mask==1).sum()
nbad = (source_mask==0).sum()
print('\n========================= Testing TRAPS direction  ========================\n')
print('Number of good TRAPS = %d' % (int(ngood)))
print('Number of bad TRAPS = %d' % (int(nbad)))

#################################################################################
#                                 Count sources                                 #
#################################################################################

# Initialize counters
LuvSrc_counts     = [0]*3
u_counts          = [0]*3
v_counts          = [0]*3
preLO_counts      = [0]*3
triv_counts       = [0]*3
mergedriv_counts  = [0]*3
LwSrc_counts      = [0]*3
mergedwwtp_counts = [0]*3
total_counts      = [0]*3

# Add expected values
LuvSrc_counts[0]    = 176
u_counts[0]         = 117
v_counts[0]         = 59
preLO_counts[0]     = 45
triv_counts[0]      = 131
mergedriv_counts[0] = 3
LwSrc_counts[0]     = 87
mergedwwtp_counts[0]= 0
total_counts[0]     = 263

# Add actual value
LuvSrc_counts[1] = df['dir'].value_counts()[0] + df['dir'].value_counts()[1]
u_counts[1] = df['dir'].value_counts()[0]
v_counts[1] = df['dir'].value_counts()[1]
preLO_counts[1] = df['name'].str.islower().sum()
triv_counts[1] = (~df['name'].str.islower() & ~df['dir'].__eq__(2)).sum()
mergedriv_counts[1] = ((df['name'].str.contains('\+') & df['dir'].__eq__(0)).sum() + 
                       (df['name'].str.contains('\+') & df['dir'].__eq__(1)).sum())
LwSrc_counts[1] = df['dir'].value_counts()[2]
mergedwwtp_counts[1] = (df['name'].str.contains('\+') & df['dir'].__eq__(2)).sum()
total_counts[1] = len(df.name.values.tolist())

# compare expected and actual values
counts = [LuvSrc_counts, u_counts, v_counts, preLO_counts, triv_counts,
         mergedriv_counts, LwSrc_counts, mergedwwtp_counts, total_counts]
for count in counts:
    if count[0] == count[1]:
        count[2] = 'Good'
    else:
        count[2] = 'Bad'

# create table
table = [LuvSrc_counts, u_counts, v_counts, preLO_counts, triv_counts,
         mergedriv_counts, LwSrc_counts, mergedwwtp_counts, total_counts]
srccount_df = pd.DataFrame(table, columns = ['Expectation', 'Actual', 'Result'],
                           index=['Total LuvSrc (Rivers)', '  u-src', '  v-src',
                                  '  pre-existing LO', '  tiny river', '  merged rivers',
                                  'Total LwSrc (WWTPs)', '  merged wwtps', 'Total Sources'])


print('\n======================== Testing number of sources ========================\n')
print(srccount_df)


#################################################################################
#      Get values in rivers.nc (based on LO/extract/river/extract_river.py)     #
#################################################################################

print('\nplotting data...')

# list of variables to extract
vn_list = ['transport', 'salt', 'temp','NO3',
           'NH4', 'Oxyg', 'TIC', 'TAlk']


dt0 = datetime.strptime(dstr0, Lfun.ds_fmt)
dt1 = datetime.strptime(dstr1, Lfun.ds_fmt)
ndays = (dt1-dt0).days + 1

# make mds_list: list of datestrings (e.g. 2021.01.01) to loop over
mds_list = []
mdt = dt0
while mdt <= dt1:
    mds_list.append(datetime.strftime(mdt, Lfun.ds_fmt))
    mdt = mdt + timedelta(days=1)

# get list of river names
# (this is a bit titchy because of NetCDF 3 limitations on strings, forcing them
# to be arrays of characters)
mds = mds_list[0]
fn = Ldir['LOo'] / 'forcing' / gridname / ('f' + mds) / frc / 'rivers.nc'
ds = xr.open_dataset(fn)
rn = ds['river_name'].values
NR = rn.shape[0]
riv_name_list = list(rn)

NT = len(mds_list)

# get state variable values
nanmat = np.nan * np.ones((NT, NR))
v_dict = dict()
for vn in vn_list:
    v_dict[vn] = nanmat.copy()
for tt,mds in enumerate(mds_list):
    this_dt = datetime.strptime(mds, Lfun.ds_fmt)
    if this_dt.day == 1 and this_dt.month == 1:
        print(' Year = %d' % (this_dt.year))
    fn = Ldir['LOo'] / 'forcing' / gridname / ('f' + mds) / frc / 'rivers.nc'
    ds = xr.open_dataset(fn)
    # The river transport is given at noon of a number of days surrounding the forcing date.
    # Here we find the index of the time for the day "mds".
    RT = pd.to_datetime(ds['river_time'].values)
    mdt = this_dt + timedelta(hours=12)
    mask = RT == mdt
    for vn in vn_list:
        if vn == 'transport':
            v_dict[vn][tt,:] = ds['river_' + vn][mask,:]
        else:
            # the rest of the variables allow for depth variation, but we
            # don't use this, so, just use the bottom value
            v_dict[vn][tt,:] = ds['river_' + vn][mask,0,:]
    ds.close()

# make transport positive
v_dict['transport'] = np.abs(v_dict['transport'])

# store output in an xarray Dataset
mdt_list = [(datetime.strptime(item, Lfun.ds_fmt) + timedelta(hours=12)) for item in mds_list]
times = pd.Index(mdt_list)

# create dataset
x = xr.Dataset(coords={'time': times,'riv': riv_name_list})

# add state variables
for vn in vn_list:
    v = v_dict[vn]
    x[vn] = (('time','riv'), v)

#################################################################################
#                               Prepare for plotting                            #
#################################################################################

# Prepare for plotting
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']
vn_list = ['transport', 'Oxyg', 'salt', 'temp',
           'NO3', 'NH4', 'TIC', 'TAlk']
var = ['Discharge [m3/s]', 'DO [mg/L]', 'Salt [psu]','Temp [C]',
       'NO3 [mmol/m3]', 'NH4 [mmol/m3]','TIC [mmol/m3]','TAlk [meq/m3]']
times = x.time.values

#################################################################################
#                              PLot TRAPS inputs                               #
#################################################################################

source_types = ['Pre-exsiting LO river', 'Tiny river', 'Point source']

for source_type in source_types:

    # get source specific information
    if source_type == 'Pre-exsiting LO river':
        # get subsets of dataset
        minindex = 0
        maxindex = preLO_counts[1]-1
        plot_name = 'LOrivInputs'
    elif source_type == 'Tiny river':
        # get subsets of dataset
        minindex = preLO_counts[1]
        maxindex = preLO_counts[1] + triv_counts[1]-1
        plot_name = 'TinyrivInputs'
    elif source_type == 'Point source':
        # get subsets of dataset
        minindex = total_counts[1] - LwSrc_counts[1]
        maxindex = total_counts[1]-1
        plot_name = 'PointSourceInputs'

    fig, axes = plt.subplots(4,2, figsize=(16, 9), sharex=True)
    ax = axes.ravel()
    for j,rn in enumerate(x.riv[minindex:maxindex].values):
        for i,vn in enumerate(vn_list):

            # get state variable
            vals = x[vn][:,x.riv == rn].values

            # convert DO from mmol/m3 to mg/L for plotting
            if vn == 'Oxyg':
                scale = 1/31.26 
            else:
                scale = 1

            # Plot values
            ax[i].plot(times,vals*scale,color='cornflowerblue',linewidth=2,alpha=0.3)

            # format plot
            if j == 0:
                # plot label
                ax[i].text(0.05,0.85,letters[i]+' '+var[i],transform=ax[i].transAxes,fontsize=14)
                # average value
                avg = scale*x[vn][:,minindex:maxindex].values.mean()
                ax[i].text(0.7,0.85,'Mean = {}'.format(round(avg,1)),transform=ax[i].transAxes,fontsize=14)
                # fontsize of tick labels
                ax[i].tick_params(axis='both', which='major', labelsize=12)
                ax[i].tick_params(axis='x', which='major', rotation=30)
                # axis limits
                if vn == 'salt':
                    ax[i].set_ylim([-0.25,0.25])
                else:
                    ax[i].set_ylim([0,scale*1.3*x[vn][:,minindex:maxindex].max()])

            # Define the date format
            if i >= 6:
                date_form = mdates.DateFormatter("%b")
                ax[i].xaxis.set_major_formatter(date_form)

        # add plot title
        if j == 0:
            plt.suptitle(source_type + ' inputs to ROMS ('+dstr0 +'-'+dstr1+')',fontsize=18)

    # Save figure
    save_path = out_dir / plot_name
    fig.savefig(save_path)
    plt.close('all')
