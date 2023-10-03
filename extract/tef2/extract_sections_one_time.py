"""
Function to do the extraction of all sections for a single history file.
"""

from argparse import ArgumentParser
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from pickle import dump
from pandas import read_pickle, to_datetime

parser = ArgumentParser()
parser.add_argument('-sect_df_fn', type=str) # path to sect_df
parser.add_argument('-in_fn', type=str) # path to history file
parser.add_argument('-out_fn', type=str) # path to outfile (temp directory)
parser.add_argument('-vn_type', type=str) # 'salt' or 'bio'
args = parser.parse_args()

sect_df = read_pickle(args.sect_df_fn)
ds = open_dataset(args.in_fn, decode_times=False)
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

if args.vn_type == 'salt':
    vn_list = ['salt']
elif args.vn_type == 'bio':
    if 'NH4' in ds.data_vars:
        vn_list = ['salt', 'temp', 'oxygen',
            'NO3', 'NH4', 'phytoplankton', 'zooplankton', 'SdetritusN', 'LdetritusN',
            'TIC', 'alkalinity']
    else:    
        vn_list = ['salt', 'temp', 'oxygen',
            'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
            'TIC', 'alkalinity']
# grid info
DX = 1/ds.pm.values
DY = 1/ds.pn.values
# Get spacing on u and v grids
dxv = DX[:-1,:] + diff(DX,axis=0)/2 # DX on the v-grid
dyu = DY[:,:-1] + diff(DY,axis=1)/2 # DY on the u-grid
# separate out u and v parts of sect_df
u_df = sect_df[sect_df.uv == 'u']
v_df = sect_df[sect_df.uv == 'v']

# Fields that do not change with time
#
CC = dict() # this is for holding fields extracted on sections
# get depth at section points
h = ds.h.values
CC['h'] = (h[sect_df.jrp, sect_df.irp]  + h[sect_df.jrm, sect_df.irm])/2
# note that we are interpolating from two rho-grid points onto the u- or v-grid
# get width at section points
dxvv = dxv[v_df.j, v_df.i]
dyuu = dyu[u_df.j, u_df.i]
dd = nan * ones(CC['h'].shape)
dd[v_df.index] = dxvv
dd[u_df.index] = dyuu
CC['dd'] = dd

# Fields that do change with time
#
# First: tracers and zeta
for vn in vn_list:
    aa = ds[vn].values.squeeze()
    CC[vn] = (aa[:, sect_df.jrp, sect_df.irp]  + aa[:, sect_df.jrm, sect_df.irm])/2
aa = ds.zeta.values.squeeze()
CC['zeta'] = (aa[sect_df.jrp, sect_df.irp]  + aa[sect_df.jrm, sect_df.irm])/2
# Then: velocity
u = ds.u.values.squeeze()
v = ds.v.values.squeeze()
# the "-1" in the reshape index below means "figure it out based on the context"
uu = u[:, u_df.j, u_df.i] * u_df.pm.to_numpy().reshape(1,-1)
vv = v[:, v_df.j, v_df.i] * v_df.pm.to_numpy().reshape(1,-1)
# merge u and v parts back into one
vel = nan * ones(CC['salt'].shape)
# I love fancy indexing!
vel[:,u_df.index] = uu
vel[:,v_df.index] = vv
CC['vel'] = vel

# put these in a Dataset
NZ, NP = CC['vel'].shape
ot = ds.ocean_time.values
attrs = {'units':ds.ocean_time.units}
ds1 = Dataset()
ds1['time'] = (('time'), ot, attrs)
# alternate code to do the same thing:
# ds1 = Dataset({'time': ('time', ot, attrs)})

ds1['h'] = (('p'), CC['h'])
ds1['dd'] = (('p'), CC['dd'])
ds1['zeta'] = (('time','p'), CC['zeta'].reshape(1,NP))
for vn in CC.keys():
    if vn not in ['zeta', 'h', 'dd']:
        vv = CC[vn] # packed (z,p)
        ds1[vn] = (('time','z', 'p'), vv.reshape(1,NZ,NP))
ds1.to_netcdf(args.out_fn, unlimited_dims='time')

