"""
Function to do the extraction of all sections for a single history file.
"""

from argparse import ArgumentParser
from xarray import open_dataset
from numpy import nan, ones, diff
import pickle
from pandas import read_pickle

parser = ArgumentParser()
parser.add_argument('-sect_df_fn', type=str) # path to sect_df
parser.add_argument('-in_fn', type=str) # path to history file
parser.add_argument('-out_fn', type=str) # path to outfile (temp directory)
parser.add_argument('-test', type=str) # True or False
args = parser.parse_args()

sect_df = read_pickle(args.sect_df_fn)

vn_list = ['salt']

# grid info
ds = open_dataset(args.in_fn)
DX = 1/ds.pm.values
DY = 1/ds.pn.values
# Get spacing on u and v grids
dxv = DX[:-1,:] + diff(DX,axis=0)/2 # DX on the v-grid
dyu = DY[:,:-1] + diff(DY,axis=1)/2 # DY on the u-grid
# separate out u and v parts of sect_df
u_df = sect_df[sect_df.uv == 'u']
v_df = sect_df[sect_df.uv == 'v']

# Fields that do not change with time
C = dict() # this is for holding raw model fields
CC = dict() # this is for holding fields extracted on sections

# get depth at section points
h = ds.h.values
CC['h'] = (h[sect_df.jrp, sect_df.irp]  + h[sect_df.jrm, sect_df.irm])/2
# note that we are interpolating from two rho-grid points onto a u- or v-grid
# point

# get width at section points
dxvv = dxv[v_df.j, v_df.i]
dyuu = dyu[u_df.j, u_df.i]
dd = nan * ones(CC['h'].shape)
dd[v_df.index] = dxvv
dd[u_df.index] = dyuu
CC['dd'] = dd

# First: tracers and zeta
for vn in vn_list:
    C[vn] = ds[vn].values.squeeze()
C['zeta'] = ds.zeta.values.squeeze()
for vn in vn_list:
    CC[vn] = (C[vn][:, sect_df.jrp, sect_df.irp]  + C[vn][:, sect_df.jrm, sect_df.irm])/2
CC['zeta'] = (C['zeta'][sect_df.jrp, sect_df.irp]  + C['zeta'][sect_df.jrm, sect_df.irm])/2

# Then: velocity
u = ds.u.values.squeeze()
v = ds.v.values.squeeze()
# the "-1" in the reshape index below means "figure it out based on the context"
uu = u[:, u_df.j, u_df.i] * u_df.pm.to_numpy().reshape(1,-1)
vv = v[:, v_df.j, v_df.i] * v_df.pm.to_numpy().reshape(1,-1)
# then merge these back into one
vel = nan * ones(CC['salt'].shape)
# I love fancy indexing!
vel[:,u_df.index] = uu
vel[:,v_df.index] = vv

CC['vel'] = vel

# print(CC['vel'].shape)

pickle.dump(CC, open(args.out_fn,'wb'))