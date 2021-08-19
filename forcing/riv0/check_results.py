"""
This is to compare results of LO/forcing/riv0 for a daily forecast
to the same file from LiveOcean/forcing/riv2

RESULT: the two files appear to be quite similar.  The names and times
are identical.  The tranports are very close.  The temperatures are
sometimes colder by as much as 2 degC in the new rivers.  Presumably this
is because I have re-made the climatologies and there was data from
more rivers this time.

"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

from lo_tools import Lfun
Ldir = Lfun.Lstart(gridname='cas6',tag='v3')
Ldir['date_string'] = '2021.04.20'

old_fn = Ldir['parent'] / 'LiveOcean_output' / Ldir['gtag'] / ('f' + Ldir['date_string']) / 'riv2' / 'rivers.nc'
new_fn = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / 'riv0' / 'rivers.nc'

# get old river name list
old = nc.Dataset(old_fn)
orn = old['river_name'][:]
NR = orn.shape[1]
orn_list = []
for ii in range(NR):
    a = orn[:,ii]
    if isinstance(a, np.ma.MaskedArray):
        a = a.data
    r = []
    for l in a:
        r.append(l.decode())
    rr = ''.join(r)
    orn_list.append(rr)
    
# get new river name list
new = nc.Dataset(new_fn)
nrn_list = list(new['river_name'][:])

# check that names match
both = dict(zip(orn_list, nrn_list))
they_match = True
for k in both.keys():
    if k != both[k]:
        print(k + ' = unmatched river name')
        
# check that times match
ot = old['river_time'][:].data
nt = new['river_time'][:].data
if ((ot - nt) != 0).any():
    print('times different')
    
# compare flow at t=0
oQ = old['river_transport'][0,:].data
nQ = new['river_transport'][0,:].data

# compare temperature at t=0
oT = old['river_temp'][0,:].data
nT = new['river_temp'][0,:].data

# plotting
plt.close('all')
fig = plt.figure(figsize=(16,8))

ax = fig.add_subplot(121)
ax.plot(oQ,nQ, 'og')
ax.axis('square')
aa = ax.axis()
ax.plot(aa[:2],aa[:2], '-b')
ax.grid(True)
ax.set_title('river_transport')
ax.set_xlabel('Old')
ax.set_ylabel('New')

ax = fig.add_subplot(122)
ax.plot(oT,nT, 'og')
#ax.axis('square')
aa = ax.axis([5,15,5,15])
ax.plot([5,15], [5,15], '-b')
ax.grid(True)
ax.set_title('river_temp')
ax.set_xlabel('Old')
ax.set_ylabel('New')

plt.show()