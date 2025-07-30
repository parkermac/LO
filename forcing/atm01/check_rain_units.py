"""
Code to explore rain units in WRF.
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path

d_str = '2019070400'
in_dir = Path('/Users/pm8/Documents/LO_data/wrf') / d_str

NT = 73
d2_list = []
for hr in range(NT):
    hr_str = ('0' + str(hr))[-2:]
    d2_list.append(in_dir / ('wrfout.ocean_d2.' + d_str + '.f' + hr_str + '.0000'))
    
rain_acc = []
for d2 in d2_list:
    ds = xr.open_dataset(d2)
    rain_acc.append(ds.RAINC[0, 10, 10].values + ds.RAINNC[0, 10, 10].values)

rain_acc = np.array(rain_acc)
hour = np.arange(NT)

# convert to [kg m-2 s-1]
rain = np.zeros(NT)
rain[:-1] = (rain_acc[1:] - rain_acc[:-1]) / 3600
rain[-1] = rain[-2]

plt.close('all')
fig = plt.figure(figsize=(14,7))

ax = fig.add_subplot(121)
ax.plot(hour, rain_acc, '-b')
ax.grid(True)
ax.set_title('Accumulated rain [mm]')
ax.set_xlabel('Time [hours]')

ax = fig.add_subplot(122)
ax.plot(hour, rain, '-r')
ax.grid(True)
ax.set_title('Rain rate [kg m-2 s-1]')
ax.set_xlabel('Time [hours]')

plt.show()
