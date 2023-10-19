"""
This code pulls out useful information about the Collias dataset.
"""

import pandas as pd
import numpy as np

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.options.display.width = 0 # auto-detect full display width

source = 'collias'
in_dir = Ldir['data'] / 'obs' / source

fn = in_dir / 'EIMDiscreteResults_2023May17_456416.csv'
print('\nFile Name:')
print(fn)
df = pd.read_csv(fn, low_memory=False,
    parse_dates=['Field_Collection_Start_Date_Time'])
print('\nColumn Names:')
print(df.columns)

# Check on data and units
print('\nData Names and Units:')
rpn_set = set(df['Result_Parameter_Name'].to_list())
for rpn in rpn_set:
    unit_set = set(df.loc[df['Result_Parameter_Name']==rpn,
        'Result_Value_Units'].to_list())
    for unit in unit_set:
        print('%30s [%s]' % (rpn,unit))

# Rename columns, for convenience
col_dict = {
    'Location_ID':'name',
    'Field_Collection_Start_Date_Time':'time',
    'Field_Collection_Upper_Depth':'depth',
    'Result_Parameter_Name':'vn',
    'Result_Value':'val',
    'Result_Value_Units':'units',
    'Calculated_Latitude_Decimal_Degrees_NAD83HARN':'lat',
    'Calculated_Longitude_Decimal_Degrees_NAD83HARN':'lon'
}
df = df.rename(columns=col_dict)

# Just look at one variable, to explore space and time extents
vn = 'Temperature, water'
this_df = df.loc[df.vn==vn,:]

lon = this_df.lon.to_numpy()
lat = this_df.lat.to_numpy()
z = -this_df.depth.to_numpy()
t = this_df.time

# Plotting
import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun
plt.close('all')
pfun.start_plot(figsize=(12,12))

fig = plt.figure()

ax = fig.add_subplot(211)
ax.plot(lon, lat,'.b')
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_title('Collias Dataset')
ax.set_xlabel('Longitude [deg]')
ax.set_ylabel('Latitude [deg]')

ax.axis([-125, -122, 47, 49])
ax = fig.add_subplot(212)
ax.plot(t,z,'.b')
ax.set_xlabel('Time')
ax.set_ylabel('Z [m]')
ax.grid(True)

plt.show()
plt.savefig(in_dir / 'info_plot.png')
pfun.end_plot()
