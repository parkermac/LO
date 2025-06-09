"""
Make a plot of the Ecology Stations, for reference.
"""

from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
import matplotlib.pyplot as plt
import pandas as pd

Ldir = Lfun.Lstart()

plt.close('all')
pfun.start_plot(figsize=(10,13))

in_dir0 = Ldir['LOo'] / 'obs' / 'ecology_nc'
in_dir = in_dir0 / 'ctd'
in_fn = in_dir / 'info_2019.p' # choosing one I assume is representative
out_fn = in_dir0 / 'station_map.png'

df = pd.read_pickle(in_fn)
sta_list = df.name.unique()

sdf = pd.DataFrame(index=sta_list, columns=['lon','lat'])

# pull out just a list of station names and locations, for plotting
for sta in sta_list:
    df1 = df[df.name==sta]
    s1 = df1.iloc[0,:]
    sdf.loc[sta,'lon'] = s1.lon
    sdf.loc[sta,'lat'] = s1.lat

fig = plt.figure()
ax = fig.add_subplot(111)
pfun.add_coast(ax,color='g')
for sta in sdf.index:
    lon = sdf.loc[sta,'lon']
    lat = sdf.loc[sta,'lat']
    ax.plot(lon,lat,'or',markersize=5)
    if sta in ['WPA113','HCB004']:
        ax.text(lon,lat,sta,rotation=0,fontsize=10)
    else:
        ax.text(lon,lat,sta,rotation=30,fontsize=10)

pfun.dar(ax)
ax.axis([-124.5,-122,46.25,49])

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('WA Ecology Station Map')
fig.tight_layout()

plt.show()

fig.savefig(out_fn)