"""
Code to test, graphically, the results of the new make_forcing_main.py
which writes only selected surface fields to NetCDF and Azure.

"""

# setup
import os, sys
sys.path.append(os.path.abspath('../../alpha'))
import Lfun
import netCDF4 as nc
sys.path.append(os.path.abspath('../../plotting'))
import pfun
import matplotlib.pyplot as plt

Ldir = Lfun.Lstart('cas6', 'v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'
f_string = 'f2019.07.04'
#f_string = 'f2020.07.14'
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'

# Note: the file name and variable list should relate to those in make_forcing_main.py.
if False:
    in_name = 'ocean_surface.nc'
    vn_list = ['Uwind', 'Vwind','salt', 'temp','u','v']
else:
    in_name = 'ocean_layers.nc'
    testing = False
    if testing:
        tag_list = ['surface','bottom','10']
        vn_list = ['oxygen_'+tag for tag in tag_list]
    else:
        tag_list = ['surface','bottom','10','20']
        vn_list = ['oxygen_'+tag for tag in tag_list] + ['temp_'+tag for tag in tag_list]
        

# PLOTTING
ds = nc.Dataset(in_dir + in_name)
ot = ds['ocean_time'][:]
lon = ds['lon_psi'][:]
lat = ds['lat_psi'][:]
plt.close('all')
fs=14
plt.rc('font', size=fs)
for vn in vn_list:
    # make one plot for each variable, with two panels: start and end time
    fig = plt.figure(figsize=(16,10))
    nplot = 1
    for tlev in [0, -1]:
        ax = fig.add_subplot(1,2,nplot)
        cs = ax.pcolormesh(lon, lat, ds[vn][tlev, 1:-1, 1:-1],
                           cmap='rainbow')
        try:
            tun = ds[vn].units
        except AttributeError:
            tun = ''
        ax.set_title(vn + ' (' + tun + ')')
        fig.colorbar(cs)
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        dt = Lfun.modtime_to_datetime(ot[tlev])
        ax.text(.98, .075, dt.strftime('%Y-%m-%d'),
            ha='right' , va='bottom', transform=ax.transAxes)
        ax.text(.98, .065, dt.strftime('%H:%M') + ' UTC',
            ha='right', va='top', transform=ax.transAxes)
        nplot += 1
plt.show()
plt.rcdefaults()

ds.close()
