# Code to extract mooring locations and print them to use in
# the mooring extraction code. 

# imports
from lo_tools import Lfun

import xarray as xr

Ldir = Lfun.Lstart()

# processed data location
source = 'ocnms'
otype = 'moor' 
data_dir = Ldir['LOo'] / 'obs' / source / otype

# named in this order for plotting 
sn_name_dict = {
    'MB042':'Makah Bay 42m',
    'MB015':'Makah Bay 15m',
    'CA042':'Cape Alava 42m',
    'CA015':'Cape Alava 15m',
    'TH042':'Teahwhit Head 42m',
    'TH015':'Teahwhit Head 15m',
    'KL027':'Kalaloch 27m',
    'KL015':'Kalaloch 15m',
    'CE042':'Cape Elizabeth 42m',
    'CE015':'Cape Elizabeth 15m'  
}

    # elif job_name == 'ooi':
    #     sta_dict = {
    #         'CE01':(-124.095, 44.6598), # Oregon Inshore (25 m)
    #         'CE02':(-124.304, 44.6393), # Oregon Shelf (80 m)
    #         'CE04':(-124.956, 44.3811), # Oregon Offshore (588 m)
    #         'PN01A':(-125.3983, 44.5096), # Slope Base (2905 m)
    #     }

sn_list = list(sn_name_dict.keys())

print("elif job_name == 'ocnms':")
print('    sta_dict = {')

for sn in sn_list:
    in_fn = data_dir / (sn + '_2011_2023_hourly.nc')
    
    ds = xr.open_dataset(in_fn, decode_times=True)

    print("        \'%s\':(%0.4f, %0.4f), # %s" % (sn, ds.lon, ds.lat, ds.attrs['Station Name']))

print('    }')