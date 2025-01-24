"""
Code to automate getting year-long tide height records from
a series of NOAA sites around the Salish Sea and NE Pacific
coast.
"""

import pandas as pd
import requests
from time import sleep
import sys

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# defaults for year(s), end date, and stations
year_list = [2017]
end_date = '1231'

noaa_sn_dict = {
    'Charleston': 9432780,
    'South Beach': 9435380,
    'Garibaldi': 9437540,
    'Toke Point': 9440910,
    'Westport': 9441102,
    'La Push': 9442396,
    'Neah Bay': 9443090,
    'Port Angeles': 9444090,
    'Friday Harbor': 9449880,
    'Cherry Point': 9449424,
    'Port Townsend': 9444900,
    'Seattle': 9447130,
    'Tacoma': 9446484}

# flag for testing
testing = True

if testing:
    noaa_sn_dict = {'Charleston': 9432780}

out_dir = Ldir['LOo'] / 'obs' / 'tide'
Lfun.make_dir(out_dir)

# Extract and save data
for year in year_list:

    print('Working on year = ' + str(year))

    for name in noaa_sn_dict.keys():

        sn = noaa_sn_dict[name] # station number
        out_fn = out_dir / ('tide_' + str(sn) + '_' + str(year) + '.p') # data (pickled Series)
        metadata_out_fn = out_dir / ('info_' + str(sn) + '_' + str(year) + '.csv') # station metadata
        # harmonics_out_fn = out_dir / ('h_' + str(sn) + '_' + str(year) + '.p') # harmonics (do later)
        
        print(' - NOAA station = ' + name)
        sys.stdout.flush()

        # debugging before we hide this in a function

        sn = str(sn) # station number
        year_str = str(year)
        url_str = ('https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?'
            + 'begin_date=' + year_str + '0101'
            + '&end_date=' + year_str + end_date
            + '&station=' + sn
            + '&product=hourly_height'
            + '&datum=MLLW&units=metric&time_zone=gmt'
            + '&application=University_of_Washington_LO'
            + '&format=json')

        ntries = 0
        while ntries < 3:
            got_data = True
            try:
                response = requests.get(url_str, timeout = (3, 10))
            except requests.exceptions.Timeout:
                print("The request timed out")
                got_data = False
            except requests.exceptions.RequestException as e:
                print("An error occurred: ", e)
                got_data = False
            if got_data == True: # get out of ntries loop
                break
            else:
                if response.status_code != 200:
                    print('Status code: ' + str(response.status_code))
                    print('Reason: ' + str(response.reason))
                    sys.stdout.flush()
                got_data = False
                sleep(3)
                print('ntries = ' + str(ntries))
                ntries += 1

        if got_data:
            # save the data in
            # (i) a csv with the metadata
            # (ii) a pickled pandas Series with the SSH time series
            json_dict = response.json()
            sn_info = json_dict['metadata']
            # sn_info is a dict like:
            # {'id': '9432780', 'name': 'Charleston', 'lat': '43.3450', 'lon': '-124.3220'}
            df0 = pd.json_normalize(json_dict['data']) # this makes a DataFrame
            dti = pd.DatetimeIndex(df0.t)
            sn_ser = pd.Series(index=dti, data=df0.v.to_numpy(dtype=float))
            sn_ser.to_pickle(out_fn)
            Lfun.dict_to_csv(sn_info, metadata_out_fn)
            if testing:
                sn_ser.plot()

    