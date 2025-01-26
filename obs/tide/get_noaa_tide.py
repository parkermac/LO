"""
Code to automate getting year-long tide height records from
a series of NOAA sites around the Salish Sea and NE Pacific
coast.

Useful website describing the tidesandcurrents api:
https://api.tidesandcurrents.noaa.gov/api/prod/

Here is a sample tidesandcurrents api query for hourly tide data:
https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=20200101&end_date=20201231&station=8518750&product=hourly_height&datum=MLLW&time_zone=lst&units=metric&application=DataAPI_Sample&format=json

"""

import pandas as pd
import requests
from time import sleep
import sys

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# defaults for year(s), end date, and stations
year_list = [2022]
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
testing = False

if testing:
    noaa_sn_dict = {
        'Charleston': 9432780,
        'South Beach': 9435380}

out_dir = Ldir['LOo'] / 'obs' / 'tide'
Lfun.make_dir(out_dir)

# Start a DataFrame to hold station information. The index will be
# the station numbers "sn" such as 9432780, stored as strings.
sn_df = pd.DataFrame(columns=['name','lat','lon'])

# Extract and save data
for year in year_list:

    print('Working on year = ' + str(year))

    for name in noaa_sn_dict.keys():

        sn = noaa_sn_dict[name] # station number
        out_fn = out_dir / ('tide_noaa_' + str(sn) + '_' + str(year) + '.p') # data (pickled Series)
         
        print(' - NOAA station = ' + name)
        sys.stdout.flush()

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
        got_data = False
        while ntries < 3:
            try:
                response = requests.get(url_str, timeout = (3, 10))
                got_data = True
            except requests.exceptions.Timeout:
                print("The request timed out")
                got_data = False
            except requests.exceptions.RequestException as e:
                print("An error occurred: ", e)
                got_data = False
            if response.status_code != 200:
                print('Status code: ' + str(response.status_code))
                print('Reason: ' + str(response.reason))
                sys.stdout.flush()
                got_data = False

            if got_data == True: # get out of ntries loop
                print('    got_data = ' + str(got_data))
                break
            else:
                print('    got_data = ' + str(got_data))
                print('    ntries = ' + str(ntries))
            ntries += 1

        if got_data:
            # Save the data in a pickled pandas Series with the SSH time series
            json_dict = response.json()

            df0 = pd.json_normalize(json_dict['data']) # this makes a DataFrame of the data
            dti = pd.DatetimeIndex(df0.t)
            sn_ser = pd.Series(index=dti, data=df0.v.to_numpy(dtype=float))
            sn_ser.to_pickle(out_fn)

            # Add the metadata to a DataFrame
            sn_info = json_dict['metadata']
            # sn_info is a dict like:
            # {'id': '9432780', 'name': 'Charleston', 'lat': '43.3450', 'lon': '-124.3220'}
            name = sn_info['name'].split(',')[0].strip()
            lat = float(sn_info['lat'])
            lon = float(sn_info['lon'])
            sn_df.loc[str(sn),['name','lat','lon']] = name, lat, lon

            if testing:
                sn_ser.plot()
    
    # save the sn_df DataFrame, and a csv version of it for convenience.
    sn_df.to_pickle(out_dir / ('sn_df_noaa_' + str(year) + '.p'))
    sn_df.to_csv(out_dir / ('sn_df_noaa_' + str(year) + '.csv'))
    

    