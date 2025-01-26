"""
Code to automate getting year-long tide height records from
a series of DFO sites around the Salish Sea and NE Pacific
coast.

This works, although the server response is not perfect and there are a few
gaps in the data, which are stored as nan's in the output series. Tedious
because we can only ask for a month of data at a time.

Takes about 10 minutes for a year.
"""

import pandas as pd
import requests
from time import sleep
import sys
from datetime import datetime
import numpy as np

from lo_tools import Lfun
Ldir = Lfun.Lstart()
import dfo_info_list as dil
import matplotlib.pyplot as plt

# defaults for year(s), end date, and stations
year_list = [2022]
month_list = list(range(1,13))
dfo_sn_dict = {
    'Point Atkinson': '07795',
    'Vancouver': '07735',
    'Patricia Bay': '07277',
    'Victoria Harbour': '07120',
    'Bamfield': '08545',
    'Tofino': '08615',
    'Campbell River': '08074',
    'New Westminster': '07654'}

# ingest station info
dfo_info_dict = dict()
for item in dil.dfo_info_list:
    dfo_info_dict[item["officialName"]] = item

# flag for testing
testing = False

if testing:
    dfo_sn_dict = {
        'Point Atkinson': '07795',
        'Campbell River': '08074'}
    month_list = [1,2]

out_dir = Ldir['LOo'] / 'obs' / 'tide'
Lfun.make_dir(out_dir)

# Extract and save data
for year in year_list:

    # Start a DataFrame to hold station information. The index will be
    # the station numbers "sn" such as 07795, stored as strings.
    sn_df = pd.DataFrame(columns=['name','lat','lon'])

    print('Working on year = ' + str(year))

    plt.close('all')

    for name in dfo_sn_dict.keys():

        sn = dfo_info_dict[name]['code'] # station number
        id = dfo_info_dict[name]['id'] # The magic number, like '5ddeee612a8a340001a3676f' that
            # is needed to get the data for this station.

        out_fn = out_dir / ('tide_dfo_' + str(sn) + '_' + str(year) + '.p') # data (pickled Series)
        
        print(' - DFO station = ' + name)
        sys.stdout.flush()

        sn = str(sn) # station number
        year_str = str(year)

        # Initialize a Series to hold the data.
        t = pd.Timestamp(datetime(year,month_list[-1],1))
        ndays = t.days_in_month
        dti = pd.date_range(datetime(year,1,1,0,0,0),datetime(year,month_list[-1],ndays,23,0,0),freq='h')
        sn_ser = pd.Series(index=dti, data=np.nan)
        # We do this because sometimes there are data gaps, so we cannot just concatenate months.

        for month in month_list:
            mo_str = ('00' + str(month))[-2:]
            print('  -- month = ' + mo_str)
            t = pd.Timestamp(datetime(year,month,1))
            ndays = str(t.days_in_month)

            url_str = ('https://api.iwls-sine.azure.cloud-nuage.dfo-mpo.gc.ca/api/v1/stations/'
                + id + '/data?time-series-code=wlo'
                + '&from=' + year_str + '-' + mo_str + '-01T00:00:00Z'
                + '&to=' + year_str + '-' + mo_str + '-' + ndays + 'T23:00:00Z'
                + '&station=' + sn
                + '&resolution=SIXTY_MINUTES')

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
                # save the data in
                # (i) a csv with the metadata
                # (ii) a pickled pandas Series with the SSH time series
                json_list = response.json()

                df0 = pd.json_normalize(json_list) # this makes a DataFrame
                dti = pd.DatetimeIndex(df0.eventDate).tz_localize(None)
                data=df0.value.to_numpy(dtype=float)

                sn_ser[dti] = data

                if month == 1:
                    # save metadata
                    lat = float(dfo_info_dict[name]['latitude'])
                    lon = float(dfo_info_dict[name]['longitude'])
                    sn_df.loc[str(sn),['name','lat','lon']] = name, lat, lon
                else:
                    pass

        sn_ser.to_pickle(out_fn)

        if testing:
            fig = plt.figure()
            sn_ser.plot(fig=fig,title=name)

    # save the sn_df DataFrame, and a csv version of it for convenience.
    sn_df.to_pickle(out_dir / ('sn_df_dfo_' + str(year) + '.p'))
    sn_df.to_csv(out_dir / ('sn_df_dfo_' + str(year) + '.csv'))


    