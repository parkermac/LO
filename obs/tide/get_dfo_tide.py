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
import dfo_info_list as dil
from importlib import reload
reload(dil)

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
testing = True

if testing:
    dfo_sn_dict = {'Point Atkinson': '07795'}
    month_list = [1]

out_dir = Ldir['LOo'] / 'obs' / 'tide'
Lfun.make_dir(out_dir)

# Extract and save data
for year in year_list:

    print('Working on year = ' + str(year))

    for name in dfo_sn_dict.keys():

        sn = dfo_info_dict[name]['code'] # station number
        id = dfo_info_dict[name]['id']

        out_fn = out_dir / ('tide_' + str(sn) + '_' + str(year) + '.p') # data (pickled Series)
        metadata_out_fn = out_dir / ('info_' + str(sn) + '_' + str(year) + '.csv') # station metadata
        
        print(' - DFO station = ' + name)
        sys.stdout.flush()

        # debugging before we hide this in a function

        sn = str(sn) # station number
        year_str = str(year)
        url_str = ('https://api.iwls-sine.azure.cloud-nuage.dfo-mpo.gc.ca/api/v1/stations/'
            + id + '/data?time-series-code=wlo'
            + '&from=' + year_str + '-01-01T00%3A00%3A00Z'
            + '&to=' + year_str + '-01-31T00%3A00%3A00Z'
            + '&station=' + sn
            + '&resolution=SIXTY_MINUTES')

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
            if response.status_code != 200:
                print('Status code: ' + str(response.status_code))
                print('Reason: ' + str(response.reason))
                sys.stdout.flush()
                got_data = False
            if got_data == True: # get out of ntries loop
                break
            else:
                sleep(3)
                print('ntries = ' + str(ntries))
                ntries += 1

        if got_data:
            # save the data in
            # (i) a csv with the metadata
            # (ii) a pickled pandas Series with the SSH time series
            json_list = response.json()

            # sn_info = json_dict['metadata']
            # # sn_info is a dict like:
            # # {'id': '9432780', 'name': 'Charleston', 'lat': '43.3450', 'lon': '-124.3220'}

            df0 = pd.json_normalize(json_list) # this makes a DataFrame
            dti = pd.DatetimeIndex(df0.eventDate)
            sn_ser = pd.Series(index=dti, data=df0.value.to_numpy(dtype=float))
            sn_ser.to_pickle(out_fn)
            # Lfun.dict_to_csv(sn_info, metadata_out_fn)
            if testing:
                sn_ser.plot()

    