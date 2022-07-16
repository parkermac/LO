"""
This is the main program for copying the daily forecast to Azure blob storage.
It copies all three days of the forecast, overwriting existing days 1 and 2 from the
previous forecast.

I created this as a standalone program because it is a very niche task.  It is designed as
a post job that runs with the forecast, but it can also run in backfill.

Test on mac or apogee:
run post_main.py -gtx cas6_v0_live -ro 0 -d 2019.07.04 -r backfill -job azure0 -test True

Run for real on apogee:
python post_main.py -gtx cas6_v0_u0kb -ro 0 -d [today's datestring] -r forecast -job azure0 < /dev/null > azure0.log &

Or just put it in the lineup of the driver_post jobs.

Test output on my mac:
====================== copy_to_azure =======================
2022.07.16 11:23:47
 -- Copying: ocean-his-0001.nc to cas6-v0-live-f2019-07-04
 -- URL to access file: https://pm2.blob.core.windows.net/cas6-v0-live-f2019-07-04/ocean-his-0001.nc
 -- Took 183.39 sec
"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

from time import time
from lo_tools import Lfun
from azure.storage.blob import BlobServiceClient
from azure.storage.blob import PublicAccess

ds0 = Ldir['date_string']
dt0 = datetime.strptime(ds0, Ldir['ds_fmt'])
if Ldir['run_type'] == 'backfill':
    dt1 = dt0
elif Ldir['run_type'] == 'forecast':
    ndays = Ldir['forecast_days']
    dt1 = dt0 + timedelta(days=ndays-1)

# screen output
print(' copy_to_azure '.center(60,'='))
print(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))

result = 'success'
this_dt = dt0
while this_dt <= dt1:
    this_ds = this_dt.strftime(Ldir['ds_fmt'])
    f_string = 'f' + this_ds
    fn_list = Lfun.get_fn_list('hourly', Ldir, this_ds, this_ds)

    if Ldir['testing'] == True:
        fn_list = [fn_list[0]]
        
    for fn in fn_list:
        # AZURE
        out_name = fn.name.replace('_','-')
        container_name = Ldir['gtagex'].replace('_','-') + '-' + f_string.replace('.','-')
        print(' -- Copying: %s to %s' % (out_name, container_name))
        
        tt0 = time()
        # get the blob_service
        account_fn = Ldir['data'] / 'accounts' / 'azure_pm_2015.05.25.txt'
        with account_fn.open() as aa:
            # make sure to strip newline characters
            account = aa.readline().strip()
            key = aa.readline().strip()
        connection_string = ('DefaultEndpointsProtocol=https' +
            ';AccountName=' + account +
            ';AccountKey=' + key + 
            ';EndpointSuffix=core.windows.net')
        blob_service = BlobServiceClient.from_connection_string(conn_str=connection_string)
        try: # create the container if needed
            blob_service.create_container(container_name, public_access=PublicAccess.Container)
        except:
            pass # assume error is because container exists
        # write it to Azure
        try:
            from azure.storage.blob import BlobClient
            blob = BlobClient.from_connection_string(conn_str=connection_string,
                container_name=container_name, blob_name=out_name)
            with fn.open('rb') as data:
                blob.upload_blob(data, overwrite=True)
            az_url = ('https://pm2.blob.core.windows.net/' + container_name + '/' + out_name)
            print(' -- URL to access file: %s' % (az_url))
            print(' -- Took %0.2f sec' % (time()-tt0))
        except Exception as e:
            print(' ** Copy failed: %s' % (str(e)))
            result = 'fail'

    this_dt += timedelta(days=1)

# -------------------------------------------------------

# test for success
result_dict['result'] = result

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
