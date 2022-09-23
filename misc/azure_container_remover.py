"""
Code to delete a range of containers.

This code has to be customized by hand to delete specific containers.

RESULT: this is fast and effective for things that are tedious to do from the portal.
"""

from azure.storage.blob import BlobServiceClient
from azure.storage.blob import PublicAccess
from lo_tools import Lfun
from time import time
from datetime import datetime, timedelta

Ldir = Lfun.Lstart()

# screen output
print(' azure_container_remover.py '.center(60,'='))
print(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
print('Deleting user-specified containers')

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

dt0 = datetime(2020,9,15)
# dt1 = datetime(2020,9,16)
dt1 = datetime.now()

dt = dt0
while dt <= dt1:
    container_name = 'f' + dt.strftime('%Y%m%d')
    print(container_name)
    
    try:
        result = blob_service.delete_container(container_name)
    except Exception as e:
        print(' no action: assume container did not exist')
        # print(e)
    
    dt = dt + timedelta(days=1)

print('DONE')
print('Total time %0.2f sec' % (time()-tt0))

