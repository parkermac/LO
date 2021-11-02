"""
Standalone program to copy a file to azure, and return the URL to access the file.
The idea is that you can run from the command line and escape to shell because
sometimes the write takes awhile.  Then look in the log file to get the URL for
sharing.

To test:
run copy_to_azure

To copy a file:
python copy_to_azure.py -fn [string of full path] > az.log &

Note that the last line of screen output from extract_box gives you
the [string of full path].

"""

from azure.storage.blob import BlobServiceClient
from azure.storage.blob import PublicAccess
from lo_tools import Lfun
import argparse
from time import time
from datetime import datetime
from pathlib import Path

# get and process arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fn', type=str) # full path to input file
parser.add_argument('-out_name', type=str) # name of output file
parser.add_argument('-container_name', type=str, default='pm-share')
args = parser.parse_args()
fn = args.fn
out_name = args.out_name
container_name = args.container_name
Ldir = Lfun.Lstart()
if fn == None:
    fn = Ldir['data'] / 'accounts' / 'az_testfile.txt'
else:
    fn = Path(fn)
if out_name == None:
    out_name = fn.name

# screen output
print(' copy_to_azure '.center(60,'='))
print(datetime.now().strftime('%Y.%m.%d %H:%M:%S'))
print('Copying: %s' % (fn))

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
    print('URL to access file: %s' % (az_url))
    print('Took %0.2f sec' % (time()-tt0))
except Exception as e:
    print('Copy failed: %s' % (str(e)))
