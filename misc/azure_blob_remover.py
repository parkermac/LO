#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to delete selected blobs or containers
"""

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

from datetime import datetime, timedelta

# ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime
start_time = datetime.now()

dt0 = datetime(2019,7,1) # start time
dt1 = datetime(2020,9,1) # end time

# # account name and key
# azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
# account = azu_dict['account']
# key = azu_dict['key']
# # get a handle to the account
# blob_service = BlockBlobService(account_name=account, account_key=key)

# get the blob_service
from azure.storage.blob import BlobServiceClient
from azure.storage.blob import PublicAccess
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']
connection_string = ('DefaultEndpointsProtocol=https' +
    ';AccountName=' + account +
    ';AccountKey=' + key + 
    ';EndpointSuffix=core.windows.net')
blob_service = BlobServiceClient.from_connection_string(conn_str=connection_string)



dt = dt0
while dt <= dt1:
    Ldir['date_string'] = dt.strftime('%Y.%m.%d')
    print('Deleting Azure files for ' + Ldir['date_string'])
    f_string = 'f' + Ldir['date_string']
    ff_string = f_string.replace('.','') # azure does not like dots in container names
    container_name = ff_string
    
    # create the container if needed
    try:
        blob_service.create_container(container_name, public_access=PublicAccess.Container)
        
    except:
        # assume error is because container exists
        pass
    result = blob_service.delete_container(container_name)
    print(' - result = ' + str(result))
    dt = dt + timedelta(days=1)
    
print('DONE')
