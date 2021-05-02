"""
Functions to use with HYCOM extractions.

"""

import time
from urllib.request import urlretrieve
from urllib.error import URLError
from socket import timeout

# specify the various experiments we will use
hy_dict = {
    'hy1': ('GLBu0.08', '90.9'), # daily output
    'hy2': ('GLBu0.08', '91.0'), # daily output
    'hy3': ('GLBu0.08', '91.1'), # daily output
    'hy4': ('GLBu0.08', '91.2'), # daily output
    'hy5': ('GLBu0.08', '93.0'), # 3 hour output
    'hy6': ('GLBy0.08', '93.0'), # 3 hour output and finer grid
    }

# specify the sub region of hycom to extract
aa = [-131, -121, 39, 53]
    
def get_backfill_url(hy, dt, var_list):
    
    dstr0 = dt.strftime('%Y-%m-%d-T00:00:00Z')
    
    glb = hy_dict[hy][0]
    exnum = hy_dict[hy][1]
    
    # specify spatial limits
    north = aa[3]
    south = aa[2]
    west = aa[0] + 360
    east = aa[1] + 360
    
    url = ('http://ncss.hycom.org/thredds/ncss/' + glb + '/expt_' + exnum + 
        '?var='+var_list +
        '&north='+str(north)+'&south='+str(south)+'&west='+str(west)+'&east='+str(east) +
        '&disableProjSubset=on&horizStride=1' +
        '&time_start='+dstr0+'&time_end='+dstr0+'&timeStride=8' +
        '&vertCoord=&addLatLon=true&accept=netcdf4')
        
    return url
    
def get_extraction(hy, dt, out_fn, var_list):
    
    url = get_backfill_url(hy, dt, var_list)
    
    # get the data and save as a netcdf file
    counter = 1
    got_file = False
    while (counter <= 3) and (got_file == False):
        print('  Attempting to get data, counter = ' + str(counter))
        tt0 = time.time()
        try:
            (a,b) = urlretrieve(url, out_fn)
            # a is the output file name
            # b is a message you can see with b.as_string()
        except URLError as ee:
            if hasattr(ee, 'reason'):
                print('  *We failed to reach a server.')
                print('  -Reason: ', ee.reason)
            elif hasattr(ee, 'code'):
                print('  *The server could not fulfill the request.')
                print('  -Error code: ', ee.code)
        except timeout:
            print('  *Socket timed out')
        else:
            got_file = True
            print('  Worked fine')
        print('  took %0.1f seconds' % (time.time() - tt0))
        counter += 1
        
    if got_file:
        result = 'success'
    else:
        result = 'fail'
        
    return result
