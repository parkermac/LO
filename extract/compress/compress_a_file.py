import sys
from time import time
import xarray as xr

arg_list = sys.argv
if len(arg_list) != 2:
    print('Error: compress_a_file needs just one argument')
    sys.exit()
    
fn = arg_list[1]

# compress that copy in place
tt0 = time()
ds = xr.load_dataset(fn)
enc_dict = {'zlib':True, 'complevel':1}
Enc_dict = {vn:enc_dict for vn in ds.data_vars if 'ocean_time' in ds[vn].dims}
ds.to_netcdf(fn, encoding=Enc_dict)
ds.close()
print('- %0.1f sec to compress %s/%s' % (time()-tt0, fn.split('/')[-2], fn.split('/')[-1]))
