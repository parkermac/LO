"""
This is the code for doing cast extractions.

Test on mac in ipython:

run extract_casts.py -g cas6 -t v3 -x lo8b -test True

"""

from pathlib import Path
import sys
from datetime import datetime

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun

Ldir = exfun.intro() # this handles the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import pandas as pd
import zfun
import zrfun
import subprocess
import netCDF4 as nc
import Lfun
import numpy as np

if Ldir['a1'] == 'woac2':
    pth = Ldir['parent'] / 'ptools_output' / 'woac2' / 'sta_df.p'
    if pth.is_file():
        sta_df = pd.read_pickle(pth)

if Ldir['testing']:
    
    # Set the filename (and hence the time) and location
    fn = Ldir['roms'] / 'output' / Ldir['gtagex'] / 'f2019.07.04' / 'ocean_his_0001.nc'
    lon = -123.228000
    lat = 48.240300
    
    # Find indicies nearest to the location
    G, S, T = zrfun.get_basic_info(fn)
    ix = zfun.find_nearest_ind(G['lon_rho'][0,:], lon)
    iy = zfun.find_nearest_ind(G['lat_rho'][:,0], lat)
    
    # Create the output directory and name the output file
    out_dir = Ldir['LOo'] / 'extract' / 'moor' / 'test'
    Lfun.make_dir(out_dir)
    out_fn = out_dir / 'test.nc'
    
    # Run ncks to do the extraction, overwriting any existing file
    cmd_list = ['ncks', '-d', 'xi_rho,'+str(ix), '-d', 'eta_rho,'+str(iy),
        '-v', 'AKs,salt,temp,NO3,phytoplankton,zooplankton,detritus,Ldetritus,oxygen,alkalinity,TIC,h',
        '-O', str(fn), str(out_fn)]
    # Note: 3-D variables will retain singelton dimensions
    # Note: We get AKs so that the s_w dimension is retained
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # and check on the results
    stdout, stderr = proc.communicate()
    if len(stdout) > 0:
        print('\n' + ' sdtout '.center(60,'-'))
        print(stdout.decode())
    if len(stderr) > 0:
        print('\n' + ' stderr '.center(60,'-'))
        print(stderr.decode())
    
    # Add z-coordinates to the file
    foo = nc.Dataset(out_fn, 'a')
    z_rho, z_w = zrfun.get_z(foo['h'][:], np.array([0.]), S)
    vv = foo.createVariable('z_rho', float, ('s_rho',))
    vv.long_name = 'vertical position on s_rho grid, positive up, zero at surface'
    vv.units = 'm'
    vv[:] = z_rho
    vv = foo.createVariable('z_w', float, ('s_w',))
    vv.long_name = 'vertical position on s_w grid, positive up, zero at surface'
    vv.units = 'm'
    vv[:] = z_w
    # and check on the results
    for vn in foo.variables:
        print('%14s: %s' % (vn, str(foo[vn].shape)))
    # Note that z_rho and z_w do not have singleton dimensions
    # Also add units to salt
    foo['salt'].units = 'g kg-1'
    foo.close()
    
    # plot for reality check
    if Ldir['lo_env'] == 'pm_mac':
        a = nc. Dataset(out_fn)
        import matplotlib.pyplot as plt
        plt.close('all')
        fig = plt.figure(figsize=(14,8))
        
        v_list = ['salt', 'temp', 'oxygen']
        NV = len(v_list)
        ii = 1
        for vn in v_list:
            ax = fig.add_subplot(1,NV,ii)
            ax.plot(a[vn][:].squeeze(), a['z_rho'][:], '-o')
            if ii == 1:
                ax.set_ylabel('Z [m]')
            ax.set_title('%s [%s]' % (vn, a[vn].units))
            ii += 1
        a.close()
        plt.show()
    
    
# -------------------------------------------------------

# test for success
if True:
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
