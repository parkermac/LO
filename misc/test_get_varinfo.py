"""
Code to test zrfun.get_varinfo(), which uses the new varinfo.yaml.

RESULT: it appears to work correctly in all cases so far.

"""

from lo_tools import zrfun
from importlib import reload
reload(zrfun)

to_test = 'atm'
vartype='state'

if to_test == 'atm':
    vn_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair','Uwind','Vwind']
    
elif to_test == 'riv':
    vn0_list = ['salt', 'temp', 'NO3', 'NH4', 'LDeN', 'SDeN', 'LDeC', 'SDeC',
            'TIC', 'alkalinity', 'Oxyg']
    vn_list = ['river_'+item for item in vn0_list]
    
elif to_test == 'ocn':
    vartype='climatology'
    avn_list = ['zeta', 'ubar', 'vbar', 'temp', 'salt', 'u', 'v']
    bvn_list = ['NO3', 'NH4', 'chlorophyll', 'phytoplankton', 'zooplankton',
            'LdetritusN', 'SdetritusN', 'LdetritusC', 'SdetritusC',
            'TIC', 'alkalinity', 'oxygen']
    vn_list = avn_list + bvn_list

elif to_test == 'ocn_bry':
    vartype='climatology'
    avn_list = ['zeta', 'ubar', 'vbar', 'temp', 'salt', 'u', 'v']
    bvn_list = ['NO3', 'NH4', 'chlorophyll', 'phytoplankton', 'zooplankton',
            'LdetritusN', 'SdetritusN', 'LdetritusC', 'SdetritusC',
            'TIC', 'alkalinity', 'oxygen']
    vn0_list = avn_list + bvn_list
    vn_list = []
    for D in ['north', 'south', 'east', 'west']:
        for vn in vn0_list:
        # rename variable
        # We have to deal with some inconsistency in the naming of bio boundary variables.
            bvn_dict = {'phytoplankton': 'phyt', 'zooplankton': 'zoop', 'chlorophyll': 'chlo',
                    'LdetritusN': 'LDeN', 'SdetritusN': 'SDeN',
                    'LdetritusC': 'LdeC', 'SdetritusC': 'SdeC'}
            if vn in bvn_dict.keys():
                vn_new = bvn_dict[vn]
            else:
                vn_new = vn
            Vn = vn_new + '_' + D
            vn_list.append(Vn)

vinfo_dict = dict()
for vn in vn_list:
    vinfo_dict[vn] = zrfun.get_varinfo(vn, vartype=vartype)
