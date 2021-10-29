""""
Dictionaries of defaults to be used for plotting.

"""

from cmocean import cm
cmap_dict =  {'salt': cm.haline,
             'temp': 'rainbow',# cm.thermal,
             'NO3': cm.turbid,
             'phytoplankton': 'Spectral_r',
             'zooplankton': 'Spectral_r',
             'oxygen': cm.oxy,
             'TIC': cm.matter,
             'alkalinity': cm.matter,
             'PH': cm.matter,
             'ARAG': 'coolwarm_r',
             'speed': 'YlOrBr'}
             
vlims_fac_dict =  {'salt': .5,
             'temp': 3,
             'NO3': 1,
             'phytoplankton': 3,
             'zooplankton': 3,
             'oxygen': 1,
             'TIC': 1,
             'alkalinity': 1,
             'PH': 1,
             'ARAG': 1,
             'speed': 1}

# Color limits
# If you use -avl True (the default) then the limits will be set by the first plot
# and then held constant at those levels thereafter.
vlims_dict = {'salt': (14, 35),
        'temp': (7, 18),
        'dye_01': (0,1),
        'NO3': (0, 44),
        'phytoplankton': (0,30),
        'zooplankton': (0, 4),
        'oxygen': (0, 8),
        'TIC': (2000, 2400),
        'alkalinity': (2000,2400),
        'PH': (7, 8.5),
        'ARAG': (0, 3),
        'speed': (0, 2)}

# Units (after multiplying by scaling factor)
units_dict = {'salt': '[ppt]',
             'temp': r' $[^{\circ}C]$',
             'NO3': r'$[\mu mol\ L^{-1}]$',
             'phytoplankton': r'$[mg\ Chl\ m^{-3}]$',
             'zooplankton': r'$[\mu mol\ N\ L^{-1}]$',
             'oxygen': r'$[mg\ L^{-1}]$',
             'TIC': r'$[\mu mol\ L^{-1}]$',
             'alkalinity': r'$[\mu\ equivalents\ L^{-1}]$',
             'PH': '',
             'ARAG': '',
             'speed': r'$[m\ s^{-1}]$'}

# Scaling factors
fac_dict =  {'salt': 1,
             'temp': 1,
             'NO3': 1,
             'phytoplankton': 2.5,
             'zooplankton': 1,
             'oxygen': 32/1000, # convert mmol m-3 to mg L-1
             'TIC': 1,
             'alkalinity': 1,
             'PH': 1,
             'ARAG': 1,
             'speed': 1}
             
# String form to use in titles
tstr_dict = {'salt': 'Salinity',
             'temp': 'Temperature',
             'NO3': 'Nitrate',
             'phytoplankton': 'Phytoplankton',
             'zooplankton': 'Zooplankton',
             'oxygen': 'DO',
             'TIC': 'DIC',
             'alkalinity': 'Alkalinity',
             'PH': 'pH',
             'ARAG': r'$\Omega_{arag}$',
             'speed': 'Speed'}
             
