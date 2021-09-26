""""
Dictionaries of defaults to be used for plotting.

"""

from cmocean import cm

# default figure size
figsize = (13,8) # laptop

# Color limits
# If you use () then the limits will be set by the first plot
# and then held constant at those levels thereafter.    
vlims_dict = {'salt': (14, 35),
        'temp': (7, 18),
        'dye_01': (0,1),
        'NO3': (0, 44),
        'phytoplankton': (0,30),
        'zooplankton': (0, 4),
        'oxygen': (0, 10),
        'TIC': (2000, 2400),
        'alkalinity': (2000,2400),
        'PH': (7, 8.5),
        'ARAG': (.2, 2.2),
        'detritus': (0,1),
        'Ldetritus': (0,1)}

# Colormaps (use _r for reverse)
cmap_dict = {'salt': cm.haline,#'Spectral_r',
             'temp': 'RdYlBu_r',#'bwr',
             'NO3': 'jet',
             'phytoplankton': 'ocean_r',
             'zooplankton': 'jet',
             'oxygen': cm.oxy,
             'TIC': 'rainbow',
             'alkalinity': 'rainbow',
             'PH': 'jet',
             'ARAG': 'rainbow',
             'detritus': 'rainbow',
             'Ldetritus': 'rainbow',
             'u': 'rainbow',
             'v': 'rainbow',
             'w': 'rainbow',
             'zeta': 'rainbow'}

# Units (after multiplying by scaling factor)
units_dict = {'salt': '$(g\ kg^{-1})$',
             'temp': ' $(^{\circ}C)$',
             'NO3': ' $(\mu mol\ L^{-1})$',
             'phytoplankton': ' $(mg\ Chl\ m^{-3})$',
             'zooplankton': ' $(\mu mol\ N\ L^{-1})$',
             'oxygen': ' $(mg\ L^{-1})$',
             'TIC': ' $(\mu mol\ L^{-1})$',
             'alkalinity': ' $(\mu\ equivalents\ L^{-1})$',
             'PH': '',
             'ARAG': '',
             'detritus': ' $(\mu mol\ L^{-1})$',
             'Ldetritus': ' $(\mu mol\ L^{-1})$',
             'w': ' $(m\ s^{-1})$',
             'u': ' $(m\ s^{-1})$',
             'v': ' $(m\ s^{-1})$',
             'ubar': ' $(m\ s^{-1})$',
             'vbar': ' $(m\ s^{-1})$',
             'zeta': ' (m)'}

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
             'detritus': 1,
             'Ldetritus': 1,
             'w': 1,
             'u': 1,
             'v': 1,
             'ubar': 1,
             'vbar': 1,
             'zeta': 1}
             
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
             'ARAG': '$\Omega_{arag}$',
             'detritus': 'Detritus',
             'Ldetritus': 'Ldetritus',
             'w': 'W',
             'u': 'U',
             'v': 'V',
             'ubar': 'Ubar',
             'vbar': 'Vbar',
             'zeta': 'Zeta'}
             
# this is used by plotting_functions.auto_lims to decide how many
# standard deviations +/- to set the color limits relative to the mean
range_dict =  {'salt': 1,
             'temp': 2.5,
             'NO3': 3,
             'phytoplankton': 3,
             'zooplankton': 3,
             'oxygen': 3,
             'TIC': 3,
             'alkalinity': 3,
             'PH': 3,
             'ARAG': 3,
             'detritus': 3,
             'Ldetritus': 3,
             'w': 3,
             'u': 3,
             'v': 3,
             'ubar': 3,
             'vbar': 3,
             'zeta': 3}
