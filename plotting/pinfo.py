""""
Dictionaries of defaults to be used for plotting.

"""

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
        'oxygen': (0, 8),
        'TIC': (2000, 2400),
        'alkalinity': (2000,2400),
        'PH': (7, 8.5),
        'ARAG': (.2, 2.2),
        'Ldetritus': ()}

# Colormaps (use _r for reverse)
cmap_dict = {'salt': 'Spectral_r',#'rainbow',
             'temp': 'bwr',#run 'jet',
             'NO3': 'jet',
             'phytoplankton': 'ocean_r',
             'zooplankton': 'jet',
             'oxygen': 'rainbow_r',#'jet_r',
             'TIC': 'rainbow',
             'alkalinity': 'rainbow',
             'PH': 'jet',
             'ARAG': 'rainbow',
             'Ldetritus': 'rainbow',
             'w': 'rainbow'}

# Units (after multiplying by scaling factor)
units_dict = {'salt': '',
             'temp': 'degC',#' $(^{\circ}C)$',
             'NO3': ' $(\mu mol\ L^{-1})$',
             'phytoplankton': ' $(mg\ Chl\ m^{-3})$',
             'zooplankton': ' $(\mu mol\ N\ L^{-1})$',
             'oxygen': ' $(ml\ L^{-1})$',
             'TIC': ' $(\mu mol\ L^{-1})$',
             'alkalinity': ' $(\mu\ equivalents\ L^{-1})$',
             'PH': '',
             'ARAG': '',
             'Ldetritus': '',
             'w': ' m/s',#' $(m s^{-1})$',
             'u': ' m/s',#' $(m s^{-1})$',
             'v': ' m/s',#' $(m s^{-1})$',
             'ubar': ' m/s',#' $(m s^{-1})$',
             'vbar': ' m/s',#' $(m s^{-1})$',
             'zeta': ' (m)'}

# Scaling factors
fac_dict =  {'salt': 1,
             'temp': 1,
             'NO3': 1,
             'phytoplankton': 2.5,
             'zooplankton': 1,
             'oxygen': 0.032/1.42903, # convert mmol m-3 to ml L-1
             'TIC': 1,
             'alkalinity': 1,
             'PH': 1,
             'ARAG': 1,
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
             'Ldetritus': 'Ldetritus',
             'w': 'W',
             'u': 'U',
             'v': 'V',
             'ubar': 'Ubar',
             'vbar': 'Vbar',
             'zeta': 'Zeta'}
             
