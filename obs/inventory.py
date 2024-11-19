"""
Code to go through all the observational bottle (and maybe ctd) data and
inventory what we have by year, source, and month.

I am writing this to help with the obsmod.html project, where I am worried
about some data gaps that may indicate errors in processing.
"""

import pandas as pd
from lo_tools import Lfun

Ldir = Lfun.Lstart()

otype = 'bottle'
year_list = list(range(2013,2024))
month_list = list(range(1,13))
source_list = ['dfo1', 'ecology_nc', 'nceiSalish', 'nceiCoastal',
    'LineP', 'nceiPNW', 'WOD']

for source in source_list:
    dfi = pd.DataFrame(index=year_list,columns=month_list)
    for year in year_list:
        try:
            df = pd.read_pickle(Ldir['LOo'] / 'obs' / source / otype / (str(year) + '.p'))
            for month in month_list:
                dfm = df[df['time'].dt.month==month]
                dfi.loc[year,month] = len(dfm)
        except:
            dfi.loc[year,:] = 0
    print(68*'=')
    print(source)
    print(dfi)


