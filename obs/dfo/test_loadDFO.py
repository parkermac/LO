"""
Code to test the module loadDFO_new.py
"""
from datetime import datetime, timedelta
from pathlib import Path

import loadDFO
from importlib import reload
reload(loadDFO)

from lo_tools import Lfun
Ldir = Lfun.Lstart()

basedir=str(Ldir['data'] / 'obs' / 'dfo')

datelims=[datetime(2016,1,1),datetime(2016,12,31)]
latlims=[48.9,49.5]
lonlims=[-124.2,-123.2]

# df_bot = loadDFO.loadDFO_bottle(basedir=basedir, dbname='DFO.sqlite',
#     datelims=datelims,latlims=latlims,lonlims=lonlims)
#
# df_ctd = loadDFO.loadDFO_CTD(basedir=basedir, dbname='DFO.sqlite',
#     datelims=datelims,latlims=latlims,lonlims=lonlims)
    
df_ctd2 = loadDFO.loadDFOCTD(basedir=basedir, dbname='DFO_CTD.sqlite',
        datelims=datelims)
