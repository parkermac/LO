"""
Code to test the module loadDFO_new.py
"""
from datetime import datetime, timedelta
from pathlib import Path

import loadDFO_new as loadDFO
from importlib import reload
reload(loadDFO)

from lo_tools import Lfun
Ldir = Lfun.Lstart()

basedir=str(Ldir['data'] / 'dfo')
dbname='DFO.sqlite'

datelims=(datetime(2016,3,1),datetime(2016,4,1)),
latlims=(49,49.4),
lonlims=(-123.6,-123.2)

# df_bot = loadDFO.loadDFO_bottle(basedir=basedir, dbname='DFO.sqlite',
#     datelims=datelims,latlims=latlims,lonlims=lonlims)

df_ctd = loadDFO.loadDFO_CTD(basedir=basedir, dbname='DFO_CTD.sqlite',
    datelims=datelims)
