"""
Code to make the sta_df.p file for use by extract_casts.py.
"""

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
from datetime import datetime, timedelta

in_dir = Ldir['data'] / 'obs' / 'newport_line'
out_dir = Ldir['LOo'] / 'obs' / 'newport_line'
Lfun.make_dir(out_dir)

fn_list = ['sta_df_NH2017.csv', 'sta_df_NH2018.csv', 'sta_df_NH2019.csv']
for fn in fn_list:
    in_fn = in_dir / fn
    a = pd.read_csv(in_fn)
    a = a.set_index('Station')
    
    dt_list = []
    for t in a['Datetime']:
        dt_list.append(datetime.strptime(t, '%d-%b-%Y'))
    a['Datetime'] = dt_list
    
    if fn == fn_list[0]:
        A = a.copy()
    else:
        A = pd.concat((A,a))
        
A.to_pickle(out_dir / 'sta_df.p')