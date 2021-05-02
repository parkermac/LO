"""
A convenience program for loading the flow DataFrame.
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
import matplotlib.pyplot as plt
import numpy as np

gtag = 'cas6_v3'
year0 = 1980
year1 = 2020

# location of historical data to plot
riv_dir = Ldir['LOo'] / 'pre' / 'river' / gtag / 'Data_historical'
all_df = pd.read_pickle(riv_dir / ('ALL_flow_' + str(year0) + '_' + str(year1) + '.p'))