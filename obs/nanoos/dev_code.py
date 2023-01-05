"""
Code to help with the development of the nanoos processing.

This gets all the column headers from every cruise I have data for, and then
gloms them together into a sorted list without repeats. Then it prints them
to the screen in a format that makes it easy to edit as a dictionary in
bot_fun.py.

"""

import pandas as pd

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# BOTTLE
source = 'nanoos'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source
year_list = [2017, 2021]
a = []
for year in year_list:
    ys = str(year)
    in_dirs = list(in_dir0.glob('*'+ys+'*'))
    fn_list = []
    for in_dir in in_dirs:
        fn_list.append(list(in_dir.glob('*labupcast*'))[0])
    for fn in fn_list:
        print(fn.name)
        df = pd.read_excel(fn)
        a += list(df.columns)
b = list(set(a))
b.sort()
print('')
for v in b:
    print("\'%s\':\'\'," % (v))