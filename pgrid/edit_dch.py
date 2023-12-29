"""
This is a convenience program to allow you to edit the items in the
dch "choices" dictionary, without having to start a new version of the grid.

This should be used sparingly because it circumvents the logic of making
all imortant choices up-front in start_grid.py.

You can have it work on any gridname by specifying it with a command line argument: -g [gridname].

The program works interactively. You select what you want to change and then
type in something that will evaluate to an object of that type. We do this with
a neat function ast.literal_eval(). Note that this can only work with a limited
set of types: strings, bytes, numbers, tuples, lists, dicts, sets, booleans, and None.

To be safe, it saves the changed file to "choices_new.p" and you have to rename this
yourself to "choices.p" if you want pgrid to use it.
"""

import ast
import sys
import pandas as pd
import gfun
import gfun_utility as gfu
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',type=str) # e.g. cas6
args = parser.parse_args()

if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()

Gr =gfun.gstart()

# load the default choices
dch_fn = Gr['gdir'] / 'choices.p'
dch = pd.read_pickle(dch_fn)

# organize for selection
lt_list = list(dch.keys())
Nlt = len(lt_list)
rlt = range(Nlt)
lt_dict = dict(zip(rlt, lt_list))

def print_dch():
    # utility to print current dch
    for ii in rlt:
        k = lt_dict[ii]
        v = dch[k]
        dt = type(v)
        print('%d: %s = %s (%s)' % (ii, k, v, dt))
        
# show original state
print(' Original choices.p for %s'.center(60,'=') % (Gr['gridname']))
print_dch()

# choose which item to edit
print('\n%s' % '** Choose item to change **')
my_ii = input('-- Input number -- ')

# get info for this item
k = lt_dict[int(my_ii)]
v = dch[k]
dt = type(v)

# get new value with the same type using fancy call to the ast module
new_v = ast.literal_eval(input('-- Type new value (must evaluate to correct type!) -- '))
# Note: I don't think we can use this to change Path objects

if type(new_v) != dt:
    # this may catch some errors that ast.literal_eval does not
    print('ERROR: incorrect expression for the expected type.')
    sys.exit()
    
dch[k] = new_v
print(' New choices.p for %s'.center(60,'=') % (Gr['gridname']))
print_dch()

# save the results to choices_new.p
dch_new_fn = Gr['gdir'] / 'choices_new.p'
print('\nSaving results to %s' % (dch_new_fn))
print('\n*** NOTE: to use this in pgrid you must rename it to choices.p by hand ***')
pd.to_pickle(dch, dch_new_fn)



