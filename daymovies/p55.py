"""
Plot fields in one or more history files.

Testing calls from ipython on mac:

run p5
run p5 -tracks True -mov True -lt hourly

"""

import os, sys
import argparse
import numpy as np
from datetime import datetime, timedelta

from lo_tools import Lfun

from importlib import reload
import plots; reload(plots)
import pfun; reload(pfun)
import pinfo; reload(pinfo)

import matplotlib.pyplot as plt
plt.close('all')

parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-gridname', type=str, default='cas6')
parser.add_argument('-tag', type=str, default='v3')
parser.add_argument('-ex_name', type=str, default='lo8b')
parser.add_argument('-ds0', type=str, default='2019.07.04')
parser.add_argument('-ds1', type=str, default='')
# arguments that allow you to bypass the interactive choices
parser.add_argument('-hn', type=int, default=1) # history file number
parser.add_argument('-lt', type=str, default='snapshot')
parser.add_argument('-pt', type=str, default='P1')
parser.add_argument('-vn', type=str, default='salt')
parser.add_argument('-dom', type=str, default='full')
parser.add_argument('-bot', default=False, type=Lfun.boolean_string)
parser.add_argument('-mov', default=False, type=Lfun.boolean_string)
parser.add_argument('-avl', default=True, type=Lfun.boolean_string)
parser.add_argument('-emask', default=False, type=Lfun.boolean_string)
parser.add_argument('-tracks', default=False, type=Lfun.boolean_string)
parser.add_argument('-ttag', default='base', type=str) # tag for tracking
parser.add_argument('-test', default=False, type=Lfun.boolean_string)

args = parser.parse_args()
Q = args.__dict__

if len(Q['ds1']) == 0:
    Q['ds1'] = Q['ds0']
    
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + Q['ex_name']

whichplot = getattr(plots, Q['pt'])

# get list of history files to plot
ds0 = Q['ds0']
ds1 = Q['ds1']
fn_list = Lfun.get_fn_list(Q['lt'], Ldir, ds0, ds1, his_num=Q['hn'])

pfun.get_ax_limits(Q)

if Q['avl'] == False:
    Q['vmin'] = pinfo.vlims_dict[Q['vn']][0]
    Q['vmax'] = pinfo.vlims_dict[Q['vn']][1]

M = pfun.get_moor_info(Q)
pfun.get_moor(ds0, ds1, Ldir, Q, M)

if Q['tracks']:
    pfun.get_tracks(Q, Ldir)

# PLOTTING

if len(fn_list) == 1:
    # plot a single image to screen
    fn = fn_list[0]
    Q['fn'] = fn
    Q['fn_out'] = ''
    whichplot(Q, M)
    
elif len(fn_list) > 1:
    # prepare a directory for results
    outdir0 = Ldir['LOo'] / 'p5' / Ldir['gtagex']
    Lfun.make_dir(outdir0, clean=False)
    if Q['bot'] == True:
        bot_tag = 'bot'
    else:
        bot_tag = 'top'
        
    moviename = Q['pt'] + '_' + Q['dom'] + '_' + Q['vn'] + '_' + bot_tag
    outdir = outdir0 + moviename + '/'
    Lfun.make_dir(outdir, clean=True)
    # plot to a folder of files
    jj = 0
    for fn in fn_list:
        nouts = ('0000' + str(jj))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
        if np.mod(jj,10) == 0:
            print('Plot %d out %d' % (jj, len(fn_list)))
        Q['fn'] = fn
        Q['fn_out'] = outfile
        whichplot(Q, M)
        # after the first plot we no longer change vlims
        Q['avl'] = False
        jj += 1
    # and make a movie
    if Q['mov']:
        ff_str = ("ffmpeg -r 8 -i " + 
            outdir + "plot_%04d.png -vcodec libx264 -pix_fmt yuv420p -crf 25 "
            + outdir + moviename + ".mp4")
        os.system(ff_str)
