"""
Plot fields in one or more history files.

Testing calls from ipython on mac:

run dm_plot

python dm_plot.py -vn speed -dom PS -mov True -lt hourly < /dev/null > test.log &

run dm_plot -tracks True -mov True -lt hourly

This call tests a new run with perfect restart file numbering. Using "temp" as the
variable name means it will not overwrite a forecast output file.
run dm_plot -tracks True -mov True -lt hourly -ds0 2021.07.04 -gtx cas6_v00_uu0mb -vn temp

"""

import argparse
import numpy as np
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi

from lo_tools import Lfun

from importlib import reload
import plots; reload(plots)
import dm_pfun; reload(dm_pfun)
import pinfo; reload(pinfo)

import matplotlib.pyplot as plt
plt.close('all')

parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-gtx', '--gtagex', type=str, default='cas6_v0_live')
parser.add_argument('-ro', '--roms_out_num', type=int, default=0)
parser.add_argument('-ds0', type=str, default='2019.07.04')
parser.add_argument('-ds1', type=str, default='')
# arguments that allow you to bypass the interactive choices
parser.add_argument('-hn', type=int, default=2) # history file number
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

gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in Q.keys():
    if a not in Ldir.keys():
        Ldir[a] = Q[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
    
whichplot = getattr(plots, Q['pt'])

# get list of history files to plot
ds0 = Q['ds0']
ds1 = Q['ds1']
fn_list = Lfun.get_fn_list(Q['lt'], Ldir, ds0, ds1, his_num=Q['hn'])

dm_pfun.get_ax_limits(Q)

if Q['avl'] == False:
    Q['vmin'] = pinfo.vlims_dict[Q['vn']][0]
    Q['vmax'] = pinfo.vlims_dict[Q['vn']][1]

M = dm_pfun.get_moor_info(Q)
dm_pfun.get_moor(ds0, ds1, Ldir, Q, M)

if Q['tracks']:
    dm_pfun.get_tracks(Q, Ldir)

# PLOTTING

if Q['bot'] == True:
    bot_tag = 'bot'
else:
    bot_tag = 'top'

if len(fn_list) == 1:
    # plot a single image to a file
    fn = fn_list[0]
    Q['fn'] = fn
    plotname = ('daymovie_' + Q['gtagex'] + '_' + Q['pt'] + '_'
        + Q['dom'] + '_' + Q['vn'] + '_' + bot_tag + '.png')
    Q['fn_out'] = Ldir['LOo'] / 'plots' / plotname
    whichplot(Q, M)
    
elif len(fn_list) > 1:
    # prepare a directory for results
    outdir0 = Ldir['LOo'] / 'daymovie' / Ldir['gtagex']
    Lfun.make_dir(outdir0)
    moviename = Q['pt'] + '_' + Q['dom'] + '_' + Q['vn'] + '_' + bot_tag
    outdir = outdir0 / moviename
    Lfun.make_dir(outdir, clean=True)
    # plot to a folder of files
    jj = 0
    for fn in fn_list:
        nouts = ('0000' + str(jj))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir / outname
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
        cmd_list = ['ffmpeg','-r','8','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
            '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
            '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/movie.mp4']
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if True:
            if len(stdout) > 0:
                print('\n'+stdout.decode())
            if len(stderr) > 0:
                print('\n'+stderr.decode())
        