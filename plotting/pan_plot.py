"""
Plot fields in one or more history files.

Examples:

Plot a single figure to the screen with default arguments:
run pan_plot
(it will prompt for the plot type)

Here is an example with explicit flags:
run pan_plot -gtx cas6_v0_live -ro 0 -0 2019.07.04 -lt snapshot -pt P_Chl_DO -avl False

When using the default -avl True for multiple plots (e.g. when making a movie)
the color limits will all be set to match those set by auto_lims() from the first plot.

Use -avl False to have color limits all set to match those set by pinfo.vlims_dict.

Using -test True will create the object "ds", an xarray Dataset of the most recent history file.

The new "mtag" flag allows you to append an optional tag to the end of a movie folder name,
which may be convenient when making the same movie for different time periods. Similarly
the new "mname" tag allows you to have a custom movie name.

"""
import os, sys
import argparse
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi

from lo_tools import Lfun

parser = argparse.ArgumentParser()

# which run to use
parser.add_argument('-gtx', '--gtagex', default='cas7_t0_x4b', type=str)
parser.add_argument('-ro', '--roms_out_num', default=0, type=int)
# 2 = Ldir['roms_out2'], etc.

# select time period and frequency
parser.add_argument('-0', '--ds0', default='2017.07.04', type=str)
parser.add_argument('-1', '--ds1', type=str)
parser.add_argument('-lt', '--list_type', default='snapshot', type=str)
# snapshot, hourly, daily, or lowpass

# arguments that allow you to bypass the interactive choices
parser.add_argument('-hn', '--his_num', default=2, type=int)
parser.add_argument('-pt', '--plot_type', type=str)

# arguments that influence other behavior
#  e.g. make a movie, override auto color limits
parser.add_argument('-mov', '--make_movie', default=False, type=Lfun.boolean_string)
parser.add_argument('-avl', '--auto_vlims', default=True, type=Lfun.boolean_string)
parser.add_argument('-save', '--save_plot', default=False, type=Lfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
parser.add_argument('-mtag', '--movie_tag', default='', type=str)
parser.add_argument('-mname', '--movie_name', default='movie', type=str)

# do things with the arguments
args = parser.parse_args()
argsd = args.__dict__
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set second date string if omitted (need to have Ldir['ds0'])
if Ldir['ds1'] == None:
    Ldir['ds1'] = Ldir['ds0']
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

# choose the type of list to make
if Ldir['list_type'] == None:
    print(' pan_plot '.center(60,'='))
    print('\n%s\n' % '** Choose List type (return for snapshot) **')
    lt_list = ['snapshot', 'daily', 'hourly ', 'allhours']
    Nlt = len(lt_list)
    lt_dict = dict(zip(range(Nlt), lt_list))
    for nlt in range(Nlt):
        print(str(nlt) + ': ' + lt_list[nlt])
    my_nlt = input('-- Input number -- ')
    if len(my_nlt)==0:
        Ldir['list_type'] = 'snapshot'
    else:
        Ldir['list_type'] = lt_dict[int(my_nlt)]
else:
    pass
    
# get the roms_plots module, looking first in LO_user
pth = Ldir['LO'] / 'plotting'
upth = Ldir['LOu'] / 'plotting'
if (upth / 'roms_plots.py').is_file():
    print('Importing roms_plots from LO_user')
    roms_plots = Lfun.module_from_file('roms_plots', upth / 'roms_plots.py')
else:
    print('Importing roms_plots from LO')
    roms_plots = Lfun.module_from_file('roms_plots', pth / 'roms_plots.py')

if Ldir['testing']:
    import roms_plots
    from importlib import reload
    reload(roms_plots)
    

# choose the type of plot to make
if Ldir['plot_type'] == None:
    print('\n%s\n' % '** Choose Plot type (return for P_basic) **')
    pt_list_raw = dir(roms_plots)
    pt_list = []
    for pt in pt_list_raw:
        if pt[:2] == 'P_':
            pt_list.append(pt)
    Npt = len(pt_list)
    pt_dict = dict(zip(range(Npt), pt_list))
    for npt in range(Npt):
        print(str(npt) + ': ' + pt_list[npt])
    my_npt = input('-- Input number -- ')
    if len(my_npt)==0:
        Ldir['plot_type'] = 'P_basic'
    else:
        Ldir['plot_type'] = pt_dict[int(my_npt)]
else:
    pass
    
whichplot = getattr(roms_plots, Ldir['plot_type'])

in_dict = dict()
in_dict['auto_vlims'] = Ldir['auto_vlims']
in_dict['testing'] = Ldir['testing']

# get list of history files to plot
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir,
    Ldir['ds0'], Ldir['ds1'], his_num=Ldir['his_num'])
    
if (Ldir['list_type'] == 'allhours') and Ldir['testing']:
    fn_list = fn_list[:4]

# PLOTTING
outdir0 = Ldir['LOo'] / 'plots'
Lfun.make_dir(outdir0)
if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
plt.close('all')

if len(fn_list) == 1:
    # plot a single image to screen, or to a file
    fn = fn_list[0]
    in_dict['fn'] = fn
    if Ldir['save_plot'] == True:
        in_dict['fn_out'] = outdir0 / (Ldir['list_type'] + '_'
            + Ldir['plot_type'] + '_' + Ldir['gtagex'] + '.png')
    else:
        in_dict['fn_out'] = ''
    whichplot(in_dict)
    
elif len(fn_list) > 1:
    # prepare a directory for results
    if len(Ldir['movie_tag']) > 0:
        outdir = outdir0 / (Ldir['list_type'] + '_' + Ldir['plot_type'] +
            '_' + Ldir['gtagex'] + '_' + Ldir['movie_tag'])
    else:
        outdir = outdir0 / (Ldir['list_type'] + '_' + Ldir['plot_type'] + '_' + Ldir['gtagex'])
    Lfun.make_dir(outdir, clean=True)
    # plot to a folder of files
    jj = 0
    for fn in fn_list:
        nouts = ('0000' + str(jj))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir / outname
        print('Plotting ' + str(fn))
        sys.stdout.flush()
        in_dict['fn'] = fn
        in_dict['fn_out'] = outfile
        whichplot(in_dict)
        # after the first plot we no longer change vlims
        in_dict['auto_vlims'] = False
        jj += 1
    # and make a movie
    if Ldir['make_movie']:
        cmd_list = ['ffmpeg','-r','8','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
            '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
            '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/' + Ldir['movie_name'] + '.mp4']
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print('\n'+stdout.decode())
        if len(stderr) > 0:
            print('\n'+stderr.decode())
            
if Ldir['testing']:
    import xarray as xr
    ds = xr.open_dataset(fn)
