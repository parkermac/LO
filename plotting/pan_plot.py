"""
Plot fields in one or more history files.

Examples:

Plot a single figure to the screen with default arguments
run pan_plot -gtx cas6_v3_lo8b -test True

Here is an example with explicit flags
run pan_plot -gtx cas6_v3_lo8b -ro 2 -0 2019.07.04 -lt snapshot -pt P_Chl_DO -avl False

Save multiple plots with color limits all set to match those set by
auto_lims() from the first plot

Use -avl False to have color limits all set to match those set by
pinfo.vlims_dict.

"""
import os, sys
import argparse
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from lo_tools import Lfun
import roms_plots
from importlib import reload
reload(roms_plots)

parser = argparse.ArgumentParser()

# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.

# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str) # list type: snapshot, hourly, or daily

# arguments that allow you to bypass the interactive choices
parser.add_argument('-hn', '--his_num', type=int, default=1)
parser.add_argument('-pt', '--plot_type', type=str)

# arguments that influence other behavior
#  e.g. make a movie, override auto color limits
parser.add_argument('-mov', '--make_movie', default=False, type=Lfun.boolean_string)
parser.add_argument('-avl', '--auto_vlims', default=True, type=Lfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)

# do things with the arguments
args = parser.parse_args()
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# testing
if Ldir['testing']:
    Ldir['roms_out_num'] = 2
    Ldir['ds0'] = '2019.07.04'
    Ldir['ds1'] = '2019.07.04'
    Ldir['list_type'] = 'snapshot'
    Ldir['plot_type'] = 'P_basic'
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
plt.close('all')

if len(fn_list) == 1:
    # plot a single image to screen
    fn = fn_list[0]
    in_dict['fn'] = fn
    in_dict['fn_out'] = ''
    whichplot(in_dict)
    
elif len(fn_list) > 1:
    # prepare a directory for results
    outdir0 = Ldir['LOo'] / 'plots'
    Lfun.make_dir(outdir0, clean=False)
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
        ff_str = ("ffmpeg -r 8 -i " + 
        str(outdir)+"/plot_%04d.png -vcodec libx264 -pix_fmt yuv420p -crf 25 "
        +str(outdir)+"/movie.mp4")
        os.system(ff_str)
