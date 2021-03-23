"""
This is the main program for making the TIDE forcing file.

Test on mac in ipython:

run make_forcing_main.py -g cas6 -t v3 -f tide0 -r backfill -s continuation -d 2019.07.04 -test True

"""

from pathlib import Path
import sys
from datetime import datetime

pth = Path(__file__).parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import forcing_functions as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import subprocess

cmd = Ldir['which_matlab']

date_string = Ldir['date_string']
out_pth = Path(Ldir['LOo']) / Ldir['gtag'] / ('f' + date_string) / Ldir['frc']
out_dir = str(out_pth)

func = ("make_forcing_worker(\'" +
    str(Ldir['grid']) + "\',\'" +
    str(Ldir['data']) + '/tide' + "\',\'" +
    out_dir + "\',\'" +
    date_string + "\')")
print(func)
    
cmd_list = [Ldir['which_matlab'], "-nodisplay", "-r", func, "&"]

proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = proc.communicate()

if Ldir['testing']:
    print('\n' + ' sdtout '.center(60,'-'))
    print(stdout.decode())
    print('\n' + ' stderr '.center(60,'-'))
    print(stderr.decode())
    
else:
    with open(out_pth / 'Info' / 'screen_output.txt', 'w') as fout:
        fout.write(stdout.decode())
    if len(stderr) > 0:
        with open(out_pth / 'Info' / 'subprocess_error.txt', 'w') as ffout:
            ffout.write(stderr.decode())

#make_tides(gridname, tag, date_string, run_type, outdir)


result_dict['result'] = 'success' # success or fail

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
