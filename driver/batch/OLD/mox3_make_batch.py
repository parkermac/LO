"""
This creates a batch script for LiveOcean -- using a template version of 
the batch file replacing things like $whatever$ with actual values and paths.
- Initial version: dsd, November 2015
- Modified and reorganized for LO: PM, 2021_05

"""

from pathlib import Path
# get the path to this directory
pth = Path(__file__).absolute().parent

# get command line arguments
import argparse

parser = argparse.ArgumentParser()

# positional arguments
parser.add_argument('-xd', '--roms_ex_dir', type=str, help="path to ROMS executable")
parser.add_argument('-rxn', '--roms_ex_name', type=str, help="name of the ROMS executable")
parser.add_argument('-rod', '--roms_out_dir', type=str, help="path to ROMS output")
parser.add_argument('-np', '--np_num', type=int, help="number of cores to use")
parser.add_argument('-N', '--cores_per_node', type=int, help="cores per node")
parser.add_argument('-x', '--ex_name', type=str, help="executable name")
parser.add_argument('-j', '--jobname', type=str, help="jobname")

args = parser.parse_args()
in_dict = args.__dict__
in_dict['node_num'] = max((1, int(in_dict['np_num'] / in_dict['cores_per_node'])))

## create lo_batch.sh - batch job script  ##########################

f  = open(pth / 'mox3_batch_BLANK.sh','r')
f2 = open(Path(in_dict['roms_out_dir']) / 'mox_batch.sh','w')

for line in f:
    for var in in_dict.keys():
        if '$'+var+'$' in line: 
            line2 = line.replace('$'+var+'$', str(in_dict[var]))
            line = line2
        else:
            line2 = line
    f2.write(line2)

f.close()
f2.close()
