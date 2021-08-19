"""
Code to simplify pushing a box extraction to azure.
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# choose the file
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'box'
box_name = Lfun.choose_item(in_dir, tag='.nc', exclude_tag='', itext='** Choose box extraction from list **')
box_fn = in_dir / box_name

# this will give screen output of the URL for sharing
az_dict = Lfun.copy_to_azure(box_fn, box_name, 'pm-share', Ldir)
if az_dict['result'] =='success':
    print('USE THIS URL TO ACCESS THE FILE')
    print(az_dict['az_url'])
elif az_dict['result'] =='fail':
    print('EXCEPTION')
    print(az_dict['exception'])