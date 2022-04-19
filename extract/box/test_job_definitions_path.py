"""
Code to test import path specification.
"""
import sys
from lo_tools import Lfun
Ldir = Lfun.Lstart()

import importlib.util

def module_from_file(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

pth = Ldir['LO'] / 'extract' / 'box'
upth = Ldir['LOu'] / 'extract' / 'box'

if (upth / 'job_definitions.py').is_file():
    print('Importing from LO_user')
    job_definitions = module_from_file('job_definitions', upth / 'job_definitions.py')
    aa, vn_list = job_definitions.get_box('LOu_test', 0, 0)
else:
    print('Importing from LO')
    job_definitions = module_from_file('job_definitions', pth / 'job_definitions.py')
    aa, vn_list = job_definitions.get_box('sequim0', 0, 0)
    

