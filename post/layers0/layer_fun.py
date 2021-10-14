"""
Functions for the layers.
"""
import os

def get_fn_list(in_dir, testing):
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item and '.nc' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    # shorten the list to be every 4 hours
    if testing:
        fn_list = fn_list[::8]
    else:
        fn_list = fn_list[::4]
    return fn_list