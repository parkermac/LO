"""
This is a module of functions used for budget calculations. It is a likely
place where users will need to make their own version, so in code that
uses it we add a "hook" to look for a version in LO_user.

"""
import pandas as pd
import sys

def get_sntup_list(gctag, vol_name):
    """
    The lists returned by this function give the information needed to:
    (1) Specify the open boundary sections, and
    (2) Specify all the segments landward of the open boundary.
    
    The open boundary can consist of one or more sections, since we often
    have cases with this requirement, like needing Admiralty Inlet and
    Deception Pass for a Puget Sound volume.
    
    In "sntup_list" we give the open boundary as a list of tuples, each
    with (section name, sign) where the sign tells you what to multiply
    the transport by in order to get flow INTO the volume. To do this
    right you have to know the sigen convention of each section, which
    you could find by looking at the plot created by plot_collection.py.
    So for the cas6_c0 collection, the ai1 section has positive out to
    sea, so we put -1 as its sign. It is okay for sntup_list to only have
    one tuple, but it must be a list.
    
    In "sect_base_list" we list the first part of each of the sections that
    will be part of the volume. So for example the "base" of ai1 and ai2 is
    "ai". We will use these bases to auto-generate all the possible sections
    that could be in the volume.
    
    In outer_sns_list make a list of the sect + sign of the
    OUTSIDE of the volume, so that we can exclude segments that have
    these as part of their boundaries.
    
    This has lists of sections that are the open boundaries to selected
    volumes.
    
    NOTE on naming conventions:
    "sntup" refers to a tuple like ('ai1',-1)
    "sect_name" or "sn" refers to a section name like 'ai1'
    "sns" refers to a section name + sign like 'ai1_p'
    
    """
    if gctag == 'cas6_c0':
        if vol_name == 'Puget Sound':
            sntup_list = [('ai1',-1),('dp',-1)] # the sign is for inflow direction
            sect_base_list = ['ai','dp','wb','hc','mb','tn','ss']
            outer_sns_list = ['ai1_p','dp_p'] # exclude segments that have these
            
    return sntup_list, sect_base_list, outer_sns_list
