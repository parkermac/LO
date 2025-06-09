"""
This convenience program updates the pre-processed varinfo file
LO_data/varinfo/varinfo_list.p.

This file is created from ROMS/External/varinfo.yaml by zrfun.make_varinfo_list()
and the resulting varinfo_list.p is used when creating forcing files. It ensures
that the naming and metadata are consistent with ROMS.

We do this just pre-processing step just to speed things up.

NOTE: you should run this program every time you do "git pull" in
LO_roms_source_git.
"""

from lo_tools import zrfun

zrfun.make_varinfo_list()