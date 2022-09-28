#!/bin/bash

ymd=$(date "+%Y.%m.%d")
echo "Copying movies and tracks to homer for "$ymd

scp /dat1/parker/LO_output/post/cas6_v0_u0kb/f$ymd/daymovie0/*.mp4 pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/Figs_active_forecast
scp /dat1/parker/LO_output/tracks/full_surf_forWeb/tracks_full.json pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks/
scp /dat1/parker/LO_output/tracks/PS_surf_forWeb/tracks_PS.json pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks/