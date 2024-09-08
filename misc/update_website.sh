#!/bin/bash

ymd=$(date "+%Y.%m.%d")
echo "Copying movies and tracks to homer for "$ymd

scp /dat1/parker/LO_output/post/cas7_t0_x4b/f$ymd/daymovie0/*.mp4 pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/Figs_active_forecast
scp /dat1/parker/LO_output/tracks2/cas7_t0_x4b/json_files/*.json pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks2/
scp /dat1/parker/LO_output/tracks2/wgh2_t0_xn0b/json_files/*.json pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks2/