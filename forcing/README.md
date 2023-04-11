# README for LO/forcing

## These folders contain the code to generate NetCDF forcing files for ROMS runs.

#### A run with complete forcing [frc] will need:
- Atmospheric forcing (atm)
- Ocean initial and boundary conditions (ocn)
- Tidal forcing (tide)
- River transport (riv)

#### The naming of versions is not completely rigorous, but in general:
- [frc]0 is the first working version of LO forecast forcing
- [frc]1 is me testing some improvement to frc[0] such as using more up-to-date time coordinate names
- [frc]N# is for nested runs, meaning it gets its values by interpolating from ROMS history files
- [frc]A# is for analytical runs
- [frc]## (repeated number, like ocn00) is for use with the updated ROMS, including automated coordinate naming using varinfo.yaml, and the banas-fennel biogeochemical model.

---

#### Specific Notes

**riv01** is like riv00 but it uses the new pre/river1 system, with more robust naming, for historical and climatological data.
