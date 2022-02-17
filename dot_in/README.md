# README for dot_in

## The purpose of this code is to create the dot_in files for a ROMS run. These are generally controlled by `LO/driver/driver_roms_mox.py`.

---

### Each of the folder names is the [gtagex] for a run. within a given folder:

`make_dot_in.py` accepts arguments (all handled by `LO/lo_tools/lo_tools/dot_in_argfun.py`) and then fills in fields in the templates BLANK.in and npzd2o_Banas_BLANK.in. It also uses the list of forcings to use from forcing_list.csv.

We keep forcing_list.csv as a separate file because it is also used by `driver_roms_mox.py` to know which forcing to copy from boiler to mox.

---

### shared

This contains the scripts and templates for making the sbatch scripts, e.g. for running mpirun on mox.

In the old LiveOcean framework these were kept in the LiveOcean_roms/makefiles/[exname]. The new organization here makes for much less repetition.

---

### Notes on specific [gtagex]:

- cas6_v0_uu0k is the first attempt to use the updated ROMS on LO_roms_source, with the executable in uu0k in LO_roms_user.  It is designed to run a physics-only version of current forecast (cas6_v0_u0mb).  The "k" is for klone.
- ae0_v0_uu1k is the first attempt to make an analytical run.  It does not have atm forcing, and so I set NFFILES == 0 near the end of BLANK.in.
