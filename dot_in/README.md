# README for dot_in

## The purpose of this code is to create the dot_in files for a ROMS run. These are generally controlled by `LO/driver/driver_roms_mox.py`.

---

### Each of the folder names is the [gtagex] for a run. within a given folder:

`make_dot_in.py` accepts arguments (all handled by `LO/alpha/dot_in_argfun.py`) and then fills in fields in the templates BLANK.in and npzd2o_Banas_BLANK.in. It also uses the list of forcings to use from forcing_list.csv.

We keep forcing_list.csv as a separate file because it is also used by `driver_roms_mox.py` to know which forcing to copy from boiler to mox.

---

### shared

This contains the scripts and templates for making the sbatch scripts, e.g. for running mpirun on mox.

In the old LiveOcean framework these were kept in the LiveOcean_roms/makefiles/[exname]. The new organization here makes for much less repetition.
