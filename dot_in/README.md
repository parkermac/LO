# README for dot_in

### The purpose of this code is to create the dot_in files for a ROMS run. These are generally controlled by `LO/driver/driver_roms2.py` or similar.

---

### Each of the folder names is the [gtagex] for a run. Within a given folder:

`make_dot_in.py` accepts arguments (all handled by `LO/lo_tools/lo_tools/dot_in_argfun.py`) and then fills in fields in the templates `BLANK.in` and `npzd2o_Banas_BLANK.in` or the newer `bio_Fennel_BLANK.in`. It also uses the list of forcings and open boundaries to use from `forcing_list.csv`.

We keep forcing_list.csv as a separate file because it is also used by `driver_roms2.py` and similar to know which forcing to copy from our server (e.g. apogee or perigee) to hyak (klone or mox).

---

### shared

This contains the scripts and templates for making the sbatch scripts, e.g. for running mpirun on hyak.

**OBSOLETE: replaced by `LO/driver/batch`.**

---

### Notes on specific [gtagex] cases:

**These all use the older (2016) LiveOcean ROMS source code:**

- `cas6_v0_u0kb` is the current daily forecast, run on klone.  In the [ex_name] **u0** is just a version number, **k** means it is to run on klone, and **b** means it has biology.
- `cas6_v0_u0mb` is the current daily forecast backup, run on mox.
- `ai0_v0_n0k` is a nested run of a high-resolution grid covering just Admiralty Inlet.
- `hc0_v0_n0k` is a nested run of a high-resolution grid covering just Hood Canal.
- `so0_v0_n0k` is a nested run of a high-resolution grid covering just South Sound.
- `cas6_v3t075_lo8` is a run like the current forecast but with 75% strength tides
- `cas6_v3t110_lo8` is a run like the current forecast but with 110% strength tides

**These cases all use the up-to-date ROMS code. The associated ROMS code is in the repos LO_roms_source, LO_roms_source_alt, LO_roms_users. The README in LO_roms_user has more description of how the updated ROMS code is organized.**

- `cas6_v0_uu0k` is a physics-only version of current forecast. It is the first (SUCCESSFUL!) attempt to use the up-to-date ROMS source code, instead of the 2016 version in LiveOcean.    The double letter **uu** in [ex_name] is my shorthand that this is using the new ROMS, but that is not likely to be a strict rule going forward.
- `cas6_v0_uu0kb` is like cas6_v0_uu0k but also uses the bio code (the Fennel code modified to closely resemble the Banas-Siedlecki-Davis equations).
- `ae0_v0_uu1k` is an analytical run (idealized estary-shelf). This has no atm forcing.
- `cas6_v00_uu0m` is a physics-only version of the updated model (new ROMS, perfect restart, EMINUSP, etc.) It was created as a one-off for Ramsey Harcourt
- `cas6_v01_uu0m` is just like cas6_v00_uu0m except one parameter in BLANK.in relating to the turbulence closure has been changed (GLS_Kmin set to 1e-10). RESULT: It blew up immediately and would not run.
- `cas6_v00_x0mb` is the current base case for the ROMS update project. It uses the Fennel bio code (with NH4), modified to be as consistent as possible with "BSD as written" (i.e. the equations and parameters from the Banas, Siedlecki, and Davis papers of the RISE/PNWTOX era, but accounting for Ammonium cycling). The bio code is contained in LO_roms_source_alt/npzd_banas.
- `cas6_traps2_x0mb` is exactly like `cas6_v00_x0mb` except that the river forcing comes from Aurora Leeson's traps2 code. This includes all the tiny rivers and point sources as horizontal sources, bypassing the current problem with vertical sources.
- `cas6_traps2_x1b` is like `cas6_traps2_x0mb` except I modified the fennel.h code to use "optimum uptake" in the nutrient limitation for NH4 (it was already in NO3). I expect this will significantly decrease the productivity, and hopefully solve our issue that DO is too high and DIN is too low. Also, this is the first attempt to isolate the edited bio code in the [ex_name] folder. I also dropped the "m" for mox. It should run on both klone and mox interchangeably.
