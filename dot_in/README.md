# README for dot_in

### The purpose of this code is to create the dot_in files for a ROMS run. These are generally controlled by `LO/driver/driver_roms4.py` or similar.

---

### Each of the folder names is the [gtagex] for a run. Within a given folder:

`make_dot_in.py` accepts arguments (all handled by `LO/lo_tools/lo_tools/dot_in_argfun.py`) and then fills in fields in the templates `BLANK.in` and `bio_Fennel_BLANK.in`. It also uses the list of forcings and open boundaries to use from `forcing_list.csv`.

We keep forcing_list.csv as a separate file because it is also used by `driver_roms4.py` and similar to know which forcing to copy from our server (e.g. apogee or perigee) to klone.

---

## Notes on specific [gtagex] cases. Not all instances are in the folder, so these notes reflect some development steps that may have been abandoned.

### CURRENT

#### Part of the daily forecast
- `cas7_t0_x4b` Primary forecast and long hindcast (late 2021 to present). Uses traps rivers and increased light ettenuation in the Salish Sea
- `oly1_t0_xn4b`100 m South Sound nest, with WET_DRY
- `wgh2_t0_xn0b` 200 m Willapa Bay and Grays Harbor nest, with WET_DRY

#### Current development runs
- `cas7_t1_x11ab` **Leading edge of our development**, with bio code from x10ab. Has changed TIDE_START. Uses new ROMS (4.3) presumably with Harcourt turbulence/advection improvements. Uses OMEGA_IMPLICIT. These allow DT = 60 sec instead of 40. 
- `cas7_t1_x11b` **Candidate for new primary forecast**, like cas7_t1_x11ab but with hourly saves and no averages.
- `cas7_t1_x10ab` experiment along the way to x11ab. Has 50% burial of particles (including Carbon) in the Salish Sea, 2017 run used for validation. The validation results look great, even for pCO2. I saved a version of 2017 output from  and earlier version of this that used HSIMT for bio tracers and renamed it cas7_hsimt_x10ab. (on apogee). The "a" in the executable name implies that we use daily saves for speed and low storage and create a daily **average** for validation.
- `wgh2_t0_xn4b` like wgh2_t0_xn0b but using updated bio code based on x4b. Has 2017 run done for mooring validation.

#### Runs to experiment with fixing tide phase
- `cas7_t0_x4tf` tide experiment using the tidal tractive force
- `cas7_t1_x4` tide experiment with changed TIDE_START
- `cas7_t1_x4tf` tide experiment using the tractive force with changed TIDE_START

#### Miscellany
- `cas7_t1_x4a` experiment using averages
- `cas2k_v0_x2b` 2 km grid developed for Samantha Siedlecki's use
- `cas7_t0_x4` physics-only test of the new cas7 grid, which is a slightly modified version of cas6. Includes Agate Pass and Swinomish Channel.
- `ae0_t0_xa0` analytical run
- `ae0_v0_uu1k` analytical run

---

### OLD: RECENT BUT STILL OBSOLETE

**These cases all use (at the time) up-to-date ROMS code. The associated ROMS code is in the repos LO_roms_source, LO_roms_source_alt, LO_roms_users. The README in LO_roms_user has more description of how the updated ROMS code is organized.**

- `cas6_v0_uu0k` is a physics-only version of current forecast. It is the first (SUCCESSFUL!) attempt to use the up-to-date ROMS source code, instead of the 2016 version in LiveOcean.    The double letter **uu** in [ex_name] is my shorthand that this is using the new ROMS, but that is not likely to be a strict rule going forward.
- `cas6_v0_uu0kb` is like cas6_v0_uu0k but also uses the bio code (the Fennel code modified to closely resemble the Banas-Siedlecki-Davis equations).
- `ae0_v0_uu1k` is an analytical run (idealized estary-shelf). This has no atm forcing.
- `cas6_v00_uu0m` is a physics-only version of the updated model (new ROMS, perfect restart, EMINUSP, etc.) It was created as a one-off for Ramsey Harcourt
- `cas6_v01_uu0m` is just like cas6_v00_uu0m except one parameter in BLANK.in relating to the turbulence closure has been changed (GLS_Kmin set to 1e-10). RESULT: It blew up immediately and would not run.
- `cas6_v00_x0mb` was the current base case for the ROMS update project. It uses the Fennel bio code (with NH4), modified to be as consistent as possible with "BSD as written" (i.e. the equations and parameters from the Banas, Siedlecki, and Davis papers of the RISE/PNWTOX era, but accounting for Ammonium cycling). The bio code is contained in LO_roms_source_alt/npzd_banas.
- `cas6_traps2_x0mb` is exactly like `cas6_v00_x0mb` except that the river forcing comes from Aurora Leeson's traps2 code. This includes all the tiny rivers and point sources as horizontal sources, bypassing the current problem with vertical sources.
- `cas6_traps2_x1b` is like `cas6_traps2_x0mb` except I modified the fennel.h code to use "optimum uptake" in the nutrient limitation for NH4 (it was already in NO3). I expect this will significantly decrease the productivity, and hopefully solve our issue that DO is too high and DIN is too low. Also, this is the first attempt to isolate the edited bio code in the [ex_name] folder. I also dropped the "m" for mox. It should run on both klone and mox interchangeably.

---

### VERY_OLD

**These all use the older (2016) LiveOcean ROMS source code:**

- `cas6_v0_u0kb` was the current daily forecast, run on klone.  In the [ex_name] **u0** is just a version number, **k** means it is to run on klone, and **b** means it has biology.
- `cas6_v0_u0mb` is the current daily forecast backup, run on mox.
- `ai0_v0_n0k` is a nested run of a high-resolution grid covering just Admiralty Inlet.
- `hc0_v0_n0k` is a nested run of a high-resolution grid covering just Hood Canal.
- `so0_v0_n0k` is a nested run of a high-resolution grid covering just South Sound.
- `cas6_v3t075_lo8` is a run like the current forecast but with 75% strength tides
- `cas6_v3t110_lo8` is a run like the current forecast but with 110% strength tides
- `shared` is a folder of batch script stuff that was replaced by LO/driver/batch.
