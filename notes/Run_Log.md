# Run Log

### Info about where the output files for some important runs are archived on our machines.

---

## Salish Runs [total 5 TB]

**perigee: /data2/parker/archive/salish**

- salish_2006_3 This is the run that Dave Sutherland used for his 2011 JPO "MoSSea rollout" paper. Uses Global-NCOM for its ocean bc's. History files deleted, but analysis files saved. [90 GB]
- salish_2006_4 This is just like _3 except that the open boundary conditions come from the higher resolution NCOM CCS. Also this one has (I think) eddy viscosity and diffusivity saved in all files. This was used for Jamie's work. NOTE: It appears to have zero values for sustr after the first month or so. [2.4 TB]
also _lp (low-passed) versions of 2005_1 and 2006_4
- salish_2005_1 This is there, all 8715 files. Based on their size they do not have eddy viscosity or diffusivity. [2 TB]
- ainlet_2006_3 Dave's 100 m resolution nested submodel of Admiralty Inlet. Has saves 2 to 771 (about 16 days of half-hourly saves). Each save file is 55 MB. [62 GB]

---

## PNWTOX Runs [total 9 TB]

**perigee: /data2/parker/archive/pnwtox**

- T40final_base_2004
- T40final_base_2005_noramp
- T40final_base_2006_noramp_v2
- T40final_base_2007_noramp
- T40final_base_2007_noramp_lp
- T40final_NoSalish_[2004-2007]

---

## Cascadia Runs [total 11 TB]

**perigee: /data2/parker/archive/cascadia**

- [B2002, C2003-2009]
- [B2002, C2003-2009]_lp Low-passed versions have little gaps at the year start and end because of the way the filtering was done using one-year folders.
- Cdia2005 Exactly like C2005, but with diagnostics and averages turned on. [5 TB]
- also some _lp versions

These are are an 8-year sequence of quasi-continuous runs (first PM ran each year on its own to get a good IC for the next year, then he re-ran the last 7 years). They are physics-only, no dye or NPZDO. Hourly files, about 714 GB per year.

---

## LiveOcean Runs (results stored in day folders like f2020.01.26)

**apogee: /dat1/parker/LO_roms/**

- cas7_t0_x4b This is the **CURRENT FORECAST** Late 2012-present. [## TB and growing]
- **OBSOLETE:** cas6_traps2_x2b This was the current forecast 2024.07.01-2024.08.25. It lands here before being archived into cas6_v1_live on perigee [1 TB]
- wgh2_t0_xn0b Nested model Willapa Bay and Grays Harbor, running as part of the daily forecast, 2023.08.07-now. [238 GB and growing]
- cas2k_v0_x2b A low-resolution version of cas6 (2 km grid) designed for long coastal runs. [2 TB]
- Various other nested or analytical runs: ae0, ai0, hc0, so0.

**perigee: /data1/parker/LO_roms/**

- **DELETED 2024.08.25** cas6_v1_live **CURRENT FORECAST ARCHIVE** 2023.05.21-now [9 TB and growing]
- cas6_v0_live **OLD FORECAST ARCHIVE**, 2016.12.15-2023.05.23 (6+ years). It was originally the files in cas6_v3_lo8b, and then we started appending files from cas6_v0_u0kb (or occasionally _u0mb) after the transition from the LiveOcean code to the LO code. [40 TB, compressed NetCDF]

**perigee: /data2/parker/LiveOcean_roms/output**

- cascadia1_base_lobio5 The original LiveOcean runs on the cascadia1 grid. [15 TB] 2013.01.01-2018.11.08
- aestus1_A1_ae1 An idealized estuary used (I think) in the JPO Variance paper. Three months long. Has averages and diagnostics. [297 GB] 2013.01.01-2013.03.31
- cas6_v3_lo8da A re-run of cas6_v3_lo8b for two months, but with diagnostics and averages saved, and no biology [21 TB] Downwelling: 2017.03.01-2017.04.01, Upwelling: 2018.05.15-2018.06.16
- sj0_v0_lo8nest This is a model of the San Juan Islands, physics only, 100 m resolution, one-way nested inside of cas6_v3_lo8b. It was created for the FHL 2019 Summer School. [1.2 TB] 2019.07.22-2019.12.25

**boiler: /data1/parker/LiveOcean_roms/output**
- NOTE: This machine is deprecated: don't expect any of these files to be available
- cas6_v3_lo8b This was the original cas6 forecast, 500 m grid in the Salish Sea, 45 rivers, full NPZDO+Carbon. [15 TB per year] 2021.01.01-late 2021. This was copied to cas6_v0_live on perigee, described above, and is no longer active.
- cas6_v3_lo8dye Has dye in Hood Canal [2.4 TB] 2019.07.03-2019.10.04

---
