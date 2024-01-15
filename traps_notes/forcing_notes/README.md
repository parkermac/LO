# LO/forcing/trapsF##

This folder contains forcing directories used to generate TRAPS forcing. Every forcing directory should contain some verions of the scripts detailed below.
#test
---
## Forcing folders

<details><summary><strong>trapsF01 (2024.01.05)</strong></summary>

trapsF01 is identical to trapsF00 (realistic conditions), except that there are no nutrients coming from WWTP effluent. i.e. NO3 and NH4 concentrations are set to zero for all point sources (all w-sources). This includes all WWTPs in both US and Canadian waters, and also factories as well. In all cases, the factories had zero or low DIN dischage to begin with. Thus, these changes to the factory nutrient discharges are not expected to make a significant difference relative to changes in WWTP discharge.

</details>

<details><summary><strong>trapsF00 (2023.12.27)</strong></summary>

This version of TRAPS forcing was updated as part of a larger TRAPS upgrade which made the code base more modular.

trapsF00 is the baseline version of TRAPS forcing that most accurately represents realistic conditions.

</details>

<details><summary><strong>trapsV00 (2023.09.04)</strong></summary>

This is the first, fully functional, version of TRAPS. The code implements WWTPs as vertical sources (and they no longer blow up!).

The code in trapsV00 has also been cleaned, commented, and restructured as part of a major refactoring effort.

</details>

---
## make_forcing_main.py & make_[source type]_forcing.py

These scripts generate a rivers.nc file for ROMS. The file is saved in LO_output/forcing/[gridname]/[fdate]/trapsV00.

<details><summary><strong>Summary</strong></summary>

The `make_forcing_main.py` script is based off of typical `make_forcing_main.py` scripts used in LiveOcean. They both generate flow, temperature, salt, and biogeochemistry forcing for sources. What makes the new `make_forcing_main.py` different is that the actual forcing generation for pre-existing LO rivers, tiny rivers, and point sources are all handled separately in three different helper scripts:

- make_LOriv_forcing.py
- make_triv_forcing.py
- make_wwtp_forcing.py

The scripts have been separated to improve readability. Now, `make_forcing_main.py` simply calls each of these helper scripts and concatenates their results into one dataset. The final dataset is saved as rivers.nc.

</details>

<details><summary><strong>Important details</strong></summary>

*Notes for make_LOriv_forcing.py*

- The flow and temperature data for all pre-existing LO rivers is unchanged compared to prior versions of LO.
- There are several pre-exsiting LO rivers for which Ecology also has data. The biogeochemistry variables for these duplicate pre-existing rivers are thus filled using the TRAPS climatology based on Ecology's data (LO_user/pre/trapsV00/make_climatology_LOrivbio.py).
- Some duplicate rivers have weird values in Ecology's dataset (i.e. zero DO, negative TIC, etc.). The algorithm opts to **not** use Ecology's data for these weird rivers, and instead leave these pre-existing rivers unchanged.
- Fraser river NH4 is set to a constant 4.43 mmol/m3 concentration, as recommended by Susan Allen.

*Overlapping rivers in make_triv_forcing.py & make_wwtp_forcing.py*

- Sometimes, a pair of tiny rivers or a pair of WWTPs may be mapped to the same cell on the model grid. They are 'overlapping' sources.
- To ensure that ROMS does not get confused, the forcing algorithm consolidates the overlapping sources into a single source.
- The names of the overlapping sources are combine using a '+'. For instance, the tiny rivers 'Perry Cr' and 'McLane Cr' get combined into a single river called 'Perry Cr+McLane Cr'
- The flowrate of the consolidate source is the sum of the two sources
- The other variables are consolidated using a weighted average based on flowrate

*WWTP open and close dates*

- LO_data/trapsV00/wwtp_open_close_dates.xlsx contains a list of WWTPs and their open/close years
- The `make_wwtp_forcing.py` script checks this file. If a WWTP is closed for the year in which forcing is being generated, then the discharge rate is padded with zeros

</details>

<details><summary><strong>Disabling TRAPS</strong></summary>

Users can choose to enable either tiny rivers, point sources, or both by toggling the logical switches on lines 37 and 38 of `make_forcing_main.py`.

![TRAPS-switch](https://github.com/ajleeson/LO_user/assets/15829099/734d31f4-e240-4506-a875-b25c0c2a9fe9)

*NOTE:* If you enable point sources, then you must also enable LwSrc in the corresponding BLANK.in file. LuvSrc will already be enabled by default because rivers introduce horizontal (or u- v-) momentum to the system. Point sources discharge vertically (w-momentum), so LwSrc must be set to 'T' true. Example screenshot below.

 ![enable_lwsrc](https://user-images.githubusercontent.com/15829099/209903422-4f3f238b-68f8-44e4-b31d-2448cc5d9053.png)

 </details>

---
## trapsfun.py

This script has similar functionality to rivfun.py. It contains helper functions used to generate forcing for TRAPS.

---
## rivfun.py

This script was borrowed directly from Parker's riv00 forcing. It is included in TRAPS because it contains necessary functions to generate forcing for pre-existing LO rivers.