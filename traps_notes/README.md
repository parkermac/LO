# LO_traps
 Module to add tiny rivers and point sources (TRAPS) to LiveOcean.

 ---
## Description

LO_traps is a module that contains all data and scripts required to add TRAPS to LiveOcean. Tiny rivers are defined as smaller rivers that do not already exist in LiveOcean. Point sources are anthropogenic point sources such as wastewater treatment plants and factories.

The TRAPS integration module is a work in progress. If you encounter any bugs or have general feedback, please email auroral@uw.edu. Thanks!

All data and source locations have been downloaded from Washington State Department of Ecology's [website](https://fortress.wa.gov/ecy/ezshare/EAP/SalishSea/SalishSeaModelBoundingScenarios.html). These data are also used in the Salish Sea Model.

---
## Getting set up

<details><summary><strong>Where to put which files</strong></summary>

To enable TRAPS, you will need to move the files in LO_traps to the correct directory within the LO system (i.e your LO_user).

First, clone the LO_traps repo onto your computer so you can pull updates easily. Note that you will still need to manually copy files from your instance of LO_traps into your instance of LO_user. Specifically:

- Copy the LO_traps/user/pre/trapsV00 directory into your LO_user/pre directory
- Copy the LO_traps/user/forcing/trapsV00 directory into your LO_user/forcing directory

So you should have:
- LO_user/pre/trapsV00
- LO_user/forcing/trapsV00

</details>

<details><summary><strong>Additional required files</strong></summary> 

The data used to generate TRAPS forcing is stored on Perigee.

On Perigee, copy the /data1/auroral/LO_data/trapsD00* folder into LO_data on your computer and whichever machine you will use to generate forcing (Perigee or Apogee).

Once this is complete you should have an LO_data/trapsD00 folder with the following files:
- **LiveOcean_SSM_rivers.xlsx:** Excel sheet with list of duplicate rivers in LiveOcean and the Salish Sea Model. When you create TRAPS climatology and when you generate forcing, the scripts will look at this excel sheet to determine which rivers to omit from LiveOcean. This ensures that TRAPS does not add duplicate rivers to LiveOcean.
- **wwtp_open_close_dates.xlsx:** Excel sheet with a list of WWTPs and the year that they closed or opened.
- **all_nonpoint_source_data.nc**: Ecology's timeseries data of state variables and lat/lon coordinates for all river mouths. Used in LO_traps/user/pre/trapsP## to generate climatology files.
- **all_point_source_data.nc:** Ecology's timeseries data of state variables and lat/lon coordinates for all point sources. Used in LO_traps/user/pre/trapsP## to generate climatology files.

*Note: trapsD00 is the current version of Ecology data used to generate traps climatologies and forcing. If new data becomes available, I will increment the version number to be trapsD01.

</details>

<details><summary><strong>Update get_lo_info.py</strong></summary> 

The last required step is to update your LO_user/get_lo_info.py file to specify a name for your traps code. In this repo, the default name is "trapsV00." See example of my get_lo_info.py:

![traps_name](https://github.com/ajleeson/LO_user/assets/15829099/2e18508c-c10f-4d1a-a1c4-e9e75897095f)

</details>

---
## Running TRAPS

After getting the required files, users should be able to add TRAPS to their model runs. An overview of the TRAPS workflow is shown below. Scroll past the figure for more detailed steps.

<details><summary><strong>TRAPS workflow diagram</strong></summary>

![traps-top-level-diagram-v4](https://github.com/ajleeson/LO_user/assets/15829099/610263e8-80e4-459d-bc4e-cbf69f98f918)

</details>

<details><summary><strong>Run steps</strong></summary>

<details><summary>1. Generate climatologies</summary>
    
This step generates climatology files for each of the TRAPS.
From your remote machine in LO_user/pre/trapsV00 in ipython:

```
run make_climatology_tinyrivs.py
run make_climatology_pointsources.py
run make_climatology_LOrivbio.py 
```

Climatology pickle files will be generated and saved in three folders in LO_output/pre/trapsV00:

- **point_sources:** Climatology files for point sources
- **tiny_rivers:** Climatology files for tiny rivers
- **LO_rivbio:** Climatology files for pre-existing LO rivers
  
If you want to look at climatology timeseries, run with ```-test True``` on your local machine. This option will create a subfolder in LO_output/pre/trapsV00/[source type]/lo_base/Data_historical/climatology_plots with a climatology figure for each source. An example figure for Burley Creek is shown below.

![Burley Cr](https://github.com/ajleeson/LO_user/assets/15829099/adc0456f-f855-4428-82c5-63f5aa1fa5b0)

</details>

<details><summary>2. Map TRAPS to the grid</summary>

This step uses the lat/lon coordinates of TRAPS to map each source to the nearest appropriate grid cell. Tiny rivers are mapped to the nearest coastal grid cell. Point sources are mapped to the nearest water cell. From your remote maching in LO_user/pre/trapsV00 in ipython:

```
run traps_placement.py -g [gridname]
```

Csv files with river directions and grid indices for the sources will be generated and saved in LO_data/grid/[gridname]

To look at where the TRAPS get mapped, run with run with ```-test True``` on your local machine. This option will create an interactive figure that you can zoom into. And example screenshot is shown below.

![traps-placements](https://github.com/ajleeson/LO_user/assets/15829099/9cb89ea3-1372-48e6-bddc-e0a979385b8e)

</details>

<details><summary>3. Generate TRAPS forcing</summary>

This step generates a rivers.nc files with forcing for all pre-existing LO rivers and TRAPS. It uses the climatologies generated in Step 1, and the grid indices and river directions generated in Step 2.

From your remote machine in LO/driver:

```
python driver_forcing3.py -g [gridname] -r backfill -s new -0 2017.01.01 -1 2017.01.02 -f trapsV00
```

</details>

<details><summary>4. Run the model</summary>

Before running the model, make sure that you enable vertical sources in your dot in file. To do this, update the boolean option in your dot in file so:

```
LwSrc = T
```

This change is necessary because point sources are implemented as vertical sources using the LwSrc module.

After completing this change, run the model as you normally would.

</details>
</details>

---
## Exceptions and nuances in the code

<details><summary><strong>Fraser river NO4</strong></summary>

Ammonium (NO4) climatology generated from Ecology's data for the Fraser River is a constant value of 0.074 mmol/m3. This concentration is lower than I expected. Since the Fraser River is so large, it is important to get this value right. I reached out to Susan Allen at UBC to learn what NO4 concentration her group uses for the Fraser. She recommended a constant concentration of 4.43 mmol/m3 which is the mean measurement from Environmental Canada ([Olson et al, 2020](https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2019JC015766&file=jgrc24099-sup-0001-Text_SI-S01.pdf)).

The 4.43 mmol/m3 NO4 concentration is implemented as an ```if``` statement in the depths of LO_traps/user/forcing/trapsV00/make_LOriv_forcing.py code.

![fraser-nh4-code](https://github.com/ajleeson/LO_user/assets/15829099/353472de-8444-48e6-a016-8ae12aca7b30)

</details>

<details><summary><strong>Shifted rivers in Hood Canal</strong></summary>

Several Hood Canal rivers in Ecology's data, like Union River, get their flow data from the Big Beef Creek USGS river gage. However, the Big Beef Creek gage became inactive in mid-2012. As a result, from mid-2012 through the end of 2014, river data for these Hood Canal rivers are a copy of prior year data. These copied data also appear to be shifted by 3 months.

To prevent river climatologies from being biased by these shifted, copied data, I have removed data from mid-2012 through the end of 2014 for the affected Hood Canal rivers. This "data cropping" is implemented in LO_traps/user/pre/trapsV00/make_climatology_tinyrivs.py.

An example hydrograph for Union River is shown below before and after the data were cropped.

![union-river-hydrograph](https://github.com/ajleeson/LO_user/assets/15829099/5381807c-d46b-4487-96e4-98724981f95e)

</details>

<details><summary><strong>WWTP open and close dates</strong></summary>

LO_data/trapsV00/wwtp_open_close_dates.xlsx is a user-modifiable sheet with the open and close dates of the WWTPs (with a yearly resolution). The information in this excel sheet is read by the LO_traps/user/forcing/trapsV00/make_wwtp_forcing.py script and turned into a series of ``if`` statements. When the user generates forcing for a year in which a WWTP is closed, then the scripts will still add the WWTP to the model grid. However, the script will set the discharge rate to be 0 m3/s.

</details>

<details><summary><strong>Overlapping sources and the Lake Stevens WWTPs</strong></summary>

For the cas7 grid, several pairs of tiny river and pairs of WWTPs get mapped to the same grid cell (despite having different lat/lon coordinates). These pairs of sources are called "overlapping" sources. To prevent ROMS from getting confused, the forcing scripts consolidate overlapping sources into a single source. The scripts sum the flowrates of both sources, and calculates a weighted average for the other state variables (e.g. temperature) based on flowrate. Even if users are not using the cas7 grid, the TRAPS forcing script will identify and consolidate overlapping sources. Note that this script can only consolidate a pair of overlapping tiny rivers, or a pair of overlapping WWTPs. The script is not able to identify whether a tiny river and WWTP are overlapping. Luckily, this scenario does not occur in the cas7 grid.

The Lake Stevens 001 and Lake Stevens 002 WWTPs overlap on the cas7 grid. However, these WWTPs are never open concurrently-- Lake Stevens 002 opens after Lake Stevens 001 closes. Thus, there is a conditional statement in the LO_traps/user/forcing/make_wwtp_forcing.py script that <i>un-</i>consolidates these WWTPs.

</details>

<details><summary><strong>Willamette River</strong></summary>

Willamette River is included in the Ecology data, and it is not explicitly a duplicate pre-existing LiveOcean river. However, Willamette River discharges into the Columbia River. The Columbia River was pre-existing to LiveOcean, and its USGS gauge is downstream of the Willamette River (meaning that the pre-existing Columbia River already includes contribution from the Willamette). Therefore, the TRAPS code needs to remove the Willamette River from being incorporated into LiveOcean.

This exception is handled in LO_user/pre/trapsV00/make_climatology_tinyrivs.py:

![Willamette](https://github.com/ajleeson/LO_user/assets/15829099/8271fb86-d892-4148-9cc7-8b0bfd2cdb75)

And also in LO_user/pre/trapsV00/traps_placement.py:

![remove_willamette](https://github.com/ajleeson/LO_user/assets/15829099/bbdc8f44-db1c-4734-aac6-fcd8ab4c54a0)

</details>

---
## Update Notes

<details><summary><strong>2023.11.28 update</strong></summary>

**LO integration**

Small changes to folder names and naming conventions to be consistent with the version of TRAPS that is integrated in LO.

</details>

<details><summary><strong>2023.09.04 update</strong></summary>

**Full code refactor**

Several improvements were made to the structure and clarity of the TRAPS code. These changes are intended to enhance TRAPS functionality and to make the code more accessible and readable for users.

</details>

<details><summary><strong>2023.01.05 update</strong></summary>

**Adding TRAPS climatology to pre-existing LO rivers**

Upon request, I have created and generated forcing for pre-existing LiveOcean rivers for which Ecology has data. These are all of the rivers in Ecology's dataset that are duplicates of LiveOcean rivers (and are thus not treated as a tiny river). As part of this update, I have created a new climatology script in LO_traps/pre/traps/make_climatology_LOrivbio.py to generate climatology for these duplicate rivers. I have also created a new folder LO_traps/user/forcing/traps1 with updated versions of make_forcing_main.py and trapsfun.py that use the new climatologies.

There are three confusing parts to the new code, listed below. As a user, it is not necessary to understand the details of these nuances to run the code.

1. Not all pre-existing LiveOcean rivers have a corresponding duplicate in Ecology's dataset. Thus, the file in LO_traps/data/LiveOcean_SSM_rivers.xlsx is frequently used to identify which of the pre-existing rivers do have Ecology data.
2. LiveOcean and Ecology's dataset use different names for the same rivers. Thus, there are several places in the code in which the name must be converted. When reading data and writing forcing for the pre-existing LiveOcean rivers, the LO name must be used. When generating climatology, or using climatology to create forcing, the Ecology/SSM name must be used. The helper function trapsfun.LO2SSM_name helps handle this conversion.
3. Some pre-existing LiveOcean rivers that have a corresponding duplicate in Ecology's dataset have weird values. I call them "weird rivers." Some characteristics include near-zero DO, negative TIC, and zero alkalinity. Rather than using Ecology's data for weird rivers, I have deferred to LiveOcean's default handling of these rivers. Thus, there are places within the code in which I subtract a list of "weird rivers" from the list of pre-existing, duplicate rivers.

</details>