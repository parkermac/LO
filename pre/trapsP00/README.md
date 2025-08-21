# LO/pre/trapsP00

This folder contains the scripts that generate climatologies for TRAPS and that map TRAPS to the model grid. Details of each can be found below.

---
## ecology_excel2netCDF.py

This script converts all of Ecology's raw data from excel format to netCDF format.

Note: Users do not need to run this script. It was already run during TRAPS development, and users can simply use the output .nc files provided on perigee in auroral's LO_data folder. The script is provided here for completion only.

---
## traps_data_ver.csv

This is a config file that allows for modularity. This file contains the name of the LO_data/trapsD## that the user wants to use for climatology generation.

In theory, different versions of trapsD## allows the user to pick different versions of Ecology data to use, should updates be released in the future. However, there is currently only one traps data version: **trapsD00**

---
## climatology_all_source_types.sh

This shellscript automates the climatology generation process.

To use TRAPS, users need to generate climatology for three TRAPS types: point sources, tiny rivers, and pre-existing LO rivers. Running this script will automatically step through the make_climatology scripts for each of these sources types (so users do not need to independently run three different climatology scripts).

---
## make_climatology scripts

The climatology scripts condense the 1999-2017 Ecology timeseries data into yearly climatologies for each source. This section describes how these scripts work.

<details><summary><strong>Summary</strong></summary>

There are three main climatology scripts:

- `make_climatology_pointsources.py`: Creates climatology files for all point sources using Ecology's data in LO_traps/data/all_point_source_data.nc
- `make_climatology_tinyrivers.py`: Creates climatology files for river mouths using Ecology's data in LO_traps/data/all_nonpoint_source_data.nc. This script does not generate climatology for pre-existing rivers in LiveOcean.
- `make_climatology_LOrivbio.py`: Creates biogeochemistry climatology files for all pre-existing LiveOcean rivers for which Ecology has data in LO_traps/data/all_nonpoint_source_data.nc. This script does not generate climatology for tiny rivers, nor does it generate flowrate or temperature climatology.

These scripts generate climatology pickle files in LO_output/pre/trapsP##/[source type]/lo_base/Data_historical.

Climatologies are generated for the following variables:
- flow (TRAPS only)
- temperature (TRAPS only)
- DO
- NO3
- NH4
- TIC
- TAlk

</details>

<details><summary><strong>Algorithm</strong></summary>

The structure of these scripts are all similar, so they will be explained generally. There are a few nuances in `make_climatology_tinyrivers.py` which are discussed explicitly.

1. First, raw data are read from LO_data/trapsD##/all_nonpoint_source_data.nc and LO_data/trapsD##/all_point_source_data.nc. 

>>> **River notes:** For rivers, the script also reads the list of pre-existing LO rivers from LO_data/trapsD##/LiveOcean_SSM_rivers.xlsx. The pre-existing rivers are omitted from  tinyriver climatology. The non-pre-existing rivers are omitted from LOrivbio climatology.

>>> **Tiny river notes:** In the raw Ecology data, there are several tiny rivers with unrealistic biogeochemistry parameters (i.e. zero DO, negative TIC, etc.). These "weird rivers" are temporarily removed from climatology generation. They are handled separately in Step 4.

2. Then, the script creates empty dataframes for DO, discharge, temperature, NO3, NH4, TIC, and TAlk. For every source, the script then fills these dataframes with the average yearly climatology of the full 1999-2017 timeseries from Ecology. Essentially, climatologies are the "average year" of each source. The script also calculates the standard deviations of these climatologies.

3. Next, the script plots the climatology summary statistics. For every state variable, the script calculates and plots the average climatology profile, the standard deviation, and the min and max climatology values. This plot is saved in LO_output/pre/trapsP##/[source type]/lo_base/Data_historical. An example figure is shown below for tiny rivers.<p style="text-align:center;"><img src="https://github.com/ajleeson/LO_user/assets/15829099/43092aad-d254-4e28-b63f-68c84103f53a" width="800"/><br></p>

4. **Only applies to tiny rivers.** The script then overwrites the biogeochemistry climatologies for "weird rivers" with the average climatology of other rivers calculated in Step 3.

5. Finally, the script saves climatology dataframes as pickle files.

</details>

---
## traps_placement.py & traps_helper.py

This section describes the `traps_placement.py` script. The `traps_placement.py` script and its helper functions in `traps_helper.py` determine where in a model domain TRAPS should be located given their lat/lon coordinates.

<details><summary><strong>Summary</strong></summary>
This is the main function that places TRAPS in the model domain. This script runs the placement function twice: once with an input of 'riv' for tiny rivers, and a second time with an input of 'wwtp' for point sources. The script reads lat/lon coordinates of TRAPS, then decides where to place the TRAPS in the model domain.

This function does not output anything, but it does save .csv files with TRAPS location indices in LO_data/grids/[gridname].

The following subsections provide more details about the placement algorithm and its helper functions.

</details>

<details><summary><strong>Algorithm</strong></summary>

*Tiny Rivers*

1. For each river listed in LO_data/trapsD##/all_nonpoint_source_data.nc, the algorithm first checks if the river is already pre-existing in LiveOcean. If it is pre-existing, then this function does nothing and skips to the next river. If the river is not pre-existing in LiveOcean, then this function proceeds to the next step.
2. This function then feeds the lat/lon coordinates of each river into `traps_helper.get_nearest_coastal_cell_riv` to obtain i,j-indices and direction of the placed river (See the "Tiny River Handling" section below for more details).
3. Finally, this function saves river information in LO_data/grids/[gridname]/triv_info.csv.

*Point Sources*

There are no pre-existing rivers in LiveOcean, nor are there any point sources that discharge to multiple grid cells in the SSM. Thus, point sources are easier to handle than tiny rivers.

1. First, the functions feeds each point source listed in LO_data/trapsD##/all_point_source_data.nc into `get_nearest_coastal_cell_wwtp` to obtain the i,j-indices of the places source (See the "Point Source Handling" section below for more details).
2. Then, this function saves point source information in LO_data/grids/[gridname]/wwtp_info.csv.

</details>

### Tiny River Handling

<details><summary><code><strong>get_nearest_coastal_cell_riv</strong></code></summary>
This function finds the closest coastal grid cell to a river mouth, then returns:

- indices of nearest coatal grid cell to river mouth
- river direction
- number of "rings" away the nearest coastal cell is from the river mouth

To calculate these values, this function follows the following steps:

1. Given river mouth lat/lon coordinates, the algorithm determines in which grid cell the river mouth is originally located in.<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/235255958-37e851f5-820e-4b53-aeca-85c101b7ddc8.png" width="500"/><br></p>

2. Checks whether the starting grid cell is a coastal cell by calling `get_cell_info_riv`. If the starting grid cell is a coastal grid cell, then the function returns the i,j-indices of the cell as well as river direction.

3. If the starting grid cell is not a coastal cell, then the function begins searching in a ring around the starting grid cell. For each cell in the surrouding ring, the function calls `get_cell_info_riv`. If no coastal grid cells are found in the first ring, then the function begins searching the next ring, and so on and so forth until a coastal cell is found.<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/235255959-fb10f648-0d58-4647-a8d0-2ae20e1bbb0b.png" width="500"/><br></p>

4. If one coastal cell is found in a ring, then the function records the coastal cell i,j-indices, the distance from the coastal cell to the river mouth, and the river direction (which are outputs of `get_cell_info_riv`).<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/235255962-fca53e68-2195-4d97-a66c-48962d2d491e.png" width="500"/><br></p>If more than one coastal cell is found in a ring, then information will be recorded for the coastal cell that is nearest to the river mouth.<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/235255966-a26a7d8b-b8b6-41b7-a134-3333a43241ef.png" width="500"/><br></p><br>

Note that this function always checks for a "nearest coastal cell" one ring further out than the first coastal cell-containing ring. This check is important for stretched grids. In a stretched grid, it is possible that the nearest coastal grid cell is located several rings away, even if there are coastal grid cells in closer rings.

<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/235260467-dc0a89b9-5a26-48e5-914e-dd0a87c043da.png" width="450"/><br></p>

</details>

<details><summary><code><strong>get_cell_info_riv</strong></code></summary>

A grid cell of interest is determined in `get_nearest_coatal_cell_riv` before being fed as an input to this function.
This function checks if the grid cell of interest is a coastal water cell. If it is a coastal water cell, then the function returns:

- indices of the coastal grid cell
- distance from the center of the grid cell to the river mouth
- direction of river flow, given the relative position of the nearest land cell

To calculate these values, this function follows the following steps:

1. Checks if the grid cell of interest is a coastal cell by checking whether any adjacent cells have a land mask. The figure below shows a simple domain with a land cell located to the North and East of the grid cell of interest.<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/234992835-2f83a04a-82b2-423a-ae6c-19eff040c75e.png" width="400"/><br></p>

2. If the grid cell is indeed coastal, then the distance from the river mouth to the grid cell is recorded as an output. The function then proceeds to steps 3 and 4. If the grid cell is not coastal, then the function ends and nothing is returned.<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/234995892-1907373b-d2ac-4284-82a3-6efb6d121563.png" width="400"/><br></p>

3. Then the function needs to decide from which land cell the river should flow (i.e. what direction does the river come from?)<br> First, the function calculates the distance from the river mouth to each adjacent land cell. <br> <p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/234995894-d6a13d85-23f7-4d08-ba52-d1e881511c8a.png" width="400"/><br></p> The river flow direction is set by whichever adjacent land cell is closest to the original river mouth lat/lon coordinates. In our simple example, the Eastern land cell is closest to the river. Thus, the function decides that the river mouth flows westward into the grid cell of interest from the eastern land cell. <br> <p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/234995895-3b0f6e21-c479-4e6c-bde2-c2da72f2d0a2.png" width="400"/><br></p>

4. Finally, the function outputs the indices of the coastal grid cell, the distance from the river mouth to the coastal grid cell, and the direction of river flow into the grid cell.

</details>

### Point Source Handling

<details><summary><code><strong>get_nearest_coastal_cell_wwtp</strong></code></summary>

This function is the point source equivalent of `get_nearest_coastal_cell_riv`. The main difference is that this function calls `get_cell_info_wwtp` rather than `get_cell_info_riv`.

The nearest coastal cell that this function is searching for is *any* water cell. This function does not search through rings if the starting grid cell is already a water cell.

This function only needs to search for the nearest coastal grid cell if the starting cell is a land cell.

<p style="text-align:center;"><img src="https://user-images.githubusercontent.com/15829099/235257876-aedcee38-b4b1-4899-a40f-bce06cb7c6ed.png" width="800"/><br></p>

</details>

<details><summary><code><strong>get_cell_info_wwtp</strong></code></summary>

A grid cell of interest is determined in `get_nearest_coatal_cell_wwtp` before being fed as an input to this function.

This function is the point source equivalent of `get_cell_info_riv`, except it is much simpler. In general, point sources are easier to handle than tiny rivers because point sources can be located on an water cell (including in open water), whereas rivers must be located on a land-adjacent water cell. Furthermore, rivers need an associated flow direction, but point sources do not. Thus, this function only needs to check whether the grid cell of interest is a water cell. If so, the function returns the i,j-indices of the grid cell of interest as well as the distance from the center of the grid cell to the point source. If the grid cell of interest is not a water cell, then nothing is returned.
</details><br>