# LO/trapsP01

This README describes the processing scripts which converts raw point source (WWTP) and nonpoint source (river) data into a netCDF file for easy incorporation into the LiveOcean model.

---
## Datasets

Two different datasets were combined in this Tiny River And Point Source (TRAPS) iteration.


<details><summary><strong>Wastewater Treatment Plants (WWTPs)</strong></summary>

</br>

<details><summary>Mohamedali et al. (2020)</summary>

[**Data source**](https://fortress.wa.gov/ecy/ezshare/EAP/SalishSea/SalishSeaModelBoundingScenarios.html)

- **Dataset Description**: Monthly point source discharge, nutrient loads, temperature for WWTPs (n=89) and industrial facilities (n=10) discharging to both US and Canadian marine waters. Developed by Washington State Department of Ecology for input to the Salish Sea Model.
- **Dataset Timespan**: January 1999 - July 2017
- **LiveOcean Handling**:
    - Industrial facitlies are omitted from the LiveOcean integration.
    - All WWTPs from this dataset are included in LiveOcean
    - WWTP discharge and nutrient concentrations are updated to values from Wasielewski et al. (2024), if available.

Figure 1 depicts locations of all point sources in Mohamedali et al. (2020), and the mean annual dissolved inorganic nitrogen (DIN) load of each source type.

<p style="text-align:center;"><img src="figures/moh20_all_loads_comparison.png" width="430"/><br>Fig 1. Top panel: locations of WWTPs and industrial facilities in the Mohamedali et al. (2020) dataset. Bottom panel: climatology nutrient load profiles for the sum of each type of facilitiy (e.g., pink is the climatology for the sum of all WWTPs).</p><br>

</details>

<details><summary>Wasielewski et al. (2024)</summary>

[**Data source**](https://www.sciencebase.gov/catalog/item/64762b37d34e4e58932d9d81)

- **Dataset Description**: Monthly point source nutrient discharge for WWTPs (n=97), industrial facilities (n=20), and fish hatcheries (n=47) discharging to Washington state watersheds. Developed by Washington State Department of Ecology and United States Geological Survey for input to a SPARROW watershed model.
- **Dataset Timespan**: January 2005 - December 2020
- **LiveOcean Handling**:
    - Industrial facitlies and fish hatcheries are omitted from the LiveOcean integration.
    - WWTPs are only incorporated into LiveOcean if the WWTp is also present in the Mohamedali et al. (2020) dataset. Any WWTP present in Wasielewski et al. (2024), but not Mohamedali et al. (2020), is omitted from LiveOcean.

Figure 2 depicts locations of all point sources in Wasielewski et al. (2024), and the mean annual total nitrogen (TN) load of each source type.
    
<p style="text-align:center;"><img src="figures/was24_all_loads_comparison.png" width="430"/><br>Fig 2. Top panel: locations of WWTPs, industrial facilities, and fish hatcheries in the Wasielewski et al. (2024) dataset. Bottom panel: climatology nutrient load profiles for the sum of each type of facilitiy (e.g., pink is the climatology for the sum of all WWTPs).</p><br>

</details>

<br>

Figure 3 depicts the locations of WWTPs across the two datasets. They are processed as follows:
- WWTPs in <span style="color:red">red</span> are present ONLY in Mohamedali et al. (2020), and they ARE integrated into LiveOcean.
- WWTPs in <span style="color:dodgerblue">blue</span> are present in both datasets, and they ARE integrated into LiveOcean. Their lat/lon data come from Mohamedali et al. (2020), but their discharge and nutrient concentrations come from Wasielewski et al. (2024)
- WWTPs in <span style="color:goldenrod">yellow</span> are present ONLY in Wasielewski et al. (2024), and they ARE NOT integrated into LiveOcean.

<p style="text-align:center;"><img src="figures/wwtp_locations.png" width="430"/><br>Fig 3. Locations of WWTPs across the datasets.</p><br>

</details>

<details><summary><strong>TO-DO: Rivers</strong></summary>

- [**Mohamedali et al. (2020)**](https://fortress.wa.gov/ecy/ezshare/EAP/SalishSea/SalishSeaModelBoundingScenarios.html), as mentioned above, provides data for both rivers and WWTPs discharging to both US and Canadian waters from. The data span from January 1999 - July 2017. All tiny rivers in LiveOcean use data from this source. All WWTPs from this dataset are included in LiveOcean, but the discharge and nutrient concentrations use values from Wasielewski et al. (2024), if available.

</details>

---
## Data Processing Steps

<details><summary><strong>rawdata_2netCDF.py</strong></summary>

This script compiles all of the excel files from Mohamedali et al. (2020) and the csv files from Wasielewski et al. (2024) into three netCDF files. These .nc files are used for later processing in the TRAPS integration workflow.

Inputs:
- The script reads raw data from the two datasets in:
    - LO_data/trapsD01/mohamedali_etal2020
    - LO_data/trapsD01/wasielewski_etal2024
- It also takes in metadata from the excel files located in LO_data/trapsD01:
    - **SSM_source_info.xlsx** contains metadata, and importantly lat/lon coordinates, for the rivers and point sources in Mohamedali et al. (2020)
    - **wwtp_names.xlsx** contains a list of all WWTPs in Mohamedali et al. (2020), with the corresponding names of WWTPs in the Wasielewski et al. (2024) dataset.
    - **LiveOcean_SSM_rivers.xlsx** contains a list of pre-existing rivers in LiveOcean, and their corresponding river names in the Mohamedali et al. (2020) dataset.

Outputs (which are saved in LO_data/trapsD01/processed_data):
- **river_data_mohamedali_etal_2020.nc** contains daily river data from the Mohamedali et al. (2020) dataset
- **wwtp_data_mohamedali_etal_2020.nc** contains daily WWTP data from the Mohamedali et al. (2020) dataset
- **wwtp_data_wasielewski_etal_2024.nc** contains daily WWTP data from the Wasielewski et al. (2024) dataset

Note that the WWTP data in the two .nc files are unique. This script already handles the nuances of cases in which the same WWTP is present in both datasets.

In theory, this script only needs to be run once.
Then, the netCDF files can be referenced to generate climatologies

This script takes about 15 minutes to run on my local machine.

<details><summary>Exceptions and nuances in data processing</summary>

- Mohamedali et al. (2020)
    - omitted industrial facilities
        - BP Cherry Point
        - Conoco Phillips
        - Intalco
        - Kimberly_Clark
        - Nippon Paper
        - Port Townsend Paper
        - Shell Oil
        - Tesoro
        - US Oil & Refining
        - West Rock
    - omitted WWTPs that are also listed in Wasielewski et al. (2024)
        - which are listed in LO_data/trapsD01/wwtp_names.xlsx
    - WWTP data saved in:
        - LO_data/trapsD01/processed_data/wwtp_data_mohamedali_etal_2020.nc
    - river lat/lon are averaged because some river mouths are split across two grid cells in SSM.

- Wasielewski et al. (2024)
    - omitted industrial facilities
    - omitted fish hatcheries
    - omitted WWTPs that are NOT also listed in Mohamedali et al. (2020)
    - Used lat/lon locations from Mohamedali et al. (2020)
        - Names of the same WWTP in both datasets are listed in LO_data/trapsD01/wwtp_names.xlsx
        - Special cases where Wasielewski et al. (2024) used the same name for two different WWTPs:
            - 'Everett Water Pollution Control Facility'
                - ID=WA0024490_Gardner corresponds to Moh20's 'OF-100'
                - ID=WA0024490_Snohomish corresponds to Moh20's 'Everett Snohomish'
            - 'OAK HARBOR STP':
                - ID=WA0020567-001 corresponds to Moh20's 'Oak Harbor RBC', which we omit anyways because it stopped operating in 2010
                - ID=WA0020567-002 corresponds to Moh20's 'Oak Harbor Lagoon'
            <p style="text-align:center;"><img src="figures/everett_and_oakharbor.png" width="430"/><br></p><br>

    - This dataset has flow, nitrate, and ammonium data. but not temp, DO, TIC, and alkalinity
        - used climatology of these variables from Mohamedali et al. (2020) WWTPs as inputs for these WWTPs
            - note that all WWTPs in Mohamedali et al. (2020) uses the same values for all of these variables.
            - was careful about leap years and non-leap years
    - Special case WWTPs:
        - removed Oak Harbor STP (WA0020567-001), which stopped operating in 2010
        - removed Lake Stevens Sewer Disctric (WA0020893-thru2012) and later combined flows with the newer Lake Stevens WWTP
        <p style="text-align:center;"><img src="figures/lake_stevens_handling.png" width="430"/><br></p><br>

        - padded end of Port Gamble WWTP (WA0022292) with zeros, because it was [decommisioned in 2017](https://ecology.wa.gov/blog/june-2017/around-the-sound-ongoing-and-future-restoration-r#:~:text=Decommissioning%20of%20the%20Port%20Gamble,be%20finished%20by%20March%202018.)

</details>

</details>

---
## Testing scripts

Files located in **trapsP01/testing/scripts**

*Note: these scripts are not necessary to the TRAPS workflow. However, they are here as testing and exploratory tools for the WWTPs in the two datasets.*

<details><summary><strong>explore_loading_profiles_before_trapsP01_processing.py</strong></summary>

This script was created prior to writing the processing scripts in trapsP01. The intention of this script was to explore the point source loading data in Mohamedali et al. (2020) and Wasielewski et al. (2024).

The decisions to omit and keep certain WWTPs from the different datasets were directly informed by the analysis in this script.

Output figures from this script are saved to **LO_output/loading_test/point_source_integration**

</details>