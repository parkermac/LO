# LO/trapsP01

This README describes the processing scripts which converts raw point source (WWTP) and nonpoint source (river) data into a netCDF file for easy incorporation into the LiveOcean model.

---
## Datasets

Two different datasets were combined in this Tiny River And Point Source (TRAPS) iteration.

### Wastewater Treatment Plants (WWTPs)

- [**Mohamedali et al. (2020)**](https://fortress.wa.gov/ecy/ezshare/EAP/SalishSea/SalishSeaModelBoundingScenarios.html)
    - **Dataset Description**: Monthly point source discharge, nutrient loads, temperature for WWTPs (n=89) and industrial facilities (n=10) discharging to both US and Canadian marine waters.
    - **Dataset Timespan**: January 1999 - July 2017
    - **LiveOcean Handling**:
        - Industrial facitlies are omitted from the LiveOcean integration.
        - All WWTPs from this dataset are included in LiveOcean
        - WWTP discharge and nutrient concentrations are updated to values from Wasielewski et al. (2024), if available.

![all_loads_comparison](https://github.com/user-attachments/assets/ada4d134-bec3-43a5-b257-d42628ba4511)

- [**Wasielewski et al. (2024)**](https://www.sciencebase.gov/catalog/item/64762b37d34e4e58932d9d81)
    - **Dataset Description**: Monthly point source nutrient discharge for WWTPs (n=??), industrial facilities (n=??), and fish hatcheries (n=??) discharging to Washington state watersheds.
    
    
    provides monthly data for points sources discharging into Washington state watersheds from January 2005 - December 2020. Only the WWTP data are integrated into LiveOcean, and industrial facility and fish hatchery data are omitted from LiveOcean integration. The Wasieleski et al. (2024) WWTP discharge and nutrient loads are only incorporated in LiveOcean if the WWTP is also present in the Mohamedali et al. (2020) dataset. Any WWTP present in Wasielewski et al. (2024), but not Mohamedali et al. (2020), is omitted from LiveOcean.

Figure 1 depicts the locations of WWTPs included in LiveOcean, and their associated data source.

**TO-DO:** Add figure of WWTP locations 

### Rivers

- [**Mohamedali et al. (2020)**](https://fortress.wa.gov/ecy/ezshare/EAP/SalishSea/SalishSeaModelBoundingScenarios.html), as mentioned above, provides data for both rivers and WWTPs discharging to both US and Canadian waters from. The data span from January 1999 - July 2017. All tiny rivers in LiveOcean use data from this source. All WWTPs from this dataset are included in LiveOcean, but the discharge and nutrient concentrations use values from Wasielewski et al. (2024), if available.

---
## Omitting Industrial Facilities and Fish Hatcheries