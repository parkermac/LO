# README for the Collias data

---

Run `get_initial_file_info.py` to see info about the data: variable names, units, and spatial and temporal extent.

Years with data are 1932 to 1975, with some gaps

---

Here are the unique values of 'Result_Parameter_Name' and 'Result_Value_Units':
```
                      Salinity [ppt]
    Alkalinity, Total as CaCO3 [mg/L]
               Ortho-Phosphate [mg/L]
                       Nitrate [mg/L]
            Temperature, water [deg C]
              Dissolved Oxygen [mL/L]
         Silicic acid (H4SiO4) [mg/L]
              Nitrogen dioxide [mg/L]
                       Density [mg/L]
```

---

All 89 columns in the original file, with some typical values:
NOTE: to display this first do: pd.set_option('display.max_rows', 500)
```
Study_ID                                                                                    Collias
Study_Name                                                                          Collias Dataset
Location_ID                                                                                  PSB388
Study_Specific_Location_ID                                                                   PSB388
Location_Name                                                                                PSB388
Location_Setting                                                                            Estuary
Field_Collection_Type                                                                   Measurement
Field_Collector                                                                          University
Field_Collection_Start_Date                                                              05/12/1969
Field_Collection_Start_Time                                                                11:18:00
Field_Collection_Start_Date_Time                                                1969-05-12 11:18:00
Field_Collection_End_Date                                                                05/12/1969
Field_Collection_End_Time                                                                  11:18:00
Field_Collection_End_Date_Time                                                  05/12/1969 11:18:00
Field_Collection_Comment                                                                        NaN
Field_Collection_Area                                                                           NaN
Field_Collection_Area_Units                                                                     NaN
Field_Collection_Reference_Point                                                      Water Surface
Field_Collection_Upper_Depth                                                                      0
Field_Collection_Lower_Depth                                                                      0
Field_Collection_Depth_Units                                                                      m
Upper_Depth_In_Feet                                                                             0.0
Lower_Depth_In_Feet                                                                             0.0
Well_Water_Level_Measuring_Point_Or_TOC_ID                                                      NaN
Sample_ID                                                                                  PSB388_1
Sample_Field_Replicate_ID                                                                       NaN
Sample_Replicate_Flag                                                                           NaN
Sample_Sub_ID                                                                                   NaN
Sample_Composite_Flag                                                                           NaN
Storm_Event_Qualifier                                                                           NaN
Sample_Matrix                                                                                 WATER
Sample_Source                                                                     Salt/Marine Water
Sample_Use                                                                                      NaN
Sample_Collection_Method                                                                        NaN
Sample_Collection_Method_Description                                                            NaN
Sample_Preparation_Method                                                                       NaN
Sample_Preparation_Method_Description                                                           NaN
Sample_Method_Other                                                                             NaN
Sample_Method_Other_Type                                                                        NaN
Sample_Method_Other_Description                                                                 NaN
Sample_Taxon_Name                                                                               NaN
Sample_Taxon_Common_Name                                                                        NaN
Sample_Taxon_TSN                                                                                NaN
Sample_Tissue_Type                                                                              NaN
Sample_Percent_Sorted                                                                           NaN
Result_Parameter_Name                                                            Temperature, water
Result_Parameter_CAS_Number                                                                     NaN
Lab_Analysis_Date                                                                               NaN
Lab_Analysis_Date_Accuracy                                                                      NaN
Lab_Analysis_Time                                                                               NaN
Result_Value                                                                                   9.96
Result_Value_Units                                                                            deg C
Result_Reporting_Limit                                                                          NaN
Result_Reporting_Limit_Type                                                                     NaN
Result_Detection_Limit                                                                          NaN
Result_Detection_Limit_Type                                                                     NaN
Result_Data_Qualifier                                                                           NaN
Result_Data_Qualifier_Description                                                               NaN
Result_Suspect_or_Rejected_Flag                                                                 NaN
Result_Suspect_Code                                                                             NaN
Result_Suspect_Description                                                                      NaN
Fraction_Analyzed                                                                               NaN
Field_Filtered_Flag                                                                             NaN
Result_Basis                                                                                    NaN
Digestion_Method                                                                                NaN
Water_Level_Accuracy                                                                            NaN
Result_Method                                                                                TEMPHG
Result_Method_Description                                        Temperature by Mercury Thermometer
Result_Comment                                                                                  NaN
Result_Additional_Comment                                                                       NaN
Result_Lab_Replicate_ID                                                                         NaN
Result_Lab_Name                                                                                 NaN
Result_Validation_Level                                                                         NaN
Result_Taxon_Name                                                                               NaN
Result_Taxon_Common_Name                                                                        NaN
Result_Taxon_Alias_Other                                                                        NaN
Result_Taxon_TSN                                                                                NaN
Result_Taxon_Parent_TSN                                                                         NaN
Result_Taxon_Unidentified_Species                                                               NaN
Result_Taxon_Life_Stage                                                                         NaN
Study_ID_Alias                                                                                  NaN
Study_QA_Planning_Level                                   LEVEL 1: Informal or no QA documentation.
Study_QA_Assessment_Level                         Level 1:  Data neither Verified nor Assessed f...
Study_Type                                                              General environmental study
EIM_Data_Entry_Review_Status                                                           Not Reviewed
Calculated_Latitude_Decimal_Degrees_NAD83HARN                                               47.6115
Calculated_Longitude_Decimal_Degrees_NAD83HARN                                           -122.36122
Record_Created_On                                                               6/7/2021 3:41:51 PM
Result_System_ID                                                                          214410966
```

---

#### Bad Data Issue:

In many of the years there are some stations, often in Hood Canal, where DO is not real, and instead is a duplicate of Salinity:

```
In [51]: run process_data -test True

1939
*** HCB548 There were 1 casts at this station ***
- took 0.2 sec to process
Total processing time 0.0 minutes

In [52]: A
Out[52]:
      salt (ppt)  temp (degC)  DO (mL/L)  SiO4 (mg/L)  NO2 (mg/L)    name                time         lon        lat  cid
z                                                                                                                        
-105       30.59         8.87      30.59        1.882     0.00448  HCB548 1939-04-29 14:00:00 -123.132933  47.388148    0
-75        30.57         9.10      30.57        1.826     0.00182  HCB548 1939-04-29 14:00:00 -123.132933  47.388148    0
-50        30.52         9.14      30.52        1.854     0.00070  HCB548 1939-04-29 14:00:00 -123.132933  47.388148    0
-30        30.43         9.10      30.43        1.826     0.00098  HCB548 1939-04-29 14:00:00 -123.132933  47.388148    0
-20        30.34         9.02      30.34        1.657     0.00252  HCB548 1939-04-29 14:00:00 -123.132933  47.388148    0
-10        30.17         8.78      30.17        1.629         NaN  HCB548 1939-04-29 14:00:00 -123.132933  47.388148    0
 0         29.69         9.24      29.69        1.573         NaN  HCB548 1939-04-29 14:00:00 -123.132933  47.388148    0
```

This is now handled with a mask in `process_data.py`.

NOTE: there are also some instances of way-too-high NO2 that will have to be dealt with.
