# Young_et_al_aerodynamic_resistance_analysis

**Author:** Adam M. Young, Ph.D.  
**Email:** adam.young[at]nau.edu  
**ORCID:** https://orcid.org/0000-0003-2668-2794  

The code provided in this repository was used to conduct the analysis in the following paper:

Young, A.M., et al. 2021. Seasonality in aerodynamic resistance across a range of North American ecosystems. 
Agricultural and Forest Meteorology. DOI: https://doi.org/10.1016/j.agrformet.2021.108613

### /code/ directory:
---------------------
Within this directory you will find the scripts needed to repeat the analysis. The names of the scripts are 
prefixed by a number indicating the chronological order the needed to reproduce the analysis.

There is also a sub-directory (/z_functions/) where custom written functions are stored.

### /ancillary_data/ directory:
---------------------
This folder houses four csv files that store key metadata for running the analysis.

1. mead2_mead3_crop_types.csv - The mead sites change crops about every other year
between soybean and corn. This file just identifies a given crop for a given year.

2. pheno_flux_sites_to_use.csv - This is a table that houses key metadata, such as
site names, lat/lon, elevation, PI measured canopy height, etc.

3. profile_metadata.csv - This is the site-specific metadata needed for the wind-profile
analysis, including information on the measurement height of wind speed from different 
sensors.

4. variables_to_import_for_fluxsites.csv - This file provides the variable names needed
for import from the raw data file download from the AmeriFlux webpage.
