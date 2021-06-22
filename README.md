# Remember-Your-Roots

The following files are the code I used to analyze biogeographic traits in the native ranges of invasive plant species. 
This README file is largley a duplicate of the MetadataOfRFilesforManuscript.R file. 
My apologies for the confusing file names, hopefully this helps to explain what is in each one. 

<br/><br/>
__Raw Data Files:__
1. PlantsCSVForR6.csv
    * This is a list of all species that are native to only the Lower 48 (from USDA PLANTS, updated June 2021) and each L48 county each species is native to.
    * In other words, each row is a unique species - county combination.  
2. countyoutlines.shp (& associated files)
    * This is a shapefile with the outlines of all CONUS counties.
    * Note that Louisiana has an outline here but that USDA PLANTS does not contain county-level data for Louisiana, so it is excluded from our analyses. 
3. (Rasters from OrgRasts folder)
    * These are all raw rasters downloaded from the various sources cited in the manuscript. 
4. PlCoMerge.csv (created and subsequently used in CodeForPub.R file)
5. SinglePlants.csv (created and subsequently used in CodeForPub.R file) 
    * This compiled a single value for all the variables across all the counties each species was found within. 
    * So one row per species, with all the biogeographic variables associated with that species also included in the same row.   

<br/><br/>
__R Code Files:__
1. Redoplants.R
    *  Compiles spatial ranges for each species
    *  Calculates the values of the biogeographic traits within each species' spatial range
    *  Requires PlantsCSVForR6.csv, countyoutlines2.shp, all rasters from OrgRasts 
2. CodeForPub.R
    * Calculates average values for all traits by invasion status (used in mansucript figure)
    * Runs Kruskal-Wallis and Wilcoxon test for all variables
3. 200GAMS.R
    * Runs first set of logistic regressions (GAMS)
4. Lasso.R
    * Runs second set of logistic regressions (GLMs w/ LASSO var. selection)
5. Multinom.R
    * Runs third and final set of logistic regressions (multinomial ones)
6. SVMloops.R
    * Runs all SVMs
7. RandomForest.R
    * Runs all Random Forests
8. EcoRegions_and_GrowthHabit_Groupings.R
    * Looks at pattersn in invasiveness across the individual ecoregions and growth habits
    * This information is used in 2 figures in the manuscript
  
