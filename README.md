# Remember-Your-Roots

The following files are the code I used to analyze biogeographic traits in the native ranges of invasive plant species. 
This README file is largley a duplicate of the MetadataOfRFilesforManuscript.R file. 
My apologies for the confusing file names, hopefully this helps to explain what is in each one. 

Raw Data Files:
1. PlantsCSVForR6.csv
  * This is all combinations of species - counties with respective variables. 
3. countyoutlines2.shp (& associated files)
  * This is a shapefile with the outlines of all CONUS counties and county-specific variables.
5. (rasters from OrgRasts folder)
  * These are all raw rasters downloaded from the various sources cited in the manuscript. 
7. PlCoMerge.csv (created and subsequently used in CodeForPub.R file)
8. SinglePlants.csv (created and subsequently used in CodeForPub.R file) 
  * This combined all the variables across all the counties each species was found within. So one row per species. 

R Code Files:
1. Redoplants.R
  *  Calculates the values of the biogeographic traits within each species' native range
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

