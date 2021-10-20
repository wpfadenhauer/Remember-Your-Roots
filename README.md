# Remember-Your-Roots

The following files (are going to) include the raw data and R code used in the "Remember Your Roots" manuscript. 
My apologies for the confusing file names, hopefully this helps to explain what is in each one. The model numbers referenced below correspond to Table 2 from the manuscript. 
Please reach out if you have any questions about content or organization in this repo. 

<br/><br/>
__Raw Data Files:__
1. SinglePlantsALLVARS2.csv 
    * This file has a single value for all the variables for each plant species. 
    * These values were calcualted from the native area inhabitated by each species in the lower 48 United States.  
    * So one row per species, with all the biogeographic variables associated with that species also included in the same row.
    * It also includes the nativity status and growth habit for each plant species.    

<br/><br/>
__R Code Files:__ 
1. CodeForPub2.R
    * Calculate average value for each invasion status group for each variable
    * Then test for sig. diffs. between those averages by running Kruskal-Wallis and Wilcoxon test for all variables
2. 500GAMS.R
    * Runs first set of logistic regressions (GAMS), aka models 2 & 6. 
3. Lasso500.R
    * Runs second set of logistic regressions (GLMs w/ LASSO var. selection), aka models 3 & 8. 
4. Multinom.R
    * Runs third and final set of logistic regressions (multinomial - one is a GAM, one is a GLM), aka models 10 & 11. 
5. 500SVMs.R
    * Runs all SVMs, aka models 4, 7, & 12. 
6. 500RandomForests.R
    * Runs all Random Forests, aka models 1, 5, & 9. 
7. UpdatedEcoRegions_andGH_Groups.R
    * Looks at patterns in invasiveness across the individual ecoregions and growth habits
    * This information is used in Figure 4 in the manuscript. 
8. NewSamplingRedo.R
    * This has all of the same models included in the above files, but using a different sampling method (taking 131 from all groups, then 86 of the first sample for training and the remaining 45 for testing).  

  
