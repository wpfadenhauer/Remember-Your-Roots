# Remember-Your-Roots

The following files include the raw data and R code used in the "Remember Your Roots" manuscript. 
My apologies for the confusing file names, hopefully this helps to explain what is in each one. The model numbers referenced below correspond to Appendix 1: Supplemental Table 1 from the manuscript. 
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
    * Calculate average value for each invasion status group (e.g. "native", "established", and "invasive") for each variable.
    * Then test for significant differences between those averages by running Kruskal-Wallis and Wilcoxon test for all variables.
2. NewSamplingRedo.R
    * This contains all of the models run for the manuscript (see Appendix 1: Supplemental Table 1 for list of models). They are in the order: 2, 6, 11, 10, 3, 8, 4, 7, 12, 9, 1, 5.
    * See Main Text Methods for descriptions of Random Forest models (#1, 5, and 9) and see Appendix 1: Supplemental Methods for descriptions of logistic models and Support Vector Machines (all other models). 
3. UpdatedEcoRegions_andGH_Groups.R
    * Looks at patterns in invasiveness across the individual ecoregions and growth habits. 
    * This information is used in Figure 4 in the manuscript. 
