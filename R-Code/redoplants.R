#Redo Plants
#Calculating variables within each plant's USDA Range

#This is the first R file (there was some slight data wrangling in FME prior to this).
#The SinglePlants.CSV file created at the end of this will be used in subsequent R files. 


require(rgdal)
require(data.table)
require(raster)
require(sf)
require(dplyr)
require(tmap)
require(maptools)
require(exactextractr)
require(adfExplorer)

plants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/databasemgmt_and_data/ShapefileFromFMEWorkflow/PlantsCSVforR6.csv")

#Some prep work ----
#Assigning invasion statues. 1= Native, 2= Established, 
# 3= Invasive. 
plants$Status <- as.factor(ifelse(plants$db_simplified == "NATI", "native",
                                  ifelse(plants$db_simplified == "NATU", "established",
                                         ifelse(plants$db_simplified == "NATU_INV" | plants$db_simplified == "INV", "invasive", 0))))
head(plants)


#Read in Counties Shapefile----
#This came from here: https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
counties = readOGR(file.path("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/county_boundaries/cb_2018_us_county_500k/cb_2018_us_county_500k.shp"))
counties$GEOID_char <- sprintf("%0.5d", as.integer(counties$GEOID))
counties <- st_as_sf(counties)

plants$GEOID_char <- sprintf("%0.5d", as.integer(plants$GEOID))
plants$sciname <- plants$`Scientific Name`

plant_names <- as.vector(unique(plants$sciname))

length(unique(plants$sciname))

#Dividing up the list of plants and their counties into 5 groups. 
plant_names2 <- sort(plant_names)
plants <- plants[order(plants$sciname),]
plant_names2[c(2000)] #"Castilleja genevievana"
plant_names2[c(4000)] #"Eriogonum racemosum"
plant_names2[c(6000)] #"Lotus mearnsii"
plant_names2[c(8000)] #"Piptochaetium avenacioides"
plant_names2[c(10541)] #"Zuckia brandegeei"

#Make plant range shapefiles in 5 batches since the files are big----
#Basically, in each one, I'm gathering the spatial boundaries of all 
#the counties that each speices is found within and binding them
#into a multipolygon. 


#First 20%
plants0.2 <- plants[c(1:71226),]
plant_names0.2 <- plant_names2[c(1:2000)]
plant0.2_multipolygons = lapply(plant_names0.2, function(x) {
  counties %>%
    filter(GEOID_char %in% plants0.2$GEOID_char[plants0.2$sciname == x]) %>% # filter counties
    dplyr::select(geometry) %>%
    summarise() %>% # comment out summarise to get single polygons
    mutate(sciname = x) })

plant0.2_multipolygons = do.call(rbind, plant0.2_multipolygons) # bind them into a data.frame

tm_shape(plant0.2_multipolygons) + tm_polygons() + tm_facets("sciname")


#Second 20%
plants0.4<- plants[c(71227:147503),]
plant_names0.4 <- plant_names2[c(2001:4000)]
plant0.4_multipolygons = lapply(plant_names0.4, function(x) {
  counties %>%
    filter(GEOID_char %in% plants0.4$GEOID_char[plants0.4$sciname == x]) %>% # filter counties
    dplyr::select(geometry) %>%
    summarise() %>% # comment out summarise to get single polygons
    mutate(sciname = x) })

plant0.4_multipolygons = do.call(rbind, plant0.4_multipolygons) # bind them into a data.frame

tm_shape(plant0.4_multipolygons) + tm_polygons() + tm_facets("sciname")


#Third 20%
plants0.6<- plants[c(147504:228458),]
plant_names0.6 <- plant_names2[c(4001:6000)]
plant0.6_multipolygons = lapply(plant_names0.6, function(x) {
  counties %>%
    filter(GEOID_char %in% plants0.6$GEOID_char[plants0.6$sciname == x]) %>% # filter counties
    dplyr::select(geometry) %>%
    summarise() %>% # comment out summarise to get single polygons
    mutate(sciname = x) })

plant0.6_multipolygons = do.call(rbind, plant0.6_multipolygons) # bind them into a data.frame

tm_shape(plant0.6_multipolygons) + tm_polygons() + tm_facets("sciname")


#Fourth 20%
plants0.8<- plants[c(228459:296902),]
plant_names0.8 <- plant_names2[c(6001:8000)]
plant0.8_multipolygons = lapply(plant_names0.8, function(x) {
  counties %>%
    filter(GEOID_char %in% plants0.8$GEOID_char[plants0.8$sciname == x]) %>% # filter counties
    dplyr::select(geometry) %>%
    summarise() %>% # comment out summarise to get single polygons
    mutate(sciname = x) })

plant0.8_multipolygons = do.call(rbind, plant0.8_multipolygons) # bind them into a data.frame

tm_shape(plant0.8_multipolygons) + tm_polygons() + tm_facets("sciname")


#Fifth 20%
plants1<- plants[c(296903:428747),]
plant_names1 <- plant_names2[c(8001:10541)]
plant1_multipolygons = lapply(plant_names1, function(x) {
  counties %>%
    filter(GEOID_char %in% plants1$GEOID_char[plants1$sciname == x]) %>% # filter counties
    dplyr::select(geometry) %>%
    summarise() %>% # comment out summarise to get single polygons
    mutate(sciname = x) })

plant1_multipolygons = do.call(rbind, plant1_multipolygons) # bind them into a data.frame

tm_shape(plant1_multipolygons) + tm_polygons() + tm_facets("sciname")








#Zonal statistics:----

#Bind together the 5 previously made multipolygon lists. 
full_multi_counties <- rbind(plant0.2_multipolygons, plant0.4_multipolygons, plant0.6_multipolygons, plant0.8_multipolygons, plant1_multipolygons)


#Going to export all results, then after I'm done with all the zonal stats,
#I'll import everything back in and bind them all together. But leaving those 
#results in my R environment slows everything down way too much. 

#Airports----
#Load Raster
s <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/AverageTiffsForR/Airports.tif")

#Do zonal stats
Air_Results <- exact_extract(s, full_multi_counties, 'mean')

fwrite(as.data.frame(Air_Results), "air_results.csv")

#Seaports----
#Load Raster
ports <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/AverageTiffsForR/ports.tif")

#Do zonal stats
ports_Results <- exact_extract(ports, full_multi_counties, 'mean')

fwrite(as.data.frame(ports_Results), "ports_results.csv")

#GHMI----
#Load Raster
ghmi_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/AverageTiffsForR/gHM_project.tif")

#Do zonal stats
ghmi_Results <- exact_extract(ghmi_ras, full_multi_counties, 'mean')

fwrite(as.data.frame(ghmi_Results), "ghmi_results.csv")

#OrgCont----
#Load Raster
SOC_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/AverageTiffsForR/SOC_prediction_1991_2010.tif")

#Do zonal stats
soc_Results <- exact_extract(SOC_ras, full_multi_counties, 'mean')

fwrite(as.data.frame(soc_Results), "soc_results.csv")

#AWC----
#Load Raster
AWC_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/AverageTiffsForR/awc_raster.tif")

#Do zonal stats
awc_Results <- exact_extract(AWC_ras, full_multi_counties, 'mean')

fwrite(as.data.frame(awc_Results), "awc_results.csv")



#tempmax----
#Load Raster
tempmax_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Range/tempmax.asc")

#Do zonal stats
tempmax_Results <- exact_extract(tempmax_ras, full_multi_counties, 'max')

fwrite(as.data.frame(tempmax_Results), "tempmax_results.csv")

#tempmin----
#Load Raster
tempmin_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Range/tempmin.asc")

#Do zonal stats
tempmin_Results <- exact_extract(tempmin_ras, full_multi_counties, 'min')

fwrite(as.data.frame(tempmin_Results), "tempmin_results.csv")

#precipmax----
#Load Raster
precip_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Range/precip.asc")

#Do zonal stats
precipmax_Results <- exact_extract(precip_ras, full_multi_counties, 'max')

fwrite(as.data.frame(precipmax_Results), "precipmax_results.csv")

#precipmin----

#Do zonal stats
precipmin_Results <- exact_extract(precip_ras, full_multi_counties, 'min')

fwrite(as.data.frame(precipmin_Results), "precipmin_results.csv")

#pHmax----
#Load Raster
ph_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Range/phfinalraster.tif")

#Do zonal stats
phmax_Results <- exact_extract(ph_ras, full_multi_counties, 'max')

fwrite(as.data.frame(phmax_Results), "phmax_results.csv")

#pHmin----
#Load Raster

#Do zonal stats
phmin_Results <- exact_extract(ph_ras, full_multi_counties, 'min')

fwrite(as.data.frame(phmin_Results), "phmin_results.csv")

#sandmax----
#Load Raster
sand_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Range/avg_sand_raster.tif")

#Do zonal stats
sandmax_Results <- exact_extract(sand_ras, full_multi_counties, 'max')

fwrite(as.data.frame(sandmax_Results), "sandmax_results.csv")

#sandmin----

#Do zonal stats
sandmin_Results <- exact_extract(sand_ras, full_multi_counties, 'min')

fwrite(as.data.frame(sandmin_Results), "sandmin_results.csv")

#siltmax----
#Load Raster
silt_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Range/avg_silt_raster.tif")

#Do zonal stats
siltmax_Results <- exact_extract(silt_ras, full_multi_counties, 'max')

fwrite(as.data.frame(siltmax_Results), "siltmax_results.csv")

#siltmin----

#Do zonal stats
siltmin_Results <- exact_extract(silt_ras, full_multi_counties, 'min')

fwrite(as.data.frame(siltmin_Results), "siltmin_results.csv")

#claymax----
#Load Raster
clay_ras <- raster("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Range/avg_clay_raster.tif")

#Do zonal stats
claymax_Results <- exact_extract(clay_ras, full_multi_counties, 'max')

fwrite(as.data.frame(claymax_Results), "claymax_results.csv")

#claymin----

#Do zonal stats
claymin_Results <- exact_extract(clay_ras, full_multi_counties, 'min')

fwrite(as.data.frame(claymin_Results), "claymin_results.csv")


#Post-Zonal Statistics: -----
#For the remaining variables, the zonal statistics were already 
#done in QGIS, and then the resulting raster histograms / counts are being 
#imported here for futher species-based calculations. 

#The rasters that were used for zonal statistics in QGIS can still
#be found in the OrgRasts folder, and the tool that was used in
#QGIS is quite literally called "Zonal Statistics". 











#Flood Plains ----
floodhist<-fread("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Counts/floodhist.csv")

#For this, 1= Floodplain, 0 = Not floodplain

#Bind county floodplain data to plants object:
plants <- merge(plants, floodhist, sort=FALSE)

#Sum all columns by species
flood_sums <- as.data.table (plants %>%
                  group_by(`Accepted Symbol`,`Scientific Name`) %>%
                    summarise_if(is.integer, sum))

#Remove nonsense columns
flood_sums <- flood_sums[, !"GEOID"]
flood_sums <- flood_sums[, !"HISTO_NODATA"]  

#Calc Prop. FP for each species
flood_sums$PropFP <- (flood_sums$HISTO_1/(flood_sums$HISTO_0 + flood_sums$HISTO_1))

#Remove extra columns
flood_sums <- flood_sums[, !"HISTO_0"]
flood_sums <- flood_sums[, !"HISTO_1"]

#Export
fwrite(as.data.frame(flood_sums), "flood_results.csv")



#Tornados ----
tornhist<-fread("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/CountTracks/torncounts.csv")

counties$ALAND_int <- bit64::as.integer64(counties$ALAND)

#Merge plants, counties, and tornados data
plants <- merge(plants, tornhist, sort=FALSE)
plants <- merge(plants, counties, by="GEOID_char", sort = FALSE)
plants <- plants[, !"geometry"]

#Sum number of tornados and total land area across species
torn_sums <- as.data.table (plants %>%
                               group_by(`Accepted Symbol`,`Scientific Name`) %>%
                               summarise_if(is.numeric, sum))
torn_sums <- torn_sums[, !"ALAND.x"]
torn_sums <- torn_sums[, !"ALAND_char"]

#Calculate tornados per area for each species
torn_sums$TornArea <- (torn_sums$NUMPOINTS/(torn_sums$ALAND_int))

torn_sums <- torn_sums[, !"NUMPOINTS"]
torn_sums <- torn_sums[, !"ALAND_int"]

#Export
fwrite(as.data.frame(torn_sums), "tornado_results.csv")



#Fire Regime (mode) ----
firehist<-fread("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Mode/fire.csv")

#There are five main types of fire regimes classified by LANDFIRE, so that's
#what the 1-5 here are. Note that they represent distinct regimes, so they're factors,
#not integers, hence the mode calculation. 

plants <- merge(plants, firehist, sort=FALSE)

#Sum all columns by species
fire_sums <- as.data.table (plants %>%
                               group_by(`Accepted Symbol`,`Scientific Name`) %>%
                               summarise_if(is.integer, sum))

fire_sums <- fire_sums %>%
  mutate(FRG = case_when(HISTO_5 > HISTO_4 & HISTO_5 > HISTO_3 & HISTO_5 > HISTO_2 & HISTO_5 > HISTO_1  ~ '5',
                         HISTO_4 > HISTO_3 & HISTO_4 > HISTO_2 & HISTO_4 > HISTO_1  ~ '4',
                         HISTO_3 > HISTO_2 & HISTO_3 > HISTO_1 ~ '3',
                         HISTO_2 > HISTO_1 ~ '2',
                         TRUE ~ '1'))

fire_sums <- fire_sums[, !"GEOID"]
fire_sums <- fire_sums[, !"HISTO_1"]
fire_sums <- fire_sums[, !"HISTO_2"]
fire_sums <- fire_sums[, !"HISTO_3"]
fire_sums <- fire_sums[, !"HISTO_4"]
fire_sums <- fire_sums[, !"HISTO_5"]

fwrite(as.data.frame(fire_sums), "fire_results.csv")

#Tropical Storms ----
stormhist<-fread("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/CountTracks/hurrcounts.csv")

counties$ALAND_int <- bit64::as.integer64(counties$ALAND)

plants <- merge(plants, stormhist, sort=FALSE)
plants <- merge(plants, counties, by="GEOID_char", sort = FALSE)

storm_sums <- as.data.table (plants %>%
                              group_by(`Accepted Symbol`,`Scientific Name`) %>%
                              summarise_if(is.numeric, sum))

storm_sums$StormArea <- (storm_sums$CountperGEOID/(storm_sums$ALAND_int))

storm_sums <- storm_sums[, !"CountperGEOID"]
storm_sums <- storm_sums[, !"GEOID.x"]
storm_sums <- storm_sums[, !"STATEFP.x"]
storm_sums <- storm_sums[, !"ALAND_int"]

fwrite(as.data.frame(storm_sums), "tropstorms_results.csv")


#Agricultural land use ----
luhist<-fread("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Counts/histcsv.csv")

#81 and 82 are agricultural land for this. 

plants <- merge(plants, luhist, sort=FALSE)

lu_sums <- as.data.table (plants %>%
                               group_by(`Accepted Symbol`,`Scientific Name`) %>%
                               summarise_if(is.integer, sum))

lu_sums$Ag <- ((lu_sums$HISTO_81 + lu_sums$HISTO_82)/(lu_sums$HISTO_0 + lu_sums$HISTO_11 + lu_sums$HISTO_12
                                                    + lu_sums$HISTO_21 + lu_sums$HISTO_22 + lu_sums$HISTO_23
                                                    + lu_sums$HISTO_24 + lu_sums$HISTO_31 + lu_sums$HISTO_41 +
                                                    lu_sums$HISTO_42 + lu_sums$HISTO_43 + lu_sums$HISTO_52 +
                                                      lu_sums$HISTO_71 + lu_sums$HISTO_90 + lu_sums$HISTO_95 +
                                                      lu_sums$HISTO_81 + lu_sums$HISTO_82))

#Will export with urban after next section. 

#Urban land use ----

#Use same luhist as above.
#21, 22, 23, and 24 are urban land. 

#Use same lu_sums as above

lu_sums$Urb <- ((lu_sums$HISTO_21 + lu_sums$HISTO_22 + lu_sums$HISTO_23 + lu_sums$HISTO_24)/
                                                      (lu_sums$HISTO_0 + lu_sums$HISTO_11 + lu_sums$HISTO_12
                                                      + lu_sums$HISTO_21 + lu_sums$HISTO_22 + lu_sums$HISTO_23
                                                      + lu_sums$HISTO_24 + lu_sums$HISTO_31 + lu_sums$HISTO_41 +
                                                        lu_sums$HISTO_42 + lu_sums$HISTO_43 + lu_sums$HISTO_52 +
                                                        lu_sums$HISTO_71 + lu_sums$HISTO_90 + lu_sums$HISTO_95 +
                                                        lu_sums$HISTO_81 + lu_sums$HISTO_82))

lu_sums <- lu_sums[, !"HISTO_0"]
lu_sums <- lu_sums[, !"HISTO_11"]
lu_sums <- lu_sums[, !"HISTO_12"]
lu_sums <- lu_sums[, !"HISTO_21"]
lu_sums <- lu_sums[, !"HISTO_22"]
lu_sums <- lu_sums[, !"HISTO_23"]
lu_sums <- lu_sums[, !"HISTO_24"]
lu_sums <- lu_sums[, !"HISTO_31"]
lu_sums <- lu_sums[, !"HISTO_41"]
lu_sums <- lu_sums[, !"HISTO_42"]
lu_sums <- lu_sums[, !"HISTO_43"]
lu_sums <- lu_sums[, !"HISTO_52"]
lu_sums <- lu_sums[, !"HISTO_71"]
lu_sums <- lu_sums[, !"HISTO_90"]
lu_sums <- lu_sums[, !"HISTO_95"]
lu_sums <- lu_sums[, !"HISTO_81"]
lu_sums <- lu_sums[, !"HISTO_82"]
lu_sums <- lu_sums[, !"GEOID"]

fwrite(as.data.frame(lu_sums), "urb_and_ag_results.csv")

#Native Species ----

sppercohist<-fread("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/Counts/2021Sp_Per_Co.csv")

sppercohist$GEOID <- sppercohist$FIPS_Code

plants <- merge(plants, sppercohist, sort=FALSE)
plants <- merge(plants, counties, by="GEOID_char", sort = FALSE)

spperco_sums <- as.data.table (plants %>%
                               group_by(`Accepted Symbol`,`Scientific Name`) %>%
                               summarise_if(is.numeric, sum))

spperco_sums$SpeciesArea <- (spperco_sums$count/(spperco_sums$ALAND_int))

spperco_sums <- spperco_sums[, !"GEOID.x"]
spperco_sums <- spperco_sums[, !"FIPS_Code"]
spperco_sums <- spperco_sums[, !"count"]
spperco_sums <- spperco_sums[, !"ALAND_int"]


fwrite(as.data.frame(spperco_sums), "natsps_results.csv")



#Compile all CSVs  ----

air<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/air_results.csv")
seaports <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/ports_results.csv")
ghmi <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/ghmi_results.csv")
soc <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/soc_results.csv")
awc <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/awc_results.csv")
tmax <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/tempmax_results.csv")
tmin <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/tempmin_results.csv")
pmax <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/precipmax_results.csv")
pmin <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/precipmin_results.csv")
phmax <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/phmax_results.csv")
phmin <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/phmin_results.csv")
sdmax <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/sandmax_results.csv")
sdmin <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/sandmin_results.csv")
stmax <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/siltmax_results.csv")
stmin <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/siltmin_results.csv")
cymax <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/claymax_results.csv")
cymin <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/claymin_results.csv")
fp <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/flood_results.csv")
tndo <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/tornado_results.csv")
frg <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/fire_results.csv")
hurr <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/tropstorms_results.csv")
agurb <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/urb_and_ag_results.csv")
natsps <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/natsps_results.csv")


#Need to get the species list back - in the order of full_multi_counties

SinglePlants <- cbind(full_multi_counties$sciname, air, seaports, ghmi, soc, awc, tmax, tmin,
                      pmax, pmin, phmax, phmin, sdmax, sdmin, stmax, stmin, cymax, cymin)

SinglePlants$`Scientific Name` <- SinglePlants$V1

SinglePlants <- merge(SinglePlants, fp, sort=FALSE)
SinglePlants <- merge(SinglePlants, tndo, sort=FALSE)
SinglePlants <- merge(SinglePlants, frg, sort=FALSE)
SinglePlants <- merge(SinglePlants, hurr, sort=FALSE, all.x = TRUE)
SinglePlants <- merge(SinglePlants, agurb, sort=FALSE)
SinglePlants <- merge(SinglePlants, natsps, sort=FALSE)

#need to replace nas with 0 for tropical storms (hurr)

SinglePlants$StormArea[is.na(SinglePlants$StormArea)] = 0

SinglePlants$Temp <- SinglePlants$tempmax_Results - SinglePlants$tempmin_Results
SinglePlants$Precip <- SinglePlants$precipmax_Results - SinglePlants$precipmin_Results
SinglePlants$pH <- SinglePlants$phmax_Results - SinglePlants$phmin_Results
SinglePlants$sand <- SinglePlants$sandmax_Results - SinglePlants$sandmin_Results
SinglePlants$silt <- SinglePlants$siltmax_Results - SinglePlants$siltmin_Results
SinglePlants$clay <- SinglePlants$claymax_Results - SinglePlants$claymin_Results


SinglePlants <- SinglePlants[, !"V1"]
SinglePlants <- SinglePlants[, !"tempmax_Results"]
SinglePlants <- SinglePlants[, !"tempmin_Results"]
SinglePlants <- SinglePlants[, !"precipmax_Results"]
SinglePlants <- SinglePlants[, !"precipmin_Results"]
SinglePlants <- SinglePlants[, !"phmax_Results"]
SinglePlants <- SinglePlants[, !"phmin_Results"]
SinglePlants <- SinglePlants[, !"sandmax_Results"]
SinglePlants <- SinglePlants[, !"sandmin_Results"]
SinglePlants <- SinglePlants[, !"siltmax_Results"]
SinglePlants <- SinglePlants[, !"siltmin_Results"]
SinglePlants <- SinglePlants[, !"claymax_Results"]
SinglePlants <- SinglePlants[, !"claymin_Results"]

#GH & NS ----
#Now, just need to attach growth habit and nativity status the the species, then we're good to go!

ghs<-fread("E:/UMass/CH3_NativeRangeAnalysis/OrgRasts/2021_Growth_Habits.csv")

ghs <- ghs[, !"Synonym Symbol"]
ghs <- ghs[, !"Scientific Name"]
ghs <- ghs[, !"Common Name"]


SinglePlants <- merge(SinglePlants, ghs)

#Some plants have multiple growth habits, but I just one want, so I'll keep the first one
#Essentially randomly choosing one of the growth habits for each sp. 

SinglePlants <- SinglePlants %>% distinct(`Accepted Symbol`, .keep_all = TRUE)

#Nativity Status

Status <- plants
Status <- Status[, !"db_detailed"]
Status <- Status[, !"db_simplified"]
Status <- Status[, !"county"]
Status <- Status[, !"Scientific Name"]
Status <- Status[, !"GEOID"]
Status <- Status[, !"state"]
Status <- Status[, !"StateAb"]
Status <- unique(Status)

SinglePlants <- merge(SinglePlants, Status)

#Export SinglePlants.csv
fwrite(as.data.frame(SinglePlants), "SinglePlants.csv")
