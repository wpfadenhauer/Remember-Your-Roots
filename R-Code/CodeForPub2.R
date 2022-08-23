#Code for Pub, take 2


require(data.table)
require(vctrs)
require(stringi)
require(ggpubr)


SpAvgs2 <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")


#Values for Bar Charts in Publication----
#Calculating averages for invasion statuses

#Sorry, this is called ANOVA_Avs_df even though we ended up using Kruskal Wallis instead of ANOVAs. 

#Start with empty dataframe:
ANOVA_Avs_df <- data.frame(Status=factor(c("1","2","3")),
                           GHMI=(c(1,2,3)), 
                           T_Range=(c(1,2,3)), 
                           P_Range=(c(1,2,3)),
                           pH_Range=(c(1,2,3)),
                           airport=(c(1,2,3)),
                           port=(c(1,2,3)),
                           urban=(c(1,2,3)),
                           agric=(c(1,2,3)),
                           AWC=(c(1,2,3)),
                           sand_range=(c(1,2,3)),
                           silt_Range=(c(1,2,3)),
                           clay_Range=(c(1,2,3)),
                           SOC=(c(1,2,3)),
                           Storm=(c(1,2,3)),
                           Tornado=(c(1,2,3)),
                           NatSpecs=(c(1,2,3)),
                           FloodPl=(c(1,2,3)),
                           RangeSz=(c(1,2,3)),
                           stringsAsFactors=FALSE) 

ANOVA_Avs_df$Status <- c("Native", "Established", "Invasive")
ANOVA_Avs_df

#Then fill those columns with actual averages (the V1 just prevents including the status over and over)
#Temp
ANOVA_Avs_df$T_Range <- SpAvgs2[,mean(Temp), by=Status]$V1

#Precip
ANOVA_Avs_df$P_Range <- SpAvgs2[,mean(Precip), by=Status]$V1

#Soil pH
ANOVA_Avs_df$pH_Range <- SpAvgs2[,mean(pH), by=Status]$V1

#Airport Distane
ANOVA_Avs_df$airport <- SpAvgs2[,mean(Air_Results), by=Status]$V1

#Port Distance
ANOVA_Avs_df$port <- SpAvgs2[,mean(ports_Results), by=Status]$V1

#GHMI
ANOVA_Avs_df$GHMI <-SpAvgs2[,mean(ghmi_Results), by=Status]$V1

#Urban Land Use
ANOVA_Avs_df$urban <- SpAvgs2[,mean(Urb, na.rm = TRUE), by=Status]$V1

#Agri Land Use
ANOVA_Avs_df$agric <-SpAvgs2[,mean(Ag, na.rm = TRUE), by=Status]$V1

#Proportion of Sand in Soil
ANOVA_Avs_df$sand_range <-SpAvgs2[,mean(sand, na.rm = TRUE), by=Status]$V1

#Proportion of Silt in Soil
ANOVA_Avs_df$silt_Range <-SpAvgs2[,mean(silt, na.rm = TRUE), by=Status]$V1

#Proportion of Clay in Soil
ANOVA_Avs_df$clay_Range <- SpAvgs2[,mean(clay, na.rm = TRUE), by=Status]$V1

#Soil Organic Carbon
ANOVA_Avs_df$SOC <-  SpAvgs2[,mean(soc_Results), by=Status]$V1

#Soil Available Water Content
ANOVA_Avs_df$AWC <- SpAvgs2[,mean(awc_Results), by=Status]$V1

#Species per county
ANOVA_Avs_df$NatSpecs <- SpAvgs2[,mean(SpeciesArea), by=Status]$V1

#Number of Tropical Storms
ANOVA_Avs_df$Storm <- SpAvgs2[,mean(StormArea), by=Status]$V1

#Number of Tornados
ANOVA_Avs_df$Tornado <- SpAvgs2[,mean(TornArea), by=Status]$V1

#Percent Floodplain
ANOVA_Avs_df$FloodPl <- SpAvgs2[,mean(PropFP), by=Status]$V1

#Range Size
ANOVA_Avs_df$RangeSz <- SpAvgs2[,mean(Range_Size), by=Status]$V1


#Export results for graphing! 
fwrite(ANOVA_Avs_df, "ANOVABarGraphs.csv")



#Kruskal Wallis and Wilcoxon Tests----
#Note that this does not include all of the data exploration and ANOVAs with various transformations we tried 
#before settling on Kruskal Wallis and Wilcoxon. 

#Code from this part mainly taken from: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

#Temp ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "Temp", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Temp", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(Temp ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$Temp, SpAvgs2$Status)

#Precip ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "Precip", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Precip", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(Precip ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$Precip, SpAvgs2$Status)

#Soil pH ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "pH", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "pH", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(pH ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$pH, SpAvgs2$Status)

#Airports ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "Air_Results", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Dist. to Airports", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(Air_Results ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$Air_Results, SpAvgs2$Status)

#Ports ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "ports_Results", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Dist. to Ports", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(ports_Results ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$ports_Results, SpAvgs2$Status)

#GHMI ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "ghmi_Results", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "GHMI", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(ghmi_Results ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$ghmi_Results, SpAvgs2$Status)

#Urban Land ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "Urb", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Urban Land", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(Urb ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$Urb, SpAvgs2$Status)

#Agricultural Land ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "Ag", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Agricultural Land", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(Ag ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$Ag, SpAvgs2$Status)

#Sand ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "sand", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Proportion of Sand", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(sand ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$sand, SpAvgs2$Status)

#Silt ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "silt", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Proportion of Silt", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(silt ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$silt, SpAvgs2$Status)

#Clay ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "clay", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Proportion of Clay", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(clay ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$clay, SpAvgs2$Status)

#Soil Organic Content ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "soc_Results", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Soil Organic Content", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(soc_Results ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$soc_Results, SpAvgs2$Status)


#Available Water Content ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "awc_Results", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Available Water Content", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(awc_Results ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$awc_Results, SpAvgs2$Status)
                     
#Other Native Species ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "SpeciesArea", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Other Native Species", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(SpeciesArea ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$SpeciesArea, SpAvgs2$Status)

#Tropical Storms ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "StormArea", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Tropical Storms", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(StormArea ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$StormArea, SpAvgs2$Status)

#Tornadoes ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "TornArea", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Tornadoes", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(TornArea ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$TornArea, SpAvgs2$Status)

#Floodplains ----
ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "PropFP", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Floodplains", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(PropFP ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$PropFP, SpAvgs2$Status)

#Range Size----

SpAvgs2$RS <- as.integer(SpAvgs2$Range_Size/1000000)

ggpubr::ggboxplot(SpAvgs2, x = "Status", y = "RS", 
                  color = "Status", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                  order = c("native", "established", "invasive"),
                  ylab = "Range Size", xlab = "Status")


#Since ANOVA assumptions were violated, use KW and Wilcoxon instead: 
kruskal.test(RS ~ Status, data=SpAvgs2)
pairwise.wilcox.test(SpAvgs2$RS, SpAvgs2$Status)




#Dataset Characteristics Section for Manuscript ----

length(sp_gam$Status[sp_gam$Status == "invasive"])
length(sp_gam$Status[sp_gam$Status == "established"])
length(sp_gam$Status[sp_gam$Status == "native"])

length(sp_gam$Status[sp_gam$`Growth Habit` == "Forb/herb"])
length(sp_gam$Status[sp_gam$`Growth Habit` == "Shrub"])
length(sp_gam$Status[sp_gam$`Growth Habit` == "Subshrub"])
length(sp_gam$Status[sp_gam$`Growth Habit` == "Vine"])
length(sp_gam$Status[sp_gam$`Growth Habit` == "Tree"])
length(sp_gam$Status[sp_gam$`Growth Habit` == "Graminoid"])


datachars <-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/datachars.csv")
length(unique(datachars$GEOID))
counts <- within(datachars, { count <- ave(GEOID, `Accepted Symbol`, FUN=function(x) length(unique(x)))})
counts <- counts[!duplicated(counts[,c("Accepted Symbol")])]
mean(counts$count)
median(counts$count)








