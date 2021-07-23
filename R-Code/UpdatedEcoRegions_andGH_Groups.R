
require(stats)
require(data.table)
require(dplyr)

# Ecoregion and Growth Habit groupings: 

#Separting by Ecoregion ----

Ecos <- fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")

Ecos4 <- Ecos[,-c(2:22, 114:118 )]


Ecos4 <- as.data.frame(Ecos4)

#---------------------

a <- c(3:86)
ER <- list()

for (i in a) {
  ER[[i]] <- Ecos4[Ecos4[, c(i)] == 1,]
}


head(ER[[83]])


#----------------------------
ERs <- list()

for( i in a){
  
  ERs[[i]] <- cbind(ER[[i]], "Ecoregion"=i-2)
  
}


head(ERs[[83]])
head(ERs[[32]])

#---------------------------------

Just_ERs <- lapply(ERs, function(x) x[(names(x) %in% c("Status", "sci_name", "Ecoregion"))])

head(Just_ERs)

work <- bind_rows(Just_ERs)

results <-table(work$Status, work$Ecoregion)

fwrite(results, "ER_fig.csv")


#------------------------------------------

#DO IT AGAIN WITH GROWTH HABITS


Ecos <- fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")

Ecos4 <- Ecos[,-c(2:22, 114:118 )]
Ecos4 <- Ecos4[,-c(3:86)]

Ecos4 <- as.data.frame(Ecos4)

names(Ecos4)[3] <- 1
names(Ecos4)[4] <- 2
names(Ecos4)[5] <- 3
names(Ecos4)[6] <- 4
names(Ecos4)[7] <- 5
names(Ecos4)[8] <- 6

#---------------------

a <- c(3:8)
GH <- list()

for (i in a) {
  GH[[i]] <- Ecos4[Ecos4[, c(i)] == 1,]
}


head(GH[[3]])


#----------------------------

GHs <- list()

for( i in a){
  
  GHs[[i]] <- cbind(GH[[i]], "GrowthHabit"=i-2)
  
}

head(GHs[[4]])

#------------------------------

Just_GHs <- lapply(GHs, function(x) x[(names(x) %in% c("Status", "sci_name", "GrowthHabit"))])

head(Just_GHs)

work <- bind_rows(Just_GHs)

results <-table(work$Status, work$GrowthHabit)


