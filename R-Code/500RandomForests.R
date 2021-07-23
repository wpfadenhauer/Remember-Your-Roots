#Random Forest!




require(data.table)
require(splitstackshape)
require(caTools)
require(randomForest)
require(dplyr)
require(caret)
require(popbio)



singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")
sp_rf <- singleplants

#-----------------------------------------------------------------------------

#MODEL 9
#Removed >0.70 corr vars, but still all three invasion statuses

rfss2 <- c(1:500)
rf2 <- list()
predrf2 <- list()
rf_imp2 <- list()

rf_cm2 <- list()
spec_rf2 <- list()
sens_rf2 <- list()
OvAcc_rf2 <-list()
kapp_rf2 <-list()
to_rf2 <-list() 


for (i in rfss2) {
  set.seed(i)
  
  trainingrf2 <- stratified(sp_rf, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  check_rf2 <- anti_join(sp_rf, trainingrf2, by=c("Accepted Symbol"))
  check_rf2 <- check_rf2[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  for_rf2 <- trainingrf2[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  
  use_for_rf2 <- within(for_rf2, Status <- relevel(as.factor(Status), ref ="native"))
  check_rf2 <- within(check_rf2, Status <- relevel(as.factor(Status), ref ="native"))
  
  names(use_for_rf2) <- make.names(names(use_for_rf2))
  names(check_rf2) <- make.names(names(check_rf2))
  
  rf2[[i]] <- randomForest(Status~.,
                           data=use_for_rf2,
                           importance=TRUE)
  
  
  #Tried type 2 (gini impurity) as well, but that's biased towards continous variables.
  rf_imp2[[i]] <- importance(rf2[[i]], type =1)
  
  
  
  predrf2[[i]] <- predict(rf2[[i]], newdata=check_rf2)
  
  rf_cm2[[i]] <- confusionMatrix(data = as.factor(predrf2[[i]]), as.factor(check_rf2$Status))  
  
  spec_rf2[[i]] <- rf_cm2[[i]]$byClass["Specificity"]
  sens_rf2[[i]] <- rf_cm2[[i]]$byClass["Sensitivity"]
  OvAcc_rf2[[i]] <-rf_cm2[[i]]$overall["Accuracy"]
  kapp_rf2[[i]] <-rf_cm2[[i]]$overall["Kappa"]
  to_rf2[[i]] <-confusionMatrix(data = as.factor(predrf2[[i]]), as.factor(check_rf2$Status))$table  
  
  message('Running Model ', i, ' of 500')
  
}  


results_rf2 <- data.frame(term=c("Air_Results","ports_Results", "soc_Results", "awc_Results", 
                                 "PropFP", "Ag", "Urb","SpeciesArea", "Range_Size", "Precip",
                                 "pH", "clay",  "2",
                                 "3", "5", "6", "7", "8","9","10", 
                                 "13", "14", "15", "16", "17", "18", "19", "20", "21",
                                 "23", "24", "26", "28", "29", "31",
                                 "33", "34", "35", "40", "41", "42",
                                 "43", "44", "45", "46", "47", "48", "49", "50", 
                                 "53", "54", "56", "57", "58", "59", "60", "62",
                                 "63", "64",  "68", "70", "72",
                                 "74", "75", "76", "77", "78", "79", "80", "81", "82",
                                 "83", "84", "GH_Subshrub", "GH_Graminoid", "GH_Shrub",
                                 "GH_Vine", "GH_Forb", "GH_Tree","FR - I","FR - II",
                                 "FR - III","FR - IV","FR - V"))


for(i in rfss2){
  results_rf2 <- bind_cols(results_rf2, rf_imp2[[i]])
}

results_rf2 <- results_rf2[, -1]  
results_rf2

this_rf2 <- as.list(rowMeans(results_rf2))
this_rf2
fwrite(this_rf2, "this_rf2.csv")
#fwrite(results_rf2, "listrf2.csv")

#Average Confusion Matrix for All Models
mean.list(to_rf2)

#Sensitivity and Specificity calculated by hand in 
#excel since there are so many classes and I didn't
#trust R. 

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc_rf2)))

#Average Kappa
mean(colMeans(as.data.frame(kapp_rf2)))


#----------------------------------------------------------------
#MODEL 1

rfss3 <- c(1:500)
rf3 <- list()
predrf3 <- list()
rf_imp3 <- list()

rf_cm3 <- list()
spec_rf3 <- list()
sens_rf3 <- list()
OvAcc_rf3 <-list()
kapp_rf3 <-list()
to_rf3 <-list() 


for (i in rfss3) {
  set.seed(i)
  
  trainingrf3 <- stratified(sp_rf, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  trainingrf3[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  
  check_rf3 <- anti_join(sp_rf, trainingrf3, by=c("Accepted Symbol"))
  check_rf3[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_rf3 <- check_rf3[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  for_rf3 <- trainingrf3[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  
  use_for_rf3 <- within(for_rf3, Status <- relevel(as.factor(Status), ref ="Native"))
  check_rf3 <- within(check_rf3, Status <- relevel(as.factor(Status), ref ="Native"))
  
  names(use_for_rf3) <- make.names(names(use_for_rf3))
  names(check_rf3) <- make.names(names(check_rf3))
  
  rf3[[i]] <- randomForest(Status~.,
                           data=use_for_rf3,
                           importance=TRUE)
  
  
  #Tried type 2 (gini impurity) as well, but that's biased towards continous variables.
  rf_imp3[[i]] <- importance(rf3[[i]], type =1)
  
  
  
  
  predrf3[[i]] <- predict(rf3[[i]], newdata=check_rf3)
  
  rf_cm3[[i]] <- confusionMatrix(data = as.factor(predrf3[[i]]), as.factor(check_rf3$Status))  
  
  spec_rf3[[i]] <- rf_cm3[[i]]$byClass["Specificity"]
  sens_rf3[[i]] <- rf_cm3[[i]]$byClass["Sensitivity"]
  OvAcc_rf3[[i]] <-rf_cm3[[i]]$overall["Accuracy"]
  kapp_rf3[[i]] <-rf_cm3[[i]]$overall["Kappa"]
  to_rf3[[i]] <-confusionMatrix(data = as.factor(predrf3[[i]]), as.factor(check_rf3$Status))$table  
  
  message('Running Model ', i, ' of 500')
  
}  


results_rf3 <- data.frame(term=c("Air_Results","ports_Results", "soc_Results", "awc_Results", 
                                 "PropFP", "Ag", "Urb","SpeciesArea", "Range_Size", "Precip",
                                 "pH", "clay",  "2",
                                 "3", "5", "6", "7", "8","9","10", 
                                 "13", "14", "15", "16", "17", "18", "19", "20", "21",
                                 "23", "24", "26", "28", "29", "31",
                                 "33", "34", "35", "40", "41", "42",
                                 "43", "44", "45", "46", "47", "48", "49", "50", 
                                 "53", "54", "56", "57", "58", "59", "60", "62",
                                 "63", "64",  "68", "70", "72",
                                 "74", "75", "76", "77", "78", "79", "80", "81", "82",
                                 "83", "84", "GH_Subshrub", "GH_Graminoid", "GH_Shrub",
                                 "GH_Vine", "GH_Forb", "GH_Tree","FR - I","FR - II",
                                 "FR - III","FR - IV","FR - V"))



for(i in rfss3){
  results_rf3 <- bind_cols(results_rf3, rf_imp3[[i]])
}

results_rf3 <- results_rf3[, -1]  
results_rf3

this_rf3 <- as.list(rowMeans(results_rf3))
this_rf3
fwrite(this_rf3, "this_rf3.csv")
#fwrite(results_rf2, "listrf2.csv")

#Average Confusion Matrix for All Models
mean.list(to_rf3)

#Sensitivity and Specificity calcualted by hand in 
#excel since there are so many classes and I didn't
#trust R. 

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc_rf3)))

#Average Kappa
mean(colMeans(as.data.frame(kapp_rf3)))



#-------------------------------------------------------
#MODEL 5

rfss4 <- c(1:500)
rf4 <- list()
predrf4 <- list()
rf_imp4 <- list()

rf_cm4 <- list()
spec_rf4 <- list()
sens_rf4 <- list()
OvAcc_rf4 <-list()
kapp_rf4 <-list()
to_rf4 <-list() 


sps_rf <- sp_rf[!sp_rf$Status == "native"]

for (i in rfss4) {
  set.seed(i)
  
  trainingrf4 <- stratified(sps_rf, group = "Status",
                            select = list(Status = c("established","invasive")),
                            size = c(86), replace=FALSE)
  
  trainingrf4[,Status:=factor(Status, labels = c("Established", "Invasive") )]
  
  check_rf4 <- anti_join(sps_rf, trainingrf4, by=c("Accepted Symbol"))
  check_rf4[,Status:=factor(Status, labels = c("Established", "Invasive") )]
  check_rf4 <- check_rf4[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  for_rf4 <- trainingrf4[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  use_for_rf4 <- within(for_rf4, Status <- relevel(as.factor(Status), ref ="Established"))
  check_rf4 <- within(check_rf4, Status <- relevel(as.factor(Status), ref ="Established"))
  
  names(use_for_rf4) <- make.names(names(use_for_rf4))
  names(check_rf4) <- make.names(names(check_rf4))
  
  rf4[[i]] <- randomForest(Status~.,
                           data=use_for_rf4,
                           importance=TRUE)
  
  
  #Tried type 2 (gini impurity) as well, but that's biased towards continous variables.
  rf_imp4[[i]] <- importance(rf4[[i]], type =1)
  
  
  predrf4[[i]] <- predict(rf4[[i]], newdata=check_rf4)
  
  rf_cm4[[i]] <- confusionMatrix(data = as.factor(predrf4[[i]]), as.factor(check_rf4$Status))  
  
  spec_rf4[[i]] <- rf_cm4[[i]]$byClass["Specificity"]
  sens_rf4[[i]] <- rf_cm4[[i]]$byClass["Sensitivity"]
  OvAcc_rf4[[i]] <-rf_cm4[[i]]$overall["Accuracy"]
  kapp_rf4[[i]] <-rf_cm4[[i]]$overall["Kappa"]
  to_rf4[[i]] <-confusionMatrix(data = as.factor(predrf4[[i]]), as.factor(check_rf4$Status))$table  
  
  message('Running Model ', i, ' of 500')
  
}  


results_rf4 <- data.frame(term=c("Air_Results","ports_Results", "soc_Results", "awc_Results", 
                                 "PropFP", "Ag", "Urb","SpeciesArea", "Range_Size", "Precip",
                                 "pH", "clay",  "2",
                                 "3", "5", "6", "7", "8","9","10", 
                                 "13", "14", "15", "16", "17", "18", "19", "20", "21",
                                 "23", "24", "26", "28", "29", "31",
                                 "33", "34", "35", "40", "41", "42",
                                 "43", "44", "45", "46", "47", "48", "49", "50", 
                                 "53", "54", "56", "57", "58", "59", "60", "62",
                                 "63", "64",  "68", "70", "72",
                                 "74", "75", "76", "77", "78", "79", "80", "81", "82",
                                 "83", "84", "GH_Subshrub", "GH_Graminoid", "GH_Shrub",
                                 "GH_Vine", "GH_Forb", "GH_Tree","FR - I","FR - II",
                                 "FR - III","FR - IV","FR - V"))


for(i in rfss4){
  results_rf4 <- bind_cols(results_rf4, rf_imp4[[i]])
}

results_rf4 <- results_rf4[, -1]  
results_rf4

this_rf4 <- as.list(rowMeans(results_rf4))
this_rf4
fwrite(this_rf4, "this_rf4.csv")
#fwrite(results_rf2, "listrf2.csv")

#Average Confusion Matrix for All Models
mean.list(to_rf4)

#Sensitivity and Specificity calcualted by hand in 
#excel since there are so many classes and I didn't
#trust R. 

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc_rf4)))

#Average Kappa
mean(colMeans(as.data.frame(kapp_rf4)))

