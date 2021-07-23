#SVM Loops - This is the "sibling" file for 200GAMs.R and Lasso.R (it runs the same loops,
#just for the SVM code which was originally put in the SVM.R file).


require(caret)
require(glmnet)
require(data.table)
require(penalizedSVM)
require(splitstackshape)
require(dplyr)

require(broom)
require(mgcv)
require(sjPlot)
require(kableExtra)
require(car)
require(popbio)



singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")
sp_svm <- singleplants


#-------------------------------------------------------------------------------------

#MODEL 4 


v <- c(1:500)
svm1 <- list()
work <- list()
work2 <- list()
work3 <- list()
svmpreds <- list()
svmtablecheck <- list()
resultssvm <- list()
svmtables <- list()
kapp4 <- list()


for (i in v) {
  
  set.seed(i)
  
  trainingsvm <- stratified(sp_svm, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  trainingsvm[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_svm <- anti_join(sp_svm, trainingsvm, by=c("Accepted Symbol"))
  check_svm[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_svm <- check_svm[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  for_svm <- trainingsvm[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                             53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  use_for_svm <- within(for_svm, Status <- relevel(Status, ref="Native"))
  
  fitControl <- trainControl(method = "none")
  
  svm1[[i]] <- train(Status~.,
                     data=use_for_svm,
                     method = "svmLinear",
                     trControl = fitControl,
                     preProcess= c("center", "scale"))
  
  
  
  work[[i]] <- varImp(svm1[[i]], scale = FALSE)
  work2[[i]] <- work[[i]]$importance
  work3[[i]] <- as.vector(work2[[i]]$Invasive)
  
  svmpreds[[i]] <- predict(svm1[[i]], newdata = check_svm)
  
  
  svmtablecheck[[i]] <- confusionMatrix(data = svmpreds[[i]], check_svm$Status)
  
  resultssvm[[i]] <- confusionMatrix(data = svmpreds[[i]], check_svm$Status)
  
  svmtables[[i]] <- confusionMatrix(data = svmpreds[[i]], check_svm$Status)$table
  
  kapp4[[i]] <-resultssvm[[i]]$overall["Kappa"]
  
  message('Running Model ', i, ' of 500')
}  


results_svm_df <- data.frame(term=c("Air_Results","ports_Results", "soc_Results", "awc_Results", 
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



for(i in v){
  results_svm_df <- bind_cols(results_svm_df, work3[[i]])
}

results_svm_df <- results_svm_df[, -1]  
results_svm_df

this <- as.list(rowMeans(results_svm_df))
fwrite(this, "this.csv")

mean.list(svmtables)

mean(colMeans(as.data.frame(kapp4)))

#-----------------------------------------------------------------------------------------------------
#MODEL 7

w81 <- c(1:500)
svm81 <- list()
wor81 <- list()
wor82 <- list()
wor83 <- list()
svmpreds81 <- list()
svmtablecheck81 <- list()
resultssvm81 <- list()
spec81 <- list()
sens81 <- list()
OvAcc81 <-list()
kapp81 <-list()
tablesonly81 <-list() 


sps_svm81 <- sp_svm[!sp_svm$Status == "native"]

for (i in w81) {
  
  set.seed(i)
  
  trainingsvm81 <- stratified(sps_svm81, group = "Status",
                              select = list(Status = c("established","invasive")),
                              size = c(86), replace=FALSE)
  
  trainingsvm81[,Status:=factor(Status)]
  check_svm81 <- anti_join(sps_svm81, trainingsvm81, by=c("Accepted Symbol"))
  check_svm81[,Status:=factor(Status)]
  
  check_svm81 <- check_svm81[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                                 53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  for_svm81 <- trainingsvm81[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                                 53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  
  use_for_svm81 <- within(for_svm81, Status <- relevel(Status, ref="established"))
  check_svm81 <- within(check_svm81, Status <- relevel(Status, ref="established"))
  
  fitControl <- trainControl(method = "none")
  
  svm81[[i]] <- train(Status~.,
                      data=use_for_svm81,
                      method = "svmLinear",
                      trControl = fitControl,
                      preProcess= c("center", "scale"))
  
  
  
  wor81[[i]] <- varImp(svm81[[i]], scale = FALSE)
  wor82[[i]] <- wor81[[i]]$importance
  wor83[[i]] <- as.vector(wor82[[i]]$invasive)
  
  svmpreds81[[i]] <- predict(svm81[[i]], newdata = check_svm81)
  
  
  svmtablecheck81[[i]] <- confusionMatrix(data = svmpreds81[[i]], check_svm81$Status)
  
  resultssvm81[[i]] <- confusionMatrix(data = svmpreds81[[i]], check_svm81$Status)
  
  spec81[[i]] <- resultssvm81[[i]]$byClass["Specificity"]
  sens81[[i]] <- resultssvm81[[i]]$byClass["Sensitivity"]
  OvAcc81[[i]] <-resultssvm81[[i]]$overall["Accuracy"]
  kapp81[[i]] <-resultssvm81[[i]]$overall["Kappa"]
  
  tablesonly81[[i]] <-  confusionMatrix(data = svmpreds81[[i]], check_svm81$Status)$table
  
  message('Running Model ', i, ' of 500')
}  


results_svm_df81 <- data.frame(term=c("Air_Results","ports_Results", "soc_Results", "awc_Results", 
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


for(i in w81){
  results_svm_df81 <- bind_cols(results_svm_df81, wor83[[i]])
}

results_svm_df81 <- results_svm_df81[, -1]  
results_svm_df81

this81 <- as.list(rowMeans(results_svm_df81))
fwrite(this81, "this81.csv")


#Average Confusion Matrix for All Models
mean.list(tablesonly81)

#Average Sensitivity (classes swapped)
mean(colMeans(as.data.frame(spec81)))

#Average Specificity (classes swapped)
mean(colMeans(as.data.frame(sens81)))

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc81)))

#Average Kappa
mean(colMeans(as.data.frame(kapp81)))

#-----------------------------------------------------------------------------------------------------
#MODEL 12

v91 <- c(1:500)
svm91 <- list()
wo91 <- list()
wo92 <- list()
wo93 <- list()
svmpreds91 <- list()
svmtablecheck91 <- list()
resultssvm91 <- list()
spec91 <- list()
sens91 <- list()
OvAcc91 <-list()
kapp91 <-list()
tablesonly91 <-list() 



for (i in v91) {
  
  set.seed(i)
  
  trainingsvm91 <- stratified(sp_svm, group = "Status",
                              select = list(Status = c("established","invasive", "native")),
                              size = c(86), replace=FALSE)
  
  trainingsvm91[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_svm91 <- anti_join(sp_svm, trainingsvm91, by=c("Accepted Symbol"))
  check_svm91[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  
  check_svm91 <- check_svm91[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                                 53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  for_svm91 <- trainingsvm91[,-c(1, 2, 5, 9, 10, 11, 16, 19, 20, 22, 24, 27, 34, 35, 45, 48, 50,
                                 53, 55, 59, 60, 61, 62, 74, 75, 78, 84, 88, 89, 90, 92, 94, 96)]
  
  use_for_svm91 <- within(for_svm91, Status <- relevel(Status, ref="Native"))
  
  fitControl <- trainControl(method = "none")
  
  svm91[[i]] <- train(Status~.,
                      data=use_for_svm91,
                      method = "svmLinear",
                      trControl = fitControl,
                      preProcess= c("center", "scale"))
  
  
  
  wo91[[i]] <- varImp(svm91[[i]], scale = FALSE)
  wo92[[i]] <- wo91[[i]]$importance
  wo93[[i]] <- as.vector(wo92[[i]]$Invasive)
  
  svmpreds91[[i]] <- predict(svm91[[i]], newdata = check_svm91)
  
  
  svmtablecheck91[[i]] <- confusionMatrix(data = svmpreds91[[i]], check_svm91$Status)
  
  resultssvm91[[i]] <- confusionMatrix(data = svmpreds91[[i]], check_svm91$Status)
  
  spec91[[i]] <- resultssvm91[[i]]$byClass["Specificity"]
  sens91[[i]] <- resultssvm91[[i]]$byClass["Sensitivity"]
  OvAcc91[[i]] <-resultssvm91[[i]]$overall["Accuracy"]
  kapp91[[i]] <-resultssvm91[[i]]$overall["Kappa"]
  
  tablesonly91[[i]] <-  confusionMatrix(data = svmpreds91[[i]], check_svm91$Status)$table
  
  message('Running Model ', i, ' of 500')
}  


results_svm_df91 <- data.frame(term=c("Air_Results","ports_Results", "soc_Results", "awc_Results", 
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


for(i in v91){
  results_svm_df91 <- bind_cols(results_svm_df91, wo93[[i]])
}

results_svm_df91 <- results_svm_df91[, -1]  
results_svm_df91

this91 <- as.list(rowMeans(results_svm_df91))
this91
fwrite(this91, "this91.csv")

#Average Confusion Matrix for All Models
mean.list(tablesonly91)

#Sensitivity and Specificity calcualted by hand in 
#excel since there are so many classes and I didn't
#trust R. 

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc91)))

#Average Kappa
mean(colMeans(as.data.frame(kapp91)))
