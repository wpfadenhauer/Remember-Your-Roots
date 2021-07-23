#Lasso 500


#Order of models in this file is as follows:
#1 Native vs. Established/Invasive, but REMOVED >0.75 corr. variables
#2 Established vs. Invasive, but REMOVED >0.75 corr. variables


require(caret)
require(glmnet)
require(data.table)
require(broom)
require(mgcv)
require(splitstackshape)
require(dplyr)
require(sjPlot)
require(kableExtra)
require(car)
require(popbio)

singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")
sp_lasso <- singleplants[,-c(2:22)]




#MODEL 3

#Run 500 with Native vs. Established and Invasive Combined



d <- c(1:500)
lassod <- list()
cv.lassod <- list ()
lasso_modeld <- list()
coefd<- list()
countd <- list()
resultsd <- list()
kapp2 <- list()
resultsg <- list()



for (i in d) {
  
  set.seed(i)
  
  traininglasso <- stratified(sp_lasso, group = "Status",
                              select = list(Status = c("established","invasive", "native")),
                              size = c(86), replace=FALSE)
  
  traininglasso[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_lasso <- anti_join(sp_lasso, traininglasso, by=c("Accepted Symbol"))
  check_lasso[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_lasso <- check_lasso[,-c(1,3,6,13,14,24,27,29,32,34,38:41,53,54,57,63,
                                 67,68,69,71,73,75)]
  for_lasso <- traininglasso[,-c(1,3,6,13,14,24,27,29,32,34,38:41,53,54,57,63,
                                 67,68,69,71,73,75)]
  use_for_lasso <- within(for_lasso, Status <- relevel(Status, ref="Native"))
  
  lassod[[i]] <- model.matrix(Status~., use_for_lasso)[,-1]
  
  cv.lassod[[i]] <- cv.glmnet(lassod[[i]], use_for_lasso$Status, alpha =1, family= "binomial")
  
  lambdad <- cv.lassod[[i]]
  
  lasso_modeld[[i]] <- glmnet(lassod[[i]], use_for_lasso$Status, alpha =1, family = "binomial", 
                              lambda = lambdad$lambda.1se)
  
  coefd[[i]] <- coef(lasso_modeld[[i]])
  
  check_lasso_matrix <- model.matrix(Status~., check_lasso)[,-1]
  
  probslasso4 <- as.data.frame(predict.glmnet(lasso_modeld[[i]], type="response", newx = check_lasso_matrix))
  class_predictionslasso4 <- ifelse(probslasso4 < 0.50, "Native", "Invasive")
  final_checklasso4 <- as.data.table(cbind(class_predictionslasso4, check_lasso$Status))
  final_checklasso4[,V2:=factor(V2, labels = c("Invasive", "Native") )]
  
  tablechecklasso4 <- confusionMatrix(data = as.factor(final_checklasso4$s0),
                                      reference = as.factor(final_checklasso4$V2))
  
  resultsd[[i]] <- confusionMatrix(data = as.factor(final_checklasso4$s0),
                                              reference = as.factor(final_checklasso4$V2))$table
  
  resultsg[[i]] <- confusionMatrix(data = as.factor(final_checklasso4$s0),
                                   reference = as.factor(final_checklasso4$V2))
  
  kapp2[[i]] <-resultsg[[i]]$overall["Kappa"]
  
  message('Running Model ', i, ' of 500')
}      


results_df4 <- data.frame(term=c("Intercept", "2",
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
                          



for(i in d){
  results_df4 <- bind_cols(results_df4, as.vector(coefd[[i]]))
}

rowSums(!results_df4==0.000)

plot(cv.lasso)
coef(lasso_model)

mean.list(resultsd)

mean(colMeans(as.data.frame(kapp2)))



#-----------------------------------------------------------------------------

#MODEL 8

#Run 500 with Established and Invasive


e <- c(1:500)
lassoe <- list()
cv.lassoe <- list ()
lasso_modele <- list()
coefe<- list()
counte <- list()
resultse <- list()
kapp3 <- list()
tables3 <- list()

sp_lasso2 <- sp_lasso[!sp_lasso$Status == "native"]

for (i in e) {
  
  set.seed(i)
  
  traininglasso2 <- stratified(sp_lasso2, group = "Status",
                               select = list(Status = c("established","invasive")),
                               size = c(86), replace=FALSE)
  
  traininglasso2[,Status:=factor(Status)]
  check_lasso2 <- anti_join(sp_lasso2, traininglasso2, by=c("Accepted Symbol"))
  check_lasso2[,Status:=factor(Status)]
  check_lasso2 <- check_lasso2[,-c(1,3,6,13,14,24,27,29,32,34,38:41,53,54,57,63,
                                   67,68,69,71,73,75)]
  for_lasso2 <- traininglasso2[,-c(1,3,6,13,14,24,27,29,32,34,38:41,53,54,57,63,
                                   67,68,69,71,73,75)]
  use_for_lasso2 <- within(for_lasso2, Status <- relevel(Status, ref="established"))
  check_lasso2 <- within(check_lasso2, Status <- relevel(Status, ref="established"))
  
  lassoe[[i]] <- model.matrix(Status~., use_for_lasso2)[,-1]
  
  cv.lassoe[[i]] <- cv.glmnet(lassoe[[i]], use_for_lasso2$Status, alpha =1, family= "binomial")
  
  lambdae <- cv.lassoe[[i]]
  
  lasso_modele[[i]] <- glmnet(lassoe[[i]], use_for_lasso2$Status, alpha =1, family = "binomial", 
                              lambda = lambdae$lambda.1se)
  
  coefe[[i]] <- coef(lasso_modele[[i]])
  
  check_lasso_matrix2 <- model.matrix(Status~., check_lasso2)[,-1]
  
  probslasso2 <- as.data.frame(predict.glmnet(lasso_modele[[i]], type="response", newx = check_lasso_matrix2))
  class_predictionslasso2 <- ifelse(probslasso2 < 0.50, "established", "invasive")
  final_checklasso2 <- as.data.table(cbind(class_predictionslasso2, check_lasso2$Status))
  final_checklasso2[,V2:=factor(V2, labels = c("established", "invasive") )]
  
  tablechecklasso2 <- confusionMatrix(data = as.factor(final_checklasso2$s0),
                                      reference = as.factor(final_checklasso2$V2))
  
  tables3[[i]] <- confusionMatrix(data = as.factor(final_checklasso2$s0),
                                              reference = as.factor(final_checklasso2$V2))$table
  
  resultse[[i]] <- confusionMatrix(data = as.factor(final_checklasso2$s0),
                                  reference = as.factor(final_checklasso2$V2))
  
  kapp3[[i]] <-resultse[[i]]$overall["Kappa"]
  
  message('Running Model ', i, ' of 500')
  
}



results_df5 <- data.frame(term=c("Intercept", "2",
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



for(i in e){
  results_df5 <- bind_cols(results_df5, as.vector(coefe[[i]]))
}


rowSums(!results_df5==0.000)

mean.list(tables3)


mean(colMeans(as.data.frame(kapp3)))
