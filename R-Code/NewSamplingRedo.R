#Redo everything with new groups -
#Take 131 from Native, Established, and Invasive. 

#Then use 86 of those 131 for training and remaining 45 for testing. 

#This way the results won't show huge bias due to native species. 

require(data.table)
require(splitstackshape)
require(dplyr)
require(mgcv)
require(caret)
require(car)
require(mgcv.helper)
require(popbio)
require(AER)
require(tidyr)
require(broom)


singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS2.csv")
singleplants$Range_Size <- as.numeric(singleplants$Range_Size)
singleplants <- as.data.frame(singleplants)
singleplants[c(3:9, 11:21)] <- scale(singleplants[c(3:9, 11:21)], center = TRUE, scale = TRUE)
sp_gam <- as.data.table(singleplants)
head(sp_gam)


#---------------------------------------------------------------------------------------
#MODEL 2
a <- c(1:500)
model <- list()
sum <- list ()
Chisq <- list()
results <- list()

spec1 <- list()
sens1 <- list()
OvAcc1 <-list()
kapp1 <- list()
tablesonly1 <- list()


#For-loop that selects a different group of species each time and then runs all of the 
#non-binary variables through the model to determine which ones are most often used
#to predict status. 
for (i in a) {
  
  set.seed(i)
  
  traininggam2 <- stratified(sp_gam, group = "Status",
                             select = list(Status = c("established","invasive", "native")),
                             size = c(131), replace=FALSE)
  
  traininggam <- stratified(traininggam2, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  traininggam[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_gam <- anti_join(traininggam2, traininggam, by=c("Scientific Name"))
  check_gam[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_gam <- check_gam[,-1]
  for_gam <- traininggam[,-1]
  use_for_gam <- within(for_gam, Status <- relevel(Status, ref="Native"))
  
  model[[i]] <- gam(Status~ s(Air_Results, bs="cs") + 
                      s(ports_Results, bs="cs") + 
                      s(soc_Results, bs="cs") +
                      s(awc_Results, bs="cs") + 
                      s(PropFP, bs="cs") + 
                      s(Ag, bs="cs") + 
                      s(Urb, bs="cs") + 
                      s(SpeciesArea, bs="cs") + 
                      s(Range_Size, bs="cs") + 
                      s(Precip, bs="cs")+
                      s(pH, bs="cs")+
                      s(clay, bs="cs"),
                    data = use_for_gam,
                    family = binomial,
                    method = "REML")
  
  sum[[i]] <- summary(model[[i]])
  Chisq[[i]] <- sum[[i]]$chi.sq
  
  probabilitiesgam1 <- predict.gam(model[[i]], type="response", newdata = check_gam)
  class_predictionsgam1 <- ifelse(probabilitiesgam1 < 0.50, "Native", "Invasive")
  final_checkgam1 <- cbind(check_gam, class_predictionsgam1)
  
  tablecheckgam1 <- confusionMatrix(data = as.factor(final_checkgam1$class_predictionsgam1),
                                    reference = as.factor(final_checkgam1$Status),
                                    prevalence = 2)
  
  #Note that the PPV and NPV from the confusion matrix function are wrong - need to calc by hand
  results[[i]] <- confusionMatrix(data = as.factor(final_checkgam1$class_predictionsgam1),
                                  reference = as.factor(final_checkgam1$Status),
                                  prevalence = 2)
  
  spec1[[i]] <- results[[i]]$byClass["Specificity"]
  sens1[[i]] <- results[[i]]$byClass["Sensitivity"]
  OvAcc1[[i]] <-results[[i]]$overall["Accuracy"]
  kapp1[[i]] <-results[[i]]$overall["Kappa"]
  
  tablesonly1[[i]] <- confusionMatrix(data = as.factor(final_checkgam1$class_predictionsgam1),
                                      reference = as.factor(final_checkgam1$Status),
                                      prevalence = 2)$table
  
  message('Running Model ', i, ' of 500')
}



results_df <- data.frame(term=c("Air_Results","ports_Results", "soc_Results",
                                "awc_Results", "PropFP", "Ag", "Urb", "SpeciesArea",
                                "Range_Size", "Precip", "pH", 
                                "clay"))


for(i in a){
  results_df <- bind_cols(results_df, Chisq[[i]])
}


#Average Confusion Matrix for All Models
mean.list(tablesonly1)

#Average Sensitivity (I know the formula is Spec,
#it's backwards for this cuz Est. is positive when it should be reference class)
mean(colMeans(as.data.frame(spec1)))

#Average Specificity (I know the formula is sens,
#it's backwards for this)
mean(colMeans(as.data.frame(sens1)))

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc1)))
sd(colMeans(as.data.frame(OvAcc1)))


#Average Kappa
mean(colMeans(as.data.frame(kapp1)))
sd(colMeans(as.data.frame(kapp1)))


#This is the number of times the `mgcv` package shrunk the variable 
#influence to 0 (Out of 400 models run, each with a different combination
#of species). 
modcount400 <- rowSums(results_df>0.001)
modcount400

plot(model[[5]], pages =1)

mean.list(results)#done, in Excel file


#MODEL ASSUMPTIONS: 

#1. Check concurvity (instead of VIF, since it's a GAM)
concurvity(model[[2]], full = FALSE)
#Nothing is consistently above 0.8 (observed), for multiple models, so we're ok. 

#2. Dispersion parameter
#"Overdispersion is irrelevant for models that estimate a scale parameter"
# according to Ben Bolker (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion)
#So it doesn't apply to our Binomial GAMs

#3. Linearity

resid<- residuals(model[[2]], type = "response")
linpred <-napredict(model[[2]], model[[2]]$linear.predictors)

scatter.smooth(linpred, resid)
#Smooth line is reasonably close to y=0, so we're ok. 

#4. Influential Outliers

model.data <- augment(model[[2]]) %>% 
  mutate(index = 1:n()) 

model.data$.cooksd > 0.5
#Nothing over 0.5, so we're ok

#-----------------------------------------------------
#MODEL 6

#Established Vs. Invasive

a <- c(1:500)
model2 <- list()
sum2 <- list ()
Chisq2 <- list()
results2 <- list()

spec2 <- list()
sens2 <- list()
OvAcc2 <-list()
kapp2 <- list()
tablesonly2 <- list()


sp_gam2 <- sp_gam[!sp_gam$Status == "native"]

#For-loop that selects a different group of species each time and then runs all of the 
#non-binary variables through the model to determine which ones are most often used
#to predict status. 
for (i in a) {
  
  set.seed(i)
  
  traininggam2 <- stratified(sp_gam2, group = "Status",
                             select = list(Status = c("established","invasive")),
                             size = c(131), replace=FALSE)
  
  traininggam <- stratified(traininggam2, group = "Status",
                            select = list(Status = c("established","invasive")),
                            size = c(86), replace=FALSE)
  
  traininggam[,Status:=factor(Status)]
  check_gam <- anti_join(traininggam2, traininggam, by=c("Scientific Name"))
  check_gam[,Status:=factor(Status)]
  check_gam <- check_gam[,-1]
  for_gam <- traininggam[,-1]
  use_for_gam <- within(for_gam, Status <- relevel(Status, ref="established"))
  check_gam <- within(check_gam, Status <- relevel(Status, ref="established"))
  
  model2[[i]] <- gam(Status~ s(Air_Results, bs="cs") + 
                       s(ports_Results, bs="cs") + 
                       s(soc_Results, bs="cs") +
                       s(awc_Results, bs="cs") + 
                       s(PropFP, bs="cs") + 
                       s(Ag, bs="cs") + 
                       s(Urb, bs="cs") + 
                       s(SpeciesArea, bs="cs") + 
                       s(Range_Size, bs="cs") + 
                       s(Precip, bs="cs")+
                       s(pH, bs="cs")+
                       s(clay, bs="cs"),
                     data = use_for_gam,
                     family = binomial,
                     method = "REML")
  
  sum2[[i]] <- summary(model2[[i]])
  Chisq2[[i]] <- sum2[[i]]$chi.sq
  
  probabilitiesgam2 <- predict.gam(model2[[i]], type="response", newdata = check_gam)
  class_predictionsgam2 <- ifelse(probabilitiesgam2 < 0.50, "established", "invasive")
  final_checkgam2 <- cbind(check_gam, class_predictionsgam2)
  
  tablecheckgam2 <- confusionMatrix(data = as.factor(final_checkgam2$class_predictionsgam2),
                                    reference = as.factor(final_checkgam2$Status))
  
  #Note that the PPV and NPV from the confusion matrix function are wrong - need to calc by hand
  results2[[i]] <- confusionMatrix(data = as.factor(final_checkgam2$class_predictionsgam2),
                                   reference = as.factor(final_checkgam2$Status))
  
  spec2[[i]] <- results2[[i]]$byClass["Specificity"]
  sens2[[i]] <- results2[[i]]$byClass["Sensitivity"]
  OvAcc2[[i]] <-results2[[i]]$overall["Accuracy"]
  kapp2[[i]] <-results2[[i]]$overall["Kappa"]
  
  tablesonly2[[i]] <- confusionMatrix(data = as.factor(final_checkgam2$class_predictionsgam2),
                                      reference = as.factor(final_checkgam2$Status))$table
  
  message('Running Model ', i, ' of 500')
}



results_df2 <- data.frame(term=c("Air_Results","ports_Results", "soc_Results", "awc_Results",
                                 "PropFP", "Ag", "Urb", "SpeciesArea",
                                 "Range_Size", "Precip", "pH", 
                                 "clay"))


for(i in a){
  results_df2 <- bind_cols(results_df2, Chisq2[[i]])
}


#Average Confusion Matrix for All Models
mean.list(tablesonly2)

#Average Sensitivity (I know the formula is Spec,
#it's backwards for this cuz Est. is positive when it should be reference class)
mean(colMeans(as.data.frame(spec2)))

#Average Specificity (I know the formula is sens,
#it's backwards for this)
mean(colMeans(as.data.frame(sens2)))

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc2)))
sd(colMeans(as.data.frame(OvAcc2)))

#Average Kappa
mean(colMeans(as.data.frame(kapp2)))
sd(colMeans(as.data.frame(kapp2)))

#This is the number of times the `mgcv` package shrunk the variable 
#influence to 0 (Out of 400 models run, each with a different combination
#of species). 
modcount5002 <- rowSums(results_df2>0.001)
modcount5002


plot(model2[[5]], pages =1)


#MODEL ASSUMPTIONS: 

#1. Check concurvity
concurvity(model2[[3]], full = FALSE)
#Nothing is consistently above 0.8 (observed), for multiple models, so we're ok. 


#2. Dispersion parameter
#"Overdispersion is irrelevant for models that estimate a scale parameter"
# according to Ben Bolker (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion)
#So it doesn't apply to our Binomial GAMs

#3. Linearity

resid<- residuals(model2[[3]], type = "response")
linpred <-napredict(model2[[3]], model2[[3]]$linear.predictors)

scatter.smooth(linpred, resid)
#Smooth line is reasonably close to y=0, so we're ok. 

#4. Influential Outliers

model.data <- augment(model2[[2]]) %>% 
  mutate(index = 1:n()) 

model.data$.cooksd > 0.5
#Nothing over 0.5, so we're ok

#----------------------------------------------------------------------------
# Multinomial Logistic Regressions


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
require(methods)

#-----------------------------------------------

#----------------------------------------------

#MODEL 11

singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS2.csv")
singleplants$Range_Size <- round(as.numeric(singleplants$Range_Size))
singleplants <- as.data.frame(singleplants)
singleplants[c(3:9, 11:21)] <- scale(singleplants[c(3:9, 11:21)], center = TRUE, scale = TRUE)
sp_gam <- as.data.table(singleplants)
head(sp_gam)



#Run 500 multinomial ones with continuous variables. 

a <- c(1:500)
model <- list()
sum <- list ()
Chisq <- list()
results <- list()

spec1 <- list()
sens1 <- list()
OvAcc1 <-list()
kapp1 <- list()
tablesonly1 <- list()
OA6 <- list()


for (i in a) {
  
  set.seed(i)
  
  traininggam2 <- stratified(sp_gam, group = "Status",
                             select = list(Status = c("established","invasive", "native")),
                             size = c(131), replace=FALSE)
  
  
  traininggam <- stratified(traininggam2, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  traininggam[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_gam <- anti_join(traininggam2, traininggam, by=c("Accepted Symbol"))
  check_gam[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_gam <- check_gam[,-1]
  for_gam <- traininggam[,-1]
  use_for_gam <- within(for_gam, Status <- relevel(Status, ref="Native"))
  
  model[[i]] <- gam(list((as.integer(Status)-1)~ s(Air_Results, bs="cs") + 
                           s(ports_Results, bs="cs") + 
                           s(soc_Results, bs="cs") +
                           s(awc_Results, bs="cs") + 
                           s(PropFP, bs="cs") + 
                           s(Ag, bs="cs") + 
                           s(Urb, bs="cs") + 
                           s(SpeciesArea, bs="cs") + 
                           s(Range_Size, bs="cs") + 
                           s(Precip, bs="cs")+
                           s(pH, bs="cs")+
                           s(clay, bs="cs"),
                          ~ s(Air_Results, bs="cs") + 
                           s(ports_Results, bs="cs") + 
                           s(soc_Results, bs="cs") +
                           s(awc_Results, bs="cs") + 
                           s(PropFP, bs="cs") + 
                           s(Ag, bs="cs") + 
                           s(Urb, bs="cs") + 
                           s(SpeciesArea, bs="cs") + 
                           s(Range_Size, bs="cs") + 
                           s(Precip, bs="cs")+
                           s(pH, bs="cs")+
                           s(clay, bs="cs")),
                    data = use_for_gam,
                    family = multinom(K=2),
                    method = "REML")
  
  
  
  sum[[i]] <- summary(model[[i]])
  Chisq[[i]] <- sum[[i]]$chi.sq
  
  probabilitiesgam1 <- predict.gam(model[[i]], type="response", newdata = check_gam)
  
  
  max <- max.col(probabilitiesgam1)
  final_checkgam1 <- cbind(check_gam, max)
  final_checkgam1[,max:=factor(max, labels = c("Native", "Established", "Invasive") )]
  
  results[[i]] <-confusionMatrix(data = as.factor(final_checkgam1$max),
                                 reference = as.factor(final_checkgam1$Status))
  
  #Note that the PPV and NPV from the confusion matrix function are wrong - need to calc by hand
  
  kapp1[[i]] <-results[[i]]$overall["Kappa"]
  OA6[[i]] <-results[[i]]$overall["Accuracy"]
  
  tablesonly1[[i]] <- confusionMatrix(data = as.factor(final_checkgam1$max),
                                      reference = as.factor(final_checkgam1$Status))$table
  
  message('Running Model ', i, ' of 500')
}




results_df <- data.frame(term=c("Air_Results","ports_Results", "soc_Results",
                                "awc_Results", "PropFP", "Ag", "Urb", "SpeciesArea",
                                "Range_Size", "Precip", "pH", 
                                "clay", "Air_Results2","ports_Results2", "soc_Results2",
                                "awc_Results2", "PropFP2", "Ag2", "Urb2", "SpeciesArea2",
                                "Range_Size2", "Precip2", "pH2", 
                                "clay2"))


for(i in a){
  results_df <- bind_cols(results_df, Chisq[[i]])
}

#Average Confusion Matrix for All Models
mean.list(tablesonly1)

#Average Kappa
mean(colMeans(as.data.frame(kapp1)))
sd(colMeans(as.data.frame(kapp1)))

mean(colMeans(as.data.frame(OA6)))
sd(colMeans(as.data.frame(OA6)))


#This is the number of times the `mgcv` package shrunk the variable 
#influence to 0 (Out of 400 models run, each with a different combination
#of species). 
modcount400 <- rowSums(results_df>0.001)
modcount400


#MODEL ASSUMPTIONS:

#1. Check concurvity (instead of VIF, since it's a GAM)
concurvity(model[[2]], full = FALSE)
#Nothing is consistently above 0.8 (observed), for multiple models, so we're ok. 

#2. Dispersion parameter
#"Overdispersion is irrelevant for models that estimate a scale parameter"
# according to Ben Bolker (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion)
#So it doesn't apply to our mutinom GAMs

#3. Linearity

resid<- residuals(model[[2]], type = "deviance")
linpred <-napredict(model[[2]], model[[2]]$linear.predictors)

scatter.smooth(linpred[,1], resid)
scatter.smooth(linpred[,2], resid)
#Smooth lines are reasonably close to y=0, so we're ok. 

#4. Influential Outliers
#Hat values / influence / pearson residuals cannot be calculated for
#multinational GAMs, so unable to determine outliers, but model fit looks ok. 

#--------------------------------------------------------------
#MODEL 10


#Run 500 with categorical variables. 
singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS2.csv")
sp_lasso <- singleplants[,-c(2:22)]


d <- c(1:500)
lassod <- list()
cv.lassod <- list ()
lasso_modeld <- list()
coefd<- list()
countd <- list()
resultsd <- list()
kapp2 <- list()
resultsg <- list()
OA5 <- list()


for (i in d) {
  
  set.seed(i)
  
  traininggam2 <- stratified(sp_lasso, group = "Status",
                             select = list(Status = c("established","invasive", "native")),
                             size = c(131), replace=FALSE)
  
  traininglasso <- stratified(traininggam2, group = "Status",
                              select = list(Status = c("established","invasive", "native")),
                              size = c(86), replace=FALSE)
  
  traininglasso[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_lasso <- anti_join(traininggam2, traininglasso, by=c("Accepted Symbol"))
  check_lasso[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_lasso <- check_lasso[,-c(1,3,6,13,14,24,27,29,32,34,38:41,53,54,57,63,
                                 67,68,69,71,73,75)]
  for_lasso <- traininglasso[,-c(1,3,6,13,14,24,27,29,32,34,38:41,53,54,57,63,
                                 67,68,69,71,73,75)]
  use_for_lasso <- within(for_lasso, Status <- relevel(Status, ref="Native"))
  
  lassod[[i]] <- model.matrix(Status~., use_for_lasso)[,-1]
  
  cv.lassod[[i]] <- cv.glmnet(lassod[[i]], use_for_lasso$Status, alpha =1, family= "multinomial")
  
  lambdad <- cv.lassod[[i]]
  
  lasso_modeld[[i]] <- glmnet(lassod[[i]], use_for_lasso$Status, alpha =1, family = "multinomial", 
                              lambda = lambdad$lambda.1se)
  
  coefd[[i]] <- coef(lasso_modeld[[i]])
  
  check_lasso_matrix <- model.matrix(Status~., check_lasso)[,-1]
  
  
  check<- predict(lasso_modeld[[i]], newx = check_lasso_matrix)
  
  
  class_predictionslasso4 <- as.data.frame(predict(lasso_modeld[[i]], type="class", newx = check_lasso_matrix))
  final_checklasso4 <- as.data.table(cbind(class_predictionslasso4, check_lasso$Status))
  
  
  tablechecklasso4 <- confusionMatrix(data = as.factor(final_checklasso4$s0),
                                      reference = as.factor(final_checklasso4$`check_lasso$Status`))
  
  resultsd[[i]] <- confusionMatrix(data = as.factor(final_checklasso4$s0),
                                   reference = as.factor(final_checklasso4$`check_lasso$Status`))$table
  
  resultsg[[i]] <- confusionMatrix(data = as.factor(final_checklasso4$s0),
                                   reference = as.factor(final_checklasso4$`check_lasso$Status`))
  
  kapp2[[i]] <-resultsg[[i]]$overall["Kappa"]
  OA5[[i]] <-resultsg[[i]]$overall["Accuracy"]
  
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
                                 "FR - III","FR - IV","FR - V", "85"))




for(i in d){
  results_df4 <- bind_cols(results_df4, as.vector(coefd[[i]]$Native))
}

rowSums(!results_df4==0.000)


mean.list(resultsd)

mean(colMeans(as.data.frame(kapp2)))
sd(colMeans(as.data.frame(kapp2)))

mean(colMeans(as.data.frame(OA5)))
sd(colMeans(as.data.frame(OA5)))



#MODEL ASSUMPTIONS: 

#1. Multicollinearity check not necessary, as it does not have large
#influences on LASSO GLMs.

#2. Dispersion parameter
sigma(lasso_modeld[[2]])
#reasonably close to 1, so not concerning. 

#3. Linearity
#Linearity is automatically met since all of our variables for GLMs are
#categorical.

#4. Influential Outliers
#LASSO is relatively robust to outliers, so no need to check. 

#------------------------------------------------------------------------------------

#Lasso 500


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

singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS2.csv")
sp_lasso <- singleplants[,-c(2:22)]



#-------------------------------------------------------------
#MODEL 3

#Run 500 with Native vs. Established and Invasive Combined



d <- c(1:500)
lassod <- list()
cv.lassod <- list ()
lasso_modeld <- list()
coefd<- list()
countd <- list()
resultsd <- list()
OA4 <- list()
kapp2 <- list()
resultsg <- list()



for (i in d) {
  
  set.seed(i)
  
  traininggam2 <- stratified(sp_lasso, group = "Status",
                             select = list(Status = c("established","invasive", "native")),
                             size = c(131), replace=FALSE)
  
  
  traininglasso <- stratified(traininggam2, group = "Status",
                              select = list(Status = c("established","invasive", "native")),
                              size = c(86), replace=FALSE)
  
  traininglasso[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_lasso <- anti_join(traininggam2, traininglasso, by=c("Accepted Symbol"))
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
  OA4[[i]] <-resultsg[[i]]$overall["Accuracy"]
  
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
                                 "FR - III","FR - IV","FR - V", "85"))




for(i in d){
  results_df4 <- bind_cols(results_df4, as.vector(coefd[[i]]))
}

rowSums(!results_df4==0.000)

plot(cv.lasso)
coef(lasso_model)

mean.list(resultsd)

mean(colMeans(as.data.frame(kapp2)))
sd(colMeans(as.data.frame(kapp2)))

mean(colMeans(as.data.frame(OA4)))
sd(colMeans(as.data.frame(OA4)))

#MODEL ASSUMPTIONS: 

#1. Multicollinearity check not necessary, as it does not have large
#influences on LASSO GLMs.

#2. Dispersion parameter
sigma(lasso_modeld[[2]])
#reasonably close to 1, so not concerning. 

#3. Linearity
#Linearity is automatically met since all of our variables for GLMs are
#categorical.

#4. Influential Outliers
#LASSO is relatively robust to outliers, so no need to check. 


#------------------------------------------------------------------


#MODEL 8

#Run 500 with Established and Invasive


e <- c(1:500)
lassoe <- list()
cv.lassoe <- list ()
lasso_modele <- list()
coefe<- list()
counte <- list()
resultse <- list()
OA2 <- list()
kapp3 <- list()
tables3 <- list()

sp_lasso2 <- sp_lasso[!sp_lasso$Status == "native"]

for (i in e) {
  
  set.seed(i)
  
  traininggam2 <- stratified(sp_lasso2, group = "Status",
                             select = list(Status = c("established","invasive", "native")),
                             size = c(131), replace=FALSE)
  
  traininglasso2 <- stratified(traininggam2, group = "Status",
                               select = list(Status = c("established","invasive")),
                               size = c(86), replace=FALSE)
  
  traininglasso2[,Status:=factor(Status)]
  check_lasso2 <- anti_join(traininggam2, traininglasso2, by=c("Accepted Symbol"))
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
  OA2[[i]] <-resultse[[i]]$overall["Accuracy"]
  
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
                                 "FR - III","FR - IV","FR - V", "85"))



for(i in e){
  results_df5 <- bind_cols(results_df5, as.vector(coefe[[i]]))
}


rowSums(!results_df5==0.000)

mean.list(tables3)


mean(colMeans(as.data.frame(kapp3)))
sd(colMeans(as.data.frame(kapp3)))

mean(colMeans(as.data.frame(OA2)))
sd(colMeans(as.data.frame(OA2)))


#MODEL ASSUMPTIONS: 

#1. Multicollinearity check not necessary, as it does not have large
#influences on LASSO GLMs.

#2. Dispersion parameter
sigma(lasso_modele[[4]])
#It's a little high, but this style of model was not particularly useful anyway.

#3. Linearity
#Linearity is automatically met since all of our variables for GLMs are
#categorical.

#4. Influential Outliers
#LASSO is relatively robust to outliers, so no need to check. 



#------------------------------------------------------------------

#Packages needed for SVMs (models 4, 7, 12)

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



singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS2.csv")
sp_svm <- singleplants


#----------------------
#MODEL 4

#Run 500 SVMs with NAT vs ESt + Invasive

v <- c(1:500)
svm1 <- list()
work <- list()
work2 <- list()
work3 <- list()
svmpreds <- list()
svmtablecheck <- list()
resultssvm <- list()
svmtables <- list()
OvAcc4 <- list()
kapp4 <- list()


for (i in v) {
  
  set.seed(i)

  trainingsvm2 <- stratified(sp_svm, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(131), replace=FALSE)
  
  trainingsvm <- stratified(trainingsvm2, group = "Status",
                             select = list(Status = c("established","invasive", "native")),
                             size = c(86), replace=FALSE)
  
  
  trainingsvm[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_svm <- anti_join(trainingsvm2, trainingsvm, by=c("Accepted Symbol"))
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
  
  OvAcc4[[i]] <-resultssvm[[i]]$overall["Accuracy"]
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
                                    "FR - III","FR - IV","FR - V", "85"))



for(i in v){
  results_svm_df <- bind_cols(results_svm_df, work3[[i]])
}

results_svm_df <- results_svm_df[, -1]  
results_svm_df

this <- as.list(rowMeans(results_svm_df))
fwrite(this, "this.csv")

mean.list(svmtables)

mean(colMeans(as.data.frame(kapp4)))
sd(colMeans(as.data.frame(kapp4)))

#Average Overall Accuracy
mean(colMeans(as.data.frame(OvAcc4)))
sd(colMeans(as.data.frame(OvAcc4)))

#-----------------------------------------------------------------------------------------------------

#MODEL 7

#Run 500 SVMs with EST vs. INV

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
  
  trainingsvm82 <- stratified(sps_svm81, group = "Status",
                              select = list(Status = c("established","invasive")),
                              size = c(131), replace=FALSE)
  
  trainingsvm81 <- stratified(trainingsvm82, group = "Status",
                              select = list(Status = c("established","invasive")),
                              size = c(86), replace=FALSE)
  
  trainingsvm81[,Status:=factor(Status)]
  check_svm81 <- anti_join(trainingsvm82, trainingsvm81, by=c("Accepted Symbol"))
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
                                      "FR - III","FR - IV","FR - V", "85"))


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
sd(colMeans(as.data.frame(OvAcc81)))


#Average Kappa
mean(colMeans(as.data.frame(kapp81)))
sd(colMeans(as.data.frame(kapp81)))


#------------------------------------------


#MODEL 12

#Run 500 SVMs with NAT vs. EST vs. INV


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
  
  trainingsvm92 <- stratified(sp_svm, group = "Status",
                              select = list(Status = c("established","invasive", "native")),
                              size = c(131), replace=FALSE)
  
  trainingsvm91 <- stratified(trainingsvm92, group = "Status",
                              select = list(Status = c("established","invasive", "native")),
                              size = c(86), replace=FALSE)
  
  trainingsvm91[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_svm91 <- anti_join(trainingsvm92, trainingsvm91, by=c("Accepted Symbol"))
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
                                      "FR - III","FR - IV","FR - V", "85"))


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
sd(colMeans(as.data.frame(OvAcc91)))

#Average Kappa
mean(colMeans(as.data.frame(kapp91)))
sd(colMeans(as.data.frame(kapp91)))


#-------------------------------------------

#Packages needed for Random Forests

require(data.table)
require(splitstackshape)
require(caTools)
require(randomForest)
require(dplyr)
require(caret)
require(popbio)



singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS2.csv")
sp_rf <- singleplants

#-----------------------------------------------------------------------------

#MODEL 9

#Run 500 Random Forests with NAT vs. EST vs. INV

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
  
  trainingrf1 <- stratified(sp_rf, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(131), replace=FALSE)
  
  trainingrf2 <- stratified(sp_rf, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  check_rf2 <- anti_join(trainingrf1, trainingrf2, by=c("Accepted Symbol"))
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
                                 "FR - III","FR - IV","FR - V", "85"))


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
sd(colMeans(as.data.frame(OvAcc_rf2)))

#Average Kappa
mean(colMeans(as.data.frame(kapp_rf2)))
sd(colMeans(as.data.frame(kapp_rf2)))


#-------------------------------------------------------------------------

#MODEL 1

#Run 500 Random Forests with NAT vs. EST + INV


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
  
  trainingrf4 <- stratified(sp_rf, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(131), replace=FALSE)
  
  trainingrf3 <- stratified(trainingrf4, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  trainingrf3[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  
  check_rf3 <- anti_join(trainingrf4, trainingrf3, by=c("Accepted Symbol"))
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
                                 "FR - III","FR - IV","FR - V", "85"))



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
sd(colMeans(as.data.frame(OvAcc_rf3)))

#Average Kappa
mean(colMeans(as.data.frame(kapp_rf3)))
sd(colMeans(as.data.frame(kapp_rf3)))

#------------------------------------------

#MODEL 5


#Run 500 Random Forests with EST vs. INV


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
  
  trainingrf5 <- stratified(sps_rf, group = "Status",
                            select = list(Status = c("established","invasive")),
                            size = c(131), replace=FALSE)
  
  trainingrf4 <- stratified(trainingrf5, group = "Status",
                            select = list(Status = c("established","invasive")),
                            size = c(86), replace=FALSE)
  
  trainingrf4[,Status:=factor(Status, labels = c("Established", "Invasive") )]
  
  check_rf4 <- anti_join(trainingrf5, trainingrf4, by=c("Accepted Symbol"))
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
                                 "FR - III","FR - IV","FR - V", "85"))


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
sd(colMeans(as.data.frame(OvAcc_rf4)))

#Average Kappa
mean(colMeans(as.data.frame(kapp_rf4)))
sd(colMeans(as.data.frame(kapp_rf4)))

