#Binomial GAMS
#Updated version of 200GAMS.R file


#Order of models in this file is as follows:
#1 Native vs. Established/Invasive
#2 Established vs. Invasive


require(data.table)
require(splitstackshape)
require(dplyr)
require(mgcv)
require(caret)
require(car)
require(mgcv.helper)
require(popbio)


singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")
singleplants$Range_Size <- as.numeric(singleplants$Range_Size)
singleplants <- as.data.frame(singleplants)
singleplants[c(3:9, 11:21)] <- scale(singleplants[c(3:9, 11:21)], center = TRUE, scale = TRUE)
sp_gam <- as.data.table(singleplants)
head(sp_gam)
#---------------------------------------------------------------------------------------
#MODEL 1
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
  
  traininggam <- stratified(sp_gam, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  traininggam[,Status:=factor(Status, labels = c("Invasive", "Invasive", "Native") )]
  check_gam <- anti_join(sp_gam, traininggam, by=c("Scientific Name"))
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

#Average Kappa
mean(colMeans(as.data.frame(kapp1)))


#This is the number of times the `mgcv` package shrunk the variable 
#influence to 0 (Out of 400 models run, each with a different combination
#of species). 
modcount400 <- rowSums(results_df>0.001)
modcount400

plot(model[[5]], pages =1)

#Check concurvity for issues
concurvity(model[[2]], full = FALSE)
#None of the pairwise are particularly high, so not concerned. 

mean.list(results)#done, in Excel file


#---------------------------------------------------------------------------------------
#MODEL 2

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
  
  traininggam <- stratified(sp_gam2, group = "Status",
                            select = list(Status = c("established","invasive")),
                            size = c(86), replace=FALSE)
  
  traininggam[,Status:=factor(Status)]
  check_gam <- anti_join(sp_gam2, traininggam, by=c("Scientific Name"))
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

#Average Kappa
mean(colMeans(as.data.frame(kapp2)))


#This is the number of times the `mgcv` package shrunk the variable 
#influence to 0 (Out of 400 models run, each with a different combination
#of species). 
modcount5002 <- rowSums(results_df2>0.001)
modcount5002


plot(model2[[5]], pages =1)

#Check concurvity for issues
concurvity(model2[[2]], full = FALSE)
#None of the pairwise are particularly high, so not concerned. 

mean.list(results2)#done, in Excel file

