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

#MODEL 10

singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")
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


for (i in a) {
  
  set.seed(i)
  
  traininggam <- stratified(sp_gam, group = "Status",
                            select = list(Status = c("established","invasive", "native")),
                            size = c(86), replace=FALSE)
  
  traininggam[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_gam <- anti_join(sp_gam, traininggam, by=c("Accepted Symbol"))
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


#This is the number of times the `mgcv` package shrunk the variable 
#influence to 0 (Out of 400 models run, each with a different combination
#of species). 
modcount400 <- rowSums(results_df>0.001)
modcount400




#--------------------------------------------------------------


#MODEL 11


#Run 500 with categorical variables. 
singleplants<-fread("E:/UMass/CH3_NativeRangeAnalysis/CH3_NativeRangeAnalysis/SinglePlantsALLVARS.csv")
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



for (i in d) {
  
  set.seed(i)
  
  traininglasso <- stratified(sp_lasso, group = "Status",
                              select = list(Status = c("established","invasive", "native")),
                              size = c(86), replace=FALSE)
  
  traininglasso[,Status:=factor(Status, labels = c("Established", "Invasive", "Native") )]
  check_lasso <- anti_join(sp_lasso, traininglasso, by=c("Accepted Symbol"))
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
  results_df4 <- bind_cols(results_df4, as.vector(coefd[[i]]$Native))
}

rowSums(!results_df4==0.000)


mean.list(resultsd)

mean(colMeans(as.data.frame(kapp2)))




