library(Hmisc)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(pROC)
library(ggfortify)
library(caret)
library(missForest)
library(randomForest)
library(Boruta)
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}
mean_imputation<-function(data,mean_vec,names){
  for(i in 1:ncol(data)){
    data[which(is.na(data[,i])==TRUE),i]<-mean_vec[names[i]]
  }
  return(data)
}
z_normalization<-function(data,mean_feature,std_feature){
  for(i in 1:ncol(data)){
    data[,i]=(data[,i]-mean_feature[i])/std_feature[i]
  }
  return(data)
}
numeric_column<-function(data,column){
  data[,column]<-as.numeric((data[,column]))
}
set.seed(123)
dir.create("/work/path")
setwd("/work/path")
all<-read.table("/path/to/clinvar_hgmd.xls",sep = "\t",stringsAsFactors = F,check.names = F,header = T) 
test1<-read.table("/path/to/paper_VariSNP.xls",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
test1Bs<-read.table("/path/to/VariSNP1100var.xls",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
test1B<-test1[which(test1$class=="B"),]
test1P<-test1[which(test1$class=="P"),]
load("E:/ML/train/mis/RF/version11/1100.118Ftrainvar.rdata")
load("E:/ML/train/mis/RF/version11/rf.65scalefeaturesRMpopTtotal.train1.rdata")
test1Brest<-subset(test1B,!is.element(test1B$AAChange.refGene,test1Bs$AAChange.refGene))
all$class<-ifelse(all$group=="clinvar.B","B","P")
allB<-all[which(all$class=="B"),]
allP<-all[which(all$class=="P"),]
allPs<-all[which(all$group=="clinvar.P" | all$group=="clinvar_hgmd.P.O"),]
allPrest<-subset(allP,!is.element(allP$AAChange.refGene,allPs$AAChange.refGene))
###trainset 
trainset1<-rbind(allPs,allB) 
trainset2<-rbind(test1Bs,test1P)
####Population frequency ###
pop<-trainset1[,c(11:44,157,206,207)]
pop[pop=="."]=0

###prediction ####
pred<-trainset1[,c(45,46,48,49,51,52,54,55,57,58,60,61,63,64,66,67,69,70,72,73,83:84,86:89,114:115,117:119,120,122,123,125:130,208)]
popP<-trainset1[,c(50,53,56,59,62,65,68,71,85)]
popP[popP=="D"]=2
popP[popP=="A"]=2
popP[popP=="T"]=1
popP[popP=="P"]=1
popP[popP=="B"]=1
popP[popP=="M"]=2
popP[popP=="N"]=1
popP[popP=="H"]=2
popP[popP=="L"]=1
popP[popP=="U"]=1.5
popP[popP=="."]=NA
pred[pred=="."]=NA
pred_mean<-apply(pred,2,function(x){mean(c(mean(as.numeric(x[which(trainset1$class=="B")]),na.rm=TRUE),
                                           mean(as.numeric(x[which(trainset1$class=="P")]),na.rm=TRUE)))})
prednames<-colnames(pred)
popP_mean<-apply(popP,2,function(x){mean(c(mean(as.numeric(x[which(trainset1$class=="B")]),na.rm=TRUE),
                                           mean(as.numeric(x[which(trainset1$class=="P")]),na.rm=TRUE)))})
pred<-mean_imputation(pred,pred_mean,prednames)
popPnames<-colnames(popP)
popP<-mean_imputation(popP,popP_mean,popPnames)
pred1<-pred[,c(c(10,12:18,20)*2,41)]
###evolutionary####
evo<-trainset1[,c(131,132,136:138,140,142,144,146,148)]
evo[evo=="."]=NA
evo_mean<-apply(evo,2,function(x){mean(c(mean(as.numeric(x[which(trainset1$class=="B")]),na.rm=TRUE),
                                         mean(as.numeric(x[which(trainset1$class=="P")]),na.rm=TRUE)))})
evonames<-colnames(evo)
evo<-mean_imputation(evo,evo_mean,evonames)
###gene intolerance
gi<-trainset1[,c(167:204,209:215)]
gi[gi=="."]=NA
gi_mean<-apply(gi,2,function(x){mean(c(mean(as.numeric(x[which(trainset1$class=="B")]),na.rm=TRUE),
                                       mean(as.numeric(x[which(trainset1$class=="P")]),na.rm=TRUE)))})
ginames<-colnames(gi)
gi<-mean_imputation(gi,gi_mean,ginames)
###other features
base<-trainset1[,c(1:5,7,10)]

####merge
class<-ifelse(trainset1$group=="clinvar.B","B","P")
trainsetN1<-cbind(base,pop,popP,pred1,evo,gi,class)
###trainset2 
pop2<-trainset2[,c(11:44,154,198,199)]
pop2[pop2=="."]=0
pred2<-trainset2[,c(45,46,48,49,51,52,54,55,57,58,60,61,63,64,66,67,69,70,72,73,83:84,86:89,114:115,117:120,122,123,125:130,200)]
pred2[pred2=="."]=NA
popP2<-trainset2[,c(50,53,56,59,62,65,68,71,85)]
popP2[popP2=="D"]=2
popP2[popP2=="A"]=2
popP2[popP2=="T"]=1
popP2[popP2=="P"]=1
popP2[popP2=="B"]=1
popP2[popP2=="M"]=2
popP2[popP2=="N"]=1
popP2[popP2=="H"]=2
popP2[popP2=="L"]=1
popP2[popP2=="U"]=1.5
popP2[popP2=="."]=NA
popP2_mean<-apply(popP2,2,function(x){mean(c(mean(as.numeric(x[which(trainset2$class=="B")]),na.rm=TRUE),
                                             mean(as.numeric(x[which(trainset2$class=="P")]),na.rm=TRUE)))})
popP2names<-colnames(popP2)
popP2<-mean_imputation(popP2,popP2_mean,popP2names)


pred2_mean<-apply(pred2,2,function(x){mean(c(mean(as.numeric(x[which(trainset2$class=="B")]),na.rm=TRUE),
                                             mean(as.numeric(x[which(trainset2$class=="P")]),na.rm=TRUE)))})
pred2_names<-colnames(pred2)
pred2<-mean_imputation(pred2,pred2_mean,pred2_names)
pred21<-pred2[,c(c(10,12:18,20)*2,41)]
evol2<-trainset2[,c(131,132,136:138,140,142,144,146,148)]
evol2[evol2=="."]=NA
evol2_mean<-apply(evol2,2,function(x){mean(c(mean(as.numeric(x[which(trainset2$class=="B")]),na.rm=TRUE),
                                             mean(as.numeric(x[which(trainset2$class=="P")]),na.rm=TRUE)))})
evol2names<-colnames(evol2)
evol2<-mean_imputation(evol2,evol2_mean,evol2names)
gi2<-trainset2[,c(160:197,201:207)]
gi2[gi2=="."]=NA
gi2[gi2=="N/A"]=NA
gi2_mean<-apply(gi2,2,function(x){mean(c(mean(as.numeric(x[which(trainset2$class=="B")]),na.rm=TRUE),
                                         mean(as.numeric(x[which(trainset2$class=="P")]),na.rm=TRUE)))})
gi2names<-colnames(gi2)
gi2<-mean_imputation(gi2,gi2_mean,gi2names)
base2<-trainset2[,c(1:5,7,10)]
class2<-trainset2$class
trainsetN2<-cbind(base2,pop2,popP2,pred21,evol2,gi2,class2)
colnames(trainsetN2)[length(trainsetN2)]<-"class"
####trainsetN
trainsetN<-rbind(trainsetN1,trainsetN2)

####normalization####
trainsetN[,8:118]<-sapply(8:118,function(x){numeric_column(trainsetN,x)})
mean_feature<-apply(trainsetN[,45:118],2,mean)
std_feature<-apply(trainsetN[,45:118], 2, sd)
normalization<-function(data_table,scale_mean,scale_std){
  data_table[,45:118]<-z_normalization(data_table[,45:118],scale_mean,scale_std)
  return(data_table)
}
TrainN<-normalization(trainsetN,mean_feature,std_feature)
TrainN1<-TrainN
#####remove features
rmc<-c("syn","non","splice_site","frameshift","syn_z","lof_z","exac_pLI","GERP++_RS_rankscore","integrated_fitCons_score","gnomad211_exome_AF_afr" ,"gnomad211_exome_AF_sas" ,
       "gnomad211_exome_AF_amr" ,"gnomad211_exome_AF_eas","gnomad211_exome_AF_nfe","gnomad211_exome_AF_fin" ,"gnomad211_exome_AF_asj" , "gnomad211_exome_AF_oth")
rmc1<-match(rmc,colnames(TrainN1))
TrainN1<-TrainN1[,-c(1:7,9:12,21:41,rmc1)]

#########################
#Use Boruta to select important features
TrainN1$class<-as.factor(TrainN1$class)
boruta_output <- Boruta(class ~ ., data=na.omit(TrainN1), doTrace=2, maxRuns = 100)  
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")]) 
DataFeaturesSelected <- TrainN1[,c(boruta_signif, "class")]
####PLOT boruta###

#Select the best mtry parameters
bestmtry <- tuneRF(DataFeaturesSelected[,1:(dim(DataFeaturesSelected)[2]-1)], 
                   DataFeaturesSelected[,c("class")], stepFactor=2.0, improve=1e-5, ntree=500)
bestmtry_val <- bestmtry[order(bestmtry[,2]),1][1]
max.mtry <-20
if(bestmtry_val > max.mtry){
  bestmtry_val <- max.mtry
}
control <- trainControl(method="cv", 
                        summaryFunction=twoClassSummary, 
                        classProbs=T,
                        savePredictions=T)
metric <- "Accuracy"
tunegrid <- expand.grid(.mtry=c(bestmtry_val))
#Model training
IRDfit.rf <- train(class~., data=DataFeaturesSelected, weights=rep(1.0,dim(DataFeaturesSelected)[1]),
                method="cforest",
                metric=metric, trControl=control, tuneGrid=tunegrid)

#####SAVE
save(fit.rf,trainset1,trainset2,allPrest,test1Brest,test1P,test1Bs,DataFeaturesSelected,
     TrainN,trainsetN,TrainN1,
     tools,tools2,file="IRDrf.69scalefeatures.rdata")
##########
Test######
##########
###test set1 ####
####Population frequency ###
test.part1<-allPrest
pop<-test.part1[,c(11:44,157,206,207,216)]
pop[pop=="."]=0

###prediction ####
pred<-test.part1[,c(45,46,48,49,51,52,54,55,57,58,60,61,63,64,66,67,69,70,72,73,83:84,86:89,114:115,117:120,122,123,125:130,208)]
popP<-test.part1[,c(50,53,56,59,62,65,68,71,85)]
popP[popP=="D"]=2
popP[popP=="A"]=2
popP[popP=="T"]=1
popP[popP=="P"]=1
popP[popP=="B"]=1
popP[popP=="M"]=2
popP[popP=="N"]=1
popP[popP=="H"]=2
popP[popP=="L"]=1
popP[popP=="U"]=1.5
popP[popP=="."]=NA
pred[pred=="."]=NA
pred_mean<-apply(pred,2,function(x){mean(c(mean(as.numeric(x[which(test.part1$class=="P")]),na.rm=TRUE)))})
prednames<-colnames(pred)
popP_mean<-apply(popP,2,function(x){mean(c(mean(as.numeric(x[which(test.part1$class=="P")]),na.rm=TRUE)))})
pred<-mean_imputation(pred,pred_mean,prednames)
popPnames<-colnames(popP)
popP<-mean_imputation(popP,popP_mean,popPnames)
pred1<-pred[,c(c(10,12:18,20)*2,41)]
###evolutionary####
evo<-test.part1[,c(131,132,136:138,140,142,144,146,148)]
evo[evo=="."]=NA
evo_mean<-apply(evo,2,function(x){mean(c(mean(as.numeric(x[which(test.part1$class=="P")]),na.rm=TRUE)))})
evonames<-colnames(evo)
evo<-mean_imputation(evo,evo_mean,evonames)
###gene intolerance
gi<-test.part1[,c(167:204,209:215)]
gi[gi=="."]=NA
gi_mean<-apply(gi,2,function(x){mean(c(mean(as.numeric(x[which(test.part1$class=="P")]),na.rm=TRUE)))})
ginames<-colnames(gi)
gi<-mean_imputation(gi,gi_mean,ginames)
###other features
base<-test.part1[,c(1:5,7,10)]

####merge
class<-ifelse(test.part1$group=="clinvar.B","B","P")
test.part1<-cbind(base,pop,popP,pred1,evo,gi,class)

###trainset2 
test.part2<-test1Brest
pop2<-test.part2[,c(11:44,154,198,199,208)]
pop2[pop2=="."]=0
pred2<-test.part2[,c(45,46,48,49,51,52,54,55,57,58,60,61,63,64,66,67,69,70,72,73,83:84,86:89,114:115,117:120,122,123,125:130,200)]
pred2[pred2=="."]=NA
popP2<-test.part2[,c(50,53,56,59,62,65,68,71,85)]
popP2[popP2=="D"]=2
popP2[popP2=="A"]=2
popP2[popP2=="T"]=1
popP2[popP2=="P"]=1
popP2[popP2=="B"]=1
popP2[popP2=="M"]=2
popP2[popP2=="N"]=1
popP2[popP2=="H"]=2
popP2[popP2=="L"]=1
popP2[popP2=="U"]=1.5
popP2[popP2=="."]=NA
popP2_mean<-apply(popP2,2,function(x){mean(c(mean(as.numeric(x[which(test.part2$class=="B")]),na.rm=TRUE)))})
popP2names<-colnames(popP2)
popP2<-mean_imputation(popP2,popP2_mean,popP2names)
pred2_mean<-apply(pred2,2,function(x){mean(c(mean(as.numeric(x[which(test.part2$class=="B")]),na.rm=TRUE)))})
pred2_names<-colnames(pred2)
pred2<-mean_imputation(pred2,pred2_mean,pred2_names)
pred21<-pred2[,c(c(10,12:18,20)*2,41)]
evol2<-test.part2[,c(131,132,136:138,140,142,144,146,148)]
evol2[evol2=="."]=NA
evol2_mean<-apply(evol2,2,function(x){mean(c(mean(as.numeric(x[which(test.part2$class=="B")]),na.rm=TRUE)))})
evol2names<-colnames(evol2)
evol2<-mean_imputation(evol2,evol2_mean,evol2names)
gi2<-test.part2[,c(160:197,201:207)]
gi2[gi2=="."]=NA
gi2[gi2=="N/A"]=NA
gi2_mean<-apply(gi2,2,function(x){mean(c(mean(as.numeric(x[which(test.part2$class=="B")]),na.rm=TRUE)))})
gi2names<-colnames(gi2)
gi2<-mean_imputation(gi2,gi2_mean,gi2names)
base2<-test.part2[,c(1:5,7,10)]
class2<-test.part2$class
test.part2<-cbind(base2,pop2,popP2,pred21,evol2,gi2,class2)
colnames(test.part2)[length(test.part2)]<-"class"
####testset1 
testset1<-rbind(test.part1,test.part2)
testset1<-testset1[,-grep("gnomAD_AAF",names(testset1))]

####normalization####
testset1[,8:118]<-sapply(8:118,function(x){numeric_column(testset1,x)})
mean_feature<-apply(testset1[,45:118],2,mean)
std_feature<-apply(testset1[,45:118], 2, sd)
normalization<-function(data_table,scale_mean,scale_std){
  data_table[,45:118]<-z_normalization(data_table[,45:118],scale_mean,scale_std)
  return(data_table)
}
testset1<-normalization(testset1,mean_feature,std_feature)
testset1<-testset1[colnames(DataFeaturesSelected)]
###prediction
test1_out<-predict(IRDfit.rf,newdata=testset1, type = "prob")
test1_out1<-test1_out[,2]
roc1 <- roc(testset1$class, test1_out1)

###test set2 ###
test2<-read.table("path/to/uniprot.xls",sep="\t",
                  header=T,stringsAsFactors = F,check.names = F)
test2$class<-ifelse(test2$class=="B","B","P")
pop2<-test2[,c(11:44,154,198,199)]
pop2[pop2=="."]=0
pred2<-test2[,c(45,46,48,49,51,52,54,55,57,58,60,61,63,64,66,67,69,70,72,73,83:84,86:89,114:115,117:120,122,123,125:130,200)]
pred2[pred2=="."]=NA
popP2<-test2[,c(50,53,56,59,62,65,68,71,85)]
popP2[popP2=="D"]=2
popP2[popP2=="A"]=2
popP2[popP2=="T"]=1
popP2[popP2=="P"]=1
popP2[popP2=="B"]=1
popP2[popP2=="M"]=2
popP2[popP2=="N"]=1
popP2[popP2=="H"]=2
popP2[popP2=="L"]=1
popP2[popP2=="U"]=1.5
popP2[popP2=="."]=NA
popP2_mean<-apply(popP2,2,function(x){mean(c(mean(as.numeric(x[which(test2$class=="B")]),na.rm=TRUE),
                                             mean(as.numeric(x[which(test2$class=="P")]),na.rm=TRUE)))})
popP2names<-colnames(popP2)
popP2<-mean_imputation(popP2,popP2_mean,popP2names)


pred2_mean<-apply(pred2,2,function(x){mean(c(mean(as.numeric(x[which(test2$class=="B")]),na.rm=TRUE),
                                             mean(as.numeric(x[which(test2$class=="P")]),na.rm=TRUE)))})
pred2_names<-colnames(pred2)
pred2<-mean_imputation(pred2,pred2_mean,pred2_names)
pred21<-pred2[,c(c(10,12:18,20)*2,41)]
evol2<-test2[,c(131,132,136:138,140,142,144,146,148)]
evol2[evol2=="."]=NA
evol2_mean<-apply(evol2,2,function(x){mean(c(mean(as.numeric(x[which(test2$class=="B")]),na.rm=TRUE),
                                             mean(as.numeric(x[which(test2$class=="P")]),na.rm=TRUE)))})
evol2names<-colnames(evol2)
evol2<-mean_imputation(evol2,evol2_mean,evol2names)
gi2<-test2[,c(160:197,201:207)]
gi2[gi2=="."]=NA
gi2[gi2=="N/A"]=NA
gi2_mean<-apply(gi2,2,function(x){mean(c(mean(as.numeric(x[which(test2$class=="B")]),na.rm=TRUE),
                                         mean(as.numeric(x[which(test2$class=="P")]),na.rm=TRUE)))})
gi2names<-colnames(gi2)
gi2<-mean_imputation(gi2,gi2_mean,gi2names)
base2<-test2[,c(1:5,7,10)]
class2<-test2$class
test2<-cbind(base2,pop2,popP2,pred21,evol2,gi2,class2)
colnames(test2)[length(test2)]<-"class"
test2[,8:118]<-sapply(8:118,function(x){numeric_column(test2,x)})
mean_feature<-apply(test2[,45:118],2,mean)
std_feature<-apply(test2[,45:118], 2, sd)
test2<-normalization(test2,mean_feature,std_feature)

#####remove features
testset2<-test2[,colnames(DataFeaturesSelected)]
test2_out<-predict(IRDfit.rf,newdata=testset2,type="prob")
test2_out1<-test2_out[,2]
roc2<-roc(testset2$class,test2_out1)
