#modeling for Cichorium intybus
#041216
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/MetaAnalysis")
#-------------------------------------------------------------
#load data
ModDat <- read.csv("CiMetaDat.csv")

names(ModDat)

#-------------------------------------------------------------
#check assumptions

#check normality of Scale.LS
require(car)
ModDat$Scale.LS.t <- ModDat$Scale.LS + 1
#not quite normal
qqp(ModDat$Scale.LS.t, "norm")

#---------------------------------------------------------------
#try removing isolates missing from one experiment
CiSumm <- as.data.frame(with(ModDat, table(IsolateID,Exp)))
#missing Exps: 2.04.21 (1), 1.05.22 (1), Geranium (2), Navel (2)
OgDat <- ModDat
ModDat <- subset(ModDat, IsolateID != c("01.05.22","02.04.21","Geranium","Navel"))

#run the model
library(lme4); library(car); library(lmerTest)

#rerunning with isolates w/ missing experiments removed. 06/05/16
#fails to converge
#fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat) + (1|Domest/PlantGeno/IndPlant), data = ModDat)
#unable to evaluate scaled gradient/ fails to converge when you remove just 1|Exp/Rep/Flat

#working model for isolates minus those with missing exps
Sys.time()
#removed indplant
fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat), data = ModDat)
Sys.time()
sink(file='CiFullMod_060516b.txt')
print("fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat), data = ModDat)")
Sys.time()
rand(fullmod)
Anova(fullmod, type=2)
anova(fullmod)
Sys.time()
sink()


#lsmeans
#run model per isolate WITHIN each plant genotype
#so include no species terms or plant genotype terms
attach(ModDat)
out <- split( ModDat , f = ModDat$PlantGeno)
head(out[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
#fails when including exp/rep
sink(file="output/ModelOutputs/CiLSMeans_060616.txt")
print("Lesion.lm <- lmer(Scale.LS ~ IsolateID + (1|Exp), data=out[[i]])")
for (i in c(1:12)) {
  print(unique(out[[i]]$PlantGeno))
  #this one works
  Lesion.lm <- lmer(Scale.LS ~ IsolateID + (1|Exp), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "IsolateID")
  print(Lesion.lsm)
}
sink()