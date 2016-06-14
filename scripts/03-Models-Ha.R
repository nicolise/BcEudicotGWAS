#modeling for Helianthus annuus
#041316
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/")
#-------------------------------------------------------------
#load data
ModDat <- read.csv("data/MetaAnalysis/HaMetaDat.csv")

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
HaSumm <- as.data.frame(with(ModDat, table(IsolateID,Exp)))
#missing Exps: 1.05.04 (1), KGB1 (2)
OgDat <- ModDat
#do separately due to error
ModDat <- subset(ModDat, IsolateID != ("1.05.04"))
ModDat <- subset(ModDat, IsolateID != ("KGB1"))

#run the model
library(lme4); library(car); library(lmerTest)

#this one fails
fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat) + (1|Domest/PlantGeno/IndPlant), data = ModDat)
#this too
fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Domest/PlantGeno/IndPlant), data = ModDat)

Sys.time()
#model is nearly unidentifiable
fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat), data = ModDat)
Sys.time()
sink(file='output/ModelOutputs/HcFullMod_060716a.txt')
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

#fails when including 1|Exp/Rep

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file="output/ModelOutputs/HaLSMeans061416.txt")
for (i in c(1:12)) {
  print(unique(out[[i]]$PlantGeno))
  #this one works
  Lesion.lm <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|IndPlant) + (1|Exp:IsolateID), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "IsolateID")
  print(Lesion.lsm)
}
sink()