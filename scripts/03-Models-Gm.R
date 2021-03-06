#modeling for Glycine max
#041216
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/")
#-------------------------------------------------------------
#load data
ModDat <- read.csv("data/MetaAnalysis/GmMetaDat.csv")
names(ModDat)

#-------------------------------------------------------------
#check assumptions
#check normality of Scale.LS
require(car)
ModDat$Scale.LS.t <- ModDat$Scale.LS + 1
#not quite normal
qqp(ModDat$Scale.LS.t, "norm")

#---------------------------------------------------------------
#run the model
library(lme4); library(car); library(lmerTest)

#fails - model nearly unidentifiable
#fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat) + (1|Domest/PlantGeno/IndPlant), data = ModDat)
#fails when dropping 1|Exp/Rep/Flat
#plus when dropping 1|Exp/Rep too
#so must remove indplant

Sys.time()
fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) , data = ModDat)
Sys.time()

sink(file='GmFullMod_041216.txt')
print("fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) , data = ModDat)")
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

#fails when (1|Exp/Rep)

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file="output/ModelOutputs/GmLSMeans_060716.txt")
for (i in c(1:12)) {
  print(unique(out[[i]]$PlantGeno))
  #this one works
  Lesion.lm <- lmer(Scale.LS ~ IsolateID + (1|Exp), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "IsolateID")
  print(Lesion.lsm)
}
sink()