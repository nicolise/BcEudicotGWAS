#modeling for Cichorium endivia
#041116
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/MetaAnalysis")
#-------------------------------------------------------------
#load data
ModDat <- read.csv("CeMetaDat.csv")

names(ModDat)

#-------------------------------------------------------------
#check assumptions

#check normality of Scale.LS
require(car)
ModDat$Scale.LS.t <- ModDat$Scale.LS + 1
#not quite normal
qqp(ModDat$Scale.LS.t, "norm")

#---------------------------------------------------------------
#remove isolates missing from one exp
CeSumm <- as.data.frame(with(ModDat, table(IsolateID,Rep)))
#missing Exps: 94.1 (1), UKRazz (1), 01.04.15 (2), Gallo1 (2)
OgDat <- ModDat
ModDat <- subset(ModDat, IsolateID != c("94.1", "UKRazz", "01.04.15", "Gallo1"))
#run the model
library(lme4); library(car); library(lmerTest)
# fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat) + (1|IndPlant), data = ModDat)
#fails to converge with 1|Exp/Rep/Flat

#fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Domest/PlantGeno/IndPlant), data = ModDat) 
#fails to converge

#working model
#p-vals in CeFullMod_041116.txt
##fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|IndPlant), data = ModDat)

#but, it does not make any sense to include 1|IndPlant without nesting under Domest/PlantGeno so--- removing it!
#p-vals in CeFullMod_060416.txt
#in dataframe with isolates dropped: failes with 1|exp/rep/flat and 1|domes/plantgeno/indplant
#also fails to converge with 1|exp/rep and 1|domest/plantgeno/indplant
#final working model below (06/04)
Sys.time()
fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat), data = ModDat)

Sys.time()
sink(file='CeFullMod_060416.txt')
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
#fails at  PI.273578 (plant 1) when including exp/rep AND 1|Indplant
#fails at PI.503585 (plant 4) when including exp/rep
sink(file="CeLSMeans041216.txt")
for (i in c(1:12)) {
  print(unique(out[[i]]$PlantGeno))
  #this one works
  Lesion.lm <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|IndPlant), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "IsolateID")
  print(Lesion.lsm)
}
sink()