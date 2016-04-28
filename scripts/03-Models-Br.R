#modeling for Brassica rapa
#040816
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/MetaAnalysis")
setwd("~/Documents/GitRepos/BcEudicotGWAS/data")
#-------------------------------------------------------------
#load data
ModDat <- read.csv("BrMetaDat.csv")

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
BrSumm <- as.data.frame(with(ModDat, table(IsolateID,Rep)))
#missing Exps: 01.05.22 (1), 02.04.21 (1), Geranium (2), Navel (2)
OgDat <- ModDat
ModDat <- subset(ModDat, IsolateID != c("01.05.22","02.04.21","Geranium","Navel"))

#run the model
library(lme4); library(car); library(lmerTest)
Sys.time()
sink(file='BrFullMod_042716.txt')
print("fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat) + (1|IndPlant), data = ModDat)")

fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/Flat) + (1|IndPlant), data = ModDat) 
#this one does not give p values
#dropping 1|IndPlant and 1|Exp/Rep/Flat and 1|Exp/Rep doesn't fix it either. 

#working model
#do get p-vals, BrFullMod_041116.txt
##fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp), data = ModDat) 

##sink(file='BrFullMod_041116b.txt')
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
#adding XX: error
#fails at  PI.347594 (plant 3) when including exp/rep/flat
#with 1|exp/rep, PI.508409 fails (plant 5)
sink(file="LSMeans041116.txt")
for (i in c(1:12)) {
  print(unique(out[[i]]$PlantGeno))
  Lesion.lm <- lmer(Scale.LS ~ IsolateID + (1|Exp), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "IsolateID")
  print(Lesion.lsm)
}
sink()
