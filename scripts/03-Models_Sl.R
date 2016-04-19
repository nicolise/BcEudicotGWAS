#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/MetaAnalysis")

#load data
ModDat <- read.csv("SlMetaDat.csv")

#-------------------------------------------------------------

#from other eudicots
#check assumptions

#check normality of Scale.LS
require(car)
ModDat$Scale.LS.t <- ModDat$Scale.LS + 1
#not quite normal
qqp(ModDat$Scale.LS.t, "norm")

#---------------------------------------------------------------
#run the model
library(lme4); library(car); library(lmerTest)

Sys.time()
#this model worked in my original BcSolGWAS analysis
#fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep/Flat), + (1|Exp/Rep) + (1|Species/PlGenoNm/IndPlant) + (1|Species/PlGenoNm/IndPlant/AorB), data = ModDat)
sink(file='SlFullMod_041816.txt')
print("fullmod <- lmer(Scale.LS ~ IsolateID + Domest/PlantGeno + IsolateID:Domest/PlantGeno + IsolateID:Domest + (1|Exp) + (1|Exp/Rep/Flat), + (1|Exp/Rep) + (1|Species/PlGenoNm/IndPlant) + (1|Species/PlGenoNm/IndPlant/AorB), data = ModDat)")
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
#original one from BcSolGWAS:
#Lesion.lm <- lmer(Scale.LS ~ Igeno + Plant/Leaf/AorB + (1|Exp), data=out[[i]])
sink(file="SlLSMeans041816.txt")
for (i in c(1:12)) {
  print(unique(out[[i]]$PlantGeno))
  Lesion.lm <- lmer(Scale.LS ~ IsolateID + IndPlant/Leaf/AorB + (1|Exp), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "IsolateID")
  print(Lesion.lsm)
}
sink()