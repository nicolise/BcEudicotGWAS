#Meta-analysis of all Bc x eudicot experiments 
#040516
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/MetaAnalysis")
#-------------------------------------------------------------
#first get all data frames to match up
SlModDat <- read.csv("SlModelData.csv")
BrModDat <- read.csv("BrModelData.csv")
CeModDat <- read.csv("CeModelData.csv")
CiModDat <- read.csv("CiModelData.csv")
GmModDat <- read.csv("GmModelData.csv")
HcModDat <- read.csv("HcModelData.csv")

#data to include: Exp, Rep, Flat, Domest, PlantGeno, IndPlant, IsolateID, Scale.LS
library(plyr)
#format: new = old
SlModDat <- dplyr::select(SlModDat, Exp = ExpBlock, Rep = PExpRep.x, Flat = AgFlat, Domest = Species, PlantGeno = PlGenoNm, IndPlant = IndPlant, IsolateID = Igeno, Scale.LS = Scale.LS)
SlModDat$Taxon <- "Solanum"
BrModDat <- dplyr::select(BrModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS)
BrModDat$Taxon <- "Brapa"
CeModDat <- dplyr::select(CeModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS)
CeModDat$Taxon <- "Cendivia"
CiModDat <- dplyr::select(CiModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS)
CiModDat$Taxon <- "Cintybus"
GmModDat <- dplyr::select(GmModDat, Exp = exp, Rep = rep, Flat = Image, Domest = origin, PlantGeno = genotypename, IndPlant = plantnumber, IsolateID = isolatename, Scale.LS = Scale.LS)
GmModDat$Taxon <- "Glycine"
GmModDat$IndPlant <- as.factor(GmModDat$IndPlant)
HcModDat <- dplyr::select(HcModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS)
HcModDat$Taxon <- "Helianthus"

write.csv(SlModDat, "SlMetaDat.csv")
write.csv(BrModDat, "BrMetaDat.csv")
write.csv(CeModDat, "CeMetaDat.csv")
write.csv(CiModDat, "CiMetaDat.csv")
write.csv(GmModDat, "GmMetaDat.csv")
write.csv(HcModDat, "HcMetaDat.csv")

FullMetaDat <- rbind(SlModDat, BrModDat, CeModDat, CiModDat, HcModDat, GmModDat)
write.csv(FullMetaDat, "FULLMetaDat.csv")
