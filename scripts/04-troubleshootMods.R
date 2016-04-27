#troubleshoot models for BcEudicot GWAS
#Nicole E Soltis
#04/27/16
#---------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/MetaAnalysis")

SlModDat <- read.csv("SlMetaDat.csv")
BrModDat <- read.csv("BrMetaDat.csv")
CeModDat <- read.csv("CeMetaDat.csv")
CiModDat <- read.csv("CiMetaDat.csv")
GmModDat <- read.csv("GmMetaDat.csv")
HcModDat <- read.csv("HcMetaDat.csv")

names(SlModDat)
attach(SlModDat)
interaction.plot(IsolateID, Rep, Lesion.Size)
#eliminate 96a208b?
SlSumm <- as.data.frame(with(SlModDat, table(IsolateID,Rep)))
#missing Reps: KGB2 missing from 96b208b
#missing Exps: 94.1 (2), Gallo3 (1)

BrSumm <- as.data.frame(with(BrModDat, table(IsolateID,Rep)))
#missing Reps: none else
#missing Exps: 01.05.22 (1), 02.04.21 (1), Geranium (2), Navel (2)

CeSumm <- as.data.frame(with(CeModDat, table(IsolateID,Rep)))
#missing Reps: none else
#missing Exps: 94.1 (1), UKRazz (1), 01.04.15 (2), Gallo1 (2)

CiSumm <- as.data.frame(with(CiModDat, table(IsolateID,Rep)))
#missing Reps: none else
#missing Exps: 01.04.15 (1), Gallo1 (1), 01.02.18 (2), FresaSD (2)

GmSumm <- as.data.frame(with(GmModDat, table(IsolateID,Rep)))
#missing Reps: none!
#missing Exps: none!

HcSumm <- as.data.frame(with(HcModDat, table(IsolateID,Rep)))
#missing Reps:
#missing Exps: 01.05.04 (1), KGB1 (2)

#no isolate:plant terms missing entirely!

SlPlant <- as.data.frame(with(SlModDat, table(IsolateID, PlantGeno)))
BrPlant <- as.data.frame(with(BrModDat, table(IsolateID, PlantGeno)))
CePlant <- as.data.frame(with(CeModDat, table(IsolateID, PlantGeno)))
CiPlant <- as.data.frame(with(CiModDat, table(IsolateID, PlantGeno)))
GmPlant <- as.data.frame(with(GmModDat, table(IsolateID, PlantGeno)))
HcPlant <- as.data.frame(with(HcModDat, table(IsolateID, PlantGeno)))
