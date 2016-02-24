#Preliminary data analysis Eudicot spp. x Botrytis cinerea lesions
#011916
#-------------------------------------------------------------
#clear memory
rm(list=ls())
#set working directory
#read-in your data files
#MyDat contains all lesion measurements along with plant, isolate, experiment ID etc.
#IsoNm and PlantNm contain both plant and isolate IDs and alternate names
#AddUnits contains one value per image (pixels per cm)

#for GJS Chicorium intybus
setwd("~/Projects/BcEudicotGWAS/data/GJS_files/CiAnalysis")
MyDat <- read.csv("FullCiData.csv")
IsoNm <- read.csv("IsolateIDs.csv")
PlantNm <- read.csv("PlantIDs.csv")
AddUnits <- read.csv("Scaling.csv")


names(MyDat)
#subset MyDat columns to include only phenotype of interest (e.g. lesion size)
CiDat <- MyDat[,c(1:6,147)]
#Add pixels per cm conversion
names(CiDat)
names(AddUnits)
#add unit conversion column to CiDat
CiDat2 <- merge(CiDat, AddUnits, by="Image")
#subset relevant columns
names(CiDat2)
CiDat <- CiDat2[,c(2:9)]
#Add column of lesion size in cm squared
CiDat <- transform(CiDat, Scale.LS=(Lesion.Size/(PixelsPcm^2)))

library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr)
library(dplyr)

#add original isolate names (B05.10 vs. "i03" or something)
names(IsoNm)
names(CiDat)
#add isolate names
CiDatv2 <- merge(CiDat, IsoNm, by="IsolateName")
names(CiDatv2)
CiDatv2 <-CiDatv2[,c(2:10)] #make sure this includes the column Isolate
SrtDat <- CiDatv2

#add a column for plant*iso interaction
names(SrtDat)
SrtDat$PbyI <- paste(SrtDat$PlantGeno, SrtDat$IsolateID, sep='_') 

#remove any duplicate entries
SrtDat <- unique(SrtDat)

#add a column of correct plant geno names
head(SrtDat)
head(PlantNm)
ModDat <- merge(SrtDat, PlantNm, by="PlantGeno") #losing some at this stage
names(ModDat)
 
#----------------------------------------------------------------------
#prepare data frame for statistics
names(ModDat)
#keep and rename certain columns
ModDat <- dplyr::select(ModDat, ExpBlock = Pexp, Igeno = Isolate, Pgeno = PPlant, AorB = PInLflt, Leaf = PInLeaf, Plant = PInPlant, AgFlat = PImage, Species = Domest, matches("."))
#remove any rows with NA for plant or isolate
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
SrtDat <- completeFun(SrtDat, c("PlantGeno", "IsolateID"))

write.csv(SrtDat, "CiModelData.csv")
