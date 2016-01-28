#Preliminary data analysis Eudicot spp. x Botrytis cinerea lesions
#011916
#-------------------------------------------------------------
#clear memory
rm(list=ls())
#set working directory
setwd("~/path/to/dir")
#read-in your data files
#MyDat contains all lesion measurements along with plant, isolate, experiment ID etc.
MyDat <- read.csv("file01.csv")
#IsoNm and PlantNm contain both plant and isolate IDs and alternate names
IsoNm <- read.csv("file02.csv")
PlantNm <- read.csv("file03.csv")
#AddUnits contains one value per image (pixels per cm)
AddUnits <- read.csv("file04.csv")

names(MyDat)
#subset MyDat columns to include only phenotype of interest (e.g. lesion size)
LsDat <- MyDat[,c(1:11,152)]
#Add pixels per cm conversion
names(LsDat)
names(AddUnits)
#add a column "Sort" that includes both experiment ID and image ID 
#add column to both LsDat and AddUnits
LsDat$Sort <- paste(LsDat$PExpRep, LsDat$PImage, sep='') 
AddUnits$Sort <- paste(AddUnits$PExpRep, AddUnits$Image, sep='')
#add unit conversion column to LsDat
LsDat2 <- merge(LsDat, AddUnits, by="Sort")
#subset relevant columns
LsDat <- LsDat2[,c(2:13,19)]
#Add column of lesion size in cm squared
LsDat <- transform(LsDat, Scale.LS=(Lesion.Size/(pixelsPcm^2)))

library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr)
library(dplyr)

#add original isolate names (B05.10 vs. "i03" or something)
unique(unlist(LsDat$Pexp))
IsoNm96 <- IsoNm
colnames(IsoNm96)[1] <- "Piso"
#add isolate names
LsDat96v2 <- merge(LsDat96, IsoNm96, by="Piso")
LsDat96v2 <-LsDat96v2[,c(1:14,18)] #make sure this includes the column Isolate
SrtDat <- LsDat96v2

#add a column for plant*iso interaction
SrtDat$PbyI <- paste(SrtDat$PPlant, SrtDat$Isolate, sep='') 

#remove any duplicate entries
SrtDat <- unique(SrtDat)

#----------------------------------------------------------------------
#linear model
#nesting: B within A as A/B or A + A:B
#fixed effects: PInLflt
#random effects: PPlant, Isolate, PInPlant, PInLeaf, Pexp

#add a domestication term
names(SrtDat)
unique(SrtDat$PPlant)
#FL NC TX MA KS OR are Domest
SrtDat$Domest <- ifelse(SrtDat$PPlant %in% c("FL","NC","TX","MA","KS","OR"),"Dm", "Wl")

#remove Control Isolates
unique(SrtDat$Isolate)
SrtDat <- SrtDat[SrtDat$Isolate!="Control",]

#prepare data frame for statistics
ModDat <- SrtDat
names(ModDat)
#keep and rename certain columns
ModDat <- dplyr::select(ModDat, ExpBlock = Pexp, Igeno = Isolate, Pgeno = PPlant, AorB = PInLflt, Leaf = PInLeaf, Plant = PInPlant, AgFlat = PImage, Species = Domest, matches("."))
#remove any rows with NA for plant or isolate
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
ModDat <- completeFun(ModDat, c("Pgeno", "Igeno"))

#add a column of correct plant geno names
head(ModDat)
PlGenNum <- dplyr::select(PlantNm, Pgeno = PlantID, matches("."))
PlGenNum$PlGenoNm <- paste("LA", PlGenNum$PlantGeno, sep='') 
PlGenNum <- PlGenNum[c(1:12),c("Pgeno","PlGenoNm")]
ModDat <- merge(ModDat, PlGenNum, by="Pgeno")
names(ModDat)
#Plant is coded as numeric and nested within Pgeno
ModDat$IndPlant <- paste(ModDat$PlGenoNm, ModDat$Plant, sep='.') 

write.csv(ModDat, "ModelData.csv")