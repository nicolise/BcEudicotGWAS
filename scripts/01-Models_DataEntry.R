#statistics for all Bc x eudicot experiments 
#040816
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
SlModDat <- dplyr::select(SlModDat, Exp = ExpBlock, Rep = PExpRep.x, Flat = AgFlat, Domest = Species, PlantGeno = PlGenoNm, IndPlant = IndPlant, IsolateID = Igeno, Scale.LS = Scale.LS, matches("."))
SlModDat$Taxon <- "Solanum"
BrModDat <- dplyr::select(BrModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS, matches("."))
BrModDat$Taxon <- "Brapa"
CeModDat <- dplyr::select(CeModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS, matches("."))
CeModDat$Taxon <- "Cendivia"
CiModDat <- dplyr::select(CiModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS, matches("."))
CiModDat$Taxon <- "Cintybus"
GmModDat <- dplyr::select(GmModDat, Exp = exp, Rep = rep, Flat = Image, Domest = origin, PlantGeno = genotypename, IndPlant = plantnumber, IsolateID = isolatename, Scale.LS = Scale.LS, matches("."))
GmModDat$Taxon <- "Glycine"
GmModDat$IndPlant <- as.factor(GmModDat$IndPlant)
HcModDat <- dplyr::select(HcModDat, Exp = Exp, Rep = Rep, Flat = ImageName, Domest = Domestication, PlantGeno = PlantID, IndPlant = PnumID, IsolateID = IsolateID, Scale.LS = Scale.LS, matches("."))
HcModDat$Taxon <- "Helianthus"

#need to match up isolate names
# IsoNameMatch <- data.frame(SlNames=factor(1:98),stringsAsFactors=T)
# SlNames <- (unique(SlModDat$IsolateID))
# length(SlNames) <- nrow(IsoNameMatch)
# IsoNameMatch$SlNames <- SlNames
# IsoNameMatch$BrNames <- unique(BrModDat$IsolateID)
# IsoNameMatch$CeNames <- unique(CeModDat$IsolateID)
# IsoNameMatch$CiNames <- unique(CiModDat$IsolateID)
# GmNames <- (unique(GmModDat$IsolateID))
# length(GmNames) <- nrow(IsoNameMatch)
# IsoNameMatch$GmNames <- GmNames
# HcNames <- (unique(HcModDat$IsolateID))
# length(HcNames) <- nrow(IsoNameMatch)
# IsoNameMatch$HcNames <- HcNames
# write.csv(IsoNameMatch, "IsolateNameMatch.csv")

#rename isolates
#format: old = new
BrModDat$IsolateID <- revalue(BrModDat$IsolateID, c("FresaS.D." = "FresaSD", "Katietomato" = "KatieTomato", "Mex-03" = "Mex03", "NobelRot" = "NobleRot", "Peppersub" = "PepperSub", "Tripple3(T3)" = "Triple3", "Tripple7(T7)" = "Triple7", "UKrazz" = "UKRazz")) 
CeModDat$IsolateID <- revalue(CeModDat$IsolateID, c("FresaS.D." = "FresaSD", "Katietomato" = "KatieTomato", "Mex-03" = "Mex03", "NobelRot" = "NobleRot", "Peppersub" = "PepperSub", "Tripple3(T3)" = "Triple3", "Tripple7(T7)" = "Triple7", "UKrazz" = "UKRazz"))
CiModDat$IsolateID <- revalue(CiModDat$IsolateID, c("FresaS.D." = "FresaSD", "Katietomato" = "KatieTomato", "Mex-03" = "Mex03", "NobelRot" = "NobleRot", "Peppersub" = "PepperSub", "Tripple3(T3)" = "Triple3", "Tripple7(T7)" = "Triple7", "UKrazz" = "UKRazz"))
HcModDat$IsolateID <- revalue(HcModDat$IsolateID, c("FresaS.D." = "FresaSD", "Katietomato" = "KatieTomato", "Mex-03" = "Mex03", "NobelRot" = "NobleRot", "Peppersub" = "PepperSub", "Tripple3(T3)" = "Triple3", "Tripple7(T7)" = "Triple7", "UKrazz" = "UKRazz"))
GmModDat$IsolateID <- revalue(GmModDat$IsolateID, c("1.01.6" = "1.01.06", "FresaS.D." = "FresaSD", "Peppersub" = "PepperSub", "Triple3(T3)" = "Triple3", "Triple7(T7)" = "Triple7", "UKrazz" = "UKRazz"))


write.csv(SlModDat, "SlMetaDat.csv")
write.csv(BrModDat, "BrMetaDat.csv")
write.csv(CeModDat, "CeMetaDat.csv")
write.csv(CiModDat, "CiMetaDat.csv")
write.csv(GmModDat, "GmMetaDat.csv")
write.csv(HcModDat, "HcMetaDat.csv")