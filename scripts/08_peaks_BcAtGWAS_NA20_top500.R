#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS/")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2); library(grid); library(plyr)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("04_bigRRoutput/At_LesionSizes_MAF20_NA20.HEM.PlotFormat.csv")
names(HEM.plotdata)
HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[,-c(1)]

#get threshhold values 
HEM.thresh <- read.csv("04_bigRRoutput/At_LesionSize_MAF20_NA20.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]

#take the top 200 over the threshold for each phenotype
#MAKE variables to hold threshold numbers
TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}

TH999pos <- HEM.thresh[4,]
for (i in 2:ncol(TH999pos)){
  assign(paste("TH999pos_", names(TH999pos[i]), sep=""),as.numeric(TH999pos[i]))
}

TH95pos <- HEM.thresh[1,]
for (i in 2:ncol(TH95pos)){
  assign(paste("TH95pos_", names(TH95pos[i]), sep=""),as.numeric(TH95pos[i]))
}

TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}

TH999neg <- HEM.thresh[8,]
for (i in 2:ncol(TH999neg)){
  assign(paste("TH999neg_", names(TH999neg[i]), sep=""),as.numeric(TH999neg[i]))
}

TH95neg <- HEM.thresh[5,]
for (i in 2:ncol(TH95neg)){
  assign(paste("TH95neg_", names(TH95neg[i]), sep=""),as.numeric(TH95neg[i]))
}

#get all SNPs > 99% Thr
for (i in 4:7){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,i)))
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,i)))
}

#for top 500 only
for (i in 4:7){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), head(arrange(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMpos.", names(HEM.plotdata[i]), sep=""))[,4])), n=500))
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), tail(arrange(get(paste("HEMneg.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))[,4])), n=500))
#combine pos and neg by group
  assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
  mydf <- paste("HEM.", names(HEM.plotdata[i]), sep="")
  renamedf <- get(mydf)
  colnames(renamedf)[4] <- "Effect"
  assign(mydf, renamedf)
  myblob <- names(HEM.plotdata[i])
  assign(mydf, cbind(get(mydf), Trait = myblob))
}

HEM.topSNPs <- rbind(HEM.Col0.Les, HEM.coi1.Les, HEM.npr1.Les, HEM.pad3.Les)

Top500SNP <- HEM.topSNPs

library(ggplot2)
plot1 <- ggplot(HEM.topSNPs, aes(x=Pos, y=Effect))
plot1 + geom_point(aes(color=factor(Trait)))+ theme_bw()

Top500SNP$Chrom <- gsub("Chromosome", "", Top500SNP$Chrom)
Top500SNP$Chrom <- as.numeric(as.character(Top500SNP$Chrom))
Top500SNP$Pos <- as.numeric(as.character(Top500SNP$Pos))
write.csv(HEM.topSNPs, "05_processing/Top1000SNPs_BcAtGWAS_NA20.csv")
