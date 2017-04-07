#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS/AtBcGWAS_no_ogs/")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("At_LesionSizenoogs_MAF20.HEM.PlotFormat.csv")
names(HEM.plotdata)
HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("At_LesionSizenoogs_MAF20.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]

#take the top 200 over the threshold for each phenotype
library(plyr)

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

names(HEM.plotdata)

#Col0.Cam: take all over 0.99
#take a pos set and a neg set
HEM.Col0.Cam.up <- HEM.plotdata[HEM.plotdata$Col0.Cam.HEM > get(paste("TH99pos_", "Col0.Cam.HEM", sep="")),] 
HEM.Col0.Cam.dn <- HEM.plotdata[HEM.plotdata$Col0.Cam.HEM < get(paste("TH99neg_", "Col0.Cam.HEM", sep="")),] 
HEM.Col0.Cam <- rbind(HEM.Col0.Cam.up, HEM.Col0.Cam.dn)
HEM.Col0.Cam  <- HEM.Col0.Cam[,c("Chrom", "Segment", "Pos", "Col0.Cam.HEM")]
names(HEM.Col0.Cam)
HEM.Col0.Cam <- rename(HEM.Col0.Cam, c("Col0.Cam.HEM" = "Effect"))
HEM.Col0.Cam$Pheno <- "Col0.Cam"
HEM.Col0.Cam <- head(arrange(HEM.Col0.Cam,desc(abs(Effect))), n = 200)

#but for now just the top 200 SNPs for each
names(HEM.plotdata)

#Col0.AT4G30530
count <- HEM.plotdata[HEM.plotdata$Col.0AT4G30530.HEM > get(paste("TH99pos_", "Col.0AT4G30530.HEM", sep="")),] 
HEM.Col0.AT4G30530 <- subset(HEM.plotdata, select=c(Chrom, Segment, Pos, Col.0AT4G30530.HEM))
HEM.Col0.AT4G30530 <- rename(HEM.Col0.AT4G30530, c("Col.0AT4G30530.HEM" = "Effect"))
HEM.Col0.AT4G30530$Pheno <- "Col0.AT4G30530"
HEM.Col0.AT4G30530 <- head(arrange(HEM.Col0.AT4G30530,desc(abs(Effect))), n = 200)

#Col0.AT3G26830
count <- HEM.plotdata[HEM.plotdata$Col0.AT3G26830.HEM > get(paste("TH99pos_", "Col0.AT3G26830.HEM", sep="")),] 
HEM.Col0.AT3G26830 <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,Col0.AT3G26830.HEM))
HEM.Col0.AT3G26830 <- rename(HEM.Col0.AT3G26830, c("Col0.AT3G26830.HEM" = "Effect"))
HEM.Col0.AT3G26830$Pheno <- "Col0.AT3G26830"
HEM.Col0.AT3G26830 <- head(arrange(HEM.Col0.AT3G26830,desc(abs(Effect))), n = 200)

#anac055.Cam
count <- HEM.plotdata[HEM.plotdata$anac055.Cam.HEM > get(paste("TH99pos_", "anac055.Cam.HEM", sep="")),] 
HEM.anac055.Cam <- subset(HEM.plotdata, select=c(Chrom, Segment, Pos, anac055.Cam.HEM))
HEM.anac055.Cam <- rename(HEM.anac055.Cam, c("anac055.Cam.HEM" = "Effect"))
HEM.anac055.Cam$Pheno <- "anac055.Cam"
HEM.anac055.Cam <- head(arrange(HEM.anac055.Cam,desc(abs(Effect))), n = 200)

#coi1.Cam
#95% threshold for this one
count <- HEM.plotdata[HEM.plotdata$coi1.Cam.HEM > get(paste("TH95pos_", "coi1.Cam.HEM", sep="")),] 
HEM.coi1.Cam <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,coi1.Cam.HEM))
HEM.coi1.Cam <- rename(HEM.coi1.Cam, c("coi1.Cam.HEM" = "Effect"))
HEM.coi1.Cam$Pheno <- "coi1.Cam"
HEM.coi1.Cam <- head(arrange(HEM.coi1.Cam,desc(abs(Effect))), n = 200)

#Col0.Cam
count <- HEM.plotdata[HEM.plotdata$Col0.Cam.HEM > get(paste("TH99pos_", "Col0.Cam.HEM", sep="")),] 
HEM.Col0.Cam <- subset(HEM.plotdata, select=c(Chrom, Segment, Pos, Col0.Cam.HEM))
HEM.Col0.Cam <- rename(HEM.Col0.Cam, c("Col0.Cam.HEM" = "Effect"))
HEM.Col0.Cam$Pheno <- "Col0.Cam"
HEM.Col0.Cam <- head(arrange(HEM.Col0.Cam,desc(abs(Effect))), n = 200)

#npr1.Cam
count <- HEM.plotdata[HEM.plotdata$npr1.Cam.HEM > get(paste("TH99pos_", "npr1.Cam.HEM", sep="")),] 
HEM.npr1.Cam <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,npr1.Cam.HEM))
HEM.npr1.Cam <- rename(HEM.npr1.Cam, c("npr1.Cam.HEM" = "Effect"))
HEM.npr1.Cam$Pheno <- "npr1.Cam"
HEM.npr1.Cam <- head(arrange(HEM.npr1.Cam,desc(abs(Effect))), n = 200)

#tga3.Cam
count <- HEM.plotdata[HEM.plotdata$tga3.Cam.HEM > get(paste("TH99pos_", "tga3.Cam.HEM", sep="")),] 
HEM.tga3.Cam <- subset(HEM.plotdata, select=c(Chrom, Segment, Pos, tga3.Cam.HEM))
HEM.tga3.Cam <- rename(HEM.tga3.Cam, c("tga3.Cam.HEM" = "Effect"))
HEM.tga3.Cam$Pheno <- "tga3.Cam"
HEM.tga3.Cam <- head(arrange(HEM.tga3.Cam,desc(abs(Effect))), n = 200)

#anac055.Les
count <- HEM.plotdata[HEM.plotdata$anac055.Les.HEM > get(paste("TH99pos_", "anac055.Les.HEM", sep="")),] 
HEM.anac055.Les <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos, anac055.Les.HEM))
HEM.anac055.Les <- rename(HEM.anac055.Les, c("anac055.Les.HEM" = "Effect"))
HEM.anac055.Les$Pheno <- "anac055.Les.s"
HEM.anac055.Les <- head(arrange(HEM.anac055.Les,desc(abs(Effect))), n = 200)

#Col0.Les
count <- HEM.plotdata[HEM.plotdata$Col0.Les.HEM > get(paste("TH99pos_", "Col0.Les.HEM", sep="")),] 
HEM.Col0.Les <- subset(HEM.plotdata, select=c(Chrom, Segment, Pos, Col0.Les.HEM))
HEM.Col0.Les <- rename(HEM.Col0.Les, c("Col0.Les.HEM" = "Effect"))
HEM.Col0.Les$Pheno <- "Col0.Les.s"
HEM.Col0.Les <- head(arrange(HEM.Col0.Les,desc(abs(Effect))), n = 200)

#npr1.Les
count <- HEM.plotdata[HEM.plotdata$npr1.Les > get(paste("TH99pos_", "npr1.Les", sep="")),] 
HEM.npr1.Les <- subset(HEM.plotdata, select=c(Chrom, Segment, Pos, npr1.Les))
HEM.npr1.Les <- rename(HEM.npr1.Les, c("npr1.Les" = "Effect"))
HEM.npr1.Les$Pheno <- "npr1.Les"
HEM.npr1.Les <- head(arrange(HEM.npr1.Les,desc(abs(Effect))), n = 200)

#tga3.Les
count <- HEM.plotdata[HEM.plotdata$tga3.Les > get(paste("TH99pos_", "tga3.Les", sep="")),] 
HEM.tga3.Les <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,tga3.Les))
HEM.tga3.Les <- rename(HEM.tga3.Les, c("tga3.Les" = "Effect"))
HEM.tga3.Les$Pheno <- "tga3.Les.s"
HEM.tga3.Les <- head(arrange(HEM.tga3.Les,desc(abs(Effect))), n = 200)

#pad3.Les
count <- HEM.plotdata[HEM.plotdata$pad3.Les > get(paste("TH99pos_", "pad3.Les", sep="")),] 
HEM.pad3.Les <- subset(HEM.plotdata, select=c(Chrom, Segment, Pos, pad3.Les))
HEM.pad3.Les <- rename(HEM.pad3.Les, c("pad3.Les" = "Effect"))
HEM.pad3.Les$Pheno <- "pad3.Les"
HEM.pad3.Les <- head(arrange(HEM.pad3.Les,desc(abs(Effect))), n = 200)

Top200SNP <- rbind(HEM.Col0.AT4G30530, HEM.Col0.AT3G26830, HEM.anac055.Cam, HEM.coi1.Cam, HEM.Col0.Cam, HEM.npr1.Cam, HEM.tga3.Cam, HEM.anac055.Les, HEM.Col0.Les, HEM.npr1.Les, HEM.tga3.Les, HEM.pad3.Les)

Top200SNP$Chrom <- gsub("Chromosome", "", Top200SNP$Chrom)
Top200SNP$Chrom <- as.numeric(as.character(Top200SNP$Chrom))
Top200SNP$Pos <- as.numeric(as.character(Top200SNP$Pos))
write.csv(Top200SNP, "Top200SNPs_BcAtGWAS_noogs.csv")
