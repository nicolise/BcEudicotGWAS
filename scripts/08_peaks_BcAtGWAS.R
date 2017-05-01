#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS/05_processing")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("At_LesionSize_MAF20.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("At_LesionSize_MAF20.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]

#take the top 200 over the threshold for each phenotype
library(plyr)

TH99 <- HEM.thresh[3,]
for (i in 2:ncol(TH99)){
  assign(paste("TH99_", names(TH99[i]), sep=""),as.numeric(TH99[i]))
}

TH999 <- HEM.thresh[4,]
for (i in 2:ncol(TH999)){
  assign(paste("TH999_", names(TH999[i]), sep=""),as.numeric(TH999[i]))
}

TH95 <- HEM.thresh[1,]
for (i in 2:ncol(TH95)){
  assign(paste("TH95_", names(TH95[i]), sep=""),as.numeric(TH95[i]))
}

names(HEM.plotdata)

#Col0.Cam: take all over 0.99
HEM.Col0.Cam <- subset(HEM.plotdata, abs(Col0.Cam) > get(paste("TH99_", "Col0.Cam", sep="")), 
                                select=c(Chrom,Segment, Pos,Col0.Cam))
HEM.Col0.Cam <- rename(HEM.Col0.Cam, c("Col0.Cam" = "Effect"))
HEM.Col0.Cam$Pheno <- "Col0.Cam"
HEM.Col0.Cam <- head(arrange(HEM.Col0.Cam,desc(Effect)), n = 200)

#but for now just the top 200 SNPs for each

#Col0.Cam
HEM.Col0.Cam <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,Col0.Cam))
HEM.Col0.Cam <- rename(HEM.Col0.Cam, c("Col0.Cam" = "Effect"))
HEM.Col0.Cam$Pheno <- "Col0.Cam"
HEM.Col0.Cam <- head(arrange(HEM.Col0.Cam,desc(abs(Effect))), n = 200)

#Col0.Les.s
HEM.Col0.Les.s <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,Col0.Les.s))
HEM.Col0.Les.s <- rename(HEM.Col0.Les.s, c("Col0.Les.s" = "Effect"))
HEM.Col0.Les.s$Pheno <- "Col0.Les.s"
HEM.Col0.Les.s <- head(arrange(HEM.Col0.Les.s,desc(abs(Effect))), n = 200)

#Col0.AT3G26830
HEM.Col0.AT3G26830 <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,Col0.AT3G26830))
HEM.Col0.AT3G26830 <- rename(HEM.Col0.AT3G26830, c("Col0.AT3G26830" = "Effect"))
HEM.Col0.AT3G26830$Pheno <- "Col0.AT3G26830"
HEM.Col0.AT3G26830 <- head(arrange(HEM.Col0.AT3G26830,desc(abs(Effect))), n = 200)

#Col0.AT4G30530
HEM.Col0.AT4G30530 <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,Col.0AT4G30530))
HEM.Col0.AT4G30530 <- rename(HEM.Col0.AT4G30530, c("Col.0AT4G30530" = "Effect"))
HEM.Col0.AT4G30530$Pheno <- "Col0.AT4G30530"
HEM.Col0.AT4G30530 <- head(arrange(HEM.Col0.AT4G30530,desc(abs(Effect))), n = 200)

#anac055.Cam
HEM.anac055.Cam <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,anac055.Cam))
HEM.anac055.Cam <- rename(HEM.anac055.Cam, c("anac055.Cam" = "Effect"))
HEM.anac055.Cam$Pheno <- "anac055.Cam"
HEM.anac055.Cam <- head(arrange(HEM.anac055.Cam,desc(abs(Effect))), n = 200)

#anac055.Les.s
HEM.anac055.Les.s <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,anac055.Les.s))
HEM.anac055.Les.s <- rename(HEM.anac055.Les.s, c("anac055.Les.s" = "Effect"))
HEM.anac055.Les.s$Pheno <- "anac055.Les.s"
HEM.anac055.Les.s <- head(arrange(HEM.anac055.Les.s,desc(abs(Effect))), n = 200)

#coi1.Cam
HEM.coi1.Cam <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,coi1.Cam))
HEM.coi1.Cam <- rename(HEM.coi1.Cam, c("coi1.Cam" = "Effect"))
HEM.coi1.Cam$Pheno <- "coi1.Cam"
HEM.coi1.Cam <- head(arrange(HEM.coi1.Cam,desc(abs(Effect))), n = 200)

#coi1.Les.s
HEM.coi1.Les.s <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,coi1.Les.s))
HEM.coi1.Les.s <- rename(HEM.coi1.Les.s, c("coi1.Les.s" = "Effect"))
HEM.coi1.Les.s$Pheno <- "coi1.Les.s"
HEM.coi1.Les.s <- head(arrange(HEM.coi1.Les.s,desc(abs(Effect))), n = 200)

#npr1.Cam
HEM.npr1.Cam <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,npr1.Cam))
HEM.npr1.Cam <- rename(HEM.npr1.Cam, c("npr1.Cam" = "Effect"))
HEM.npr1.Cam$Pheno <- "npr1.Cam"
HEM.npr1.Cam <- head(arrange(HEM.npr1.Cam,desc(abs(Effect))), n = 200)

#npr1.Les.s
HEM.npr1.Les.s <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,npr1.Les.s))
HEM.npr1.Les.s <- rename(HEM.npr1.Les.s, c("npr1.Les.s" = "Effect"))
HEM.npr1.Les.s$Pheno <- "npr1.Les.s"
HEM.npr1.Les.s <- head(arrange(HEM.npr1.Les.s,desc(abs(Effect))), n = 200)

#tga3.Cam
HEM.tga3.Cam <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,tga3.Cam))
HEM.tga3.Cam <- rename(HEM.tga3.Cam, c("tga3.Cam" = "Effect"))
HEM.tga3.Cam$Pheno <- "tga3.Cam"
HEM.tga3.Cam <- head(arrange(HEM.tga3.Cam,desc(abs(Effect))), n = 200)

#tga3.Les.s 
HEM.tga3.Les.s <- subset(HEM.plotdata, select=c(Chrom,Segment, Pos,tga3.Les.s))
HEM.tga3.Les.s <- rename(HEM.tga3.Les.s, c("tga3.Les.s" = "Effect"))
HEM.tga3.Les.s$Pheno <- "tga3.Les.s"
HEM.tga3.Les.s <- head(arrange(HEM.tga3.Les.s,desc(abs(Effect))), n = 200)

Top200SNP <- rbind(HEM.tga3.Cam, HEM.tga3.Les.s, HEM.npr1.Les.s, HEM.npr1.Cam, HEM.coi1.Les.s, HEM.coi1.Cam, HEM.anac055.Les.s, HEM.anac055.Cam, HEM.Col0.AT4G30530, HEM.Col0.AT3G26830, HEM.Col0.Cam, HEM.Col0.Les.s)

Top200SNP$Chrom <- gsub("Chromosome", "", Top200SNP$Chrom)
Top200SNP$Chrom <- as.numeric(as.character(Top200SNP$Chrom))
Top200SNP$Pos <- as.numeric(as.character(Top200SNP$Pos))
write.csv(Top200SNP, "Top200SNPs_BcAtGWAS_trueMAF.csv")