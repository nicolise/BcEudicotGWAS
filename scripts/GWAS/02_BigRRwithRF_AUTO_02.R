#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcEudicotGWAS/")
setwd("~/Projects/BcEudicotGWAS")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/BcAtGWAS/04_bigRRoutput/At_phenos.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("data/BcAtGWAS/04_bigRRoutput/At_phenos.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]
TH99 <- HEM.thresh[3,]
TH99_Col0.Cam <- as.numeric(TH99[2])
TH99_Les.s <- as.numeric(TH99[3])
TH99_Col0.AT3G26830 <- as.numeric(TH99[4])
TH99_Col.0AT4G30530 <- as.numeric(TH99[5])
TH99_anac055.Cam <- as.numeric(TH99[6])
TH99_anac055.Les.s <- as.numeric(TH99[7])
TH99_coi1.Cam <- as.numeric(TH99[8])
TH99_coi1.Les.s <- as.numeric(TH99[9])
TH99_npr1.Cam <- as.numeric(TH99[10])
TH99_npr1.Les.s <- as.numeric(TH99[11])
TH99_tga3.Cam <- as.numeric(TH99[12])
TH99_tga3.Les.s <- as.numeric(TH99[13])

TH999 <- HEM.thresh[4,]
TH999_Col0.Cam <- as.numeric(TH999[2])
TH999_Col0.Les.s <- as.numeric(TH999[3])
TH999_Col0.AT3G26830 <- as.numeric(TH999[4])
TH999_Col.0AT4G30530 <- as.numeric(TH999[5])
TH999_anac055.Cam <- as.numeric(TH999[6])
TH999_anac055.Les.s <- as.numeric(TH999[7])
TH999_coi1.Cam <- as.numeric(TH999[8])
TH999_coi1.Les.s <- as.numeric(TH999[9])
TH999_npr1.Cam <- as.numeric(TH999[10])
TH999_npr1.Les.s <- as.numeric(TH999[11])
TH999_tga3.Cam <- as.numeric(TH999[12])
TH999_tga3.Les.s <- as.numeric(TH999[13])

TH95 <- HEM.thresh[1,]
TH95_Col0.Cam <- as.numeric(TH95[2])
TH95_Col0.Les.s <- as.numeric(TH95[3])
TH95_Col0.AT3G26830 <- as.numeric(TH95[4])
TH95_Col.0AT4G30530 <- as.numeric(TH95[5])
TH95_anac055.Cam <- as.numeric(TH95[6])
TH95_anac055.Les.s <- as.numeric(TH95[7])
TH95_coi1.Cam <- as.numeric(TH95[8])
TH95_coi1.Les.s <- as.numeric(TH95[9])
TH95_npr1.Cam <- as.numeric(TH95[10])
TH95_npr1.Les.s <- as.numeric(TH95[11])
TH95_tga3.Cam <- as.numeric(TH95[12])
TH95_tga3.Les.s <- as.numeric(TH95[13])

TH975 <- HEM.thresh[2,]
TH975_Col0.Cam <- as.numeric(TH975[2])
TH975_Col0.Les.s <- as.numeric(TH975[3])
TH975_Col0.AT3G26830 <- as.numeric(TH975[4])
TH975_Col.0AT4G30530 <- as.numeric(TH975[5])
TH975_anac055.Cam <- as.numeric(TH975[6])
TH975_anac055.Les.s <- as.numeric(TH975[7])
TH975_coi1.Cam <- as.numeric(TH975[8])
TH975_coi1.Les.s <- as.numeric(TH975[9])
TH975_npr1.Cam <- as.numeric(TH975[10])
TH975_npr1.Les.s <- as.numeric(TH975[11])
TH975_tga3.Cam <- as.numeric(TH975[12])
TH975_tga3.Les.s <- as.numeric(TH975[13])

# #Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata2 <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Pos)), ]

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequencial		-- accurate position indexing
##RFF code isn't working... replaced "Pos" with "pos"
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

hist(HEM.plotdata$Index)

#make plots for each phenotype
jpeg("plots/BcAtGWAS/AtCol0_Cam.ManhattanPlot.jpg")
qplot(Index,abs(Col0.Cam), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "Camalexin_AtCol0", colour=factor(Chrom)) +
 geom_hline(yintercept=TH99_Col0.Cam) +
  geom_text(aes(0,TH99_Col0.Cam, label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black") +
geom_hline(yintercept=TH95_Col0.Cam, colour = "blue") +
  geom_text(aes(0,TH95_Col0.Cam, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")
dev.off()

#how many SNPs are above a certain threshhold?
topSNP <- sum(HEM.plotdata$LA0410 >= 0.001) #41
highSNP <- sum(HEM.plotdata$LA0410 >= TH999_LA0410) #spits out number
totalSNP <- sum(HEM.plotdata$LA0410 >= 0)
highSNP/totalSNP*100

#rest of plots
jpeg("plots/Practice.ManhattanPlot.jpg")
qplot(Index,abs(tga3.Les.s), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "tga3.Les.s", colour=factor(Chrom)) +
  geom_hline(yintercept=TH999_tga3.Les.s) +
  geom_text(aes(0,TH999_tga3.Les.s, label = ".999 Threshold", vjust = 1.5, hjust = .05), col = "black") 
#+
#  geom_hline(yintercept=TH95_tga3.Cam, colour = "blue") +
#  geom_text(aes(0,TH95_tga3.Cam, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")
dev.off()

