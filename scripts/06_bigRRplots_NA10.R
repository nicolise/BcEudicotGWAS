#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2); library(grid); library(plyr)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/BcAtGWAS/04_bigRRoutput/At_LesionSizes_MAF20_NA10.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[,-c(1)]

#get threshhold values 
HEM.thresh <- read.csv("data/BcAtGWAS/04_bigRRoutput/At_LesionSize_MAF20_NA10.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]

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

#Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Segment <- as.numeric(as.character(HEM.plotdata$Segment))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Segment, Pos)), ]

#now make segments line up consecutively
HEM.plotdata$Chrom.Seg <- paste(HEM.plotdata$Chrom, HEM.plotdata$Segment, sep=".")
HEM.plotdata$Chrom.Seg <- as.numeric(HEM.plotdata$Chrom.Seg)

#let's try making the chrom.seg integers so that R isn't confused
unique(HEM.plotdata$Chrom.Seg)
HEM.plotdata$Chrom.Seg.F <- as.factor(HEM.plotdata$Chrom.Seg)
unique(HEM.plotdata$Chrom.Seg.F)
recode.vars <- data.frame(OGvals = c(1,1.1,2,2.1,2.2,2.3,2.4,2.5, 3, 3.1, 3.2, 4,  4.1, 5,  5.1, 6,  6.1, 6.2, 6.3, 7,  7.1, 7.2, 8,  8.1, 8.2, 9,  9.1, 10, 10.1, 11, 11.1, 12, 12.1, 13, 13.1, 13.2, 14, 14.1, 14.2, 15, 15.1, 15.2, 15.3, 15.4, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
print(recode.vars)

HEM.plotdata$Chrom.Seg.Int <- recode.vars$newvals[match(HEM.plotdata$Chrom.Seg.F, recode.vars$OGvals)]
unique(HEM.plotdata$Chrom.Seg.Int)

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing

for (i in unique(HEM.plotdata$Chrom.Seg.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom.Seg.Int==i-1)$Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Pos+lastbase
  }
  #set ticks to be a list of existing ticks, plus the current Index
  #floor rounds it down to the nearest whole number
  # ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
  
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index)/2)+1])
}
#ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

#create a custom color scale
myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#make plots for each phenotype
names(HEM.plotdata)
#without the loop
for (i in c(4:7)){
  #jpeg(paste("paper/plots/ActualPaper/bw_Sl_LesionSize_trueMAF20_NA10_lowTR_", names(HEM.plotdata[7]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
  jpeg(paste("plots/At_LesionSize_trueMAF20_NA10_lowTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
  plot(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
         theme_bw()+
         colScale+
         geom_point(aes(color = factor(Chrom)))+
         labs(list(y="SNP Effect Estimate", title=paste("Phenotype ", names(HEM.plotdata[i]))))+
         guides(col = guide_legend(nrow = 8, title="Chromosome"))+
         geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
         geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
         geom_text(aes(0,get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", vjust = 1.2, hjust=.05), col = "black")+
         geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_text(aes(0,get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", vjust = 1.5, hjust=.05), col = "black")+
         theme(legend.position="none")+
         #NA20 chromosomes
         scale_x_continuous(name="Chromosome", breaks = c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682,  38915229), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
         #NA10 chromosomes
         # scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0))
  dev.off()
}

#highTR BW
#plot 2706 as example, [9]
#with labeling removed for paper
#jpeg(paste("paper/plots/ActualPaper/FigR5/Routs/pm_BW_Sl_MAF20_highTR_", names(HEM.plotdata[8]), ".ManhattanPlot.jpg", sep=""), width=7.5, height=4, units='in', res=600)
#4 to 15
for (i in c(4:7)){
  jpeg(paste("plots/At_LesionSize_trueMAF20_NA10_hiTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
  plot(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
         theme_bw()+
         colScale+
         geom_point(aes(color = factor(Chrom)))+
         labs(list(y="SNP Effect Estimate", title=paste("Phenotype ", names(HEM.plotdata[i]))))+
         guides(col = guide_legend(nrow = 8, title="Chromosome"))+
         geom_hline(yintercept=get(paste("TH999pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_hline(yintercept=get(paste("TH999neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_text(aes(0,get(paste("TH999neg_", names(HEM.plotdata[i]), sep="")), label = "99.9% Threshold", vjust = 1.2, hjust = .05), col = "black")+
         theme(legend.position="none")+
         scale_x_continuous(name="Chromosome", breaks = c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682,  38915229), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0))
  dev.off()
}