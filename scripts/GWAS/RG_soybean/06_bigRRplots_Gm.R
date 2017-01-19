#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcEudicotGWAS/data/RG_files")
setwd("~/Projects/BcEudicotGWAS/data/RG_files")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("04_bigRRoutput/Gm_Les.Pheno_HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("04_bigRRoutput/Gm_Les.Pheno_HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]
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

# #Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Segment <- as.numeric(as.character(HEM.plotdata$Segment))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata2 <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Segment, Pos)), ]

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#want to figure out where to add +500 to draw breaks between chromosomes
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=+lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

colnames(HEM.plotdata)[16] <- "myIndex"

#make plots for each phenotype
#it isn't working with the loop for some reason
#for (m in 4:15){ #[15]
  jpeg(paste("05_plots/Gm_LesionSize_MAF20_lowTR_", names(HEM.plotdata[11]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
  ggplot(HEM.plotdata, aes(x=myIndex, y=abs(HEM.plotdata[11])))+
    theme_bw()+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[11]))))+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))+
    geom_hline(yintercept=get(paste("TH95_", names(HEM.plotdata[11]), sep="")), colour = "blue") +
    geom_text(aes(0,get(paste("TH95_", names(HEM.plotdata[11]), sep="")), label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")+
#    geom_hline(yintercept=get(paste("TH99_", names(HEM.plotdata[11]), sep="")), colour = "black") +
#    geom_text(aes(0,get(paste("TH99_", names(HEM.plotdata[11]), sep="")), label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black")+
   expand_limits(y=0)
dev.off()
#}
#[15]
#just need: column 5, 6, 12, 14, 15
#[6]
jpeg(paste("05_plots/Gm_LesionSize_MAF20_highTR_", names(HEM.plotdata[6]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=myIndex, y=abs(HEM.plotdata[6])))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[6]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH999_", names(HEM.plotdata[6]), sep=""))) +
  geom_text(aes(0,get(paste("TH999_", names(HEM.plotdata[6]), sep="")), label =
  ".999 Threshold", vjust = 1.5, hjust = .05), col = "black")+
  expand_limits(y=-0.001)
dev.off()
#}
#stop [6]

#or just do them each by hand
jpeg("plots/Sp_LesionSize_greyIntx.jpeg", width=6, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=Index, y=abs(LA4345)))+
  labs(list(y="SNP Effect Estimate", title = "Lesion Size on LA4345", x="Chromosome position"))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=TH95_LA4345, colour = "blue") +
  geom_text(aes(0,TH95_LA4345, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")+
  geom_hline(yintercept=TH99_LA4345) +
  geom_text(aes(0,TH99_LA4345, label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black")
dev.off()

#geom_hline(yintercept=TH999_LA3475) +
#  geom_text(aes(0,TH999_LA3475, label = ".999 Threshold", vjust = 1.5, hjust = .05), col = "black")

#how many SNPs are above a certain threshhold?
topSNP <- sum(HEM.plotdata$LA0410 >= 0.001) #41
highSNP <- sum(HEM.plotdata$LA0410 >= TH999_LA0410) #spits out number
totalSNP <- sum(HEM.plotdata$LA0410 >= 0)
highSNP/totalSNP*100