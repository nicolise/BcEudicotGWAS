#Meta-analysis of all Bc x eudicot experiments 
#040516
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/")
#-------------------------------------------------------------
#load data
MyData <- read.csv("data/MetaAnalysis/FULLMetaDat.csv")
#rename so that Domest is Domesticated vs. Wild vs. landraces
MyData$Domest[MyData$Domest == 'D'] <- 'Domesticated'
MyData$Domest[MyData$Domest == 'Dm'] <- 'Domesticated'
MyData$Domest[MyData$Domest == 'W'] <- 'Wild'
MyData$Domest[MyData$Domest == 'Wl'] <- 'Wild'

#-------------------------------------------------------------
#violin plot of lesion size by plant species
attach(MyData)
library(ggplot2)
p <- ggplot(MyData, aes(factor(Taxon), Scale.LS, fill=Domest))

tiff("plots/Eudicots_LesionSize_beanplots.tiff", width=10, height=4, units='in', res=600)
p + geom_violin(trim=T)+
  ylim(0,5) +
  theme_bw()+
  theme(text = element_text(size=14), axis.title.x = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14, angle = 30, hjust = 0.65, vjust = 0.8)) +
  scale_x_discrete(labels=c("Brassica rapa", "Cichorium endivia", "Cichorium intybus", "Glycine max", "Helianthus annuus", "Solanum spp."))+
  ylab("Scaled Lesion Size")+
 geom_boxplot(width=0.1, position = position_dodge(width = 0.9))+
  labs(fill="")
dev.off()


#get mean lesion size per isolate
library(plyr)
means <- ddply(MyData, "IsolateID", summarise, meanLS = mean(Scale.LS))
write.csv(means, "output/MeanLsAllTaxa.csv")


#-------------------------------------------------------------
##scatterplot by plant species
names(MyData)
attach(MyData)
#add PlantNum as an integer sorted by mean lesion size
library(plyr)
#1 value per plant:isolate combination
FigDat <- ddply(MyData, c("PlantGeno", "IsolateID", "Domest", "Taxon"), summarise,
                 meanBYall   = mean(Scale.LS))
#1 value per isolate
FigDat$meanBYiso <- ave(FigDat$meanBYall, FigDat$IsolateID, FUN = mean)
#1 value per isolate, per taxon
FigDat$meanBYisoTax <- ave(FigDat$meanBYiso, FigDat$Taxon, FUN = mean)
#1 value per isolate, per species / domest
FigDat$meanBYisoSp <- ave(FigDat$meanBYisoTax, FigDat$Domest, FUN = mean)
#1 value per taxon to scale lesion size
FigDat$meanBYtax <- ave(FigDat$meanBYall, FigDat$Taxon, FUN = mean)
#add a vector of 1:6 per domest in order of meanBYplant
MDmeans <- ddply(MyData, c("PlantGeno","Domest", "Taxon"), summarise, meanBYplant=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Taxon, MDmeans$Domest, MDmeans$meanBYplant),] 
MDmeans$PlNum <- (rep(1:6, times=12, each=1))
#combine DFs
FigDat <- merge(FigDat, MDmeans[,c(1,4:5)], by="PlantGeno")
library(ggplot2)
attach(FigDat)
names(FigDat)

#scale LS by Taxon
FigDat$scaleLS <- (FigDat$meanBYall / FigDat$meanBYtax)

#draw the figure
attach(FigDat)
ggplot(FigDat, aes(x = PlNum, y = scaleLS))+
  geom_point()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(meanBYiso), group=factor(IsolateID)), show.legend=F)+
  facet_grid(~Taxon*Domest, scales="fixed", space="free_x")+
  theme(text = element_text(size=18))+ #, axis.text.x = element_text(angle = 45, hjust = 1)
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())

#draw a figure highlighting saprophytes
#find 5 lowest lesion sizes
sort(FigDat$meanBYiso, TRUE)
hist(FigDat$meanBYiso, breaks=100)
MyMins <- c(0.1658313, 0.2960654, 0.3719376, 0.4307290, 0.4507375)
#Davis Navel, 
sort(FigDat$meanBYiso, F)
MyMaxs <- c(1.3252362, 1.2791752, 1.2661299, 1.2450447, 1.2349104)

#and for scaled
sort(FigDat$scaleLS, TRUE)

#list isolates by mean 
IsoVals <- FigDat[,c("IsolateID", "meanBYiso")]
IsoVals <- unique(IsoVals)
#LowLes is only lowest 5 isolates, HiLes is highest 5
FigDat$LowLes <- FigDat$scaleLS
FigDat$LowLes[FigDat$meanBYiso > 0.45074] <- NA
FigDat$LowLes[FigDat$meanBYiso < 0.45072] <- NA
FigDat$HiLes <- FigDat$scaleLS
FigDat$HiLes[FigDat$meanBYiso < 1.245] <- NA
FigDat$HiLes[FigDat$meanBYiso > 1.246] <- NA

#dimensions: 20 by 2
#draw grey figure
ggplot(FigDat, aes(x = PlNum, y = scaleLS))+
  theme_bw()+
  geom_line(size=1, color = "grey54", aes(group=factor(IsolateID)), show.legend=F, alpha=0.4)+
  geom_line(size=1, aes(x =  PlNum, y = HiLes, group=factor(IsolateID)), color = "black")+
  facet_grid(~Taxon*Domest, scales="fixed", space="free_x")+
  theme(text = element_text(size=18))+ #, axis.text.x = element_text(angle = 45, hjust = 1)
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())

#----------------------------------------------------------------

