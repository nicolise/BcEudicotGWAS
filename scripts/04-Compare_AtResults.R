#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/path/to/dir")
AtDat <- read.csv("96BcRAWlesiondata.csv")
ModDat <- read.csv("ModelData.csv")
names(AtDat)

library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr)
library(dplyr)

#subset to include only lesion size
AtDat <- AtDat[,c(1:10,151)]

#add a column of lesion size in cm squared for AtDat
AtDat <- transform(AtDat, Scale.LS=(Lesion.Size/(Scale..Pix.cm.^2)))

#-----------------------------------------------------------------
#Now to combine data sets
names(ModDat)
List1 <- unique(ModDat$Igeno)
List2 <- unique(AtDat$Isolate)
AtDat <- AtDat[AtDat$Isolate!="Ctrl",]
df1 <- data.frame(matrix(unlist(List1), nrow=100, byrow=T),stringsAsFactors=FALSE)
df2 <- data.frame(matrix(unlist(List2), nrow=98, byrow=T),stringsAsFactors=FALSE)

#remove unmatched from ModDat
#AtDatOG <- AtDat
AtDat<- AtDatOG
#ModDat <- ModDat[ModDat$Igeno!=c("01.02.05","94.1"),]
#AtDat <- AtDat[AtDat$Isolate!="01.02.13", "01.02.19", ,]
ModDat$Igeno <- revalue(ModDat$Igeno, c("MEAP6G" = "MEAPGG"))
AtDat$Isolate <- revalue(AtDat$Isolate, c("Apple 404" = "Apple404", "Triple 3 (T3)" = "Triple3", "Davis Navel" = "DavisNavel", "Gallo 2" = "Gallo2", "Kern A2" = "KernA2", "Kern B2" = "KernB2", "Kern B1" = "KernB1", "Apple 517" = "Apple517", "UK razz" = "UKRazz", "Triple 7 (T7)" = "Triple7", "Philo Menlo" = "PhiloMenlo", "Pepper sub" = "PepperSub", "Noble Rot" = "NobleRot", "Mex-03" = "Mex03", "MEA PGG" = "MEAPGG", "Katie tomato" = "KatieTomato", "Gallo 1" = "Gallo1", "Gallo 2" = "Gallo2", "Fresa S.D." = "FresaSD", "Fresa 525" = "Fresa525", "Esparato Fresa" = "EsparatoFresa", "BPA 1 " = "BPA1"))

#keep only: Experiment, Plant, Isolate, Scale.LS
#and ExpBlock, Igeno, Pgeno, Species, Scale.LS
AtDat2 <- AtDat[,c(1,5,7,12)]
ModDat2 <- ModDat[,c(1,2,3,8,16)]
names(AtDat2)
names(ModDat2)
ModDat2 <- dplyr::select(ModDat2, Experiment = ExpBlock, Plant = Pgeno, Isolate = Igeno, matches("."))
AtDat2$Species <- "At"
AtDat2$Experiment <- as.factor(AtDat2$Experiment)

FullDat <- rbind(ModDat2, AtDat2)

#means by isolate*plant
#violin plot for Wild, Domest, At
#add "Plant" to list of Les.means factors if want to look within plant genos
names(FullDat)
Les.means <- ddply(FullDat, c("Isolate", "Species"), summarise, mean=mean(Scale.LS))
Les.means$NumSp <- ifelse(Les.means$Species == "At", 3, ifelse (Les.means$Species == "Dm", 1, 2))

#add Les.means scaled to mean of each
tapply(Les.means$mean, Les.means$Species, mean)
#At = 0.2649255
#Dm = 0.7309409
#Wl = 0.6008727
Les.means$GpMean <- ifelse(Les.means$Species == "At", 0.2649255, ifelse (Les.means$Species == "Dm", 0.7309409, 0.6008727))
Les.means <- transform(Les.means, RelLes=(mean/GpMean))

#add a factor for species x exp
Les.means2 <- ddply(FullDat, c("Isolate", "Species", "Experiment"), summarise, mean=mean(Scale.LS))
Les.means2$NumSp <- ifelse(Les.means2$Species == "At", 3, ifelse (Les.means2$Species == "Dm", 1, 2))
tapply(Les.means2$mean, Les.means2$Species, mean)
#At = 0.2631917
#Dm = 0.7366425
#Wl = 0.620993
Les.means2$GpMean <- ifelse(Les.means2$Species == "At", 0.2631917, ifelse (Les.means2$Species == "Dm", 0.7366425, 0.620993))
Les.means2 <- transform(Les.means2, RelLes=(mean/GpMean))
Les.means2$ExpbySp <- paste(Les.means2$NumSp, Les.means2$Experiment, sep='')


#-------------------------------------------------------------------------
#scatterplot
myitalic <- element_text(face="italic")
names(Les.means)
#add a column of mRL (mean of Relative Lesion Size) per isolate
#sort dataframe by mRL
#then color by the new factor mRL
Les.means$mRL<- ave(Les.means$RelLes, Les.means$Isolate)
attach(Les.means)
Les.means <- Les.means[order(mRL),]
attach(Les.means)
ggplot(Les.means, aes(x = NumSp, y = RelLes))+
  geom_point()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(mRL), group=factor(Isolate)), show_guide=F)+
  scale_x_discrete(breaks=c("1","2","3"),
                   labels=c(expression(paste(italic("S. lycopersicum"))), expression(paste(italic("S. pimpinellifolium"))), expression(paste(italic("A. thaliana")))))+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 15, hjust = 1, vjust=1))+
  labs(y=expression(Relative ~ Lesion ~ Area), x=element_blank())
#-----------------------------------------------
#scatterplot by experiment
attach(Les.means2)
Les.means3 <- Les.means2[Les.means2$Species!="At",]
names(Les.means3)
Les.means3$ExpbyIso <- paste(Les.means3$Experiment, Les.means3$Isolate, sep='')
file <- ddply(Les.means3, c("Experiment", "Isolate"), summarise, Newmean = mean(mean))
Les.means4 <- Les.means3[Les.means3$Isolate!="EsparatoFresa",]

ggplot(Les.means4, aes(x = ExpbySp, y = mean))+
  geom_point()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(Isolate), group=factor(Isolate)), show_guide=F)+
  scale_x_discrete(breaks=c("196a","196b","296a","296b"),
                   labels=c(expression(paste(italic("S. lycopersicum"))~ E1), expression(paste(italic("S. lycopersicum"))~ E2), expression(paste(italic("S. pimpinellifolium"))~E1), expression(paste(italic("S. pimpinellifolium"))~E2)))+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 15, hjust = 1, vjust=1))+
  labs(y=expression(Mean ~ Lesion ~ Area~(cm^2)), x=element_blank())
#--------------------------------------------------------------------------------------
#violin plot

library(ggplot2)
attach(Les.means)
ggplot (data = Les.means, 
  aes(x=factor(NumSp), y=mean))+
  geom_violin(adjust = 0.7, scale = "width", fill="mediumspringgreen")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 15, hjust = 1, vjust=1),
        text = element_text(size=24))+
  scale_x_discrete(breaks=c("1","2","3"),
                   labels=c(expression(paste(italic("S. lycopersicum"))), expression(paste(italic("S. pimpinellifolium"))), expression(paste(italic("A. thaliana")))))+
  geom_boxplot(width=0.1)+
  labs(y=expression(Mean~Lesion~Area~(cm^2)), x=NULL)
