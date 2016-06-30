rm(list=ls())
setwd("~/Projects")

tomdat <- read.csv("BcSolGWAS/data/GWAS_files/03_bigRRinput/Sl_Pheno_bigRR.csv")
tomdat <- tomdat[,2:14]

tomdat2 <- reshape(tomdat, varying = c(2:13),idvar="Igeno", direction="long", sep="")
tomdat2 <- tomdat2
library("plyr")
tomdat2 <- rename(tomdat2, c("time"="Plant", "LA"="LSMeans"))
tomdat2$Plant <- as.factor(tomdat2$Plant)
tomdat2$Plant <- revalue(tomdat2$Plant, c("4355" = "LA4355", "4345" = "LA4345", "3475" = "LA3475", "3008" = "LA3008", "2706" = "LA2706", "2176" = "LA2176", "2093" = "LA2093", "1684" = "LA1684", "1589" = "LA1589", "1547" = "LA1547", "480" = "LA0480", "410" = "LA0410"))
tomdat3 <- ddply(tomdat2, "Plant", summarize,
      mean = mean(LSMeans), 
      CV = sd(LSMeans)/mean(LSMeans))

#soy data (had trouble with reshape)
soydat <- read.csv("BcEudicotGWAS/data/Thresholding/RG_soyphenos.csv")
soymeans <- as.data.frame(apply(soydat[,2:13], 2, mean))
soymeansSD <- as.data.frame(apply(soydat[,2:13], 2, sd))
soymeans$SD <- soymeansSD
soymeans <- rename(soymeans, c("apply(soydat[, 2:13], 2, mean)"="mean"))
soymeans$CV <- soymeans$SD / soymeans$mean
#write.csv(soymeans, "BcEudicotGWAS/data/Thresholding/SoyCV.csv")

#now At data
Atdat<- read.csv("BcEudicotGWAS/data/Thresholding/At_Pheno_bigRR.csv")
Atdat <- Atdat[,3:17]
Atdat <- rename(Atdat, c("Rename" = "Isolate"))
Atmeans <- as.data.frame(apply(Atdat[,2:15], 2, mean))
AtSD <- as.data.frame(apply(Atdat[,2:13], 2, sd, na.rm=T))
Atmeans$SD <- AtSD

#get CV and mean for lesion size for each tomato genotype
soyCV <- read.csv("BcEudicotGWAS/data/Thresholding/SoyCV.csv")
AtCV <- read.csv("BcEudicotGWAS/data/Thresholding/AtCV.csv")
tomCV <- tomdat3
write.csv(tomCV, "BcEudicotGWAS/data/Thresholding/TomCV.csv")

#fill those in to the threshold data
thresh <- read.csv("BcEudicotGWAS/data/Thresholding/AllThresholds.csv")
allCV <- read.csv("BcEudicotGWAS/data/Thresholding/AllCV.csv")
mydata <- merge(thresh, allCV, by="Plant")

mydata_Th99 <- subset(mydata, Threshold=="0.99 Thresh")
mydata_Th95 <- subset(mydata, Threshold=="0.95 Thresh")
mydata_Th999 <- subset(mydata, Threshold=="0.999 Thresh")

library(ggplot2)
plot1 <- ggplot(mydata_Th95, aes(x=Fx, y=NormMean))
plot1 + geom_point(aes(color=factor(Exp))) + 
  scale_x_continuous(limits=c(0,5e-05))

plot2 <- ggplot(mydata_Th99, aes(x=Fx, y=NormMean))
plot2 + geom_point(aes(color=factor(Exp)))

plot3 <- ggplot(mydata_Th999, aes(x=Fx, y=NormMean))
plot3 + geom_point(aes(color=factor(Exp)))+ 
  theme_bw() +
  scale_x_continuous(limits=c(0,0.00075)) +
  scale_y_continuous(limits=c(0.6, 1.4))

plot4 <- ggplot(mydata_Th95, aes(x=Fx, y=CV))
plot4 + geom_point(aes(color=factor(Exp))) + 
  theme_bw() +
  scale_x_continuous(limits=c(0,4e-05)) +
  scale_y_continuous(limits=c(0.3, 0.6))

plot5 <- ggplot(mydata_Th99, aes(x=Fx, y=CV))
plot5 + geom_point(aes(color=factor(Exp)))

plot6 <- ggplot(mydata_Th999, aes(x=Fx, y=CV))
plot6 + geom_point(aes(color=factor(Exp))) + 
  theme_bw() +
  scale_x_continuous(limits=c(0,.00075)) +
  scale_y_continuous(limits=c(0.3, 0.55))

