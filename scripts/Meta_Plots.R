#Meta-analysis of all Bc x eudicot experiments 
#040516
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcEudicotGWAS/data/MetaAnalysis")
#-------------------------------------------------------------
#load data
MyData <- read.csv("FULLMetaDat.csv")
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
p + geom_violin(trim=T)+
  ylim(0,5) +
  theme( axis.title.x = element_blank()) +
  ylab("Scaled Lesion Size")+
 geom_boxplot(width=0.1, position = position_dodge(width = 0.9))

#-------------------------------------------------------------
##scatterplot by plant species
names(MyData)
attach(MyData)

MyData$avgLS <- ave(MyData$Scale.LS, MyData$IsolateID)
#add PlantNum as an integer sorted by mean lesion size
library(plyr)
FigDat <- ddply(MyData, c("PlantGeno", "IsolateID", "Domest", "Taxon", "avgLS"), summarise,
                 mLS   = mean(Scale.LS))
#"mean" is mean lesion size by plant genotype.
MDmeans <- ddply(MyData, c("PlantGeno","Domest", "Taxon"), summarise, mean=mean(Scale.LS))
#add TxMean to scale lesion size by taxon
MDmeans$TxMean <- ave(MDmeans$mean, MDmeans$Taxon)
MDmeans <- MDmeans[order(MDmeans$Taxon, MDmeans$Domest, MDmeans$mean),] 
MDmeans$PlNum <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)
FigDat <- merge(FigDat, MDmeans[,c(1,4:6)], by="PlantGeno")
FigDat <- dplyr::select(FigDat, PlOrder = mean, matches("."))
FigDat$PlOrder <- as.numeric(FigDat$PlOrder)
library(ggplot2)
attach(FigDat)
names(FigDat)

#scale mLs by Taxon
FigDat$mLS <- (FigDat$mLS / FigDat$TxMean)

#add a column of mmLS (mean of mean lesion size) per isolate
#sort dataframe by mmLS 
#then color by the new factor mmLS
FigDat$mmLS <- ave(FigDat$mLS, FigDat$IsolateID)
attach(FigDat)
FigDat <- FigDat[order(mmLS),]
ggplot(FigDat, aes(x = PlNum, y = mLS))+
  geom_point()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(mmLS), group=factor(IsolateID)), show.legend=F)+
  #scale_x_discrete(breaks=c("1","2","3", "4", "5", "6", "1", "2", "3", "4", "5", "6"),
                   #labels=c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", "LA0410", "LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))+
  facet_grid(~Taxon*Domest, scales="fixed", space="free_x")+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
#+geom_smooth(aes(group = 2), size = 2, method = "lm", se = T)
#----------------------------------------------------------------

