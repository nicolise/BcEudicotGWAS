#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcEudicotGWAS/")

#########################
# This makes the bigRR_update run through the GPU
# You need to do this first to mask the native 'bigRR_update' in the bigRR package
# one alternative to family = gaussian(link = identity) is family = poisson(link = log)
bigRR_update <- function (obj, Z, family = gaussian(link = identity), tol.err = 1e-06, 
                          tol.conv = 1e-08) 
{
  w <- as.numeric(obj$u^2/(1 - obj$leverage))
  w[w < tol.err] <- tol.err
  bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, 
        tol.err = tol.err, tol.conv = tol.conv, GPU = TRUE )
}
########################
#NOTE1 FROM RACHEL:  we need bigRR1.3-9 to get GPU option
# must download R-Forge version bigRR1.3-9tar.gz and manually install
# https://r-forge.r-project.org/R/?group_id=1301
# install.packages("bigRR", repos="http://R-Forge.R-project.org")

#NOTE2 FROM RACHEL: need package 'gputools' but CRAN version fails to install
# must first install Nvidia's CUDA toolkit -current version is 7.5
# installed from developer.nvidia.com/cuda-downloads

library(bigRR) #check if version is 1.3-9

#Get genotype data
SNPs <- read.csv("data/BcAtGWAS/03_bigRRinput/binSNP_bigRR_MAF20hp.csv", row.names = 1)
FullSNPs <- SNPs
SNPs <- FullSNPs
#add a column with position as chr.base
SNPs$Chr.Base <- do.call(paste, c(SNPs[c("X.CHROM","POS")], sep="."))

#Rachel's attempt------------------=====================================

rownames(SNPs) <- SNPs[,93] #set the new column of chrom.base as rownames - this could maybe be written as: rownames(SNPs) <- SNPs$Chr.Base?
SNPs <- SNPs[,4:92] #take out first three cols (X.CHROM, POS, REF) and new last col (Chr.Base). dim(SNPs) should now be [345485, 91], colnames(SNPs) are all Bc Isolates, rownames(SNPs) are all Chr.Base

#troubleshoot: try removing X1.02.19 because it wasn't in tomato GWAS
#SNPs <- SNPs[,-16]

SNPdf <- SNPs
#Rachel's attempt ----FIN-----============================================

#what does this do? 
#loads SNP dataframe as a matrix?
#makes SNP states numeric (also transposes SNP matrix)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
  SNPs[i,] <- as.numeric(SNPs[i,])
}

Phenos <- read.csv("data/BcAtGWAS/03_bigRRinput/At_Pheno_bigRR.csv", row.names = 1)
#note: in data frame replaced one NA (npr1.Cam, 01.02.02) with 7.91 (average for npr1)

#remove duplicate row 77 (KGB1 duplicate)
Phenos <- Phenos[-77,]

#add constant to remove negative values
Phenos$Col0.Cam <- Phenos$Col0.Cam + 9
Phenos$anac055.Cam <- Phenos$anac055.Cam + 4
Phenos$coi1.Cam <- Phenos$coi1.Cam + 2
Phenos$tga3.Cam <- Phenos$tga3.Cam + 3

#try standardizing Lesion Size Data (large values ~ 1000s)
Phenos$Col0.Les.s <- (Phenos$Col0.Les - mean(Phenos$Col0.Les))/sd(Phenos$Col0.Les) + 2
Phenos$anac055.Les.s <- (Phenos$anac055.Les - mean(Phenos$anac055.Les))/sd(Phenos$anac055.Les) +3
Phenos$coi1.Les.s <- (Phenos$coi1.Les - mean(Phenos$coi1.Les))/sd(Phenos$coi1.Les) +3
Phenos$npr1.Les.s <- (Phenos$npr1.Les - mean(Phenos$npr1.Les))/sd(Phenos$npr1.Les) +3
Phenos$pad3.Les.s <- (Phenos$pad3.Les - mean(Phenos$pad3.Les))/sd(Phenos$pad3.Les) +3
Phenos$tga3.Les.s <- (Phenos$tga3.Les - mean(Phenos$tga3.Les))/sd(Phenos$coi1.Les) +2

#only keep columns with standardized phenotypes
attach(Phenos)
Phenos <- Phenos[,c("Rename","Col0.Cam","Col0.Les.s","Col0.AT3G26830","Col0.AT2G30770" ,"Col.0AT4G30530","anac055.Cam","anac055.Les.s","coi1.Cam","coi1.Les.s","npr1.Cam","npr1.Les.s","tga3.Cam","tga3.Les.s")]

#troubleshoot: try removing X1.02.19 since it wasn't used in tomato GWAS
#Phenos <- Phenos[-16,]

#tried skipping first few phenos (6:14)
#complete is: 2:14
#error without permutations for "Col0.AT2G30770" (5)
dat <- as.data.frame((Phenos[,c(2:4,6:14)]))  #INSERT PHENOTYPE COLUMNS HERE
write.csv(dat, "mydata.csv")
dat <- read.csv("mydata.csv")
dat <- dat[,-1]


outpt.BLUP <- colnames(SNPs)
outpt.HEM <- colnames(SNPs)
thresh.BLUP <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)
thresh.HEM <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)

#Calculate BLUP and HEMs for all phenotypes
Sys.time()
#should be for (i in 1:dim(dat)[2])  but only trying with 1:3
for(i in 1:dim(dat)[2]) { #i will be each plant genotype
  print(colnames(dat)[i])
  MyX <- matrix(1, dim(dat)[1], 1)
  
  Pheno.BLUP.result <- bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE)
  Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)
  
  ##outpt.BLUP <- cbind(outpt.BLUP, Pheno.BLUP.result$u)
  outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
  
  #Permute Thresholds for Phenos - this is what takes forever
  ## troubleshoot: it's breaking here
  ##perm.u.BLUP <- vector()
  perm.u.HEM <- vector()
  #should be p in 1:1000
  for(p in 1:1000) {  
    if(p %% 10 == 0) {print(paste("Thresh sample:", p, "--", Sys.time()))}
    temp.Pheno <- sample(dat[,i], length(dat[,i]), replace = FALSE) ## skip to troubleshoot
    print(try(temp.BLUP  <- bigRR(y = temp.Pheno, X = MyX, Z = SNPs, GPU = TRUE),silent = TRUE))
    #works up to here with 3 i's and 3 p's
    try(temp.HEM <- bigRR_update(temp.BLUP, SNPs)) #REF change- was bigRR_update(Pheno.BLUP.result...
    ##perm.u.BLUP <- c(perm.u.BLUP, temp.BLUP$u) #REF change- ...c(perm.u.HEM...)
    perm.u.HEM <- c(perm.u.HEM, temp.HEM$u)
    
  }
  #write.csv(perm.u.HEM, paste("PermEffects_",colnames(dat)[i],".csv",sep=""))
  ##thresh.BLUP$"0.95Thresh"[i] <- quantile(perm.u.BLUP,0.95)
  ##thresh.BLUP$"0.975Thresh"[i] <- quantile(perm.u.BLUP,0.975)
  ##thresh.BLUP$"0.99Thresh"[i] <- quantile(perm.u.BLUP,0.99)
  ##thresh.BLUP$"0.999Thresh"[i] <- quantile(perm.u.BLUP,0.999)
  thresh.HEM$"0.95Thresh"[i] <- quantile(abs(perm.u.HEM),0.95)
  thresh.HEM$"0.975Thresh"[i] <- quantile(abs(perm.u.HEM),0.975)
  thresh.HEM$"0.99Thresh"[i] <- quantile(abs(perm.u.HEM),0.99)
  thresh.HEM$"0.999Thresh"[i] <- quantile(abs(perm.u.HEM),0.999)
  ##colnames(outpt.BLUP)[i+1] <- paste(colnames(dat)[i],"BLUP",sep=".")
  colnames(outpt.HEM)[i+1] <- paste(colnames(dat)[i],"HEM",sep=".")
}

#Give column names to the thresholds from the HEM list
for(j in 1:length(thresh.HEM)) {
  names(thresh.HEM[[j]]) <- colnames(dat)
}

#RF-give row names to thresh.HEM and thresh.BLUP so that threshhold values will line up correctly with phenotypes, and you can see which threshold value is displayed
##thresh.BLUP$"0.95Thresh" <- c("0.95 Thresh", thresh.BLUP$"0.95Thresh")
##thresh.BLUP$"0.975Thresh" <- c("0.975 Thresh", thresh.BLUP$"0.975Thresh")
##thresh.BLUP$"0.99Thresh" <- c("0.99 Thresh", thresh.BLUP$"0.99Thresh")
##thresh.BLUP$"0.999Thresh" <- c("0.999 Thresh", thresh.BLUP$"0.999Thresh")
thresh.HEM$"0.95Thresh" <- c("0.95 Thresh", thresh.HEM$"0.95Thresh")
thresh.HEM$"0.975Thresh" <- c("0.975 Thresh", thresh.HEM$"0.975Thresh")
thresh.HEM$"0.99Thresh" <- c("0.99 Thresh", thresh.HEM$"0.99Thresh")
thresh.HEM$"0.999Thresh" <- c("0.999 Thresh", thresh.HEM$"0.999Thresh")

#Write results to output
##write.csv(rbind(thresh.BLUP$"0.95Thresh",thresh.BLUP$"0.975Thresh",thresh.BLUP$"0.99Thresh",thresh.BLUP$"0.999Thresh",outpt.BLUP),"AthalianaLesion.BLUP.csv")
write.csv(rbind(thresh.HEM$"0.95Thresh",thresh.HEM$"0.975Thresh",thresh.HEM$"0.99Thresh",thresh.HEM$"0.999Thresh",outpt.HEM),"data/BcAtGWAS/04_bigRRoutput/AthalianaLesion.HEM.csv")
Sys.time()

#write.csv(outpt.HEM,"data/BcAtGWAS/04_bigRRoutput/Troubleshoot/AthalianaLesion_2.HEM.csv")
Sys.time()


### can stop here

# #This part failed for me -- NES
# used to store upper percentile of SNPS
# #Write just the positive positions (RF- effect size greater than .99 thresh)
# sig.HEM <- data.frame()
# for(i in 1:dim(outpt.HEM)[1]) {
#   if(i %% 1000 == 0) {print(paste(i, "--", Sys.time()))}
#   if(any(abs(as.numeric(outpt.HEM[i,2:5]))-abs(as.numeric(thresh.HEM$"0.99Thresh"[2:5]))>0)) { 
#     # change output.HEM[1,2:5] to appropriate column dims
#     sig.HEM <- unname(as.matrix(rbind(sig.HEM,outpt.HEM[i,])))
#   }
# }
# colnames(sig.HEM) <- colnames(outpt.HEM)
# write.csv(rbind(thresh.HEM$"0.99Thresh",sig.HEM),"SolanumLesionSizePoisson.HEM.99Sig.csv")

# Phenos[13,] -> Col0.LSMeans
# Phenos[-13,] -> Phenos
# 
# MyX <- matrix(1, dim(Phenos)[1], 1)
# 
# Pheno.BLUP.result <- bigRR(y = Phenos[,1], X = MyX, Z = SNPs)
# 
# Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)

#RF the code below generates plots comparing .BLUP and .HEM methods - not necessary

split.screen(c(1, 2))
split.screen(c(2, 1), screen = 1)
screen(3); plot(abs(Pheno.BLUP.result$u), cex = .6, col = 'slateblue')
screen(4); plot(abs(Pheno.HEM.result$u), cex = .6, col = 'olivedrab')
screen(2); plot(abs(Pheno.BLUP.result$u), abs(Pheno.HEM.result$u), cex = .6, pch = 19, col = 'darkmagenta')

split.screen(c(1, 2))
split.screen(c(2, 1), screen = 1)
screen(3); plot(abs(Pheno.BLUP.result$leverage), cex = .6, col = 'slateblue')
screen(4); plot(abs(Pheno.HEM.result$leverage), cex = .6, col = 'olivedrab')