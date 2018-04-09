#Nicole E Soltis
#04/07/18
#A01_TABtoPEDnMAP.R

#---------------------------------------------------------------------------
rm(list=ls())

#using same genotype input from T4 BcAtGWAS bigRR 
setwd("~/Projects/BcGenome/data")
setwd("~/Documents/GitRepos/BcGenome/data")
SNPnames <- read.csv("Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]
names(SNPnames)[1]<- "Isolate"

#read in edited files: T4 from BcAt bigRR
#see script BcAt_RNAGWAS/scripts/bigRRGWAS/allreads/02_PrepGenos_allreads_T4.R
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/")
mySNPs <- read.csv("T4/01_prepFiles/hp_binMAF20_20NA.csv")
SNPs_renamed <- mySNPs

setwd("~/Projects/BcEudicotGWAS")
#here are the phenotypes prior to lsmeans
Phenos <- read.csv("data/BcAtGWAS/03_bigRRinput/At_Pheno_bigRR.csv", row.names = 1)
Phenos <- Phenos[,-c(1)]

#change names from genotype file to match phenotype file
#File SNPs_renamed has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_renamed ) <- SNPnames[ match( names( SNPs_renamed ) , SNPnames[ , 'SNPname' ] ) ,'Isolate' ] 

## now only keep genotypes and phenotypes that match
#only keep phenotype rows that match SNP names
IsosfromSNP <- as.data.frame(names(SNPs_renamed))
names(IsosfromSNP)[1] <- "IsoList"
matchedPhenos <- Phenos
intersect(matchedPhenos$Isolate, IsosfromSNP$IsoList)
setdiff(matchedPhenos$Isolate, IsosfromSNP$IsoList)
sort(unique(IsosfromSNP$IsoList))
#have no SNPs for 1.02.13, that drops
#and rename matchedPhenos$Isolate MEAPGG to MEAP6G
levels(matchedPhenos$Isolate) <- c(levels(matchedPhenos$Isolate),"MEAP6G")
matchedPhenos$Isolate[matchedPhenos$Isolate=="MEAPGG"] <- "MEAP6G"
matchedPhenos <- matchedPhenos[matchedPhenos$Isolate %in% IsosfromSNP$IsoList, ]
#yay, only drops 01.02.13


#only keep SNP rows that match phenotype names
IsosfromPheno <- as.data.frame(matchedPhenos[,1])
matchedSNPs <- SNPs_renamed
mySNPs <- mySNPs[,-c(1)]
SNPs3 <- mySNPs[,c(1:3)]
matchedSNPs <- matchedSNPs[names(matchedSNPs) %in% (matchedPhenos$"Isolate")]
matchedSNPs <- matchedSNPs[, order(names(matchedSNPs))]
matchedSNPs2 <- cbind(SNPs3, matchedSNPs)

#sort pheno match
matchedPhenos2 <- matchedPhenos[order(matchedPhenos$Isolate),] 

#check for matching names between SNPMatch2 and PhenoMatch2
setdiff((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(4:99)]))
intersect((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(4:99)])) #good

write.csv(matchedSNPs2, "T4/GEMMA/")
write.csv(matchedPhenos2, "T4/")

#and now for making PED format for PLINK!
#do not need positional info: just SNP states for PED
#turn df sideways (individuals as rows, SNPs as columns)
#split each genotype into 2 identical columns (PED assumes diploid)
#add a first column: FAM1 (no info on isolate families)
#second column: isolate ID
#third column: father ID (a column of zeros)
#fourth column: mother ID (a column of zeros)
#fifth column: individual sex = 1 (all assumed same)
#sixth  column: binary  phenotype (all = 1)
#fix column order



#get rid of Chrom, Contig, Pos
mySNPs2 <- matchedSNPs2[,-c(1:3)]

#for PED, NA must be replaced with 0 for genotypes, else NA will be read as an allele
#so first, set all genotypes = 0 to =2
mySNPs2[mySNPs2==0] <- 2
mySNPs2[is.na(mySNPs2)] <- 0

#turn all SNPs to "diploid"
#haha, it takes 4 days to do this as a "for" loop (for each row, rbind twice)
#because is.na <-0 before this step, there should be NO heterozygous SNP calls
#this is super fast:
mySNPs3 <- mySNPs2[rep(1:nrow(mySNPs2),each=2),] 

#transpose and format for PED
mySNPs4 <- as.data.frame(t(mySNPs3))
#add binary phenotype = 1 (6)
mySNPs4 <- cbind("Pheno" = 1, mySNPs4)
#add individual sex = 1 (5)
mySNPs4 <- cbind("sex" = 1, mySNPs4)
#add Mother = 0 (4)
mySNPs4 <- cbind("Mother" = 0, mySNPs4)
#add Father = 0 (3)
mySNPs4 <- cbind("Father" = 0, mySNPs4)
#turn row names into column 2
mySNPs4 <- cbind(rownames(mySNPs4), mySNPs4)
colnames(mySNPs4)[1] <- 'Isolate'
#add the fam column (1)
mySNPs4 <- cbind("FAM" = "FAM1", mySNPs4)
myPED <- mySNPs4

#add a phenotype for PED? 
#NA is fine for missing phenotypes
#since many phenotypes, just add as consecutive columns to *.fam, and run GEMMA in a loop over phenotypes

#make a MAP file for plink (need it to make the bed (binary ped) file from ped)
myMAP <- mySNPs[,c("Chrom","Pos")]
myMAP2 <- myMAP
myMAP2$SNPID <- paste("SNP",myMAP2$Pos, sep="")
myMAP2$SNPcM <- 0
myMAP2 <- myMAP2[,c(1,3,4,2)]
setwd("~/Documents/GitRepos/BcSolGWAS/")
#MAP2 still has chromosome 1:18
write.table(myMAP2, "data/GEMMA_files/D_01_PLINK/dpbinMAF20NA10.map", row.names=FALSE, col.names=FALSE)
#add a column of "SNP identifiers" in excel?

write.csv(mySNPs3, "data/GEMMA_files/D_01_PLINK/dp_binMAF20_10NA.csv")
write.csv(mySNPs, "data/GEMMA_files/D_01_PLINK/hp_binMAF20_10NA.csv")
Sys.time()
write.table(myPED, "data/GEMMA_files/D_01_PLINK/dpbinMAF20NA10.ped", row.names=FALSE, col.names=FALSE)
Sys.time()
