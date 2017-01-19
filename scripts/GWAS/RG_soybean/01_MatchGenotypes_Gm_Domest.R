#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in Gm phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcEudicotGWAS/data/RG_soybean/")
SNPs <- read.csv("02_csvprep/Domest/hp_binaryMAF20.csv", row.names = 1)
SNPs_rename <- SNPs

SNPnames <- read.csv("02_csvprep/Domest/Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,4)]
#read in model data with 2 trays dropped
Pheno.D <- read.csv("02_csvprep/Domest/RG_phenoDandW.csv")

#change names from genotype file to isolate names
#File SNPs_rename has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

#change names from phenotype file to isolate names
Pheno.D$REF <- SNPnames[ match( Pheno.D$REF , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

library(plyr)
unique(Pheno.D$REF)
unique(SNPnames$GenoRename)

## now only keep genotypes and phenotypes that match

#save practice file
#miniSNPs <- SNPs_rename[c(1:3),]
#write.csv(miniSNPs, "miniSNP_practice.csv")

#only keep phenotype rows that match SNP names
SNPMt <- as.data.frame(names(SNPs_rename))
PhenoMatch.D <- Pheno.D
PhenoMatch.D <- PhenoMatch.D[PhenoMatch.D$REF%in% SNPMt$"names(SNPs_rename)", ]

#only keep SNP rows that match phenotype names
#select column with isolate names
PhenoMt.D <- as.data.frame(PhenoMatch.D[,1])
SNPMatch.D <- SNPs_rename

#save first 3 columns: won't match Pheno names
SNPs3 <- SNPs_rename[,c(1:3)]
SNPMatch.D <- SNPMatch.D[names(SNPMatch.D) %in% (PhenoMt.D[,1])]
SNPMatch.D <- SNPMatch.D[ , order(names(SNPMatch.D))]
SNPMatch2.D <- cbind(SNPs3,SNPMatch.D)
#remove 01.01.06 duplicate
SNPMatch2.D <- SNPMatch2.D[,-9]

#sort pheno match
PhenoMatch2.D <- PhenoMatch.D[order(PhenoMatch.D$REF),] 
#remove 01.01.06 duplicate
PhenoMatch2.D <- PhenoMatch2.D[-6,]

#troubleshoot matching
genosList <- names(SNPMatch2.D)
removeList <- c("X.CHROM", "POS", "REF")
genosList <- setdiff(genosList,removeList)
myCheck <- data.frame(PhenoMatch2.D$REF, genosList)
setdiff(genosList,PhenoMatch2.D$REF)
setdiff(names(SNPMatch2.D), PhenoMatch2.D$REF)
#looks fine


#save the files
write.csv(SNPMatch2.D, "03_bigRRinput/Domest/binSNP_bigRR_MAF20hp.D.csv")
write.csv(PhenoMatch2.D, "03_bigRRinput/Domest/Gm_Les.Pheno_bigRR.D.csv")
