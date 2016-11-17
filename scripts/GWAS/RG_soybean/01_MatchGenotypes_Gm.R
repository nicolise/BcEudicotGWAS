#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in Gm phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcEudicotGWAS/data/RG_soybean/")
SNPs <- read.csv("02_csvprep/hp_binaryMAF20.csv", row.names = 1)
SNPs_rename <- SNPs

SNPnames <- read.csv("02_csvprep/Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,4)]
#read in model data with 2 trays dropped
Pheno.Les <- read.csv("02_csvprep/phenolesion.csv")
Pheno.OA <- read.csv("02_csvprep/pheno_OA.csv")

#change names from genotype file to isolate names
#File SNPs_rename has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

#change names from phenotype file to isolate names
Pheno.Les$REF <- SNPnames[ match( Pheno.Les$REF , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 
Pheno.OA$REF <- SNPnames[ match( Pheno.OA$REF , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

library(plyr)
unique(Pheno.Les$REF)
unique(Pheno.OA$REF)
unique(SNPnames$GenoRename)

## now only keep genotypes and phenotypes that match

#save practice file
#miniSNPs <- SNPs_rename[c(1:3),]
#write.csv(miniSNPs, "miniSNP_practice.csv")

#only keep phenotype rows that match SNP names
SNPMt <- as.data.frame(names(SNPs_rename))
PhenoMatch.Les <- Pheno.Les
PhenoMatch.OA <- Pheno.OA
PhenoMatch.Les <- PhenoMatch.Les[PhenoMatch.Les$REF %in% SNPMt$"names(SNPs_rename)", ]
PhenoMatch.OA <- PhenoMatch.OA[PhenoMatch.OA$REF%in% SNPMt$"names(SNPs_rename)", ]

#only keep SNP rows that match phenotype names
#select column with isolate names
PhenoMt.Les <- as.data.frame(PhenoMatch.Les[,1])
PhenoMt.OA <- as.data.frame(PhenoMatch.OA[,1])
SNPMatch.Les <- SNPs_rename
SNPMatch.OA <- SNPs_rename

#save first 3 columns: won't match Pheno names
SNPs3 <- SNPs_rename[,c(1:3)]
SNPMatch.Les <- SNPMatch.Les[names(SNPMatch.Les) %in% (PhenoMt.Les[,1])]
SNPMatch.OA <- SNPMatch.OA[names(SNPMatch.OA) %in% (PhenoMt.OA[,1])]
SNPMatch.Les <- SNPMatch.Les[ , order(names(SNPMatch.Les))]
SNPMatch.OA <- SNPMatch.OA[ , order(names(SNPMatch.OA))]
SNPMatch2.Les <- cbind(SNPs3,SNPMatch.Les)
SNPMatch2.OA <- cbind(SNPs3,SNPMatch.OA)

#sort pheno match
PhenoMatch2.Les <- PhenoMatch.Les[order(PhenoMatch.Les$REF),] 
PhenoMatch2.OA <- PhenoMatch.OA[order(PhenoMatch.OA$REF),] 

#save the files
write.csv(SNPMatch2.OA, "03_bigRRinput/binSNP_bigRR_MAF20hp.OA.csv")
write.csv(SNPMatch2.Les, "03_bigRRinput/binSNP_bigRR_MAF20hp.Les.csv")
write.csv(PhenoMatch2.Les, "03_bigRRinput/Gm_Les.Pheno_bigRR.csv")
write.csv(PhenoMatch2.OA, "03_bigRRinput/Gm_OA.Pheno_bigRR.csv")
