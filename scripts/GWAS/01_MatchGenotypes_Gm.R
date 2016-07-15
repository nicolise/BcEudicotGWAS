#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in Gm phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcEudicotGWAS/data")
SNPs <- read.csv("SNP_files/hp_binaryMAF20.csv", row.names = 1)
#SNPs <- read.csv("miniSNP_practice.csv") 
#SNPsDF <- SNPs
#SNPsDF <- SNPsDF[c(1:2),]
#write.csv(SNPsDF, "SNPgenos.csv")
SNPs_rename <- SNPs

SNPnames <- read.csv("SNP_files/Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]
#read in model data with 2 trays dropped
Phenos <- read.csv("GmModelDataDT2.csv")

#change names from genotype file to match phenotype file
#File SNPs_rename has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

#change isolate names to match my coding
library(plyr)
unique(Phenos$isolatename)
unique(SNPnames$GenoRename)
#format: old = new
Phenos$isolatename <- revalue(Phenos$isolatename, c("MEA P6G" = "MEAP6G", "Kern A2" = "KernA2", "Pepper sub" = "PepperSub", "Fresa 525" = "Fresa525", "Gallo 2" = "Gallo2", "Apple 404" = "Apple404", "Fresa S.D." = "FresaSD", "BPA 1" = "BPA1", "Noble Rot" = "NobleRot", "Katie Tomato" = "KatieTomato", "Davis Navel" = "DavisNavel", "Apple 517" = "Apple517", "Kern B2" = "KernB2", "Triple 3 (T3)" = "T3", "Philo Menlo" = "PhiloMenlo", "Kern B1" = "KernB1", "Triple 7 (T7)" = "T7", "Esparato Fresa" = "EsparatoFresa", "UK razz " = "UKRazz", "Gallo 1" = "Gallo1"))

## now only keep genotypes and phenotypes that match

#save practice file
#miniSNPs <- SNPs_rename[c(1:3),]
#write.csv(miniSNPs, "miniSNP_practice.csv")

#only keep phenotype rows that match SNP names
SNPMt <- as.data.frame(names(SNPs_rename))
PhenoMatch <- Phenos
PhenoMatch <- PhenoMatch[PhenoMatch$isolatename %in% SNPMt$"names(SNPs_rename)", ]

#only keep SNP rows that match phenotype names
#select column with isolate names
PhenoMt <- as.data.frame(PhenoMatch[,14])
SNPMatch <- SNPs_rename
#save first 3 columns: won't match Pheno names
SNPs3 <- SNPs_rename[,c(1:3)]
SNPMatch <- SNPMatch[names(SNPMatch) %in% (PhenoMt$"PhenoMatch[, 14]")]
SNPMatch <- SNPMatch[ , order(names(SNPMatch))]
SNPMatch2 <- cbind(SNPs3,SNPMatch)

#sort pheno match
PhenoMatch2 <- PhenoMatch[order(PhenoMatch$isolatename),] 

#save the files
write.csv(SNPMatch2, "binSNP_bigRR_MAF20hp.csv")
write.csv(PhenoMatch2, "Sl_Pheno_bigRR.csv")
