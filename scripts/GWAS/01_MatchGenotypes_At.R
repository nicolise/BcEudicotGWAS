#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files/csvPrep")
setwd("~/Projects/BcSolGWAS/data/GWAS_files/02_csvprep")
SNPs <- read.csv("hp_binaryMAF20.csv", row.names = 1)
SNPs_rename <- SNPs

SNPnames <- read.csv("Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]

#change wd to get the phenotype file
setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS")
Phenos <- read.csv("LSMeanCamLes4Map_FIN.csv")
#this is the same as LSMeanCamLes4Map.csv except for a column added with renamed isolates

#change names from genotype file to match phenotype file
#File SNPs_rename has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

## now only keep genotypes and phenotypes that match

#only keep phenotype rows that match SNP names
SNPMt <- as.data.frame(names(SNPs_rename))
PhenoMatch <- Phenos
#90 isolates match between the genotyping and the BcAtGWAS phenotypes
PhenoMatch <- PhenoMatch[PhenoMatch$Rename %in% SNPMt$"names(SNPs_rename)", ]

#only keep SNP rows that match phenotype names
PhenoMt <- as.data.frame(PhenoMatch[,2])
SNPMatch <- SNPs_rename
SNPs3 <- SNPs_rename[,c(1:3)]
SNPMatch <- SNPMatch[names(SNPMatch) %in% (PhenoMt$"PhenoMatch[, 2]")]
SNPMatch <- SNPMatch[ , order(names(SNPMatch))]
SNPMatch2 <- cbind(SNPs3,SNPMatch)

#REMOVE genotype column with 01.01.06.1
#remove SNP column "X1.01.06.1"
SNPMatch2 <- SNPMatch2[,-9]

#sort pheno match
PhenoMatch2 <- PhenoMatch[order(PhenoMatch$Rename),] 

#save them files
write.csv(SNPMatch2, "binSNP_bigRR_MAF20hp.csv")
write.csv(PhenoMatch2, "At_Pheno_bigRR.csv")
#------------------------------------------------------------------------------
#extra things
#miniSNPs <- as.data.frame(t(miniSNPs))
#miniPhenos <- subset(Phenos, Igeno %in% SNPs_rename[0,])
