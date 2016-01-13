# GAPIT - Genomic Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# Written by Zhiwu Zhang, Alex Lipka and Feng Tian
# Last update: June 29, 2012

#Step 0: Set directory to load GAPIT source files (this step is omitted for using package)
#######################################################################################
library('MASS')
library(multtest)
library(gplots)
library(compiler) #for cmpfun


#Import library (each time to start R)
library(multtest)
library("gplots")
library("scatterplot3d")

wd <- getwd()
source(paste("/home/DNA/share_files/GWAS_NOV_2014/", "emma.txt", sep="/"))
source(paste("/home/DNA/share_files/GWAS_NOV_2014/", "gapit_functions.txt", sep="/"))
#source("http://www.zzlab.net/GAPIT/emma.txt")
#source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
#Create GAPIT directory . Set it as working directory
setwd(paste(wd, "gapit_results", sep="/"))

##
options <- commandArgs(trailingOnly=T)
pheno_file <- "/mnt/data2/1k_exome/gapit/extracted_phenotype.csv.foramtted"
pheno_index <- as.numeric(options[1])
##geno_file <- paste(wd, "merged_imputed.hmp.txt", sep="/")
geno_file_numeric <- "/mnt/data2/1k_exome/gapit/genotype.csv_formatted"
snp_info = "/mnt/data2/1k_exome/gapit/snp_info.txt"


#############################################################################################
#Import GAPIT


#Tutorial 4: Genome Prediction
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
#myY  <- read.table("traits_for_gapit.txt", head = TRUE)
myY  <- read.table(pheno_file, head = TRUE)
myY <- myY[, c(1, pheno_index)];
myGD <- read.table(geno_file_numeric, sep="\t", header=T)
myGM <- read.table(snp_info, sep="\t", header=T)
#Step 2: Run GAPIT
myGAPIT <- GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
Model.selection = TRUE,
SNP.test=TRUE
)

#save.image("gapit.RData")
