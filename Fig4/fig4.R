
#################
#################
# 10-07-2024
# Written by Abebe Fola 
#################
#################

#################
#################
# This is the code to generate combined Fig4 
#################
#################


getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Bioinformatics/Abefola_github/EPHI_SentinelSites_DR_project/Fig4") # change this to your working directory


##!1) read in the VCF file with vcfR and check the file
# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(adegenet)
library(ade4)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(ggplot2)
library(ggpubr)

IBC2COREL_vcf<- read.vcfR("borkesamples_biallelicsnp05q30sample05.recode.vcf", verbose = FALSE)

# get plotting data-------------------
samplessiBokre<- read.csv("Borkre_sites_map_coordinatesmod.csv", header=TRUE, sep=",")

head(samplessiBokre$Longitude)

### column names: 
#[1] "Nr"                         "R_EA_H_V_BARCODE_NAME"      "barcode"                   
#[4] "R_EA_H_V_BARCODE"           "plate"                      "row"                       
#[7] "col"                        "gbs_library_id"             "HHID_Surveyb"              
#[10] "hc_order"                   "sort_order"                 "cluster_groups"            
#[13] "counts"                     "sample"                     "variety"                   
#[16] "cluster"                    "matched_released_varieties" "status"                    
#[19] "lat"                   "long"                  "VarCountsbyHHID"

### critical ones are latitude and logitude plus whatever category you want to color by.


####################################################
####################################################
# Part II:  R script to check VCF file quality,read depth, to calculate MAF, SNP and Sample missingness rate
####################################################

##  Please install the R packages first
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("vcfR")
BiocManager::install("SeqArray")
BiocManager::install("adegenet")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")
BiocManager::install("SeqVarTools")

# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(adegenet)
library(ade4)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(ggplot2)
library(ggpubr)

## This is the version info of the R packages used by Fola et al. 
#adegenet_2.1.3 ade4_1.7-16   SeqArray_1.26.2    vcfR_1.12.0    tidyr_1.1.2   devtools_2.3.2  ggplot2_3.3.3 dplyr_1.0.3  SNPRelate_1.20.1 gdsfmt_1.22.0       

## The following are some sources as reference:    
# Use vcfR package for downstream analyis https://knausb.github.io/vcfR_documentation/index.html
# other packages for  population genomic anaysis in R https://grunwaldlab.github.io/Population_Genetics_in_R/Getting_ready_to_use_R.html


#############################
############################
# COI vs Site (Fig4A)
#############################
############################

coiBokreEPHIAshe <- read.csv("borkesuccseq_metadatacoi.csv")


ggplot(coiBokreEPHIAshe, aes(x = region, y = coimean, fill= region)) + 
  geom_violin(trim=FALSE)+
  #scale_fill_manual(values=site)+
  scale_fill_manual(values = c( "#009E73", "#0072B2", "#e6194b", "#f58231", "#CC79A7")) +
  ylab("COI") +
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("Complexity of infection Regional Level")




IBC2COREL_vcf<- read.vcfR("borkesamples_biallelicsnp05q30sample05.recode.vcf", verbose = FALSE)


#############################
############################
## PCA analysis - Fig4B
#############################
############################


input='borkesamples_biallelicsnp05q30sample05.recode.vcf'


#) Please change the name of the prefix of output files if needed
name = "pf_ephi"

# Define pops 



Bokre.site_col1 = c("#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                             "#e6194b",  "#fabebe", "darkred", "#46f0f0", "#3cb44b", "gray", "#008080")  
                             

site = c("#e6194b", "tomato4",  "#f58231", "#911eb4",    
                  "#f032e6",  "#fabebe", "#ffe119","#46f0f0","#4363d8","#008080", "#3cb44b", "gray")     
                  
site2 =  c( "#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                     "#e6194b",  "#fabebe", "darkred",  "#46f0f0", "#3cb44b", "gray", "#008080")
                     
shape.mutation <- c(5, 19)

# Pop1 Location
pop1 <-read.table("hc_modified.tsv", sep="\t")
A1<-colnames(pop1) <- c("sample.id",	"geo")

#Pop2 Year

pop2 <-read.table("mutation.tsv", sep="\t")

A1<-colnames(pop2) <- c("sample.id",	"mutation")


#Colors and shape 

# define vcf file 
vcf.fn <-input

ephi = paste('ephi_', name, '.gds', sep = "")
ephi

# Read the VCF file and save it as GDS format file
snpgdsVCF2GDS(vcf.fn, ephi,  method="biallelic.only")

## Open the SNP GDS file
genofile <- snpgdsOpen(ephi)

#calculate the eigenvectors and eigenvalues for principal component analysis.
ephi_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

#View(ephi_pca)
ephi_pca$varprop[1:10]*100


pc.percent <- ephi_pca$varprop*100
head(round(pc.percent, 2))

# eigenvalues
EV=plot(ephi_pca$eigenval[1:10], type = "o", col = "red", xlab = "PCA", ylab = "Eigenvalues",main = "Eigenvalues for EPHI data analysis")

par(mfrow=c(1, 1))
EV
# variance proportion (%)
plot(ephi_pca$varprop[1:10]*100, type = "o", 
     col = "black", siz=1, xlab = "PC", ylab = "(%)",
     main = "Bokre samples Variance Explained")

barplot(ephi_pca$varprop[1:10]*100, type = "o", col = "darkblue", xlab = "PC", ylab = "Percentage of variance explained",
        main = " variance explained ")

## Get data "sample.id" from a GDS node
## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

EV1 = ephi_pca$eigenvect[,1]    # the first eigenvector
EV2 = ephi_pca$eigenvect[,2]    # the second eigenvector
EV3 = ephi_pca$eigenvect[,3]    # the 3rd eigenvector
EV4 = ephi_pca$eigenvect[,4]    # the 4th eigenvector


## Get a list of pops following the order of the sample IDs and plot pca1 and 2

## Merge pop with you pca file

# Site

pop1_metadata<-read.table("hc_modified.tsv", sep="\t")
A1<-colnames(pop1_metadata) <- c("sample.id", "geo")

pop1 = factor(pop1_metadata$geo)[match(ephi_pca$sample.id, 
                                       pop1_metadata$sample.id)]
pop2 <-read.table("mutation.tsv", sep="\t")

A1<-colnames(pop2) <- c("sample.id",	"mutation")

pop2 = factor(pop2$mutation)[match(ephi_pca$sample.id, 
                                   pop2$sample.id)]

# Write pca file.tsv format
tab <- data.frame(ephi_pca$sample.id,pop1, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA1_2Bokre_final.tsv', sep=""), 
            quote = F, row.names = F, sep="\t")



# Larger Geographic areas


largersite = c( "#009E73", "#0072B2", "#e6194b", "#f58231", "#CC79A7") 

pop3 <-read.table("region_modified.tsv", sep="\t")
A1<-colnames(pop3) <- c("sample.id", "geo")

pop3 = factor(pop3$geo)[match(ephi_pca$sample.id, 
                              pop3$sample.id)] 

plot(EV1, EV2,xlab="PC1(6.3%)", ylab="PC2(3.8)", col=largersite[as.integer(pop3)],
     pch=19,cex=1.5)
legend("topright", legend=levels(pop3), bg="transparent",pch=16, cex=1,col=largersite,text.col=largersite)

abline(v=0.00, h=0.00, col="black", lwd=2, lty=2)

dev.off()


# Region and mutation combination 
# Write pca file.tsv format
tab <- data.frame(ephi_pca$sample.id,pop3, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA1_2Bokre_pop3final.tsv', sep=""), 
            quote = F, row.names = F, sep="/t")

plot(EV1, EV2,xlab="PC1(6.3%)", ylab="PC2(3.8)", col=largersite[as.integer(pop3)],
     pch=shape.mutation[as.integer(pop2)],cex=1.5)
legend("topright", legend=levels(pop3), bg="transparent",pch=16, cex=1,col=largersite,text.col=largersite)
legend("topleft", legend=levels(pop2), bg="transparent",pch=shape.mutation, cex=1)
abline(v=0, h=0, col="black", lwd=1.5, lty=2)








