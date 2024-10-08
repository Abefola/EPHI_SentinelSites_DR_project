

#################
#################
# 10-07-2024
# Written by Abebe Fola 
#################
#################

#################
#################
# This is the code to generate combined Fig2 
#################
#################


getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Bioinformatics/Abefola_github/EPHI_SentinelSites_DR_project/Fig3") # change this to your working directory

#################
# load all required libraries 
################
# UPset plot

#Install Library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)
library(UpSetR)

#Load the data 

MDR1_622I <- read.csv("MDR1_622I_genotypes.csv", header=TRUE, sep="," )

#Columns contain
# Sample Id, mutant names, and Wild genotype (if it contains wild alleles for all loci )

### Define colors 

#bar_cols <- c("darkcyan", "tomato4")

bar_cols3 <- c("#0072B2","darkcyan", "black", "grey","orange", "#CC79A7" ,"tomato4", "coral")


upset(MDR1_622I, nsets = 9, nintersects = 30, mb.ratio = c(0.5, 0.5), 
      order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE), keep.order = TRUE, main.bar.color = bar_cols3)

# Save plot 
ggsave("MDR1_622I.svg", dpi=600, width=7.5, height=7)
ggsave("MDR1_622I.pdf", dpi=600, width=7.5, height=7)

#################
#################
# k13_mdr1upsetplotfile
#################
#################


#Load the data 

k13_mdr1upsetplot <- read.csv("MDR1_675V_genotypes.csv", header=TRUE, sep="," )


#Columns contain
# Sample Id, mutant names, and Wild genotype (if it contains wild alleles for all loci )

### Define colors 

#bar_cols <- c("darkcyan", "tomato4")

#bar_cols <- c("#0072B2","darkgreen","darkcyan", "grey", "black", "red","tomato4","orange")

bar_cols2 <- c("#4363d8",  "#008080","#46f0f0","#f58231", "#911eb4",  "#f032e6",    
                        "#e6194b" )  
                        
upset(k13_mdr1upsetplot, nsets = 9, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE), keep.order = TRUE, main.bar.color = bar_cols2)


ggsave("MDR1_675V.svg", dpi=600, width=7.5, height=7)
ggsave("MDR1_675V.pdf", dpi=600, width=7.5, height=7)
