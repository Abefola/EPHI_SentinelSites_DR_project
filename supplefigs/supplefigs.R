
#################
#################
# 10-07-2024
# Written by Abebe Fola 
#################
#################

#################
#################
# This is the code to generate all supplementary figs
#################
#################



#################
#################

#Figure S1. Sample Collection, Processing and Data Analysis Workflow for P. falciparum Drug Resistance Study in Ethiopia
# This figure was generated using adobe and powerpoint
#################
#################

#################
#################
# Fig S2 retained Sample coverage and SNP distribution  
#################
#################



getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Bioinformatics/Abefola_github/EPHI_SentinelSites_DR_project/supplefigs")

# sample distributions

##################################################
# google map method (ggmap)
##################################################

getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/2024_focus-areas/EPHI_622I_project_phaseone/repooledIBC2CORE")
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

IBC2COREL_vcf<- read.vcfR("borkesamples_biallelicsnp05q30sample05.recode.vcf", verbose = FALSE)

##!2) vcf Summarization and check metadata

head(IBC2COREL_vcf) 
queryMETA(IBC2COREL_vcf) 
queryMETA(IBC2COREL_vcf, element = 'DP')
head(is.polymorphic(IBC2COREL_vcf, na.omit = TRUE))
head(is.biallelic(IBC2COREL_vcf))
queryMETA(IBC2COREL_vcf, element = 'FORMAT=<ID=DP')
strwrap(IBC2COREL_vcf@meta[1:7])

# If you needed you can Subset samples
IBC2COREL_vcf[,1:10]

# The fix region
# The fix region contains information for each variant which is sometimes summarized over all samples. The first eight columns of the fixed region and are titled CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO. This is per variant information which is ‘fixed’, or the same, over all samples. The first two columns indicate the location of the variant by chromosome and position within that chromosome. Here, the ID field has not been used, so it consists of missing data (NA). The REF and ALT columns indicate the reference and alternate allelic staBokre. When multiple alternate allelic staBokre are present they are delimited with commas. The QUAL column attempts to summarize the quality of each variant over all samples. The FILTER field is not used here but could contain information on whether a variant has passed some form of quality assessment.
tail(getFIX(IBC2COREL_vcf))

# The gt region
#The gt (genotype) region contains information about each variant for each sample. The values for each variant and each sample are colon delimited. Multiple types of data for each genotype may be stored in this manner. The format of the data is specified by the FORMAT column (column nine). Here we see that we have information for GT, AD, DP, GQ and PL. The definition of these acronyms can be referenced by querying the the meta region, as demonstrated previously. Every variant does not necessarily have the same information (e.g., SNPs and indels may be handled differently), so the rows are best treated independently. Different variant callers may include different information in this region.
IBC2COREL_vcf@gt[1:6, 1:4]

##!4)Creating chromR objects and plots

#Then plot important statistics summed over entire VCF
chrom_Bokre_samples <- create.chromR(vcf=IBC2COREL_vcf)



#quick check cumulative read depth distribution across samples - Fig S2A


DP_samples <- extract.gt(IBC2COREL_vcf, element='DP', as.numeric=TRUE)
rownames(DP_samples ) <- 1:nrow(DP_samples )

head(DP_samples )

heatmap.bp(DP_samples)

is.na(DP_samples [na.omit(DP_samples  == 0)]) <- TRUE

heatmap.bp(log(DP_samples), cbarplot = F, rbarplot = F, min= 0.2, max =0.5, legend = TRUE, clabels = F, rlabels = TRUE, na.rm = TRUE,
           scale = c("column"),
           #col.ramp = viridisLite::viridis(n = 100, alpha = 1))
           #col.ramp = colorRampPalette(c("yellow", "orange", "red"))(100))
           
           col.ramp = colorRampPalette(c("grey", "darkred"))(110))

par(mar=c(5,4,4,2))
boxplot(DP_samples , col=2:8, las=3)
title(ylab = " Depth (DP)")




# Load SNP density formatted data _ Fig S2B


samples855_snp_density<- read.csv("bokresamplessnp_coordinates.csv", header = T) # Look the attached file. This one of intermediate file for selection analysis. It contains four columns Snp_coordinate, Chro_no, SNP_position and P_value (not required for SNP density plot)

#plot SNP-density plot


library(CMplot)
# Sample plot, make sure `samples855_snp_density` is correctly formatted
CMplot(samples855_snp_density, 
       plot.type = "d", 
       bin.size = 50000, 
       chr.den.col = c("darkgreen", "yellow", "red"), 
       file = "jpg", 
       file.output = TRUE, 
       verbose = TRUE, 
       width = 9, 
       height = 6, 
       dpi = 300)

CMplot(samples855_snp_density, 
       plot.type = "d", 
       bin.size = 25000, 
       chr.den.col = c( "#008080", "#008080", "#008080"), 
       file = "jpg",                # File type (pdf)
       file.output = TRUE, 
       verbose = TRUE, 
       width = 9, 
       height = 6, 
       dpi = 300)


##############
##############
# Fig S3 - Probe performance
##############
##############
data <- read.csv("sequence perforemancemod.csv")

# Calculate mean coverage for each probe
results <- data %>%
  select(-samples) %>%  # Exclude the sample column
  summarise(across(everything(), list(
    mean = ~mean(.)
  ), .names = "{col}_{fn}"))

# Reshape the results for better presentation
results_long <- results %>%
  pivot_longer(cols = everything(), 
               names_to = c("Probe", "Statistic"), 
               names_pattern = "(.*)_(.*)")


# Filter the results to only keep the 'mean' statistic
mean_data <- results_long %>% filter(Statistic == "mean")

# Categorize the mean coverage into the specified ranges
mean_data <- mean_data %>%
  mutate(Coverage_Category = case_when(
    value < 10 ~ "< 10",
    value >= 10 & value < 50 ~ "10 to 50",
    value >= 50 & value < 100 ~ "50 to 100",
    value >= 100 ~ "> 100"
  ))

# Count the number of probes and calculate the proportion in each category
category_stats <- mean_data %>%
  group_by(Coverage_Category) %>%
  summarise(
    Count = n(),
    Proportion = n() / nrow(mean_data)
  ) %>%
  mutate(Label = paste(Count, "(", round(Proportion * 100, 1), "%)", sep = ""))

# Merge the category statistics into the main data
mean_data <- mean_data %>%
  left_join(category_stats, by = "Coverage_Category")

# Color Universal Design (CUD) color palette (color-blind friendly)
colors <- c("< 10" = "#E69F00", 
            "10 to 50" = "#56B4E9", 
            "50 to 100" = "#009E73", 
            "> 100" = "#F0E442")

# Plot the Mean Coverage with color-blind friendly colors and custom legend
ggplot(mean_data, aes(x = Probe, y = value, fill = Coverage_Category)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = colors) +
  theme_grey() +
  labs(title = "Mean Coverage by Probe",
       x = "Probe",
       y = "Mean Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = NULL) +  # Optional: remove x-axis labels if too cluttered
  guides(
    fill = guide_legend(title = "Coverage Category", 
                        labels = category_stats$Label)  # Custom legend labels
  )



##############
##############
# Figure S4. WHO validated K13 mutations prevalence health facility level 
##############
##############

df_freq_comb_year <- read.csv("allefre_Ethiopia_2019_2023.csv")


cbp3 <- c("#46f0f0","#008080","#fabebe", "#4363d8")

k13mutationcomb <- read.csv("k13mutationcombmod.csv")

# plot
ggplot(k13mutationcomb,aes(x=rehc, y=prevalence, fill = marker))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Mutation name", y="Prevalence (%)") +
  geom_col(width = 0.5, position = position_dodge(0.7))+
  scale_fill_discrete(name="Markers") +
  theme(axis.text.x=element_text(size=rel(1.2), angle=270))+
  #theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #scale_fill_jco()+
  scale_fill_manual(values=cbp3)+
  ggtitle("K13mutations per HC")


##############
##############
#Figure S5. Prevalence of non-synonymous mutations across the K13 gene, coloured according to amino-acid residues within beta-propeller domain where validated resistance mutations are located or not. 
##############
##############

FigS5 <- read.csv("K13_mutations.csv", header=TRUE, sep="," )

# plot
cbp2 <- c("darkred", "#0072B2") # my color choice 

kelch13 <- ggplot(FigS5, aes(x=reorder(Mutation, AA_coden_position), y=prevalence, fill = WHO_Category
))
kelch13  + 
  geom_bar(stat = "identity")+
  labs(x="Mutation name", y="Prevalence (%)") +
  scale_fill_discrete(name="WHO Category") +
  theme(axis.text.x=element_text(size=rel(1.2), angle=270))+
  #theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #scale_fill_jco()+
  scale_fill_manual(values=cbp2)
ggtitle("NS mutations across Kelch13 gene in Ethiopia")



##############
##############
# Figure S6. Network analysis showing highly related parasite pairs (IBD≥0.99) of 622I mutants (A) regional level and (B) health facility (district) level
##############
##############

onlymutants<-read.csv("final_ibd_mle_long_onlymutants.csv", header=TRUE)
set.seed(100)

dbs2 <- onlymutants

dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.90)
relations <- dbs3 %>%
  dplyr::select(p1,p2,deme_p1)

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,deme_p1, deme_p2)


dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.90)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,IBD,deme_p1)

for.checks <- dbs3 %>%
  dplyr:: select(p1,p2,IBD,deme_p1,deme_p2)

relations <- dbs3 %>%
  dplyr::select(p1,p2,IBD,deme_p1)

colnames(relations) <- c("from","to","F","deme_p1")

relations <- relations[order(relations$deme_p1),]
relations$deme_p1 <- as.factor(relations$deme_p1)


dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.90)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,IBD,deme_p1)

colnames(relations) <- c("target","source","F","deme_p1")
relations <- relations[order(relations$deme_p1),]

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,IBD,deme_p1,deme_p2)
colnames(for.checks) <- c("target","source","F","deme_p1","deme_p2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

links90MUTANTS <- relations %>% dplyr::select(target, source, F)
write.csv(links90MUTANTS, "links90MUTANTSprec.csv" )

#links90<- read.csv("links95mod1prec.csv", header = T)
seq1 <- for.checks %>% dplyr::select(target,deme_p1) 
seq2 <- for.checks %>% dplyr::select(source,deme_p2) %>% rename(target=source,deme_p1=deme_p2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,deme_p1)
nodes <- nodes[order(nodes$deme_p1),]
nodes <- droplevels(nodes)

metadcomplete <- read.csv("final_reference_Bokre_metadatacoimod.csv")

ta90MUTANT <- metadcomplete %>% 
  dplyr::filter(metadcomplete$p1 %in% nodes$target) %>%
  dplyr::select(p1, deme_p1)


write.csv(ta90MUTANT, "ta90MUTANTprec.csv") # modify and order per pop
#ta95<-read.csv("ta95msmtprec.csv", header = T)

#ta2 <- ta %>%
# select(seqID, hrp23_status_f)


IBDNW90MUTANT <- graph_from_data_frame(d=links90MUTANTS, vertices=ta90MUTANT, directed=F) 


#Figure S6A
colorsregion <- c(rep("#0072B2",2),   rep("#009E73", 42), rep("#e6194b", 2),rep("#009E73", 41))

plot(IBDNW90MUTANT, vertex.label=NA,vertex.size=7, main="IBD>=0.99, n=500", vertex.color=colorsregion)

legend(x="topleft", legend=c("Gamnbela","Amhara","Oromia", "Amhara"), col=c( "#0072B2", "#009E73", "#e6194b", "#009E73"), cex=1.5, pch=c(19))


#Figure S6B
colorsHC <- c(rep("#f032e6",2),   rep("#4363d8", 42), rep("#f58231", 2),rep("#fabebe", 27), rep("#008080", 14))

plot(IBDNW90MUTANT, vertex.label=NA,vertex.size=8, main="IBD>=0.99, n=500", vertex.color=colorsHC)

legend(x="topleft", legend=c("Abol","Andansa","Asendabo", "Jiga", "Woreta"), col=c( "#f032e6", "#4363d8", "#f58231",
                                                                                             "#fabebe",  "#008080"), cex=1.5, pch=c(19))



####################

# Figure S7. Network analysis showing relatedness of mutant vs wildtype parasite pairs with different IBD thresholds 
#(IBD≥0.75,  A , and IBD≥0.50, B. 
####################

inf.m7<-read.csv("final_ibd_mle_long_622I.csv", header=TRUE)
set.seed(100)

#######################
#Figure S7A - 0.75
######################

dbs2 <- inf.m7

dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.75)
relations <- dbs3 %>%
  dplyr::select(p1,p2,k13_Arg622Ile_p1)

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,k13_Arg622Ile_p1, k13_Arg622Ile_p2)


dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.75)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1)

for.checks <- dbs3 %>%
  dplyr:: select(p1,p2,IBD,k13_Arg622Ile_p1,k13_Arg622Ile_p2)

relations <- dbs3 %>%
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1)

colnames(relations) <- c("from","to","F","k13_Arg622Ile_p1")

relations <- relations[order(relations$k13_Arg622Ile_p1),]
relations$k13_Arg622Ile_p1 <- as.factor(relations$k13_Arg622Ile_p1)


dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.75)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1)

colnames(relations) <- c("target","source","F","k13_Arg622Ile_p1")
relations <- relations[order(relations$k13_Arg622Ile_p1),]

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1,k13_Arg622Ile_p2)
colnames(for.checks) <- c("target","source","F","k13_Arg622Ile_p1","k13_Arg622Ile_p2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

links75 <- relations %>% dplyr::select(target, source, F)
#write.csv(links75, "links75prec.csv" )

#links90<- read.csv("links95mod1prec.csv", header = T)
seq1 <- for.checks %>% dplyr::select(target,k13_Arg622Ile_p1)
seq2 <- for.checks %>% dplyr::select(source,k13_Arg622Ile_p2) %>% rename(target=source,k13_Arg622Ile_p1=k13_Arg622Ile_p2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,k13_Arg622Ile_p1)
nodes <- nodes[order(nodes$k13_Arg622Ile_p1),]
nodes <- droplevels(nodes)

metadcomplete <- read.csv("final_Bokre_k13_metadata.csv")

ta75 <- metadcomplete %>%
  dplyr::filter(metadcomplete$p1 %in% nodes$target) %>%
  dplyr::select(p1, k13_Arg622Ile_p1)



#write.csv(ta75, "ta75prec.csv") # modify and order per pop
#ta95<-read.csv("ta95msmtprec.csv", header = T)

#ta2 <- ta %>%
# select(seqID, hrp23_status_f)


IBDNW75 <- graph_from_data_frame(d=links75, vertices=ta75, directed=F)

colors <- c(rep("#008080",488),   rep("#CC79A7", 102))

plot(IBDNW75, vertex.label=NA,vertex.size=7, main="IBD>=0.75", vertex.color=colors)

legend(x="topleft", legend=c("Wildtype","Mutant"), col=c( "#008080", "#CC79A7"), cex=0.7, pch=c(19))


#######################
#Figure S7B - 0.50
######################

dbs2 <- inf.m7

dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.50)
relations <- dbs3 %>%
  dplyr::select(p1,p2,k13_Arg622Ile_p1)

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,k13_Arg622Ile_p1, k13_Arg622Ile_p2)


dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.50)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1)

for.checks <- dbs3 %>%
  dplyr:: select(p1,p2,IBD,k13_Arg622Ile_p1,k13_Arg622Ile_p2)

relations <- dbs3 %>%
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1)

colnames(relations) <- c("from","to","F","k13_Arg622Ile_p1")

relations <- relations[order(relations$k13_Arg622Ile_p1),]
relations$k13_Arg622Ile_p1 <- as.factor(relations$k13_Arg622Ile_p1)


dbs3 <- dbs2 %>%
  dplyr:: filter(IBD >= 0.50)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1)

colnames(relations) <- c("target","source","F","k13_Arg622Ile_p1")
relations <- relations[order(relations$k13_Arg622Ile_p1),]

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,IBD,k13_Arg622Ile_p1,k13_Arg622Ile_p2)
colnames(for.checks) <- c("target","source","F","k13_Arg622Ile_p1","k13_Arg622Ile_p2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

links05 <- relations %>% dplyr::select(target, source, F)
write.csv(links05, "links05prec.csv" )

#links90<- read.csv("links95mod1prec.csv", header = T)
seq1 <- for.checks %>% dplyr::select(target,k13_Arg622Ile_p1)
seq2 <- for.checks %>% dplyr::select(source,k13_Arg622Ile_p2) %>% rename(target=source,k13_Arg622Ile_p1=k13_Arg622Ile_p2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,k13_Arg622Ile_p1)
nodes <- nodes[order(nodes$k13_Arg622Ile_p1),]
nodes <- droplevels(nodes)

metadcomplete <- read.csv("final_Bokre_k13_metadata.csv")

ta05 <- metadcomplete %>%
  dplyr::filter(metadcomplete$p1 %in% nodes$target) %>%
  dplyr::select(p1, k13_Arg622Ile_p1)
IBDNW05 <- graph_from_data_frame(d=links05, vertices=ta05, directed=F)


colors <- c(rep("#008080",448),   rep("#CC79A7", 89))

plot(IBDNW05, vertex.label=NA,vertex.size=7, main="IBD>=0.50", vertex.color=colors)

legend(x="topleft", legend=c("Wildtype","Mutant"), col=c( "#008080", "#CC79A7"), cex=0.7, pch=c(19))




####################

# Figures S8. Pairwise IBD at health center level. Boxes indicate the interquartile range, the line indicates the median, the whiskers show the 95% confidence intervals, and the dots show outlier values. 

####################

#1. load the file 
data <- read.csv("final_ibd_mle_long_622I.csv", stringsAsFactors = FALSE)


# 2. Filter for IBD > 0 and where deme_p1 == deme_p2 (same health centre)
data_filtered <- inf.m8 %>%
  filter(IBD > 0, deme_p1 == deme_p2)

# 3. Create a new 'mutation_status' column reflecting the mutation status of each pair
data_filtered$mutation_status <- dplyr::case_when(
  data_filtered$deme_p1 == "mutant" & data_filtered$deme_p1 == "mutant" ~ "mutant_mutant",   # Both are mutant
  data_filtered$deme_p1 == "wildtype" & data_filtered$deme_p1 == "wildtype" ~ "wildtype_wildtype",   # Both are wildtype
  TRUE ~ "wildtype_mutant"  # One is wildtype and the other is mutant
)

# 4. Visualize the pairwise IBD sharing for IBD > 0, stratified by mutation status within the same health centre
ggplot(data_filtered, aes(x = mutation_status, y = IBD, fill = mutation_status)) +
  geom_boxplot() +
  facet_wrap(~ deme_p1, scales = "free") +  # Facet the plot by deme_p1 (health centre), scales = "free" allows each health centre to have its own axis range
  labs(
    title = "Pairwise IBD Sharing by Mutation Status Within Health Centres (IBD > 0)",
    x = "Mutation Status",
    y = "IBD Sharing"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability


# 4. Custom colors for visualization
custom_colors <- c("mutant_mutant" = "#CC79A7",  # Teal for mutant_mutant
                   "wildtype_wildtype" = "#008080",  # Blue for wildtype_wildtype
                   "wildtype_mutant" = "#E69F00")  # Pink for wildtype_mutant

# 5. Visualize the pairwise IBD sharing for IBD > 0, stratified by mutation status within the same health centre
ggplot(data_filtered, aes(x = mutation_status, y = IBD, fill = mutation_status)) +
  geom_boxplot() +
  facet_wrap(~ deme_p1, scales = "free") +  # Facet the plot by deme_p1 (health centre), scales = "free" allows each health centre to have its own axis range
  labs(
    title = "Pairwise IBD Sharing by Mutation Status Within Health Centres (IBD > 0)",
    x = "Mutation Status",
    y = "IBD Sharing"
  ) +
  scale_fill_manual(values = custom_colors) +  # Apply the custom color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
    legend.position = "none"  # Remove legend since mutation_status is already labeled on the x-axis
  )





################

# Fig S9 A, B AND C 

#############

####
#PLOT > Zero per site
####

site = c( "#CC79A7", "#008080")

ibd_mle_long_mtdt_622I <- read.csv("FigS9.csv")

Totalpairsgreater0 <- ibd_mle_long_mtdt_622I %>%
  filter (malecotf>0.0)



my_comparisons1 <- list( c("mutant","wildtype"))
ibdany<- ggboxplot(Totalpairsgreater0, x = "k13_Ala675Val_p1", y = "malecotf",
                   color = "k13_Ala675Val_p1") # Note the same code works fo =r 441l 574L 
ibdany+
  stat_compare_means()+
  #scale_y_continuous(trans = 'log2')+
  theme_bw () +
  labs(x="675V", y="IBD") +
  #scale_color_brewer(palette="Set2")+
  scale_color_manual(values = site)+
  theme(axis.text.x = element_text(angle =90))+ 
  stat_compare_means(comparisons=my_comparisons1 )+ 
  # Default method = "kruskal.test" for multiple groups # Global p-value
  ggtitle("IBD vs 675V")

# Save plot (change this for each mutations)
ggsave("FigS9A.svg", dpi=600, width=7.5, height=7)
ggsave("FigS9A.pdf", dpi=600, width=7.5, height=7)



#####################
#Figure S10 Haplotype Distribution and Prevalence Across Regions and Overall
###################

genotype_metadata <- read.csv("SPHAPLOTYPES.csv")

# Define colors for the plot
Bokre.site_col1 = c("#f032e6", "#4363d8", "#f58231", "#ffe119", "#911eb4", 
                             "#e6194b", "#fabebe", "darkred", "#46f0f0", "#3cb44b", "gray", "#008080")         
                             
# Data processing
HCcrt1 <- genotype_metadata %>% 
  filter(!is.na(IRNGE)) %>%
  #group_by(Region) %>% 
  summarise(total_ALT = sum(IRNGE>= 1, na.rm = TRUE),
            total_REF = sum(IRNGE == 0, na.rm = TRUE),
            N = total_ALT + total_REF,
            PREVALENCE = total_ALT / N * 100) %>% 
  mutate(ci.lower = PREVALENCE - 1.96 * sqrt(PREVALENCE * (100 - PREVALENCE) / N),
         ci.upper = PREVALENCE + 1.96 * sqrt(PREVALENCE * (100 - PREVALENCE) / N))


print(HCcrt1)
# Plotting
ggplot(data = HCcrt,
       aes(y = PREVALENCE,
           x = reorder(Region, Region),  # This will sort the Health_facility alphabetically
           fill = Region,
           ymin = ci.lower,
           ymax = ci.upper)) +
  theme(axis.text.x = element_text(angle = 90)) +  # Rotate x-axis labels for better visibility
  scale_fill_manual(values = Bokre.site_col1) +  # Apply the custom color palette
  geom_bar(stat = "identity") +  # Create the bars
  geom_errorbar(position = "dodge") +  # Add error bars
  ggtitle("IRNGE Prevalence per Region in Ethiopia")  # Title


write.csv(HCcrt, "IRNGEGRegion.csv", row.names = F)

###
#FigS11 A
###
## SP markers

library(miplicorn)
library(MIPanalyzer)
library(wesanderson)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(tidyverse)
library(dplyr)
library(tidyr)
library(viridis)
library(leaflet)
library(rhandsontable)
library(sp)
library(rgeos)

#FILES PATH
ref_file <- ("reference_AA_table.csv") 
alt_file <- ("alternate_AA_table.csv")
cov_file <- ("coverage_AA_table.csv")


relevant_mutations <- c("crt", "dhfr-ts", "dhps", "k13", "mdr1")
all_regions1<- read_tbl_ref_alt_cov(ref_file,
                                    alt_file,
                                    cov_file,
                                    gene == "dhfr-ts" | gene == "dhps"
)

prevalenceDHPS_DHFR_mutations<-mutation_prevalence(all_regions1, 5)

plot(prevalenceDHPS_DHFR_mutations, "Prevalence of SP resistance markers ")



###
#FigS11B
###

Bokre.site_col1 = c("#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                             "#e6194b",  "#fabebe", "darkred", "#46f0f0", "#3cb44b", "gray", "#008080")         
                             
# Region

genotype_metadata <-read.csv("genotype_metadata.csv")

HCcrt<- genotype_metadata %>% 
  filter(!is.na(crt_Lys76Thr)) %>%
  group_by(Health_facility,) %>% 
  summarise(total_ALT = sum(crt_Lys76Thr >= 1, na.rm = TRUE),
            total_REF = sum(crt_Lys76Thr == 0, na.rm = TRUE),
            N = total_ALT + total_REF,
            PREVALENCE = total_ALT / N*100) %>% 
  mutate(ci.lower = PREVALENCE - 1.96 * sqrt(PREVALENCE * (100 - PREVALENCE) / N),
         ci.upper = PREVALENCE + 1.96 * sqrt(PREVALENCE * (100 - PREVALENCE) / N))

#write.csv(HCcrt, "HC_crt_Lys76Thrpregion.csv", row.names = F)
ggplot(data = HCcrt,
       aes (
         y = PREVALENCE,
         x = reorder(Health_facility,PREVALENCE), fill=Health_facility,
         ymin = ci.lower,
         ymax = ci.upper)) +
  theme(axis.text.x = element_text(angle = 90))+
  
  scale_fill_manual(values=Bokre.site_col1)+
  geom_bar(stat = "identity") +
  geom_errorbar(position = "dodge") +
  ggtitle("76T  prevalence per HF Ethiopia")





mdr1haps <- read.csv("CRTHAP.csv")


# Load your haplotype data (assuming the CSV file is named 'haplotype_file.csv')
#haplotype_data <- read.csv("haplotype_file.csv")  # Modify the file path if necessary

# Calculate haplotype frequencies at the regional level (Region)
regional_haplotype_freq <- mdr1haps %>%
  group_by(Region, CTHAP) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100)


regional_haplotype_freq <- mdr1haps %>%
  group_by(Region, CTHAP) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Region) %>%  # Group by Region to calculate percentages within each region
  mutate(percentage = count / sum(count) * 100)  # Calculate percentage within each region


color_blind_friendly_paletteSP <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                                             "#0072B2", "#D55E00", "#CC79A7", "#999999")
                                             
# Create a pie chart for the overall haplotype distribution with confidence intervals
ggplot(overall_haplotype_freq, aes(x = "", y = percentage, fill = CTHAP)) +
  geom_bar(stat = "identity", width = 1) +  # Bar chart to create pie chart
  coord_polar(theta = "y") +  # Convert bar chart to pie chart
  theme_void() +  # Remove axis and background grid
  labs(title = "Overall Haplotype Distribution and Prevalence") +
  theme(legend.title = element_blank()) +  # Remove legend title
  # Text for the "Other" haplotype inside the pie chart
  geom_text(aes(label = ifelse(CTHAP == "Other", paste0(CTHAP, ": ", round(percentage, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), size = 5) +  # Inside pie text for "Other"
  # Text for all other haplotypes outside the pie chart
  geom_text(aes(label = ifelse(CTHAP != "Other", paste0(CTHAP, ": ", round(percentage, 1), "%"), "")), 
            size = 5, color = "black", 
            check_overlap = TRUE, nudge_y = 0.1) +  # Outside pie text for other haplotypes
  scale_fill_manual(values = color_blind_friendly_palette)  # Apply colorblind-friendly palette




################
################
# Figure S12A. Complexity of infections and population structure of P. falciparum populations. 

################
################

coiBokreEPHIAshe<- read.csv("FigS12.csv")


# COI over all 
ggplot(coiBokreEPHIAshe, aes(coimean)) + geom_histogram(binwidth = 0.5, fill = "grey30")+
  labs(y="Number of Samples", x="COI") 

coiBokreEPHIAshe %>% 
  group_by(coimean) %>% 
  tally()  %>% 
  arrange(n)


################
################
# FigS12B
################
################


IBC2COREL_vcf<- read.vcfR("borkesamples_biallelicsnp05q30sample05.recode.vcf", verbose = FALSE)

input='borkesamples_biallelicsnp05q30sample05.recode.vcf'


#) Please change the name of the prefix of output files if needed
name = "pf_ephi"


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


# plot pca as pdf 

figpdf = paste('PCA1_2BokreEPHIAshe_final.pdf', sep="")
pdf(file = figpdf)

plot(EV1, EV2,xlab="PC1(6.3%)", ylab="PC2(3.8)", col=site2[as.integer(pop1)],
     pch=19,cex=1.5)
legend("topright", legend=levels(pop1), bg="transparent",pch=16, cex=1,col=site2,text.col=site2)

abline(v=0.00, h=0.00, col="black", lwd=2, lty=2)

dev.off()

#plot pca site and mutation as pdf 

plot(EV1, EV2,xlab="PC1(6.3%)", ylab="PC2(3.8%)", col=site2[as.integer(pop1)],
     pch=shape.mutation[as.integer(pop2)],cex=1.2)
legend("topright", legend=levels(pop1), bg="transparent",pch=16, cex=1,col=site2,text.col=site2)
legend("topleft", legend=levels(pop2), bg="transparent",pch=shape.mutation, cex=1)
abline(v=0, h=0, col="black", lwd=1.5, lty=2)




###################
#Figure S13.  Heatmap showing pairwise IBD sharing within and between health facilities
###################

inf.m8<-read.csv("final_ibd_mle_long_EPHIBokreamples.csv", header=TRUE)
#metadcomplete



# Summarize the data by calculating the mean IBD for each combination of mutation_status and deme_p1
library(dplyr)


data_filtered <- inf.m8 %>%
  filter(IBD > 0, deme_p1 == deme_p2)

heatmap_data <- data_filtered %>%
  group_by(deme_p1, deme_p2) %>%
  summarise(mean_IBD = mean(IBD), .groups = 'drop')

# Create the heatmap

# Create a new variable to distinguish within-deme vs between-deme pairs
data_filtered <- data_filtered %>%
  mutate(pair_type = ifelse(deme_p1 == lag(deme_p1), "Within Deme", "Between Deme"))

# Summarize data by calculating the mean IBD for each combination of mutation_status, deme_p1, and pair_type
heatmap_data <- data_filtered %>%
  group_by(deme_p1, deme_p2, pair_type) %>%
  summarise(mean_IBD = mean(IBD), .groups = 'drop')

# Create the heatmap
ggplot(heatmap_data, aes(x = deme_p1, y = deme_p1, fill = mean_IBD)) +
  geom_tile() +
  facet_wrap(~ pair_type) +  # Facet by within or between deme pairs
  labs(
    title = "Heatmap of Mean IBD Sharing by Mutation Status and Pair Type (Within vs Between Deme)",
    x = "Mutation Status",
    y = "Health Centre (Deme)",
    fill = "Mean IBD"
  ) +
  scale_fill_gradient(low = "blue", high = "red") +  # Color gradient for mean IBD
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
    axis.text.y = element_text(angle = 0, hjust = 1)  # Ensure y-axis labels are readable
  )



inf.m7<-read.csv("final_ibd_mle_long_EPHIBokreamples.csv", header=TRUE)
set.seed(100)

dbs2 <- inf.m7

dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf >= 0.75)
relations <- dbs3 %>%
  dplyr::select(p1,p2,deme_p1)

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,deme_p1, deme_p2)


dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf >= 0.75)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,malecotf,deme_p1)

for.checks <- dbs3 %>%
  dplyr:: select(p1,p2,malecotf,deme_p1,deme_p2)

relations <- dbs3 %>%
  dplyr::select(p1,p2,malecotf,deme_p1)

colnames(relations) <- c("from","to","F","deme_p1")

relations <- relations[order(relations$deme_p1),]
relations$deme_p1 <- as.factor(relations$deme_p1)


dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf >= 0.75)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,malecotf,deme_p1)

colnames(relations) <- c("target","source","F","deme_p1")
relations <- relations[order(relations$deme_p1),]

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,malecotf,deme_p1,deme_p2)
colnames(for.checks) <- c("target","source","F","deme_p1","deme_p2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

links75 <- relations %>% dplyr::select(target, source, F)
write.csv(links75, "links99prec.csv" )

#links90<- read.csv("links95mod1prec.csv", header = T)
seq1 <- for.checks %>% dplyr::select(target,deme_p1) 
seq2 <- for.checks %>% dplyr::select(source,deme_p2) %>% rename(target=source,deme_p1=deme_p2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,deme_p1)
nodes <- nodes[order(nodes$deme_p1),]
nodes <- droplevels(nodes)

metadcomplete <- read.csv("final_reference_Bokre_metadatacoimod.csv")

ta75 <- metadcomplete %>% 
  dplyr::filter(metadcomplete$p1 %in% nodes$target) %>%
  dplyr::select(p1, deme_p1)


write.csv(ta75, "ta99msmtprec.csv") # modify and order per pop
#ta95<-read.csv("ta95msmtprec.csv", header = T)

#ta2 <- ta %>%
# select(seqID, hrp23_status_f)


IBDNW75 <- graph_from_data_frame(d=links75, vertices=ta75, directed=F) 


colors<-  c( "#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                      "#e6194b",  "#fabebe", "darkred",  "#46f0f0", "#3cb44b", "gray", "#008080")
                      


colors <- c(rep("#f032e6",19),   rep("#4363d8", 61), rep("#f58231", 13), rep("#ffe119", 22), rep( "#911eb4", 1),rep("#e6194b", 26),rep("#fabebe", 80), rep("darkred", 74), rep("#46f0f0", 125), rep("#3cb44b", 10), rep("gray", 57), rep("#008080", 35))

plot(IBDNW99, vertex.label=NA,vertex.size=7, main="IBD>=0.75", vertex.color=colors)

legend(x="topleft", legend=c("Abol","Andansa","Asendabo","Batu","Dila" ,"Erer", "Jiga", "Lante", "Metehara" , "Secha", "Shile", "Woreta"), col=c( "#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                                                                                                                                                           "#e6194b",  "#fabebe", "darkred",  "#46f0f0", "#3cb44b", "gray", "#008080"), cex=0.7, pch=c(19))
