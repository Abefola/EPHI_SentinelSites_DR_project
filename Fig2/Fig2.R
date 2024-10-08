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

setwd("C:/Users/afola/Desktop/Bioinformatics/Abefola_github/EPHI_SentinelSites_DR_project/Fig2") # change this to your working directory

#################
# load all required libraries 
#################

#devtools::install_github("ropensci/rnaturalearth")
#remove.packages (rnaturalearthdata, lib)


#Install Library

library(ggspatial)
library(scatterpie)
library(scales)
library (ggforce)
library (ggmap)
# run the next line if you already have rstan installed
# remove.packages(c("StanHeaders", "rstan"))

#install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(rstan)
library(gt)
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
library(sf)
library(rgeos)
library(ggpubr)
library(scales)
library(gridExtra)
library(scatterplot3d)
library(knitr)
library(ggExtra)
library(GGally)
library(viridis)
library("ggsci")
library(rnaturalearth)
library(rnaturalearthdata)
library(UpSetR)

#readings
#https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/

############
# load base map layers 
############


# new approach

df_freq <- read.table("allefre_Ethiopia_2019_2023.txt")

admin8 <- ne_download(scale="large", type = "admin_1_states_provinces_lines",
                      category = "cultural", returnclass = "sf")
rivers8 <- ne_download(scale = 10, type = 'rivers_lake_centerlines',
                       category = 'physical', returnclass = "sf")
lakes8 <- ne_download(scale = "large", type = 'lakes',
                      category = 'physical', returnclass = "sf")
sov110 <- ne_download(scale="medium", type = "sovereignty",
                      category = "cultural", returnclass = "sf")
admin110 <- ne_download(scale="large", type = "populated_places",
                        category = "cultural", returnclass = "sf")
#AddisAbaba <- data.frame(name="Addis Ababa", long=32.5833, lat=0.31666)


df_freq_comb <- read.csv("allefre_Ethiopia_2019_2023mod1.csv")

#view(df_freq_comb)

pie_data <- df_freq_comb %>% 

  dplyr:: select(District, year, Locus, freq, lat, long)

pie_data <- pie_data %>%
  pivot_wider(values_from = freq, names_from = Locus) %>%
  #mutate(WT = 1-PfK13) %>%
  mutate(radius = Wildtype) %>%
  filter(radius >=0)

#rescale radius
minr<-min(pie_data$radius)  # min in the data original scale
maxr<-max(pie_data$radius)  # max in the data original scale
min_scale<-0.1  ## min in the new scale
max_scale<-0.3 ## max in the new scale
pie_data$radius_adj <- min_scale + (max_scale-min_scale)*pie_data$radius/(maxr-minr)

pie_chart_df <- pie_data %>% filter(year == 2022)

#View(pie_chart_df)
ggplot() +
  geom_sf(data=sov110, color='black', size=0.8, fill = ifelse(sov110$ADMIN == "Ethiopia", 'grey98', 'grey90')) +
  # geom_sf(data=rivers8, color="cyan4", size=0.5, alpha=0.5) +
  geom_sf(data=lakes8, color="grey30", fill ="aliceblue", size= 0.5) +
  geom_sf(data=admin8, color="grey80", size= 0.6) +
  geom_sf(data=admin110, color="grey80", size= 0.4) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(which_north = "true",  height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         pad_x = unit(0.85, "cm"),
                         pad_y = unit(0.6, "cm"))+
  annotate("text", x = 34.5, y = 13, label = "Sudan",
           color="grey3", size=5 , fontface="italic") +
  annotate("text", x = 32.5, y = 9.5, label = "South Sudan", angle= 90,
           color="grey3", size=5 , fontface="italic") +
  annotate("text", x = 39, y = 15, label = "Eritrea",
           color="grey3", size=5 , fontface="italic") +
  annotate("text", x = 44.5, y=4, label = "Somalia", angle= 90,
           color="grey3", size=5 , fontface="italic") +
  annotate("text", x = 37.5, y = 3.5, label = "Kenya", angle= 90,
           color="grey3", size=5 , fontface="italic") +
  
  coord_sf(xlim = c(31, 46), ylim = c(3,15), expand = T) +
  theme_void() +
  #coord_sf(xlim = c(32, 40), ylim = c(7.5,15), expand = T)
  #coord_sf(xlim = c(29.5, 40), ylim = c(-1,-12), expand = TRUE) + # orginal Tanzi data
  
  geom_scatterpie(data = pie_chart_df, 
                  aes(x=long, y=lat, group=District, r = radius_adj), 
                  cols=c("R622I", "A675V", "P441L", "P574L", "Wildtype")) +
  geom_scatterpie_legend(pie_chart_df$radius_adj,
                         x=28.7,
                         y=0.7,
                         labeller = function(x) round((x*0.5441068-0.1)/0.35 *100),
                         #labeller = function(x) round(((0.639899*(x-0.1))/0.35)*100), #convert adjusted radius size back to the mutant/genotyped
                         n=4)+
  theme(legend.text=element_text(size=10), 
        axis.text = element_text(size = 10), 
        legend.position = c(0.085,0.84)) +
  scale_fill_manual(values = c("#fabebe", "#4363d8",  "#46f0f0","#008080", "white"),
                    labels = c("R622I", "A675V", "P441L", "P574L", "Wildtype"),
                    name = "Mutations")


# Save plot 
ggsave("Fig2A.svg", dpi=600, width=7.5, height=7)
ggsave("Fig2A.pdf", dpi=600, width=7.5, height=7)




################

# Fig2B

#############

####
#PLOT > Zero per site
####

site = c( "#CC79A7", "#008080")

ibd_mle_long_mtdt_622I <- read.csv("final_ibd_mle_long_622I.csv")

Totalpairsgreater0 <- ibd_mle_long_mtdt_622I %>%
  filter (malecotf>0.0)



my_comparisons1 <- list( c("mutant","wildtype"))
ibdany<- ggboxplot(Totalpairsgreater0, x = "k13_Arg622Ile_p1", y = "malecotf",
                   color = "k13_Arg622Ile_p1")
ibdany+
  stat_compare_means()+
  #scale_y_continuous(trans = 'log2')+
  theme_bw () +
  labs(x="622I", y="IBD") +
  #scale_color_brewer(palette="Set2")+
  scale_color_manual(values = site)+
  theme(axis.text.x = element_text(angle =90))+ 
  stat_compare_means(comparisons=my_comparisons1 )+ 
  # Default method = "kruskal.test" for multiple groups # Global p-value
  ggtitle("IBD vs 622I")

# Save plot 
ggsave("Fig2B.svg", dpi=600, width=7.5, height=7)
ggsave("Fig2B.pdf", dpi=600, width=7.5, height=7)




###############

# Fig2C network IBD sharing mutant vs wildtype

###############


inf.IBD<-read.csv("final_ibd_mle_long_622I.csv", header=TRUE)
set.seed(100)

dbs2 <- inf.IBD

dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf > 0.99)
relations <- dbs3 %>%
  dplyr::select(p1,p2,k13_Arg622Ile_p1)

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,k13_Arg622Ile_p1, k13_Arg622Ile_p2)


dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf > 0.99)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,malecotf,k13_Arg622Ile_p1)

for.checks <- dbs3 %>%
  dplyr:: select(p1,p2,malecotf,k13_Arg622Ile_p1,k13_Arg622Ile_p2)

relations <- dbs3 %>%
  dplyr::select(p1,p2,malecotf,k13_Arg622Ile_p1)

colnames(relations) <- c("from","to","F","k13_Arg622Ile_p1")

relations <- relations[order(relations$k13_Arg622Ile_p1),]
relations$k13_Arg622Ile_p1 <- as.factor(relations$k13_Arg622Ile_p1)


dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf > 0.99)
relations <- dbs3 %>%
  
  dplyr::select(p1,p2,malecotf,k13_Arg622Ile_p1)

colnames(relations) <- c("target","source","F","k13_Arg622Ile_p1")
relations <- relations[order(relations$k13_Arg622Ile_p1),]

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,malecotf,k13_Arg622Ile_p1,k13_Arg622Ile_p2)
colnames(for.checks) <- c("target","source","F","k13_Arg622Ile_p1","k13_Arg622Ile_p2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

links99 <- relations %>% dplyr::select(target, source, F)
write.csv(links99, "links622IIBD1prec.csv" )

#links90<- read.csv("links95mod1prec.csv", header = T)
seq1 <- for.checks %>% dplyr::select(target,k13_Arg622Ile_p1) 
seq2 <- for.checks %>% dplyr::select(source,k13_Arg622Ile_p2) %>% rename(target=source,k13_Arg622Ile_p1=k13_Arg622Ile_p2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,k13_Arg622Ile_p1)
nodes <- nodes[order(nodes$k13_Arg622Ile_p1),]
nodes <- droplevels(nodes)

metadcomplete <- read.csv("final_Bokre_k13_metadata.csv")

ta99 <- metadcomplete %>% 
  dplyr::filter(metadcomplete$p1 %in% nodes$target) %>%
  dplyr::select(p1, k13_Arg622Ile_p1)


write.csv(ta99, "ta99622IIBDprec.csv") # modify and order per pop
#ta99<-read.csv("ta99622Iprec.csv", header = T)

#ta2 <- ta %>%
# select(seqID, hrp23_status_f)


IBDNW99 <- graph_from_data_frame(d=links99, vertices=ta99, directed=F) 


#colors<-  c( "#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                      #"#e6194b",  "#fabebe", "darkred",  "#46f0f0", "#3cb44b", "gray", "#008080")
                      


colors <- c(rep("#008080",410),   rep("#CC79A7", 83))

plot(IBDNW99, vertex.label=NA,vertex.size=7, main="IBD>=0.99, n=493", vertex.color=colors)

legend(x="topleft", legend=c("Wildtype","Mutant"), col=c( "#008080", "#CC79A7"), cex=0.7, pch=c(19))

# Save plot 
ggsave("Fig2C.svg", dpi=600, width=7.5, height=7)
ggsave("Fig2C.pdf", dpi=600, width=7.5, height=7)





