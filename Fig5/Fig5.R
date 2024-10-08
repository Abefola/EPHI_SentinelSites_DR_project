

#################
#################
# 10-07-2024
# Written by Abebe Fola 
#################
#################

#################
#################
# This is the code to generate combined Fig5 
#################
#################

#####
# Fig5A
#####
ibd_mle_long<- read.csv("final_ibd_Bokresamples.csv")


#IBD distribution

mainplot <- ibd_mle_long %>%
  ggplot() +
  geom_histogram(aes(x=malecotf, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "grey45") +
  xlab("IBD") + ylab("frequency (%)") +
  theme_classic()

insetplot <- ibd_mle_long %>%
  ggplot() +
  geom_histogram(aes(x=malecotf, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "grey45") +
  xlab("IBD>0.5") + ylab("frequency (%)") +
  theme_classic() +
  coord_cartesian(xlim = c(0.5,1), ylim = c(0,5)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))
# cowplot
cowplot::ggdraw() +
  cowplot::draw_plot(mainplot, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(insetplot, x = 0.5, y= 0.3, width = 0.4, height = 0.4)


#####
# Fig5B
#####

colors =  c( "#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                      "#e6194b",  "#fabebe", "darkred",  "#46f0f0", "#3cb44b", "gray", "#008080")
                      

ibd_final <- read.csv("final_ibd_mle_long_EPHIBokreamples.csv") # modified to add deme_p1 (within vs between pop )


Totalpairsgreater0 <- ibd_final %>% # half sibling and above 
  filter (malecotf>0.0)

ggplot(Totalpairsgreater0) +
  aes(x = deme_p1, y = malecotf, color = deme_p1 ) +
  geom_boxplot(stat = "boxplot",
               position = "dodge2",
               outlier.stroke = 0.05) +
  theme(legend.position = "dodge2") +
  ylab("Pairwise IBD") +
  scale_color_manual(values = colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("IBD sharing (IBD>0.0) within Sites") 





#####
# Fig5C
# part 1
#......................
# read in map base
#####

library(tidyverse) # ggplot2, dplyr, tidyr, readr, purrr, tibble
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(scatterpie)
library(scales)
library (ggforce)
library (ggmap)
library (ggplot2)

library(ggridges)

# make plot

############
# load base map layers 
############

#load base map layers
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

coinot_filtered<-read.csv("final_ibd_mle_long_EPHIBokreamples.csv")

# GT edit: generate mock mean IBD values for each unique site
set.seed(123) # For reproducibility
unique_sites <- unique(coinot_filtered$deme_p1)
mock_mean_ibd <- data.frame(
  deme_p1 = unique_sites,
  mean_ibd = runif(length(unique_sites), min = 0, max = 1) # Random values between 0 and 1
)

withinpairs90 <- coinot_filtered %>%
  filter(IBD >= 0.90) %>%
  filter(deme_p1 != deme_p2) %>%
  left_join(mock_mean_ibd, by = "deme_p1")

ggplot() +
  geom_sf(data = sov110, color = 'black', size = 0.3, fill = ifelse(sov110$ADMIN == "Ethiopia", 'grey98', 'grey90')) +
  geom_sf(data = lakes8, color = "grey40", fill = "aliceblue", size = 0.5) +
  geom_sf(data = admin8, color = "grey80", size = 0.6) +
  geom_sf(data = admin110, color = "grey80", size = 0.4) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.65, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotate("text", x = 34.5, y = 13, label = "Sudan",
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 32.9, y = 9.5, label = "South Sudan", angle = 90,
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 39, y = 15, label = "Eritrea",
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 44.5, y = 4, label = "Somalia", angle = 45,
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 37.5, y = 3.3, label = "Kenya", angle = 45,
           color = "grey3", size = 5, fontface = "italic") +
  coord_sf(xlim = c(33.3, 46), ylim = c(1, 16.2), expand = T) +  
  geom_curve(data = withinpairs90,
             alpha = 0.5, size = 0.7,
             aes(x = longnum_p1, y = latnum_p1,
                 xend = longnum_p2, yend = latnum_p2,
                 color = IBD), show.legend = T) +
  geom_point(data = withinpairs90, aes(x = longnum_p1, y = latnum_p1, fill = mean_ibd),
             shape = 21, size = 4, color = "black") +
  scale_color_gradient("Pairwise IBD", low = "lightpink", high = "red") +
  scale_fill_gradient("Mean Site IBD", low = "lightblue", high = "darkblue") +
  theme_minimal()


# Add major motorways

# define the bounding box for Ethiopia
bbox <- c(xmin = 32.99, ymin = 3.4, xmax = 47.98, ymax = 14.85)

# get major highways from OpenStreetMap
ethiopia_highways <- opq(bbox = bbox) %>%
  add_osm_feature(key = 'highway', value = c('motorway', 'trunk', 'primary')) %>%
  osmdata_sf()

roads <- ethiopia_highways$osm_lines

ggplot() +
  geom_sf(data = sov110, color = 'black', size = 0.3, fill = ifelse(sov110$ADMIN == "Ethiopia", 'grey98', 'grey90')) +
  geom_sf(data = lakes8, color = "grey40", fill = "aliceblue", size = 0.5) +
  geom_sf(data = admin8, color = "grey80", size = 0.6) +
  geom_sf(data = admin110, color = "grey80", size = 0.4) +
  geom_sf(data = roads, color = "green", size = 0.5, linetype = "solid") +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.65, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotate("text", x = 34.5, y = 13, label = "Sudan",
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 32.9, y = 9.5, label = "South Sudan", angle = 90,
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 39, y = 15, label = "Eritrea",
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 44.5, y = 4, label = "Somalia", angle = 45,
           color = "grey3", size = 5, fontface = "italic") +
  annotate("text", x = 37.5, y = 3.3, label = "Kenya", angle = 45,
           color = "grey3", size = 5, fontface = "italic") +
  coord_sf(xlim = c(33.3, 46), ylim = c(1, 16.2), expand = TRUE) +  
  geom_curve(data = withinpairs90,
             alpha = 0.5, size = 0.7,
             aes(x = longnum_p1, y = latnum_p1,
                 xend = longnum_p2, yend = latnum_p2,
                 color = IBD), show.legend = TRUE) +
  geom_point(data = withinpairs90, aes(x = longnum_p1, y = latnum_p1, fill = mean_ibd),
             shape = 21, size = 4, color = "black") +
  scale_color_gradient("Pairwise IBD", low = "lightpink", high = "red") +
  scale_fill_gradient("Mean Site IBD", low = "lightblue", high = "darkblue") +
  theme_minimal()




##########

# Fig5D

##########
inf.m7<-read.csv("final_ibd_mle_long_EPHIBokreamples.csv", header=TRUE)
set.seed(100)

dbs2 <- inf.m7

dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf >= 0.99)
relations <- dbs3 %>%
  dplyr::select(p1,p2,deme_p1)

for.checks <- dbs3 %>%
  dplyr::select(p1,p2,deme_p1, deme_p2)


dbs3 <- dbs2 %>%
  dplyr:: filter(malecotf >= 0.99)
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
  dplyr:: filter(malecotf >= 0.99)
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

links99 <- relations %>% dplyr::select(target, source, F)
write.csv(links99, "links99prec.csv" )

#links90<- read.csv("links95mod1prec.csv", header = T)
seq1 <- for.checks %>% dplyr::select(target,deme_p1) 
seq2 <- for.checks %>% dplyr::select(source,deme_p2) %>% rename(target=source,deme_p1=deme_p2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,deme_p1)
nodes <- nodes[order(nodes$deme_p1),]
nodes <- droplevels(nodes)

metadcomplete <- read.csv("final_reference_Bokre_metadatacoimod.csv")

ta99 <- metadcomplete %>% 
  dplyr::filter(metadcomplete$p1 %in% nodes$target) %>%
  dplyr::select(p1, deme_p1)


write.csv(ta99, "ta99msmtprec.csv") # modify and order per pop
#ta95<-read.csv("ta95msmtprec.csv", header = T)

#ta2 <- ta %>%
# select(seqID, hrp23_status_f)


IBDNW99 <- graph_from_data_frame(d=links99, vertices=ta99, directed=F) 


colors<-  c( "#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                      "#e6194b",  "#fabebe", "darkred",  "#46f0f0", "#3cb44b", "gray", "#008080")
                      


colors <- c(rep("#f032e6",19),   rep("#4363d8", 61), rep("#f58231", 13), rep("#ffe119", 22), rep( "#911eb4", 1),rep("#e6194b", 26),rep("#fabebe", 80), rep("darkred", 74), rep("#46f0f0", 125), rep("#3cb44b", 10), rep("gray", 57), rep("#008080", 35))

plot(IBDNW99, vertex.label=NA,vertex.size=7, main="IBD>=0.99, n=500", vertex.color=colors)

legend(x="topleft", legend=c("Abol","Andansa","Asendabo","Batu","Dila" ,"Erer", "Jiga", "Lante", "Metehara" , "Secha", "Shile", "Woreta"), col=c( "#f032e6", "#4363d8", "#f58231", "#ffe119",   "#911eb4", 
                                                                                                                                                           "#e6194b",  "#fabebe", "darkred",  "#46f0f0", "#3cb44b", "gray", "#008080"), cex=0.7, pch=c(19))