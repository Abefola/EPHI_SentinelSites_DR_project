#################
#################
# 10-07-2024
# Written by Abebe Fola 
#################
#################

#################
#################
# This is the code to generate combined Fig
#################
#################



#Set your directory 

setwd("C:/Users/afola/Desktop/Bioinformatics/Abefola_github/EPHI_SentinelSites_DR_project/Fig1") # change this to your working directory

# To check the current dir

getwd() 

#################
# load all required libraries 
#################

#devtools::install_github("ropensci/rnaturalearth")
#remove.packages (rnaturalearthdata, lib)

# get plotting data-------------------
samplessites<- read.csv("Fig1A_file.csv", header=TRUE, sep=",")

view(samplessites)

### column names: 
#[1] "Nr"                         "R_EA_H_V_BARCODE_NAME"      "barcode"                   
#[4] "R_EA_H_V_BARCODE"           "plate"                      "row"                       
#[7] "col"                        "gbs_library_id"             "HHID_Surveyb"              
#[10] "hc_order"                   "sort_order"                 "cluster_groups"            
#[13] "counts"                     "sample"                     "variety"                   
#[16] "cluster"                    "matched_released_varieties" "status"                    
#[19] "lat"                   "long"                  "VarCountsbyHHID"

### critical ones are latitude and logitude plus whatever category you want to color by.

#===================================================================#
#           load packages, download maps (GADM)                     #
#===================================================================#
# read in shapefile and .RData formats from Global Administrative Areas (gadm.org), 
library(ggplot2)
library(maptools)
library(rgdal)
library(ggmap) 
library(mapproj) 
library(raster) 
library(maps)
library(sp)
library(RColorBrewer)
library(tidyverse) # ggplot2, dplyr, tidyr, readr, purrr, tibble
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(scatterpie)
library(scales)
library (ggforce)
library (ggmap)
library(dplyr)
library (ggplot2)
library(ggpubr)

# global administrative district maps----------
adm1 <- readRDS("ETH_adm1.rds")

plot(adm1)

fadm1 = fortify(adm1)

class(fadm1)

head(fadm1)
dim(fadm1)

ggplot(data=fadm1,aes(x=long,y=lat,  border="gray10", fill=TRUE, bg="gray30", group=group)) + geom_path() #+theme_bw()
       
ggplot(data=fadm1,aes(x=long,y=lat,  border="gray10", fill=TRUE, bg="gray30", group=group)) + geom_path() + theme_bw()+ geom_point(x=samplessites$Longitude, y=samplessites$Latitude, col="orange",cex=1, pch=19)

ggplot(data=fadm1,aes(x=long,y=lat,  border="gray10", fill=TRUE, bg="gray30", group=group)) + geom_path() + geom_point(x=samplessites$Longitude, y=samplessites$Latitude)


#+theme_bw()

plot(adm1)

# Add a point on the map for each airport:
points(x=samplessites$Longitude, y=samplessites$Latitude, col="red",cex=1, pch=19)



# global administrative district maps----------
adm1 <- getData('GADM', country='ETH', level=0)
adm2 <- getData('GADM', country='ETH', level=1)
adm3 <- getData('GADM', country='ETH', level=2)

fadm1 = fortify(adm1)
fadm2 = fortify(adm2)
fadm3 = fortify(adm3)
view (fadm2)

#===================================================================#
#                       Prepare names of regions                    #
#===================================================================#
### Go to http://gadm.org/country, select Country: "ETH" and File format: "R (SpatialPolygonsDataFrame)"
### save in the working directory and load as below
GADM1 <- readRDS("ETH_adm1.rds")
ETH.adm1.spdf<-get("GADM1") # complex list

plot(GADM1)


### Get centroids of spatialPolygonDataFrame and convert to dataframe ----------- 
ETH.adm1.centroids.df <- data.frame(long=coordinates(ETH.adm1.spdf)[, 1], lat=coordinates(ETH.adm1.spdf)[, 2])

### Get names and id numbers corresponding to administrative areas-----------
ETH.adm1.centroids.df[, 'ID_1'] <- ETH.adm1.spdf@data[,'ID_1']
ETH.adm1.centroids.df[, 'NAME_1'] <- ETH.adm1.spdf@data[,'NAME_1']
save(ETH.adm1.centroids.df, file = "ETH.adm1.centroids.df.RData")


#===================================================================#
# Plot sample sites on ETH map, color by HC
#===================================================================#
 
ggplot(samplessites, aes(x=Longitude, y=Latitude)) +
 
  
  geom_polygon(data=fadm2, aes(x=long, y=lat, group = group), fill = "white", colour="grey0", size = 0.2) +
  geom_jitter(data=samplessites, position=position_jitter(width=0.2, height=0.2), size = 5, 
              # aes(x=lat, y=long),alpha = 6/10)+ #, color=variety, shape=factor(status)
              # aes(x=lat, y=long, alpha = 6/10, color=variety, shape=factor(status))) + ### if you want to color and shape
              aes(x=Longitude, y=Latitude, color=HC )) +
  # scale_color_brewer(palette = "",type = "qual") + 
 scale_color_manual(values = c( "#f032e6","#4363d8","#f58231","#ffe119", "#911eb4","#e6194b","#fabebe","tomato4","#46f0f0","#3cb44b","gray50","#008080")) +
  geom_text(data = ETH.adm1.centroids.df, aes(label = NAME_1, x = long, y = lat, group = NAME_1), size = 3)+
  coord_equal(ratio=1) +
  labs(title = "Sentinel site Health Centres", x = "Longitude", y = "Latitude") + 
   annotate("text", x = 34.5, y = 13, label = "Sudan",
            color="grey3", size=5 , fontface="italic") +
   annotate("text", x = 32.9, y = 9.5, label = "South Sudan", angle= 90,
            color="grey3", size=5 , fontface="italic") +
   annotate("text", x = 39, y = 15, label = "Eritrea",
            color="grey3", size=5 , fontface="italic") +
   annotate("text", x = 44.5, y=4, label = "Somalia", angle= 45,
            color="grey3", size=5 , fontface="italic") +
   annotate("text", x = 37.5, y = 3.3, label = "Kenya", angle= 45,
            color="grey3", size=5 , fontface="italic")+
   
  theme_grey(base_size = 10)+
coord_sf(xlim = c(29.3, 45.5), ylim = c(3,14.7), expand = T)



###################
#Africa map - fig1A_insert
###################

#readings
#https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/

############
# load base map layers 
############

admin8 <- ne_download(scale="large", type = "admin_1_states_provinces_lines",
                      category = "cultural", returnclass = "sf")
rivers8 <- ne_download(scale = 10, type = 'rivers_lake_centerlines',
                       category = 'physical', returnclass = "sf")
lakes8 <- ne_download(scale = "large", type = 'lakes',
                      category = 'physical', returnclass = "sf")
sov10 <- ne_download(scale="large", type = "sovereignty",
                     category = "cultural", returnclass = "sf")
admin10<- ne_download(scale="large", type = "populated_places",
                      category = "cultural", returnclass = "sf")



#ETH_regions <- dplyr::right_join(ETH_regions,keydrugresistance_prev, by = "NAME_1" )

ETHO<-st_read("gadm36_ETH_0.shp") # Source https://www.diva-gis.org/datadown
ETHO_regions<-st_read("gadm36_ETH_1.shp")
View(ETHO_regions)
ETHO_districts<-st_read("gadm36_ETH_2.shp") # name NAME_2 for combining
#View(ETHO_districts)
#Set coordinate reference system for ETHOia
ETHO_crs<-st_crs(ETHO)
africa <- ne_countries(scale="medium", type = "sovereignty", continent = "Africa", returnclass = "sf", )




ggplot() +
  geom_sf(data=africa, fill="gray90")+
  geom_sf(data = ETHO, fill="white", lwd=0.5) +
  geom_sf(data = ETHO_regions, fill= "white",) +
  theme_bw()

#Fig1B

# Number of samples in each Region:

dbs_metadata<- read.csv("Fig1B_file.csv.csv", header=T,sep=",")

dbs_metadata %>% 
  group_by(Region) %>% 
  tally()  %>% 
  arrange(n)
#Bokre.site_col = c( "#e6194b", "tomato4",  "#f58231", "#911eb4",    
                               #"#f032e6",  "#fabebe", "#ffe119","#46f0f0","#4363d8","#008080", "#3cb44b", "gray")                 
                            
Bokre.site_col = c( "#f032e6","#4363d8","#f58231","#ffe119", "#911eb4","#e6194b","#fabebe","tomato4","#46f0f0","#3cb44b","gray","#008080")



## count number of  samples at HC

dbs_metadata %>% 
  group_by(Health_facility) %>% 
  tally()  %>% 
  arrange(n)

barplotHC <-  ggplot(dbs_metadata, aes(x = reorder(Health_facility,Region), fill= Health_facility))
barplotHC + 
  geom_bar(aes(fill= Health_facility)) +
 theme_bw() +
  labs(x="Health_facility", y="Number of samples") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill=guide_legend(title="Health Centre"))+
  scale_fill_manual(values=Bokre.site_col)+
  #scale_fill_jco()+
  #theme() +
  coord_flip() +
  ggtitle("Samples Size Distribution per Health centre")










