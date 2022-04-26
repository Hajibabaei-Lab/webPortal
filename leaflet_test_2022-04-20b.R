# Teresita M. Porter, April 20, 2022
# WWF watershed shapefiles from Mike Wright
# Data cleaned by Artin Mashayekhi

# read in shapefile
library(sf)
library(rgdal)
library(dplyr)

# plotting & data flipping
library(tidyverse)
library(leaflet)
library(geojsonio)
library(htmlwidgets)
library(htmltools)
library(spatialEco)

# automatically extract points in a polygon
library(raster)
library(sp)

###############################
# read in a shapefile
system.time(wwf_read <- st_read("Subwatersheds/WSC_subwatersheds.shp", 
                                quiet = TRUE))

# # a slightly slower way to read in a shapefile
# canada_raw <- readOGR(dsn = "data/gcd_000b11a_e", layer = "gcd_000b11a_e",
#                       use_iconv=TRUE, encoding="CP1250")

# check out the geometry (already WGS84)
wwf_geom <- st_geometry(wwf_read)

# switch off spherical geometry (may not be needed)
sf::sf_use_s2(FALSE)

# id should really be a key to map watersheds from WWF to our samples (done below, remove id step here)
wwf_wgs84 <- wwf_read %>%
  st_buffer(0) %>% # Make invalid geometries valid
  st_transform(crs = 4326) %>% # Convert coordinates to WGS84 already in WGS84
  mutate(id = c(1:length(rownames(wwf_read)))) # Add column with id to each site

# simplify so that map loads faster
simplified <- rmapshaper::ms_simplify(wwf_wgs84)



###############################
# read in site info
s <- read.csv("sites.csv", header=TRUE, stringsAsFactors = FALSE)

# get unique lat, long, and label (Watershed or Site.Sample)
s2 <- unique(s[,c("Watershed","Lat","Long")])

# convert coordinate to numeric
s2$Lat <- as.numeric(s2$Lat)
s2$Long <- as.numeric(s2$Long)

# drop row if NA, data is incomplete
s3 <- na.omit(s2) 

# palette, too many watersheds to get a distinct unique color for each
# pal <- colorFactor(
#   palette = 'Dark1',
#   domain = s3$Watershed
# )

# plot sites using circles (or with clustering for easier viewing)
leaflet(s3) %>%
  addTiles() %>%
  addPolygons(data=simplified, color = "darkgrey", fillColor="white",weight = 0.5, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              popup = simplified$WSCSDA) %>%
  # cluster points that are close to each other when zoomed out
  # click to zoom in
  # current labels are watershed but could be anything else
  #addCircleMarkers(~Long, ~Lat, popup=~htmlEscape(Watershed),
  #                  clusterOptions = markerClusterOptions()
  #                 ) 
  # add points as circles
  addCircles(~Long, ~Lat, popup=~htmlEscape(Watershed),
             color = "#0033CC") 



###############################
# get a unique Site_Sample -> WSCSDA map 
# automatically extract points in a polygon
# get unique coordinates for each Site_Sample
s4 <- s[,c("Site_Sample","Lat","Long")]

# remove any rows with NAs or that are empty 
s5 <- na.omit(s4) 
s5 <- s5[!s5$Site_Sample=="",]

# turn into sf obj (combines lat and long into a single field)
s5 <- st_as_sf(s5, coords = c("Long", "Lat"), crs = 4326) 

# now find which points map to which polygons (watersheds)
new_shape <- data.frame(point.in.poly(s5, simplified[,c("WSCSDA", "geometry")]))

# full map
pts_poly_map <- unique(new_shape[,c(1:2)])
# [1] 523   2

# map with rows with NAs dropped (fix this root problem in excel)
pts_poly_map2 <- na.omit(pts_poly_map) 
# [1] 503   2

# find which Site_Samples can't be assigned so this can be fixed manually
unassigned3 <- pts_poly_map[is.na(pts_poly_map$WSCSDA),]
