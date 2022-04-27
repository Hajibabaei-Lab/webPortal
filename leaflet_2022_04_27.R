# Teresita M. Porter, April 20, 2022
# WWF watershed .shp files from Mike Wright
# Data cleaned by Artin Mashayekhi

# read in .shp file
library(sf)
library(dplyr)

# plotting & data flipping
library(leaflet)
library(htmltools)
library(spatialEco)



###############################
# read in a .shp file
system.time(wwf_read <- st_read("Subwatersheds/WSC_subwatersheds.shp", quiet = TRUE))

# id should really be a key to map watersheds from WWF to our samples (done below, remove id step here)
wwf_wgs84 <- wwf_read %>%
  st_buffer(0) %>% # Make invalid geometries valid
  st_transform(crs = 4326) %>% # Convert coordinates to WGS84 already in WGS84
  mutate(id = c(1:length(rownames(wwf_read)))) # Add column with id to each site

# simplify so that map loads faster
simplified <- rmapshaper::ms_simplify(wwf_wgs84)



###############################
# read in metadata file containing sites data with site sample, latitude, and longitude information. 
s <- read.csv("metadata.csv", header=TRUE, stringsAsFactors = FALSE)

# get unique latitude, longitude, and label
s2 <- unique(s[,c("Watershed","Lat","Long")])

# drop rows that have NA values
s3 <- na.omit(s2) 

# plot sites using circles
leaflet(s3) %>%
  addTiles() %>%
  addPolygons(data=simplified, color = "darkgrey", fillColor="white",weight = 0.5, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              popup = simplified$WSCSDA) %>%
  # alternative plotting using clustered points
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
##omit the WSCSDA entry if your metadata doesn't have the values beforehand
s4 <- s[,c("Site_Sample", "Lat", "Long", "WSCSDA")]

# remove any rows with NAs or that are empty 
s5 <- na.omit(s4) 
s5 <- s5[!s5$Site_Sample=="",]

# turn into sf obj (combines lat and long into a single field)
s5 <- st_as_sf(s5, coords = c("Long", "Lat"), crs = 4326) 

# now find which points map to which polygons (watersheds)
##new_shape <- data.frame(point.in.poly(s5, simplified[,c("WSCSDA", "geometry")]))
  
# full map
pts_poly_map <- unique(s5[,c(1:2)])
##for mapping new_shape when you don't have the WSCSDA values beforehand
###pts_poly_map <- unique(new_shape[,c(1:2)])

# map with rows with NAs dropped
pts_poly_map2 <- na.omit(pts_poly_map)

# use this to find which Site_Samples can't be assigned properly, so they can be fixed manually
##unassigned3 <- pts_poly_map[is.na(pts_poly_map$WSCSDA),]
