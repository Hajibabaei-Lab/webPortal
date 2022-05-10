# Teresita M. Porter, May 3, 2002
# WWF watershed .shp files from Mike Wright
# Data cleaned by Artin Mashayekhi
# plot richness as cloropleth

###############################
# read in .shp file
library(sf)

# plotting & data flipping
library(leaflet)
library(leaflet.minicharts)
library(htmltools)
library(dplyr)
library(reshape2)


###############################
# read in a .shp file
wwf_read <- st_read("Subwatersheds/WSC_subwatersheds.shp", quiet = TRUE)

# fix geometry data
wwf_wgs84 <- wwf_read %>%
  st_buffer(0) %>% # Make invalid geometries valid
  st_transform(crs = 4326) # Convert coordinates to WGS84

# Simplify so that map loads faster
simplified <- rmapshaper::ms_simplify(wwf_wgs84)


###############################
# Read in taxonomy file
# STREAM_2018-2020_taxonomy.csv
tax <- read.csv("taxonomy.csv", header = TRUE, stringsAsFactors = FALSE)
# 275525     33

# Keep taxonomy File_Names as is, it's the metadata filenames that are mangled
# Rename Sample_Name to File_Name
tax$File_Name <- tax$SampleName

# Replace any dashes with underscores for separation
tax$File_Name <- gsub("-", "_", tax$File_Name)

# Split File_Number into different columns using the underscore as the separator
tax <- cbind(tax, data.frame(do.call('rbind', strsplit(as.character(tax$File_Name), '_', fixed = TRUE))))
names(tax)[34:43] <- c("Project", "Group", "Substrate", "Site", "Sample", "Replicate", "Date", "Marker", "IlluminaRun", "Amplicon")

# Put the field back together to just keep key fields
tax$File_Name <- paste(tax$Project, tax$Group, tax$Substrate, tax$Site, tax$Sample, tax$Replicate, tax$Date, tax$Marker, sep="-")

# Remove this sample with missing lat/lon info from taxonomy file right away, can't work with this sequence data
tax <- tax[!tax$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]
# 275206     43

# Create taxon field for assignments
## For 200bp COI query, 95% confidence at species rank (sBP >= 0.70), 
## 99% confidence genus (gBP >= 0.30), family (fBP >= 0.20), no cutoff needed for order+
tax$taxon <- ifelse(tax$sBP >= 0.70, paste(tax$Phylum, tax$Class, tax$Species, sep=";"),
                    ifelse(tax$gBP >= 0.30, paste(tax$Phylum, tax$Class, tax$Genus, sep=";"),
                           ifelse(tax$fBP >= 0.20, paste(tax$Phylum, tax$Class, tax$Family, sep=";"),
                                  paste(tax$Phylum, tax$Class, tax$Order, sep=";"))))


###############################
# Read in metadata
s <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# Get unique latitude, longitude, file name, and WSCSDA
s2 <- unique(s[,c("File_Name", "WSCSDA", "Lat","Long")])

# Drop row(s) if NA
s3 <- na.omit(s2)

# Get unique taxon and site per File_Name
tax.stats <- unique(tax[,c("File_Name", "taxon", "Site")])

# Map taxon, samples, and sites to s3 by File_Name
tax.meta <- merge(tax.stats, s3, by = "File_Name", all.x = TRUE)

# # File_Name in taxonomy file that don't match anything in metadata
# setdiff(tax.stats$File_Name, s3$File_Name)

# # Fix this when the dataset entry is completed
# tax.meta <- tax.meta[!tax.meta$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]

# Group by WSCSDA and unique taxon counts, sites, and File_Name samples
tax.richness <- data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(richness = n_distinct(taxon)))
tax.sites <- merge(data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(sites = n_distinct(Site))), tax.richness, by = "WSCSDA", all.x = TRUE)
tax.final <- merge(data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(samples = n_distinct(File_Name))), tax.sites, by = "WSCSDA", all.x = TRUE)

# Add richness, samples, and sites to the .shp file
shp <- merge(simplified, tax.final, by = "WSCSDA", all.x = TRUE)
# Set NAs to zero, gives warning but still works
shp[, 6:8][is.na(shp[, 6:8])] <- 0

# Create color palette
mybins <- c(0,1,500,1000,1500,Inf)
mypalette <- colorBin( palette="YlOrBr", domain=shp$richness, na.color="transparent", bins=mybins)
mylabels <- c("0", "<=500", "<=1000", "<=1500", ">1500")

# pal.tmp <- RColorBrewer::brewer.pal(5, "YlOrBr")
# pal <- colorNumeric(palette = pal.tmp, domain = shp$richness)

# Basic chloropleth with leaflet
m <- leaflet(simplified) %>% 
  addTiles()  %>% 
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data=shp, color= "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~mypalette(richness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$richness, "taxa<br>", 
                            "Sites: ", shp$sites, "<br>", 
                            "Samples: ", shp$samples, "<br>"))  %>%
  addLegend(data=shp, pal=mypalette, values=~richness, 
            opacity=1, title = "Richness", position = "bottomleft",
            labFormat = function(type, cuts, p) {paste0(mylabels)})
# Display chloropleth
m

# Simplified but with centroid coordinates, included watershed labels too for plotting
simplified.centroid <- st_centroid(simplified)
centroid.coords <- data.frame(st_coordinates(simplified.centroid))
simplified.centroid <- data.frame(simplified.centroid$WSCSDA, centroid.coords, simplified.centroid$WSCSDA_EN)
colnames(simplified.centroid) <- c("WSCSDA", "Long", "Lat", "WSCSDA_EN")


###############################
# Major Taxa
major.order <- c("Ephemeroptera", "Plecoptera_Insecta", "Trichoptera", "Odonata")
major.family <-c("Chironomidae")
# Get the File_Name, Order, and Family
tax.major.order <- tax[c("File_Name", "Order")]
tax.major.family <- tax[c("File_Name", "Family")]
# Get rid of unwanted orders and families
tax.major.order2 <- tax.major.order[grepl(paste(major.order, collapse = "|"), tax.major.order$Order),]
tax.major.family2 <- tax.major.family[grepl(paste(major.family, collapse = "|"), tax.major.family$Family),]

# Keep key fields, File_Name, Order, Family, taxon, site
tax.orderstats <- tax[, c("File_Name", "Order", "taxon", "Site")]
tax.familystats <- tax[, c("File_Name", "Family", "taxon", "Site")]

# Map taxon, samples, and sites to s3 by File_Name
tax.ordermeta <- merge(tax.orderstats, s3, by = "File_Name", all.x = TRUE)
tax.familymeta <- merge(tax.familystats, s3, by = "File_Name", all.x = TRUE)

# Fix this later
tax.ordermeta <- tax.ordermeta[!tax.ordermeta$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]
tax.familymeta <- tax.familymeta[!tax.familymeta$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]

# Group by WSCSDA and unique taxon counts, sites, and File_Name samples
tax.order <- data.frame(tax.ordermeta %>% group_by(WSCSDA, Order) %>% summarise(taxon = n_distinct(taxon)))
tax.family <- data.frame(tax.familymeta %>% group_by(WSCSDA, Family) %>% summarise(taxon = n_distinct(taxon)))

# Convert long to wide format and remove NAs
tax.order.wide <- dcast(tax.order, WSCSDA ~ Order, value.var = "taxon")
tax.order.wide[is.na(tax.order.wide)] <- 0
tax.family.wide <- dcast(tax.family, WSCSDA ~ Family, value.var = "taxon")
tax.family.wide[is.na(tax.family.wide)] <- 0

# Keep relevant columns
tax.order.wide2 <- tax.order.wide[, names(tax.order.wide) %in% major.order]
tax.order.wide2$WSCSDA <- tax.order.wide$WSCSDA
tax.family.wide2 <- tax.family.wide[, names(tax.family.wide) %in% major.family]

# Combine the 5 major taxa
tax.major.combo <- tax.order.wide2
tax.major.combo$Chironomidae <- tax.family.wide2

# Combine major taxa, WSCSDA, and coordinates
tax.major.final <- merge(tax.major.combo, simplified.centroid, by = "WSCSDA", all.x = TRUE)
names(tax.major.final)[names(tax.major.final) == 'Plecoptera_Insecta'] <- 'Plecoptera'

# # display the usual data in the popup
# my_popups <- simplified %>% 
#   group_by(WSCSDA_EN) %>% 
#     mutate(popup = paste0(WSCSDA_EN, "<br>")) %>% 
#   pull(popup)

my_popups <- paste("Watershed: ", tax.major.final$WSCSDA_EN, "<br>",
                   "Ephemeroptera: ", tax.major.final$Ephemeroptera, "<br>", 
                   "Odonata: ", tax.major.final$Odonata, "<br>", 
                   "Plecoptera: ", tax.major.final$Plecoptera, "<br>",
                   "Trichoptera: ", tax.major.final$Trichoptera, "<br>",
                   "Chironomidae: ", tax.major.final$Chironomidae, "<br>")

# Add color
pal <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')

# Generate leaflet map with pie charts
a <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data=simplified, color= "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = tax.major.final$Long, lat = tax.major.final$Lat, type = "pie",
                chartdata = tax.major.final[, c("Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera", "Chironomidae")],
                popup = popupArgs(html = my_popups))
# Display pie chart map
a


###############################
# Get the File_Name and Phylum
tax.ranks <- tax[c("File_Name", "Phylum")]
# Count taxa per phylum
phylum <- tax %>% count(tax$Phylum)
# Find the top 5 taxa for all the phyla
phylum <- head(phylum[order(-phylum$n),], 5)
topPhyla <- phylum$`tax$Phylum`

# Keep key fields, File_Name, Phylum, taxon, Site
tax.stats2 <- tax[, c("File_Name", "Phylum", "taxon", "Site")]

# Map taxon, samples, and sites to s3 by File_Name
tax.meta2 <- merge(tax.stats2, s3, by = "File_Name", all.x = TRUE)
# Fix this later
tax.meta2 <- tax.meta2[!tax.meta2$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]

# Group by WSCSDA and unique taxon counts, sites, and File_Name samples
tax.phylum <- data.frame(tax.meta2 %>% group_by(WSCSDA, Phylum) %>% summarise(taxon = n_distinct(taxon)))

# convert long to wide format
tax.phylum.wide <- dcast(tax.phylum, WSCSDA ~ Phylum, value.var = "taxon")
# set NA to 0
tax.phylum.wide[is.na(tax.phylum.wide)] <- 0

# only keep cols if they are in the target list (and WSCSDA key field for combining with shp file)
tax.phylum.wide2 <- tax.phylum.wide[,names(tax.phylum.wide) %in% topPhyla]
tax.phylum.wide2$WSCSDA <- tax.phylum.wide$WSCSDA

# Combine top 5 phyla, WSCSDA, and coordinates
tax.phylum.final <- merge(tax.phylum.wide2, simplified.centroid, by = "WSCSDA", all.x = TRUE)

#popups
my_popups2 <- paste("Watershed: ", tax.phylum.final$WSCSDA_EN, "<br>",
                   "Annelida: ", tax.phylum.final$Annelida, "<br>", 
                   "Arthropoda: ", tax.phylum.final$Arthropoda, "<br>", 
                   "Chordata: ", tax.phylum.final$Chordata, "<br>",
                   "Cnidaria: ", tax.phylum.final$Cnidaria, "<br>",
                   "Mollusca: ", tax.phylum.final$Mollusca, "<br>")

# Add color
pal <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')

# Generate leaflet map with pie charts
p <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data=simplified, color= "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = tax.phylum.final$Long, lat = tax.phylum.final$Lat, type = "pie", 
                chartdata = tax.phylum.final[, c("Annelida", "Arthropoda", "Chordata", "Cnidaria", "Mollusca")],
                popup = popupArgs(html = my_popups2))
# Display pie chart map
p


###############################
# Take user input on resolution
choice <- readline(prompt = "Select Resolution (GlobalESV, Genus, Family, Species, Phylum): ")

# Function for generating pie chart map
fx <- function(choice)  {
  # Find top 5 of resolution
  one <- tax %>% count(tax[choice])
  one <- head(one[order(-one$n), ], 5)
  top.five <- one[, 1]
  
  # Create data of relevant information by resolution and remove incomplete entry
  df <- tax[, c("File_Name", choice, "taxon", "Site")]
  df.meta <- merge(df, s3, by = "File_Name", all.x = TRUE)
  df.meta <- df.meta[!df.meta$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI", ]
  
  # Group by WSCSDA and resolution choice with File_Name and Sites. Also convert to Wide Format
  df.choice <- data.frame(df.meta %>% group_by(WSCSDA, df.meta[choice]) %>% summarise(taxon = n_distinct(taxon)))
  df.choice.wide <- dcast(df.choice, WSCSDA ~ df.choice[, choice], value.var = "taxon")
  df.choice.wide[is.na(df.choice.wide)] <- 0
  
  # Remove entries in resolution not in top 5
  df.choice.wide2 <- df.choice.wide[names(df.choice.wide) %in% top.five]
  df.choice.wide2$WSCSDA <- df.choice.wide$WSCSDA
  
  # Merge wide resolution with .shp data
  df.final <- merge(df.choice.wide2, simplified.centroid, by = "WSCSDA", all.x = TRUE)
  
  # Create map popup
  df.popups <- paste("Watershed: ", df.final$WSCSDA_EN, "<br>",
                      colnames(df.final[2]), ":", df.final[, 2], "<br>", 
                     colnames(df.final[3]), ":", df.final[, 3], "<br>", 
                     colnames(df.final[4]), ":", df.final[, 4], "<br>",
                     colnames(df.final[5]), ":", df.final[, 5], "<br>",
                     colnames(df.final[6]), ":", df.final[, 6], "<br>")

  # Add color
  pal <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')
  
  # Generate leaflet map with pie charts
  a <- leaflet(simplified) %>%
    addTiles() %>%
    # addProviderTiles(providers$Stamen.Terrain) %>%
    setView(lat = 60, lng = -95, zoom = 3) %>%
    addPolygons(data=simplified, color= "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
                fillColor = pal, fillOpacity = 0.5
    )  %>%
    addMinicharts(lng = df.final$Long, lat = df.final$Lat, type = "pie", 
                  chartdata = df.final[, 2:6],
                  popup = popupArgs(html = df.popups))
  return(a)
}

fx(choice)
