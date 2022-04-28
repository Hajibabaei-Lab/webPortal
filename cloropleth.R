# Teresita M. Porter, April 26, 2022
# WWF watershed .shp files from Mike Wright
# Data cleaned by Artin Mashayekhi
# plot richness as cloropleth

# read in .shp file
library(sf)

# plotting & data flipping
library(leaflet)
library(htmltools)
library(dplyr)

###############################
# read in a .shp file
system.time(wwf_read <- st_read("Subwatersheds/WSC_subwatersheds.shp", quiet = TRUE))

# fix geometry data
wwf_wgs84 <- wwf_read %>%
  st_buffer(0) %>% # Make invalid geometries valid
  st_transform(crs = 4326) # Convert coordinates to WGS84

# Simplify so that map loads faster
simplified <- rmapshaper::ms_simplify(wwf_wgs84)


###############################
# Read in taxonomy file
tax <- read.csv("taxonomy.csv", header = TRUE, stringsAsFactors = FALSE)

# Rename Sample_Name to File_Name
tax$File_Name <- tax$Sample_Name

# Replace any dashes with underscores for separation
tax$File_Name <- gsub("-", "_", tax$File_Name)

# Split File_Number into different columns using the underscore as the separator
tax <- cbind(tax, data.frame(do.call('rbind', strsplit(as.character(tax$File_Name), '_', fixed = TRUE))))
names(tax)[34:43] <- c("Project", "Group", "Substrate", "Site", "Sample", "Replicate", "Date", "Marker", "IlluminaRun", "Amplicon")

# Put the field back together to match File_name in s
tax$File_Name <- paste(tax$Project, tax$Group, tax$Substrate, tax$Site, tax$Sample, tax$Replicate, tax$Date, tax$Marker, sep="-")

# Create taxon field for assignments
## For 200bp COI query, 95% confidence at species rank (sBP >= 0.70), 
## 99% confidence genus (gBP >= 0.30), family (fBP >= 0.20), no cutoff needed for order+
tax$taxon <- ifelse(tax$sBP >= 0.70, paste(tax$Phylum, tax$Class, tax$Species, sep=";"),
                    ifelse(tax$gBP >= 0.30, paste(tax$Phylum, tax$Class, tax$Genus, sep=";"),
                           ifelse(tax$fBP >= 0.20, paste(tax$Phylum, tax$Class, tax$Family, sep=";"),
                                  paste(tax$Phylum, tax$Class, tax$Order, sep=";"))))

# Give taxon and site per File_Name
tax.stats <- tax[,c("File_Name", "taxon", "Site")]


###############################
# Read in metadata
s <- read.csv("metadata.csv", header=TRUE, stringsAsFactors = FALSE)

# Get unique latitude, longitude, file name, and WSCSDA
s2 <- unique(s[,c("File_Name", "WSCSDA", "Lat","Long")])

# Drop row(s) if NA
s3 <- na.omit(s2)

# Map taxon, samples, and sites to s3 by File_Name
tax.meta <- merge(tax.stats, s3, by = "File_Name", all.x = TRUE)

# Fix this when the dataset entry is completed
tax.meta <- tax.meta[!tax.meta$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]

# Group by WSCSDA and unique taxon counts, sites, and File_Name samples
tax.richness <- data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(richness = n_distinct(taxon)))
tax.sites <- merge(data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(sites = n_distinct(Site))), tax.richness, all.x = TRUE)
tax.stats <- merge(data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(samples = n_distinct(File_Name))), tax.sites, all.x = TRUE)

# Add richness, samples, and sites to the .shp file
shp <- merge(simplified, tax.stats, by = "WSCSDA", all.x = TRUE)
# Set NAs to zero
shp[, 6:8][is.na(shp[, 6:8])] <- 0

# Create color palette
pal.tmp <- RColorBrewer::brewer.pal(3, "BuGn")
pal <- colorNumeric(palette = pal.tmp, domain = shp$richness)

# Basic choropleth with leaflet
m <- leaflet(simplified) %>% 
  addTiles()  %>% 
  setView( lat=10, lng=0 , zoom=2) %>%
  addPolygons(data=shp, color= "white", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~pal(richness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$richness, "taxa<br>", 
                            "Sites: ", shp$sites, "<br>", 
                            "Samples: ", shp$samples, "<br>"))
# Display chloropleth
m


#######################
# adapt to summarize major taxa such as Ephemeroptera, Plecoptera, Trichoptera, Odonata, Chironomidae, Other
# adapt to summarize number of taxa in each phylum
# ranks for these taxa can be found at ncbi taxonomy page, ex. https://www.ncbi.nlm.nih.gov/taxonomy/?term=Ephemroptera
# leaflet can be used to display this info as piecharts
# put together a preliminary report in R markdown to get feedback from Mehrdad
# https://rmarkdown.rstudio.com/