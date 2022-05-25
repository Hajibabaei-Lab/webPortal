# Teresita M. Porter, May 10, 2022
# WWF watershed .shp files from Mike Wright
# Data cleaned by Artin Mashayekhi
# plot richness as chloropleth

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
tax <- read.csv("taxonomy.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
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
tax$File_Name <- paste(tax$Project, tax$Group, tax$Substrate, tax$Site, tax$Sample, tax$Replicate, tax$Date, tax$Marker, sep = "-")

# Remove this sample with missing lat/lon info from taxonomy file right away, can't work with this sequence data
tax <- tax[!tax$File_Name == "STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]
# 275206     43

# Create taxon field for assignments
## For 200bp COI query, 95% confidence at species rank (sBP >= 0.70), 
## 99% confidence genus (gBP >= 0.30), family (fBP >= 0.20), no cutoff needed for order+
tax$taxon <- ifelse(tax$sBP >= 0.70, paste(tax$Phylum, tax$Class, tax$Species, sep = ";"),
                    ifelse(tax$gBP >= 0.30, paste(tax$Phylum, tax$Class, tax$Genus, sep = ";"),
                           ifelse(tax$fBP >= 0.20, paste(tax$Phylum, tax$Class, tax$Family, sep = ";"),
                                  paste(tax$Phylum, tax$Class, tax$Order, sep = ";"))))


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

# Group by WSCSDA and unique taxon counts, sites, and File_Name samples
tax.richness <- data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(richness = n_distinct(taxon)))
tax.sites <- merge(data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(sites = n_distinct(Site))), tax.richness, by = "WSCSDA", all.x = TRUE)
tax.final <- merge(data.frame(tax.meta %>% group_by(WSCSDA) %>% summarise(samples = n_distinct(File_Name))), tax.sites, by = "WSCSDA", all.x = TRUE)

# Add richness, samples, and sites to the .shp file
shp <- merge(simplified, tax.final, by = "WSCSDA", all.x = TRUE)
names(shp)[8] <- "TaxonRichness"

# calc esv.richness
# Get unique taxon and site per File_Name
esv.stats <- unique(tax[, c("File_Name", "GlobalESV", "Site")])
# Map taxon, samples, and sites to s3 by File_Name
esv.meta <- merge(esv.stats, s3, by = "File_Name", all.x = TRUE)
# calc richness
esv.richness <- data.frame(esv.meta %>% group_by(WSCSDA) %>% summarise(richness = n_distinct(GlobalESV)))

# Add richness, samples, and sites to the .shp file
shp <- merge(shp, esv.richness, by = "WSCSDA", all.x = TRUE)
names(shp)[9] <- "ESVrichness"

# calc species.richness see https://github.com/terrimporter/CO1Classifier for cutoffs
# filter to only keep species id's, 95% confidence, 200bp, COI, sBP >= 0.70
tax.sp <- tax[tax$sBP >= 0.70,]
# Get unique taxon and site per File_Name
species.stats <- unique(tax.sp[, c("File_Name", "Species", "Site")])
# Map taxon, samples, and sites to s3 by File_Name
species.meta <- merge(species.stats, s3, by = "File_Name", all.x = TRUE)
# calc richness
species.richness <- data.frame(species.meta %>% group_by(WSCSDA) %>% summarise(richness = n_distinct(Species)))

# Add richness, samples, and sites to the .shp file
shp <- merge(shp, species.richness, by = "WSCSDA", all.x = TRUE)
names(shp)[10] <- "SpeciesRichness"

# calc genus.richness
# filter to only keep species id's, 99% confidence, 200bp, COI, gBP >= 0.30
tax.g <- tax[tax$gBP >= 0.30,]
# Get unique taxon and site per File_Name
genus.stats <- unique(tax.g[, c("File_Name", "Genus", "Site")])
# Map taxon, samples, and sites to s3 by File_Name
genus.meta <- merge(genus.stats, s3, by = "File_Name", all.x = TRUE)
# calc richness
genus.richness <- data.frame(genus.meta %>% group_by(WSCSDA) %>% summarise(richness = n_distinct(Genus)))

# Add richness, samples, and sites to the .shp file
shp <- merge(shp, genus.richness, by = "WSCSDA", all.x = TRUE)
names(shp)[11] <- "GenusRichness"

# calc family.richness
# filter to only keep species id's, 99% confidence, 200bp, COI, fBP >= 0.20
tax.f <- tax[tax$fBP >= 0.20,]
# Get unique taxon and site per File_Name
family.stats <- unique(tax.f[, c("File_Name", "Family", "Site")])
# Map taxon, samples, and sites to s3 by File_Name
family.meta <- merge(family.stats, s3, by = "File_Name", all.x = TRUE)
# calc richness
family.richness <- data.frame(family.meta %>% group_by(WSCSDA) %>% summarise(richness = n_distinct(Family)))

# Add richness, samples, and sites to the .shp file
shp <- merge(shp, family.richness, by = "WSCSDA", all.x = TRUE)
names(shp)[12] <- "FamilyRichness"

# Set NAs to zero, gives warning but still works
shp[, 6:12][is.na(shp[, 6:12])] <- 0

###############################
# create plots for each tab in Rmd

# ESV tab
# Create color palette
#range(shp$ESVrichness)
var_cut <- cut(shp$ESVrichness, breaks = 5)
#levels(var_cut)
# [1] "(-59.4,1.19e+04]"    "(1.19e+04,2.38e+04]" "(2.38e+04,3.57e+04]"
# [4] "(3.57e+04,4.75e+04]" "(4.75e+04,5.95e+04]"
mybins <- c(0, 1, 11900, 23800, 35700, 47500, 59500)
mypalette <- colorBin(palette = "YlOrBr", domain = shp$ESVrichness, na.color = "transparent", bins = mybins)
mylabels <- c("0", "<=11,900", "<=23,800", "<=35,700", "<=47,500", "<=59,500")

# ESV richness
m.esv <- leaflet(simplified) %>% 
  addTiles()  %>% 
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = shp, color = "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~mypalette(ESVrichness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$ESVrichness, "ESVs<br>", 
                            "Sites: ", shp$sites, "<br>", 
                            "Samples: ", shp$samples, "<br>"))  %>%
  addLegend(data = shp, pal = mypalette, values = ~ESVrichness, 
            opacity = 1, title = "ESV Richness", position = "bottomleft",
            labFormat = function(type, cuts, p){paste0(mylabels)})
# Display chloropleth
m.esv

# Species tab
# Create color palette
#range(shp$SpeciesRichness)
var_cut <- cut(shp$SpeciesRichness, breaks = 5)
#levels(var_cut)
# [1] "(-0.782,156]" "(156,313]"    "(313,469]"    "(469,626]"   
# [5] "(626,783]" 
mybins <- c(0, 1, 156, 313, 469, 626, 783)
mypalette <- colorBin(palette = "YlOrBr", domain = shp$SpeciesRichness, na.color = "transparent", bins = mybins)
mylabels <- c("0", "<=156", "<=313", "<=469", "<=626", "<=783")

# species richness
m.species <- leaflet(simplified) %>% 
  addTiles()  %>% 
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = shp, color = "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~mypalette(SpeciesRichness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$SpeciesRichness, "species<br>", 
                            "Sites: ", shp$sites, "<br>", 
                            "Samples: ", shp$samples, "<br>"))  %>%
  addLegend(data = shp, pal = mypalette, values = ~SpeciesRichness, 
            opacity = 1, title = "Species Richness", position = "bottomleft",
            labFormat = function(type, cuts, p){paste0(mylabels)})
# Display chloropleth
m.species

# Genus tab
# Create color palette
#range(shp$GenusRichness)
var_cut <- cut(shp$GenusRichness, breaks = 5)
#levels(var_cut)
# [1] "(-0.681,136]" "(136,272]"    "(272,409]"    "(409,545]"   
# [5] "(545,682]" 
mybins <- c(0, 1, 136, 272, 409, 545, 682)
mypalette <- colorBin(palette = "YlOrBr", domain = shp$GenusRichness, na.color = "transparent", bins = mybins)
mylabels <- c("0", "<=136", "<=272", "<=409", "<=545", "<=682")

# genus richness
m.genus <- leaflet(simplified) %>% 
  addTiles()  %>% 
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = shp, color = "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~mypalette(GenusRichness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$GenusRichness, "genera<br>", 
                            "Sites: ", shp$sites, "<br>", 
                            "Samples: ", shp$samples, "<br>"))  %>%
  addLegend(data = shp, pal = mypalette, values = ~GenusRichness, 
            opacity = 1, title = "Genus Richness", position = "bottomleft",
            labFormat = function(type, cuts, p){paste0(mylabels)})
# Display chloropleth
m.genus

# Family tab
# Create color palette
#range(shp$FamilyRichness)
var_cut <- cut(shp$FamilyRichness, breaks = 5)
#levels(var_cut)
# [1] "(-0.375,75]" "(75,150]"    "(150,225]"   "(225,300]"   "(300,375]" 
mybins <- c(0, 1, 75, 150, 225, 300, 375)
mypalette <- colorBin(palette = "YlOrBr", domain = shp$FamilyRichness, na.color = "transparent", bins = mybins)
mylabels <- c("0", "<=75", "<=150", "<=225", "<=300", "<=375")

# family richness
m.family <- leaflet(simplified) %>% 
  addTiles()  %>% 
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = shp, color = "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~mypalette(FamilyRichness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$FamilyRichness, "families<br>", 
                            "Sites: ", shp$sites, "<br>", 
                            "Samples: ", shp$samples, "<br>"))  %>%
  addLegend(data = shp, pal = mypalette, values = ~FamilyRichness, 
            opacity = 1, title = "Family Richness", position = "bottomleft",
            labFormat = function(type, cuts, p){paste0(mylabels)})
# Display chloropleth
m.family

# Taxon tab (finest taxonomic level of resolution with good confidence)
# Create color palette
#range(shp$TaxonRichness)
var_cut <- cut(shp$TaxonRichness, breaks = 5)
#levels(var_cut)
# [1] "(-1.6,320]"         "(320,640]"          "(640,960]"         
# [4] "(960,1.28e+03]"     "(1.28e+03,1.6e+03]"
mybins <- c(0, 1, 320, 640, 960, 1280, 1600)
mypalette <- colorBin(palette = "YlOrBr", domain = shp$TaxonRichness, na.color = "transparent", bins = mybins)
mylabels <- c("0", "<=320", "<=640", "<=960", "<=1280", "<=1600")

# family richness
m.taxon <- leaflet(simplified) %>% 
  addTiles()  %>% 
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = shp, color = "darkgrey", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~mypalette(TaxonRichness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$TaxonRichness, "taxa<br>", 
                            "Sites: ", shp$sites, "<br>", 
                            "Samples: ", shp$samples, "<br>"))  %>%
  addLegend(data = shp, pal = mypalette, values = ~TaxonRichness, 
            opacity = 1, title = "Taxon Richness", position = "bottomleft",
            labFormat = function(type, cuts, p){paste0(mylabels)})
# Display chloropleth
m.taxon

# Simplified but with centroid coordinates, included watershed labels too for plotting
simplified.centroid <- st_centroid(simplified)
centroid.coords <- data.frame(st_coordinates(simplified.centroid))
simplified.centroid <- data.frame(simplified.centroid$WSCSDA, centroid.coords, simplified.centroid$WSCSDA_EN)
colnames(simplified.centroid) <- c("WSCSDA", "Long", "Lat", "WSCSDA_EN")


###############################
# Major Taxa
target_orders <- c("Ephemeroptera", "Plecoptera_Insecta", "Trichoptera", "Odonata")
target_family <-c("Chironomidae")

# filter by target taxa
# COI, 200bp, 99% confidence for any oBP >= 0
tax.major.orders <- tax[tax$Order %in% target_orders,]
# COI, 200bp, 99% confidence if fBP >= 0.20
tax.major.family <- tax[tax$fBP >= 0.20 & 
                          tax$Family %in% target_family,]

# put together filtered records
tax.major <- rbind(tax.major.orders, tax.major.family)

# simplify Plecoptera_Insecta
tax.major$Order <- gsub("Plecoptera_Insecta", "Plecoptera", tax.major$Order)

# function to create df for plotting minicharts
make_minichart_df <- function(df, rank){
  # df options: tax.major
  # rank options: GlobalESV, Species, Genus, Family, taxon
  # cutoff options: NA, species sBP >= 0.70, genus gBP >= 0.30, species sBP >= 0.20, NA
  
  # Keep key fields, File_Name, Order, Family, rank, site
  # filter by bootstrap values where needed
  if (rank == "GlobalESV") {
    df <- df[, c("File_Name","Order","Family","GlobalESV","Site")]
  } else if (rank == "Species") {
    df <- df[df$sBP >= 0.70,]
    df <- df[, c("File_Name","Order","Family","Species","Site")]
  } else if (rank == "Genus") {
    df <- df[df$gBP >= 0.30,]
    df <- df[, c("File_Name","Order","Family","Genus","Site")]
  } else if (rank == "Family") {
    df <- df[df$fBP >= 0.20,]
    df <- df[, c("File_Name","Order","Family","Family","Site")]
  } else {
    df <- df[, c("File_Name","Order","Family","taxon","Site")]
  }

  # Map taxon, samples, and sites to s3 by File_Name
  tax.meta <- merge(unique(df), unique(s3), by = "File_Name", all.x = TRUE)
  
  # Group by WSCSDA and unique taxon counts, sites, and File_Name samples
  length_unique_rank <- paste0('length(unique(', rank, '))')
  tax.rank <- data.frame(tax.meta %>% 
                           group_by(WSCSDA, Order) %>% 
                           dplyr::summarize(rank = n_distinct(!!as.symbol(rank))))
  
  # Convert long to wide format and remove NAs
  tax.rank.wide <- dcast(tax.rank, WSCSDA ~ Order, value.var = "rank")
  tax.rank.wide[is.na(tax.rank.wide)] <- 0

  # Combine major taxa, WSCSDA, and coordinates
  tax.major.final <- merge(tax.rank.wide, simplified.centroid, by = "WSCSDA", all.x = TRUE)

}

# now get df at correct resolution for each plot of major taxa
df.esv <- make_minichart_df(tax.major, "GlobalESV")
df.species <- make_minichart_df(tax.major, "Species")
df.genus <- make_minichart_df(tax.major, "Genus")
df.family <- make_minichart_df(tax.major, "Family")
df.taxon <- make_minichart_df(tax.major, "taxon")

# ESV tab
my_popups <- paste("Watershed: ", prettyNum(df.esv$WSCSDA_EN, big.mark = ","), "<br>",
                   "Ephemeroptera: ", prettyNum(df.esv$Ephemeroptera, big.mark = ","), "<br>", 
                   "Odonata: ", prettyNum(df.esv$Odonata, big.mark = ","), "<br>", 
                   "Plecoptera: ", prettyNum(df.esv$Plecoptera, big.mark = ","), "<br>",
                   "Trichoptera: ", prettyNum(df.esv$Trichoptera, big.mark = ","), "<br>",
                   "Chironomidae: ", prettyNum(df.esv$Diptera, big.mark = ","), "<br>")

# Add color
pal <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')

# Generate leaflet map with pie charts
a.esv <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.esv$Long, lat = df.esv$Lat, type = "pie",
                chartdata = df.esv[, c("Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera", "Diptera")],
                popup = popupArgs(html = my_popups))
# Display pie chart map
a.esv

# Species tab
my_popups <- paste("Watershed: ", df.species$WSCSDA_EN, "<br>",
                   "Ephemeroptera: ", df.species$Ephemeroptera, "<br>", 
                   "Odonata: ", df.species$Odonata, "<br>", 
                   "Plecoptera: ", df.species$Plecoptera, "<br>",
                   "Trichoptera: ", df.species$Trichoptera, "<br>",
                   "Chironomidae: ", df.species$Diptera, "<br>")

# Generate leaflet map with pie charts
a.species <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.species$Long, lat = df.species$Lat, type = "pie",
                chartdata = df.species[, c("Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera", "Diptera")],
                popup = popupArgs(html = my_popups))
# Display pie chart map
a.species

# Genus Tab
my_popups <- paste("Watershed: ", df.genus$WSCSDA_EN, "<br>",
                   "Ephemeroptera: ", df.genus$Ephemeroptera, "<br>", 
                   "Odonata: ", df.genus$Odonata, "<br>", 
                   "Plecoptera: ", df.genus$Plecoptera, "<br>",
                   "Trichoptera: ", df.genus$Trichoptera, "<br>",
                   "Chironomidae: ", df.genus$Diptera, "<br>")

# Generate leaflet map with pie charts
a.genus <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.genus$Long, lat = df.genus$Lat, type = "pie",
                chartdata = df.genus[, c("Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera", "Diptera")],
                popup = popupArgs(html = my_popups))
# Display pie chart map
a.genus

# Family Tab
my_popups <- paste("Watershed: ", df.family$WSCSDA_EN, "<br>",
                   "Ephemeroptera: ", df.family$Ephemeroptera, "<br>", 
                   "Odonata: ", df.family$Odonata, "<br>", 
                   "Plecoptera: ", df.family$Plecoptera, "<br>",
                   "Trichoptera: ", df.family$Trichoptera, "<br>",
                   "Chironomidae: ", df.family$Diptera, "<br>")

# Generate leaflet map with pie charts
a.family <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.family$Long, lat = df.family$Lat, type = "pie",
                chartdata = df.family[, c("Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera", "Diptera")],
                popup = popupArgs(html = my_popups))
# Display pie chart map
a.family

# Taxon Tab
my_popups <- paste("Watershed: ",df.taxon$WSCSDA_EN, "<br>",
                   "Ephemeroptera: ", df.taxon$Ephemeroptera, "<br>", 
                   "Odonata: ", df.taxon$Odonata, "<br>", 
                   "Plecoptera: ", df.taxon$Plecoptera, "<br>",
                   "Trichoptera: ", df.taxon$Trichoptera, "<br>",
                   "Chironomidae: ", df.taxon$Diptera, "<br>")

# Generate leaflet map with pie charts
a.taxon <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.taxon$Long, lat = df.taxon$Lat, type = "pie",
                chartdata = df.taxon[, c("Ephemeroptera", "Odonata", "Plecoptera", "Trichoptera", "Diptera")],
                popup = popupArgs(html = my_popups))
# Display pie chart map
a.taxon


###############################
# Get the File_Name and Phylum
tax.ranks <- tax[c("File_Name", "Phylum")]
# Count taxa per phylum
phylum <- tax %>% count(tax$Phylum)
# Find the top 5 taxa for all the phyla
phylum <- head(phylum[order(-phylum$n),], 5)
topPhyla <- phylum$`tax$Phylum`

tax.top.phylum <- tax[tax$Phylum %in% topPhyla,]

make_minichart_df <- function(df, rank){
  # df options: tax.major
  # rank options: GlobalESV, Species, Genus, Family, taxon
  # cutoff options: NA, species sBP >= 0.70, genus gBP >= 0.30, species sBP >= 0.20, NA
  
  # Keep key fields, File_Name, Order, Family, rank, site
  # filter by bootstrap values where needed
  if (rank == "GlobalESV") {
    df <- df[, c("File_Name","Phylum","GlobalESV","Site")]
  } else if (rank == "Species") {
    df <- df[df$sBP >= 0.70,]
    df <- df[, c("File_Name","Phylum","Species","Site")]
  } else if (rank == "Genus") {
    df <- df[df$gBP >= 0.30,]
    df <- df[, c("File_Name","Phylum","Genus","Site")]
  } else if (rank == "Family") {
    df <- df[df$fBP >= 0.20,]
    df <- df[, c("File_Name","Phylum","Family","Site")]
  } else {
    df <- df[, c("File_Name","Phylum","taxon","Site")]
  }
  
  # Map taxon, samples, and sites to s3 by File_Name
  tax.meta <- merge(unique(df), unique(s3), by = "File_Name", all.x = TRUE)
  
  # Group by WSCSDA and unique taxon counts, sites, and File_Name samples
  length_unique_rank <- paste0('length(unique(', rank, '))')
  tax.rank <- data.frame(tax.meta %>% 
                           group_by(WSCSDA, Phylum) %>% 
                           dplyr::summarize(rank = n_distinct(!!as.symbol(rank))))
  
  # Convert long to wide format and remove NAs
  tax.rank.wide <- dcast(tax.rank, WSCSDA ~ Phylum, value.var = "rank")
  tax.rank.wide[is.na(tax.rank.wide)] <- 0
  
  # Combine major taxa, WSCSDA, and coordinates
  tax.major.final <- merge(tax.rank.wide, simplified.centroid, by = "WSCSDA", all.x = TRUE)
  
}

df.esv <- make_minichart_df(tax.top.phylum, "GlobalESV")
df.species <- make_minichart_df(tax.top.phylum, "Species")
df.genus <- make_minichart_df(tax.top.phylum, "Genus")
df.family <- make_minichart_df(tax.top.phylum, "Family")
df.taxon <- make_minichart_df(tax.top.phylum, "taxon")

# ESV tab
my_popups <- paste("Watershed: ", df.esv$WSCSDA_EN, "<br>",
                   paste(names(df.esv[2]), ": ", sep = ""), prettyNum(df.esv[, 2], big.mark = ","), "<br>", 
                   paste(names(df.esv[3]), ": ", sep = ""), prettyNum(df.esv[, 3], big.mark = ","), "<br>", 
                   paste(names(df.esv[4]), ": ", sep = ""), prettyNum(df.esv[, 4], big.mark = ","), "<br>",
                   paste(names(df.esv[5]), ": ", sep = ""), prettyNum(df.esv[, 5], big.mark = ","), "<br>",
                   paste(names(df.esv[6]), ": ", sep = ""), prettyNum(df.esv[, 6], big.mark = ","), "<br>")

# Generate leaflet map with pie charts
p.esv <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.esv$Long, lat = df.esv$Lat, type = "pie",
                chartdata = df.esv[2:6],
                popup = popupArgs(html = my_popups))
# Display pie chart map
p.esv

# Species tab
my_popups <- paste("Watershed: ", df.species$WSCSDA_EN, "<br>",
                   paste(names(df.species[2]), ": ", sep = ""), df.species[, 2], "<br>", 
                   paste(names(df.species[3]), ": ", sep = ""), df.species[, 3], "<br>", 
                   paste(names(df.species[4]), ": ", sep = ""), df.species[, 4], "<br>",
                   paste(names(df.species[5]), ": ", sep = ""), df.species[, 5], "<br>",
                   paste(names(df.species[6]), ": ", sep = ""), df.species[, 6], "<br>")

# Generate leaflet map with pie charts
p.species <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.species$Long, lat = df.species$Lat, type = "pie",
                chartdata = df.species[2:6],
                popup = popupArgs(html = my_popups))
# Display pie chart map
p.species

# Genus tab
my_popups <- paste("Watershed: ", df.genus$WSCSDA_EN, "<br>",
                   paste(names(df.genus[2]), ": ", sep = ""), df.genus[, 2], "<br>", 
                   paste(names(df.genus[3]), ": ", sep = ""), df.genus[, 3], "<br>", 
                   paste(names(df.genus[4]), ": ", sep = ""), df.genus[, 4], "<br>",
                   paste(names(df.genus[5]), ": ", sep = ""), df.genus[, 5], "<br>",
                   paste(names(df.genus[6]), ": ", sep = ""), df.genus[, 6], "<br>")

# Generate leaflet map with pie charts
p.genus <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.genus$Long, lat = df.genus$Lat, type = "pie",
                chartdata = df.genus[2:6],
                popup = popupArgs(html = my_popups))
# Display pie chart map
p.genus

# Family tab
my_popups <- paste("Watershed: ", df.family$WSCSDA_EN, "<br>",
                   paste(names(df.family[2]), ": ", sep = ""), df.family[, 2], "<br>", 
                   paste(names(df.family[3]), ": ", sep = ""), df.family[, 3], "<br>", 
                   paste(names(df.family[4]), ": ", sep = ""), df.family[, 4], "<br>",
                   paste(names(df.family[5]), ": ", sep = ""), df.family[, 5], "<br>",
                   paste(names(df.family[6]), ": ", sep = ""), df.family[, 6], "<br>")

# Generate leaflet map with pie charts
p.family <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.family$Long, lat = df.family$Lat, type = "pie",
                chartdata = df.family[2:6],
                popup = popupArgs(html = my_popups))
# Display pie chart map
p.family

# Taxon tab
my_popups <- paste("Watershed: ", df.taxon$WSCSDA_EN, "<br>",
                   paste(names(df.taxon[2]), ": ", sep = ""), df.taxon[, 2], "<br>", 
                   paste(names(df.taxon[3]), ": ", sep = ""), df.taxon[, 3], "<br>", 
                   paste(names(df.taxon[4]), ": ", sep = ""), df.taxon[, 4], "<br>",
                   paste(names(df.taxon[5]), ": ", sep = ""), df.taxon[, 5], "<br>",
                   paste(names(df.taxon[6]), ": ", sep = ""), df.taxon[, 6], "<br>")

# Generate leaflet map with pie charts
p.taxon <- leaflet(simplified) %>%
  addTiles() %>%
  # addProviderTiles(providers$Stamen.Terrain) %>%
  setView(lat = 60, lng = -95, zoom = 3) %>%
  addPolygons(data = simplified, color = "darkgrey", stroke = TRUE, 
              weight = 1, smoothFactor = 0.2,
              fillColor = pal, fillOpacity = 0.5
  )  %>%
  addMinicharts(lng = df.taxon$Long, lat = df.taxon$Lat, type = "pie",
                chartdata = df.taxon[2:6], popup = popupArgs(html = my_popups))
# Display pie chart map
p.taxon

