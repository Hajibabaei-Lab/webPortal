# Teresita M. Porter, April 26, 2022
# WWF watershed shapefiles from Mike Wright
# Data cleaned by Artin Mashayekhi
# plot richness as cloropleth

# read in shapefile
library(sf)

# plotting & data flipping
library(leaflet)
library(htmltools)
library(dplyr)

###############################
# read in a shapefile
system.time(wwf_read <- st_read("WWF_Subwatersheds/WSC_subwatersheds.shp", 
                                quiet = TRUE))

# id should really be a key to map watersheds from WWF to our samples (done below, remove id step here)
wwf_wgs84 <- wwf_read %>%
  st_buffer(0) %>% # Make invalid geometries valid
  st_transform(crs = 4326)  # Convert coordinates to WGS84 already in WGS84

# simplify so that map loads faster
simplified <- rmapshaper::ms_simplify(wwf_wgs84)




###############################
# read in taxonomy file
tax <- read.csv("STREAM_2018-2020_taxonomy.csv", header=TRUE, stringsAsFactors = FALSE)

# reformat SampleName in tax to match File_Name in s (below)
# replace "_" with "-" in tax
# remove run number and amplicon
# call it File_Name
tax$File_Name <- tax$SampleName
tax$File_Name <- gsub("-", "_", tax$File_Name)
tax <- cbind(tax, data.frame(do.call('rbind', strsplit(as.character(tax$File_Name),'_',fixed=TRUE))))
names(tax)[34:43] <- c("Project", "Group", "Substrate", "Site", "Sample", "Replicate", "Date", "Marker", "IlluminaRun", "Amplicon")
# put the field back together to match File_name in s (metadata file below)
tax$File_Name <- paste(tax$Project, tax$Group, tax$Substrate, tax$Site, tax$Sample, tax$Replicate, tax$Date, tax$Marker, sep="-")
# some File_Name records in metadata do not have Replicate field, 
# but these haven't been sequenced yet so don't worry about it for now

# create taxon field for assignments we can be confident in
# for 200bp COI query, 95% confidence at species rank (sBP >= 0.70), 
# 99% confidence genus (gBP >= 0.30), family (fBP >= 0.20), order+ no cutoff needed
# This set contains
tax$taxon <- ifelse(tax$sBP >= 0.70, paste(tax$Phylum, tax$Class, tax$Species, sep=";"),
                    ifelse(tax$gBP >= 0.30, paste(tax$Phylum, tax$Class, tax$Genus, sep=""),
                          ifelse(tax$fBP >= 0.20, paste(tax$Phylum, tax$Class, tax$Family, sep=";"),
                                 paste(tax$Phylum, tax$Class, tax$Order, sep=";"))))

# create a data frame that can be mapped to the shapefile via metadata
# ex. here index by File_Name, summarize by Site_Sample (pool across replicates), summarize by WSCSDA to get total richness per watershed

# this gives richness per file_name
tax.taxon <- tax[,c("File_Name", "taxon")]



###############################
# read in Artin's cleaned up metadata
# need to mapp Site.Sample to WSCSDA field from simplified shapefile (see below)
s <- read.csv("STREAM_Metadata_20Apr2022_AM_TP.csv", header=TRUE, stringsAsFactors = FALSE)
# manually fix all the File_Name without the Amplicon

# get unique lat long and label (Watershed or Site.Sample)
# Include File_Name to map taxonomy
# Include WSCSDA to map to shapefile
s2 <- unique(s[,c("File_Name", "WSCSDA", "Lat","Long")])

# drop row if NA, data is incomplete
# need to fill in the blanks
# EVERY sample MUST have lat & long, fix root problem in Excel with Mike's help
#s3 <- na.omit(s2) 
s3 <- s2[!is.na(s2$Lat),]
s3 <- s3[!is.na(s3$Long),]

# now map richness to s3 by File_Name
tax.meta <- merge(tax.richness, s3, by="File_Name", all.x=TRUE)

# # find all the malformed File_Name records and manually fix in STREAM_Metadata_20Apr2022_AM_TP.csv
# tax.meta[is.na(tax.meta$WSCSDA),]

# STREAM-DFONLX-B-LALD-000X-X-20201011-COI in metadata is missing info, but need it because we have sequences for it
# Artin to get this info from Mike, until then, drop
tax.meta <- tax.meta[!tax.meta$File_Name=="STREAM-DFONLX-B-LALD-000X-X-20201011-COI",]

# now group by WSCSDA and count unique taxa
tax.richness <- data.frame(tax.meta %>% group_by(WSCSDA) %>% dplyr::summarise(richness = n_distinct(taxon)))

# add richness to the shapefile
shp <- merge(simplified, tax.richness, by="WSCSDA", all.x=TRUE)
# set NA to zero
shp$richness[is.na(shp$richness)] <- 0



# create color palette
pal.tmp <- brewer.pal(3, "BuGn")

pal <- colorNumeric(
  palette = pal.tmp,
  domain = shp$richness)




# Basic choropleth with leaflet
m <- leaflet(simplified) %>% 
  addTiles()  %>% 
  setView( lat=10, lng=0 , zoom=2) %>%
  addPolygons(data=shp, color= "white", stroke = TRUE, weight = 1, smoothFactor = 0.2,
              fillColor = ~pal(richness), fillOpacity = 1,
              popup = paste("Watershed: ", shp$WSCSDA_EN, "<br>",
                            "Richness: ", shp$richness, "taxa<br>"))

m




# adapt this script to also plot number of "Sites: and Samples:" 
# per watershed in the popup box
# Sites, number of unique sites
# Samples, equivalent to number of unique filenames

# adapt to summarize major taxa such as Ephemeroptera, Plecoptera, Trichoptera, Odonata, Chironomidae, Other
# adapt to summarize number of taxa in each phylum
# ranks for these taxa can be found at ncbi taxonomy page, ex. https://www.ncbi.nlm.nih.gov/taxonomy/?term=Ephemroptera
# leaflet can be used to display this info as piecharts
https://cran.r-project.org/web/packages/leaflet.minicharts/vignettes/introduction.html

# put together a preliminary report in R markdown to get feedback from Mehrdad
https://rmarkdown.rstudio.com/


