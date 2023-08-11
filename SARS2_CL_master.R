################################################################################
########### SARS-CoV-2 genomic surveillance in Chile ###########################
################################################################################

############## Rhys Inward, Bernardo Gutierrez #################################

## Package load
library(readxl) # to read xlsx file with specific sheet 
library(stringi) # to remove accents from strings
library(viridis) # for nice ggplot colors
library(wesanderson) # ggplot color help
library(tidyverse)
library(safejoin)
library(data.table)
library(dplyr)
library(zoo) 
library(countrycode) 
library(purrr) 
library(readr) 
library(stringr)   
library(tidyr) 
library(usdata)
library(seqinr)
library(lubridate)
library(ggplot2)
library(ape)
library(patchwork)

## Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

## Set working directory
setwd("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/")

## Create folder path for results
folres <- (path = "/Users/user/Documents/SARS2_Chile_local/Results/")


## Source main functions to run scripts
files.sources =
  list.files(path = "/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/Main")

for (i in 1:length(files.sources)) {
  source(paste0(c("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/Main/",
                  files.sources[i]), collapse = ''))
}


#### PROCESS CHILE ISP METADATA ################################################

# Load Chile metadata
# (Updated data; all pre-processing has been done)
# Community and airport surveillance data sets loaded separately
community_metadata <-
  as.data.frame(fread('Data_rhys/community_ISP_missing_updated_n7957.qced.v1.tsv',
                      drop = c("gender")))
community_metadata$surveillance <- 'Community'

airport_metadata <-
  as.data.frame(fread('Data_rhys/airport_ISP_missing_updated_n788.qced.v1.tsv'))
airport_metadata$surveillance <- 'Airport'

# Bind airport and community data sets
chile_metadata_combined <- rbind(community_metadata,airport_metadata)

# Process Chile metadata
chile_isp_metadata <- filter(GISAID.processed,`Accession ID` %in%
                               chile_metadata_combined$accession)

# Add column to indicate where sequences were collected 
chile_isp_metadata <- left_join(chile_isp_metadata, chile_metadata_combined,
                                by = c('Accession ID' = 'accession'))
chile_isp_metadata <- dplyr::select(chile_isp_metadata,c(1:6,10,14))

# Remove duplicates
chile_isp_metadata <-
  chile_isp_metadata[!duplicated(chile_isp_metadata$`Accession ID`),]

# Call sequence Variants from Pango lineage designation 
chile_isp_metadata <-
  chile_isp_metadata %>%
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(chile_isp_metadata$`Pango lineage`[
                                 grep("AY.",
                                      chile_isp_metadata$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(chile_isp_metadata$`Pango lineage`[
                                 grep("P.1.",
                                      chile_isp_metadata$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

# Add column specifying epiweek
chile_isp_metadata <- chile_isp_metadata %>%
  mutate(Date = as.Date(Date, "%Y-%m-%d"),
         epiweek = epiweek(Date))

# Add column specifying cumulative epiweek
chile_isp_metadata <- chile_isp_metadata %>%
  mutate(c_epiweek = epiweek(Date),
         c_epiweek = ifelse(year(Date) == 2021 &
                            epiweek < 53, 53 + epiweek, epiweek),
         c_epiweek = ifelse(year(Date) == 2022 &
                            epiweek < 53, 106 + epiweek, epiweek),
         c_epiweek = ifelse(year(Date) == 2022 &
                            epiweek == 158, -52 + epiweek, epiweek))

# Add column specifying epiweek start date
chile_isp_metadata <- chile_isp_metadata %>%
  mutate(epiyear = epiyear(Date))

tmp <- vector()
for(i in 1:nrow(chile_isp_metadata)){
  ifelse(chile_isp_metadata$epiyear[i] == "2020",
         d <- as.Date("2019-12-29") + chile_isp_metadata$epiweek[i]*7,
         d <- as.Date("2021-01-03") + chile_isp_metadata$epiweek[i]*7)
  tmp <- c(as.Date(tmp), d)
}
chile_isp_metadata$epiweek_start <- tmp
rm(d)

# Remove accents and clean up region names
chile_isp_metadata$State <- stri_trans_general(chile_isp_metadata$State,
                                               "Latin-ASCII")
chile_isp_metadata$State[
  grep("Higgins", chile_isp_metadata$State)] <- "O Higgins"
chile_isp_metadata$State[
  grep("Magallanes", chile_isp_metadata$State)] <- "Magallanes"
chile_isp_metadata$State[
  grep("Parinacota", chile_isp_metadata$State)] <- "Arica y Parinacota"
chile_isp_metadata$State[
  grep("Araucania", chile_isp_metadata$State)] <- "Araucania"
chile_isp_metadata$State[
  grep("Antofagasta", chile_isp_metadata$State)] <- "Antofagasta"
chile_isp_metadata$State[
  grep("Maule", chile_isp_metadata$State)] <- "El Maule"

# Create Variant-specific data sets
chile_metadata_alpha <- chile_isp_metadata[chile_isp_metadata$Variant=="Alpha",]
chile_metadata_gamma <- chile_isp_metadata[chile_isp_metadata$Variant=="Gamma",]
chile_metadata_lambda <- chile_isp_metadata[chile_isp_metadata$Variant=="Lambda",]
chile_metadata_mu <- chile_isp_metadata[chile_isp_metadata$Variant=="Mu",]
chile_metadata_delta <- chile_isp_metadata[chile_isp_metadata$Variant=="Delta",]


#### PROCESS GISAID METADATA ###################################################

# Load GISAID metadata
# Data set downloaded on 2022-06-30
metadata <- as.data.frame(fread('Data_rhys/metadata_2022_06_30.tsv',
                                drop = c("Type", "Sequence length","Host",
                                         "Patient age", "Gender","Clade",
                                         "Pangolin version", "Variant",
                                         "AA Substitutions", "Submission date",
                                         "Is reference?", "Is complete?",
                                         "Is high coverage?", "Is low coverage?",
                                         "N-Content", "GC-Content",
                                         "Additional location information")))


# Process metadata using custom function
GISAID.processed <- GISAID_process(metadata)

# Define date range from Chile ISP metadata
first_CL_seq <- as.Date(min(chile_isp_metadata$Date))
last_CL_seq <- as.Date(max(chile_isp_metadata$Date))

# Filter GISAID metadata to only contain sequences from ISP study period
GISAID_final <- GISAID.processed %>%
  filter(Date >= first_CL_seq & Date <= last_CL_seq)

GISAID_final <- GISAID_final %>%
  mutate(Date = as.Date(Date, "%Y-%m-%d"),
         epiweek = epiweek(Date),
         epiweek = ifelse(year(Date) == 2021 &
                            epiweek < 53, 53 + epiweek, epiweek),
         epiweek = ifelse(year(Date) == 2022 &
                            epiweek < 53, 106 + epiweek, epiweek),
         epiweek = ifelse(year(Date) == 2022 &
                            epiweek == 158, -52 + epiweek, epiweek))


#### PROCESS CHILE GISAID METADATA #############################################

# Filter global metadata to match Chile metadata dates and
# include column specifying VOI/VOC variant designation
GISAID_CL <- GISAID_final[GISAID_final$Country == "Chile",] %>%
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(GISAID_final$`Pango lineage`[
                                 grep("AY.", GISAID_final$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(GISAID_final$`Pango lineage`[
                                 grep("P.1.", GISAID_final$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

# Remove accents and clean up regions
GISAID_CL$State <-stri_trans_general(GISAID_CL$State, "Latin-ASCII")
GISAID_CL$State[grep("Higgins", GISAID_CL$State)] <- "O Higgins"
GISAID_CL$State[grep("Magallanes", GISAID_CL$State)] <- "Magallanes"
GISAID_CL$State[grep("Parinacota", GISAID_CL$State)] <- "Arica y Parinacota"
GISAID_CL$State[grep("Antofagasta", GISAID_CL$State)] <- "Antofagasta"
GISAID_CL$State[grep("Araucania", GISAID_CL$State)] <- "Araucania"
GISAID_CL$State[grep("Maule", GISAID_CL$State)] <- "El Maule"

# Create Variant-specific data sets
GISAID_CL_alpha <- GISAID_CL[GISAID_CL$Variant=="Alpha",]
GISAID_CL_gamma <- GISAID_CL[GISAID_CL$Variant=="Gamma",]
GISAID_CL_lambda <- GISAID_CL[GISAID_CL$Variant=="Lambda",]
GISAID_CL_mu <- GISAID_CL[GISAID_CL$Variant=="Mu",]
GISAID_CL_delta <- GISAID_CL[GISAID_CL$Variant=="Delta",]


#### PROCESS PRELIMINARY TRANSMISSION LINEAGES #################################

# Load Transmission lineage data
# Data set version 30 June (single tree, algorithmic matching)
tl <- read.csv("30June_TLs/30June_chile_TL_summaries_v2.csv") %>%
  select(-taxa)
tl$tmrca <- ymd(tl$tmrca)
tl$first_sample_date <- ymd(tl$first_sample_date)
tl$last_sample_date <- ymd(tl$last_sample_date)

tl_nosingletons <- tl[tl$singleton=="False",] %>% select(-singleton)
tl_nosingletons$earliest_taxa_state_comp[
  tl_nosingletons$earliest_taxa_state_comp == "{'airport': 1, 'community': 1}"] <-
  "Undetermined"
tl_nosingletons$earliest_taxa_state_comp[
  grep("airport", tl_nosingletons$earliest_taxa_state_comp)] <- "Airport"
tl_nosingletons$earliest_taxa_state_comp[
  grep("community", tl_nosingletons$earliest_taxa_state_comp)] <- "Community"


#### UTILITIES #################################################################

### Create VOC colour palette
vocs_colors <- c("Alpha" = "#A8201A", "Delta" = "#CA5D22",  "Gamma" = "#EC9A29",
                 "Lambda" = "#0F8B8D", "Mu" = "#143642", "Other" = "#A9ABB3")

### Create South American countries colour palette
latam <- c("Argentina", "Brazil", "Peru", "Uruguay", "Paraguay",
           "Ecuador", "Colombia", "Venezuela", "Bolivia")
latam_palette <- NatParksPalettes$Acadia[1] %>% unlist()
names(latam_palette) <- latam


