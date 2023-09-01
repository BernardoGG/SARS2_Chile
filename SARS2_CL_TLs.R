################################################################################
####### SARS-CoV-2 Chile transmission lineage extraction and analysis ##########
################################################################################

# Source scripts and data
source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")

########################### Bernardo Gutierrez #################################

library(lme4)
library(janitor)
library(lubridate)
library(stringr)
library(stringi) # to remove accents from strings
library(seqinr)
library(ggplot2)
library(scales)
library(patchwork)
library(NatParksPalettes)
library(tidyverse)

############################### WRANGLING ######################################
### Import TL data ####
#### Files generated using fertree
## Import TLs function
read_TLs <- function(x){
  read.table(x, sep = "\t") |>
    row_to_names(row_number = 1) |>
    select(-tree, -source) |>
    mutate(ntaxa = as.numeric(ntaxa)) |>
    mutate(tmrca = as.numeric(tmrca)) |>
    mutate(ptmrca = as.numeric(ptmrca)) |>
    mutate(first_seen = as.numeric(first_seen)) |>
    mutate(last_seen = as.numeric(last_seen)) |>
    mutate(tmrca = decimal_date(last_CL_seq) - tmrca) |>
    mutate(ptmrca = decimal_date(last_CL_seq) - ptmrca) |>
    mutate(first_seen = decimal_date(last_CL_seq) - first_seen) |>
    mutate(last_seen = decimal_date(last_CL_seq) - last_seen) |>
    mutate(tmrca_date = as.Date(date_decimal(tmrca))) |>
    mutate(ptmrca_date = as.Date(date_decimal(ptmrca))) |>
    mutate(first_seen_date = as.Date(date_decimal(first_seen))) |>
    mutate(last_seen_date = as.Date(date_decimal(last_seen))) |>
    mutate(lineage = c(1, 1:tail(lineage, n = 1)+ 1))
}

## Alpha
alpha_tl <-
  read_TLs("Phylogenetics/DTA_Alpha_ISP/DTA_Alpha_ISP_transmission_lineages.tsv") |>
  mutate(lineage = paste0("Alpha_TL_", lineage))

## Gamma
gamma_tl <-
  read_TLs("Phylogenetics/DTA_Gamma_ISP/DTA_Gamma_ISP_transmission_lineages.tsv") |>
  mutate(lineage = paste0("Gamma_TL_", lineage))

## Lambda
lambda_tl <-
  read_TLs("Phylogenetics/DTA_Lambda_ISP/DTA_Lambda_ISP_transmission_lineages.tsv") |>
  mutate(lineage = paste0("Lambda_TL_", lineage))

## Mu
mu_tl <-
  read_TLs("Phylogenetics/DTA_Mu_ISP/DTA_Mu_ISP_transmission_lineages.tsv") |>
  mutate(lineage = paste0("Mu_TL_", lineage))

## Delta
delta_tl <-
  read_TLs("Phylogenetics/DTA_Delta_ISP/DTA_Delta_ISP_transmission_lineages.tsv") |>
  mutate(lineage = paste0("Delta_TL_", lineage))

#### Unify data sets
## Full data set
CL_fertree <- data.frame(variant = c(rep("Alpha", nrow(alpha_tl)), rep("Gamma", nrow(gamma_tl)),
                                     rep("Lambda", nrow(lambda_tl)), rep("Mu", nrow(mu_tl)),
                                     rep("Delta", nrow(delta_tl))),
                            lineage = c(alpha_tl$lineage, gamma_tl$lineage,
                                        lambda_tl$lineage, mu_tl$lineage, delta_tl$lineage),
                            tmrca = c(alpha_tl$tmrca, gamma_tl$tmrca,
                                      lambda_tl$tmrca, mu_tl$tmrca, delta_tl$tmrca),
                            tmrca_date = c(alpha_tl$tmrca_date, gamma_tl$tmrca_date,
                                           lambda_tl$tmrca_date, mu_tl$tmrca_date, delta_tl$tmrca_date),
                            ptmrca = c(alpha_tl$ptmrca, gamma_tl$ptmrca,
                                       lambda_tl$ptmrca, mu_tl$ptmrca, delta_tl$ptmrca),
                            ptmrca_date = c(alpha_tl$ptmrca_date, gamma_tl$ptmrca_date,
                                            lambda_tl$ptmrca_date, mu_tl$ptmrca_date, delta_tl$ptmrca_date),
                            first_seen = c(alpha_tl$first_seen, gamma_tl$first_seen,
                                           lambda_tl$first_seen, mu_tl$first_seen, delta_tl$first_seen),
                            first_seen_date = c(alpha_tl$first_seen_date, gamma_tl$first_seen_date,
                                                lambda_tl$first_seen_date, mu_tl$first_seen_date, delta_tl$first_seen_date),
                            last_seen = c(alpha_tl$last_seen, gamma_tl$last_seen,
                                          lambda_tl$last_seen, mu_tl$last_seen, delta_tl$last_seen),
                            last_seen_date = c(alpha_tl$last_seen_date, gamma_tl$last_seen_date,
                                               lambda_tl$last_seen_date, mu_tl$last_seen_date, delta_tl$last_seen_date),
                            taxa = c(alpha_tl$taxa, gamma_tl$taxa, lambda_tl$taxa, mu_tl$taxa, delta_tl$taxa),
                            ntaxa = c(alpha_tl$ntaxa, gamma_tl$ntaxa, lambda_tl$ntaxa, mu_tl$ntaxa, delta_tl$ntaxa))


## Singletons
CL_singletons <- CL_fertree |> filter(ntaxa == 1) |> mutate(seen = first_seen) |> mutate(seen_date = first_seen_date) |>
  select(-ntaxa, -first_seen, -first_seen_date, -last_seen, -last_seen_date)

## Transmission lineages
CL_TLs <- CL_fertree |> filter(ntaxa != 1) |> select(-taxa) |> mutate(lineage = as.factor(lineage))

#### Annotate with airport metadata
## Transmission lineages/singletons containing airport sequences
airport_TLs <- vector()
for(i in 1:nrow(airport_metadata)){
  n <- CL_fertree$lineage[grepl(airport_metadata$accession[i],
                                CL_fertree$taxa, fixed = TRUE)]
  airport_TLs <- c(airport_TLs, n) 
}
airport_TLs <- unique(airport_TLs)

## Transmission lineages/singletons containing airport sequences
CL_TLs$airport <- CL_TLs$lineage %in% airport_TLs
CL_singletons$airport <- CL_singletons$lineage %in% airport_TLs

## Annotate sequences found in TLs (not singletons)
# Airport surveillance
airport_circ <- vector()
for(i in 1:nrow(airport_metadata)){
  n <- ifelse(TRUE %in% str_detect(airport_metadata$accession[i],
                                   CL_fertree$taxa[CL_fertree$ntaxa > 1]),
              "TL",
              "singleton")
  airport_circ <- c(airport_circ, n) 
}
airport_metadata$circulation <- airport_circ


# Community surveillance
community_circ <- vector()
for(i in 1:nrow(community_metadata)){
  n <- ifelse(TRUE %in% str_detect(community_metadata$accession[i],
                                   CL_fertree$taxa[CL_fertree$ntaxa > 1]),
              "TL",
              "singleton")
  community_circ <- c(community_circ, n) 
}
community_metadata$circulation <- community_circ

## Add columns that reference epidemic weeks
airport_metadata$epiyear <- epiyear(airport_metadata$collection_date)
airport_metadata$epiweek <- epiweek(airport_metadata$collection_date)

community_metadata$epiyear <- epiyear(community_metadata$collection_date)
community_metadata$epiweek <- epiweek(community_metadata$collection_date)

CL_TLs$epiyear_tmrca <- epiyear(CL_TLs$tmrca_date)
CL_TLs$epiweek_tmrca <- epiweek(CL_TLs$tmrca_date)
CL_TLs$epiyear_first_seen <- epiyear(CL_TLs$first_seen_date)
CL_TLs$epiweek_first_seen <- epiweek(CL_TLs$first_seen_date)
CL_TLs$epiyear_ptmrca <- epiyear(CL_TLs$ptmrca_date)
CL_TLs$epiweek_ptmrca <- epiweek(CL_TLs$ptmrca_date)

CL_singletons$epiyear_tmrca <- epiyear(CL_singletons$tmrca_date)
CL_singletons$epiweek_tmrca <- epiweek(CL_singletons$tmrca_date)
CL_singletons$epiyear_collection <- epiyear(CL_singletons$seen_date)
CL_singletons$epiweek_collection <- epiweek(CL_singletons$seen_date)
CL_singletons$epiyear_ptmrca <- epiyear(CL_singletons$ptmrca_date)
CL_singletons$epiweek_ptmrca <- epiweek(CL_singletons$ptmrca_date)

# Estimate the date in which every epiweek starts - airport data
tmp <- vector()
for(i in 1:nrow(airport_metadata)){
  ifelse(airport_metadata$epiyear[i] == "2020",
         d <- as.Date("2019-12-29") + airport_metadata$epiweek[i]*7,
         d <- as.Date("2021-01-03") + airport_metadata$epiweek[i]*7)
  tmp <- c(as.Date(tmp), d)
}
airport_metadata$epiweek_start <- tmp
rm(d)

# Estimate the date in which every epiweek starts - community data
tmp <- vector()
for(i in 1:nrow(community_metadata)){
  ifelse(community_metadata$epiyear[i] == "2020",
         d <- as.Date("2019-12-29") + community_metadata$epiweek[i]*7,
         d <- as.Date("2021-01-03") + community_metadata$epiweek[i]*7)
  tmp <- c(as.Date(tmp), d)
}
community_metadata$epiweek_start <- tmp
rm(d)

# Estimate the date in which every epiweek starts - Chile TLs data
tmp <- vector()
for(i in 1:nrow(CL_TLs)){
  ifelse(CL_TLs$epiyear_tmrca[i] == "2020",
         d <- as.Date("2019-12-29") + CL_TLs$epiweek_tmrca[i]*7,
         d <- as.Date("2021-01-03") + CL_TLs$epiweek_tmrca[i]*7)
  tmp <- c(as.Date(tmp), d)
}
CL_TLs$epiweek_tmrca_start <- tmp
rm(d)

tmp <- vector()
for(i in 1:nrow(CL_TLs)){
  ifelse(CL_TLs$epiyear_ptmrca[i] == "2020",
         d <- as.Date("2019-12-29") + CL_TLs$epiweek_ptmrca[i]*7,
         d <- as.Date("2021-01-03") + CL_TLs$epiweek_ptmrca[i]*7)
  tmp <- c(as.Date(tmp), d)
}
CL_TLs$epiweek_ptmrca_start <- tmp
rm(d)

tmp <- vector()
for(i in 1:nrow(CL_TLs)){
  ifelse(CL_TLs$epiyear_first_seen[i] == "2020",
         d <- as.Date("2019-12-29") + CL_TLs$epiweek_first_seen[i]*7,
         d <- as.Date("2021-01-03") + CL_TLs$epiweek_first_seen[i]*7)
  tmp <- c(as.Date(tmp), d)
}
CL_TLs$epiweek_first_seen_start <- tmp
rm(d)

# Estimate the date in which every epiweek starts - Chile singletons data
tmp <- vector()
for(i in 1:nrow(CL_singletons)){
  ifelse(CL_singletons$epiyear_tmrca[i] == "2020",
         d <- as.Date("2019-12-29") + CL_singletons$epiweek_tmrca[i]*7,
         d <- as.Date("2021-01-03") + CL_singletons$epiweek_tmrca[i]*7)
  tmp <- c(as.Date(tmp), d)
}
CL_singletons$epiweek_tmrca_start <- tmp
rm(d)

tmp <- vector()
for(i in 1:nrow(CL_singletons)){
  ifelse(CL_singletons$epiyear_ptmrca[i] == "2020",
         d <- as.Date("2019-12-29") + CL_singletons$epiweek_ptmrca[i]*7,
         d <- as.Date("2021-01-03") + CL_singletons$epiweek_ptmrca[i]*7)
  tmp <- c(as.Date(tmp), d)
}
CL_singletons$epiweek_ptmrca_start <- tmp
rm(d)

tmp <- vector()
for(i in 1:nrow(CL_singletons)){
  ifelse(CL_singletons$epiyear_collection[i] == "2020",
         d <- as.Date("2019-12-29") + CL_singletons$epiweek_collection[i]*7,
         d <- as.Date("2021-01-03") + CL_singletons$epiweek_collection[i]*7)
  tmp <- c(as.Date(tmp), d)
}
CL_singletons$epiweek_collection_start <- tmp
rm(d, tmp)

### Sequence geocoding for ISP data set ####
## Read geocoding data
geocode_CL_full <- as.data.frame(fread('Data/CL_geocoded_ISP_seq.csv')) |>
  select(accession, collection_date, comuna, variant, surveillance, X, Y) |>
  filter(accession %in% chile_isp_metadata$`Accession ID`) |>
  rename(long = X, lat = Y)

## Generate geocoding file with tip names
geocode_CL_full <- geocode_CL_full |> arrange(accession)
name_lib <- chile_isp_metadata |> arrange(`Accession ID`)

geocode_CL_full$virus_name <- name_lib$`Virus name`[
  geocode_CL_full$accession %in% name_lib$`Accession ID`
]

geocode_CL_full$taxa <- str_c(geocode_CL_full$virus_name, "|",
                              geocode_CL_full$accession, "|",
                              geocode_CL_full$collection_date)

chile_isp_metadata$Taxa <- str_c(chile_isp_metadata$`Virus name`, "|",
                                 chile_isp_metadata$`Accession ID`, "|",
                                 chile_isp_metadata$Date)

## Remove missing data
geocode_CL <- geocode_CL_full |> drop_na(long, lat)

## Create vector of accession numbers with missing data
geocode_CL_missing <- geocode_CL_full$accession[geocode_CL_full$comuna=="desconocido"]

## Filter down to just include transmission lineages
taxaunlist <- vector()
for(i in 1:length(taxalist)){
  x <- unlist(taxalist[[i]])
  taxaunlist <- c(taxaunlist, x)
}
taxaunlist <- data.frame(taxa = taxaunlist) |> arrange(taxa)

taxa_df <- taxaunlist
taxaunlist <- str_split_fixed(taxaunlist$taxa, fixed("|"), 3) |>
  as.data.frame() |>
  rename(virus_name = V1, accession = V2, collection_date = V3)
taxaunlist <- bind_cols(taxa_df, taxaunlist)

geocode_TLs <- geocode_CL |>
  filter(accession %in% taxaunlist$accession)

geocode_TLs_file <- geocode_TLs |> select(taxa, lat, long)

## Save geocoding metadata file
write.table(geocode_TLs_file,
            file = "Data/SARS2_CL_TLs_geocode.tsv",
            row.names = FALSE,
            col.names = c("taxa", "lat", "long"),
            sep = "\t",
            quote = FALSE)


### Extract TL subtrees ####
## Transmission lineages >= 60 sequences extracted
CL_large_TLs <- CL_fertree |> filter(ntaxa >= 60)  # <- Expect 20 TLs

## Extract lists of taxa per each TL
taxalist <- list()
for(i in 1:nrow(CL_large_TLs)){
  taxa = strsplit(CL_large_TLs$taxa[i], ";") |>
    unlist() |>
    str_trim()
  taxalist[[i]] = taxa
}
names(taxalist) <- c(CL_large_TLs$lineage)

lapply(1:length(taxalist),
       function(i) write.table(taxalist[[i]],
                             file = paste0(names(taxalist[i]), ".tsv"),
                             row.names = FALSE,
                             sep = "\t",
                             quote = FALSE)) # <-- Export TL tips to extract with fertree

## Extract sequences of taxa per each TL
# List with all VOC sequences
VOCs_fasta <- c(Alpha_ISP_fasta, Gamma_ISP_fasta, Lambda_ISP_fasta,
                Mu_ISP_fasta, Delta_ISP_fasta)

# Object with TL names, sequence names and sequence EPI_ISL
TL_names <- names(taxalist) |> paste0("_aln")
seq_names <- names(VOCs_fasta)
seq_acce <- str_split_fixed(names(VOCs_fasta), fixed("|"),3) |>
  as.data.frame() |>
  select(V2) |> rename(accession = V2)

# Rename sequences as EPI_ISL
names(VOCs_fasta) <- seq_acce$accession
accelist <- list()
for(i in 1:nrow(CL_large_TLs)){
  taxa = strsplit(CL_large_TLs$taxa[i], ";") |>
    unlist() |>
    str_trim() |>
    str_split_fixed(fixed("|"), 3) |>
    as.data.frame() |>
    select(V2)
  accelist[[i]] = taxa
}
names(accelist) <- c(CL_large_TLs$lineage)

# Create object containing all ISP sequences
TL_seqs <- list()
for(i in 1:length(accelist)){
  TL_seqs[[i]] <- VOCs_fasta[names(VOCs_fasta) %in% unlist(accelist[[i]]) &
                              names(VOCs_fasta) %in% geocode_CL$accession]
}
names(TL_seqs) <- TL_names

x <- VOCs_fasta[names(VOCs_fasta) %in% unlist(accelist[[2]]) &
             names(VOCs_fasta) %in% geocode_CL$accession]

# Export FASTA files for each TL
# Only use EPI-ISL codes as sequence names for consistency
lapply(1:length(TL_seqs),
       function(i) write.fasta(TL_seqs[[i]],
                               names = names(TL_seqs[[i]]),
                               file.out = paste0("Data/Transmission_lineages/",
                                             names(TL_seqs[i]),
                                             ".fasta")))

# Export taxa names to prune trees in fertree
lapply(1:length(TL_seqs),
       function(i) write.table(
         taxaunlist$taxa[taxaunlist$accession %in% names(TL_seqs[[i]])],
         file = paste0("Phylogenetics/SC2_CL_TLs/",
                       names(accelist[i]), "_pruned.tsv"),
         row.names = FALSE,
         col.names = FALSE,
         sep = "\t",
         quote = FALSE))

# Export lat-long for individual transmission lineages
geocode_TLs_extract <- geocode_TLs_file |>
  mutate(taxa = str_extract(taxa, "(EPI_ISL_\\d+)"))

lapply(1:length(TL_seqs),
       function(i) write.table(
         geocode_TLs_extract[geocode_TLs_extract$taxa %in%
                               names(TL_seqs[[i]]),1:3],
         file = paste0("Phylogenetics/SC2_CL_TLs/",
                       names(accelist[i]), "_geocode.tsv"),
         row.names = FALSE,
         col.names = c("taxa", "lat", "long"),
         sep = "\t",
         quote = FALSE))

# Export collection dates for sequences in individual transmission lineages
lapply(1:length(TL_seqs),
       function(i) write.table(
         chile_isp_metadata[chile_isp_metadata$`Accession ID` %in% names(TL_seqs[[i]]),c(2,6)],
         file = paste0("Phylogenetics/SC2_CL_TLs/",
                       names(accelist[i]), "_dates.tsv"),
         row.names = FALSE,
         col.names = c("taxa", "date"),
         sep = "\t",
         quote = FALSE))

### Import domestic spread data frames ####
#### Files generated using custom scripts (Python)
## Import invasion trees data frame function
invasions_folder <- "Phylogenetics/SC2_CL_TLs/BEAST_TL_continuous_phylogeo/"
read_invasions <- function(variant, tl){
  read.table(paste0(invasions_folder, variant, "_TL_", tl, ".tsv"), sep = "\t") |>
    row_to_names(row_number = 1) |>
    select(-head_node, -tail_node, -length, -geo_distance,
           -head_comuna_amb, -tail_comuna_amb) |>
    mutate(lineage =  paste0(variant, "_TL_", tl)) |>
    mutate(variant =  variant) |>
    mutate(head_lat_3395 = as.numeric(head_lat_3395)) |>
    mutate(head_long_3395 = as.numeric(head_long_3395)) |>
    mutate(tail_lat_3395 = as.numeric(tail_lat_3395)) |>
    mutate(tail_long_3395 = as.numeric(tail_long_3395)) |>
    mutate(head_lat_4326 = as.numeric(head_lat_4326)) |>
    mutate(head_long_4326 = as.numeric(head_long_4326)) |>
    mutate(tail_lat_4326 = as.numeric(tail_lat_4326)) |>
    mutate(tail_long_4326 = as.numeric(tail_long_4326)) |>
    mutate(head_dec_date = as.numeric(head_dec_date)) |>
    mutate(tail_dec_date = as.numeric(tail_dec_date)) |>
    mutate(head_date = as.Date(head_date)) |>
    mutate(tail_date = as.Date(tail_date))
}

## Alpha
invasion_a13 <- read_invasions("Alpha", "13")

## Gamma
invasion_g35 <- read_invasions("Gamma", "35")
invasion_g55 <- read_invasions("Gamma", "55")

## Lambda
invasion_l90 <- read_invasions("Lambda", "90")
invasion_l95 <- read_invasions("Lambda", "95")
invasion_l103 <- read_invasions("Lambda", "103")
invasion_l104 <- read_invasions("Lambda", "104")
invasion_l107 <- read_invasions("Lambda", "107")

## Mu
invasion_m15 <- read_invasions("Mu", "15")
invasion_m26 <- read_invasions("Mu", "26")
invasion_m61 <- read_invasions("Mu", "61")
invasion_m73 <- read_invasions("Mu", "73")
invasion_m82 <- read_invasions("Mu", "82")

## Delta
invasion_d10 <- read_invasions("Delta", "10")
invasion_d64 <- read_invasions("Delta", "64")
invasion_d166 <- read_invasions("Delta", "166")
invasion_d179 <- read_invasions("Delta", "179")
invasion_d186 <- read_invasions("Delta", "186")
invasion_d191 <- read_invasions("Delta", "191")
invasion_d254 <- read_invasions("Delta", "254")

invasions_list <- list(invasion_a13, invasion_g55, invasion_g35, invasion_l103,
                       invasion_l104, invasion_l107, invasion_l90, invasion_l95,
                       invasion_m15, invasion_m26, invasion_m61, invasion_m73,
                       invasion_m82, invasion_d10, invasion_d179, invasion_d186,
                       invasion_d191, invasion_d254, invasion_d64, invasion_d166)

## Create a single data frame with all invasion trees
CL_invasions <- Reduce(
  function(x, y, ...) bind_rows(x, y, ...), 
  invasions_list
)

## Create a single data frame with all invasion trees
comunas_to_fix_head <- data.frame(
  to_fix = unique(
    CL_invasions$head_comuna[
      !CL_invasions$head_comuna %in% unique(
        lockdowns$comuna_residencia)]) |> sort(),
  fixed = c("Aysen", "Casablanca", "Coyhaique", "Puyehue", "La Union",
            "Llay-Llay", "Llanquihue", "Canela", "Quilaco", "Quilleco",
            "Ranquil", "Ranquil", "Requinoa", "San Rosendo",
            "San Jose de Maipo", "Til Til", "Yerbas Buenas"))

## Clean-up comuna names to match lockdown data
comunas_to_fix_tail <- data.frame(
  to_fix = unique(
    CL_invasions$tail_comuna[
      !CL_invasions$tail_comuna %in% unique(
        lockdowns$comuna_residencia)]) |> sort(),
  fixed = c("Aysen", "Casablanca", "Coyhaique", "Puyehue", "La Union",
            "La Union", "Llay-Llay", "Llanquihue", "Canela", "Quilaco",
            "Quilleco", "Ranquil", "Ranquil", "Requinoa", "San Rosendo",
            "San Jose de Maipo", "Til Til", "Yerbas Buenas"))

for(i in 1:nrow(CL_invasions)){
  for(j in 1:nrow(comunas_to_fix_head)){
    CL_invasions$head_comuna[i] <-
      ifelse(CL_invasions$head_comuna[i] %in% comunas_to_fix_head$to_fix[j],
             comunas_to_fix_head$fixed[j],
             CL_invasions$head_comuna[i])
  }
}

for(i in 1:nrow(CL_invasions)){
  for(j in 1:nrow(comunas_to_fix_tail)){
    CL_invasions$tail_comuna[i] <-
      ifelse(CL_invasions$tail_comuna[i] %in% comunas_to_fix_tail$to_fix[j],
             comunas_to_fix_tail$fixed[j],
             CL_invasions$tail_comuna[i])
  }
}

rm(comunas_to_fix_head, comunas_to_fix_tail)

## Add column with lockdown tier over time period of transition
CL_invasions <- CL_invasions |>
  left_join(lockdowns, by = c("head_comuna" = "comuna_residencia",
                              "head_date" = "Fecha")) |>
  select(-region_residencia, head_lockdown_tier = Paso)

CL_invasions <- CL_invasions |>
  left_join(lockdowns, by = c("tail_comuna" = "comuna_residencia",
                              "tail_date" = "Fecha")) |>
  select(-region_residencia, tail_lockdown_tier = Paso)

## Add column with recoded lockdown tier
CL_invasions <- CL_invasions |>
  mutate(head_lockdown_tier_recoded = case_when(
    head_lockdown_tier == 1 ~ "Full lockdown",
    head_lockdown_tier == 2 ~ "Weekend lockdown",
    head_lockdown_tier == 3 | head_lockdown_tier == 4 ~ "No lockdown")) |>
  mutate(tail_lockdown_tier_recoded = case_when(
    tail_lockdown_tier == 1 ~ "Full lockdown",
    tail_lockdown_tier == 2 ~ "Weekend lockdown",
    tail_lockdown_tier == 3 | tail_lockdown_tier == 4 ~ "No lockdown")) |>
  mutate(head_lockdown_tier_recoded = as.factor(head_lockdown_tier_recoded)) |>
  mutate(tail_lockdown_tier_recoded = as.factor(tail_lockdown_tier_recoded))
CL_invasions$head_lockdown_tier_recoded <- factor(
  CL_invasions$head_lockdown_tier_recoded,
  levels = c("Full lockdown", "Weekend lockdown", "No lockdown"))
CL_invasions$tail_lockdown_tier_recoded <- factor(
  CL_invasions$tail_lockdown_tier_recoded,
  levels = c("Full lockdown", "Weekend lockdown", "No lockdown"))

################################# PLOTS ########################################
### Plot TLs in Chile ####
# Singleton collection times per variant
ggplot(CL_singletons, aes(x = seen_date, fill = variant)) +
  geom_histogram(binwidth = 7, color = 'black', size = 0.2) +
  labs(x = "Sample collection date",
       y = "No. of singletons observed per week",
       fill = "Variant") +
  theme_minimal() + scale_fill_manual(values = vocs_colors) +
  xlim(min(CL_TLs$tmrca_date), max(CL_TLs$last_seen_date)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# TL TMRCAs per variant
CL_TLs


# TL sizes per variant
CL_TLs |>
  ggplot(aes(x = fct_reorder(CL_TLs$lineage, CL_TLs$ntaxa, .desc = TRUE),
             y = ntaxa, fill = variant)) + geom_col() +
  labs(x = "Transmission lineages", y = "No. of sequences per TL",
       fill = "Variant") +
  theme_bw() + scale_fill_manual(values = vocs_colors) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l", outside = TRUE) +
  coord_cartesian(clip = "off") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.x = element_blank())

# TL persistence in time colored by variant
CL_TLs |> ggplot() +
  geom_linerange(aes(x = fct_reorder(lineage, first_seen_date, .desc = FALSE),
                     ymin = first_seen_date, ymax = last_seen_date, color = variant),
                 linewidth = 0.7, lineend='round') +
  geom_linerange(aes(x = fct_reorder(lineage, first_seen_date, .desc = FALSE),
                     ymin = tmrca_date, ymax = first_seen_date, color = variant),
                 linewidth = 0.4, alpha = 0.3) +
  geom_point(aes(x = fct_reorder(lineage, first_seen_date, .desc = FALSE),
                 y = last_seen_date, color = variant, size = log10(ntaxa)), alpha = 0.3) +
  geom_point(aes(x = fct_reorder(lineage, first_seen_date, .desc = FALSE),
                 y = tmrca_date, color = variant), shape = 18) +
  labs(x = "Transmission lineages", y = "Persistence over time",
       color = "Variant", size = "Number of genomes\n(Log10)") +
  theme_bw() + scale_color_manual(values = vocs_colors) +
  coord_flip() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.y = element_blank())

ggsave("Figures/Fig_TLs_persistence.pdf", dpi = 300,
       height = 12, width = 10, units = "in", bg = "white")

ggsave("Figures/Fig_TLs_persistence.png", dpi = 300,
       height = 12, width = 10, units = "in", bg = "white")

# TLs detected at the airport by TMRCA date
ggplot(CL_TLs[CL_TLs$airport == TRUE,]) +
  geom_histogram(aes(x = tmrca_date, fill = variant),
                 binwidth = 7, color = 'black', size = 0.2) +
  labs(x = "TL TMRCA date",
       y = "No. of TL introductions detected during \nairport surveillance per week",
       fill = "Variant") +
  theme_minimal() + scale_fill_manual(values = vocs_colors) +
  xlim(min(CL_TLs$tmrca_date), max(CL_TLs$last_seen_date)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggplot(CL_singletons[CL_singletons$airport == TRUE,]) +
  geom_histogram(aes(x = tmrca_date, fill = variant),
                 binwidth = 7, color = 'black', size = 0.2) +
  labs(x = "Singleton TMRCA date",
       y = "No. of singleton introductions detected \nduring airport surveillance per week",
       fill = "Variant") +
  theme_minimal() + scale_fill_manual(values = vocs_colors) +
  xlim(min(CL_TLs$tmrca_date), max(CL_TLs$last_seen_date)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggplot(CL_singletons[CL_singletons$airport == TRUE,]) +
  geom_histogram(aes(x = seen_date, fill = variant),
                 binwidth = 7, color = 'black', size = 0.2) +
  labs(x = "Singleton sample collection date",
       y = "No. of singleton introductions detected \nduring airport surveillance per week",
       fill = "Variant") +
  theme_minimal() + scale_fill_manual(values = vocs_colors) +
  xlim(min(CL_TLs$tmrca_date), max(CL_TLs$last_seen_date)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Sequences associated to TLs or as singletons per surveillance scheme
ggplot(airport_metadata[airport_metadata$variant == "alpha" |
                          airport_metadata$variant == "gamma" |
                          airport_metadata$variant == "lambda" |
                          airport_metadata$variant == "mu" |
                          airport_metadata$variant == "delta",]) +
  geom_freqpoly(aes(x = collection_date, color = circulation),
                 binwidth = 7, size = 1.2) +
  theme_minimal() +
  xlim(min(CL_TLs$tmrca_date), max(CL_TLs$last_seen_date)) +
  labs(
    title = "Airport surveillance sequences associated to \ntransmission lineages and singletons",
    x = "Sequence collection date",
    y = "Count")

ggplot(community_metadata[community_metadata$variant == "alpha" |
                            community_metadata$variant == "gamma" |
                            community_metadata$variant == "lambda" |
                            community_metadata$variant == "mu" |
                            community_metadata$variant == "delta",]) +
  geom_freqpoly(aes(x = collection_date, color = circulation),
                binwidth = 7, size = 1.2) +
  theme_minimal() +
  xlim(min(CL_TLs$tmrca_date), max(CL_TLs$last_seen_date)) +
  labs(
    title = "Community surveillance sequences associated to \ntransmission lineages and singletons",
    x = "Sequence collection date",
    y = "Count")

ggplot(CL_singletons) +
  geom_histogram(aes(x = tmrca_date, fill = variant))

### Manuscript plots
a <- CL_fertree |> select(-taxa) |>
  ggplot(aes(x = tmrca_date, fill = variant, color = variant)) +
  geom_density(linewidth = 0.6, alpha = 0.05) +
  labs(x = "TL TMRCA",
       y = "Proportion of introductions\nover time per VOC",
       color = "Variant") +
  theme_minimal() + scale_color_manual(values = vocs_colors) +
  scale_fill_manual(values = vocs_colors) +
  scale_x_date(breaks = "1 month", labels = NULL,
               limits = c(range(SARS2_CL_epiweeks$epiweek_start)[1],
                          range(SARS2_CL_epiweeks$epiweek_start)[2])) +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(fill = "none")

ggsave("Figures/CL_introductions_density_plot.png", dpi = 300,
       height = 2, width = 7, units = "in", bg = "white")

ggsave("Figures/CL_introductions_density_plot.pdf", dpi = 300,
       height = 2, width = 7, units = "in", bg = "white")


### Plot invasions according to lockdown tiers ####
# Normalising factor by transmission lineage
norm_tl <- CL_invasions$lineage |> table()

# Normalising factor by variant
norm_voc <- CL_invasions$variant |> table()

# Normalising factor by lockdown tier
norm_tier_head <- CL_invasions$head_lockdown_tier_recoded |> table()
norm_tier_tail <- CL_invasions$tail_lockdown_tier_recoded |> table()

# Normalising factor by comuna
norm_comuna_head <- CL_invasions$head_comuna |> table()
norm_comuna_tail <- CL_invasions$tail_comuna |> table()

# Normalising factor by comuna, imports vs internal counts
norm_comuna_tail_within <- CL_invasions |>
  filter(head_comuna == tail_comuna) |> select(tail_comuna) |> table()
norm_comuna_tail_between <- CL_invasions |>
  filter(head_comuna != tail_comuna) |> select(tail_comuna) |> table()

# Viral movements (normalised) by lockdown tier
within_h <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(variant, lineage, head_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_voc[
    as.character(variant)]) |>
  ggplot(aes(x = head_lockdown_tier_recoded,
             y = norm_count, fill = variant)) +
  geom_boxplot(position = position_dodge(0.6), alpha = 0.3,
               linewidth = 0.3, width = 0.4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(0.6), stroke = 0.5,
               stackratio = 0.6, dotsize = 0.4) +
  scale_fill_manual(values = vocs_colors) +
  labs(title = "Outgoing viral movements within comunas", x = "",
       y = "Proportion of movements out of all\ntransitions inferred for each variant") +
  theme_classic() + theme(axis.text.x = element_blank(), legend.position = "none")

within_t <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(variant, lineage, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_voc[
    as.character(variant)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = norm_count, fill = variant)) +
  geom_boxplot(position = position_dodge(0.6), alpha = 0.3,
               linewidth = 0.3, width = 0.4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(0.6), stroke = 0.5,
               stackratio = 0.6, dotsize = 0.4) +
  scale_fill_manual(values = vocs_colors) +
  labs(title = "Incoming viral movements within comunas", x = "",
       y = "") +
  theme_classic() + theme(axis.text.x = element_blank())

between_h <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(variant, lineage, head_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_voc[
    as.character(variant)]) |>
  ggplot(aes(x = head_lockdown_tier_recoded,
             y = norm_count, fill = variant)) +
  geom_boxplot(position = position_dodge(0.6), alpha = 0.3,
               linewidth = 0.3, width = 0.4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(0.6), stroke = 0.5,
               stackratio = 0.6, dotsize = 0.4) +
  scale_fill_manual(values = vocs_colors) +
  labs(title = "Outgoing viral movements between comunas", x = "",
       y = "Proportion of movements out of all\ntransitions inferred for each variant") +
  theme_classic() + theme(legend.position = "none")

between_t <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(variant, lineage, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_voc[
    as.character(variant)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = norm_count, fill = variant)) +
  geom_boxplot(position = position_dodge(0.6), alpha = 0.3,
               linewidth = 0.3, width = 0.4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(0.6), stroke = 0.5,
               stackratio = 0.6, dotsize = 0.4) +
  scale_fill_manual(values = vocs_colors) +
  labs(title = "Incoming viral movements between comunas", x = "",
       y = "") +
  theme_classic()

(within_h | within_t) / (between_h | between_t)


# Viral imports (normalised) into comunas by lockdown tier
# Generate data frame
CL_invasions

# Plot
within_t2 <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(variant, lineage, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(variant = as.factor(variant)) |>
  mutate(variant = fct_relevel(variant, "Delta", after = Inf)) |>
  mutate(norm_count = count / norm_tier_tail[
    as.character(tail_lockdown_tier_recoded)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = norm_count, fill = variant)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(0.8), stroke = 0.5,
               stackratio = 0.8, dotsize = 0.8) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.2,
               linewidth = 0.3, width = 0.6, colour = "#555555") +
  stat_summary(fun = median, geom = 'line',
              aes(group = variant, color = variant),
              position = position_dodge(width = 0.8), alpha = 0.1,
              linewidth = 1.5) +
  scale_fill_manual(values = vocs_colors) +
  scale_colour_manual(values = vocs_colors) +
  scale_y_continuous(trans = 'log10') +
  labs(title = "Viral movements within comunas", x = "",
       y = "Proportion of inferred\nviral movements (Log10)") +
  theme_classic() + theme(axis.text.x = element_blank(),
                          legend.position = "none")

between_t2 <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(variant, lineage, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(variant = as.factor(variant)) |>
  mutate(variant = fct_relevel(variant, "Delta", after = Inf)) |>
  mutate(norm_count = count / norm_tier_tail[
    as.character(tail_lockdown_tier_recoded)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = norm_count, fill = variant)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(0.8), stroke = 0.5,
               stackratio = 0.8, dotsize = 0.8) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.2,
               linewidth = 0.3, width = 0.6, colour = "#555555") +
  stat_summary(fun = median, geom = 'line',
               aes(group = variant, color = variant),
               position = position_dodge(width = 0.8), alpha = 0.1,
               linewidth = 1.5) +
  scale_fill_manual(values = vocs_colors) +
  scale_colour_manual(values = vocs_colors) +
  scale_y_continuous(trans = 'log10') +
  labs(title = "Viral movements between comunas", x = "", colour = "",
       y = "Proportion of inferred\nviral movements (Log10)", fill = "Variant") +
  theme_classic() + theme(axis.text.x = element_text(size = 10),
                          legend.position = "bottom") + guides(color = "none")

within_t2 / between_t2

ggsave("Figures/SC2_CL_lockdown_tiers_transitions_imports.pdf", dpi = 300,
       height = 7.15, width = 6, bg = "white")
ggsave("Figures/SC2_CL_lockdown_tiers_transitions_imports.png", dpi = 300,
       height = 7.15, width = 6, bg = "white")


## Viral movements into individual comunas by lockdown tier
# Imports
imports_c <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_comuna_tail_between[
    as.character(tail_comuna)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = count, fill = tail_lockdown_tier_recoded)) +
  geom_violin(alpha = 0.5, linewidth = 0.6, width = 0.2, colour = "#555555",
               outlier.shape = NA) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', stroke = 0.5,
               stackratio = 0.4, dotsize = 0.6, alpha = 0.1,
               fill = "#555555", position = position_jitter(
                 width = 0.2, height = 0.05)) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "", colour = "",
       y = "Counts of viral importations") +
  theme_classic() + theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          legend.position = "none")

# Movements within comunas
movements_c <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_comuna_tail_within[
    as.character(tail_comuna)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = count, fill = tail_lockdown_tier_recoded)) +
  geom_violin(alpha = 0.5, linewidth = 0.6, width = 0.2, colour = "#555555",
               outlier.shape = NA) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', stroke = 0.5,
               stackratio = 0.4, dotsize = 0.6, alpha = 0.1,
               fill = "#555555", position = position_jitter(
                 width = 0.2, height = 0.05)) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "", colour = "",
       y = "Counts of viral movements\nwithin comunas") +
  theme_classic() + theme(axis.text.x = element_text(size = 12),
                          legend.position = "none")

imports_c / movements_c

ggsave("Figures/SC2_CL_lockdown_tiers_by_comuna.pdf", dpi = 300,
       height = 7.15, width = 6, bg = "white")
ggsave("Figures/SC2_CL_lockdown_tiers_by_comuna.png", dpi = 300,
       height = 7.15, width = 6, bg = "white")


################################ EXPORTS #######################################
### Export TL summaries and FASTA files ####
## Singletons
write.table(CL_singletons)

## Transmission lineages
write.table(CL_TLs)

## FASTA files
write.fasta(Alpha_TL_13_fasta, names = names(Alpha_TL_13_fasta),
            file.out = "Data/Transmission_lineages/Alpha_TL_13.fasta")










chile_isp_metadata[chile_isp_metadata$Taxa=="hCoV-19/Chile/RM-266418/2020|EPI_ISL_1167921|2020-12-30",]
geocode_CL_full[geocode_CL_full$taxa=="hCoV-19/Chile/RM-266418/2020|EPI_ISL_1167921|2020-12-30",]

chile_isp_metadata[chile_isp_metadata$Taxa=="hCoV-19/Chile/LR-262486/2020|EPI_ISL_1167923|2020-12-22",]
geocode_CL_full[geocode_CL_full$taxa=="hCoV-19/Chile/LR-262486/2020|EPI_ISL_1167923|2020-12-22",]




################################ SANDBOX #######################################
## Viral importations
CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = count, fill = tail_lockdown_tier_recoded)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.6, width = 0.2, colour = "#555555",
               outlier.shape = NA) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', stroke = 0.5,
               stackratio = 0.2, dotsize = 0.6, alpha = 0.1,
               fill = "#555555", position = position_jitter(
                 width = 0.2, height = 0.1)) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "", colour = "",
       y = "Inferred viral importations (Log10)") +
  theme_classic() + theme(axis.text.x = element_text(size = 10),
                          legend.position = "none")

CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_comuna_tail[
    as.character(tail_comuna)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = norm_count, fill = tail_lockdown_tier_recoded)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.6, width = 0.2, colour = "#555555",
               outlier.shape = NA) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', stroke = 0.5,
               stackratio = 0.2, dotsize = 0.6, alpha = 0.1,
               fill = "#555555", position = position_jitter(
                 width = 0.2, height = 0.1)) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "", colour = "",
       y = "Proportion of inferred viral importations relative\nto all importations per comuna (Log10)") +
  theme_classic() + theme(axis.text.x = element_text(size = 10),
                          legend.position = "none")

## Movements within comunas
CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = count, fill = tail_lockdown_tier_recoded)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.6, width = 0.2, colour = "#555555",
               outlier.shape = NA) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', stroke = 0.5,
               stackratio = 0.2, dotsize = 0.6, alpha = 0.1,
               fill = "#555555", position = position_jitter(
                 width = 0.2, height = 0.1)) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "", colour = "",
       y = "Inferred viral movements\nwithin comunas (Log10)") +
  theme_classic() + theme(axis.text.x = element_text(size = 10),
                          legend.position = "none")

CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(norm_count = count / norm_comuna_tail[
    as.character(tail_comuna)]) |>
  ggplot(aes(x = tail_lockdown_tier_recoded,
             y = norm_count, fill = tail_lockdown_tier_recoded)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.6, width = 0.2, colour = "#555555",
               outlier.shape = NA) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', stroke = 0.5,
               stackratio = 0.2, dotsize = 0.6, alpha = 0.1,
               fill = "#555555", position = position_jitter(
                 width = 0.2, height = 0.1)) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "", colour = "",
       y = "Proportion of inferred viral movements relative\nto all movements per comuna (Log10)") +
  theme_classic() + theme(axis.text.x = element_text(size = 10),
                          legend.position = "none")


CL_invasions_MEM_within <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_MEM_between <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))


anova(lm(count ~ tail_lockdown_tier_recoded, CL_invasions_MEM_within))
anova(lm(count ~ tail_lockdown_tier_recoded, CL_invasions_MEM_between))
m1 <-lm(count ~ tail_lockdown_tier_recoded, data = CL_invasions_MEM_within)
summary(m1)

anova(lm(norm_count ~ tail_lockdown_tier_recoded, CL_invasions_MEM_within))
anova(lm(norm_count ~ tail_lockdown_tier_recoded, CL_invasions_MEM_between))


(lmer(count ~ tail_lockdown_tier_recoded + (tail_lockdown_tier_recoded | tail_comuna),
     CL_invasions_MEM_within)

ggplot(data = CL_invasions_MEM_between,
       aes(x = tail_lockdown_tier_recoded, y = norm_count,
           col = tail_comuna, group = tail_comuna)) + #to add the colours for different classes
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter") + #to add some random noise for plotting purposes
  theme_minimal() +
  theme(legend.position = "none") +
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8) + # to add regression line
  labs(title    = "Viral movements vs. Lockdown tier",
       subtitle = "add colours for different comunas and regression lines")


CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna) |>
  summarise(count = n()) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)])) |>
  select(norm_count) |> hist(breaks = 20)






colnames(CL_invasions)

ggplot(CL_invasions) + geom_histogram(aes(x = tail_date - head_date)) +
  labs(x = "Branch lenghts (days)", y = "") +
  theme_minimal()

lockdowns |>
  group_by(comuna_residencia) |>
  arrange(Fecha)


lockdowns$lockdown_tier[lockdowns$comuna_residencia == "Santiago"]

CL_invasions |> arrange(tail_date) |>
  filter(tail_comuna == "Santiago") |>
  select(tail_lockdown_tier_recoded) |>
  as.vector()

lockdowns$