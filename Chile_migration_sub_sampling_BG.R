################################################################################
########### SARS-CoV-2 Chile genomic subsampling ###############################
################################################################################

############## Rhys Inward, Bernardo Gutierrez #################################

# Clear working environment
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#set WD
setwd("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/")

# Folder path for results
folres <- (path = "/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/Results/")


# Main functions to run script
files.sources = list.files(path = "/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/Main")
for (i in 1:length(files.sources)) {
  source(paste0(c("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/Main/",
                  files.sources[i]), collapse = ''))
  }

#Load packages

library(tidyverse)
library(stringi)
library(safejoin)
library(data.table)
library(lubridate)
library(dplyr)
library(zoo) 
library(countrycode) 
library(purrr) 
library(readr) 
library(stringr)   
library(tidyr) 
library(usdata)
library(seqinr)
library(ggplot2)
library(ape)

#### MIGRATION INFORMED SUB-SAMPLING ###########################################

# Estimate data set sizes per variant
numseq_ISP_alpha <- nrow(chile_metadata_alpha)
numseq_GISAIDCL_alpha <- nrow(GISAID_CL_alpha)

numseq_ISP_gamma <- nrow(chile_metadata_gamma)
numseq_GISAIDCL_gamma <- nrow(GISAID_CL_gamma)

numseq_ISP_lambda <- nrow(chile_metadata_lambda)
numseq_GISAIDCL_lambda <- nrow(GISAID_CL_lambda)

numseq_ISP_mu <- nrow(chile_metadata_mu)
numseq_GISAIDCL_mu <- nrow(GISAID_CL_mu)

numseq_ISP_delta <- nrow(chile_metadata_delta)
numseq_GISAIDCL_delta <- nrow(GISAID_CL_delta)

# We aim for a 1:1:1 ratio between Chile sequences, G-data informed sequences 
# and sequences from countries outside LATAM with direct flights to Chile.
# G-data version: **CHECK**

# Load migration to Chile data
migration_to_chile <-
  as.data.frame(fread('Data_rhys/total_movements_to_chile.csv'))

# Scale the number of target sequences by the exportation intensity
# Estimate the required number of sequences from each location
migration_to_chile$num_alpha_ISP <-
  migration_to_chile$exportation_intensity * round((numseq_ISP_alpha))
migration_to_chile$num_gamma_ISP <-
  migration_to_chile$exportation_intensity * round((numseq_ISP_gamma))
migration_to_chile$num_lambda_ISP <-
  migration_to_chile$exportation_intensity * round((numseq_ISP_lambda))
migration_to_chile$num_mu_ISP <-
  migration_to_chile$exportation_intensity * round((numseq_ISP_mu))
migration_to_chile$num_delta_ISP <-
  migration_to_chile$exportation_intensity * round((numseq_ISP_delta))

migration_to_chile$num_alpha_GISAID <-
  migration_to_chile$exportation_intensity * round((numseq_GISAIDCL_alpha))
migration_to_chile$num_gamma_GISAID <-
  migration_to_chile$exportation_intensity * round((numseq_GISAIDCL_gamma))
migration_to_chile$num_lambda_GISAID <-
  migration_to_chile$exportation_intensity * round((numseq_GISAIDCL_lambda))
migration_to_chile$num_mu_GISAID <-
  migration_to_chile$exportation_intensity * round((numseq_GISAIDCL_mu))
migration_to_chile$num_delta_GISAID <-
  migration_to_chile$exportation_intensity * round((numseq_GISAIDCL_delta))

# Add the expected sequences by country over the variant circulation time span 
# to get a time-invariant subsample
CL_migration_alpha_ISP <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_alpha_ISP = sum(num_alpha_ISP))

CL_migration_gamma_ISP <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_gamma_ISP = sum(num_gamma_ISP))

CL_migration_lambda_ISP <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_lambda_ISP = sum(num_lambda_ISP))

CL_migration_mu_ISP <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_mu_ISP = sum(num_mu_ISP))

CL_migration_delta_ISP <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_delta_ISP = sum(num_delta_ISP))


CL_migration_alpha_GISAID <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_alpha_GISAID = sum(num_alpha_GISAID))

CL_migration_gamma_GISAID <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_gamma_GISAID = sum(num_gamma_GISAID))

CL_migration_lambda_GISAID <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_lambda_GISAID = sum(num_lambda_GISAID))

CL_migration_mu_GISAID <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_mu_GISAID = sum(num_mu_GISAID))

CL_migration_delta_GISAID <- migration_to_chile |>
  group_by(origin_admin_0) |>
  dplyr :: summarise(num_delta_GISAID = sum(num_delta_GISAID))

# Round to the nearest whole number
CL_migration_alpha_ISP$num_alpha_ISP <- ceiling(CL_migration_alpha_ISP$num_alpha_ISP)
CL_migration_gamma_ISP$num_gamma_ISP <- ceiling(CL_migration_gamma_ISP$num_gamma_ISP)
CL_migration_lambda_ISP$num_lambda_ISP <- ceiling(CL_migration_lambda_ISP$num_lambda_ISP)
CL_migration_mu_ISP$num_mu_ISP <- ceiling(CL_migration_mu_ISP$num_mu_ISP)
CL_migration_delta_ISP$num_delta_ISP <- ceiling(CL_migration_delta_ISP$num_delta_ISP)

CL_migration_alpha_GISAID$num_alpha_GISAID <- ceiling(CL_migration_alpha_GISAID$num_alpha_GISAID)
CL_migration_gamma_GISAID$num_gamma_GISAID <- ceiling(CL_migration_gamma_GISAID$num_gamma_GISAID)
CL_migration_lambda_GISAID$num_lambda_GISAID <- ceiling(CL_migration_lambda_GISAID$num_lambda_GISAID)
CL_migration_mu_GISAID$num_mu_GISAID <- ceiling(CL_migration_mu_GISAID$num_mu_GISAID)
CL_migration_delta_GISAID$num_delta_GISAID <- ceiling(CL_migration_delta_GISAID$num_delta_GISAID)


# Select country-specific metadata from top five high-flow countries
Argentina_sequences <- filter (GISAID_final, Country == 'Argentina')
Bolivia_sequences <- filter (GISAID_final, Country == 'Bolivia')
Brazil_sequences <- filter (GISAID_final, Country == 'Brazil')
Colombia_sequences <- filter (GISAID_final, Country == 'Colombia')
Peru_sequences <- filter (GISAID_final, Country == 'Peru')

# Add Variant identification in country-specific data sets
Argentina_sequences <- Argentina_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Argentina_sequences$`Pango lineage`[
                                 grep("AY.", Argentina_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Argentina_sequences$`Pango lineage`[
                                 grep("P.1.", Argentina_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

Bolivia_sequences <- Bolivia_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Bolivia_sequences$`Pango lineage`[
                                 grep("AY.", Bolivia_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Bolivia_sequences$`Pango lineage`[
                                 grep("P.1.", Bolivia_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

Brazil_sequences <- Brazil_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Brazil_sequences$`Pango lineage`[
                                 grep("AY.", Brazil_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Brazil_sequences$`Pango lineage`[
                                 grep("P.1.", Brazil_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))


Colombia_sequences <- Colombia_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Colombia_sequences$`Pango lineage`[
                                 grep("AY.", Colombia_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Colombia_sequences$`Pango lineage`[
                                 grep("P.1.", Colombia_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

Peru_sequences <- Peru_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Peru_sequences$`Pango lineage`[
                                 grep("AY.", Peru_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Peru_sequences$`Pango lineage`[
                                 grep("P.1.", Peru_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

# Sample uniformly across time, based on G-data estimates
# Alpha: sampling range of 36 weeks
wk_alpha <- ceiling(as.numeric(max(chile_metadata_alpha$Date) -
                                 min(chile_metadata_alpha$Date))/7)

# Gamma: sampling range of 37 weeks
wk_gamma <- ceiling(as.numeric(max(chile_metadata_gamma$Date) -
                                 min(chile_metadata_gamma$Date))/7)

# Lambda: sampling range of 36 weeks
wk_lambda <- ceiling(as.numeric(max(chile_metadata_lambda$Date) -
                                 min(chile_metadata_lambda$Date))/7)

# Mu: sampling range of 26 weeks
wk_mu <- ceiling(as.numeric(max(chile_metadata_mu$Date) -
                                 min(chile_metadata_mu$Date))/7)

# Delta: sampling range of 17 weeks
wk_delta <- ceiling(as.numeric(max(chile_metadata_delta$Date) -
                                 min(chile_metadata_delta$Date))/7)


### Sampling per variant per country
set.seed(5)

## Alpha
# Sample Argentina, Alpha, ISP
sample <- Argentina_sequences[Argentina_sequences$Variant=="Alpha",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_alpha_ISP[CL_migration_alpha_ISP$origin_admin_0=="Argentina",2])/wk_alpha)
sample$count <- ifelse(sample$count>max, max, sample$count)

ARG_ISP_alpha <- as.data.frame(Argentina_sequences[Argentina_sequences$Variant=="Alpha",] |>
                           group_by(epiweek) |>
                           slice(sample(min(max, n()))))

# Sample Bolivia, Alpha => No Alpha sequences from Bolivia available

# Sample Brazil, Alpha, ISP
sample <- Brazil_sequences[Brazil_sequences$Variant=="Alpha",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_alpha_ISP[CL_migration_alpha_ISP$origin_admin_0=="Brazil",2])/wk_alpha)
sample$count <- ifelse(sample$count>max, max, sample$count)

BRA_ISP_alpha <- as.data.frame(Brazil_sequences[Brazil_sequences$Variant=="Alpha",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n())))) # <-- Brazil is over-represented with 1 seq per epiweek

# Sample Colombia, Alpha, ISP
sample <- Colombia_sequences[Colombia_sequences$Variant=="Alpha",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_alpha_ISP[CL_migration_alpha_ISP$origin_admin_0=="Colombia",2])/wk_alpha)
sample$count <- ifelse(sample$count>max, max, sample$count)

COL_ISP_alpha <- as.data.frame(Colombia_sequences[Colombia_sequences$Variant=="Alpha",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Peru, Alpha, ISP => All Alpha sequences from Peru included to match expected number
PER_ISP_alpha <- Peru_sequences[Peru_sequences$Variant=="Alpha",]


# Sample Argentina, Alpha, GISAID
sample <- Argentina_sequences[Argentina_sequences$Variant=="Alpha",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_alpha_GISAID[CL_migration_alpha_GISAID$origin_admin_0=="Argentina",2])/wk_alpha)
sample$count <- ifelse(sample$count>max, max, sample$count)

ARG_GISAID_alpha <- as.data.frame(Argentina_sequences[Argentina_sequences$Variant=="Alpha",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Bolivia, Alpha => No Alpha sequences from Bolivia available

# Sample Brazil, Alpha, GISAID
sample <- Brazil_sequences[Brazil_sequences$Variant=="Alpha",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_alpha_GISAID[CL_migration_alpha_GISAID$origin_admin_0=="Brazil",2])/wk_alpha)
sample$count <- ifelse(sample$count>max, max, sample$count)

BRA_GISAID_alpha <- as.data.frame(Brazil_sequences[Brazil_sequences$Variant=="Alpha",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n())))) # <-- Brazil is over-represented with 1 seq per epiweek

# Sample Colombia, Alpha, GISAID
sample <- Colombia_sequences[Colombia_sequences$Variant=="Alpha",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_alpha_GISAID[CL_migration_alpha_GISAID$origin_admin_0=="Colombia",2])/wk_alpha)
sample$count <- ifelse(sample$count>max, max, sample$count)

COL_GISAID_alpha <- as.data.frame(Colombia_sequences[Colombia_sequences$Variant=="Alpha",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Peru, Alpha, GISAID => All Alpha sequences from Peru included to match expected number
PER_GISAID_alpha <- Peru_sequences[Peru_sequences$Variant=="Alpha",]


## Gamma
# Sample Argentina, Gamma, ISP
sample <- Argentina_sequences[Argentina_sequences$Variant=="Gamma",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_gamma_ISP[CL_migration_gamma_ISP$origin_admin_0=="Argentina",2])/wk_gamma)
sample$count <- ifelse(sample$count>max, max, sample$count)

ARG_ISP_gamma <- as.data.frame(Argentina_sequences[Argentina_sequences$Variant=="Gamma",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Bolivia, Gamma => All Gamma sequences from Bolivia used
BOL_ISP_gamma <- Bolivia_sequences[Bolivia_sequences$Variant=="Gamma",]

# Sample Brazil, Gamma, ISP
sample <- Brazil_sequences[Brazil_sequences$Variant=="Gamma",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_gamma_ISP[CL_migration_gamma_ISP$origin_admin_0=="Brazil",2])/wk_gamma)
sample$count <- ifelse(sample$count>max, max, sample$count)

BRA_ISP_gamma <- as.data.frame(Brazil_sequences[Brazil_sequences$Variant=="Gamma",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n())))) # <-- Brazil is over-represented

# Sample Colombia, Gamma, ISP => All Gamma sequences from Peru included to match expected number
COL_ISP_gamma <- Colombia_sequences[Colombia_sequences$Variant=="Gamma",]

# Sample Peru, Gamma, ISP
sample <- Peru_sequences[Peru_sequences$Variant=="Gamma",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_gamma_ISP[CL_migration_gamma_ISP$origin_admin_0=="Peru",2])/wk_gamma)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_ISP_gamma <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Gamma",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))


# Sample Argentina, Gamma, GISAID
sample <- Argentina_sequences[Argentina_sequences$Variant=="Gamma",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_gamma_GISAID[CL_migration_gamma_GISAID$origin_admin_0=="Argentina",2])/wk_gamma)
sample$count <- ifelse(sample$count>max, max, sample$count)

ARG_GISAID_gamma <- as.data.frame(Argentina_sequences[Argentina_sequences$Variant=="Gamma",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Bolivia, Gamma => All Gamma sequences from Bolivia used
BOL_GISAID_gamma <- Bolivia_sequences[Bolivia_sequences$Variant=="Gamma",]

# Sample Brazil, Gamma, ISP
sample <- Brazil_sequences[Brazil_sequences$Variant=="Gamma",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_gamma_GISAID[CL_migration_gamma_GISAID$origin_admin_0=="Brazil",2])/wk_gamma)
sample$count <- ifelse(sample$count>max, max, sample$count)

BRA_GISAID_gamma <- as.data.frame(Brazil_sequences[Brazil_sequences$Variant=="Gamma",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n())))) # <-- Brazil is over-represented

# Sample Colombia, Gamma, ISP => All Gamma sequences from Peru included to match expected number
COL_GISAID_gamma <- Colombia_sequences[Colombia_sequences$Variant=="Gamma",]

# Sample Peru, Gamma, ISP
sample <- Peru_sequences[Peru_sequences$Variant=="Gamma",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_gamma_GISAID[CL_migration_gamma_GISAID$origin_admin_0=="Peru",2])/wk_gamma)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_GISAID_gamma <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Gamma",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))


## Lambda
# Sample Argentina, Lambda, ISP
sample <- Argentina_sequences[Argentina_sequences$Variant=="Lambda",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_lambda_ISP[CL_migration_lambda_ISP$origin_admin_0=="Argentina",2])/wk_lambda)
sample$count <- ifelse(sample$count>max, max, sample$count)

ARG_ISP_lambda <- as.data.frame(Argentina_sequences[Argentina_sequences$Variant=="Lambda",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Bolivia, Lambda => All Lambda sequences from Bolivia used
BOL_ISP_lambda <- Bolivia_sequences[Bolivia_sequences$Variant=="Lambda",]

# Sample Brazil, Lambda => All Lambda sequences from Bolivia used
BRA_ISP_lambda <- Brazil_sequences[Brazil_sequences$Variant=="Lambda",]

# Sample Colombia, Lambda, ISP => All Lambda sequences from Colombia included to match expected number
COL_ISP_lambda <- Colombia_sequences[Colombia_sequences$Variant=="Lambda",]

# Sample Peru, Lambda, ISP
sample <- Peru_sequences[Peru_sequences$Variant=="Lambda",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_lambda_ISP[CL_migration_lambda_ISP$origin_admin_0=="Peru",2])/wk_lambda)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_ISP_lambda <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Lambda",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))


# Sample Argentina, Lambda, GISAID
sample <- Argentina_sequences[Argentina_sequences$Variant=="Lambda",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_lambda_GISAID[CL_migration_lambda_GISAID$origin_admin_0=="Argentina",2])/wk_lambda)
sample$count <- ifelse(sample$count>max, max, sample$count)

ARG_GISAID_lambda <- as.data.frame(Argentina_sequences[Argentina_sequences$Variant=="Lambda",] |>
                                  group_by(epiweek) |>
                                  slice(sample(min(max, n()))))

# Sample Bolivia, Lambda => All Lambda sequences from Bolivia used
BOL_GISAID_lambda <- Bolivia_sequences[Bolivia_sequences$Variant=="Lambda",]

# Sample Brazil, Lambda => All Lambda sequences from Bolivia used
BRA_GISAID_lambda <- Brazil_sequences[Brazil_sequences$Variant=="Lambda",]

# Sample Colombia, Lambda, ISP => All Lambda sequences from Colombia included to match expected number
COL_GISAID_lambda <- Colombia_sequences[Colombia_sequences$Variant=="Lambda",]

# Sample Peru, Lambda, ISP
sample <- Peru_sequences[Peru_sequences$Variant=="Lambda",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_lambda_GISAID[CL_migration_lambda_GISAID$origin_admin_0=="Peru",2])/wk_lambda)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_GISAID_lambda <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Lambda",] |>
                                  group_by(epiweek) |>
                                  slice(sample(min(max, n()))))


## Mu
# Sample Argentina, Mu, ISP => All Mu sequences from Argentina included to match expected number
ARG_ISP_mu <- Argentina_sequences[Argentina_sequences$Variant=="Mu",]

# Sample Bolivia, Mu => All Mu sequences from Bolivia used
BOL_ISP_mu <- Bolivia_sequences[Bolivia_sequences$Variant=="Mu",]

# Sample Brazil, Mu => All Mu sequences from Brazil used
BRA_ISP_mu <- Brazil_sequences[Brazil_sequences$Variant=="Mu",]

# Sample Colombia, Mu, ISP
sample <- Colombia_sequences[Colombia_sequences$Variant=="Mu",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_mu_ISP[CL_migration_mu_ISP$origin_admin_0=="Colombia",2])/wk_mu)
sample$count <- ifelse(sample$count>max, max, sample$count)

COL_ISP_mu <- as.data.frame(Colombia_sequences[Colombia_sequences$Variant=="Mu",] |>
                              group_by(epiweek) |>
                              slice(sample(min(max, n()))))

# Sample Peru, Mu, ISP
sample <- Peru_sequences[Peru_sequences$Variant=="Mu",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_mu_ISP[CL_migration_mu_ISP$origin_admin_0=="Peru",2])/wk_mu)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_ISP_mu <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Mu",] |>
                                  group_by(epiweek) |>
                                  slice(sample(min(max, n()))))

# Sample Argentina, Mu => All Mu sequences from Argentina included to match expected number
ARG_GISAID_mu <- Argentina_sequences[Argentina_sequences$Variant=="Mu",]

# Sample Bolivia, Mu => All Mu sequences from Bolivia used
BOL_GISAID_mu <- Bolivia_sequences[Bolivia_sequences$Variant=="Mu",]

# Sample Brazil, Mu => All Mu sequences from Brazil used
BRA_GISAID_mu <- Brazil_sequences[Brazil_sequences$Variant=="Mu",]

# Sample Colombia, Mu, ISP
sample <- Colombia_sequences[Colombia_sequences$Variant=="Mu",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_mu_GISAID[CL_migration_mu_GISAID$origin_admin_0=="Colombia",2])/wk_mu)
sample$count <- ifelse(sample$count>max, max, sample$count)

COL_GISAID_mu <- as.data.frame(Colombia_sequences[Colombia_sequences$Variant=="Mu",] |>
                              group_by(epiweek) |>
                              slice(sample(min(max, n()))))

# Sample Peru, Mu, ISP
sample <- Peru_sequences[Peru_sequences$Variant=="Mu",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_mu_GISAID[CL_migration_mu_GISAID$origin_admin_0=="Peru",2])/wk_mu)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_GISAID_mu <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Mu",] |>
                              group_by(epiweek) |>
                              slice(sample(min(max, n()))))

## Delta
# Sample Argentina, Delta, ISP
sample <- Argentina_sequences[Argentina_sequences$Variant=="Delta",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_delta_ISP[CL_migration_delta_ISP$origin_admin_0=="Argentina",2])/wk_delta)
sample$count <- ifelse(sample$count>max, max, sample$count)

ARG_ISP_delta <- as.data.frame(Argentina_sequences[Argentina_sequences$Variant=="Delta",] |>
                                  group_by(epiweek) |>
                                  slice(sample(min(max, n()))))

# Sample Bolivia, Delta => All Delta sequences from Bolivia used
BOL_ISP_delta <- Bolivia_sequences[Bolivia_sequences$Variant=="Delta",]

# Sample Brazil, Delta => All Delta sequences from Brazil used
sample <- Brazil_sequences[Brazil_sequences$Variant=="Delta",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_delta_ISP[CL_migration_delta_ISP$origin_admin_0=="Brazil",2])/wk_delta)
sample$count <- ifelse(sample$count>max, max, sample$count)

BRA_ISP_delta <- as.data.frame(Brazil_sequences[Brazil_sequences$Variant=="Delta",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Colombia, Delta => All Delta sequences from Colombia included to match expected number
sample <- Colombia_sequences[Colombia_sequences$Variant=="Delta",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_delta_ISP[CL_migration_delta_ISP$origin_admin_0=="Colombia",2])/wk_delta)
sample$count <- ifelse(sample$count>max, max, sample$count)

COL_ISP_delta <- as.data.frame(Colombia_sequences[Colombia_sequences$Variant=="Delta",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Peru, Delta, ISP
sample <- Peru_sequences[Peru_sequences$Variant=="Delta",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_delta_ISP[CL_migration_delta_ISP$origin_admin_0=="Peru",2])/wk_delta)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_ISP_delta <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Delta",] |>
                                  group_by(epiweek) |>
                                  slice(sample(min(max, n()))))

# Sample Argentina, Delta, GISAID => All Delta sequences from Argentina included to match expected number
ARG_GISAID_delta <- Argentina_sequences[Argentina_sequences$Variant=="Delta",]

# Sample Bolivia, Delta => All Delta sequences from Bolivia used
BOL_GISAID_delta <- Bolivia_sequences[Bolivia_sequences$Variant=="Delta",]

# Sample Brazil, Delta, GISAID
sample <- Brazil_sequences[Brazil_sequences$Variant=="Delta",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_delta_GISAID[CL_migration_delta_GISAID$origin_admin_0=="Brazil",2])/wk_delta)
sample$count <- ifelse(sample$count>max, max, sample$count)

BRA_GISAID_delta <- as.data.frame(Brazil_sequences[Brazil_sequences$Variant=="Delta",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Colombia, Delta, GISAID
sample <- Colombia_sequences[Colombia_sequences$Variant=="Delta",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_delta_GISAID[CL_migration_delta_GISAID$origin_admin_0=="Colombia",2])/wk_delta)
sample$count <- ifelse(sample$count>max, max, sample$count)

COL_GISAID_delta <- as.data.frame(Colombia_sequences[Colombia_sequences$Variant=="Delta",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))

# Sample Peru, Delta, GISAID
sample <- Peru_sequences[Peru_sequences$Variant=="Delta",] |>
  group_by(epiweek) |>
  summarise(count = n()) |>
  as.data.frame()
max <- ceiling(as.numeric(CL_migration_delta_GISAID[CL_migration_delta_GISAID$origin_admin_0=="Peru",2])/wk_delta)
sample$count <- ifelse(sample$count>max, max, sample$count)

PER_GISAID_delta <- as.data.frame(Peru_sequences[Peru_sequences$Variant=="Delta",] |>
                                 group_by(epiweek) |>
                                 slice(sample(min(max, n()))))


## Sampling seven countries with direct flights to Chile
# Select countries
Mexico_sequences <- filter (GISAID_final, Country == 'Mexico')
Canada_sequences <- filter (GISAID_final, Country == 'Canada')
France_sequences <- filter (GISAID_final, Country == 'France')
Panama_sequences <- filter (GISAID_final, Country == 'Panama')
Spain_sequences <- filter (GISAID_final, Country == 'Spain')
USA_sequences <- filter (GISAID_final, Country == 'USA')
UK_sequences <- filter (GISAID_final, Country == 'United Kingdom')

# Add Variant column
Mexico_sequences <- Mexico_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Mexico_sequences$`Pango lineage`[
                                 grep("AY.", Mexico_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Mexico_sequences$`Pango lineage`[
                                 grep("P.1.", Mexico_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))


Canada_sequences <- Canada_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Canada_sequences$`Pango lineage`[
                                 grep("AY.", Canada_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Canada_sequences$`Pango lineage`[
                                 grep("P.1.", Canada_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

France_sequences <- France_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(France_sequences$`Pango lineage`[
                                 grep("AY.", France_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(France_sequences$`Pango lineage`[
                                 grep("P.1.", France_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

Panama_sequences <- Panama_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Panama_sequences$`Pango lineage`[
                                 grep("AY.", Panama_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Panama_sequences$`Pango lineage`[
                                 grep("P.1.", Panama_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

Spain_sequences <- Spain_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(Spain_sequences$`Pango lineage`[
                                 grep("AY.", Spain_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(Spain_sequences$`Pango lineage`[
                                 grep("P.1.", Spain_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

USA_sequences <- USA_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(USA_sequences$`Pango lineage`[
                                 grep("AY.", USA_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(USA_sequences$`Pango lineage`[
                                 grep("P.1.", USA_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))

UK_sequences <- UK_sequences |>
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(UK_sequences$`Pango lineage`[
                                 grep("AY.", UK_sequences$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(UK_sequences$`Pango lineage`[
                                 grep("P.1.", UK_sequences$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             TRUE ~ "Other"))


# Select a uniform sample of these countries over time
# Numbers of sequences per data set
alpha_ISP_ext <- ceiling((numseq_ISP_alpha/7)/wk_alpha)
gamma_ISP_ext <- ceiling((numseq_ISP_gamma/7)/wk_gamma)
lambda_ISP_ext <- ceiling((numseq_ISP_lambda/7)/wk_lambda)
mu_ISP_ext <- ceiling((numseq_ISP_mu/7)/wk_mu)
delta_ISP_ext <- ceiling((numseq_ISP_delta/7)/wk_delta)

alpha_GISAID_ext <- ceiling((numseq_GISAIDCL_alpha/7)/wk_alpha)
gamma_GISAID_ext <- ceiling((numseq_GISAIDCL_gamma/7)/wk_gamma)
lambda_GISAID_ext <- ceiling((numseq_GISAIDCL_lambda/7)/wk_lambda)
mu_GISAID_ext <- ceiling((numseq_GISAIDCL_mu/7)/wk_mu)
delta_GISAID_ext <- ceiling((numseq_GISAIDCL_delta/7)/wk_delta)

# Generate sampled data sets
Mexico_sequences_alpha_ISP <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_ISP_ext, n()))))
Mexico_sequences_gamma_ISP <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_ISP_ext, n()))))
Mexico_sequences_lambda_ISP <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Lambda",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(lambda_ISP_ext, n()))))
Mexico_sequences_mu_ISP <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Mu",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(mu_ISP_ext, n()))))
Mexico_sequences_delta_ISP <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_ISP_ext, n()))))

Mexico_sequences_alpha_GISAID <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_GISAID_ext, n()))))
Mexico_sequences_gamma_GISAID <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_GISAID_ext, n()))))
Mexico_sequences_lambda_GISAID <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Lambda",] |>
                                               group_by(epiweek) |>
                                               slice(sample(min(lambda_GISAID_ext, n()))))
Mexico_sequences_mu_GISAID <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Mu",] |>
                                           group_by(epiweek) |>
                                           slice(sample(min(mu_GISAID_ext, n()))))
Mexico_sequences_delta_GISAID <- as.data.frame(Mexico_sequences[Mexico_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_GISAID_ext, n()))))

Canada_sequences_alpha_ISP <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_ISP_ext, n()))))
Canada_sequences_gamma_ISP <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_ISP_ext, n()))))
Canada_sequences_lambda_ISP <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Lambda",] |>
                                               group_by(epiweek) |>
                                               slice(sample(min(lambda_ISP_ext, n()))))
Canada_sequences_mu_ISP <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Mu",] |>
                                           group_by(epiweek) |>
                                           slice(sample(min(mu_ISP_ext, n()))))
Canada_sequences_delta_ISP <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_ISP_ext, n()))))

Canada_sequences_alpha_GISAID <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Alpha",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(alpha_GISAID_ext, n()))))
Canada_sequences_gamma_GISAID <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Gamma",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(gamma_GISAID_ext, n()))))
Canada_sequences_lambda_GISAID <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Lambda",] |>
                                                  group_by(epiweek) |>
                                                  slice(sample(min(lambda_GISAID_ext, n()))))
Canada_sequences_mu_GISAID <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Mu",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(mu_GISAID_ext, n()))))
Canada_sequences_delta_GISAID <- as.data.frame(Canada_sequences[Canada_sequences$Variant=="Delta",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(delta_GISAID_ext, n()))))

France_sequences_alpha_ISP <- as.data.frame(France_sequences[France_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_ISP_ext, n()))))
France_sequences_gamma_ISP <- as.data.frame(France_sequences[France_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_ISP_ext, n()))))
France_sequences_lambda_ISP <- as.data.frame(France_sequences[France_sequences$Variant=="Lambda",] |>
                                               group_by(epiweek) |>
                                               slice(sample(min(lambda_ISP_ext, n()))))
France_sequences_mu_ISP <- as.data.frame(France_sequences[France_sequences$Variant=="Mu",] |>
                                           group_by(epiweek) |>
                                           slice(sample(min(mu_ISP_ext, n()))))
France_sequences_delta_ISP <- as.data.frame(France_sequences[France_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_ISP_ext, n()))))

France_sequences_alpha_GISAID <- as.data.frame(France_sequences[France_sequences$Variant=="Alpha",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(alpha_GISAID_ext, n()))))
France_sequences_gamma_GISAID <- as.data.frame(France_sequences[France_sequences$Variant=="Gamma",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(gamma_GISAID_ext, n()))))
France_sequences_lambda_GISAID <- as.data.frame(France_sequences[France_sequences$Variant=="Lambda",] |>
                                                  group_by(epiweek) |>
                                                  slice(sample(min(lambda_GISAID_ext, n()))))
France_sequences_mu_GISAID <- as.data.frame(France_sequences[France_sequences$Variant=="Mu",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(mu_GISAID_ext, n()))))
France_sequences_delta_GISAID <- as.data.frame(France_sequences[France_sequences$Variant=="Delta",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(delta_GISAID_ext, n()))))


Panama_sequences_alpha_ISP <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_ISP_ext, n()))))
Panama_sequences_gamma_ISP <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_ISP_ext, n()))))
Panama_sequences_lambda_ISP <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Lambda",] |>
                                               group_by(epiweek) |>
                                               slice(sample(min(lambda_ISP_ext, n()))))
Panama_sequences_mu_ISP <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Mu",] |>
                                           group_by(epiweek) |>
                                           slice(sample(min(mu_ISP_ext, n()))))
Panama_sequences_delta_ISP <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_ISP_ext, n()))))

Panama_sequences_alpha_GISAID <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Alpha",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(alpha_GISAID_ext, n()))))
Panama_sequences_gamma_GISAID <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Gamma",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(gamma_GISAID_ext, n()))))
Panama_sequences_lambda_GISAID <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Lambda",] |>
                                                  group_by(epiweek) |>
                                                  slice(sample(min(lambda_GISAID_ext, n()))))
Panama_sequences_mu_GISAID <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Mu",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(mu_GISAID_ext, n()))))
Panama_sequences_delta_GISAID <- as.data.frame(Panama_sequences[Panama_sequences$Variant=="Delta",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(delta_GISAID_ext, n()))))


Spain_sequences_alpha_ISP <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_ISP_ext, n()))))
Spain_sequences_gamma_ISP <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_ISP_ext, n()))))
Spain_sequences_lambda_ISP <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Lambda",] |>
                                               group_by(epiweek) |>
                                               slice(sample(min(lambda_ISP_ext, n()))))
Spain_sequences_mu_ISP <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Mu",] |>
                                           group_by(epiweek) |>
                                           slice(sample(min(mu_ISP_ext, n()))))
Spain_sequences_delta_ISP <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_ISP_ext, n()))))

Spain_sequences_alpha_GISAID <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Alpha",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(alpha_GISAID_ext, n()))))
Spain_sequences_gamma_GISAID <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Gamma",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(gamma_GISAID_ext, n()))))
Spain_sequences_lambda_GISAID <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Lambda",] |>
                                                  group_by(epiweek) |>
                                                  slice(sample(min(lambda_GISAID_ext, n()))))
Spain_sequences_mu_GISAID <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Mu",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(mu_GISAID_ext, n()))))
Spain_sequences_delta_GISAID <- as.data.frame(Spain_sequences[Spain_sequences$Variant=="Delta",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(delta_GISAID_ext, n()))))


USA_sequences_alpha_ISP <- as.data.frame(USA_sequences[USA_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_ISP_ext, n()))))
USA_sequences_gamma_ISP <- as.data.frame(USA_sequences[USA_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_ISP_ext, n()))))
USA_sequences_lambda_ISP <- as.data.frame(USA_sequences[USA_sequences$Variant=="Lambda",] |>
                                               group_by(epiweek) |>
                                               slice(sample(min(lambda_ISP_ext, n()))))
USA_sequences_mu_ISP <- as.data.frame(USA_sequences[USA_sequences$Variant=="Mu",] |>
                                           group_by(epiweek) |>
                                           slice(sample(min(mu_ISP_ext, n()))))
USA_sequences_delta_ISP <- as.data.frame(USA_sequences[USA_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_ISP_ext, n()))))

USA_sequences_alpha_GISAID <- as.data.frame(USA_sequences[USA_sequences$Variant=="Alpha",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(alpha_GISAID_ext, n()))))
USA_sequences_gamma_GISAID <- as.data.frame(USA_sequences[USA_sequences$Variant=="Gamma",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(gamma_GISAID_ext, n()))))
USA_sequences_lambda_GISAID <- as.data.frame(USA_sequences[USA_sequences$Variant=="Lambda",] |>
                                                  group_by(epiweek) |>
                                                  slice(sample(min(lambda_GISAID_ext, n()))))
USA_sequences_mu_GISAID <- as.data.frame(USA_sequences[USA_sequences$Variant=="Mu",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(mu_GISAID_ext, n()))))
USA_sequences_delta_GISAID <- as.data.frame(USA_sequences[USA_sequences$Variant=="Delta",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(delta_GISAID_ext, n()))))


UK_sequences_alpha_ISP <- as.data.frame(UK_sequences[UK_sequences$Variant=="Alpha",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(alpha_ISP_ext, n()))))
UK_sequences_gamma_ISP <- as.data.frame(UK_sequences[UK_sequences$Variant=="Gamma",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(gamma_ISP_ext, n()))))
UK_sequences_lambda_ISP <- as.data.frame(UK_sequences[UK_sequences$Variant=="Lambda",] |>
                                               group_by(epiweek) |>
                                               slice(sample(min(lambda_ISP_ext, n()))))
UK_sequences_mu_ISP <- as.data.frame(UK_sequences[UK_sequences$Variant=="Mu",] |>
                                           group_by(epiweek) |>
                                           slice(sample(min(mu_ISP_ext, n()))))
UK_sequences_delta_ISP <- as.data.frame(UK_sequences[UK_sequences$Variant=="Delta",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(delta_ISP_ext, n()))))

UK_sequences_alpha_GISAID <- as.data.frame(UK_sequences[UK_sequences$Variant=="Alpha",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(alpha_GISAID_ext, n()))))
UK_sequences_gamma_GISAID <- as.data.frame(UK_sequences[UK_sequences$Variant=="Gamma",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(gamma_GISAID_ext, n()))))
UK_sequences_lambda_GISAID <- as.data.frame(UK_sequences[UK_sequences$Variant=="Lambda",] |>
                                                  group_by(epiweek) |>
                                                  slice(sample(min(lambda_GISAID_ext, n()))))
UK_sequences_mu_GISAID <- as.data.frame(UK_sequences[UK_sequences$Variant=="Mu",] |>
                                              group_by(epiweek) |>
                                              slice(sample(min(mu_GISAID_ext, n()))))
UK_sequences_delta_GISAID <- as.data.frame(UK_sequences[UK_sequences$Variant=="Delta",] |>
                                                 group_by(epiweek) |>
                                                 slice(sample(min(delta_GISAID_ext, n()))))


### Merge all non-Chile sequences into variant-specific data sets
## ISP data
Alpha_ISP <- as.data.frame(rbind(ARG_ISP_alpha, BRA_ISP_alpha, COL_ISP_alpha, PER_ISP_alpha,
                                 Mexico_sequences_alpha_ISP, Canada_sequences_alpha_ISP,
                                 France_sequences_alpha_ISP, Panama_sequences_alpha_ISP,
                                 Spain_sequences_alpha_ISP, USA_sequences_alpha_ISP,
                                 UK_sequences_alpha_ISP))

Alpha_ISP$surveillance <- "background_sequences"

Gamma_ISP <- as.data.frame(rbind(ARG_ISP_gamma, BRA_ISP_gamma, COL_ISP_gamma, PER_ISP_gamma,
                                 Mexico_sequences_gamma_ISP, Canada_sequences_gamma_ISP,
                                 France_sequences_gamma_ISP, Panama_sequences_gamma_ISP,
                                 Spain_sequences_gamma_ISP, USA_sequences_gamma_ISP,
                                 UK_sequences_gamma_ISP))

Gamma_ISP$surveillance <- "background_sequences"

Lambda_ISP <- as.data.frame(rbind(ARG_ISP_lambda, BRA_ISP_lambda, COL_ISP_lambda, PER_ISP_lambda,
                                 Mexico_sequences_lambda_ISP, Canada_sequences_lambda_ISP,
                                 France_sequences_lambda_ISP, Panama_sequences_lambda_ISP,
                                 Spain_sequences_lambda_ISP, USA_sequences_lambda_ISP,
                                 UK_sequences_lambda_ISP))

Lambda_ISP$surveillance <- "background_sequences"

Mu_ISP <- as.data.frame(rbind(ARG_ISP_mu, BRA_ISP_mu, COL_ISP_mu, PER_ISP_mu,
                                  Mexico_sequences_mu_ISP, Canada_sequences_mu_ISP,
                                  France_sequences_mu_ISP, Panama_sequences_mu_ISP,
                                  Spain_sequences_mu_ISP, USA_sequences_mu_ISP,
                                  UK_sequences_mu_ISP))

Mu_ISP$surveillance <- "background_sequences"

Delta_ISP <- as.data.frame(rbind(ARG_ISP_delta, BRA_ISP_delta, COL_ISP_delta, PER_ISP_delta,
                                  Mexico_sequences_delta_ISP, Canada_sequences_delta_ISP,
                                  France_sequences_delta_ISP, Panama_sequences_delta_ISP,
                                  Spain_sequences_delta_ISP, USA_sequences_delta_ISP,
                                  UK_sequences_delta_ISP))

Delta_ISP$surveillance <- "background_sequences"

## GISAID data
Alpha_GISAID <- as.data.frame(rbind(ARG_GISAID_alpha, BRA_GISAID_alpha, COL_GISAID_alpha, PER_GISAID_alpha,
                                 Mexico_sequences_alpha_GISAID, Canada_sequences_alpha_GISAID,
                                 France_sequences_alpha_GISAID, Panama_sequences_alpha_GISAID,
                                 Spain_sequences_alpha_GISAID, USA_sequences_alpha_GISAID,
                                 UK_sequences_alpha_GISAID))

Alpha_GISAID$surveillance <- "background_sequences"

Gamma_GISAID <- as.data.frame(rbind(ARG_GISAID_gamma, BRA_GISAID_gamma, COL_GISAID_gamma, PER_GISAID_gamma,
                                 Mexico_sequences_gamma_GISAID, Canada_sequences_gamma_GISAID,
                                 France_sequences_gamma_GISAID, Panama_sequences_gamma_GISAID,
                                 Spain_sequences_gamma_GISAID, USA_sequences_gamma_GISAID,
                                 UK_sequences_gamma_GISAID))

Gamma_GISAID$surveillance <- "background_sequences"

Lambda_GISAID <- as.data.frame(rbind(ARG_GISAID_lambda, BRA_GISAID_lambda, COL_GISAID_lambda, PER_GISAID_lambda,
                                  Mexico_sequences_lambda_GISAID, Canada_sequences_lambda_GISAID,
                                  France_sequences_lambda_GISAID, Panama_sequences_lambda_GISAID,
                                  Spain_sequences_lambda_GISAID, USA_sequences_lambda_GISAID,
                                  UK_sequences_lambda_GISAID))

Lambda_GISAID$surveillance <- "background_sequences"

Mu_GISAID <- as.data.frame(rbind(ARG_GISAID_mu, BRA_GISAID_mu, COL_GISAID_mu, PER_GISAID_mu,
                              Mexico_sequences_mu_GISAID, Canada_sequences_mu_GISAID,
                              France_sequences_mu_GISAID, Panama_sequences_mu_GISAID,
                              Spain_sequences_mu_GISAID, USA_sequences_mu_GISAID,
                              UK_sequences_mu_GISAID))

Mu_GISAID$surveillance <- "background_sequences"

Delta_GISAID <- as.data.frame(rbind(ARG_GISAID_delta, BRA_GISAID_delta, COL_GISAID_delta, PER_GISAID_delta,
                                 Mexico_sequences_delta_GISAID, Canada_sequences_delta_GISAID,
                                 France_sequences_delta_GISAID, Panama_sequences_delta_GISAID,
                                 Spain_sequences_delta_GISAID, USA_sequences_delta_GISAID,
                                 UK_sequences_delta_GISAID))

Delta_GISAID$surveillance <- "background_sequences"

### Merge non-Chile and Chile metadata to create complete data sets
## ISP data
Meta_Alpha_ISP <- rbind(select(Alpha_ISP, colnames(chile_metadata_alpha)),
                   chile_metadata_alpha)
Meta_Gamma_ISP <- rbind(select(Gamma_ISP, colnames(chile_metadata_gamma)),
                        chile_metadata_gamma)
Meta_Lambda_ISP <- rbind(select(Lambda_ISP, colnames(chile_metadata_lambda)),
                        chile_metadata_lambda)
Meta_Mu_ISP <- rbind(select(Mu_ISP, colnames(chile_metadata_mu)),
                        chile_metadata_mu)
Meta_Delta_ISP <- rbind(select(Delta_ISP, colnames(chile_metadata_delta)),
                        chile_metadata_delta)

Meta_Alpha_GISAID <- rbind(select(Alpha_GISAID, colnames(chile_metadata_alpha)),
                        chile_metadata_alpha)
Meta_Gamma_GISAID <- rbind(select(Gamma_GISAID, colnames(chile_metadata_gamma)),
                        chile_metadata_gamma)
Meta_Lambda_GISAID <- rbind(select(Lambda_GISAID, colnames(chile_metadata_lambda)),
                         chile_metadata_lambda)
Meta_Mu_GISAID <- rbind(select(Mu_GISAID, colnames(chile_metadata_mu)),
                     chile_metadata_mu)
Meta_Delta_GISAID <- rbind(select(Delta_GISAID, colnames(chile_metadata_delta)),
                        chile_metadata_delta)

# Generate tip names to match phylogenetic trees
Meta_Alpha_ISP$`Virus name` <- str_c(Meta_Alpha_ISP$`Virus name`, "|",
                                     Meta_Alpha_ISP$`Accession ID`, "|", Meta_Alpha_ISP$Date)
Meta_Gamma_ISP$`Virus name` <-  str_c(Meta_Gamma_ISP$`Virus name`, "|",
                                      Meta_Gamma_ISP$`Accession ID`, "|", Meta_Gamma_ISP$Date)
Meta_Lambda_ISP$`Virus name` <-  str_c(Meta_Lambda_ISP$`Virus name`, "|",
                                       Meta_Lambda_ISP$`Accession ID`, "|", Meta_Lambda_ISP$Date)
Meta_Mu_ISP$`Virus name` <-  str_c(Meta_Mu_ISP$`Virus name`, "|",
                                   Meta_Mu_ISP$`Accession ID`, "|", Meta_Mu_ISP$Date)
Meta_Delta_ISP$`Virus name` <-  str_c(Meta_Delta_ISP$`Virus name`, "|",
                                      Meta_Delta_ISP$`Accession ID`, "|", Meta_Delta_ISP$Date)

Meta_Alpha_GISAID$`Virus name` <- str_c(Meta_Alpha_GISAID$`Virus name`, "|",
                                        Meta_Alpha_GISAID$`Accession ID`, "|", Meta_Alpha_GISAID$Date)
Meta_Gamma_GISAID$`Virus name` <-  str_c(Meta_Gamma_GISAID$`Virus name`, "|",
                                         Meta_Gamma_GISAID$`Accession ID`, "|", Meta_Gamma_GISAID$Date)
Meta_Lambda_GISAID$`Virus name` <-  str_c(Meta_Lambda_GISAID$`Virus name`, "|",
                                          Meta_Lambda_GISAID$`Accession ID`, "|", Meta_Lambda_GISAID$Date)
Meta_Mu_GISAID$`Virus name` <-  str_c(Meta_Mu_GISAID$`Virus name`, "|",
                                      Meta_Mu_GISAID$`Accession ID`, "|", Meta_Mu_GISAID$Date)
Meta_Delta_GISAID$`Virus name` <-  str_c(Meta_Delta_GISAID$`Virus name`, "|",
                                         Meta_Delta_GISAID$`Accession ID`, "|", Meta_Delta_GISAID$Date)


# Write metadata files
write.table(Meta_Alpha_ISP, file = "Data/SARS2_CL_Alpha_ISP.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Gamma_ISP, file = "Data/SARS2_CL_Gamma_ISP.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Lambda_ISP, file = "Data/SARS2_CL_Lambda_ISP.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Mu_ISP, file = "Data/SARS2_CL_Mu_ISP.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Delta_ISP, file = "Data/SARS2_CL_Delta_ISP.csv",
            row.names = FALSE, sep = ",", quote = FALSE)

write.table(Meta_Alpha_GISAID, file = "Data/SARS2_CL_Alpha_GISAID.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Gamma_GISAID, file = "Data/SARS2_CL_Gamma_GISAID.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Lambda_GISAID, file = "Data/SARS2_CL_Lambda_GISAID.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Mu_GISAID, file = "Data/SARS2_CL_Mu_GISAID.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
write.table(Meta_Delta_GISAID, file = "Data/SARS2_CL_Delta_GISAID.csv",
            row.names = FALSE, sep = ",", quote = FALSE)
