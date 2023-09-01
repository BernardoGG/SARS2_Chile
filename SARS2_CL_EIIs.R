################################################################################
######## SARS-CoV-2 Chile EII comparisons with importation dynamics ############
################################################################################

########################### Bernardo Gutierrez #################################

library(tidyverse)
library(patchwork)
library(NatParksPalettes)
library(scales)

########################## Weekly phylogeo counts ##############################
## Create data frame with numbers of sequences and introductions per week
# Time series for epiweek cumulative cases per data set
CL_imports <- bind_rows(CL_TLs, CL_singletons) |>
  select(-taxa, -seen, -seen_date, -epiyear_collection,
         -epiweek_collection, -epiweek_collection_start) |>
  mutate(epiweek_tmrca_start = as.Date(epiweek_tmrca_start)) |>
  mutate(epiweek_ptmrca_start = as.Date(epiweek_ptmrca_start))

a <- airport_metadata |>
  group_by(epiweek_start) |>
  summarise(n = n())
b <- community_metadata |>
  group_by(epiweek_start) |>
  summarise(n = n())
c <- CL_TLs |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
d <- CL_singletons |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
e <- CL_TLs[CL_TLs$airport == TRUE,] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
f <- CL_singletons[CL_singletons$airport == TRUE,] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
g <- CL_TLs |>
  group_by(epiweek_first_seen_start) |>
  summarise(n = n())
h <- CL_singletons |>
  group_by(epiweek_collection_start) |>
  summarise(n = n())
i <- CL_TLs[CL_TLs$airport == TRUE,] |>
  group_by(epiweek_first_seen_start) |>
  summarise(n = n())
j <- CL_singletons[CL_singletons$airport == TRUE,] |>
  group_by(epiweek_collection_start) |>
  summarise(n = n())
k <- CL_imports[CL_imports$variant == "Alpha",] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
l <- CL_imports[CL_imports$variant == "Alpha" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
m <- CL_imports[CL_imports$variant == "Gamma",] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
n <- CL_imports[CL_imports$variant == "Gamma" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
o <- CL_imports[CL_imports$variant == "Lambda",] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
p <- CL_imports[CL_imports$variant == "Lambda" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
q <- CL_imports[CL_imports$variant == "Mu",] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
r <- CL_imports[CL_imports$variant == "Mu" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
s <- CL_imports[CL_imports$variant == "Delta",] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
t <- CL_imports[CL_imports$variant == "Delta" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_tmrca_start) |>
  summarise(n = n())
u <- CL_TLs |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
v <- CL_singletons |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
w <- CL_TLs[CL_TLs$airport == TRUE,] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
x <- CL_singletons[CL_singletons$airport == TRUE,] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
y <- CL_imports[CL_imports$variant == "Alpha",] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
z <- CL_imports[CL_imports$variant == "Alpha" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
aa <- CL_imports[CL_imports$variant == "Gamma",] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
ab <- CL_imports[CL_imports$variant == "Gamma" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
ac <- CL_imports[CL_imports$variant == "Lambda",] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
ad <- CL_imports[CL_imports$variant == "Lambda" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
ae <- CL_imports[CL_imports$variant == "Mu",] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
af <- CL_imports[CL_imports$variant == "Mu" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
ag <- CL_imports[CL_imports$variant == "Delta",] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())
ah <- CL_imports[CL_imports$variant == "Delta" &
                  CL_imports$airport == TRUE,] |>
  group_by(epiweek_ptmrca_start) |>
  summarise(n = n())

colnames(a) <- c("epiweek", "n")
colnames(b) <- c("epiweek", "n")
colnames(c) <- c("epiweek", "n")
colnames(d) <- c("epiweek", "n")
colnames(e) <- c("epiweek", "n")
colnames(f) <- c("epiweek", "n")
colnames(g) <- c("epiweek", "n")
colnames(h) <- c("epiweek", "n")
colnames(i) <- c("epiweek", "n")
colnames(j) <- c("epiweek", "n")
colnames(k) <- c("epiweek", "n")
colnames(l) <- c("epiweek", "n")
colnames(m) <- c("epiweek", "n")
colnames(n) <- c("epiweek", "n")
colnames(o) <- c("epiweek", "n")
colnames(p) <- c("epiweek", "n")
colnames(q) <- c("epiweek", "n")
colnames(r) <- c("epiweek", "n")
colnames(s) <- c("epiweek", "n")
colnames(t) <- c("epiweek", "n")
colnames(u) <- c("epiweek", "n")
colnames(v) <- c("epiweek", "n")
colnames(w) <- c("epiweek", "n")
colnames(x) <- c("epiweek", "n")
colnames(y) <- c("epiweek", "n")
colnames(z) <- c("epiweek", "n")
colnames(aa) <- c("epiweek", "n")
colnames(ab) <- c("epiweek", "n")
colnames(ac) <- c("epiweek", "n")
colnames(ad) <- c("epiweek", "n")
colnames(ae) <- c("epiweek", "n")
colnames(af) <- c("epiweek", "n")
colnames(ag) <- c("epiweek", "n")
colnames(ah) <- c("epiweek", "n")

SARS2_CL_epiweeks <- full_join(a, b, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, c, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, d, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, e, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, f, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, g, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, h, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, i, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, j, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, k, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, l, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, m, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, n, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, o, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, p, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, q, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, r, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, s, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, t, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, u, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, v, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, w, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, x, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, y, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, z, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, aa, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, ab, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, ac, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, ad, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, ae, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, af, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, ag, by = "epiweek")
SARS2_CL_epiweeks <- full_join(SARS2_CL_epiweeks, ah, by = "epiweek")

colnames(SARS2_CL_epiweeks) <- c("epiweek_start", "airport_seqs",
                                 "community_seqs", "TLs_tmrca",
                                 "singletons_tmrca", "TLs_airport_tmrca",
                                 "singletons_airport_tmrca", "TLs_first_seen",
                                 "singletons_collection",
                                 "TLs_airport_first_seen",
                                 "singletons_airport_collection",
                                 "all_alpha", "airport_alpha",
                                 "all_gamma", "airport_gamma",
                                 "all_lambda", "airport_lambda",
                                 "all_mu", "airport_mu",
                                 "all_delta", "airport_delta",
                                 "TLs_ptmrca", "singletons_ptmrca",
                                 "TLs_airport_ptmrca",
                                 "singletons_airport_ptmrca",
                                 "all_alpha_ptmrca", "airport_alpha_ptmrca",
                                 "all_gamma_ptmrca", "airport_gamma_ptmrca",
                                 "all_lambda_ptmrca", "airport_lambda_ptmrca",
                                 "all_mu_ptmrca", "airport_mu_ptmrca",
                                 "all_delta_ptmrca", "airport_delta_ptmrca")

SARS2_CL_epiweeks$all_airport_tmrca <- SARS2_CL_epiweeks$TLs_airport_tmrca +
  SARS2_CL_epiweeks$singletons_airport_tmrca
SARS2_CL_epiweeks$all_airport_ptmrca <- SARS2_CL_epiweeks$TLs_airport_ptmrca +
  SARS2_CL_epiweeks$singletons_airport_ptmrca
SARS2_CL_epiweeks$all_airport_observed <- SARS2_CL_epiweeks$TLs_airport_first_seen +
  SARS2_CL_epiweeks$singletons_airport_collection

SARS2_CL_epiweeks <- SARS2_CL_epiweeks |> replace(is.na(.), 0)

SARS2_CL_epiweeks$all_imports <- SARS2_CL_epiweeks$TLs_tmrca +
  SARS2_CL_epiweeks$singletons_tmrca
SARS2_CL_epiweeks$all_imports_ptmrca <- SARS2_CL_epiweeks$TLs_ptmrca +
  SARS2_CL_epiweeks$singletons_ptmrca

SARS2_CL_epiweeks <- arrange(SARS2_CL_epiweeks, epiweek_start)

rm(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,ag,ah)


################### Import EII estimates and components ########################
#### Files generated by Joseph Tsui and Rosario Evans-Pena

## Relevant states for Brazil and the USA
br_states <- c("Rio de Janeiro", "Santa Catarina", "Sao Paulo", "Bahia",
               "Goias", "Parana", "Rio Grande do Sul")
us_states <- c("Florida", "Georgia", "New York", "Texas", "California")

### EIIs
## Import EIIs functions
read_weekly_EIIs <- function(x){
  read.csv(x, sep = ",") |>
    mutate(week = as.Date(week)) |>
    mutate(EII = as.numeric(EII))
}

read_monthly_EIIs <- function(x){
  read.csv(x, sep = ",") |>
    mutate(month = as.Date(month)) |>
    mutate(EII = as.numeric(EII))
}

read_weekly_l_EIIs <- function(x){
  read.csv(x, sep = ",") |>
    mutate(week = as.Date(week)) |>
    mutate(Alpha.EII = as.numeric(Alpha.EII)) |>
    mutate(Gamma.EII = as.numeric(Gamma.EII)) |>
    mutate(Lambda.EII = as.numeric(Lambda.EII)) |>
    mutate(Mu.EII = as.numeric(Mu.EII)) |>
    mutate(Delta.EII = as.numeric(Delta.EII)) |>
    mutate(Others.EII = as.numeric(Others.EII))
}

## Air-based EIIs
# Alpha
eii_alpha <- read_weekly_EIIs("EIIs_export/Alpha_weekly_EII_v2.csv") |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start) |>
  filter(Chile_airport == "SCL")

# Gamma
eii_gamma <- read_weekly_EIIs("EIIs_export/Gamma_weekly_EII_v2.csv") |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start) |>
  filter(Chile_airport == "SCL")

# Lambda
eii_lambda <- read_weekly_EIIs("EIIs_export/Lambda_weekly_EII_v2.csv") |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start) |>
  filter(Chile_airport == "SCL")

# Mu
eii_mu <- read_weekly_EIIs("EIIs_export/Mu_weekly_EII_v2.csv") |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start) |>
  filter(Chile_airport == "SCL")

# Delta
eii_delta <- read_weekly_EIIs("EIIs_export/Delta_weekly_EII_v2.csv") |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start) |>
  filter(Chile_airport == "SCL")

## Land-based EIIs
# Argentina
l_eii_arg <-
  read_weekly_l_EIIs(
    "EIIs_export/Argentina_Bolivia_Peru_land-based_EIIs/Argentina_weekly_border-crossing_EII.csv"
    ) |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start)

# Bolivia
l_eii_bol <-
  read_weekly_l_EIIs(
    "EIIs_export/Argentina_Bolivia_Peru_land-based_EIIs/Bolivia_weekly_border-crossing_EII.csv"
  ) |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start)

# Peru
l_eii_per <-
  read_weekly_l_EIIs(
    "EIIs_export/Argentina_Bolivia_Peru_land-based_EIIs/Peru_weekly_border-crossing_EII.csv"
  ) |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start)

# Alpha (land)
l_eii_alpha <- bind_rows(l_eii_arg |> select(week, Alpha.EII) |>
                           mutate(country = "Argentina"),
                         l_eii_bol |> select(week, Alpha.EII) |>
                           mutate(country = "Bolivia"),
                         l_eii_per |> select(week, Alpha.EII) |>
                           mutate(country = "Peru")) |>
  rename(EII = Alpha.EII) |> replace(is.na(.), 0)

# Gamma (land)
l_eii_gamma <- bind_rows(l_eii_arg |> select(week, Gamma.EII) |>
                           mutate(country = "Argentina"),
                         l_eii_bol |> select(week, Gamma.EII) |>
                           mutate(country = "Bolivia"),
                         l_eii_per |> select(week, Gamma.EII) |>
                           mutate(country = "Peru")) |>
  rename(EII = Gamma.EII) |> replace(is.na(.), 0)

# Lambda (land)
l_eii_lambda <- bind_rows(l_eii_arg |> select(week, Lambda.EII) |>
                           mutate(country = "Argentina"),
                         l_eii_bol |> select(week, Lambda.EII) |>
                           mutate(country = "Bolivia"),
                         l_eii_per |> select(week, Lambda.EII) |>
                           mutate(country = "Peru")) |>
  rename(EII = Lambda.EII) |> replace(is.na(.), 0)

# Mu (land)
l_eii_mu <- bind_rows(l_eii_arg |> select(week, Mu.EII) |>
                           mutate(country = "Argentina"),
                         l_eii_bol |> select(week, Mu.EII) |>
                           mutate(country = "Bolivia"),
                         l_eii_per |> select(week, Mu.EII) |>
                           mutate(country = "Peru")) |>
  rename(EII = Mu.EII) |> replace(is.na(.), 0)

# Delta (land)
l_eii_delta <- bind_rows(l_eii_arg |> select(week, Delta.EII) |>
                           mutate(country = "Argentina"),
                         l_eii_bol |> select(week, Delta.EII) |>
                           mutate(country = "Bolivia"),
                         l_eii_per |> select(week, Delta.EII) |>
                           mutate(country = "Peru")) |>
  rename(EII = Delta.EII) |> replace(is.na(.), 0)

### EII components
## Import total case counts
weekly_totalcasecounts <- bind_rows(
  bind_rows(
    read.csv(
      "EIIs_export/processed_data_export/OWID_world_weekly_new_cases.csv",
      sep = ",") |>
      select(-iso_code),
    read.csv(
      "EIIs_export/processed_data_export/USA_state_weekly_new_cases.csv",
      sep = ",") |>
      rename(location = State)
  ), read.csv(
    "EIIs_export/processed_data_export/Brazil_state_weekly_new_cases.csv",
    sep = ",") |>
    rename(location = State)
) |>
  mutate(week = as.Date(week)) |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start) |>
  arrange(location)

## Import variant case counts
weekly_voccasecounts <- bind_rows(
  bind_rows(
    read.csv(
      "EIIs_export/processed_data_export/OWID_world_weekly_VOC_counts_props.csv",
      sep = ",") |>
      rename(location = Country),
    read.csv(
      "EIIs_export/processed_data_export/USA_state_weekly_VOC_counts_props.csv",
      sep = ",") |>
      rename(location = State)
  ), read.csv(
    "EIIs_export/processed_data_export/Brazil_state_weekly_VOC_counts_props.csv",
    sep = ",") |>
    rename(location = State)
) |>
  mutate(week = as.Date(week)) |>
  filter(week %in% SARS2_CL_epiweeks$epiweek_start) |>
  arrange(location)

## Import population sizes
popsizes <- bind_rows(
  bind_rows(
    read.csv(
      "EIIs_export/processed_data_export/OWID_world_populations.csv",
      sep = ",") |>
      select(-iso3) |>
      rename(location = country),
    read.csv(
      "EIIs_export/processed_data_export/USA_state_populations.csv",
      sep = ",") |>
      rename(location = State)
  ), read.csv(
    "EIIs_export/processed_data_export/Brazil_state_populations.csv",
    sep = ",") |>
    rename(location = State)
) |>
  filter(location %in% weekly_totalcasecounts$location) |>
  arrange(location)

## Create sequencing intensity data frame
seq_intensity <- data.frame(
  week = weekly_voccasecounts$week,
  location = weekly_voccasecounts$location,
  genomes = weekly_voccasecounts$total.count,
  cases = weekly_totalcasecounts$new_cases,
  seq_proportion = (weekly_voccasecounts$total.count/
    weekly_totalcasecounts$new_cases)) |>
  mutate(Source = ifelse(location %in% br_states, "Brazil",
                         ifelse(location %in% us_states, "USA",
                                ifelse(location %in% countries, location,
                                       "Other"))))


################ Set data frames for comparisons by variant ####################
## EII vs phylogenetic importation data sets

# Create vectors with likely importation sources per variant
# These are defined as places that at any point in the time series show an EII
# equal to or greater than 1.
sources_alpha <- unique(eii_alpha$source[eii_alpha$EII >= 1])
sources_gamma <- unique(eii_gamma$source[eii_gamma$EII >= 1])
sources_lambda <- unique(eii_lambda$source[eii_lambda$EII >= 1])
sources_mu <- unique(eii_mu$source[eii_mu$EII >= 1])
sources_delta <- unique(eii_delta$source[eii_delta$EII >= 1])
sources_neighbours <- unique(neighbours$source)

# Create grouped data frames for each likely importation source and numbers of
# imports estimated from phylodynamic analyses

# Alpha
alpha_comp <- eii_alpha |>
  mutate(Source = ifelse(source %in% sources_alpha |
                           source %in% sources_neighbours,
                         country, "Other")) |>
  mutate(all_imports = rep(SARS2_CL_epiweeks$all_alpha,
                           length(unique(eii_alpha$source)))) |>
  mutate(airport_imports = rep(SARS2_CL_epiweeks$airport_alpha,
                               length(unique(eii_alpha$source)))) |>
  select(-Chile_airport)
alpha_comp <- left_join(alpha_comp, l_eii_alpha, by = c("week", "country")) |>
  rename(EII = EII.x) |>
  rename(l_EII = EII.y)

# Gamma
gamma_comp <- eii_gamma |>
  filter(Chile_airport == "SCL") |>
  mutate(Source = ifelse(source %in% sources_gamma |
                           source %in% sources_neighbours,
                         country, "Other")) |>
  mutate(all_imports = rep(SARS2_CL_epiweeks$all_gamma,
                           length(unique(eii_gamma$source)))) |>
  mutate(airport_imports = rep(SARS2_CL_epiweeks$airport_gamma,
                               length(unique(eii_gamma$source)))) |>
  select(-Chile_airport)
gamma_comp <- left_join(gamma_comp, l_eii_gamma, by = c("week", "country")) |>
  rename(EII = EII.x) |>
  rename(l_EII = EII.y)

# Lambda
lambda_comp <- eii_lambda |>
  filter(Chile_airport == "SCL") |>
  mutate(Source = ifelse(source %in% sources_lambda |
                           source %in% sources_neighbours, 
                         country, "Other")) |>
  mutate(all_imports = rep(SARS2_CL_epiweeks$all_lambda,
                           length(unique(eii_lambda$source)))) |>
  mutate(airport_imports = rep(SARS2_CL_epiweeks$airport_lambda,
                               length(unique(eii_lambda$source)))) |>
  select(-Chile_airport)
lambda_comp <- left_join(lambda_comp, l_eii_lambda, by = c("week", "country")) |>
  rename(EII = EII.x) |>
  rename(l_EII = EII.y)

# Mu
mu_comp <- eii_mu |>
  filter(Chile_airport == "SCL") |>
  mutate(Source = ifelse(source %in% sources_mu |
                           source %in% sources_neighbours,
                         country, "Other")) |>
  mutate(all_imports = rep(SARS2_CL_epiweeks$all_mu,
                           length(unique(eii_mu$source)))) |>
  mutate(airport_imports = rep(SARS2_CL_epiweeks$airport_mu,
                               length(unique(eii_mu$source)))) |>
  select(-Chile_airport)
mu_comp <- left_join(mu_comp, l_eii_mu, by = c("week", "country")) |>
  rename(EII = EII.x) |>
  rename(l_EII = EII.y)

# Delta
delta_comp <- eii_delta |>
  filter(Chile_airport == "SCL") |>
  mutate(Source = ifelse(source %in% sources_delta |
                           source %in% sources_neighbours,
                         country, "Other")) |>
  mutate(all_imports = rep(SARS2_CL_epiweeks$all_delta,
                           length(unique(eii_delta$source)))) |>
  mutate(airport_imports = rep(SARS2_CL_epiweeks$airport_delta,
                               length(unique(eii_delta$source)))) |>
  select(-Chile_airport)
delta_comp <- left_join(delta_comp, l_eii_delta, by = c("week", "country")) |>
  rename(EII = EII.x) |>
  rename(l_EII = EII.y)


#################### Granger causality tests per variant #######################
#### Alpha ####
alpha_florida <- alpha_comp |>
  filter(source == sources_alpha[1])
alpha_france <- alpha_comp |>
  filter(source == sources_alpha[2])
alpha_spain <- alpha_comp |>
  filter(source == sources_alpha[3])
alpha_arg <- alpha_comp |>
  filter(source == sources_neighbours[1])
alpha_bol <- alpha_comp |>
  filter(source == sources_neighbours[2])
alpha_per <- alpha_comp |>
  filter(source == sources_neighbours[3])
alpha_all <- alpha_comp |>
  group_by(week) |>
  summarise(all_imports = mean(all_imports),
            airport_imports = mean(airport_imports),
            EII = sum(EII),
            l_EII = sum(l_EII, na.rm = TRUE))

### Full Alpha, all imports
lag <- alpha_all |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = alpha_all)
lmtest::grangertest(EII ~ all_imports, order = n, data = alpha_all)

lag <- alpha_all |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = alpha_all)

lag <- alpha_all |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = alpha_all)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = alpha_all)

### Florida (US), all imports
lag <- alpha_florida |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = alpha_florida)
lmtest::grangertest(EII ~ all_imports, order = n, data = alpha_florida)

lag <- alpha_florida |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = alpha_florida)

### France, all imports
lag <- alpha_france |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = alpha_france)
lmtest::grangertest(EII ~ all_imports, order = n, data = alpha_france)

lag <- alpha_france |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = alpha_france)

### Spain, all imports
lag <- alpha_spain |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = alpha_spain)
lmtest::grangertest(EII ~ all_imports, order = n, data = alpha_spain)

lag <- alpha_spain |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = alpha_spain)

### Argentina, all imports
lag <- alpha_arg |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = alpha_arg)
lmtest::grangertest(EII ~ all_imports, order = n, data = alpha_arg)

lag <- alpha_arg |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = alpha_arg)

lag <- alpha_arg |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = alpha_arg)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = alpha_arg)

### Bolivia, all imports
ifelse(any((alpha_bol$EII == 0) == TRUE),
       noquote("No EII estimates for this country"),
       noquote("Proceed with test"))

### Peru, all imports
lag <- alpha_per |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = alpha_per)
lmtest::grangertest(EII ~ all_imports, order = n, data = alpha_per)

lag <- alpha_per |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = alpha_per)

lag <- alpha_per |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = alpha_per)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = alpha_per)


#### Gamma ####
gamma_rio <- gamma_comp |>
  filter(source == sources_gamma[1])
gamma_catarina <- gamma_comp |>
  filter(source == sources_gamma[2])
gamma_sp <- gamma_comp |>
  filter(source == sources_gamma[3])
gamma_col <- gamma_comp |>
  filter(source == sources_gamma[4])
gamma_arg <- gamma_comp |>
  filter(source == sources_neighbours[1])
gamma_bol <- gamma_comp |>
  filter(source == sources_neighbours[2])
gamma_per <- gamma_comp |>
  filter(source == sources_gamma[5])
gamma_all <- gamma_comp |>
  group_by(week) |>
  summarise(all_imports = mean(all_imports),
            airport_imports = mean(airport_imports),
            EII = sum(EII),
            l_EII = sum(l_EII, na.rm = TRUE))

### Full Gamma, all imports
lag <- gamma_all |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_all)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_all)

lag <- gamma_all |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_all)

lag <- gamma_all |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = gamma_all)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = gamma_all)

### Rio de Janeiro (BR), all imports
lag <- gamma_rio |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_rio)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_rio)

lag <- gamma_rio |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_rio)

### Santa Catarina (BR), all imports
lag <- gamma_catarina |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_catarina)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_catarina)

lag <- gamma_catarina |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_catarina)

### Sao Paulo (BR), all imports
lag <- gamma_sp |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_sp)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_sp)

lag <- gamma_sp |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_sp)

### Colombia, all imports
lag <- gamma_col |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_col)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_col)

lag <- gamma_col |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_col)

### Argentina, all imports
lag <- gamma_arg |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_arg)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_arg)

lag <- gamma_arg |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_arg)

lag <- gamma_arg |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = gamma_arg)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = gamma_arg)

### Bolivia, all imports
lag <- gamma_bol |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_bol)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_bol)

lag <- gamma_bol |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_bol)

lag <- gamma_bol |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = gamma_bol)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = gamma_bol)

### Peru, all imports
lag <- gamma_per |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = gamma_per)
lmtest::grangertest(EII ~ all_imports, order = n, data = gamma_per)

lag <- gamma_per |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = gamma_per)

lag <- gamma_per |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = gamma_per)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = gamma_per)


#### Lambda ####
lambda_arg <- lambda_comp |>
  filter(source == sources_neighbours[1])
lambda_bol <- lambda_comp |>
  filter(source == sources_neighbours[2])
lambda_per <- lambda_comp |>
  filter(source == sources_lambda[1])
lambda_all <- lambda_comp |>
  group_by(week) |>
  summarise(all_imports = mean(all_imports),
            airport_imports = mean(airport_imports),
            EII = sum(EII),
            l_EII = sum(l_EII, na.rm = TRUE))

### Full Lambda, all imports
lag <- lambda_all |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = lambda_all)
lmtest::grangertest(EII ~ all_imports, order = n, data = lambda_all)

lag <- lambda_all |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = lambda_all)

lag <- lambda_all |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = lambda_all)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = lambda_all)

### Argentina, all imports
lag <- lambda_arg |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = lambda_arg)
lmtest::grangertest(EII ~ all_imports, order = n, data = lambda_arg)

lag <- lambda_arg |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = lambda_arg)

lag <- lambda_arg |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = lambda_arg)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = lambda_arg)

### Bolivia, all imports
lag <- lambda_bol |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = lambda_bol)
lmtest::grangertest(EII ~ all_imports, order = n, data = lambda_bol)

lag <- lambda_bol |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = lambda_bol)

lag <- lambda_bol |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = lambda_bol)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = lambda_bol)

### Peru, all imports
lag <- lambda_per |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = lambda_per)
lmtest::grangertest(EII ~ all_imports, order = n, data = lambda_per)

lag <- lambda_per |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = lambda_per)

lag <- lambda_per |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = lambda_per)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = lambda_per)


#### Mu ####
mu_col <- mu_comp |>
  filter(source == sources_mu[1])
mu_arg <- mu_comp |>
  filter(source == sources_neighbours[1])
mu_bol <- mu_comp |>
  filter(source == sources_neighbours[2])
mu_per <- mu_comp |>
  filter(source == sources_neighbours[3])
mu_all <- mu_comp |>
  group_by(week) |>
  summarise(all_imports = mean(all_imports),
            airport_imports = mean(airport_imports),
            EII = sum(EII),
            l_EII = sum(l_EII, na.rm = TRUE))


### Full Mu, all imports
lag <- mu_all |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = mu_all)
lmtest::grangertest(EII ~ all_imports, order = n, data = mu_all)

lag <- mu_all |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = mu_all)

lag <- mu_all |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = mu_all)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = mu_all)

### Colombia, all imports
lag <- mu_col |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = mu_col)
lmtest::grangertest(EII ~ all_imports, order = n, data = mu_col)

lag <- mu_col |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = mu_col)

### Argentina, all imports
lag <- mu_arg |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = mu_arg)
lmtest::grangertest(EII ~ all_imports, order = n, data = mu_arg)

lag <- mu_arg |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = mu_arg)

lag <- mu_arg |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = mu_arg)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = mu_arg)

### Bolivia, all imports
lag <- mu_bol |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = mu_bol)
lmtest::grangertest(EII ~ all_imports, order = n, data = mu_bol)

lag <- mu_bol |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = mu_bol)

lag <- mu_bol |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = mu_bol)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = mu_bol)

### Peru, all imports
lag <- mu_per |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = mu_per)
lmtest::grangertest(EII ~ all_imports, order = n, data = mu_per)

lag <- mu_per |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = mu_per)

lag <- mu_per |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = mu_per)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = mu_per)


#### Delta ####
delta_sp <- delta_comp |>
  filter(source == sources_delta[1])
delta_fl <- delta_comp |>
  filter(source == sources_delta[2])
delta_ga <- delta_comp |>
  filter(source == sources_delta[3])
delta_ny <- delta_comp |>
  filter(source == sources_delta[4])
delta_tx <- delta_comp |>
  filter(source == sources_delta[5])
delta_fra <- delta_comp |>
  filter(source == sources_delta[6])
delta_nld <- delta_comp |>
  filter(source == sources_delta[7])
delta_esp <- delta_comp |>
  filter(source == sources_delta[8])
delta_arg <- delta_comp |>
  filter(source == sources_neighbours[1])
delta_bol <- delta_comp |>
  filter(source == sources_neighbours[2])
delta_per <- delta_comp |>
  filter(source == sources_neighbours[3])
delta_all <- delta_comp |>
  group_by(week) |>
  summarise(all_imports = mean(all_imports),
            airport_imports = mean(airport_imports),
            EII = sum(EII),
            l_EII = sum(l_EII, na.rm = TRUE))

### Full Delta, all imports
lag <- delta_all |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_all)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_all)

lag <- delta_all |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_all)

lag <- delta_all |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = delta_all)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = delta_all)

### Sao Paulo (BRA), all imports
lag <- delta_sp |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_sp)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_sp)

lag <- delta_sp |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_sp)

### Florida (US), all imports
lag <- delta_fl |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_fl)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_fl)

lag <- delta_fl |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_fl)

### Georgia (US), all imports
lag <- delta_ga |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_ga)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_ga)

lag <- delta_ga |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_ga)

### New York (US), all imports
lag <- delta_ny |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_ny)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_ny)

lag <- delta_ny |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_ny)

### Texas (US), all imports
lag <- delta_tx |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_tx)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_tx)

lag <- delta_tx |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_tx)

### France, all imports
lag <- delta_fra |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_fra)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_fra)

lag <- delta_fra |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_fra)

### Netherlands, all imports
lag <- delta_nld |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_nld)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_nld)

lag <- delta_nld |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_nld)

### Spain, all imports
lag <- delta_esp |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_esp)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_esp)

lag <- delta_esp |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_esp)

### Argentina, all imports
lag <- delta_arg |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_arg)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_arg)

lag <- delta_arg |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_arg)

lag <- delta_arg |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = delta_arg)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = delta_arg)

### Bolivia, all imports
lag <- delta_bol |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_bol)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_bol)

lag <- delta_bol |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_bol)

lag <- delta_bol |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = delta_bol)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = delta_bol)

### Peru, all imports
lag <- delta_per |>
  select(EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ EII, order = n, data = delta_per)
lmtest::grangertest(EII ~ all_imports, order = n, data = delta_per)

lag <- delta_per |>
  select(EII, all_imports) |>
  mutate(EII = c(0, diff(EII))) |>
  mutate(all_imports = c(0, diff(all_imports))) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(diff(all_imports) ~ diff(EII), order = n, data = delta_per)

lag <- delta_per |>
  select(l_EII, all_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(all_imports ~ l_EII, order = n, data = delta_per)
lmtest::grangertest(l_EII ~ all_imports, order = n, data = delta_per)


################################ Plots #########################################
# Colors for country sets
sources <- unique(c(unique(sources_alpha), unique(sources_gamma),
                    unique(sources_lambda), unique(sources_mu),
                    unique(sources_delta)))
countries <- unique(c(unique(alpha_comp$country), unique(gamma_comp$country), 
                      unique(lambda_comp$country), unique(mu_comp$country), 
                      unique(delta_comp$country)))
s_palette <- c(NatParksPalettes$BryceCanyon[1] |> unlist(), "lightgrey")
names(s_palette) <- c(countries, "Other")
s_palette_2 <- c(s_palette,
                 latam_palette[names(latam_palette) %in% unique(neighbours$source)])

# Alpha
a <- SARS2_CL_epiweeks |>
  select(epiweek_start, airport_alpha, all_alpha) |>
  pivot_longer(!epiweek_start, names_to = "Importation_route",
               values_to = "Count") |>
  as.data.frame()

imp_palette <- (NatParksPalettes$KingsCanyon[1] |> unlist())[3:4]
names(imp_palette) <- unique(a$Importation_route)

w <- ggplot(a, aes(x = epiweek_start, y = Count, fill = Importation_route)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = imp_palette) +
  labs(title = "Alpha", x = "", y = "New introductions\n(weekly)") +
  theme_minimal()

x <- ggplot() +
  geom_line(data = alpha_comp,
            aes(x = week,y = EII, group = source, color = Source)) +
  geom_point(data = alpha_comp,
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette) +
  geom_line(data = alpha_all,
            aes(x = week, y = EII/2), alpha = 0.2, linewidth = 2,
            color = "darkred") +
  scale_y_continuous(name = "EII by location",
                     sec.axis = sec_axis(trans = ~.*2,
                                         name = "EII global (faded red)")) +
  theme(legend.position = "top") +
  labs(x = "") + theme_minimal()

b <- weekly_voccasecounts$Alpha.prop[
  weekly_voccasecounts$location %in% unique(sources_alpha)]

y <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_alpha)) |>
  mutate(Source = ifelse(location != "Florida", location, "USA")) |>
  ggplot(aes(x = week,
             y = new_cases * b,
             group = location, color = Source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Estimated new\nAlpha cases") + theme_minimal()

z <- ggplot(seq_intensity[seq_intensity$location %in% sources_alpha,]) +
  geom_line(aes(x = week, y = seq_proportion, color = Source)) +
  xlim(c(as.Date(min(SARS2_CL_epiweeks$epiweek_start)),
         as.Date(max(SARS2_CL_epiweeks$epiweek_start)))) +
  scale_y_continuous(labels = percent) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Percentage of\nsequenced cases") + theme_minimal()

c <- w / x / y / z
ggsave("Figures/Alpha_imports_plot.png")

# Gamma
a <- SARS2_CL_epiweeks |>
  select(epiweek_start, airport_gamma, all_gamma) |>
  pivot_longer(!epiweek_start, names_to = "Importation_route",
               values_to = "Count") |>
  as.data.frame()

imp_palette <- (NatParksPalettes$KingsCanyon[1] |> unlist())[3:4]
names(imp_palette) <- unique(a$Importation_route)

w <- ggplot(a, aes(x = epiweek_start, y = Count, fill = Importation_route)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = imp_palette) +
  labs(title = "Gamma", x = "", y = "New introductions\n(weekly)") +
  theme_minimal()

x <- ggplot() +
  geom_line(data = gamma_comp,
            aes(x = week,y = EII, group = source, color = Source)) +
  geom_point(data = gamma_comp,
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette) +
  geom_line(data = gamma_all,
            aes(x = week, y = EII/2), alpha = 0.2, linewidth = 2,
            color = "darkred") +
  scale_y_continuous(name = "EII by location",
                     sec.axis = sec_axis(trans = ~.*2,
                                         name = "EII global (faded red)")) +
  theme(legend.position = "top") +
  labs(x = "") + theme_minimal()

b <- weekly_voccasecounts$Gamma.prop[
  weekly_voccasecounts$location %in% unique(sources_gamma)]

y <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_gamma)) |>
  mutate(Source = ifelse(location == "Rio de Janeiro" |
                           location == "Santa Catarina" |
                           location == "Sao Paulo",
                         "Brazil", location)) |>
  ggplot(aes(x = week,
             y = new_cases * b,
             group = location, color = Source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Estimated new\nGamma cases") + theme_minimal()

z <- ggplot(seq_intensity[seq_intensity$location %in% sources_gamma,]) +
  geom_line(aes(x = week, y = seq_proportion, group = location,
                colour = Source)) +
  xlim(c(as.Date(min(SARS2_CL_epiweeks$epiweek_start)),
         as.Date(max(SARS2_CL_epiweeks$epiweek_start)))) +
  scale_y_continuous(labels = percent) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Percentage of\nsequenced cases") + theme_minimal()

c <- w / x / y / z
ggsave("Figures/Gamma_imports_plot.png")

# Gamma(by pTMRCA)
a <- SARS2_CL_epiweeks |>
  select(epiweek_start, airport_gamma_ptmrca, all_gamma_ptmrca) |>
  pivot_longer(!epiweek_start, names_to = "Importation_route",
               values_to = "Count") |>
  as.data.frame()

imp_palette <- (NatParksPalettes$KingsCanyon[1] |> unlist())[3:4]
names(imp_palette) <- unique(a$Importation_route)

w <- ggplot(a, aes(x = epiweek_start, y = Count, fill = Importation_route)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = imp_palette) +
  labs(title = "Gamma (imports by pTMRCA)", x = "",
       y = "New introductions\n(weekly)") +
  theme_minimal()

x <- ggplot() +
  geom_line(data = gamma_comp,
            aes(x = week,y = EII, group = source, color = Source)) +
  geom_point(data = gamma_comp,
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette) +
  geom_line(data = gamma_all,
            aes(x = week, y = EII/2), alpha = 0.2, linewidth = 2,
            color = "darkred") +
  scale_y_continuous(name = "EII by location",
                     sec.axis = sec_axis(trans = ~.*2,
                                         name = "EII global (faded red)")) +
  theme(legend.position = "top") +
  labs(x = "") + theme_minimal()

b <- weekly_voccasecounts$Gamma.prop[
  weekly_voccasecounts$location %in% unique(sources_gamma)]

y <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_gamma)) |>
  mutate(Source = ifelse(location == "Rio de Janeiro" |
                           location == "Santa Catarina" |
                           location == "Sao Paulo",
                         "Brazil", location)) |>
  ggplot(aes(x = week,
             y = new_cases * b,
             group = location, color = Source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Estimated new\nGamma cases") + theme_minimal()

z <- ggplot(seq_intensity[seq_intensity$location %in% sources_gamma,]) +
  geom_line(aes(x = week, y = log(seq_proportion), group = location,
                colour = location)) +
  xlim(c(as.Date(min(SARS2_CL_epiweeks$epiweek_start)),
         as.Date(max(SARS2_CL_epiweeks$epiweek_start)))) +
  scale_y_continuous(labels = percent) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Percentage of\nsequenced cases") + theme_minimal()

c <- w / x / y / z
ggsave("Figures/Gamma_imports_pTMRCA_plot.png")

# Lambda
a <- SARS2_CL_epiweeks |>
  select(epiweek_start, airport_lambda, all_lambda) |>
  pivot_longer(!epiweek_start, names_to = "Importation_route",
               values_to = "Count") |>
  as.data.frame()

imp_palette <- (NatParksPalettes$KingsCanyon[1] |> unlist())[3:4]
names(imp_palette) <- unique(a$Importation_route)

w <- ggplot(a, aes(x = epiweek_start, y = Count, fill = Importation_route)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = imp_palette) +
  labs(title = "Lambda", x = "", y = "New introductions\n(weekly)") +
  theme_minimal()

x <- ggplot() +
  geom_line(data = lambda_comp,
            aes(x = week,y = EII, group = source, color = Source)) +
  geom_point(data = lambda_comp,
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette) +
  geom_line(data = lambda_all,
            aes(x = week, y = EII/2), alpha = 0.2, linewidth = 2,
            color = "darkred") +
  scale_y_continuous(name = "EII by location",
                     sec.axis = sec_axis(trans = ~.*2,
                                         name = "EII global (faded red)")) +
  theme(legend.position = "top") +
  labs(x = "") + theme_minimal()

b <- weekly_voccasecounts$Lambda.prop[
  weekly_voccasecounts$location %in% unique(sources_lambda)]

y <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_lambda)) |>
  mutate(Source = location) |>
  ggplot(aes(x = week,
             y = new_cases * b,
             group = location, color = Source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Estimated new\nLambda cases") + theme_minimal()

z <- ggplot(seq_intensity[seq_intensity$location %in% sources_lambda,]) +
  geom_line(aes(x = week, y = seq_proportion, group = location,
                colour = Source)) +
  xlim(c(as.Date(min(SARS2_CL_epiweeks$epiweek_start)),
         as.Date(max(SARS2_CL_epiweeks$epiweek_start)))) +
  scale_y_continuous(labels = percent) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Percentage of\nsequenced cases") + theme_minimal()

c <- w / x / y / z
ggsave("Figures/Lambda_imports_plot.png")

# Mu
a <- SARS2_CL_epiweeks |>
  select(epiweek_start, airport_mu, all_mu) |>
  pivot_longer(!epiweek_start, names_to = "Importation_route",
               values_to = "Count") |>
  as.data.frame()

imp_palette <- (NatParksPalettes$KingsCanyon[1] |> unlist())[3:4]
names(imp_palette) <- unique(a$Importation_route)

w <- ggplot(a, aes(x = epiweek_start, y = Count, fill = Importation_route)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = imp_palette) +
  labs(title = "Mu", x = "", y = "New introductions\n(weekly)") +
  theme_minimal()

x <- ggplot() +
  geom_line(data = mu_comp,
            aes(x = week,y = EII, group = source, color = Source)) +
  geom_point(data = mu_comp,
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette) +
  geom_line(data = mu_all,
            aes(x = week, y = EII/2), alpha = 0.2, linewidth = 2,
            color = "darkred") +
  scale_y_continuous(name = "EII by location",
                     sec.axis = sec_axis(trans = ~.*2,
                                         name = "EII global (faded red)")) +
  theme(legend.position = "top") +
  labs(x = "") + theme_minimal()

b <- weekly_voccasecounts$Mu.prop[
  weekly_voccasecounts$location %in% unique(sources_mu)]

y <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_mu)) |>
  mutate(Source = location) |>
  ggplot(aes(x = week,
             y = new_cases * b,
             group = location, color = Source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Estimated new\nMu cases") + theme_minimal()

z <- ggplot(seq_intensity[seq_intensity$location %in% sources_mu,]) +
  geom_line(aes(x = week, y = seq_proportion, group = location,
                colour = Source)) +
  xlim(c(as.Date(min(SARS2_CL_epiweeks$epiweek_start)),
         as.Date(max(SARS2_CL_epiweeks$epiweek_start)))) +
  scale_y_continuous(labels = percent) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Percentage of\nsequenced cases") + theme_minimal()

c <- w / x / y / z
ggsave("Figures/Mu_imports_plot.png")

# Delta
a <- SARS2_CL_epiweeks |>
  select(epiweek_start, airport_delta, all_delta) |>
  pivot_longer(!epiweek_start, names_to = "Importation_route",
               values_to = "Count") |>
  as.data.frame()

imp_palette <- (NatParksPalettes$KingsCanyon[1] |> unlist())[3:4]
names(imp_palette) <- unique(a$Importation_route)

w <- ggplot(a, aes(x = epiweek_start, y = Count, fill = Importation_route)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = imp_palette) +
  labs(title = "Mu", x = "", y = "New introductions\n(weekly)") +
  theme_minimal()

x <- ggplot() +
  geom_line(data = delta_comp,
            aes(x = week,y = EII, group = source, color = Source)) +
  geom_point(data = delta_comp,
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette) +
  geom_line(data = delta_all,
            aes(x = week, y = EII/2), alpha = 0.2, linewidth = 2,
            color = "darkred") +
  scale_y_continuous(name = "EII by location",
                     sec.axis = sec_axis(trans = ~.*2,
                                         name = "EII global (faded red)")) +
  theme(legend.position = "top") +
  labs(x = "") + theme_minimal()

b <- weekly_voccasecounts$Delta.prop[
  weekly_voccasecounts$location %in% unique(sources_delta)]

y <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_delta)) |>
    mutate(Source = ifelse(location == "Florida" |
                           location == "Georgia" |
                           location == "New York" |
                           location == "Texas", "USA",
                         ifelse(location == "Sao Paulo",
                                "Brazil", location))) |>
  ggplot(aes(x = week,
             y = new_cases * b,
             group = location, color = Source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Estimated new Delta cases") + theme_minimal()

z <- ggplot(seq_intensity[seq_intensity$location %in% sources_delta,]) +
  geom_line(aes(x = week, y = seq_proportion, group = location,
                colour = location)) +
  xlim(c(as.Date(min(SARS2_CL_epiweeks$epiweek_start)),
         as.Date(max(SARS2_CL_epiweeks$epiweek_start)))) +
  scale_y_continuous(labels = percent) +
  scale_colour_manual(values = s_palette) + 
  labs(x = "", y = "Percentage of\nsequenced cases") + theme_minimal()

c <- w / x / y / z
ggsave("Figures/Delta_imports_plot.png")


### Manuscript figure plots
## VOC EIIs
a <- ggplot() +
  geom_line(data = alpha_comp[alpha_comp$Source != "Other",],
            aes(x = week, y = EII, group = source, color = Source)) +
  geom_point(data = alpha_comp[alpha_comp$Source != "Other",],
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette_2[
    names(s_palette_2) %in% alpha_comp$Source[alpha_comp$Source != "Other"]]) +
  geom_line(data = alpha_comp,
            aes(x = week, y = l_EII, color = country),
            linetype = "dashed") +
  geom_point(data = alpha_comp,
             aes(x = week, y = l_EII, color = country)) +
  scale_y_continuous(name = element_blank()) +
  scale_x_date(breaks = "1 month", labels = NULL) +
  labs(title = "Alpha", x = "", color = NULL) + theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8))

b <- ggplot() +
  geom_line(data = gamma_comp[gamma_comp$Source != "Other",],
            aes(x = week, y = EII, group = source, color = Source)) +
  geom_point(data = gamma_comp[gamma_comp$Source != "Other",],
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values =s_palette_2[
    names(s_palette_2) %in% gamma_comp$Source[gamma_comp$Source != "Other"]]) +
  geom_line(data = gamma_comp,
            aes(x = week, y = l_EII, color = country),
            linetype = "dashed") +
  geom_point(data = gamma_comp,
             aes(x = week, y = l_EII, color = country)) +
  scale_y_continuous(name = element_blank()) +
  scale_x_date(breaks = "1 month", labels = NULL) +
  labs(title = "Gamma", x = "", color = NULL) + theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8))

c <- ggplot() +
  geom_line(data = lambda_comp[lambda_comp$Source != "Other",],
            aes(x = week, y = EII, group = source, color = Source)) +
  geom_point(data = lambda_comp[lambda_comp$Source != "Other",],
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette_2[
    names(s_palette_2) %in% lambda_comp$Source[lambda_comp$Source != "Other"]]) +
  geom_line(data = lambda_comp,
            aes(x = week, y = l_EII, color = country),
            linetype = "dashed") +
  geom_point(data = lambda_comp,
             aes(x = week, y = l_EII, color = country)) +
  scale_y_continuous(name = "EII by location") +
  scale_x_date(breaks = "1 month", labels = NULL) +
  labs(title = "Lambda", x = "", color = NULL) + theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 8))

d <- ggplot() +
  geom_line(data = mu_comp[mu_comp$Source != "Other",],
            aes(x = week, y = EII, group = source, color = Source)) +
  geom_point(data = mu_comp[mu_comp$Source != "Other",],
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette_2[
    names(s_palette_2) %in% mu_comp$Source[mu_comp$Source != "Other"]]) +
  geom_line(data = mu_comp,
            aes(x = week, y = l_EII, color = country),
            linetype = "dashed") +
  geom_point(data = mu_comp,
             aes(x = week, y = l_EII, color = country)) +
  scale_y_continuous(name = element_blank()) +
  scale_x_date(breaks = "1 month", labels = NULL) +
  labs(title = "Mu", x = "", color = NULL) + theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8))

e <- ggplot() +
  geom_line(data = delta_comp[delta_comp$Source != "Other",],
            aes(x = week, y = EII, group = source, color = Source)) +
  geom_point(data = delta_comp[delta_comp$Source != "Other",],
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_colour_manual(values = s_palette_2[
    names(s_palette_2) %in% delta_comp$Source[delta_comp$Source != "Other"]]) +
  geom_line(data = delta_comp,
            aes(x = week, y = l_EII, color = country),
            linetype = "dashed") +
  geom_point(data = delta_comp,
             aes(x = week, y = l_EII, color = country)) +
  scale_y_continuous(name = element_blank()) +
  scale_x_date(breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Delta", x = "", color = NULL) + theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8))

fig2b <- a / b / c / d / e
ggsave("Figures/all_EIIs_plot.png", height = 10,
       width = 7, dpi = 300, units = "in")

ggsave("Figures/all_EIIs_plot.pdf", height = 10,
       width = 7, dpi = 300, units = "in")

########################### Mixed effects model ################################
anova(lm(alpha_all$all_imports ~ alpha_all$EII))
anova(lm(alpha_florida$all_imports ~ alpha_florida$EII))
anova(lm(alpha_france$all_imports ~ alpha_france$EII))
anova(lm(alpha_spain$all_imports ~ alpha_spain$EII))


anova(lm(lambda_all$all_imports ~ lambda_all$EII))
anova(lm(lambda_per$all_imports ~ lambda_per$EII))

plot(lambda_per$all_imports, lambda_per$EII)
hist(lambda_per$EII)

################################# Sandbox ######################################

ggplot(SARS2_CL_epiweeks) +
  geom_line(aes(x = epiweek_start, y = airport_seqs), color = "red") +
  geom_col(aes(x = epiweek_start, y = all_airport_tmrca), fill = "blue") +
  theme_minimal() +
  xlim(min(CL_TLs$tmrca_date), max(CL_TLs$last_seen_date)) +
  labs(x = "Epidemiological week",
       y = "Count")

ggplot(SARS2_CL_epiweeks) +
  geom_line(aes(x = epiweek_start,
                y = all_airport_tmrca),
            color = "darkred",
            linewidth = 1.5) +
  geom_line(aes(x = epiweek_start,
                y = all_imports),
            color = "darkgreen",
            linewidth = 1.5) +
  theme_minimal()


ggplot(SARS2_CL_epiweeks) +
  geom_col(aes(x = epiweek_start, y = all_airport_tmrca/all_imports)) +
  theme_minimal()

ggplot(SARS2_CL_epiweeks) +
  geom_col(aes(x = epiweek_start, y = all_airport_observed/airport_seqs))




a <- SARS2_CL_epiweeks |>
  select(epiweek_start, airport_gamma, all_gamma) |>
  pivot_longer(!epiweek_start, names_to = "Importation_route",
               values_to = "Count") |>
  as.data.frame()

imp_palette <- (NatParksPalettes$KingsCanyon[1] |> unlist())[3:4]
names(imp_palette) <- unique(a$Importation_route)

y <- ggplot(a, aes(x = epiweek_start, y = Count, fill = Importation_route)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = imp_palette) +
  labs(x = "", y = "New introductions (weekly)") + theme_minimal()


b <- weekly_voccasecounts$Gamma.prop[
  weekly_voccasecounts$location %in% unique(sources_gamma) |
    weekly_voccasecounts$location == "Argentina"]

z <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_gamma) |
           location == "Argentina") |>
  mutate(Source = ifelse(location == "Rio de Janeiro" |
                           location == "Santa Catarina" |
                           location == "Sao Paulo",
                         "Brazil", location)) |>
  mutate(Transparency = ifelse(location == "Argentina", "Yes", "No")) |>
  ggplot() +
  geom_line(aes(x = week,
                y = new_cases * b,
                group = location,
                color = Source,
                alpha = Transparency),
            linewidth = 1) +
  scale_alpha_discrete(range = c(0.2, 1), guide = FALSE) +
  labs(x = "", y = "Estimated new Gamma cases") + theme_minimal()

y / z



y <- ggplot() +
  geom_line(data = gamma_comp_neigh,
            aes(x = week,y = EII, group = source, color = Source)) +
  geom_point(data = gamma_comp_neigh,
             aes(x = week, y = EII, group = source, color = Source)) +
  scale_y_continuous(name = "EII by location") +
  theme(legend.position = "top") +
  labs(title = "Gamma (with neighbouring countries)", x = "") + theme_minimal()

b <- weekly_voccasecounts$Gamma.prop[
  weekly_voccasecounts$location %in% unique(sources_gamma) |
    weekly_voccasecounts$location == "Argentina" |
    weekly_voccasecounts$location == "Bolivia"]

z <- weekly_totalcasecounts |>
  filter(location %in% unique(sources_gamma) |
           location == "Argentina" |
           location == "Bolivia") |>
  mutate(Source = ifelse(location == "Rio de Janeiro" |
                           location == "Santa Catarina" |
                           location == "Sao Paulo",
                         "Brazil", location)) |>
  ggplot(aes(x = week,
             y = new_cases * b,
             group = location, color = Source)) +
  geom_line(linewidth = 1) +
  labs(x = "", y = "Estimated new Gamma cases") + theme_minimal()

y / z



## Full Alpha, airport imports
lag <- alpha_all |>
  select(EII, airport_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(airport_imports ~ EII, order = n, data = alpha_all)
lmtest::grangertest(EII ~ airport_imports, order = n, data = alpha_all)

## Florida (US), airport imports
lag <- alpha_florida |>
  select(EII, airport_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(airport_imports ~ EII, order = n, data = alpha_florida)
lmtest::grangertest(EII ~ airport_imports, order = n, data = alpha_florida)

## France, airport imports
lag <- alpha_france |>
  select(EII, airport_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(airport_imports ~ EII, order = n, data = alpha_france)
lmtest::grangertest(EII ~ airport_imports, order = n, data = alpha_france)

## Spain, airport imports
lag <- alpha_spain |>
  select(EII, airport_imports) |>
  vars::VARselect(lag.max = 12, type = "both")
n <- min(lag$selection)

lmtest::grangertest(airport_imports ~ EII, order = n, data = alpha_spain)
lmtest::grangertest(EII ~ airport_imports, order = n, data = alpha_spain)

