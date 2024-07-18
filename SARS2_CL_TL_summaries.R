################################################################################
######### SARS-CoV-2 Chile transmission lineage phylogeo summaries #############
################################################################################

########################### Bernardo Gutierrez #################################

library(tidyverse)
library(patchwork)

######################### Import TL phylogeo table #############################
#### Files generated using SARS2_CL_merge_TL_phylogeo.R script

phylogeo_CL <- read.csv(
  "/Phylogenetics/SC2_CL_TLs/BEAST_TL_continuous_phylogeo/CL_transmission_lineages_phylogeo.csv",
  sep = ",")

##################### Summaries for onward transmission ########################
CL_lineages <- bind_rows(CL_TLs |> 
                           mutate(singleton = rep("No", nrow(CL_TLs))) |>
                           select(variant, ntaxa, airport, singleton),
                         CL_singletons |>
                           mutate(ntaxa = rep(1, nrow(CL_singletons)),
                                  singleton = rep("Yes", nrow(CL_singletons))) |>
                           select(variant, ntaxa, airport, singleton))

ggplot(CL_lineages |> filter(singleton == "No")) +
  geom_density(aes(x = ntaxa, group = variant)) +
  facet_wrap(vars(variant)) +
  theme_minimal()

CL_lineages_summary <- full_join(
  full_join(CL_lineages |> group_by(variant) |>
    summarize(n.intro = n(), n.taxa = sum(ntaxa)),
    CL_lineages |> filter(ntaxa != 1) |> group_by(variant) |>
      summarize(n.intro.onwards.transmission = n(), n.desc.circ.lin = sum(ntaxa),
                median.desc.circ.lin = median(ntaxa), mean.desc.circ.lin = mean(ntaxa),
                max.desc.circ.lin = max(ntaxa)) |>
      mutate(prop.desc.circ = n.intro.onwards.transmission / n.desc.circ.lin)) |>
    mutate(n.intro.hpd.high = c(43,344,83,124,92),
           n.intro.hpd.low = c(35,327,64,95,79),
           prop.intro.onwards = n.intro.onwards.transmission / n.intro,
           prop.intro.onwards.hpd.high = n.intro.onwards.transmission / n.intro.hpd.high,
           prop.intro.onwards.hpd.low = n.intro.onwards.transmission / n.intro.hpd.low),
  voc_cases_cl |> filter(Variant != "Other") |> rename("variant" = "Variant") |>
    group_by(variant) |> summarise(total.cases.voc = sum(cases_voc))) |>
  mutate(variant = factor(
    variant, levels = c("Alpha", "Gamma", "Lambda", "Mu", "Delta")))

upper <- ggplot(CL_lineages_summary) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", size = 0.8, color = "grey") +
  geom_linerange(aes(x = variant, ymin = prop.intro.onwards.hpd.low,
                     ymax = prop.intro.onwards.hpd.high, colour = variant), linewidth = 0.5) +
  geom_point(aes(x = variant, y = prop.intro.onwards,
                 size = n.intro, color = variant), alpha = 0.3) +
  geom_point(aes(x = variant, y = prop.intro.onwards,
                 color = variant)) +
  scale_size(range = c(10, 25)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = vocs_colors) +
  labs(x = element_blank(),
       y = "Proportion of viral introductions\nthat show onward transmission",
       size = "Number of viral introductions") +
  guides(size = guide_legend(override.aes = list(shape = 21,
                                                 fill = "transparent",
                                                 color = "darkgrey")),
         color = "none") + theme_minimal()

CL_lineages_timeline <- left_join(
  voc_cases_cl |> filter(Variant != "Other") |>
    rename("variant" = "Variant") |>
    select(-percentage, -new_cases_smoothed),
  CL_TLs |> group_by(variant, epiweek_tmrca_start) |>
    rename("epiweek_start" = "epiweek_tmrca_start") |>
    summarise(n.intro = n(), mean.n.desc = mean(ntaxa),
              median.n.desc = median(ntaxa)))

lower <- ggplot(CL_lineages_timeline) +
  geom_smooth(aes(x = median.n.desc, y = cases_voc, color = variant), method = "lm",
              formula = y ~ x, alpha = 0.05, size = 0.1, se = FALSE) +
  #  geom_smooth(aes(x = mean.n.desc, y = cases_voc), method = "lm",
  #              formula = y ~ x, se = FALSE, size = 1, color = "lightgrey") +
  geom_point(aes(x = median.n.desc, y = cases_voc,
                 size = n.intro, color = variant), alpha = 0.3) +
  scale_size(range = c(3, 7)) +
  scale_colour_manual(values = vocs_colors) +
  labs(y = "Total weekly cases\nestimated per variant",
       x = "Median number of descendants for viral introductions\nthat shows onward transmission",
       color = "Variant",
       size = "Number of viral introductions\nthat show onward transmission") +
  guides(size = guide_legend(
    override.aes = list(shape = 21, fill = "transparent", color = "darkgrey"))) +
  theme_minimal()

upper + lower + plot_layout(heights = c(6, 10))

ggsave("Figures/Chile_VOC_introductions_onward_transmission.png",
       dpi = 300, height = 10, width = 10, bg = "white")

###################### First detection of TLs in Chile #########################
## Create vector with EPI_ISL code fr earliest sequence(s) per TL
## Note that some TLs have multiple sequences detected on the earliest date for
## that TL; in such cases, all of them are included.

x <- vector()
for(i in 1:nrow(CL_fertree)){
  x[i] <- str_extract_all(CL_fertree$taxa[i],
                          regex("((19|2[0-9])[0-9]{2})-(0[1-9]|1[012])-(0[1-9]|[12][0-9]|3[01])")) |>
    unlist() |> as.Date() |> min() |> as.character()
}

first_taxa <- NULL
for(i in 1:nrow(CL_fertree)){
  y <- CL_fertree$taxa[i] |>
    str_extract_all(paste0("\\bEPI_ISL_\\d{7,8}\\b(?=[^>|]*?\\|", x[i], ")"),
                    simplify = T) |> as.vector()
  first_taxa <- c(y, first_taxa)
}


## Extract provinces and comunas where each VOC was first detected following
## importation.

first_comuna_alpha <- data.frame(
  Province = chile_metadata_alpha$State[
    chile_metadata_alpha$`Accession ID` %in% first_taxa],
  Comuna = chile_metadata_alpha$comuna[
    chile_metadata_alpha$`Accession ID` %in% first_taxa] |>
    str_to_title())

first_comuna_gamma <- data.frame(
  Province = chile_metadata_gamma$State[
    chile_metadata_gamma$`Accession ID` %in% first_taxa],
  Comuna = chile_metadata_gamma$comuna[
    chile_metadata_gamma$`Accession ID` %in% first_taxa] |>
    str_to_title())

first_comuna_lambda <- data.frame(
  Province = chile_metadata_lambda$State[
    chile_metadata_lambda$`Accession ID` %in% first_taxa],
  Comuna = chile_metadata_lambda$comuna[
    chile_metadata_lambda$`Accession ID` %in% first_taxa] |>
    str_to_title())

first_comuna_mu <- data.frame(
  Province = chile_metadata_mu$State[
    chile_metadata_mu$`Accession ID` %in% first_taxa],
  Comuna = chile_metadata_mu$comuna[
    chile_metadata_mu$`Accession ID` %in% first_taxa] |>
    str_to_title())

first_comuna_delta <- data.frame(
  Province = chile_metadata_delta$State[
    chile_metadata_delta$`Accession ID` %in% first_taxa],
  Comuna = chile_metadata_delta$comuna[
    chile_metadata_delta$`Accession ID` %in% first_taxa] |>
    str_to_title())

## Clean up comuna names
first_comuna_alpha$Comuna <- stri_trans_general(first_comuna_alpha$Comuna,
                                                "Latin-ASCII")
first_comuna_gamma$Comuna <- stri_trans_general(first_comuna_gamma$Comuna,
                                                "Latin-ASCII")
first_comuna_lambda$Comuna <- stri_trans_general(first_comuna_lambda$Comuna,
                                                "Latin-ASCII")
first_comuna_mu$Comuna <- stri_trans_general(first_comuna_mu$Comuna,
                                                "Latin-ASCII")
first_comuna_delta$Comuna <- stri_trans_general(first_comuna_delta$Comuna,
                                                "Latin-ASCII")

first_comuna_alpha$Comuna[
  grep("Con Con", first_comuna_alpha$Comuna)] <- "Concon"
first_comuna_alpha$Comuna[
  grep("Desconocido", first_comuna_alpha$Comuna)] <- NA

first_comuna_gamma$Comuna[
  grep("Coyhaique", first_comuna_gamma$Comuna)] <- "Coihaique"
first_comuna_gamma$Comuna[
  grep("Alto Hospicio", first_comuna_gamma$Comuna)] <- "Iquique"
first_comuna_gamma$Comuna[
  grep("Desconocido", first_comuna_gamma$Comuna)] <- NA

first_comuna_lambda$Comuna[
  grep("Coyhaique", first_comuna_lambda$Comuna)] <- "Coihaique"
first_comuna_lambda$Comuna[
  grep("Alto Hospicio", first_comuna_lambda$Comuna)] <- "Iquique"
first_comuna_lambda$Comuna[
  grep("Desconocido", first_comuna_lambda$Comuna)] <- NA

first_comuna_mu$Comuna[
  grep("Alto Hospicio", first_comuna_mu$Comuna)] <- "Iquique"
first_comuna_mu$Comuna[
  grep("Desconocido", first_comuna_mu$Comuna)] <- NA

first_comuna_delta$Comuna[
  grep("Con Con", first_comuna_delta$Comuna)] <- "Concon"
first_comuna_delta$Comuna[
  grep("Alto Hospicio", first_comuna_delta$Comuna)] <- "Iquique"
first_comuna_delta$Comuna[
  grep("Desconocido", first_comuna_delta$Comuna)] <- NA


######################### Mapping first detections #############################
# Read map
cl_map <- sf::st_read("Data/Comunas_de_Chile.geojson") |>
  sf::st_make_valid()

# Add centroids
cl_map <- cl_map  |>
  mutate(lon = map_dbl(geometry, ~sf::st_centroid(.x)[[1]]),
         lat = map_dbl(geometry, ~sf::st_centroid(.x)[[2]]))

# Comuna name cleanup
cl_map$comuna <- str_to_title(cl_map$comuna)
cl_map$comuna <- stringi::stri_trans_general(cl_map$comuna, "Latin-ASCII")
cl_map$comuna <- str_replace(cl_map$comuna, " De ", " de ")

# Add columns showing number of 'first detections' per comuna per VOC
alpha_seeds <- table(first_comuna_alpha$Comuna) |> as.data.frame() |>
  rename(comuna = Var1)
gamma_seeds <- table(first_comuna_gamma$Comuna) |> as.data.frame() |>
  rename(comuna = Var1)
lambda_seeds <- table(first_comuna_lambda$Comuna) |> as.data.frame() |>
  rename(comuna = Var1)
mu_seeds <- table(first_comuna_mu$Comuna) |> as.data.frame() |>
  rename(comuna = Var1)
delta_seeds <- table(first_comuna_delta$Comuna) |> as.data.frame() |>
  rename(comuna = Var1)

cl_map <- left_join(cl_map, alpha_seeds, by = "comuna") |>
  rename(alpha_seeds = Freq)
cl_map <- left_join(cl_map, gamma_seeds, by = "comuna") |>
  rename(gamma_seeds = Freq)
cl_map <- left_join(cl_map, lambda_seeds, by = "comuna") |>
  rename(lambda_seeds = Freq)
cl_map <- left_join(cl_map, mu_seeds, by = "comuna") |>
  rename(mu_seeds = Freq)
cl_map <- left_join(cl_map, delta_seeds, by = "comuna") |>
  rename(delta_seeds = Freq)

# Add marker column for border crossings
# Data from 'SARS2_CL_EIIs.R'
border_crossings <- crossings |>
  group_by(commune) |>
  summarise(total_crossings = sum(value)) |>
  rename(comuna = commune)

cl_map <- left_join(cl_map, border_crossings, by = "comuna")

# Add lat/long for border crossings
for(i in 1:nrow(cl_map)){
  cl_map$lon_bc[i] <- ifelse(is.na(cl_map$total_crossings[i]),
                             NA,
                             cl_map$lon[i])}

for(i in 1:nrow(cl_map)){
  cl_map$lat_bc[i] <- ifelse(is.na(cl_map$total_crossings[i]),
                             NA,
                             cl_map$lat[i])}

## Plot maps
a <- ggplot(cl_map) +
  geom_sf(aes(fill = alpha_seeds), colour = "#7f7f7f", linewidth = 0.06) +
  scale_fill_gradient(low = "#F7CCCA", high = "#A8201A",
                      na.value = "white") +
  geom_point(aes(x = lon_bc, y = lat_bc,
                 size = total_crossings/max(
                   total_crossings[is.na(
                     total_crossings) == FALSE])), alpha = 0.5) +
  theme_void() + labs(alpha = "Total border\ncrossings", fill = "Alpha",
                      size = NULL)

b <- ggplot(cl_map) +
  geom_sf(aes(fill = gamma_seeds), colour = "#7f7f7f", linewidth = 0.06) +
  scale_fill_gradient(low = "#F8DCB4", high = "#EC9A29",
                      na.value = "white") +
  geom_point(aes(x = lon_bc, y = lat_bc,
                 size = total_crossings/max(
                   total_crossings[is.na(
                     total_crossings) == FALSE])), alpha = 0.5) +
  theme_void() + labs(alpha = "Total border\ncrossings", fill = "Gamma",
                      size = NULL)

c <- ggplot(cl_map) +
  geom_sf(aes(fill = lambda_seeds), colour = "#7f7f7f", linewidth = 0.06) +
  scale_fill_gradient(low = "#C8F8F9", high = "#0F8B8D",
                      na.value = "white") +
  geom_point(aes(x = lon_bc, y = lat_bc,
                 size = total_crossings/max(
                   total_crossings[is.na(
                     total_crossings) == FALSE])), alpha = 0.5) +
  theme_void() + labs(alpha = "Total border\ncrossings", fill = "Lambda",
                      size = NULL)

d <- ggplot(cl_map) +
  geom_sf(aes(fill = mu_seeds), colour = "#7f7f7f", linewidth = 0.06) +
  scale_fill_gradient(low = "#D0E8F1", high = "#143642",
                      na.value = "white") +
  geom_point(aes(x = lon_bc, y = lat_bc,
                 size = total_crossings/max(
                   total_crossings[is.na(
                     total_crossings) == FALSE])), alpha = 0.5) +
  theme_void() + labs(alpha = "Total border\ncrossings", fill = "Mu",
                      size = NULL)

e <- ggplot(cl_map) +
  geom_sf(aes(fill = delta_seeds), colour = "#7f7f7f", linewidth = 0.06) +
  scale_fill_gradient(low = "#F6DACB", high = "#CA5D22",
                      na.value = "white") +
  geom_point(aes(x = lon_bc, y = lat_bc,
                 size = total_crossings/max(
                   total_crossings[is.na(
                     total_crossings) == FALSE])), alpha = 0.5) +
  theme_void() + labs(alpha = "Total border\ncrossings", fill = "Delta",
                      size = NULL)

f <- a | b | c | d | e

ggsave("Figures/CL_first detections_plot.png", dpi = 300,
       height = 7, width = 10, units = "in", bg = "white")
ggsave("Figures/CL_first detections_plot.pdf", dpi = 300,
       height = 7, width = 10, units = "in", bg = "white")
rm(a,b,c,d,e,f)
