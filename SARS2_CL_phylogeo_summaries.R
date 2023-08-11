################################################################################
###### SARS-CoV-2 Chile top transmission lineage phylogeographic analyses ######
################################################################################

########################### Bernardo Gutierrez #################################

library(tidyverse)
library(janitor)
library(lubridate)
library(stringi)
library(ggplot2)
library(scales)
library(vegan)
library(sf)
library(abjutils)
library(seraphim)
library(patchwork)
library(NatParksPalettes)

#################### Source phylogeography summaries ###########################

source('SARS2_CL_merge_TL_phylogeo.R')


######################## Utilities for analyses ################################

## Santiago Metropolitan comuna grouping
# Santiago Metropolitan Region (admin1)
santiago_reg <- c("Cerillos", "Cerro Navia", "Conchali", "El Bosque",
                  "Estacion Central", "Huechuraba", "Independencia",
                  "La Cisterna", "La Florida", "La Granja", "La Pintana",
                  "La Reina", "Las Condes", "Lo Barnechea", "Lo Espejo",
                  "Lo Prado", "Macul", "Maipu", "Nunoa", "Pedro Aguirre Cerda",
                  "Penalolen", "Providencia", "Pudahuel", "Quilicura",
                  "Quinta Normal", "Recoleta", "Renca", "San Joaquin",
                  "San Miguel", "San Ramon", "Santiago", "Vitacura",
                  "Puente Alto", "San Jose de Maipo", "Pirque", "Colina",
                  "Lampa", "Til Til", "Buin", "Calera de Tango", "Paine",
                  "San Bernardo", "Alhue", "Curacavi", "Maria Pinto",
                  "Melipilla", "San Pedro", "El Monte", "Isla de Maipo",
                  "Padre Hurtado", "Penaflor", "Talagante")

# Santiago Province - urban (admin2)
santiago_prov <- c("Cerillos", "Cerro Navia", "Conchali", "El Bosque",
                   "Estacion Central", "Huechuraba", "Independencia",
                   "La Cisterna", "La Florida", "La Granja", "La Pintana",
                   "La Reina", "Las Condes", "Lo Barnechea", "Lo Espejo",
                   "Lo Prado", "Macul", "Maipu", "Nunoa", "Pedro Aguirre Cerda",
                   "Penalolen", "Providencia", "Pudahuel", "Quilicura",
                   "Quinta Normal", "Recoleta", "Renca", "San Joaquin",
                   "San Miguel", "San Ramon", "Santiago", "Vitacura")

# Santiago metropolitan area (continuous urban area)
santiago_metro <- c("Cerrillos", "Cerro Navia", "Conchali", "El Bosque",
                    "Estacion Central", "Huechuraba", "Independencia",
                    "La Cisterna", "La Florida", "La Granja", "La Pintana",
                    "La Reina", "Las Condes", "Lo Barnechea", "Lo Espejo",
                    "Lo Prado", "Macul", "Maipu", "Nunoa", "Padre Hurtado",
                    "Pedro Aguirre Cerda", "Penalolen", "Providencia",
                    "Pudahuel", "Puente Alto", "Quilicura", "Quinta Normal",
                    "Recoleta", "Renca", "San Bernardo", "San Joaquin",
                    "San Miguel", "San Ramon", "Santiago", "Vitacura")

phylogeo_CL$head_comuna[phylogeo_CL$head_comuna=="San Jode de Maipo"] <- "San Jose de Maipo"
phylogeo_CL$tail_comuna[phylogeo_CL$tail_comuna=="San Jode de Maipo"] <- "San Jose de Maipo"
phylogeo_CL$head_comuna <- rm_accent(phylogeo_CL$head_comuna)
phylogeo_CL$tail_comuna <- rm_accent(phylogeo_CL$tail_comuna)

phylogeo_CL$head_santiago_reg <- replace(phylogeo_CL$head_comuna, phylogeo_CL$head_comuna %in% santiago_reg, "Santiago Metropolitan Region")
phylogeo_CL$head_santiago_prov <- replace(phylogeo_CL$head_comuna, phylogeo_CL$head_comuna %in% santiago_prov, "Santiago Province")
phylogeo_CL$head_santiago_metro <- replace(phylogeo_CL$head_comuna, phylogeo_CL$head_comuna %in% santiago_metro, "Santiago Metropolitan Area")
phylogeo_CL$tail_santiago_reg <- replace(phylogeo_CL$tail_comuna, phylogeo_CL$tail_comuna %in% santiago_reg, "Santiago Metropolitan Region")
phylogeo_CL$taiol_santiago_prov <- replace(phylogeo_CL$tail_comuna, phylogeo_CL$tail_comuna %in% santiago_prov, "Santiago Province")
phylogeo_CL$tail_santiago_metro <- replace(phylogeo_CL$tail_comuna, phylogeo_CL$tail_comuna %in% santiago_metro, "Santiago Metropolitan Area")


#################### Import and process Chile map ##############################

## Import map
map_CL <- rgdal::readOGR(dsn = "Data/gadm36_CHL_3/", layer = "gadm36_CHL_3")
map_CL <- world[world$name_long == "Chile", ]

map_CL <- sf::read_sf("Data/gadm36_CHL_3/") %>% filter(NAME_3 != "Ocean Islands")
tmap::tm_shape(map_CL) + tmap::tm_borders()

unique(map_CL_df$id)
unique(phylogeo_CL$head_comuna)

#ggplot() +
#  geom_polygon(data = map_CL, aes(x = long, y = lat, group = group), colour = "grey") +
#  theme_void() + theme(legend.position = "none")

################## Probability of exporting/importing ##########################

## Number of outgoing viral movements
source_count <- phylogeo_CL %>% count(head_santiago_metro, variant)

## Frequency of outgoing viral movements
sink_count <- phylogeo_CL %>% count(tail_santiago_metro, variant)

ggplot(source_count[source_count$variant=="Alpha",]) + geom_col(aes(head_santiago_metro, n))
ggplot(sink_count) + geom_col(aes(tail_santiago_metro, n, group = variant))


####################### Shannon Index per comuna ###############################

## Global Shannon index per comuna
count_internal <-
  phylogeo_CL[phylogeo_CL$head_comuna==phylogeo_CL$tail_comuna,] %>%
  count(invasion_tree, tail_comuna)

count_external <-
  phylogeo_CL[phylogeo_CL$head_comuna!=phylogeo_CL$tail_comuna,] %>%
  count(invasion_tree, tail_comuna)


data.frame(comuna = count_external %>%
             group_by(tail_comuna) %>%
             summarise(),
           shannon_index = count_external %>%
             group_by(tail_comuna) %>%
             group_map(~ diversity(.x$n)) %>%
             unlist())

data.frame(comuna = count_internal %>%
             group_by(tail_comuna) %>%
             summarise(),
           shannon_index = count_internal %>%
             group_by(tail_comuna) %>%
             group_map(~ diversity(.x$n)) %>%
             unlist())

## Weekly Shannon index per comuna
x <- estShannon(shannon_count$n) %>% unlist()
diversity(shannon_count$n)

##################### Spatial diffusion summaries ##############################
## Setting function parameters for all trees
localTreesDirectory <-
  "Phylogenetics/Transmission_lineage_processing/diffusion_statistics"
randomSampling = FALSE
coordinateAttribiuteName = "location"

## Create function to extract features of individual TLs


allTrees <- scan(file = "xx.trees", what = "", sep = "\n", quiet = TRUE)
burnin <- 0
nberOfTreesToSample = X
mostRecentSamplingDatum = X

