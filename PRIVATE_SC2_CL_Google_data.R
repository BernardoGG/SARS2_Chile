################################################################################
################## ***** ANALYSES WITH PRIVATE DATA ***** ######################
################################################################################
################################################################################
################## SARS-CoV-2 Chile Google mobility trends #####################
################################################################################

########################### Bernardo Gutierrez #################################

library(tidyverse)

library(patchwork)
library(NatParksPalettes)
library(scales)

################################ Import data ###################################

# Google mobility data **PRIVATE DATA**
google_raw <- read.csv("Data/2022-08-05_Chile.csv",
                       sep = ",") %>%
  select(-X) %>%
  mutate(week_start = as.Date(week_start)) %>%
  mutate(week_end = as.Date(week_end)) %>%
  filter(week_start >= as.Date("2020-11-08") &
           week_start <= as.Date("2021-10-17"))

google_raw <- read.csv("Data/G_data_CL.csv",
                       sep = ",") %>%
  select(-X) %>%
  mutate(week_start = as.Date(week_start)) %>%
  mutate(week_end = as.Date(week_end)) %>%
  filter(week_start >= as.Date("2020-11-08") &
           week_start <= as.Date("2021-10-17"))

# Google lookup table
google_lookup <- read.csv("Data/CL_mobility_direct_flights_countries_lookup.csv",
                          sep = ",") %>% select(id, name)

# Add country names to Google data table
for(i in 1:nrow(google_raw)){
  google_raw$origin_location[i] <-
    google_lookup$name[google_lookup$id == google_raw$origin[i]]
}

for(i in 1:nrow(google_raw)){
  google_raw$destination_location[i] <-
    google_lookup$name[google_lookup$id == google_raw$destination[i]]
}

# Subset incoming movements into Chile
google_CL <- google_raw %>%
  filter(destination_location == "Chile") %>%
  select(-origin, -destination, -destination_location, -week_end)



ggplot(google_CL[google_CL$origin_location %in% latam,]) +
  geom_line(aes(x = week_start, y = movement, colour = origin_location)) +
  theme_minimal()
