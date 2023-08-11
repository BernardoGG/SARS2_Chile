################################################################################
############# SARS-CoV-2 Chile Estimated importation indices ###################
################################################################################

# Source scripts and data
source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")

########################### Bernardo Gutierrez #################################

library(tidyverse)
library(janitor)
library(lubridate)
library(stringr)
library(seqinr)
library(scales)
library(patchwork)
library(NatParksPalettes)

############################## Data wrangling ##################################
## Import flight data (no. of passengers per month) into Chile ####
## Obtained from the Ministry of Transportation of Chile
flights <- read.csv("Data/EII_data/CL_flights_Sep2020-Dec2021_incoming.csv",
                    sep =",")

## Create column for USA and Brazil Adm1 identifiers
flights$source <- ifelse(flights$Adm0 != "United States of America" &
                           flights$Adm0 != "Brazil",
                         flights$Adm0,
                         flights$Adm1)

## Create data sets for the Santiago de Chile airport and all other airports
flights_scl <- flights %>% filter(Airport == "Santiago") %>%
  select(source, Adm0, Sep_2020:Dec_2021) %>%
  gather(month, passenger_volume, Sep_2020:Dec_2021)

flights_scl$month <- factor(flights_scl$month,
                            levels = colnames(flights)[7:22])
flights_scl$date <- rep(seq(as.Date("2020/09/01"), by = "month",
                            length.out = length(unique(flights_scl$month))),
                        each = nrow(flights_scl[
                          flights_scl$month == "Nov_2020",]))
flights_scl <- flights_scl %>% replace(is.na(.), 0)

flights_other <- flights %>% filter(Airport != "Santiago") %>%
  select(source, Adm0, Sep_2020:Dec_2021) %>%
  gather(month, passenger_volume, Sep_2020:Dec_2021)
flights_other$month <- factor(flights_other$month,
                            levels = colnames(flights)[7:22])
flights_other$date <- rep(seq(as.Date("2020/09/01"), by = "month",
                            length.out = length(unique(flights_other$month))),
                        each = nrow(flights_other[
                          flights_other$month == "Nov_2020",]))
flights_other <- flights_other %>% replace(is.na(.), 0)

flights_cl <- bind_rows(flights_scl, flights_other)
flights_cl$airport <- factor(c(rep("SCL", nrow(flights_scl),),
                        rep("Other", nrow(flights_other))),
                        levels = c("SCL", "Other"))


## Import land border crossings (no. of people per month) into Chile ####
# Obtained from the Information Transparency Agency of Chile
# Request placed by Leo Ferres
crossings <- read.csv("Data/EII_data/crossings.csv",
                    sep =",") %>%
  mutate(value =  as.numeric(gsub(",","",value))) %>%
  mutate(commune = stringi::stri_trans_general(commune, "Latin-ASCII")) %>%
  mutate(region = stringi::stri_trans_general(region, "Latin-ASCII")) %>%
  mutate(region = ifelse(region == "Arica y Parinacota (XV)",
                         "Arica y Parinacota",
                         region)) %>%
  mutate(region = ifelse(region == "Antofagasta (II)",
                         "Antofagasta",
                         region)) %>%
  mutate(date = case_when(month == "enero" ~ as.Date("2021-01-01"),
                          month == "febrero" ~ as.Date("2021-02-01"),
                          month == "marzo" ~ as.Date("2021-03-01"),
                          month == "abril" ~ as.Date("2021-04-01"),
                          month == "mayo" ~ as.Date("2021-05-01"),
                          month == "junio" ~ as.Date("2021-06-01"),
                          month == " julio" ~ as.Date("2021-07-01"),
                          month == "agosto" ~ as.Date("2021-08-01"),
                          month == "septiembre" ~ as.Date("2021-09-01"),
                          month == "octubre" ~ as.Date("2021-10-01"),
                          month == "noviembre" ~ as.Date("2021-11-01"),
                          month == "diciembre" ~ as.Date("2021-12-01"))) %>%
  filter(type == "entry") %>%
  filter(by == "land") %>%
  select(-type, -by, -year, -month)

# Create column with cross-border country
for(i in 1:nrow(crossings)){
  crossings$border_with[i] =
    if(crossings$region[i] == "Arica y Parinacota"){
      if(crossings$commune[i] == "Arica"){
        "Peru"
      } else {
        "Bolivia"
      }
    } else if(crossings$region[i] == "Tarapaca"){
      "Bolivia"
    } else if(crossings$region[i] == "Antofagasta"){
      if(crossings$commune[i] == "Ollague"){
        "Bolivia"
      } else {
        "Argentina"
      }
    } else {
      "Argentina"
    }
}

## Combined data set for land and air imports for neighboring countries ####
neighbours <- flights_cl[flights_cl$Adm0 == "Argentina" |
                           flights_cl$Adm0 == "Bolivia" |
                           flights_cl$Adm0 == "Peru",] %>%
  select(-Adm0) %>%
  group_by(source, date) %>%
  summarise(traveler_volume = sum(passenger_volume)) %>%
  ungroup() %>%
  mutate(port_of_entry = rep("airport", length(unique(flights_cl$date))*3)) %>%
  bind_rows(crossings %>%
              select(border_with, date, value) %>%
              group_by(border_with, date) %>%
              summarise(traveler_volume = sum(value)) %>%
              ungroup() %>%
              mutate(port_of_entry = rep("border crossing",
                                         length(unique(crossings$date))*3)) %>%
              rename(source = border_with))


############################## Summary plots ###################################
## Summary of number of international travelers by location
ggplot(flights) +
  geom_col(aes(x = Adm0, y = TOTAL, fill = Adm1)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Summaries of incoming flights into Santiago
# Flights by country
ggplot(flights_scl) +
  geom_col(aes(x = month, y = passenger_volume, fill = Adm0)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =
  "Numbers of passengers arrving to \nSantiago de Chile International Airport")


# South American flights by month
ggplot(flights_scl[flights_scl$Adm0 %in% latam,] %>%
         filter(month != "Sep_2020") %>%
         filter(month != "Oct_2020") %>%
         filter(month != "Nov_2020") %>%
         filter(month != "Nov_2021") %>%
         filter(month != "Dec_2021")) +
  geom_col(aes(x = month,
               y = passenger_volume,
               fill = Adm0),
           position = 'dodge') +
  theme_light() +
  scale_fill_manual(values = latam_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =
  "Numbers of passengers arrving to \nSantiago de Chile International Airport",
       x = "", y = "Monthly arrivals")

## Flights vs border crossings
# Combined
ggplot(neighbours) +
  geom_col(aes(x = date,
               y = traveler_volume,
               fill = source),
           position = 'dodge') +
  facet_wrap(~port_of_entry, ncol = 1) +
  theme_light() +
  scale_fill_manual(values = latam_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlim(as.Date("2020-12-15"), as.Date("2021-10-15")) +
  ylim(0, 60000) +
  labs(x = "", y = "Monthly incoming travelers \ninto Chile")

ggplot(neighbours) +
  geom_col(aes(x = date,
               y = traveler_volume,
               fill = port_of_entry)) +
  facet_wrap(~source, ncol = 1) +
  theme_light() +
  scale_fill_manual(values = c("#00AACC", "#8B5D33")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlim(as.Date("2020-12-15"), as.Date("2021-10-15")) +
  labs(x = "", y = "Monthly incoming travelers \ninto Chile")


# South American flights by month (2021) - stacked
a <- ggplot(flights_scl[flights_scl$Adm0 %in% latam,] %>%
              filter(month != "Sep_2020") %>%
              filter(month != "Oct_2020") %>%
              filter(month != "Nov_2020") %>%
              filter(month != "Dec_2020") %>%
              filter(month != "Nov_2021") %>%
              filter(month != "Dec_2021")) +
  geom_col(aes(x = date,
               y = passenger_volume,
               fill = Adm0)) +
  theme_light() +
  scale_fill_manual(values = latam_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlim(as.Date("2020-12-15"), as.Date("2021-10-15")) +
  ylim(c(0, 60000)) +
  labs(x = "",
       y = "Monthly arrivals to Santiago de \nChile International Airport")

# Border crossings by month - stacked
b <- ggplot(crossings) +
  geom_col(aes(x = date,
               y = value,
               fill = border_with)) +
  theme_light() +
  scale_fill_manual(values = latam_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlim(as.Date("2020-12-15"), as.Date("2021-10-15")) +
  ylim(c(0, 60000)) +
  labs(x = "", y = "Monthly incoming border \ncrossings into Chile")

a / b

ggsave("Figures/CL_intl_travel_stacked.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")


# South American flights by month (2021)
a <- ggplot(flights_scl[flights_scl$Adm0 %in% latam,] %>%
              filter(month != "Sep_2020") %>%
              filter(month != "Oct_2020") %>%
              filter(month != "Nov_2020") %>%
              filter(month != "Dec_2020") %>%
              filter(month != "Nov_2021") %>%
              filter(month != "Dec_2021")) +
  geom_col(aes(x = date,
               y = passenger_volume,
               fill = Adm0),
           position = 'dodge') +
  theme_light() +
  scale_fill_manual(values = latam_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlim(as.Date("2020-12-15"), as.Date("2021-10-15")) +
  ylim(c(0, 28000)) +
  labs(x = "",
       y = "Monthly arrivals to Santiago de \nChile International Airport")

# Border crossings by month
b <- ggplot(crossings) +
  geom_col(aes(x = date,
               y = value,
               fill = border_with),
           position = 'dodge',
           width = 10) +
  theme_light() +
  scale_fill_manual(values = latam_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlim(as.Date("2020-12-15"), as.Date("2021-10-15")) +
  ylim(c(0, 28000)) +
  labs(x = "", y = "Monthly incoming border \ncrossings into Chile")

a / b

ggsave("Figures/CL_intl_travel.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")



tmp <- aggregate(flights_scl$passenger_volume[flights_scl$Adm0 %in% latam],
                 by = list(flights_scl$month[flights_scl$Adm0 %in% latam]),
                 FUN = sum) %>% data.frame()

props <- flights_scl[flights_scl$Adm0 %in% latam,]
for(i in 1:nrow(props)){
  props$passenger_proportion[i] =
    props$passenger_volume[i] / tmp$x[tmp$Group.1 == props$month[i]]
}

ggplot(props) +
  geom_col(aes(x = month,
               y = passenger_proportion,
               fill = Adm0),
           position = 'dodge') +
  theme_light() +
  scale_fill_manual(values = latam_palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =
         "Proportion of passengers arrving to \nSantiago de Chile International Airport",
       x = "", y = "Monthly proportions")




################################# Sandbox ######################################
# Total flights (coloured by country)
ggplot(flights_scl) +
  geom_col(aes(x = month, y = passenger_volume)) +
  facet_wrap(~ Adm0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =
  "Numbers of passengers arrving to \nSantiago de Chile International Airport")



## Summaries of incoming flights into Chile (excl. Santiago)
# Flights by country
ggplot(flights_other) +
  geom_col(aes(x = month, y = passenger_volume, fill = Adm0)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =


# Total flights (coloured by country)
ggplot(flights_other) +
  geom_col(aes(x = month, y = passenger_volume)) +
  facet_wrap(~ Adm0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =
  "Numbers of passengers arrving to \nother airports in Chile")

ggplot(flights_cl) +
  geom_col(aes(x = month, y = passenger_volume, fill = airport)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title =
  "Numbers of passengers arrving to \ninternational airports in Chile")

## Summaries of incoming flights into Santiago
# Flights by country
ggplot(flights_other) +
  geom_col(aes(x = month, y = passenger_volume)) +
  facet_wrap(~ Adm0)

############################### Data export ####################################
## Processed flight data
write.csv(flights_cl, "Data/EII_data/CL_flights_Sep2020-Dec2021_processed.csv",
          row.names = FALSE)

## Processed combined flights and land movement data (neighboring countries)
write.csv(neighbours, "Data/EII_data/CL_border_crossings_2021_processed.csv",
          row.names = FALSE, quote = FALSE)
