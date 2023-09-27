################################################################################
####################### SARS-CoV-2 Chile lockdowns #############################
################################################################################

# Source scripts and data
source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")

########################### Bernardo Gutierrez #################################

library(tidyverse)
library(janitor)
library(lubridate)
library(stringr)
library(seqinr)
library(stringi) # to remove accents from strings
library(scales)
library(zoo)
library(patchwork)
library(NatParksPalettes)

############################## Data wrangling ##################################
## Import case counts ####
comuna_cases <- read.table("Data/MSC_Epi/data_clean.csv", sep = ",") |>
  row_to_names(row_number = 1) |> drop_na(`Casos Confirmados`, Region) |>
  filter(`Casos Confirmados` != "-" & Region != "") |>
  mutate(Comuna = stri_trans_general(Comuna, "Latin-ASCII")) |>
  mutate(`Casos Confirmados` = as.numeric(`Casos Confirmados`)) |>
  mutate(Fecha = ymd(Fecha)) |>
  drop_na(`Casos Confirmados`) |>
  select(Fecha, `Casos Confirmados`, Comuna, Poblacion) |>
  group_by(Comuna) |>
  arrange(Fecha) |>
  mutate(`Casos nuevos` = `Casos Confirmados` -
           lag(`Casos Confirmados`, default = first(`Casos Confirmados`))) |>
  arrange(Comuna) |>
  filter(Fecha <= last_CL_seq & Fecha >= as.Date("2020-08-09")) |>
  mutate(`Casos nuevos` = case_when(`Casos nuevos` < 0 ~ 0,
                                    TRUE ~ `Casos nuevos`))

comuna_cases$Comuna[comuna_cases$Comuna == "Calera"] <- "Caldera"
comuna_cases$Comuna[comuna_cases$Comuna == "Coihaique"] <- "Coyhaique"
comuna_cases$Comuna[comuna_cases$Comuna == "Llaillay"] <- "Llay-Llay"
comuna_cases$Comuna[comuna_cases$Comuna == "OHiggins"] <- "O'Higgins"
comuna_cases$Comuna[comuna_cases$Comuna == "Paiguano"] <- "Paihuano"
comuna_cases$Comuna[comuna_cases$Comuna == "Tiltil"] <- "Til Til"
comuna_cases$Comuna[comuna_cases$Comuna == "Treguaco"] <- "Trehuaco"

## Import 'Paso a paso' lockdown data ####
## Obtained from the Ministry of Science data portal
lockdowns <- read.csv("Data/paso_a_paso_std.csv",
                    sep =",") |>
  select(region_residencia, comuna_residencia, Fecha, Paso) |>
  filter(Fecha <= last_CL_seq & Fecha >= as.Date("2020-08-09")) |>
  mutate(Fecha = as.Date(Fecha)) |>
  mutate(Paso = as.numeric(Paso)) |>
  mutate(comuna_residencia = stri_trans_general(comuna_residencia,
                                                "Latin-ASCII")) |>
  mutate(region_residencia = stri_trans_general(region_residencia,
                                                "Latin-ASCII"))
lockdowns$region_residencia <-
  factor(lockdowns$region_residencia,
         levels = c("Arica y Parinacota", "Tarapaca", "Antofagasta", "Atacama",
                    "Coquimbo", "Valparaiso", "Metropolitana", "O'Higgins",
                    "Maule", "Nuble", "Biobio", "La Araucania", "Los Rios",
                    "Los Lagos", "Aysen", "Magallanes"))

# Remove duplicates
# Some comunas have records for different tiers for the same date, this function
# only keeps the lower value for every comuna for every date
lockdowns <- lockdowns |>
  group_by(comuna_residencia, Fecha) |>
  filter(Paso == min(Paso)) |>
  ungroup()
lockdowns <- lockdowns[!duplicated(lockdowns),]

# Add column with custom lockdown tier
lockdowns <- lockdowns |>
  mutate(lockdown_tier = case_when(
    Paso == 1 ~ "Full lockdown",
    Paso == 2 ~ "Weekend lockdown",
    Paso == 3 | Paso == 4 ~ "No lockdown",
    .default = "Opening")) |>
  mutate(lockdown_tier = as.factor(lockdown_tier))

lockdowns$lockdown_tier <- factor(lockdowns$lockdown_tier,
  levels = c("Full lockdown", "Weekend lockdown", "No lockdown", "Opening"))

# Add column to lockdown df with comuna population
lockdowns <- lockdowns |>
  left_join(comuna_cases |> select(Comuna, Poblacion) |> unique(),
            by = c("comuna_residencia" = "Comuna")) |>
  mutate(Poblacion = as.double(Poblacion)) |>
  rename("population" = "Poblacion")

# Add column to lockdown df with new cases
lockdowns <- lockdowns |> rename("fecha" = "Fecha")
lockdowns <- lockdowns |>
  left_join(comuna_cases |> select(Comuna, Fecha, `Casos nuevos`),
            by = c("comuna_residencia" = "Comuna", "fecha" = "Fecha"))
lockdowns <- lockdowns |> rename("Fecha" = "fecha")
lockdowns <- lockdowns |> rename("casos_nuevos" = "Casos nuevos")


## Create df for time series of lockdowns in Chile
lockdowns_timeseries <- lockdowns |>
  filter(lockdown_tier == "Full lockdown") |>
  group_by(Fecha) |>
  summarise(lockdown_comunas = n(),
            lockdown_population = sum(population, na.rm = TRUE),
            cases = sum(casos_nuevos, na.rm = TRUE)) |>
  mutate(lockdown_comunas_percentage =
           lockdown_comunas / (length(unique(lockdowns$comuna_residencia)))*100) |>
  mutate(lockdown_population_percentage = 
           lockdown_population / (sum(unique(lockdowns$population), na.rm = TRUE))*100) |>
  mutate(cases_7day_avg = rollmean(cases, k = 7, fill = NA, align = 'right'))

cases_timeseries <- read.table("Data/owid-covid-data_20211205.csv", sep = ",") |>
  row_to_names(row_number = 1) |>
  select(location, date, new_cases, new_cases_smoothed, new_cases_per_million,
         new_cases_smoothed_per_million) |> filter(location == "Chile") |>
  select(-location) |>
  mutate(date = ymd(date))

cases_timeseries_red <- cases_timeseries |>
  filter(date >= min(chile_isp_metadata$epiweek_start) &
         date <= max(chile_isp_metadata$epiweek_start))

lockdowns_timeseries <- lockdowns_timeseries |>
  full_join(cases_timeseries_red |> select(date, new_cases_per_million,
                                        new_cases_smoothed_per_million),
            by = c("Fecha" = "date")) |>
  mutate(new_cases_per_million = as.double(new_cases_per_million)) |>
  mutate(new_cases_smoothed_per_million = as.double(new_cases_smoothed_per_million))

## Create df for lockdown durations in Chile
lockdowns_duration <- lockdowns |>
  dplyr::select(comuna_residencia, Fecha, lockdown_tier) |>
  group_by(comuna_residencia, lockdown_tier) |>
  arrange(Fecha) |>
  mutate(stringency_change = lockdown_tier !=
           lag(lockdown_tier, default = first(lockdown_tier)) |
           as.integer(Fecha - lag(Fecha, default = first(Fecha))) > 1,
         stringency_group = cumsum(stringency_change)) |>
  group_by(comuna_residencia, lockdown_tier, stringency_group) |>
  summarise(start_date = min(Fecha), end_date = max(Fecha), duration = n()) |>
  dplyr::select(-stringency_group) |>
  arrange(lockdown_tier, comuna_residencia, start_date)


################################# Plots ########################################
## Plot distribution of lockdown duration times for each stringency tier ####
lockdowns_duration |>
  filter(lockdown_tier != "Opening") |>
  ggplot() +
  geom_histogram(aes(x = round(duration / 7, digits = 0),
                     fill = lockdown_tier,
                     group = lockdown_tier),
                 bins = max(round(lockdowns_duration$duration / 7, digits = 0)),
                 alpha = 0.5) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_x_continuous(limits = c(0, 63), breaks = seq(0, 64, 2)) +
  theme_minimal() + labs(x = "Lockdown duration (weeks)", y = "") +
  theme(legend.position = "none", strip.text.x = element_text(size = 12)) +
  facet_wrap(vars(lockdown_tier), 3, 1)

ggsave("../Figures/CL_lockdown_durations.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")
ggsave("../Figures/CL_lockdown_durations.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")


## Plot case counts and population under full lockdown ####
ggplot(lockdowns_timeseries) +
  geom_col(aes(x = Fecha,
               y = lockdown_comunas_percentage),
           alpha = 0.6, fill = "orange") +
  geom_col(aes(x = Fecha,
               y = lockdown_population_percentage),
           alpha = 0.4, fill = "lightblue") +
  geom_line(aes(x = Fecha, y = new_cases_smoothed_per_million / 4), linewidth = 0.5,
            colour = "darkred") +
  geom_vline(xintercept = as.Date("2021-07-23"), linetype = "dashed", # 70% vaccination rate
             linewidth = 0.6, colour = "grey") + 
  geom_vline(xintercept = as.Date("2021-09-20"), linetype = "dashed", # 75% vaccination rate
             linewidth = 0.8, colour = "darkgrey") +
  scale_x_date(limits = c(min(chile_isp_metadata$epiweek_start),
                          max(chile_isp_metadata$epiweek_start))) +
  scale_y_continuous(name = "Percentage",
                     sec.axis = sec_axis(
                       trans = ~.*4,
                       name =
                  "New cases per million habitants (7-day rolling average)")) +
  theme_minimal() + labs(x = "") +
  theme(axis.title.y.right = element_text(color = "darkred"),
        axis.text.y.right = element_text(color = "darkred"))

ggsave("Figures/CL_cases_lockdowns.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")
ggsave("Figures/CL_cases_lockdowns.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")


ggplot(cases_timeseries) +
  geom_line(aes(x = date, y = as.double(new_cases_smoothed_per_million) / 4), linewidth = 0.5,
            colour = "darkred") +
  scale_y_continuous(name =
                         "New cases per million habitants (7-day rolling average)") +
  scale_x_date(limits = c(min(lockdowns$Fecha),
                          max(lockdowns$Fecha))) +
  theme_minimal() +
  theme(axis.title.y = element_text(color = "darkred"),
        axis.text.y = element_text(color = "darkred"))

ggsave("Figures/CL_cases.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")


## Plots for selected comunas over specific time intervals ####
## Alpha context
lockdowns |>
  filter(Fecha >= as.Date("2020-11-29") & Fecha <= as.Date("2021-05-02")) |>
  filter(comuna_residencia == "Antofagasta" |
           comuna_residencia == "La Florida" |
           comuna_residencia == "Maipu" |
           comuna_residencia == "Puente Alto" |
           comuna_residencia == "Santiago") |>
  ggplot(aes(x = Fecha, y = 6-Paso)) +
  geom_line(size = 2, colour = "lightblue") +
  facet_wrap(vars(comuna_residencia), nrow = 5, ncol = 1)

lockdowns |>
  filter(Fecha >= as.Date("2020-11-29") & Fecha <= as.Date("2021-05-02")) |>
  filter(comuna_residencia == "Antofagasta" |
           comuna_residencia == "La Florida" |
           comuna_residencia == "Maipu" |
           comuna_residencia == "Puente Alto" |
           comuna_residencia == "Santiago") |>
  ggplot(aes(x = Fecha, y = comuna_residencia, fill = lockdown_tier)) +
  geom_tile(alpha = 0.2) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  labs(x = "Date", y = "Comuna", fill = "Lockdown stringency") +
  theme_minimal()


## Chile context
lockdowns |>
  filter(region_residencia == "Metropolitana") |>
  ggplot(aes(x = Fecha, y = 6-Paso)) +
  geom_line(size = 1, colour = "lightblue") +
  labs(title = "Santiago Metropolitan", x = "Date", y = "Lockdown stringency") +
  facet_wrap(vars(comuna_residencia), ncol = 4)

lockdowns |>
  filter(region_residencia == "Metropolitana") |>
  ggplot(aes(x = Fecha, y = comuna_residencia, fill = lockdown_tier)) +
  geom_tile(alpha = 0.2) +
  scale_x_date(limits = c(as.Date("2020-09-28"), as.Date("2021-09-13"))) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  labs(x = "Date", y = "Comuna", fill = "Lockdown stringency") +
  theme_minimal()

lockdowns |>
  filter(comuna_residencia == "Santiago") |>
  ggplot(aes(x = Fecha, y = comuna_residencia, fill = lockdown_tier)) +
  geom_tile(alpha = 0.2) +
  scale_x_date(limits = c(as.Date("2020-09-28"), as.Date("2021-09-13"))) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  labs(x = "Date", y = "Comuna", fill = "Lockdown stringency") +
  theme_minimal()


lockdowns |>
  ggplot(aes(x = Fecha, y = comuna_residencia, fill = lockdown_tier)) +
  geom_tile(alpha = 0.8) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450", "#d3d3d3")) +
  labs(x = "Date", y = "Comuna", fill = "Lockdown stringency") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 4),
        legend.position = "top") +
  facet_grid(rows = vars(region_residencia), space = "free", scales = "free_y")

ggsave("Figures/CL_lockdown_tiers.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white") # Region letters to be edited


# Plot new cases in the same style as the lockdown tier plot
lockdowns |>
  filter(is.na(casos_nuevos) == FALSE) |>
  mutate(casos_nuevos = na.locf(casos_nuevos)) |>
  ggplot(aes(x = as.factor(Fecha), y = comuna_residencia, fill = log(casos_nuevos))) +
  geom_tile() +
  viridis::scale_fill_viridis() +
  labs(x = "Date", y = "Comuna", fill = "New cases (log10)") +
  theme_void() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 4),
        legend.position = "top") +
  facet_grid(rows = vars(region_residencia), space = "free", scales = "free_y")
