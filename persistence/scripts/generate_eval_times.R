library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(purrr)
library(patchwork)

last_sampling_dates <- read_tsv("./metadata/TL_last_sampling_date.tsv")
lockdown_dates <- read_tsv("./metadata/CL_lockdowns.tsv")

detect_change_points <- function(data_vector) {
  change_points <- c(T) # Include first date by default
  
  for (i in 2:(length(data_vector) -1)) {
    if (data_vector[i] != data_vector[i - 1]) {
      change_points <- c(change_points, T)
    } else {
      change_points <- c(change_points, F)
    }
  }
  
  change_points <- c(change_points, T) # Include last date also by default
  return(change_points)
}

change_point_dates <- lockdown_dates %>%
  filter(!is.na(lockdown_tier)) %>% # TODO: Check why there are NAs here
  group_by(comuna) %>%
  group_modify(~{
    .x %>%
      arrange(date) %>%
      mutate(
        change_points = detect_change_points(lockdown_tier)
      ) %>%
      filter(change_points)
  })

anc_eval_times <- change_point_dates %>%
  group_by(comuna) %>%
  group_modify(~{
    tmp_dates <- .x %>%
      arrange(date) %>%
      pull(date)
    map2_dfr(tmp_dates[1:length(tmp_dates)-1], tmp_dates[2:length(tmp_dates)],~{
      tmp_eval <- seq(.x, .y, by="1 week")
      if(!(.y %in% tmp_eval)) {
        tmp_eval <- c(tmp_eval, .y)
      }
      data.frame(
        evalTimes = tmp_eval[2:length(tmp_eval)], # First value will be ancTime
        ancTime = .x
      )
    })
  })

# Plot to confirm
select_locations <- c(
  "Pirque", 
  "Puente Alto", 
  "Santiago", 
  "San Jose de Maipo", 
  "San Bernardo", 
  "Arica",
  "Los Angeles", 
  "Requinoa", 
  "Copiapo", 
  "Iquique", 
  "Valparaiso", 
  "Antofagasta", 
  "Concepcion", 
  "Puerto Montt"
)

location_plots <- lockdown_dates %>%
  filter(comuna %in% select_locations) %>%
  group_by(comuna) %>%
  group_map(~{
    select_times <- anc_eval_times %>% filter(comuna == .y)
    .x %>% 
      ggplot(aes(date, stringency)) + 
        geom_line() + 
        theme_bw() +
        geom_vline(aes(xintercept = evalTimes), data = select_times, color="indianred", linetype="dashed", linewidth = 0.25) +
        geom_vline(aes(xintercept = ancTime), data = select_times, color="steelblue", linetype="dashed" , linewidth = 0.5) +
        ggtitle(.y)
    
  })

wrap_plots(location_plots)
ggsave("./plots/select_location_anc_eval_times.pdf", w = 15, h = 15)
  

#Convert to decimal dates
last_sampling_dates <- last_sampling_dates %>%
  mutate(
    last_sampling_date = decimal_date(last_sampling_date)
  )

anc_eval_times <- anc_eval_times %>%
  mutate(
    evalTimes = decimal_date(evalTimes),
    ancTime = decimal_date(ancTime)
  )

# Check all lineages with tree files
lineages <- list.files("./tree_files/", pattern="*_minimal$") %>% str_replace("_minimal", "")

select_last_sampling_dates <- last_sampling_dates %>%
  filter(transmission_lineage %in% lineages)

# For each lineage we need to run persistence for lock down dates of each comuna
# lsd <- select_last_sampling_dates %>% pull(last_sampling_date) %>% first()

anc_eval_times_df <- select_last_sampling_dates %>%
  group_by(transmission_lineage, last_sampling_date) %>%
  group_modify(~{
    lsd <- .y$last_sampling_date
    anc_eval_times %>% 
      filter(evalTimes < lsd) %>%
      mutate(
        evalTimes= lsd - evalTimes,
        ancTime = lsd - ancTime
      ) %>%
      group_by(comuna) %>%
      group_modify(~{
        data.frame(
          evalTimes = paste(.x$evalTimes, collapse=","),
          ancTimes = paste(.x$ancTime, collapse=",")  
        )
      })
  })

anc_eval_times_df %>% 
  filter(
    comuna %in% select_locations
  ) %>%
  mutate(
    comuna = str_replace_all(comuna, " ", "_")
  ) %>%
  group_by(transmission_lineage) %>%
  group_walk(~{
      lineage <- as.character(.y)
      print(lineage)
      .x %>%
        mutate(
          transmission_lineage = lineage
        ) %>% # Add linegae back in 
        relocate(transmission_lineage) %>% # Push to front for bash script
        write_tsv(paste0("./args_files/", lineage,"_args.tsv"))
  })
