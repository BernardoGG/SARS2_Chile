library(readr)
library(dplyr)
library(stringr)
library(ape)
library(purrr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tidyr)
library(coda)
library(khroma)

args <- commandArgs(trailingOnly = TRUE)

# location <- args[1]
lineage <- args[1]
# lineage <- "Delta_TL_10"

print(paste0("Running ", lineage))

out_folder <- paste0("./output/", lineage)

out_files <- list.files(out_folder, pattern="*.tsv", full.names=TRUE, recursive = TRUE)

df <- map_dfr(out_files, ~{
  treeIdfromFile <- str_extract(.x, "_[0-9]+.tsv") %>% str_replace_all("_|.tsv", "")
  location <- str_split(.x, "/") %>% .[[1]] %>% .[4]
  read_csv(.x, show_col_types=FALSE) %>%
    mutate(
      treeId = treeIdfromFile,
      location = location
    )
}) %>%
  mutate(
    stateAtEvaluationTime = stringi::stri_trans_general(stateAtEvaluationTime, "latin-ascii") # Remove non-ascii characters
  )

# Filter to branches that are in the state of interest
df_filtered <- df %>%
  filter(stateAtEvaluationTime == str_replace_all(location, "_", " "))

# Read in most recent sampling date
mrsd_df <- read_tsv("./metadata/TL_last_sampling_date.tsv")
mrsd <- mrsd_df %>%
  filter(transmission_lineage == lineage) %>%
  pull(last_sampling_date) %>%
  first() %>%
  decimal_date()

# Read in lockdown dates
lockdown_dates <- read_tsv("./metadata/CL_lockdowns.tsv")

# Filter to unique lineages
df_prop <- df_filtered %>%
  group_by(treeId, evaluationTime, ancestralTime, stateAtEvaluationTime) %>%
  summarise(
    anc_eval_diff = first(ancestralTime - evaluationTime), # evaluationtTime and ancestralTime should be same for all rows
    # Number of persistent lineages at eval time
    persistentLineagesAtEval = sum(persistenceTime > anc_eval_diff),
    # Number of lineages from new introductions at eval time
    introductionLineagesAtEval = sum(persistenceTime <= anc_eval_diff),
    # independenceTime used to filter to "unique" lineages
    persistentsFromUnique = sum(
      (independenceTime > anc_eval_diff) &
        (persistenceTime > anc_eval_diff)
    ),
    introductionsFromUnique = sum(
      (independenceTime > anc_eval_diff) &
        (persistenceTime <= anc_eval_diff)
    ),
    propPersistentFromUnique = if_else((persistentsFromUnique + introductionsFromUnique) > 0, persistentsFromUnique/(persistentsFromUnique + introductionsFromUnique), 0)
  ) %>%
  mutate(
    evaluationTime = as.Date(date_decimal(mrsd - evaluationTime)),
    ancestralTime = as.Date(date_decimal(mrsd - ancestralTime)),
    treeId = as.factor(treeId)
  )


# Get lockdown tiers for unique ancestral and evalTimes
anc_eval_times_tier <- df_prop %>% 
  ungroup() %>%
  distinct(stateAtEvaluationTime, ancestralTime, evaluationTime) 

anc_eval_times_tier_tmp <- anc_eval_times_tier %>%
  group_by(stateAtEvaluationTime, ancestralTime, evaluationTime) %>%
  group_modify(~{
    .x %>%
      mutate(
        lockdown_tier = lockdown_dates %>% 
          filter(comuna == .y$stateAtEvaluationTime & 
                   date <= .y$evaluationTime & 
                   date >= .y$ancestralTime) %>% 
          count(lockdown_tier) %>% 
          arrange(-n) %>% 
          pull(lockdown_tier) %>%
          first() %>%
          as.character() 
      )
    }
  )

# Join with lockdown tiers
df_prop_tiered <- df_prop %>%
  left_join(anc_eval_times_tier_tmp, by = c("stateAtEvaluationTime", "ancestralTime", "evaluationTime"))

# Get binned weeks from ancestral to evaluationTime
df_prop_tiered <- df_prop_tiered %>%
  arrange(stateAtEvaluationTime, ancestralTime, evaluationTime) %>%
  mutate(
    diff_weeks = (evaluationTime - ancestralTime) %>% as.numeric("weeks") %>% round()
  )

# Export to summarize across all lineages using other script
summarized_out_folder <- paste0("./summarized_outputs/", lineage, "/")
summarized_out_name <- paste0(summarized_out_folder, lineage, "_summarized.tsv")
dir.create(summarized_out_folder)
df_prop_tiered %>%
  write_tsv(summarized_out_name)
print(paste0("Exported to ", summarized_out_name))

summarise_hpd_lower <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x))[1])
}

summarise_hpd_upper <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x))[2])
}

df_summarized_all  <- df_prop_tiered %>%
  group_by(lockdown_tier, diff_weeks) %>%
  summarise(
    # Proportion of persistent lineages
    propPersistentFromUnique_lower = summarise_hpd_lower(propPersistentFromUnique),
    propPersistentFromUnique_upper = summarise_hpd_upper(propPersistentFromUnique),
    propPersistentFromUnique_median = median(propPersistentFromUnique),
    # Number of persistent lineages
    persistentsFromUnique_lower = summarise_hpd_lower(persistentsFromUnique),
    persistentsFromUnique_upper = summarise_hpd_upper(persistentsFromUnique),
    persistentsFromUnique_median = median(persistentsFromUnique),
    # Number of unique introductions
    introductionsFromUnique_lower = summarise_hpd_lower(introductionsFromUnique),
    introductionsFromUnique_upper = summarise_hpd_upper(introductionsFromUnique),
    introductionsFromUnique_median = median(introductionsFromUnique)
  )



p_prop <- df_summarized_all %>%
  ggplot() +
  geom_pointrange(aes(x = diff_weeks, propPersistentFromUnique_median, ymin=propPersistentFromUnique_lower, ymax = propPersistentFromUnique_upper)) +
  facet_grid(lockdown_tier ~ ., scales="free") + theme_bw() + xlab("Number of weeks since ancestral time")

p_num_persistent <- df_summarized_all %>%
  ggplot() +
  geom_pointrange(aes(x = diff_weeks, persistentsFromUnique_median, ymin=persistentsFromUnique_lower, ymax = persistentsFromUnique_upper)) +
  facet_grid(lockdown_tier ~ ., scales="free") + theme_bw() + xlab("Number of weeks since ancestral time")

p_num_unique <- df_summarized_all %>%
  ggplot() +
  geom_pointrange(aes(x = diff_weeks, introductionsFromUnique_median, ymin=introductionsFromUnique_lower, ymax = introductionsFromUnique_upper)) +
  facet_grid(lockdown_tier ~ ., scales="free") + theme_bw() + xlab("Number of weeks since ancestral time")

p_prop + p_num_persistent + p_num_unique
ggsave(paste0("./plots/", lineage, "_binned_by_weeks.pdf"), w = 15, h = 7.5)

# Proportion of lineages at evalTime
df_prop_descendants_summarised <- df_prop_tiered %>%
  mutate(    
    propDescendantsFromIntroductionsAtEval = introductionLineagesAtEval/(introductionLineagesAtEval + persistentLineagesAtEval)
  ) %>%
  group_by(lockdown_tier, diff_weeks) %>%
  summarise(
    propDescendantsFromIntroductionsAtEval_lower = summarise_hpd_lower(propDescendantsFromIntroductionsAtEval),
    propDescendantsFromIntroductionsAtEval_upper = summarise_hpd_upper(propDescendantsFromIntroductionsAtEval),
    propDescendantsFromIntroductionsAtEval_median = median(propDescendantsFromIntroductionsAtEval)
  )

df_prop_descendants_summarised %>%
  ggplot() +
  geom_pointrange(aes(x = diff_weeks, propDescendantsFromIntroductionsAtEval_median, ymin=propDescendantsFromIntroductionsAtEval_lower, ymax = propDescendantsFromIntroductionsAtEval_upper)) +
  facet_grid(lockdown_tier ~ ., scales="free") + theme_bw() + xlab("Number of weeks since ancestral time")
ggsave(paste0("./plots/",lineage,"_number_of_descendants_introductions_binned_by_week.pdf"), w = 7.5, h = 5)
