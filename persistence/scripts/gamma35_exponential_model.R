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

## Import Gamma TL 35 summary file
gamma35_persistence <- read.csv("persistence/summarized_outputs/Gamma_TL_35/Gamma_TL_35_summarized.tsv",
                                sep = "\t") |>
  group_by(lockdown_tier, diff_weeks) |>
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

# Add very small constant to avoid zero values
small_constant <- 0.001

gamma35_persistence$propPersistentFromUnique_median_corr <-
  gamma35_persistence$propPersistentFromUnique_median + small_constant

gamma35_persistence$propPersistentFromUnique_upper_corr <-
  gamma35_persistence$propPersistentFromUnique_upper + small_constant

gamma35_persistence$propPersistentFromUnique_lower_corr <-
  gamma35_persistence$propPersistentFromUnique_lower + small_constant

## Fit exponential function to each lockdown tier
# Fit to media values
g35_expfit_full <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      data = gamma35_persistence |>
                        filter(lockdown_tier == "Full lockdown"))

g35_expfit_wknd <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      data = gamma35_persistence |>
                        filter(lockdown_tier == "Weekend lockdown"))

g35_expfit_no <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                    data = gamma35_persistence |>
                      filter(lockdown_tier == "No lockdown"))

# Weigh by HPD
g35_expfit_hpd_full <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      weights = 1 / (propPersistentFromUnique_upper -
                                       propPersistentFromUnique_lower +
                                       small_constant),
                      data = gamma35_persistence |>
                        filter(lockdown_tier == "Full lockdown"))

g35_expfit_hpd_wknd <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      weights = 1 / (propPersistentFromUnique_upper -
                                       propPersistentFromUnique_lower +
                                       small_constant),
                      data = gamma35_persistence |>
                        filter(lockdown_tier == "Weekend lockdown"))

g35_expfit_hpd_no <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                    weights = 1 / (propPersistentFromUnique_upper -
                                     propPersistentFromUnique_lower +
                                     small_constant),
                      data = gamma35_persistence |>
                        filter(lockdown_tier == "No lockdown"))

# Estimate decay rate per lockdown tier
slope_full <- coef(g35_expfit_full)[2]
slope_wknd <- coef(g35_expfit_wknd)[2]
slope_no <- coef(g35_expfit_no)[2]

rate_full <- exp(slope_full) -1
rate_wknd <- exp(slope_wknd) -1
rate_no <- exp(slope_no) -1


## Plot Gamma TL 35 proportion of lineages over time per lockdown tier
g35_persist_prop <- gamma35_persistence |>
  ggplot() +
  geom_point(aes(x = diff_weeks,
                 y = propPersistentFromUnique_median), size = 3) +
  geom_line(aes(x = diff_weeks,
                y = propPersistentFromUnique_median)) +
  geom_ribbon(aes(x = diff_weeks,
                ymin = propPersistentFromUnique_lower,
                ymax = propPersistentFromUnique_upper), alpha = 0.1) +
  facet_grid(lockdown_tier ~ ., scales="free") + theme_bw() + xlab("Number of weeks since ancestral time")
