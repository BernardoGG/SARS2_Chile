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

# Split data frames by lockdown tier
gamma35_persistence_full <- gamma35_persistence |>
  filter(lockdown_tier == "Full lockdown")

gamma35_persistence_wknd <- gamma35_persistence |>
  filter(lockdown_tier == "Weekend lockdown")

gamma35_persistence_no <- gamma35_persistence |>
  filter(lockdown_tier == "No lockdown")

## Fit exponential function to each lockdown tier
# Fit to media values
g35_expfit_full <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      data = gamma35_persistence_full)

g35_expfit_wknd <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      data = gamma35_persistence_wknd)

g35_expfit_no <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                    data = gamma35_persistence_no)

# Weighted by HPD
g35_expfit_hpd_full <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      weights = 1 / (propPersistentFromUnique_upper -
                                       propPersistentFromUnique_lower +
                                       small_constant),
                      data = gamma35_persistence_full)

g35_expfit_hpd_wknd <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                      weights = 1 / (propPersistentFromUnique_upper -
                                       propPersistentFromUnique_lower +
                                       small_constant),
                      data = gamma35_persistence_wknd)

g35_expfit_hpd_no <- lm(log(propPersistentFromUnique_median_corr) ~ diff_weeks,
                    weights = 1 / (propPersistentFromUnique_upper -
                                     propPersistentFromUnique_lower +
                                     small_constant),
                      data = gamma35_persistence_no)

# Estimate decay rate per lockdown tier
slope_full <- coef(g35_expfit_full)[2]
slope_wknd <- coef(g35_expfit_wknd)[2]
slope_no <- coef(g35_expfit_no)[2]

rate_full <- exp(slope_full) - 1
rate_wknd <- exp(slope_wknd) - 1
rate_no <- exp(slope_no) - 1

slope_hpd_full <- coef(g35_expfit_hpd_full)[2]
slope_hpd_wknd <- coef(g35_expfit_hpd_wknd)[2]
slope_hpd_no <- coef(g35_expfit_hpd_no)[2]

rate_hpd_full <- exp(slope_hpd_full) - 1
rate_hpd_wknd <- exp(slope_hpd_wknd) - 1
rate_hpd_no <- exp(slope_hpd_no) - 1


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
  facet_grid(lockdown_tier ~ ., scales="free") + theme_bw() +
  xlab("Number of weeks since ancestral time")



## Plot the original data and the fitted curve
a_full <- exp(coef(g35_expfit_full)[1])
a_wknd <- exp(coef(g35_expfit_wknd)[1])
a_no <- exp(coef(g35_expfit_no)[1])

# Create a sequence of time values for the fitted curves
fitted_time_full <- seq(min(gamma35_persistence_full$diff_weeks),
                        max(gamma35_persistence_full$diff_weeks),
                        length.out = 100)
fitted_time_wknd <- seq(min(gamma35_persistence_wknd$diff_weeks),
                        max(gamma35_persistence_wknd$diff_weeks),
                        length.out = 100)
fitted_time_no <- seq(min(gamma35_persistence_no$diff_weeks),
                        max(gamma35_persistence_no$diff_weeks),
                        length.out = 100)

# Calculate the fitted values using the exponential model
fitted_values_full <- a_full * exp(slope_full * fitted_time_full)
fitted_values_wknd <- a_full * exp(slope_wknd * fitted_time_wknd)
fitted_values_no <- a_full * exp(slope_no * fitted_time_no)

# Fitted values data frames
g35_fitted_full <- data.frame(time = fitted_time_full,
                              fit_p = fitted_values_full)
g35_fitted_wknd <- data.frame(time = fitted_time_wknd,
                              fit_p = fitted_values_wknd)
g35_fitted_no <- data.frame(time = fitted_time_no,
                              fit_p = fitted_values_no)

## Plot Gamma 35 decay rates per lockdown tier
(ggplot(gamma35_persistence_full) +
    geom_point(aes(x = diff_weeks,
                   y = propPersistentFromUnique_median), size = 3.5) +
    geom_line(aes(x = diff_weeks,
                  y = propPersistentFromUnique_median), linewidth = 1) +
    geom_ribbon(aes(x = diff_weeks,
                    ymin = propPersistentFromUnique_lower,
                    ymax = propPersistentFromUnique_upper), alpha = 0.1) +
    geom_line(data = g35_fitted_full, aes(x = time, y = fit_p),
              colour = "darkred", alpha = 0.3, linewidth = 2) +
    xlim(0,15) + theme_minimal() + theme(axis.title.x = element_blank(),
                                         axis.title.y = element_blank())) /
  (ggplot(gamma35_persistence_wknd) +
     geom_point(aes(x = diff_weeks,
                    y = propPersistentFromUnique_median), size = 3.5) +
     geom_line(aes(x = diff_weeks,
                   y = propPersistentFromUnique_median), linewidth = 1) +
     geom_ribbon(aes(x = diff_weeks,
                     ymin = propPersistentFromUnique_lower,
                     ymax = propPersistentFromUnique_upper), alpha = 0.1) +
     geom_line(data = g35_fitted_wknd, aes(x = time, y = fit_p),
               colour = "darkred", alpha = 0.3, linewidth = 2) +
     xlim(0,15) + labs(y = "Proportion of persisting lineages") +
     theme_minimal() + theme(axis.title.x = element_blank())) /
  (ggplot(gamma35_persistence_no) +
     geom_point(aes(x = diff_weeks,
                    y = propPersistentFromUnique_median), size = 3.5) +
     geom_line(aes(x = diff_weeks,
                   y = propPersistentFromUnique_median), linewidth = 1) +
     geom_ribbon(aes(x = diff_weeks,
                     ymin = propPersistentFromUnique_lower,
                     ymax = propPersistentFromUnique_upper), alpha = 0.1) +
     geom_line(data = g35_fitted_no, aes(x = time, y = fit_p),
               colour = "darkred", alpha = 0.3, linewidth = 2) +
     xlim(0,15) + labs(x = "Weeks after lockdown implementation") +
     theme_minimal() + theme(axis.title.y = element_blank()))

ggsave("../Figures/Gamma_35_persitence_decay.pdf", dpi = 300,
       height = 10, width = 7.15, bg = "white")
ggsave("../Figures/Gamma_35_persitence_decay.png", dpi = 300,
       height = 10, width = 7.15, bg = "white")
