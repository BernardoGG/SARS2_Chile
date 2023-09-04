################################################################################
################ SARS-CoV-2 Chile lockdown effects modelling ###################
################################################################################

# Source scripts, packages and data
source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")
library(lme4)
library(MASS)
library(pscl)
library(glmmTMB)
library(rcompanion)
library(tidyverse)
library(lubridate)
library(patchwork)

########################### Bernardo Gutierrez #################################
### Set data frames ####
## Add column with numbers of new cases per report date
CL_invasions <- CL_invasions |>
  left_join(comuna_cases, by = c("head_comuna" = "Comuna",
                              "head_date" = "Fecha")) |>
  select(-`Casos Confirmados`, head_date_new_cases = `Casos nuevos`)

CL_invasions <- CL_invasions |>
  left_join(comuna_cases, by = c("tail_comuna" = "Comuna",
                                 "tail_date" = "Fecha")) |>
  select(-`Casos Confirmados`, tail_date_new_cases = `Casos nuevos`)


## Internal movement and importation tables
# Include counts and total new cases per comuna per lockdown tier
CL_invasions_model_within <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_model_between <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

# Include counts, total new cases per comuna per lockdown tier and dates
CL_invasions_model_time_within <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded, tail_date) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_model_time_between <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded, tail_date) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

# Include counts and total new cases per comuna per lockdown tier before and
# after implementation of the mobility pass (May 26)
CL_invasions_no_MobPass_within <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  filter(tail_date < as.Date("2021-05-26")) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_MobPass_within <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  filter(tail_date >= as.Date("2021-05-26")) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_no_MobPass_between <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  filter(tail_date < as.Date("2021-05-26")) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_MobPass_between <- CL_invasions |>
  filter(head_comuna != tail_comuna) |>
  filter(tail_date >= as.Date("2021-05-26")) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

# Include counts and total new cases per comuna per lockdown tier before and
# after implementation of the mobility pass (May 26)
CL_invasions_MobPassVar_within <- CL_invasions |>
  mutate(MobPass = ifelse(tail_date >= as.Date("2021-05-26"), "Yes", "No")) |>
  filter(head_comuna == tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded, MobPass) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_MobPassVar_between <- CL_invasions |>
  mutate(MobPass = ifelse(tail_date >= as.Date("2021-05-26"), "Yes", "No")) |>
  filter(head_comuna != tail_comuna) |>
  group_by(tail_comuna, tail_lockdown_tier_recoded, MobPass) |>
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(count = as.numeric(count)) |>
  mutate(tail_comuna = as.factor(tail_comuna)) |>
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))


### Model setup ####
####### Model 1 #######
## LM of within-comuna movements, accounting for new cases
m1_within <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
        data = CL_invasions_model_within)
# Observed vs predicted values of within-comuna movements
m1_within_pred <- predict(m1_within)

## LM of cross-comuna movements, accounting for new cases
m1_between <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
               data = CL_invasions_model_between)
# Observed vs predicted values of cross-comuna movements
m1_between_pred <- predict(m1_between)

####### Model 2 #######
## Poisson GLM of within-comuna movements, accounting for new cases
m2_within <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
                 data = CL_invasions_model_within,
          family = "poisson")
# Observed vs predicted values of within-comuna movements
m2_within_pred <- predict(m2_within)

## Poisson GLM of cross-comuna movements, accounting for new cases
m2_between <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
                 data = CL_invasions_model_between,
                 family = "poisson")
# Observed vs predicted values of cross-comuna movements
m2_between_pred <- predict(m2_between)

####### Model 3 #######
## Negative binomial GLM of within-comuna movements, accounting for new cases
m3_within <- glm.nb(count ~ tail_lockdown_tier_recoded + total_cases,
                     data = CL_invasions_model_time_within)
# Observed vs predicted values of cross-comuna movements
m3_within_pred <- predict(m3_within)

## Negative binomial GLM of cross-comuna movements, accounting for new cases
m3_between <- glm.nb(count ~ tail_lockdown_tier_recoded + total_cases,
                  data = CL_invasions_model_time_between)
# Observed vs predicted values of cross-comuna movements
m3_between_pred <- predict(m3_between)

####### Model 4 #######
## LM of within-comuna movements, accounting for new cases and Mobility Pass
m4_within <- glm(count ~ tail_lockdown_tier_recoded + total_cases + MobPass,
                 data = CL_invasions_MobPassVar_within)
# Observed vs predicted values of within-comuna movements
m4_within_pred <- predict(m4_within)

## LM of cross-comuna movements, accounting for new cases and Mobility Pass
m4_between <- glm(count ~ tail_lockdown_tier_recoded + total_cases + MobPass,
                  data = CL_invasions_MobPassVar_between)
# Observed vs predicted values of cross-comuna movements
m4_between_pred <- predict(m4_between)

####### Model 5 #######
## Poisson GLM of within-comuna movements, accounting for new cases and Mobility Pass
m5_within <- glm(count ~ tail_lockdown_tier_recoded + total_cases + MobPass,
                 data = CL_invasions_MobPassVar_between,
                 family = "poisson")
# Observed vs predicted values of within-comuna movements
m5_within_pred <- predict(m5_within)

## Poisson GLM of cross-comuna movements, accounting for new cases and Mobility Pass
m5_between <- glm(count ~ tail_lockdown_tier_recoded + total_cases + MobPass,
                  data = CL_invasions_MobPassVar_between,
                  family = "poisson")
# Observed vs predicted values of cross-comuna movements
m5_between_pred <- predict(m5_between)

####### Model 6 #######
## Negative binomial GLM of within-comuna movements, accounting for new cases and date
m6_within <- glm.nb(count ~ tail_lockdown_tier_recoded + total_cases + MobPass,
                    data = CL_invasions_MobPassVar_within, maxit = 1000)
# Observed vs predicted values of cross-comuna movements
m6_within_pred <- predict(m6_within)

## Negative binomial GLM of cross-comuna movements, accounting for new cases and date
m6_between <- glm.nb(count ~ tail_lockdown_tier_recoded + total_cases + MobPass,
                     data = CL_invasions_MobPassVar_between, maxit = 1000)
# Observed vs predicted values of cross-comuna movements
m6_between_pred <- predict(m6_between)


### Model test and comparisons ####
####### Within comunas #######
summary(m1_within)
summary(m2_within)
summary(m3_within)
summary(m4_within)
summary(m5_within)
summary(m6_within)

####### Between comunas #######
summary(m1_between)
summary(m2_between)
summary(m3_between)
summary(m4_between)
summary(m5_between)
summary(m6_between)

####### Model tests #######
# Compare models with and without Mobility Pass variable
compareGLM(m1_within, m4_within)
compareGLM(m2_within, m5_within)
compareGLM(m3_within, m6_within)

compareGLM(m1_between, m4_between)
compareGLM(m2_between, m5_between)
compareGLM(m3_between, m6_between)

## LRT between Poisson and NB models to test goodness of fit based on dispersion
## parameter
pchisq(2 * (logLik(m3_within) - logLik(m2_within)), df = 1, lower.tail = FALSE)
pchisq(2 * (logLik(m3_between) - logLik(m2_between)), df = 1, lower.tail = FALSE)

pchisq(2 * (logLik(m6_within) - logLik(m5_within)), df = 1, lower.tail = FALSE)
pchisq(2 * (logLik(m6_between) - logLik(m5_between)), df = 1, lower.tail = FALSE)

(est_within <- cbind(Estimate = coef(m3_within), confint(m3_within)))
(est_between <- cbind(Estimate = coef(m3_between), confint(m3_between)))
exp(est_within)
exp(est_between)

(est_within <- cbind(Estimate = coef(m6_within), confint(m6_within)))
(est_between <- cbind(Estimate = coef(m6_between), confint(m6_between)))
exp(est_within)
exp(est_between)

## Plot deviance residuals against predicted values to evaluate model fit
# NOTE: Dev. resid. expected between -2 and 2, tighter values show better model
# fit.
m1_within_resids <- data.frame(predicted_values = m1_within_pred,
                               deviance_residuals = resid(m1_within))
m2_within_resids <- data.frame(predicted_values = m2_within_pred,
                               deviance_residuals = resid(m2_within))
m3_within_resids <- data.frame(predicted_values = m3_within_pred,
                               deviance_residuals = resid(m3_within))
m4_within_resids <- data.frame(predicted_values = m4_within_pred,
                               deviance_residuals = resid(m4_within))
m5_within_resids <- data.frame(predicted_values = m5_within_pred,
                               deviance_residuals = resid(m5_within))
m6_within_resids <- data.frame(predicted_values = m6_within_pred,
                               deviance_residuals = resid(m6_within))

m1_between_resids <- data.frame(predicted_values = m1_between_pred,
                               deviance_residuals = resid(m1_between))
m2_between_resids <- data.frame(predicted_values = m2_between_pred,
                               deviance_residuals = resid(m2_between))
m3_between_resids <- data.frame(predicted_values = m3_between_pred,
                               deviance_residuals = resid(m3_between))
m4_between_resids <- data.frame(predicted_values = m4_between_pred,
                                deviance_residuals = resid(m4_between))
m5_between_resids <- data.frame(predicted_values = m5_between_pred,
                                deviance_residuals = resid(m5_between))
m6_between_resids <- data.frame(predicted_values = m6_between_pred,
                                deviance_residuals = resid(m6_between))

(ggplot(m1_within_resids, aes(x = predicted_values, y = deviance_residuals)) +
  geom_point(colour = 'gray') + labs(title = "Linear", x = "", y = "") +
    theme_minimal() |
ggplot(m1_between_resids, aes(x = predicted_values, y = deviance_residuals)) +
  geom_point(colour = 'gray') + labs(title = "", x = "", y = "") +
    theme_minimal()) /

(ggplot(m2_within_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
  geom_point(colour = 'gray') + labs(title = "Poisson", x = "",
                                     y = "Deviance residuals") +
   theme_minimal() |
ggplot(m2_between_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
  geom_point(colour = 'gray') + labs(title = "", x = "", y = "") +
  theme_minimal()) /

(ggplot(m3_within_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
  geom_point(colour = 'gray') + labs(title = "Negative binomial",
                                     x = "Predicted values", y = "") +
   theme_minimal() |
ggplot(m3_between_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
  geom_point(colour = 'gray') + labs(title = "",
                                     x = "Predicted values", y = "") +
  theme_minimal())

(ggplot(m4_within_resids, aes(x = predicted_values, y = deviance_residuals)) +
    geom_point(colour = 'gray') + labs(title = "Linear", x = "", y = "") +
    theme_minimal() |
    ggplot(m4_between_resids, aes(x = predicted_values, y = deviance_residuals)) +
    geom_point(colour = 'gray') + labs(title = "", x = "", y = "") +
    theme_minimal()) /
  
  (ggplot(m5_within_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
     geom_point(colour = 'gray') + labs(title = "Poisson", x = "",
                                        y = "Deviance residuals") +
     theme_minimal() |
     ggplot(m5_between_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
     geom_point(colour = 'gray') + labs(title = "", x = "", y = "") +
     theme_minimal()) /
  
  (ggplot(m6_within_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
     geom_point(colour = 'gray') + labs(title = "Negative binomial",
                                        x = "Predicted values", y = "") +
     theme_minimal() |
     ggplot(m6_between_resids, aes(x = log(predicted_values), y = deviance_residuals)) +
     geom_point(colour = 'gray') + labs(title = "",
                                        x = "Predicted values", y = "") +
     theme_minimal())

## Test effects of lockdown tier by comparing model which includes and excludes
## this variable - applicable for model that best fits the data (NB)
m3_within_noLD <- update(m3_within, . ~ . - tail_lockdown_tier_recoded)
m3_between_noLD <- update(m3_between, . ~ . - tail_lockdown_tier_recoded)

anova(m3_within, m3_within_noLD)
anova(m3_between, m3_between_noLD)

## Test effects of mobility pass by comparing model which includes and excludes
## this variable - applicable for model that best fits the data (NB)
m6_within_noMP <- update(m6_within, . ~ . - MobPass)
m6_between_noMP <- update(m6_between, . ~ . - MobPass)

anova(m6_within, m6_within_noMP)
anova(m6_between, m6_between_noMP)

## Test full model versus minimal model
anova(m3_within_noLD, m6_within)
anova(m3_between_noLD, m6_between)


## Get NB estimates for data generated under those models holding numbers of
## cases at their mean.
new_m3_within <-
  data.frame(total_cases = mean(CL_invasions_model_within$total_cases),
             tail_lockdown_tier_recoded = factor(1:3, levels = 1:3,
                                    labels = levels(
                                      CL_invasions_model_within$tail_lockdown_tier_recoded)))
new_m3_within$phat <- predict(m3_within, new_m3_within, type = "response")
new_m3_within

new_m3_between <-
  data.frame(total_cases = mean(CL_invasions_model_between$total_cases),
             tail_lockdown_tier_recoded = factor(1:3, levels = 1:3,
                                                 labels = levels(
                                                   CL_invasions_model_between$tail_lockdown_tier_recoded)))
new_m3_between$phat <- predict(m3_between, new_m3_between, type = "response")
new_m3_between

## Plot predicted viral movement estimates over the range of possible new cases
## by lockdown tier
new_m3_within_plot <- data.frame(
  total_cases = rep(seq(from = min(CL_invasions_model_within$total_cases),
                        to = max(CL_invasions_model_within$total_cases),
                        length.out = 100), 3),
  tail_lockdown_tier_recoded = factor(rep(1:3, each = 100),
                                          levels = 1:3,
                                          labels = levels(
                       CL_invasions_model_within$tail_lockdown_tier_recoded)))

new_m3_within_plot <- cbind(new_m3_within_plot,
                            predict(m3_within, new_m3_within_plot, type = "link",
                              se.fit = TRUE))

new_m3_within_plot <- within(new_m3_within_plot, {
  count = exp(fit)
  LL <- exp(fit - 1.96*se.fit)
  UL <- exp(fit + 1.96*se.fit)
})

new_m3_between_plot <- data.frame(
  total_cases = rep(seq(from = min(CL_invasions_model_between$total_cases),
                        to = max(CL_invasions_model_between$total_cases),
                        length.out = 100), 3),
  tail_lockdown_tier_recoded = factor(rep(1:3, each = 100),
                                          levels = 1:3,
                                          labels = levels(
                       CL_invasions_model_between$tail_lockdown_tier_recoded)))

new_m3_between_plot <- cbind(new_m3_between_plot,
                            predict(m3_between, new_m3_between_plot, type = "link",
                              se.fit = TRUE))

new_m3_between_plot <- within(new_m3_between_plot, {
  count = exp(fit)
  LL <- exp(fit - 1.96*se.fit)
  UL <- exp(fit + 1.96*se.fit)
})

# Plots
ggplot(new_m3_within_plot, aes(total_cases, count)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = tail_lockdown_tier_recoded),
              alpha = .25) +
  geom_line(aes(colour = tail_lockdown_tier_recoded), size = 2) +
  labs(x = "Total cases per comuna", y = "Predicted viral movements")

ggplot(new_m3_between_plot, aes(total_cases, count)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = tail_lockdown_tier_recoded),
              alpha = .25) +
  geom_line(aes(colour = tail_lockdown_tier_recoded), size = 2) +
  labs(x = "Total cases per comuna", y = "Predicted viral movements")


### Model setup - within viral movements as branching events ####
## Estimate frequency of branching events under different lockdown tiers
# Counts of unique occurrences of comunas as 'head' nodes
# Corresponds to all branching event taking place per comuna
ancestor_total_counts <- CL_invasions |>
  group_by(head_comuna) |>
  summarise(count = n())

# Counts of unique occurrences of comunas as 'head' nodes with the same 'tail' node
# Corresponds to all branching events resulting in descendants within the same comuna
ancestor_local_counts <- CL_invasions |>
  filter(head_comuna == tail_comuna) |>
  group_by(head_comuna) |>
  summarise(count = n())

# All LD tiers
counts_all_tiers <- CL_invasions |>
  group_by(head_comuna, tail_comuna) |>
  summarise(count = n())
ancestor_counts_all_tiers <- counts_all_tiers |>
  filter(head_comuna == tail_comuna) |>
  group_by(head_comuna) |>
  summarise(branching = sum(count, na.rm = TRUE))

# Full lockdown
counts_full_ld <- CL_invasions |>
  filter(head_lockdown_tier_recoded == "Full lockdown") |>
  group_by(head_comuna, tail_comuna) |>
  summarise(count = n())
ancestor_counts_full_ld <- counts_full_ld |>
  filter(head_comuna == tail_comuna) |>
  group_by(head_comuna) |>
  summarise(branching = sum(count, na.rm = TRUE))

# Weekend lockdown
counts_weekend_ld <- CL_invasions |>
  filter(head_lockdown_tier_recoded == "Weekend lockdown") |>
  group_by(head_comuna, tail_comuna) |>
  summarise(count = n())
ancestor_counts_weekend_ld <- counts_weekend_ld |>
  filter(head_comuna == tail_comuna) |>
  group_by(head_comuna) |>
  summarise(branching = sum(count, na.rm = TRUE))

# No lockdown
counts_no_ld <- CL_invasions |>
  filter(head_lockdown_tier_recoded == "No lockdown") |>
  group_by(head_comuna, tail_comuna) |>
  summarise(count = n())
ancestor_counts_no_ld <- counts_no_ld |>
  filter(head_comuna == tail_comuna) |>
  group_by(head_comuna) |>
  summarise(branching = sum(count, na.rm = TRUE))

branching_all_tiers <- left_join(ancestor_local_counts, ancestor_counts_all_tiers) |>
  mutate(branching = ifelse(is.na(branching), 0, branching)) |>
  mutate(branching_freq = branching / count)

branching_full_ld <- left_join(ancestor_local_counts, ancestor_counts_full_ld) |>
  mutate(branching = ifelse(is.na(branching), 0, branching)) |>
  mutate(branching_freq = branching / count)

branching_weekend_ld <- left_join(ancestor_local_counts, ancestor_counts_weekend_ld) |>
  mutate(branching = ifelse(is.na(branching), 0, branching)) |>
  mutate(branching_freq = branching / count)

branching_no_ld <- left_join(ancestor_local_counts, ancestor_counts_no_ld) |>
  mutate(branching = ifelse(is.na(branching), 0, branching)) |>
  mutate(branching_freq = branching / count)

## Plot frequencies of branching events under different lockdown tiers
# Frequencies of branching events
full_hist <- ggplot(branching_full_ld) +
  geom_histogram(aes(x = branching_freq), fill = "#B81D13", bins = 21) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(x = "", y = "Full lockdown") + theme_minimal()

weekend_hist <- ggplot(branching_weekend_ld) +
  geom_histogram(aes(x = branching_freq), fill = "#EFB700", bins = 21) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(x = "", y = "Weekend lockdown") + theme_minimal()

no_hist <- ggplot(branching_no_ld) +
  geom_histogram(aes(x = branching_freq), fill = "#008450", bins = 21) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(x = "Within-comuna branching frequency", y = "No lockdown") + theme_minimal()

full_hist / weekend_hist / no_hist

ggsave("../Figures/CL_branching_frequencies_hist.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")
ggsave("../Figures/CL_branching_frequencies_hist.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")


full_dens <- ggplot(branching_full_ld) +
  geom_density(aes(x = branching_freq), kernel = "gaussian", fill = "#B81D13",
               alpha = 0.5) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x = "", y = "Full lockdown") + theme_minimal()

weekend_dens <- ggplot(branching_weekend_ld) +
  geom_density(aes(x = branching_freq), kernel = "gaussian", fill = "#EFB700",
               alpha = 0.5) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x = "", y = "Weekend lockdown") + theme_minimal()

no_dens <- ggplot(branching_no_ld) +
  geom_density(aes(x = branching_freq), kernel = "gaussian", fill = "#008450",
               alpha = 0.5) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x = "Within-comuna branching frequency", y = "No lockdown") + theme_minimal()

full_dens / weekend_dens / no_dens

ggsave("../Figures/CL_branching_frequencies_dens.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")
ggsave("../Figures/CL_branching_frequencies_dens.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")

# Counts of branching events
ref_c <- ggplot(ancestor_total_counts) +
  geom_histogram(aes(x = count), bins = 10)

full_hist_c <- ggplot(branching_full_ld) +
  geom_histogram(aes(x = branching), fill = "#B81D13", bins = 36) +
  scale_x_continuous(limits = c(-6, 400)) +
  scale_y_continuous(limits = c(0, 180)) +
  labs(x = "", y = "Full lockdown") + theme_minimal()

weekend_hist_c <- ggplot(branching_weekend_ld) +
  geom_histogram(aes(x = branching), fill = "#EFB700", bins = 36) +
  scale_x_continuous(limits = c(-6, 400)) +
  scale_y_continuous(limits = c(0, 180)) +
  labs(x = "", y = "Weekend lockdown") + theme_minimal()

no_hist_c <- ggplot(branching_no_ld) +
  geom_histogram(aes(x = branching), fill = "#008450", bins = 36) +
  scale_x_continuous(limits = c(-6, 400)) +
  scale_y_continuous(limits = c(0, 180)) +
  labs(x = "Within-comuna branching events", y = "No lockdown") + theme_minimal()

full_hist_c / weekend_hist_c / no_hist_c

# Comparison of branching frequencies
left_join(
  left_join(
    branching_full_ld |> rename(branching_freq_full = branching_freq) |>
      dplyr::select(head_comuna, branching_freq_full),
    branching_weekend_ld |> rename(branching_freq_wkd = branching_freq) |>
      dplyr::select(head_comuna, branching_freq_wkd)
  ),
  branching_no_ld |> rename(branching_freq_no = branching_freq) |>
    dplyr::select(head_comuna, branching_freq_no)
) |>
  pivot_longer(
    !head_comuna, names_to = "lockdown_tier", values_to = "branching_freq"
  ) |>
  ggplot(aes(x = lockdown_tier,
             y = branching_freq, fill = lockdown_tier)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.6, width = 0.2, colour = "#555555",
              outlier.shape = NA) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', stroke = 0.5,
               stackratio = 0.1, dotsize = 0.6, alpha = 0.1,
               position = position_jitter(
                 width = 0.2, height = 0.05)) +
  scale_fill_manual(values = c("#B81D13", "#EFB700", "#008450")) +
  scale_x_discrete(labels = c("Full lockdown", "Weekend lockdown", "No lockdown")) +
  labs(x = "", colour = "",
       y = "Within-comuna branching frequency") +
  theme_classic() + theme(axis.text.x = element_text(size = 15),
                          axis.text.y = element_text(size = 10),
                          axis.title.y = element_text(size = 15),
                          axis.ticks.x = element_blank(),
                          legend.position = "none")







############# Sandbox ##############
## Set model test for branching events within comunas under lockdown tiers
CL_invasions_branching <- left_join(
  left_join(
    branching_full_ld |> rename(`Full lockdown` = branching) |>
      dplyr::select(head_comuna, `Full lockdown`),
    branching_weekend_ld |> rename(`Weekend lockdown` = branching) |>
      dplyr::select(head_comuna, `Weekend lockdown`)
  ),
  branching_no_ld |> rename(`No lockdown` = branching) |>
    dplyr::select(head_comuna, `No lockdown`)
) |>
  pivot_longer(
    !head_comuna, names_to = "lockdown_tier", values_to = "branching_count"
  ) |>
  left_join(CL_invasions |>
              dplyr::select(head_comuna, head_lockdown_tier_recoded, head_date_new_cases) |>
              group_by(head_comuna, head_lockdown_tier_recoded) |>
              summarise(total_cases = sum(head_date_new_cases)),
            by = c("head_comuna" = "head_comuna",
                   "lockdown_tier" = "head_lockdown_tier_recoded"))






## Model 4: LM of within-comuna branching events, accounting for new cases
m4_branching <- glm(branching_count ~ lockdown_tier + total_cases,
                    data = CL_invasions_branching)

# Observed vs predicted values of within-comuna branching events
m4_branching_pred <- predict(m4_branching)
plot(m4_branching_pred, CL_invasions_branching$branching_count)


## Model 5: Poisson GLM of within-comuna branching events, accounting for new cases
m5_branching <- glm(branching_count ~ lockdown_tier + total_cases,
                    data = CL_invasions_branching,
                    family = "poisson")

# Observed vs predicted values of within-comuna branching events
m5_branching_pred <- predict(m5_branching)
plot(m5_branching_pred, CL_invasions_branching$branching_count)


## Model 6: Negative binomial GLM of within-comuna branching events, accounting for new cases
m6_branching <- glm.nb(branching_count ~ lockdown_tier + total_cases,
                       data = CL_invasions_branching)

# Observed vs predicted values of within-comuna branching events
m6_branching_pred <- predict(m6_branching)
plot(m6_branching_pred, CL_invasions_branching$branching_count)


## Model 7: Zero-inflated Poisson model of within-comuna branching events, accounting for new cases
m7_branching <- zeroinfl(branching_count ~
                           lockdown_tier + total_cases | lockdown_tier + total_cases,
                         data = CL_invasions_branching,
                         dist = "poisson")
m7_branching_pred <- predict(m7_branching)

# Observed vs predicted values of within-comuna branching events
m7_branching_pred <- predict(m7_branching)
plot(m7_branching_pred, CL_invasions_branching$branching_count)


### Model test and comparisons
summary(m4_branching)
summary(m5_branching)
summary(m6_branching)
summary(m7_branching)






summary(glm.nb(count ~ MobPass + tail_lockdown_tier_recoded + total_cases,
                data = CL_invasions_MobPassVar_within, maxit = 1000))
exp(m4_between$coefficients)

summary(glm.nb(count ~ MobPass + tail_lockdown_tier_recoded + total_cases,
               data = CL_invasions_MobPassVar_between, maxit = 1000))


summary(glm.nb(count ~ MobPass + tail_lockdown_tier_recoded + total_cases +
                 MobPass*tail_lockdown_tier_recoded,
               data = CL_invasions_MobPassVar_within, maxit = 1000))

summary(glm.nb(count ~ MobPass + tail_lockdown_tier_recoded + total_cases +
                 MobPass*tail_lockdown_tier_recoded,
               data = CL_invasions_MobPassVar_between, maxit = 1000))




# Mixed effects negative binomial models
summary(glmmTMB(count ~ MobPass + tail_lockdown_tier_recoded +
                   (1|total_cases),
               data = CL_invasions_MobPassVar_within,
               family = nbinom2))

summary(glmmTMB(count ~ MobPass + tail_lockdown_tier_recoded +
                  (1|total_cases),
                data = CL_invasions_MobPassVar_between,
                family = nbinom2))



summary(glmmTMB(count ~ tail_lockdown_tier_recoded +
                  (1|total_cases),
                data = CL_invasions_model_between,
                family = nbinom2))

summary(glmmTMB(branching_count ~ lockdown_tier + (1|total_cases),
                data = CL_invasions_branching,
                family = nbinom2))


# Plot distribution of viral moivements over time, mark Mobility Pass date
ggplot(CL_invasions |> filter(head_comuna != tail_comuna)) +
  geom_histogram(aes(x = tail_date, fill = as.factor(variant))) +
  geom_vline(xintercept = as.Date("2021-05-26"), linetype = "dashed") +
  labs(x = "Importation into new comuna estimated date",
       y = "Count", fill = "Variant") +
  scale_fill_manual(values = vocs_colors) +
  theme_minimal()

ggplot(CL_invasions |> filter(head_comuna == tail_comuna)) +
  geom_histogram(aes(x = tail_date, fill = as.factor(variant))) +
  geom_vline(xintercept = as.Date("2021-05-26"), linetype = "dashed") +
  labs(x = "Within-comuna viral movement estimated date",
       y = "Count", fill = "Variant") +
  scale_fill_manual(values = vocs_colors) +
  theme_minimal()

