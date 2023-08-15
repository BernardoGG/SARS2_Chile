################################################################################
################ SARS-CoV-2 Chile lockdown effects modelling ###################
################################################################################

# Source scripts, packages and data
source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")
library(MASS)
library(rcompanion)
library(tidyverse)
library(lubridate)
library(patchwork)

########################### Bernardo Gutierrez #################################
### Set data frames ####
## Add column with numbers of new cases per report date
CL_invasions <- CL_invasions %>%
  left_join(comuna_cases, by = c("head_comuna" = "Comuna",
                              "head_date" = "Fecha")) %>%
  select(-`Casos Confirmados`, head_date_new_cases = `Casos nuevos`)

CL_invasions <- CL_invasions %>%
  left_join(comuna_cases, by = c("tail_comuna" = "Comuna",
                                 "tail_date" = "Fecha")) %>%
  select(-`Casos Confirmados`, tail_date_new_cases = `Casos nuevos`)


## Internal movement and importation tables
# Include counts and total new cases per comuna per lockdown tier
CL_invasions_MEM_within <- CL_invasions %>%
  filter(head_comuna == tail_comuna) %>%
  group_by(tail_comuna, tail_lockdown_tier_recoded) %>%
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) %>%
  as.data.frame() %>%
  mutate(count = as.numeric(count)) %>%
  mutate(tail_comuna = as.factor(tail_comuna)) %>%
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_MEM_between <- CL_invasions %>%
  filter(head_comuna != tail_comuna) %>%
  group_by(tail_comuna, tail_lockdown_tier_recoded) %>%
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) %>%
  as.data.frame() %>%
  mutate(count = as.numeric(count)) %>%
  mutate(tail_comuna = as.factor(tail_comuna)) %>%
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

# Include counts, total new cases per comuna per lockdown tier and dates
CL_invasions_MEM_time_within <- CL_invasions %>%
  filter(head_comuna == tail_comuna) %>%
  group_by(tail_comuna, tail_lockdown_tier_recoded, tail_date) %>%
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) %>%
  as.data.frame() %>%
  mutate(count = as.numeric(count)) %>%
  mutate(tail_comuna = as.factor(tail_comuna)) %>%
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))

CL_invasions_MEM_time_between <- CL_invasions %>%
  filter(head_comuna != tail_comuna) %>%
  group_by(tail_comuna, tail_lockdown_tier_recoded, tail_date) %>%
  summarise(count = n(),
            total_cases = sum(tail_date_new_cases, na.rm = TRUE)) %>%
  as.data.frame() %>%
  mutate(count = as.numeric(count)) %>%
  mutate(tail_comuna = as.factor(tail_comuna)) %>%
  mutate(norm_count = as.numeric(count / norm_comuna_tail[
    as.character(tail_comuna)]))


### Models and statistical testing ####
## Model 1: LM of within-comuna movements, accounting for new cases
m1_within <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
        data = CL_invasions_MEM_within)

# Observed vs predicted values of within-comuna movements
m1_within_pred <- predict(m1_within)
plot(m1_within_pred, CL_invasions_MEM_within$count)

## Model 1: LM of cross-comuna movements, accounting for new cases
m1_between <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
               data = CL_invasions_MEM_between)

# Observed vs predicted values of cross-comuna movements
m1_between_pred <- predict(m1_between)
plot(m1_between_pred, CL_invasions_MEM_between$count)


## Model 2: Poisson GLM of within-comuna movements, accounting for new cases
m2_within <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
                 data = CL_invasions_MEM_within,
          family = "poisson")

# Observed vs predicted values of within-comuna movements
m2_within_pred <- predict(m2_within)
plot(m2_within_pred, CL_invasions_MEM_within$count)


## Model 2: Poisson GLM of cross-comuna movements, accounting for new cases
m2_between <- glm(count ~ tail_lockdown_tier_recoded + total_cases,
                 data = CL_invasions_MEM_between,
                 family = "poisson")

# Observed vs predicted values of cross-comuna movements
m2_between_pred <- predict(m2_between)
plot(m2_between_pred, CL_invasions_MEM_between$count)


## Model 3: Negative binomial GLM of within-comuna movements, accounting for new cases and date
m3_within <- glm.nb(count ~ tail_lockdown_tier_recoded + total_cases,
                     data = CL_invasions_MEM_time_within)

# Observed vs predicted values of cross-comuna movements
m3_within_pred <- predict(m3_within)
plot(m3_within_pred, CL_invasions_MEM_time_within$count)


## Model 3: Negative binomial GLM of cross-comuna movements, accounting for new cases and date
m3_between <- glm.nb(count ~ tail_lockdown_tier_recoded + total_cases,
                  data = CL_invasions_MEM_time_between)

# Observed vs predicted values of cross-comuna movements
m3_between_pred <- predict(m3_between)
plot(m3_between_pred, CL_invasions_MEM_time_between$count)


### Model test and comparisons

## Within comunas
summary(m1_within)
summary(m2_within)
summary(m3_within)

# Compare models under different distn families
compareGLM(m1_within, m2_within, m3_within)

## Between comunas
summary(m1_between)
summary(m2_between)
summary(m3_between)

# Compare models under different distn families
compareGLM(m1_between, m2_between, m3_between)

## LRT between Poisson and NB models to test goodness of fit based on dispersion
## parameter
pchisq(2 * (logLik(m3_within) - logLik(m2_within)), df = 1, lower.tail = FALSE)
pchisq(2 * (logLik(m3_between) - logLik(m2_between)), df = 1, lower.tail = FALSE)

(est_within <- cbind(Estimate = coef(m3_within), confint(m3_within)))
exp(est_within)
(est_between <- cbind(Estimate = coef(m3_between), confint(m3_between)))
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
m1_between_resids <- data.frame(predicted_values = m1_between_pred,
                               deviance_residuals = resid(m1_between))
m2_between_resids <- data.frame(predicted_values = m2_between_pred,
                               deviance_residuals = resid(m2_between))
m3_between_resids <- data.frame(predicted_values = m3_between_pred,
                               deviance_residuals = resid(m3_between))

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

## Test effects of lockdown tier by comparing model which includes and excludes
## this variable - applicable for model that best fits the data (NB)
m3_within_noLD <- update(m3_within, . ~ . - tail_lockdown_tier_recoded)
m3_between_noLD <- update(m3_between, . ~ . - tail_lockdown_tier_recoded)

anova(m3_within, m3_within_noLD)
anova(m3_between, m3_between_noLD)

## Get NB estimates for data generated under that model holding numbers of cases
## at their mean.
new_m3_within <-
  data.frame(total_cases = mean(CL_invasions_MEM_within$total_cases),
             tail_lockdown_tier_recoded = factor(1:3, levels = 1:3,
                                    labels = levels(
                                      CL_invasions_MEM_within$tail_lockdown_tier_recoded)))
new_m3_within$phat <- predict(m3_within, new_m3_within, type = "response")
new_m3_within

new_m3_between <-
  data.frame(total_cases = mean(CL_invasions_MEM_between$total_cases),
             tail_lockdown_tier_recoded = factor(1:3, levels = 1:3,
                                                 labels = levels(
                                                   CL_invasions_MEM_between$tail_lockdown_tier_recoded)))
new_m3_between$phat <- predict(m3_between, new_m3_between, type = "response")
new_m3_between


## Estimate frequency of branching events under different lockdown tiers
# Counts of unique occurrences of comunas as 'head' nodes
# Corresponds to all branching event taking place per comuna
ancestor_total_counts <- CL_invasions %>%
  group_by(head_comuna) %>%
  summarise(count = n())

# Counts of unique occurrences of comunas as 'head' nodes with the same 'tail' node
# Corresponds to all branching events resulting in descendants within the same comuna
ancestor_local_counts <- CL_invasions %>%
  filter(head_comuna == tail_comuna) %>%
  group_by(head_comuna) %>%
  summarise(count = n())

# All LD tiers
counts_all_tiers <- CL_invasions %>%
  group_by(head_comuna, tail_comuna) %>%
  summarise(count = n())
ancestor_counts_all_tiers <- counts_all_tiers %>%
  filter(head_comuna == tail_comuna) %>%
  group_by(head_comuna) %>%
  summarise(branching = sum(count, na.rm = TRUE))

# Full lockdown
counts_full_ld <- CL_invasions %>%
  filter(head_lockdown_tier_recoded == "Full lockdown") %>%
  group_by(head_comuna, tail_comuna) %>%
  summarise(count = n())
ancestor_counts_full_ld <- counts_full_ld %>%
  filter(head_comuna == tail_comuna) %>%
  group_by(head_comuna) %>%
  summarise(branching = sum(count, na.rm = TRUE))

# Weekend lockdown
counts_weekend_ld <- CL_invasions %>%
  filter(head_lockdown_tier_recoded == "Weekend lockdown") %>%
  group_by(head_comuna, tail_comuna) %>%
  summarise(count = n())
ancestor_counts_weekend_ld <- counts_weekend_ld %>%
  filter(head_comuna == tail_comuna) %>%
  group_by(head_comuna) %>%
  summarise(branching = sum(count, na.rm = TRUE))

# No lockdown
counts_no_ld <- CL_invasions %>%
  filter(head_lockdown_tier_recoded == "No lockdown") %>%
  group_by(head_comuna, tail_comuna) %>%
  summarise(count = n())
ancestor_counts_no_ld <- counts_no_ld %>%
  filter(head_comuna == tail_comuna) %>%
  group_by(head_comuna) %>%
  summarise(branching = sum(count, na.rm = TRUE))

branching_all_tiers <- left_join(ancestor_local_counts, ancestor_counts_all_tiers) %>%
  mutate(branching = ifelse(is.na(branching), 0, branching)) %>%
  mutate(branching_freq = branching / count)

branching_full_ld <- left_join(ancestor_local_counts, ancestor_counts_full_ld) %>%
  mutate(branching = ifelse(is.na(branching), 0, branching)) %>%
  mutate(branching_freq = branching / count)

branching_weekend_ld <- left_join(ancestor_local_counts, ancestor_counts_weekend_ld) %>%
  mutate(branching = ifelse(is.na(branching), 0, branching)) %>%
  mutate(branching_freq = branching / count)

branching_no_ld <- left_join(ancestor_local_counts, ancestor_counts_no_ld) %>%
  mutate(branching = ifelse(is.na(branching), 0, branching)) %>%
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
    branching_full_ld %>% rename(branching_freq_full = branching_freq) %>%
      select(head_comuna, branching_freq_full),
    branching_weekend_ld %>% rename(branching_freq_wkd = branching_freq) %>%
      select(head_comuna, branching_freq_wkd)
  ),
  branching_no_ld %>% rename(branching_freq_no = branching_freq) %>%
    select(head_comuna, branching_freq_no)
) %>%
  pivot_longer(
    !head_comuna, names_to = "lockdown_tier", values_to = "branching_freq"
  ) %>%
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

### Sandbox ####
## Model 3: LM of within-comuna movements, accounting for new cases and dates
CL_invasions_MEM_time_within$tail_date_min <- min(CL_invasions_MEM_time_within$tail_date)
CL_invasions_MEM_time_within <- CL_invasions_MEM_time_within %>%
  mutate(tail_date_n = as.numeric(tail_date - tail_date_min))

m3_within <-lm(count ~ tail_lockdown_tier_recoded + total_cases + tail_date_n,
               data = CL_invasions_MEM_time_within %>% filter())
summary(m3_within)

  
plot(CL_invasions_MEM_time_within$tail_date, CL_invasions_MEM_time_within$count)

# Observed vs predicted values of within-comuna movements
m3_within_pred <- predict(m3_within)
plot(CL_invasions_MEM_time_within$count, m3_within_pred)

## Model 3: LM of cross-comuna movements, accounting for new cases and dates
m3_between <-lm(count ~ tail_lockdown_tier_recoded + total_cases + tail_date,
                data = CL_invasions_MEM_time_between)
summary(m3_between)

# Observed vs predicted values of cross-comuna movements
m3_between_pred <- predict(m3_between)
plot(CL_invasions_MEM_time_between$count, m3_between_pred)


## Model 4: Poisson GLM of within-comuna movements, accounting for new cases and date
m4_within <- glm(count ~ tail_lockdown_tier_recoded + total_cases + tail_date,
                 data = CL_invasions_MEM_time_within,
                 family = "poisson")
summary(m4_within)
exp(m4_within$coefficients)

# Observed vs predicted values of within-comuna movements
m4_within_pred <- predict(m4_within)
plot(CL_invasions_MEM_time_within$count, m4_within_pred)


## Model 4: Poisson GLM of cross-comuna movements, accounting for new cases and date
m4_between <- glm(count ~ tail_lockdown_tier_recoded + total_cases + tail_date,
                  data = CL_invasions_MEM_time_between,
                  family = "poisson")
summary(m4_between)
exp(m4_between$coefficients)

# Observed vs predicted values of cross-comuna movements
m4_between_pred <- predict(m4_between)
plot(CL_invasions_MEM_time_between$count, m4_between_pred)
