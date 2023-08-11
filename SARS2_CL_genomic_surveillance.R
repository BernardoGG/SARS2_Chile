################################################################################
########### SARS-CoV-2 genomic surveillance in Chile ###########################
################################################################################

#################### Bernardo Gutierrez ########################################

# Source scripts and data
source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")

#### Plot variants from various surveillance contexts ###
## Contextualize ISP databases versus available genomic data on GISAID ####
# GISAID data
lins_cl_gisaid <- GISAID_CL %>%
  arrange(epiweek, Variant) %>%
  group_by(epiweek, Variant) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count))


gisaid_panel <- ggplot(data = lins_cl_gisaid , aes(x = epiweek, y = percentage, fill = as.factor(Variant))) +
  geom_bar(stat = "identity", colour = "white") +
  labs(x = element_blank(),
       y = "GISAID",
       fill = "Lineage / VOC") + scale_fill_manual(values = vocs_colors) + theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# ISP surveillance databases
lins_cl <- chile_metadata %>%
  arrange(epiweek, Variant) %>%
  group_by(epiweek, Variant) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count))


isp_panel <- ggplot(data = lins_cl , aes(x = epiweek, y = percentage, fill = as.factor(Variant))) +
  geom_bar(stat = "identity", colour = "white") +
  labs(x = "Epidemiological week",
       y = "ISP database",
       fill = "Lineage / VOC") + scale_fill_manual(values = vocs_colors) + theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

gisaid_panel / isp_panel


## Compare airport versus community surveillance data sets
# Community surveillance
lins_comm <- chile_isp_metadata[chile_isp_metadata$surveillance == "Community",] %>%
  arrange(epiweek_start, Variant) %>%
  group_by(epiweek_start, Variant) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count))


comm_panel <- ggplot(data = lins_comm , aes(x = epiweek_start, y = percentage,
                                            fill = as.factor(Variant))) +
  geom_bar(stat = "identity", colour = "white", show.legend = FALSE) +
  labs(x = element_blank(),
       y = "Community") + scale_fill_manual(values = vocs_colors) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


# Airport surveillance
lins_airp <- chile_isp_metadata[chile_isp_metadata$surveillance == "Airport",] %>%
  arrange(epiweek_start, Variant) %>%
  group_by(epiweek_start, Variant) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count))

airp_panel <- ggplot(data = lins_airp , aes(x = epiweek_start, y = percentage,
                                            fill = as.factor(Variant))) +
  geom_bar(stat = "identity", colour = "white") +
  labs(x = "Epidemiological week start date",
       y = "Airport",
       fill = "VOC") + scale_fill_manual(values = vocs_colors) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

surv_props <- comm_panel / airp_panel +
  plot_layout(guides = 'collect')

ggsave("Figures/CL_surveillance_props.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")

ggsave("Figures/CL_surveillance_props.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")

### Evaluate sampling intensity across time ####
# Create data frame containing number of sequences per epiweek per data set
surv_intensity_epiweeks <- full_join(as.data.frame(table(GISAID_CL$epiweek)),
                                     as.data.frame(table(chile_metadata$epiweek)),
                                     by = "Var1")
colnames(surv_intensity_epiweeks) <- c("Epiweek", "GISAID", "ISP")

gisaid_vars <- group_by(GISAID_CL, epiweek, Variant) %>%
  summarize(count = n())
isp_vars <- group_by(chile_metadata, epiweek, Variant) %>%
  summarize(count = n())
surv_intensity_epiweeks_vars <- full_join(gisaid_vars, isp_vars,
                                          by = c("epiweek", "Variant"))
colnames(surv_intensity_epiweeks_vars) <- c("Epiweek", "Variant", "GISAID", "ISP")
surv_intensity_epiweeks_vars$ISP[is.na(surv_intensity_epiweeks_vars$ISP)] <- 0

# Epiweeks plots of coverage between data sets
epiweeks_plot <- surv_intensity_epiweeks %>%
  pivot_longer(c(ISP, GISAID), "Surveillance")

ggplot(data = epiweeks_plot,
       aes(x = Epiweek, y = value, fill = as.factor(Surveillance))) +
  geom_bar(stat = "identity", colour = "white") +
  labs(x = element_blank(),
       y = "No. of sequences",
       fill = "ISP Surveillance") +
  scale_fill_manual(values = c("#DAD2D8", "#126168")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Epiweeks plots of coverage between data sets per variant
epiweeks_plot_vars <- surv_intensity_epiweeks_vars %>%
  pivot_longer(c(ISP, GISAID), "Surveillance")

v_plot <- "Mu"
ggplot(data = epiweeks_plot_vars[epiweeks_plot_vars$Variant==v_plot,],
       aes(x = Epiweek, y = value, fill = as.factor(Surveillance))) +
  geom_bar(stat = "identity", colour = "white") +
  xlim(54, 94) +
  labs(title = v_plot,
       x = element_blank(),
       y = "No. of sequences",
       fill = "ISP Surveillance") +
  scale_fill_manual(values = c("#DAD2D8", "#126168")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Sampling over time per variant
ggplot(data = epiweeks_plot_vars,
       aes(x = Epiweek, y = value, color = Variant)) +
  geom_line() +
  xlim(54, 94) +
  labs(x = element_blank(),
       y = "No. of sequences",
       color = "VOI/VOC") +
  scale_color_manual(values = vocs_colors) +
  facet_wrap(~Surveillance, ncol = 1) +
  theme_minimal()


### Evaluate sampling intensity across regions ####
# Create data frame containing number of sequences per region per data set
surv_intensity_cl <- full_join(as.data.frame(table(GISAID_CL$State)),
                               as.data.frame(table(chile_metadata$State)),
                               by = "Var1") %>% slice(-14)
colnames(surv_intensity_cl) <- c("Region", "GISAID", "ISP")
surv_intensity_cl$ISP[is.na(surv_intensity_cl$ISP)] <- 0

# Semi-log scales plot
plot_ly(surv_intensity_cl, x = ~GISAID, y = ~ISP) %>% add_markers() %>%
  layout(fig, xaxis = list(type = "log"), yaxis = list(type = "log"))

# Untransformed scales plot
surv_plot <- surv_intensity_cl %>%
  pivot_longer(c(ISP, GISAID), "Surveillance")

ggplot(data = surv_plot , aes(x = Region, y = value, fill = as.factor(Surveillance))) +
  geom_bar(stat = "identity", colour = "white") +
  labs(x = element_blank(),
       y = "No. of sequences",
       fill = "ISP Surveillance") +
  scale_fill_manual(values = c("#DAD2D8", "#126168")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Plot transmission lineage summaries #####
## Importation intensity proxies (TMRCAs) across time
ggplot(tl_nosingletons, aes(x = tmrca, fill = earliest_taxa_state_comp)) +
  geom_histogram(bins = (as.numeric(max(tl_nosingletons$tmrca) -
                                      min(tl_nosingletons$tmrca)))/7) +
  labs(x = "Transmission lineage TMRCA date",
       y = "Introductions per week",
       fill = "Earliest sequence detection") +
  scale_fill_manual(values = c("#EC9A29", "#126168", "#DAD2D8")) +
  theme_minimal()
