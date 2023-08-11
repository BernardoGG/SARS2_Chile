################################################################################
################## SARS-CoV-2 VOCs in South America plots ######################
################################################################################

########################### Bernardo Gutierrez #################################

# Source scripts and data
#source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")

# Packages
library(tidyverse)
library(janitor)
library(devtools)
#install_github("cran/GADMTools")
library(GADMTools)
library(patchwork)
library(NatParksPalettes)


####################### Shapefiles for GIS mapping #############################
## Import maps for countries in the region
southamerica <- gadm_sp_loadCountries(fileNames = c('ARG', 'BOL', 'BRA', 'CHL',
                                                    'COL', 'ECU', 'GUF', 'GUY',
                                                    'PAN', 'PRY', 'PER', 'SUR',
                                                    'URY', 'VEN'),
                                          level = 0, basefile = './shapefiles',
                                      simplify = 0.1)

gadm_sa <- fortify(southamerica[[2]])
gadm_sa$country <- as.factor(c(rep("Argentina", length(grep("ARG_", gadm_sa$id))),
                               rep("Bolivia", length(grep("BOL_", gadm_sa$id))),
                               rep("Brazil", length(grep("BRA_", gadm_sa$id))),
                               rep("Chile", length(grep("CHL_", gadm_sa$id))),
                               rep("Colombia", length(grep("COL_", gadm_sa$id))),
                               rep("Ecuador", length(grep("ECU_", gadm_sa$id))),
                               rep("", length(grep("GUF_", gadm_sa$id))),
                               rep("", length(grep("GUY_", gadm_sa$id))),
                               rep("", length(grep("PAN_", gadm_sa$id))),
                               rep("Paraguay", length(grep("PRY_", gadm_sa$id))),
                               rep("Peru", length(grep("PER_", gadm_sa$id))),
                               rep("", length(grep("SUR_", gadm_sa$id))),
                               rep("Uruguay", length(grep("URY_", gadm_sa$id))),
                               rep("Venezuela", length(grep("VEN_", gadm_sa$id)))))

############################## GISAID data #####################################
# Filter global metadata to match Chile metadata dates and
# include column specifying VOI/VOC variant designation
GISAID_SA <- GISAID_final[GISAID_final$Country=="Argentina" |
                            GISAID_final$Country=="Bolivia" |
                            GISAID_final$Country=="Brazil" |
                            GISAID_final$Country=="Chile" |
                            GISAID_final$Country=="Colombia" |
                            GISAID_final$Country=="Ecuador" |
                            GISAID_final$Country=="Paraguay" |
                            GISAID_final$Country=="Peru" |
                            GISAID_final$Country=="Uruguay" |
                            GISAID_final$Country=="Venezuela",] %>%
  mutate(Variant = case_when(`Pango lineage` == "B.1.1.7" |
                               `Pango lineage` == "Q.1" ~ "Alpha",
                             `Pango lineage` %in%
                               unique(GISAID_final$`Pango lineage`[
                                 grep("AY.", GISAID_final$`Pango lineage`)]) ~
                               "Delta",
                             `Pango lineage` %in%
                               unique(GISAID_final$`Pango lineage`[
                                 grep("P.1", GISAID_final$`Pango lineage`)]) ~
                               "Gamma",
                             `Pango lineage` == "C.37" |
                               `Pango lineage` == "C.37.1" ~ "Lambda",
                             `Pango lineage` == "B.1.621" |
                               `Pango lineage` == "B.1.621.1" |
                               `Pango lineage` == "B.1.621.2" |
                               `Pango lineage` == "BB.2" ~ "Mu",
                             `Pango lineage` == "B.1.526"  ~ "Iota",
                             TRUE ~ "Other"))

# Generate counts of sequences per ~73-day epochs per country
VOCs_T1 <- as.data.frame(table(GISAID_SA$Variant[
  GISAID_SA$Date < as.Date("2021-03-06")],
  GISAID_SA$Country[GISAID_SA$Date < as.Date("2021-03-06")]))
colnames(VOCs_T1) <- c("VOC", "Country", "Count")

VOCs_T2 <- as.data.frame(table(GISAID_SA$Variant[
  GISAID_SA$Date < as.Date("2021-05-18") &
    GISAID_SA$Date >= as.Date("2021-03-06")],
  GISAID_SA$Country[GISAID_SA$Date < as.Date("2021-05-18") &
                      GISAID_SA$Date >= as.Date("2021-03-06")]))
colnames(VOCs_T2) <- c("VOC", "Country", "Count")

VOCs_T3 <- as.data.frame(table(GISAID_SA$Variant[
  GISAID_SA$Date < as.Date("2021-07-30") &
    GISAID_SA$Date >= as.Date("2021-05-18")],
  GISAID_SA$Country[GISAID_SA$Date < as.Date("2021-07-30") &
                      GISAID_SA$Date >= as.Date("2021-05-18")]))
colnames(VOCs_T3) <- c("VOC", "Country", "Count")

VOCs_T4 <- as.data.frame(table(GISAID_SA$Variant[
  GISAID_SA$Date >= as.Date("2021-07-30")],
  GISAID_SA$Country[GISAID_SA$Date >= as.Date("2021-07-30")]))
colnames(VOCs_T4) <- c("VOC", "Country", "Count")

# Add column with total number of sequences per country for that epoch
# *NOTE*: 'Tot' is sensitive to order of functions being executed
Total_T1 <- VOCs_T1 %>% group_by(Country) %>% summarise(sum(Count))
Total_T2 <- VOCs_T2 %>% group_by(Country) %>% summarise(sum(Count))
Total_T3 <- VOCs_T3 %>% group_by(Country) %>% summarise(sum(Count))
Total_T4 <- VOCs_T4 %>% group_by(Country) %>% summarise(sum(Count))

Tot <- vector()
for(i in 1:nrow(Total_T1)){
  x <- rep(Total_T1$`sum(Count)`[i], 6)
  Tot <- c(Tot, x)}
VOCs_T1$Total <- Tot

Tot <- vector()
for(i in 1:nrow(Total_T2)){
  x <- rep(Total_T2$`sum(Count)`[i], 6)
  Tot <- c(Tot, x)}
VOCs_T2$Total <- Tot

Tot <- vector()
for(i in 1:nrow(Total_T3)){
  x <- rep(Total_T3$`sum(Count)`[i], 6)
  Tot <- c(Tot, x)}
VOCs_T3$Total <- Tot

Tot <- vector()
for(i in 1:nrow(Total_T4)){
  x <- rep(Total_T4$`sum(Count)`[i], 6)
  Tot <- c(Tot, x)}
VOCs_T4$Total <- Tot

# Create data frame with VOC proportions over four epochs
voc_props <- data.frame(
  country = unique(GISAID_SA$Country) %>% sort(),
  gamma_t1 = VOCs_T1[VOCs_T1$VOC=="Gamma",]$Count/
    VOCs_T1[VOCs_T1$VOC=="Gamma",]$Total,
  gamma_t2 = VOCs_T2[VOCs_T2$VOC=="Gamma",]$Count/
    VOCs_T2[VOCs_T2$VOC=="Gamma",]$Total,
  gamma_t3 = VOCs_T3[VOCs_T3$VOC=="Gamma",]$Count/
    VOCs_T3[VOCs_T3$VOC=="Gamma",]$Total,
  gamma_t4 = VOCs_T4[VOCs_T4$VOC=="Gamma",]$Count/
    VOCs_T4[VOCs_T4$VOC=="Gamma",]$Total,
  lambda_t1 = VOCs_T1[VOCs_T1$VOC=="Lambda",]$Count/
    VOCs_T1[VOCs_T1$VOC=="Lambda",]$Total,
  lambda_t2 = VOCs_T2[VOCs_T2$VOC=="Lambda",]$Count/
    VOCs_T2[VOCs_T2$VOC=="Lambda",]$Total,
  lambda_t3 = VOCs_T3[VOCs_T3$VOC=="Lambda",]$Count/
    VOCs_T3[VOCs_T3$VOC=="Lambda",]$Total,
  lambda_t4 = VOCs_T4[VOCs_T4$VOC=="Lambda",]$Count/
    VOCs_T4[VOCs_T4$VOC=="Lambda",]$Total,
  mu_t1 = VOCs_T1[VOCs_T1$VOC=="Mu",]$Count/
    VOCs_T1[VOCs_T1$VOC=="Mu",]$Total,
  mu_t2 = VOCs_T2[VOCs_T2$VOC=="Mu",]$Count/
    VOCs_T2[VOCs_T2$VOC=="Mu",]$Total,
  mu_t3 = VOCs_T3[VOCs_T3$VOC=="Mu",]$Count/
    VOCs_T3[VOCs_T3$VOC=="Mu",]$Total,
  mu_t4 = VOCs_T4[VOCs_T4$VOC=="Mu",]$Count/
    VOCs_T4[VOCs_T4$VOC=="Mu",]$Total,
  delta_t1 = VOCs_T1[VOCs_T1$VOC=="Delta",]$Count/
    VOCs_T1[VOCs_T1$VOC=="Delta",]$Total,
  delta_t2 = VOCs_T2[VOCs_T2$VOC=="Delta",]$Count/
    VOCs_T2[VOCs_T2$VOC=="Delta",]$Total,
  delta_t3 = VOCs_T3[VOCs_T3$VOC=="Delta",]$Count/
    VOCs_T3[VOCs_T3$VOC=="Delta",]$Total,
  delta_t4 = VOCs_T4[VOCs_T4$VOC=="Delta",]$Count/
    VOCs_T4[VOCs_T4$VOC=="Delta",]$Total)


########################## VOC proportion maps #################################
# Add columns with proportions of individual VOCs to map file
gadm_sa$gamma_t1 <- voc_props$gamma_t1[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$gamma_t2 <- voc_props$gamma_t2[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$gamma_t3 <- voc_props$gamma_t3[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$gamma_t4 <- voc_props$gamma_t4[match(unlist(gadm_sa$country),
                                             voc_props$country)]

gadm_sa$lambda_t1 <- voc_props$lambda_t1[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$lambda_t2 <- voc_props$lambda_t2[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$lambda_t3 <- voc_props$lambda_t3[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$lambda_t4 <- voc_props$lambda_t4[match(unlist(gadm_sa$country),
                                             voc_props$country)]

gadm_sa$mu_t1 <- voc_props$mu_t1[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$mu_t2 <- voc_props$mu_t2[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$mu_t3 <- voc_props$mu_t3[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$mu_t4 <- voc_props$mu_t4[match(unlist(gadm_sa$country),
                                             voc_props$country)]

gadm_sa$delta_t1 <- voc_props$delta_t1[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$delta_t2 <- voc_props$delta_t2[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$delta_t3 <- voc_props$delta_t3[match(unlist(gadm_sa$country),
                                             voc_props$country)]
gadm_sa$delta_t4 <- voc_props$delta_t4[match(unlist(gadm_sa$country),
                                             voc_props$country)]

# Plot maps for proportions of each VOC per epoch
g1 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = gamma_t1),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#EC9A29", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  labs(y = "Gamma") +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

g2 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = gamma_t2),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#EC9A29", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

g3 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = gamma_t3),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#EC9A29", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

g4 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = gamma_t4),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#EC9A29", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank())


l1 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = lambda_t1),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#0F8B8D", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  labs(y = "Lambda") +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

l2 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = lambda_t2),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#0F8B8D", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

l3 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = lambda_t3),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#0F8B8D", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

l4 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = lambda_t4),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#0F8B8D", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank())


m1 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = mu_t1),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#143642", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  labs(y = "Mu") +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

m2 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = mu_t2),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#143642", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

m3 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = mu_t3),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#143642", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

m4 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = mu_t4),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#143642", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank())


d1 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = delta_t1),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#CA5D22", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  labs(x = "2020-12-22 to\n2021-03-06", y = "Delta") +
  theme(line = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

d2 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = delta_t2),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#CA5D22", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  labs(x = "2021-03-07 to\n2021-05-18") +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

d3 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = delta_t3),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#CA5D22", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  labs(x = "2021-05-19 to\n2021-07-30") +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

d4 <- ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = delta_t4),
               color = 'grey20', size = 0.25) +
  scale_fill_gradient(low = 'white', high = "#CA5D22", limits = c(0, 1)) +
  annotate("rect", xmin = -83, xmax = -78, ymin = -40, ymax = -20,
           color = "white", fill = "white") +
  annotate("rect", xmin = -82, xmax = -80, ymin = 3, ymax = 4,
           color = "white", fill = "white") +
  annotate("rect", xmin = -83, xmax = -78, ymin = 10.5, ymax = 12.5,
           color = "white", fill = "white") +
  coord_fixed() + theme_minimal() + xlim(c(-84, -34)) + ylim(-60, 12.5) +
  labs(x = "2021-07-31 to\n2021-10-12") +
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
  

# Plot all maps
plot <- (g1 | g2 | g3 | g4) / (l1 | l2 | l3 | l4) /
  (m1 | m2 | m3 | m4) /  (d1 | d2 | d3 | d4)

ggsave(plot = plot, "Figures/SouthAm_VOCs_map.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")

ggsave(plot = plot, "Figures/SouthAm_VOCs_map.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")


####################### VOC proportion timelines ###############################
# Generate data frame counting sequences by country and epiweek
lins_sa_gisaid <- GISAID_SA %>%
  arrange(epiweek, Variant) %>%
  group_by(Country, epiweek, Variant) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count))

# Plot proportions
plot <- ggplot(data = lins_sa_gisaid ,
               aes(x = epiweek, y = percentage, fill = as.factor(Variant))) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  labs(x = "GISAID epiweek",
       y = element_blank(),
       fill = "Variant") +
  scale_fill_manual(values = c(vocs_colors[1:3], "#7b6d8d", vocs_colors[4:6])) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(vars(Country), 5, 2)

ggsave(plot = plot, "Figures/SouthAm_VOCs_timeline.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")

ggsave(plot = plot, "Figures/SouthAm_VOCs_timeline.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")


####### Sandbox ########

ggplot() +
  geom_polygon(data = gadm_sa,
               aes(x = long, y = lat, group = group, fill = country),
               color = 'grey20', size = 0.4) +
  coord_fixed() + xlim(c(-84, -34)) + ylim(-60, 12.5) + theme_void()
