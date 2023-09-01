################################################################################
########## SARS-CoV-2 Chile phylogenetic analysis metadata creator #############
################################################################################

########################### Bernardo Gutierrez #################################

library(coda)
library(ggmcmc)
library(tidyverse)
library(janitor)
library(patchwork)
library(NatParksPalettes)

########################### Renaming function ##################################

#### Rename parameters - airport v community
transitions.renaming <- function(X) {X |>
  mutate(Count = case_when(Parameter == 'c_airport_background.count[1]' ~
                             'Background to Airport',
                           Parameter == 'c_airport_community.count[1]' ~
                             'Community to Airport',
                           Parameter == 'c_background_airport.count[1]' ~
                             'Airport to Background',
                           Parameter == 'c_background_community.count[1]' ~
                             'Community to Background',
                           Parameter == 'c_community_airport.count[1]' ~
                             'Airport to Community',
                           Parameter == 'c_community_background.count[1]' ~
                             'Background to Community')) |>
  mutate(Description = case_when(Parameter == 'c_airport_background.count[1]' ~
                                   'Imports to airport',
                                 Parameter == 'c_airport_community.count[1]' ~
                                   'Community/Airport',
                                 Parameter == 'c_background_airport.count[1]' ~
                                   'Airport/Background',
                                 Parameter == 'c_background_community.count[1]' ~
                                   'Exports',
                                 Parameter == 'c_community_airport.count[1]' ~
                                   'Imports through airport to community',
                                 Parameter == 'c_community_background.count[1]' ~
                                   'Imports directly to community'))
}

transitions.renaming.2 <- function(X) {X |>
    mutate(Count = case_when(Parameter == 'c_import.count[1]' ~
                               'Introductions',
                             Parameter == 'c_export.count[1]' ~
                               'Exports'))
}

#### Read BEAST log - airport v community
read_BEAST_counts <- function(x){
  read.table(x, sep = "\t") |>
    row_to_names(row_number = 1) |>
    mutate_if(is.character,as.numeric) |>
    select(`c_airport_community.count[1]`, `c_airport_background.count[1]`,
           `c_background_community.count[1]`, `c_background_airport.count[1]`,
           `c_community_airport.count[1]`, `c_community_background.count[1]`) |>
    mcmc() |> ggs() |> select(-Chain) |> filter(Iteration < 500) |>
    transitions.renaming() |> mutate(across(c(Count, Description), factor))
}

#### Read BEAST log - importations v exportations
read_BEAST_imports <- function(x){
  read.table(x, sep = "\t") |>
    row_to_names(row_number = 1) |>
    mutate_if(is.character,as.numeric) |>
    select(`c_import.count[1]`, `c_export.count[1]`) |>
    mcmc() |> ggs() |> select(-Chain) |> filter(Iteration < 500) |>
    transitions.renaming.2() |> mutate(across(Count, factor))
}

################ Metadata files from filtered phylogenies ######################
#### Import BEAST DTA log files
# Airport vs Community transitions
alpha_dta_avc <- 
  read_BEAST_counts("Phylogenetics/DTA_Alpha_ISP/DTA_SARS2_CL_Alpha_ISP.log")

gamma_dta_avc <-
  read_BEAST_counts("Phylogenetics/DTA_Gamma_ISP/DTA_SARS2_CL_Gamma_ISP.log")

lambda_dta_avc <-
  read_BEAST_counts("Phylogenetics/DTA_Lambda_ISP/DTA_SARS2_CL_Lambda_ISP.log")

mu_dta_avc <-
  read_BEAST_counts("Phylogenetics/DTA_Mu_ISP/DTA_SARS2_CL_Mu_ISP.log")

delta_dta_avc <-
  read_BEAST_counts("Phylogenetics/DTA_Delta_ISP/DTA_SARS2_CL_Delta_ISP.log")


# Imports vs Exports
alpha_dta_ive <- 
  read_BEAST_imports("Phylogenetics/DTA_Alpha_ISP/DTA_SARS2_CL_Alpha_ISP.log")

gamma_dta_ive <-
  read_BEAST_imports("Phylogenetics/DTA_Gamma_ISP/DTA_SARS2_CL_Gamma_ISP.log")

lambda_dta_ive <-
  read_BEAST_imports("Phylogenetics/DTA_Lambda_ISP/DTA_SARS2_CL_Lambda_ISP.log")

mu_dta_ive <-
  read_BEAST_imports("Phylogenetics/DTA_Mu_ISP/DTA_SARS2_CL_Mu_ISP.log")

delta_dta_ive <-
  read_BEAST_imports("Phylogenetics/DTA_Delta_ISP/DTA_SARS2_CL_Delta_ISP.log")


################################# Plots ########################################
#### Plot objects for posterior density of Markov Jumps
# Airport v Community
alpha_avc <- ggplot(alpha_dta_avc[alpha_dta_avc$Count== "Background to Airport"|
                                alpha_dta_avc$Count== "Background to Community"|
                                alpha_dta_avc$Count== "Airport to Community",],
                    aes(x = fct_relevel(Count, "Background to Airport",
                                        "Airport to Community"),
                        y = value, fill = Count)) +
  geom_violin(scale = "width") + geom_boxplot(fill = "white", width = 0.1,
                                              outlier.shape = NA) +
  theme_minimal() + scale_fill_manual(values = natparks.pals("Yellowstone", 3)) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8), axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, linewidth = 0.8))

gamma_avc <- ggplot(gamma_dta_avc[gamma_dta_avc$Count== "Background to Airport"|
                                gamma_dta_avc$Count== "Background to Community"|
                                gamma_dta_avc$Count== "Airport to Community",],
                    aes(x = fct_relevel(Count, "Background to Airport",
                                        "Airport to Community"),
                        y = value, fill = Count)) +
  geom_violin(scale = "width") + geom_boxplot(fill = "white", width = 0.1,
                                              outlier.shape = NA) +
  theme_minimal() + scale_fill_manual(values = natparks.pals("Yellowstone", 3)) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8), axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, linewidth = 0.8))

lambda_avc <- ggplot(lambda_dta_avc[lambda_dta_avc$Count== "Background to Airport"|
                                  lambda_dta_avc$Count== "Background to Community"|
                                  lambda_dta_avc$Count== "Airport to Community",],
                     aes(x = fct_relevel(Count, "Background to Airport",
                                         "Airport to Community"),
                         y = value, fill = Count)) +
  geom_violin(scale = "width") + geom_boxplot(fill = "white", width = 0.1,
                                              outlier.shape = NA) +
  theme_minimal() + scale_fill_manual(values = natparks.pals("Yellowstone", 3)) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8), axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, linewidth = 0.8))

mu_avc <- ggplot(mu_dta_avc[mu_dta_avc$Count== "Background to Airport"|
                          mu_dta_avc$Count== "Background to Community"|
                          mu_dta_avc$Count== "Airport to Community",],
                 aes(x = fct_relevel(Count, "Background to Airport",
                                     "Airport to Community"),
                     y = value, fill = Count)) +
  geom_violin(scale = "width") + geom_boxplot(fill = "white", width = 0.1,
                                              outlier.shape = NA) +
  theme_minimal() + scale_fill_manual(values = natparks.pals("Yellowstone", 3)) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8), axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, linewidth = 0.8))


delta_avc <- ggplot(delta_dta_avc[delta_dta_avc$Count== "Background to Airport"|
                                delta_dta_avc$Count== "Background to Community"|
                                delta_dta_avc$Count== "Airport to Community",],
                    aes(x = fct_relevel(Count, "Background to Airport",
                                        "Airport to Community"),
                        y = value, fill = Count)) +
  geom_violin(scale = "width") + geom_boxplot(fill = "white", width = 0.1,
                                              outlier.shape = NA) +
  theme_minimal() + scale_fill_manual(values = natparks.pals("Yellowstone", 3)) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.x = element_text(size = 8),
        panel.border = element_rect(colour = "black",
                                    fill = NA, linewidth = 0.8))

improutes <- alpha_avc / gamma_avc / lambda_avc / mu_avc / delta_avc

ggsave("Figures/phylogeo_importroutes_violinplots.png", height = 10,
       width = 5, dpi = 300, units = "in")
ggsave("Figures/phylogeo_importroutes_violinplots.pdf", height = 10,
       width = 5, dpi = 300, units = "in")

# Imports per VOC
imports <- data.frame(Alpha = alpha_dta_ive$value[alpha_dta_ive$Count=="Introductions"],
                      Gamma = gamma_dta_ive$value[gamma_dta_ive$Count=="Introductions"],
                      Lambda = lambda_dta_ive$value[lambda_dta_ive$Count=="Introductions"],
                      Mu = mu_dta_ive$value[mu_dta_ive$Count=="Introductions"],
                      Delta = delta_dta_ive$value[delta_dta_ive$Count=="Introductions"]) |>
  gather(key = "VOC", value = "Introductions")

impvocs <- ggplot(imports, aes(x = fct_relevel(VOC, "Alpha", "Gamma", "Lambda", "Mu"),
                               y = Introductions, fill = VOC)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_boxplot(fill = "white", width = 0.2, linewidth = 0.2, outlier.shape = NA) +
  theme_minimal() + scale_fill_manual(values = vocs_colors) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10))

ggsave("Figures/phylogeo_importVOCs_violinplots.png", height = 3.5,
       width = 6, dpi = 300, units = "in", bg = "white")

ggsave("Figures/phylogeo_importVOCs_violinplots.pdf", height = 3.5,
       width = 6, dpi = 300, units = "in", bg = "white")
