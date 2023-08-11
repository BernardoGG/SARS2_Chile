################################################################################
####### SARS-CoV-2 Chile transmission lineage phylogeo file synthesis ##########
################################################################################

########################### Bernardo Gutierrez #################################

####################### Import TL phylogeo tables #############################
#### Files generated using Python script

phylogeo_alpha13 <- read.csv("Phylogenetics/Transmission_lineage_processing/Alpha_TL_13.tsv",
                      sep = "\t") #TODO sort out folders
phylogeo_delta10 <- read.csv("Phylogenetics/Transmission_lineage_processing/Delta_TL_10.tsv",
                               sep = "\t")
phylogeo_delta64 <- read.csv("Phylogenetics/Transmission_lineage_processing/Delta_TL_64.tsv",
                               sep = "\t")
phylogeo_delta166 <- read.csv("Phylogenetics/Transmission_lineage_processing/Delta_TL_166.tsv",
                               sep = "\t")
phylogeo_delta179 <- read.csv("Phylogenetics/Transmission_lineage_processing/Delta_TL_179.tsv",
                               sep = "\t")
phylogeo_delta186 <- read.csv("Phylogenetics/Transmission_lineage_processing/Delta_TL_186.tsv",
                               sep = "\t")
phylogeo_delta191 <- read.csv("Phylogenetics/Transmission_lineage_processing/Delta_TL_191.tsv",
                               sep = "\t")
phylogeo_delta254 <- read.csv("Phylogenetics/Transmission_lineage_processing/Delta_TL_254.tsv",
                               sep = "\t")
phylogeo_gamma35 <- read.csv("Phylogenetics/Transmission_lineage_processing/Gamma_TL_35.tsv",
                               sep = "\t")
phylogeo_gamma55 <- read.csv("Phylogenetics/Transmission_lineage_processing/Gamma_TL_55.tsv",
                               sep = "\t")
phylogeo_lambda90 <- read.csv("Phylogenetics/Transmission_lineage_processing/Lambda_TL_90.tsv",
                               sep = "\t")
phylogeo_lambda95 <- read.csv("Phylogenetics/Transmission_lineage_processing/Lambda_TL_95.tsv",
                                sep = "\t")
phylogeo_lambda103 <- read.csv("Phylogenetics/Transmission_lineage_processing/Lambda_TL_103.tsv",
                                sep = "\t")
phylogeo_lambda104 <- read.csv("Phylogenetics/Transmission_lineage_processing/Lambda_TL_104.tsv",
                                sep = "\t")
phylogeo_lambda107 <- read.csv("Phylogenetics/Transmission_lineage_processing/Lambda_TL_107.tsv",
                                sep = "\t")
phylogeo_mu15 <- read.csv("Phylogenetics/Transmission_lineage_processing/Mu_TL_15.tsv",
                                sep = "\t")
phylogeo_mu26 <- read.csv("Phylogenetics/Transmission_lineage_processing/Mu_TL_26.tsv",
                            sep = "\t")
phylogeo_mu61 <- read.csv("Phylogenetics/Transmission_lineage_processing/Mu_TL_61.tsv",
                            sep = "\t")
phylogeo_mu73 <- read.csv("Phylogenetics/Transmission_lineage_processing/Mu_TL_73.tsv",
                            sep = "\t")
phylogeo_mu82 <- read.csv("Phylogenetics/Transmission_lineage_processing/Mu_TL_82.tsv",
                            sep = "\t")

###################### Add columns to data frames #############################

phylogeo_alpha13$variant <- rep("Alpha", nrow(phylogeo_alpha13))
phylogeo_delta10$variant <- rep("Delta", nrow(phylogeo_delta10))
phylogeo_delta64$variant <- rep("Delta", nrow(phylogeo_delta64))
phylogeo_delta166$variant <- rep("Delta", nrow(phylogeo_delta166))
phylogeo_delta179$variant <- rep("Delta", nrow(phylogeo_delta179))
phylogeo_delta186$variant <- rep("Delta", nrow(phylogeo_delta186))
phylogeo_delta191$variant <- rep("Delta", nrow(phylogeo_delta191))
phylogeo_delta254$variant <- rep("Delta", nrow(phylogeo_delta254))
phylogeo_gamma35$variant <- rep("Gamma", nrow(phylogeo_gamma35))
phylogeo_gamma55$variant <- rep("Gamma", nrow(phylogeo_gamma55))
phylogeo_lambda90$variant <- rep("Lambda", nrow(phylogeo_lambda90))
phylogeo_lambda95$variant <- rep("Lambda", nrow(phylogeo_lambda95))
phylogeo_lambda103$variant <- rep("Lambda", nrow(phylogeo_lambda103))
phylogeo_lambda104$variant <- rep("Lambda", nrow(phylogeo_lambda104))
phylogeo_lambda107$variant <- rep("Lambda", nrow(phylogeo_lambda107))
phylogeo_mu15$variant <- rep("Mu", nrow(phylogeo_mu15))
phylogeo_mu26$variant <- rep("Mu", nrow(phylogeo_mu26))
phylogeo_mu61$variant <- rep("Mu", nrow(phylogeo_mu61))
phylogeo_mu73$variant <- rep("Mu", nrow(phylogeo_mu73))
phylogeo_mu82$variant <- rep("Mu", nrow(phylogeo_mu82))

phylogeo_alpha13$invasion_tree <- rep("Alpha_TL_13", nrow(phylogeo_alpha13))
phylogeo_delta10$invasion_tree <- rep("Delta_TL_10", nrow(phylogeo_delta10))
phylogeo_delta64$invasion_tree <- rep("Delta_TL_64", nrow(phylogeo_delta64))
phylogeo_delta166$invasion_tree <- rep("Delta_TL_66", nrow(phylogeo_delta166))
phylogeo_delta179$invasion_tree <- rep("Delta_TL_179", nrow(phylogeo_delta179))
phylogeo_delta186$invasion_tree <- rep("Delta_TL_186", nrow(phylogeo_delta186))
phylogeo_delta191$invasion_tree <- rep("Delta_TL_191", nrow(phylogeo_delta191))
phylogeo_delta254$invasion_tree <- rep("Delta_TL_254", nrow(phylogeo_delta254))
phylogeo_gamma35$invasion_tree <- rep("Gamma_TL_35", nrow(phylogeo_gamma35))
phylogeo_gamma55$invasion_tree <- rep("Gamma_TL_55", nrow(phylogeo_gamma55))
phylogeo_lambda90$invasion_tree <- rep("Lambda_TL_90", nrow(phylogeo_lambda90))
phylogeo_lambda95$invasion_tree <- rep("Lambda_TL_95", nrow(phylogeo_lambda95))
phylogeo_lambda103$invasion_tree <- rep("Lambda_TL_103", nrow(phylogeo_lambda103))
phylogeo_lambda104$invasion_tree <- rep("Lambda_TL_104", nrow(phylogeo_lambda104))
phylogeo_lambda107$invasion_tree <- rep("Lambda_TL_107", nrow(phylogeo_lambda107))
phylogeo_mu15$invasion_tree <- rep("Mu_TL_15", nrow(phylogeo_mu15))
phylogeo_mu26$invasion_tree <- rep("Mu_TL_26", nrow(phylogeo_mu26))
phylogeo_mu61$invasion_tree <- rep("Mu_TL_61", nrow(phylogeo_mu61))
phylogeo_mu73$invasion_tree <- rep("Mu_TL_73", nrow(phylogeo_mu73))
phylogeo_mu82$invasion_tree <- rep("Mu_TL_82", nrow(phylogeo_mu82))


########################## Combine data frames ################################

phylogeo_CL <-rbind(phylogeo_alpha13, phylogeo_delta10, phylogeo_delta64,
                    phylogeo_delta166, phylogeo_delta179, phylogeo_delta186,
                    phylogeo_delta191, phylogeo_delta254, phylogeo_gamma35,
                    phylogeo_gamma55, phylogeo_lambda90, phylogeo_lambda95,
                    phylogeo_lambda103, phylogeo_lambda104, phylogeo_lambda107,
                    phylogeo_mu15, phylogeo_mu26, phylogeo_mu61, phylogeo_mu73,
                    phylogeo_mu82)

############################### Write file #####################################

write.table(phylogeo_CL,
            file = "Phylogenetics/Transmission_lineage_processing/CL_transmission_lineages_phylogeo.tsv",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
