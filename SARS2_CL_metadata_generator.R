################################################################################
########## SARS-CoV-2 Chile phylogenetic analysis metadata creator #############
################################################################################

########################### Bernardo Gutierrez #################################

# Source scripts and data
source("/Users/user/Documents/SARS2_Chile_local/SARS2_Chile/SARS2_CL_master.R")


################ Metadata files from filtered phylogenies ######################

#### Importing data
# Import TempEst-filtered phylogenetic trees
phylo_alpha_gisaid <-
  read.nexus("Phylogenetics/Alpha/TempEst_SC2_CL_Alpha_GISAID_tempestfiltered_corr_root.tree")
phylo_alpha_isp <-
  read.nexus("Phylogenetics/Alpha/TempEst_SARS2_CL_Alpha_ISP_v2_subtree1_heur.tree")

phylo_gamma_gisaid <-
  read.nexus("Phylogenetics/Gamma/TempEst_SARS2_CL_Gamma_GISAID_2022_08_09_18_v2_noroot.tree")
phylo_gamma_isp <-
  read.nexus("Phylogenetics/Gamma/TempEst_SARS2_CL_Gamma_ISP_2022_08_09_18_v2_heur.tree")

phylo_lambda_gisaid <-
  read.nexus("Phylogenetics/Lambda/TempEst_SARS2_CL_Lambda_GISAID_2022_08_09_18_v3_heur.tree")
phylo_lambda_isp <-
  read.nexus("Phylogenetics/Lambda/TempEst_SARS2_CL_Lambda_ISP_2022_08_09_18_v3_heur.tree")

phylo_mu_gisaid <-
  read.nexus("Phylogenetics/Mu/TempEst_SARS2_CL_Mu_GISAID_2022_08_09_17_v2_heuristic.tree")
phylo_mu_isp <-
  read.nexus("Phylogenetics/Mu/TempEst_SARS2_CL_Mu_ISP_2022_08_09_v2_heuristic.tree")

phylo_delta_gisaid <-
  read.nexus("Phylogenetics/Delta/TempEst_SARS2_CL_Delta_GISAID_v2_subtree1_heur.tree")
phylo_delta_isp <-
  read.nexus("Phylogenetics/Delta/TempEst_SARS2_CL_Delta_ISP_2022_08_09_17_v2_heur.tree")


# Extract tip accession numbers
phylo_alpha_gisaid_id <- phylo_alpha_gisaid$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)
phylo_alpha_isp_id <- phylo_alpha_isp$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)

phylo_gamma_gisaid_id <- phylo_gamma_gisaid$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)
phylo_gamma_isp_id <- phylo_gamma_isp$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)

phylo_lambda_gisaid_id <- phylo_lambda_gisaid$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)
phylo_lambda_isp_id <- phylo_lambda_isp$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)

phylo_mu_gisaid_id <- phylo_mu_gisaid$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)
phylo_mu_isp_id <- phylo_mu_isp$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)

phylo_delta_gisaid_id <- phylo_delta_gisaid$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)
phylo_delta_isp_id <- phylo_delta_isp$tip.label |> str_split_fixed(pattern = fixed("|"), 3) |>
  as.data.frame() |> select(V2) |> rename(`Accession ID` = V2)


#### Subset files to match tips in TempEst-filtered phylogenies
# Metadata files
Meta_Alpha_GISAID_filt <- Meta_Alpha_GISAID |> filter(Meta_Alpha_GISAID$`Accession ID` %in% phylo_alpha_gisaid_id$`Accession ID`)
Meta_Alpha_ISP_filt <- Meta_Alpha_ISP |> filter(Meta_Alpha_ISP$`Accession ID` %in% phylo_alpha_isp_id$`Accession ID`)

Meta_Gamma_GISAID_filt <- Meta_Gamma_GISAID |> filter(Meta_Gamma_GISAID$`Accession ID` %in% phylo_gamma_gisaid_id$`Accession ID`)
Meta_Gamma_ISP_filt <- Meta_Gamma_ISP |> filter(Meta_Gamma_ISP$`Accession ID` %in% phylo_gamma_isp_id$`Accession ID`)

Meta_Lambda_GISAID_filt <- Meta_Lambda_GISAID |> filter(Meta_Lambda_GISAID$`Accession ID` %in% phylo_lambda_gisaid_id$`Accession ID`)
Meta_Lambda_ISP_filt <- Meta_Lambda_ISP |> filter(Meta_Lambda_ISP$`Accession ID` %in% phylo_lambda_isp_id$`Accession ID`)

Meta_Mu_GISAID_filt <- Meta_Mu_GISAID |> filter(Meta_Mu_GISAID$`Accession ID` %in% phylo_mu_gisaid_id$`Accession ID`)
Meta_Mu_ISP_filt <- Meta_Mu_ISP |> filter(Meta_Mu_ISP$`Accession ID` %in% phylo_mu_isp_id$`Accession ID`)

Meta_Delta_GISAID_filt <- Meta_Delta_GISAID |> filter(Meta_Delta_GISAID$`Accession ID` %in% phylo_delta_gisaid_id$`Accession ID`)
Meta_Delta_ISP_filt <- Meta_Delta_ISP |> filter(Meta_Delta_ISP$`Accession ID` %in% phylo_delta_isp_id$`Accession ID`)


#### Data export
# Export tsv files with metadata
write.table(Meta_Alpha_GISAID_filt, file = "Data/SC2_CL_Alpha_GISAID_aln.fasta/SARS2_CL_Alpha_GISAID.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Meta_Alpha_ISP_filt, file = "Data/SC2_CL_Alpha_ISP_aln/SARS2_CL_Alpha_ISP.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(Meta_Gamma_GISAID_filt, file = "Data/SC2_CL_Gamma_GISAID_aln.fasta/SARS2_CL_Gamma_GISAID.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Meta_Gamma_ISP_filt, file = "Data/SC2_CL_Gamma_ISP_aln.fasta/SARS2_CL_Gamma_ISP.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(Meta_Lambda_GISAID_filt, file = "Data/SC2_CL_Lambda_GISAID_aln.fasta/SARS2_CL_Lambda_GISAID.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Meta_Lambda_ISP_filt, file = "Data/SC2_CL_Lambda_ISP_aln.fasta/SARS2_CL_Lambda_ISP.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(Meta_Mu_GISAID_filt, file = "Data/SC2_CL_Mu_GISAID_aln.fasta/SARS2_CL_Mu_GISAID.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Meta_Mu_ISP_filt, file = "Data/SC2_CL_Mu_ISP_aln.fasta/SARS2_CL_Mu_ISP.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(Meta_Delta_GISAID_filt, file = "Data/SC2_CL_Delta_GISAID_aln/SARS2_CL_Delta_GISAID.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Meta_Delta_ISP_filt, file = "Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP.tsv",
            row.names = FALSE, sep = "\t", quote = FALSE)


# Export tsv files with collection dates exclusively
write.table(Meta_Alpha_GISAID_filt[,c(1,6)], file = "Data/SC2_CL_Alpha_GISAID_aln.fasta/SARS2_CL_Alpha_GISAID_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)
write.table(Meta_Alpha_ISP_filt[,c(1,6)], file = "Data/SC2_CL_Alpha_ISP_aln/SARS2_CL_Alpha_ISP_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)

write.table(Meta_Gamma_GISAID_filt[,c(1,6)], file = "Data/SC2_CL_Gamma_GISAID_aln.fasta/SARS2_CL_Gamma_GISAID_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)
write.table(Meta_Gamma_ISP_filt[,c(1,6)], file = "Data/SC2_CL_Gamma_ISP_aln.fasta/SARS2_CL_Gamma_ISP_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)

write.table(Meta_Lambda_GISAID_filt[,c(1,6)], file = "Data/SC2_CL_Lambda_GISAID_aln.fasta/SARS2_CL_Lambda_GISAID_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)
write.table(Meta_Lambda_ISP_filt[,c(1,6)], file = "Data/SC2_CL_Lambda_ISP_aln.fasta/SARS2_CL_Lambda_ISP_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)

write.table(Meta_Mu_GISAID_filt[,c(1,6)], file = "Data/SC2_CL_Mu_GISAID_aln.fasta/SARS2_CL_Mu_GISAID_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)
write.table(Meta_Mu_ISP_filt[,c(1,6)], file = "Data/SC2_CL_Mu_ISP_aln.fasta/SARS2_CL_Mu_ISP_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)

write.table(Meta_Delta_GISAID_filt[,c(1,6)], file = "Data/SC2_CL_Delta_GISAID_aln/SARS2_CL_Delta_GISAID_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)
write.table(Meta_Delta_ISP_filt[,c(1,6)], file = "Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP_dates.tsv",
            row.names = FALSE, col.names = c("id", "date"), sep = "\t", quote = FALSE)


### ISP metadata
# Export ISP tsv files for importations (CL v International)
imp_alpha <- Meta_Alpha_ISP_filt |> select(`Virus name`, surveillance) |>
  mutate(surveillance = case_when(surveillance == 'background_sequences' ~ 'non-Chile', TRUE ~ 'Chile'))
write.table(imp_alpha, file = "Data/SC2_CL_Alpha_ISP_aln/SARS2_CL_Alpha_ISP_importation.tsv",
            row.names = FALSE, col.names = c("id", "location"), sep = "\t", quote = FALSE)

imp_gamma <- Meta_Gamma_ISP_filt |> select(`Virus name`, surveillance) |>
  mutate(surveillance = case_when(surveillance == 'background_sequences' ~ 'non-Chile', TRUE ~ 'Chile'))
write.table(imp_gamma, file = "Data/SC2_CL_Gamma_ISP_aln.fasta/SARS2_CL_Gamma_ISP_importation.tsv",
            row.names = FALSE, col.names = c("id", "location"), sep = "\t", quote = FALSE)

imp_lambda <- Meta_Lambda_ISP_filt |> select(`Virus name`, surveillance) |>
  mutate(surveillance = case_when(surveillance == 'background_sequences' ~ 'non-Chile', TRUE ~ 'Chile'))
write.table(imp_lambda, file = "Data/SC2_CL_Lambda_ISP_aln.fasta/SARS2_CL_Lambda_ISP_importation.tsv",
            row.names = FALSE, col.names = c("id", "location"), sep = "\t", quote = FALSE)

imp_mu <- Meta_Mu_ISP_filt |> select(`Virus name`, surveillance) |>
  mutate(surveillance = case_when(surveillance == 'background_sequences' ~ 'non-Chile', TRUE ~ 'Chile'))
write.table(imp_mu, file = "Data/SC2_CL_Mu_ISP_aln.fasta/SARS2_CL_Mu_ISP_importation.tsv",
            row.names = FALSE, col.names = c("id", "location"), sep = "\t", quote = FALSE)

imp_delta <- Meta_Delta_ISP_filt |> select(`Virus name`, surveillance) |>
  mutate(surveillance = case_when(surveillance == 'background_sequences' ~ 'non-Chile', TRUE ~ 'Chile'))
write.table(imp_delta, file = "Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP_importation.tsv",
            row.names = FALSE, col.names = c("id", "location"), sep = "\t", quote = FALSE)


# Export ISP tsv files for surveillance (Airport v Community v Background)
surv_alpha <- Meta_Alpha_ISP_filt |> select(`Virus name`, surveillance)
write.table(surv_alpha, file = "Data/SC2_CL_Alpha_ISP_aln/SARS2_CL_Alpha_ISP_surveillance.tsv",
            row.names = FALSE, col.names = c("id", "surveillance"), sep = "\t", quote = FALSE)

surv_gamma <- Meta_Gamma_ISP_filt |> select(`Virus name`, surveillance)
write.table(surv_gamma, file = "Data/SC2_CL_Gamma_ISP_aln.fasta/SARS2_CL_Gamma_ISP_surveillance.tsv",
            row.names = FALSE, col.names = c("id", "surveillance"), sep = "\t", quote = FALSE)

surv_lambda <- Meta_Lambda_ISP_filt |> select(`Virus name`, surveillance)
write.table(surv_lambda, file = "Data/SC2_CL_Lambda_ISP_aln.fasta/SARS2_CL_Lambda_ISP_surveillance.tsv",
            row.names = FALSE, col.names = c("id", "surveillance"), sep = "\t", quote = FALSE)

surv_mu <- Meta_Mu_ISP_filt |> select(`Virus name`, surveillance)
write.table(surv_mu, file = "Data/SC2_CL_Mu_ISP_aln.fasta/SARS2_CL_Mu_ISP_surveillance.tsv",
            row.names = FALSE, col.names = c("id", "surveillance"), sep = "\t", quote = FALSE)

surv_delta <- Meta_Delta_ISP_filt |> select(`Virus name`, surveillance)
write.table(surv_delta, file = "Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP_surveillance.tsv",
            row.names = FALSE, col.names = c("id", "surveillance"), sep = "\t", quote = FALSE)


# Export ISP tsv files for domestic spread (Regions v Background)
spread_alpha <- Meta_Alpha_ISP_filt |> select(`Virus name`, State)
spread_alpha$State[spread_alpha$`Virus name` %in%
                     Meta_Alpha_ISP_filt$`Virus name`[
                       Meta_Alpha_ISP_filt$surveillance == 'background_sequences']] <- "Background"
write.table(spread_alpha, file = "Data/SC2_CL_Alpha_ISP_aln/SARS2_CL_Alpha_ISP_regions.tsv",
            row.names = FALSE, col.names = c("id", "region"), sep = "\t", quote = FALSE)


spread_gamma <- Meta_Gamma_ISP_filt |> select(`Virus name`, State)
spread_gamma$State[spread_gamma$`Virus name` %in%
                     Meta_Gamma_ISP_filt$`Virus name`[
                       Meta_Gamma_ISP_filt$surveillance == 'background_sequences']] <- "Background"
write.table(spread_gamma, file = "Data/SC2_CL_Gamma_ISP_aln.fasta/SARS2_CL_Gamma_ISP_regions.tsv",
            row.names = FALSE, col.names = c("id", "region"), sep = "\t", quote = FALSE)

spread_lambda <- Meta_Lambda_ISP_filt |> select(`Virus name`, State)
spread_lambda$State[spread_lambda$`Virus name` %in%
                     Meta_Lambda_ISP_filt$`Virus name`[
                       Meta_Lambda_ISP_filt$surveillance == 'background_sequences']] <- "Background"
write.table(spread_lambda, file = "Data/SC2_CL_Lambda_ISP_aln.fasta/SARS2_CL_Lambda_ISP_regions.tsv",
            row.names = FALSE, col.names = c("id", "region"), sep = "\t", quote = FALSE)

spread_mu <- Meta_Mu_ISP_filt |> select(`Virus name`, State)
spread_mu$State[spread_mu$`Virus name` %in%
                     Meta_Mu_ISP_filt$`Virus name`[
                       Meta_Mu_ISP_filt$surveillance == 'background_sequences']] <- "Background"
write.table(spread_mu, file = "Data/SC2_CL_Mu_ISP_aln.fasta/SARS2_CL_Mu_ISP_regions.tsv",
            row.names = FALSE, col.names = c("id", "region"), sep = "\t", quote = FALSE)

spread_delta <- Meta_Delta_ISP_filt |> select(`Virus name`, State)
spread_delta$State[spread_delta$`Virus name` %in%
                     Meta_Delta_ISP_filt$`Virus name`[
                       Meta_Delta_ISP_filt$surveillance == 'background_sequences']] <- "Background"
write.table(spread_delta, file = "Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP_regions.tsv",
            row.names = FALSE, col.names = c("id", "region"), sep = "\t", quote = FALSE)





######## Filter sequence alignments files from filtered phylogenies ############

# Import alignments
Alpha_GISAID_fasta <- read.fasta("Data/SC2_CL_Alpha_GISAID_aln.fasta/SARS2_CL_Alpha_GISAID_2022_08_09_18.fasta")
Alpha_ISP_fasta <- read.fasta("Data/SC2_CL_Alpha_ISP_aln/SARS2_CL_Alpha_ISP_2022_08_09_18_v2.fasta")

Gamma_ISP_fasta <- read.fasta("Data/SC2_CL_Gamma_ISP_aln.fasta/SARS2_CL_Gamma_ISP_2022_08_09_18_v2.fasta")

Lambda_ISP_fasta <- read.fasta("Data/SC2_CL_Lambda_ISP_aln.fasta/SARS2_CL_Lambda_ISP_2022_08_09_18_v2.fasta")

Mu_ISP_fasta <- read.fasta("Data/SC2_CL_Mu_ISP_aln.fasta/SARS2_CL_Mu_ISP_2022_08_09_18_v2.fasta")

Delta_ISP_fasta <- read.fasta("Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP_2022_08_09_17_v2.fasta")

# Filter sequence alignments
Alpha_GISAID_fasta_filt <- Alpha_GISAID_fasta |> filter(names(Alpha_GISAID_fasta) %in% phylo_alpha_gisaid_id$`Accession ID`)
Alpha_ISP_fasta_filt <- Alpha_ISP_fasta |> filter(Meta_Alpha_ISP$`Accession ID` %in% phylo_alpha_isp_id$`Accession ID`)

Meta_Lambda_GISAID_filt <- Meta_Lambda_GISAID |> filter(Meta_Lambda_GISAID$`Accession ID` %in% phylo_lambda_gisaid_id$`Accession ID`)
Meta_Lambda_ISP_filt <- Meta_Lambda_ISP |> filter(Meta_Lambda_ISP$`Accession ID` %in% phylo_lambda_isp_id$`Accession ID`)

Meta_Mu_GISAID_filt <- Meta_Mu_GISAID |> filter(Meta_Mu_GISAID$`Accession ID` %in% phylo_mu_gisaid_id$`Accession ID`)
Meta_Mu_ISP_filt <- Meta_Mu_ISP |> filter(Meta_Mu_ISP$`Accession ID` %in% phylo_mu_isp_id$`Accession ID`)


######## Wee filtering sequence alignments files from filtered phylogenies for BEAST ############

# Import alignment
Delta_ISP_fasta <- read.fasta("Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP_2022_08_09_17_v2.fasta")

# Import TreeTime
treetime_delta_isp <- read.tree("Phylogenetics/TreeTime_Delta_ISP/TreeTime_SC2_CL_Delta_ISP_v2.nwk")

# Extract tip accession numbers
treetime_delta_isp_id <- treetime_delta_isp$tip.label |> as.vector()  |> str_remove_all(pattern = "'")

# Subset fasta file to match tips in TimeTree
Delta_ISP_fasta_BEASTinput <- Delta_ISP_fasta[names(Delta_ISP_fasta) %in% treetime_delta_isp_id]

# Write FASTA file
write.fasta(Delta_ISP_fasta_BEASTinput, names = names(Delta_ISP_fasta_BEASTinput),
            file.out = "Data/SC2_CL_Delta_ISP_aln.fasta/SARS2_CL_Delta_ISP_BEASTinput.fasta")
