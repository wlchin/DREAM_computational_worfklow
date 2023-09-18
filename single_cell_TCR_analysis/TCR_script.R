# DREAM scTCRseq
library(tidyverse)
library(readxl)
library(data.table)
library(vegan)
library(igraph)
library(ggraph)
library(cowplot)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(scales)
library(pheatmap)
library(gridExtra)
library(gghighlight)
library(ggrepel)

#1- Get data
#####
## make a list of directories for each csv to input into map function
#unfiltered
path <- "/Users/nicolaprincipe/OneDrive - The University of Western Australia/02_Code/1239-2/DREAM/DREAM_scTCRseq2022_all"
IDs <- dir(path)
path_end <- "outs/all_contig_annotations.csv"
files <- paste(path, IDs, path_end, sep = "/")

data <- map(files, read_csv, id = 'data_path') #sample names in order to ensure list is the same order of sample names

##ID each input file as a column in the data frame
meta_data <- read_csv("/Users/nicolaprincipe/OneDrive - The University of Western Australia/02_Code/1239-2/DREAM/2022_TCR_Metdata_CORRECT.csv")[1:100,1:6] %>%
  mutate(`External ID` = str_replace(`External ID`, 'C2D2', 'C2D1')) %>%
  separate(`GEX ID`, into = c('tmp1', 'tmp2'), sep = 2) %>%
  unite('GEX ID', tmp1, tmp2, sep = '_') %>%
  mutate(`GEX ID` = str_replace(`GEX ID`, '__', '_'))

phenotypes <- read_csv("/Users/nicolaprincipe/OneDrive - The University of Western Australia/02_Code/1239-2/DREAM/metadata_with_umap.csv") %>%
  select('barcode' = ...1, 'GEX ID' = 'LibraryName', predicted.celltype.l1, predicted.celltype.l2, predicted.celltype.l3, 
         Timepoint, PFS_6M, X_axis_umap, Y_axis_umap) %>%
  separate(barcode, into = c('barcode', 'tmp'), sep = '_') %>% 
  separate(`GEX ID`, into = c('tmp1', 'tmp2'), sep = 2) %>% 
  unite('GEX ID', tmp1, tmp2, sep = '_') %>% 
  mutate(`GEX ID` = str_replace(`GEX ID`, '__', '_'))

table(phenotypes$`GEX ID`, phenotypes$predicted.celltype.l2, useNA = 'ifany')

names(data) <- IDs 

unique(phenotypes$predicted.celltype.l2)


#####
#2- Make the data usable
#####
data_long <- bind_rows(data, .id = "sample") %>%
  left_join(meta_data, by = c("sample" = "Sample")) %>%
  separate(col = "External ID", into = c("Patient", "Timepoint"), sep = "_" ) %>%
  mutate(Patient = str_remove(Patient, '[pP]'),
         Timepoint = if_else(Timepoint == 'Bl', 'BL', Timepoint)) %>%
  left_join(phenotypes, by = c('GEX ID', 'barcode'))

table(data_long$predicted.celltype.l1)
table(data_long$`GEX ID`, data_long$predicted.celltype.l1, useNA = 'ifany')

PFS <- data_long %>%
  select('GEX ID', Patient, Timepoint.x, PFS_6M) %>%
  distinct() %>%
  filter(!is.na(PFS_6M)) %>%
  write_csv('metadata_PFS.csv')

PFS_short <- PFS %>%
  select(Patient, PFS_6M) %>%
  distinct()

#make the data usable
data_summary <- data_long %>%
  group_by(Patient, Timepoint.x) %>%
  summarise(n_cells = n_distinct(barcode)) #counting the number of cells (barcodes), and making it a column in the df


#####
#3- Defining a clone - bipartite graph
#####
#filter out non CDR3s/TCRs
reject_vector <- str_detect(data_long$cdr3, "_") | 
  str_detect(data_long$cdr3, "\\*") | 
  str_length(data_long$cdr3) > 20 | 
  str_length(data_long$cdr3) < 8

rejects <- data_long %>% #more a check to see what gets deleted
  filter(reject_vector)

data_long <- data_long %>%
  filter(!reject_vector) %>%
  filter(chain %in% c('TRA', 'TRB')) %>%
  filter(!is.na(cdr3)) 

#make dataframe required for the network
data_edge <- data_long %>% #unnest_longer isnt working
  select(barcode, cdr3, chain, `GEX ID`, Patient, Timepoint.x) %>%
  group_by(barcode, `GEX ID`, Patient, Timepoint.x) %>%
  pivot_wider(names_from = chain, values_from = cdr3) %>%
  mutate(check_col = length(unlist(TRA)) != length(unlist(TRB)) & all(c(length(unlist(TRA)) > 0, length(unlist(TRB)) > 0)), # make it so it only changes the columns that need to be changed
         TRA_length = length(unlist(TRA)), # get the number of CDR3s for each chain
         TRB_length = length(unlist(TRB)),
         TRA = if_else(check_col, list(rep(unlist(TRA), TRB_length)), TRA), # replicates the list of CDR3s to match the length of the other list
         TRB = if_else(check_col, list(rep(unlist(TRB), TRA_length)), TRB)) %>%
  unnest_longer(c(TRA,TRB)) %>% #unnest the nested dataframe - multiple TRA cdr3s or TRB cdr3s in the one cell, 
  ##unnest separates into multiple lines, chosen what columns you want to unnest
  #filter(!is.na(TRA), !is.na(TRB)) %>% #only take cells with a TRA TRB pairs, so remove cells that have only an TRA or only a TRB
  ungroup %>% 
  unite("barcode", barcode, `GEX ID`, sep = 'x') %>%
  select(TRA, TRB, barcode) %>% #graph_from_data_frame function needs edges as the first two columns of the dataframe
  group_by(TRA, TRB) %>% #nest cells with the same clones on the same line
  summarise(barcodes = paste0(barcode, collapse = ':'),
            n_barcodes = n_distinct(barcode)) %>%
  pivot_longer(c(TRA, TRB), names_to = 'chain', values_to = 'cdr3') %>%
  distinct() %>% #gets rid of duplicate rows
  filter(!is.na(cdr3))

data_network <- full_join(data_edge, data_edge, by = c('chain', 'cdr3')) %>%
  filter(barcodes.x != barcodes.y) %>% #puts exact matches i.e. from the same cell in the dataframe as well, so need to delete them
  mutate(smallest = pmin(barcodes.x, barcodes.y), #pmin works on whatever you need to find the min & max of i.e. for letters, A is start of the alphabet so smaller than B
         largest = pmax(barcodes.x, barcodes.y), #look at notebook to see explanation of why we did this
         temp = 1) %>%
  group_by(smallest, largest) %>% 
  mutate(temp2 = cumsum(temp)) %>%
  filter(temp2 == 1) %>%
  ungroup() %>%
  select(barcodes.x, barcodes.y, chain, cdr3) 

data_vertices <- data_edge %>% #need to include a vertices data frame to make sure the cells with only a TRA or TRB get included (won't have any edges)
  group_by(barcodes) %>%
  arrange(chain, cdr3) %>%
  summarise(chains = paste0(chain, collapse = ':'),
            cdr3s = paste0(cdr3, collapse = ':'))

g <- graph_from_data_frame(data_network, directed = F, vertices = data_vertices)

g_components <- decompose(g) #each subnetwork is called a component, each component is going to be a clone (TRA and TRB shared)

V(g_components[[1]]) # V means vertices (another way of saying nodes)
length(V(g_components[[1]]))

g_vertices <- map(g_components,V) #list of all the barcodes (vertices) in each component

clone_names <- paste('clone', as.character(1:length(g_vertices)), sep = '_') #just giving each clone a name from 1 to 92005
# length(g_vertices) gives the total number of vertices in that list

names(g_vertices) <- clone_names #names each list item the clone_names data frame

#get out all the barcodes that go with each clone in a data frame 
clones_barcodes <- data.frame(clone.barcode = names(unlist(g_vertices))) %>%
  separate(clone.barcode, into = c('clone', 'barcodes'), sep = '\\.') %>%
  mutate(barcode = str_split(barcodes, ':')) %>%
  unnest(cols = c(barcode)) %>%
  separate(barcode, into = c('barcode', 'GEX ID'), sep = 'x')

g_vertices_lengths_list <- map(g_components, ~ length(V(.))) #map function - list, ~ function
## want the number of nodes (length(V()) for each g_component - the '.' means it will put in each g_component

g_vertices_lengths <- map_dbl(g_components, ~ length(V(.))) #map out as a double vector, not as a data frame, put as a vector in the values section so we can use the table functions

component_sizes <- as.data.frame.table(table(g_vertices_lengths)) # this table is tallying the number of vertices in each graph component

ggplot(component_sizes, aes(g_vertices_lengths, Freq)) +
  geom_bar(stat = 'identity', fill = "steelblue") +
  theme_classic() +
  labs(x = 'Number of TCRa, TCRB or TCRa:TCRB in a clone', y = 'Number of clones') +
  theme(axis.text=element_text(colour = "black")) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Number of TCRa, TCRB or TCRaB pairs in each clone")


#####
#4- Defining expanded and singlets
#####
#Expanded vs singlets for each predicted cell type l2
expanded <- data_clones %>%
  group_by(predicted.celltype.l2, `GEX ID`, clone, Timepoint.x, PFS_6M) %>%
  summarize(count = n()) %>%
  filter(count > 1) %>%
  group_by(predicted.celltype.l2, Timepoint.x, PFS_6M) %>%
  summarize(num_repeated_clones = n())

singlets <- data_clones %>%
  group_by(predicted.celltype.l2, `GEX ID`, clone, Timepoint.x, PFS_6M) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  group_by(predicted.celltype.l2, Timepoint.x, PFS_6M) %>%
  summarize(num_unique_clones = n())

celltype <- full_join(expanded,singlets) %>%
  ungroup()

celltype$num_repeated_clones <- coalesce(celltype$num_repeated_clones, 0)
celltype$num_unique_clones <- coalesce(celltype$num_unique_clones, 0)

Tcells <- celltype %>%
  filter(predicted.celltype.l2 %in% c('CD4 TCM', 'Treg', 'CD8 TEM', 'CD4 Naive',
                                      'CD8 Naive', 'CD4 Proliferating', 'dnT', 'CD8 Proliferating',
                                      'CD4 TEM', 'CD4 CTL', 'gdT', 'CD8 TCM', 'MAIT', 'NK Proliferating',
                                      'NK_CD56bright', 'NK')) %>%
  group_by(Timepoint.x, PFS_6M) %>%
  mutate(prop_expanded = num_repeated_clones/sum(num_repeated_clones),
         prop_singlets = num_unique_clones/sum(num_unique_clones))

write.csv(Tcells, "Tcells_expandedsinglets.csv")

#Expanded vs singlets for each patient
pt_expanded <- data_clones %>%
  group_by(predicted.celltype.l2, `GEX ID`, Patient, clone, Timepoint.x, PFS_6M) %>%
  summarize(count = n()) %>%
  filter(count > 1) %>%
  group_by(`GEX ID`, Patient, predicted.celltype.l2, Timepoint.x, PFS_6M) %>%
  summarize(num_repeated_clones = n())

pt_singlets <- data_clones %>%
  group_by(predicted.celltype.l2, `GEX ID`, Patient, clone, Timepoint.x, PFS_6M) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  group_by(`GEX ID`, Patient, predicted.celltype.l2, Timepoint.x, PFS_6M) %>%
  summarize(num_unique_clones = n())

pt_celltype <- full_join(pt_expanded,pt_singlets) %>%
  ungroup()

pt_celltype$num_repeated_clones <- coalesce(pt_celltype$num_repeated_clones, 0)
pt_celltype$num_unique_clones <- coalesce(pt_celltype$num_unique_clones, 0)

#Labelling clones as 'expanded' or 'singlets' in data_clones frame with cell barcodes
clonetypes <- data_clones %>%
  group_by(predicted.celltype.l2, `GEX ID`, Patient, clone, Timepoint.x, PFS_6M) %>%
  summarize(count = n())
head(clonetypes)

clonetypes <- clonetypes %>%
  mutate(Clone_Type = ifelse(count > 1, 'Expanded', 'Singlet'))

data_clones_labelled <- data_clones %>%
  left_join(clonetypes[, c("predicted.celltype.l2", "GEX ID", "Patient", "clone", "Timepoint.x", "PFS_6M", "count", "Clone_Type")], 
            by = c("predicted.celltype.l2", "GEX ID", "Patient", "clone", "Timepoint.x", "PFS_6M"))

write.csv(data_clones_labelled, 'data_clones_expanded.csv')


#######
#5 - Logistic regression on expanded and singlets with two sided t-tests from COVID19 nature paper
######
#Just T cells
pt_Tcells <- pt_celltype %>%
  filter(predicted.celltype.l2 %in% c('CD4 TCM', 'Treg', 'CD8 TEM', 'CD4 Naive',
                                      'CD8 Naive', 'CD4 Proliferating', 'dnT', 'CD8 Proliferating',
                                      'CD4 TEM', 'CD4 CTL', 'gdT', 'CD8 TCM', 'MAIT', 'NK Proliferating',
                                      'NK_CD56bright', 'NK')) %>%
  group_by(`GEX ID`) %>%
  mutate(prop_expanded = num_repeated_clones/sum(num_repeated_clones),
         prop_singlets = num_unique_clones/sum(num_unique_clones))

write.csv(pt_Tcells, "Tcells_expandedsingletsperpatient.csv")

pt_Tcells2 <- pt_celltype %>%
  filter(predicted.celltype.l2 %in% c('CD4 TCM', 'Treg', 'CD8 TEM', 'CD4 Naive',
                                      'CD8 Naive', 'CD4 Proliferating', 'dnT', 'CD8 Proliferating',
                                      'CD4 TEM', 'CD4 CTL', 'gdT', 'CD8 TCM', 'MAIT', 'NK Proliferating',
                                      'NK_CD56bright', 'NK')) %>%
  group_by(`GEX ID`) %>%
  mutate(prop_expanded = num_repeated_clones/(sum(num_repeated_clones) + sum(num_unique_clones)),
         prop_singlets = num_unique_clones/sum(num_unique_clones)) %>%
  pivot_longer(cols = c(prop_expanded, prop_singlets), names_to = 'logic_choice', values_to = 'proportion') %>%
  mutate(logic_choice = if_else(logic_choice == 'prop_expanded', 1, 0)) %>% 
  #filter(logic_choice %in% c('prop_expanded')) %>%
  filter(!is.na(predicted.celltype.l2)) # removed NA cell type, when there are fewer p values to adjust, there is less adjustment

pt_Tcell_list <- pt_Tcells2 %>%
  group_by(predicted.celltype.l2, Timepoint.x) %>% #add response, timepoint
  group_split()

comparison_names <- map(pt_Tcell_list, ~ .x[c('predicted.celltype.l2', 'Timepoint.x')]) %>% # get all the grouping values together
  map(distinct) %>% 
  bind_rows()

fitted_models <- map(pt_Tcell_list, ~ glm(proportion ~ logic_choice+PFS_6M, # fit the logistic regression
                                          data = .x, 
                                          family = binomial("logit"))) %>% 
  map(summary) %>% 
  map('coefficients') %>% 
  map(~ .x[2,4]) %>% 
  unlist()

padjusted <- p.adjust(fitted_models, method = 'BH') # adjust p value with Benjamani Hochberg

comparison_names$padjusted <- padjusted # add p values to the comparisons data frame

write.csv(comparison_names, "logisiticregression_Tcellpertimepoint_clonetype_PFS_6M_10Feb.csv")



#####
#6 - Defining volatile vs nonvolatile clones
#####
data_clones <- full_join(data_long, clones_barcodes) %>%
  dplyr::select(barcode, clone, predicted.celltype.l1, clone, predicted.celltype.l2, Timepoint.x, Patient, PFS_6M) %>%
  distinct() %>% #gets rid of duplicate rows so when joining the two data frames together, data_long has a row for each TCR sequenced, rather than row per barcode, want row per barcode to count them later - summarise function
  group_by(Patient) %>%
  fill(PFS_6M, .direction = 'updown') %>% #fill in NAs with PFS_6M score from the patient
  group_by(Patient, Timepoint.x, PFS_6M, clone) %>%
  summarise(count = n_distinct(barcode)) %>%
  group_by(Patient) %>%
  mutate(Timepoint.y = factor(Timepoint.x, levels = c('BL', 'C2D1', 'C3D1')) %>% as.numeric()) %>% #transformed timepoint character string into numeric to use min max normalisation
  mutate(Timepoint.y = (Timepoint.y - min(Timepoint.y))/(max(Timepoint.y) - min(Timepoint.y)), #normalise timepoints between 0 - 1 = max-min normalisation
         n_timepoints = n_distinct(Timepoint.y)) %>% #some patients do not have all 3 timepoints
  filter(n_timepoints == 3) %>% #want to take the patients that have the full 3 timepoints - 31 patients
  select(-Timepoint.x) %>%
  arrange(desc(count))

table(data_clones$PFS_6M, useNA = 'ifany')

n_duplicates <- rep(1:nrow(data_clones), data_clones$count) #ecdf (cumulative distribution) needs each clone to be repeated for how many times it appeared,
# at the moment the dataframe has a count for how many times that clone has appeared. So we are duplicating each clone by how many counts it has.
# The 'dense' method in frank will bring all these clones together to work out the ranking per sample

## below is the correct one for data_clones_ranked with correct frank and ecdf
data_clones_ranked <- data_clones %>%
  ungroup() %>%
  slice(n_duplicates) %>% #adds in the duplicate rows for each clone
  group_by(Patient, Timepoint.y) %>% # important for ranking, make a rank of each clone, for each patient per timepoint, also need this later to add in phenotype dynamics
  mutate(clone_rank_first = frank(count, ties.method = 'dense'), #top clone, largest number to be rank 1
         clone_rank = ecdf(clone_rank_first)(clone_rank_first)) %>% #ecdf input needs the rank (given by frank), and a name for each rank. So instead of giving it a name, just use the ranking as the name
  dplyr::select(-count, -clone_rank_first) %>%
  distinct() %>%
  pivot_wider(names_from = Timepoint.y, values_from = clone_rank, values_fill = 0) %>% #need to make sure every clone has a value at each timepoint - so if a clone doesn't appear at a timepoint, in a sample, it gets a 0
  pivot_longer(c(`0`, `0.5`, `1`), names_to = 'Timepoint.y', values_to = 'clone_rank') %>% #transform dataframe back to og with the addition values
  mutate(Timepoint.y = as.numeric(Timepoint.y))

# explain clone_rank_first and clone_rank
# above two lines make rank 1-22 like the Carlisle, JITC, 2022 paper, but then transforms it into a rank from 0 to 1. 
# this is different to a min-max normalisation, normalises the ranks to take into account the distribution of the clones e.g. one sample will have heaps of the top 50 clones, and then the rest of the clones have a count of 1.
# if we did min-max normalisation, these will be treated as equals
# clone_rank_first - frank, ties.method = 'dense' ranks the top clone as the highest number (not 1) e.g. if there are 22 clones, the top clone i.e. most count, will be rank 22 - Carlisle paper did this the other way round (top clone = rank 1)
# clone_rank - does the R = F(Y) equation using ecdf() to rank them 0 to 1

data_clones_integral_approx <- data_clones_ranked %>%
  group_by(Patient, clone, PFS_6M) %>%
  summarise(integral_approx = integrate(approxfun(Timepoint.y, clone_rank), lower = 0, upper = 1)$value) #integral is area under curve to get a value for clones changing over time

data_clones_ranked_volatility <- data_clones_ranked %>%
  left_join(data_clones_integral_approx) %>%
  mutate(volatility = (clone_rank - integral_approx)**2) %>% #the higher the volatility, the more the clone is changing rank - high rank at some timepoints, low rank at other timepoints
  group_by(Patient, clone, PFS_6M, integral_approx) %>%
  summarise(volatility_integral_approx = integrate(approxfun(Timepoint.y, volatility), lower = 0, upper = 1)$value) #volatility_integral_approx - value of t

ggplot(data_clones_ranked_volatility, aes(integral_approx, volatility_integral_approx)) +
  geom_point() +
  #scale_fill_manual(values = c("pink", "skyblue")) +
  geom_segment(aes(x = 0.5, xend = 1, y = 0.025, yend = 0.025)) +
  geom_segment(aes(x = 0.5, xend = 0.5, y = 0, yend = 0.025)) +
  theme_classic() +
  theme(axis.text=element_text(colour = "black")) +
  labs(x = "Integral", y = "Rank - volatility_integral")

data_clones <- full_join(data_long, clones_barcodes) %>%
  select(barcode, `GEX ID`, clone, predicted.celltype.l1, predicted.celltype.l2) %>%
  distinct() %>%
  merge(PFS, by = 'GEX ID')

nonvolatileclones <- data_clones_ranked_volatility %>%
  filter(volatility_integral_approx < 0.025) %>%
  filter(integral_approx > 0.5)

nonvolatile <- data.frame(nonvolatileclones$Patient, nonvolatileclones$clone)
colnames(nonvolatile)[1] <- "Patient"
colnames(nonvolatile)[2] <- "clone"

data_nonvoltatileclones <- data_clones %>% 
  semi_join(nonvolatile, by = c("Patient", "clone"))


######
#7 - Tracking TCR clones on CD8TEM seurat clusters 
######
clones_CD8TEMclusters <- read.csv("clones_with_clusters.csv") #adding Melvin's seurat clustering of CD8 TEMs
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 0] <- "CD8_TEM_1"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 1] <- "CD8_TEM_2"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 2] <- "CD8_TEM_3"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 4] <- "CD8_TEM_4"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 5] <- "CD8_TEM_5"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 7] <- "CD8_TEM_6"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 8] <- "CD8_TEM_7"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 3] <- "NA"
clones_CD8TEMclusters$seurat_cluster[clones_CD8TEMclusters$seurat_cluster == 6] <- "NA"
clones_CD8TEMclusters <- clones_CD8TEMclusters %>%
  filter(seurat_cluster %in% c("CD8_TEM_1", "CD8_TEM_2", "CD8_TEM_3", "CD8_TEM_4", "CD8_TEM_5", "CD8_TEM_6", "CD8_TEM_7"))

CD8TEM_clusters2 <- left_join(data_clones_labelled, clones_CD8TEMclusters, by = "barcode")
CD8TEM_clusters <- CD8TEM_clusters2 %>%
  filter(predicted.celltype.l2 %in% c("CD8 TEM")) %>%
  select(`GEX ID`, barcode, clone, predicted.celltype.l2, Patient, Timepoint.x, PFS_6M.x, Clone_Type, seurat_cluster) %>%
  filter(seurat_cluster %in% c("CD8_TEM_1", "CD8_TEM_2", "CD8_TEM_3", "CD8_TEM_4", "CD8_TEM_5", "CD8_TEM_6", "CD8_TEM_7"))

# Merging seurat clusters to nonvolatile clones dataframe
CD8TEM_nonvolatile <- data_nonvoltatileclones %>%
  filter(predicted.celltype.l2 %in% c('CD8 TEM')) %>%
  left_join(clones_CD8TEMclusters, by = "barcode") %>%
  select(`GEX ID`, barcode, clone, predicted.celltype.l2, Patient, Timepoint.x, PFS_6M.x,seurat_cluster)

pt_CD8TEM_nonvolatile <- CD8TEM_nonvolatile %>%
  group_by(`GEX ID`, Patient, Timepoint.x, PFS_6M.x, clone, seurat_cluster) %>%
  summarize(num_nonvolatile_clones = n())

#get clone names not at 3 timepoints per patient
trackCD8TEMclones <- CD8TEM_nonvolatile %>%
  group_by(`GEX ID`, Patient, Timepoint.x, PFS_6M.x, clone, seurat_cluster) %>%
  summarize(num_nonvolatile_clones = n()) %>%
  ungroup() %>%
  filter(num_nonvolatile_clones > 1) %>%
  subset(select = c(Patient, PFS_6M.x, Timepoint.x, clone, seurat_cluster)) %>%
  group_by(Patient, clone) %>%
  tidyr::pivot_wider(names_from = Timepoint.x, values_from = seurat_cluster, values_fn = toString) %>%
  separate_rows(BL, sep=",") %>%
  separate_rows(C2D1, sep=",") %>%
  separate_rows(C3D1, sep=",") %>%
  mutate(across(c("BL", "C2D1", "C3D1"), ~str_trim(.x, side = "both"))) %>%
  separate_rows(BL, sep=",") %>%
  separate_rows(C2D1, sep=",") %>%
  separate_rows(C3D1, sep=",") %>%
  na.omit() %>%
  filter(!rowSums(across(c("BL", "C2D1", "C3D1"), ~(.x == "NA")))) %>%
  group_by(clone) %>%
  mutate(cloneX = if(n()>1) paste(clone, "-", row_number(), sep="") else clone) %>%
  ungroup()

write.csv(trackCD8TEMclones, "trackCD8TEMclones.csv")
write.csv(trackCD8TEMclones, "trackCD8TEMclones_string.csv") #makes CD8TEM1-3shared.csv that includes only persistent clones, 
                                                              #and if they are the same or different phenotype across all timepoints 

clones_CD8TEM_shared <- read.csv("CD8TEM1-3_shared.csv")

# Get barcode from RNA for Melvin
CD8TEM_clusters3 <- CD8TEM_clusters2 %>%
  filter(predicted.celltype.l2 %in% c("CD8 TEM")) %>%
  select(`GEX ID`, barcode_RNA, barcode, clone, predicted.celltype.l2, Patient, Timepoint.x, PFS_6M.x, Clone_Type, seurat_cluster) %>%
  filter(seurat_cluster %in% c("CD8_TEM_1", "CD8_TEM_2", "CD8_TEM_3", "CD8_TEM_4", "CD8_TEM_5", "CD8_TEM_6", "CD8_TEM_7"))

# Merge CD8TEM_clusters with clones_CD8TEM_shared on the 'clone' and 'Patient' columns
merged_data <- merge(CD8TEM_clusters3, clones_CD8TEM_shared[,c("Patient","clone","persistent")], by = c("clone", "Patient"), all.x = TRUE)

# Replace NA's in the 'persistent' column with "unique"
merged_data$persistent[is.na(merged_data$persistent)] <- "unique"

write.csv(merged_data, "CD8TEM_clusters_TCRlabels.csv") #double check in excel that everything is labelled correctly 

# frequency of shared CD8TEM1-3 out of all CD8 T cells
CD8all_clusters <- CD8TEM_clusters2 %>%
  filter(predicted.celltype.l2 %in% c("CD8 Naive","CD8 Proliferating","CD8 TCM","CD8 TEM")) %>%
  select(`GEX ID`, barcode, clone, predicted.celltype.l2, Patient, Timepoint.x, PFS_6M.x, Clone_Type, seurat_cluster)

persistentlabels_revised <- read.csv("CD8TEM_clusters_TCRlabels.csv")
colnames(CD8all_clusters)[1] ="GEX.ID"

CD8all_clusters_persistentlabels <- full_join(CD8all_clusters, persitentlabels_revised, by = c("GEX.ID", "barcode",
                                                                                               "Timepoint.x", "PFS_6M.x")) %>%
  select(GEX.ID, barcode, clone.x, predicted.celltype.l2.x, Patient.x, Timepoint.x, PFS_6M.x, Clone_Type.x, seurat_cluster.x, persistent)

CD8all_shared_freqperpatient <- CD8all_clusters_persistentlabels %>%
  group_by(GEX.ID, Patient.x, Timepoint.x, PFS_6M.x, predicted.celltype.l2.x, seurat_cluster.x, persistent) %>%
  summarise(num_cells = n()) %>%
  ungroup() %>%
  group_by(GEX.ID, Patient.x, Timepoint.x, PFS_6M.x) %>% 
  mutate(total_cells = sum(num_cells)) %>%
  ungroup() %>%
  group_by(GEX.ID, Patient.x, Timepoint.x, PFS_6M.x, predicted.celltype.l2.x, seurat_cluster.x, persistent) %>%
  summarise(freq = num_cells/total_cells)

write.csv(CD8all_shared_freqperpatient, "persistent_freqofallCD8s_perpatient.csv")


#########
#8 - Inverse shannon's on CD8TEM clusters
########
CD8TEM_diversity <- CD8TEM_clusters %>%
  group_by(`GEX ID`, Patient, Timepoint.x, PFS_6M.x, seurat_cluster, clone) %>%
  dplyr::summarise(n_barcodes = n_distinct(barcode)) %>%
  dplyr::summarise(H = vegan::diversity(n_barcodes), # base e
                   Hmax = log(n()), # base e
                   Hnorm = H/Hmax,
                   n_unique_clones = n_distinct(clone),
                   total_clones = sum(n_barcodes)) %>%
  write.csv("basicTCRstats_CD8TEM_RvsNR.csv")

CD8TEM_diversity2 <- CD8TEM_clusters %>%
  group_by(`GEX ID`, Patient, Timepoint.x, seurat_cluster, clone) %>%
  dplyr::summarise(n_barcodes = n_distinct(barcode)) %>%
  dplyr::summarise(H = vegan::diversity(n_barcodes), # base e
                   Hmax = log(n()), # base e
                   Hnorm = H/Hmax,
                   n_unique_clones = n_distinct(clone),
                   total_clones = sum(n_barcodes)) %>%
  write.csv("basicTCRstats_CD8TEM.csv")


######
#9 - UMAP on CD8TEM
######
seurat_object <- readRDS("UMAP/UMAP_CD8TEM.rds")
Nicola_df_TCR <- read.csv("CD8TEM_clusters_TCRlabels_revised.csv")

data <- seurat_object %>%
  filter(predicted.celltype.l1 %in% c('CD8 T')) %>%
  select(X, predicted.celltype.l2, seurat_clusters, UMAP_1, UMAP_2)

data <- data %>%
  rename(barcode_RNA = X)
data$seurat_clusters <- as.character(data$seurat_clusters)

data$seurat_clusters[data$seurat_clusters == 0] <- "CD8_TEM_1"
data$seurat_clusters[data$seurat_clusters == 1] <- "CD8_TEM_2"
data$seurat_clusters[data$seurat_clusters == 2] <- "CD8_TEM_3"
data$seurat_clusters[data$seurat_clusters == 4] <- "CD8_TEM_4"
data$seurat_clusters[data$seurat_clusters == 5] <- "CD8_TEM_5"
data$seurat_clusters[data$seurat_clusters == 7] <- "CD8_TEM_6"
data$seurat_clusters[data$seurat_clusters == 8] <- "CD8_TEM_7"
data$seurat_clusters[data$seurat_clusters == 3] <- "NA"
data$seurat_clusters[data$seurat_clusters == 6] <- "NA"
data_clean <- data %>%
  filter(seurat_clusters %in% c("CD8_TEM_1", "CD8_TEM_2", "CD8_TEM_3", "CD8_TEM_4", "CD8_TEM_5", "CD8_TEM_6", "CD8_TEM_7"))

data_labels <- inner_join(data_clean, Nicola_df_TCR, by = "barcode_RNA")

# UMAP of seurat clusters
data_labels_clusters <- data_labels %>%
  group_by(seurat_clusters) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))

ggplot(data = data_labels, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(colour = seurat_clusters)) +
  gghighlight() +
  geom_text_repel(
    data = data_labels_clusters, 
    aes(label = seurat_clusters),
    size = 5,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_classic() +
  theme(legend.position = "right") +
  labs(colour = expression(paste("CD8", " ", T[EM], " ", "cluster"))) +
  theme(axis.text=element_text(colour = "black"),
        axis.line.x = element_line(size = 0.5),  
        axis.text.x = element_text(size = 12, vjust = -0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.line.y = element_line(size = 0.5),
        plot.title = element_text(size = 16, hjust = 0.5))


# Expanded vs Singlets
ggplot(data = data_labels, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(colour = Clone_Type)) +
  geom_point(data = subset(data_labels, Clone_Type == "Singlet"), aes(colour = Clone_Type)) +
  geom_point(data = subset(data_labels, Clone_Type == "Expanded"), aes(colour = Clone_Type)) +
  scale_color_manual(values = c("Expanded" = "#E76254", "Singlet" = "#ADD8E6")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(colour = "Clone Type") +
  theme(axis.text=element_text(colour = "black"),
        axis.line.x = element_line(size = 0.5),  
        axis.text.x = element_text(size = 12, vjust = -0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.line.y = element_line(size = 0.5),
        plot.title = element_text(size = 16, hjust = 0.5))

# Convert 'persistent' to a factor if it isn't already
data_labels$persistent <- as.factor(data_labels$persistent)

# Reorder levels
data_labels$persistent <- factor(data_labels$persistent, levels = c("persistent_different", "unique", "persistent_same"))

ggplot(data = data_labels, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = subset(data_labels, persistent == "unique"), aes(colour = persistent)) +
  geom_point(data = subset(data_labels, persistent == "persistent_different"), aes(colour = persistent)) +
  geom_point(data = subset(data_labels, persistent == "persistent_same"), aes(colour = persistent)) +
  scale_color_manual(values = c("persistent_same" = "#E76254", "persistent_different" = "#72BCD5", "unique" = "#bcbcbc")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(colour = "Clone Type") +
  theme(axis.text=element_text(colour = "black"),
        axis.line.x = element_line(size = 0.5),  
        axis.text.x = element_text(size = 12, vjust = -0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.line.y = element_line(size = 0.5),
        plot.title = element_text(size = 16, hjust = 0.5))