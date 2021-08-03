library(dplyr)
library(ggpubr)
library(plyr)
library(ggplot2)
library(scales)
library(tidyr)
library(stringr)
library(tibble)
library(patchwork)
library(circlize)
library("ggplot2")
library("ggdendro")
library("reshape2")
library(grid)
library(plot.matrix)
library(pheatmap)
library(RColorBrewer)
library(xlsx)
library("ggpubr")


# 1. Load Data 
means <- read.table('~/significant_means.txt', header = T, sep = '\t')
pvalue <- read.table('~/pvalues.txt', header = T, sep = '\t')
complex <- read.table('~/deconvoluted.txt', header = T, sep = '\t')

# 2. Filter Data

## Filtering for interactions between cells that come from the same TME_pattern patient group

exclude <- means %>% dplyr::select(contains("_A") & contains("_B")) %>% colnames() # for means
microenv <- means %>% dplyr::select(!contains(exclude)) 
microenv <-  microenv %>% mutate(interacting_pair = str_replace(interacting_pair, " ", "")) 


pmicroenv <- pvalue %>% dplyr::select(!contains(pexclude)) # for pvalues
pexclude <- pvalue %>% dplyr::select(contains("_A") & contains("_B")) %>% colnames()
pmicroenv <-  pmicroenv %>% mutate(interacting_pair = str_replace(interacting_pair, " ", "")) 
 


## Rearrange data and add TME_pattern, A = N3MC, B = CP2E

exB <- microenv %>% dplyr::select(contains("_B.")) %>% colnames() 
groupA <- microenv %>% dplyr::select(!contains(exB)) 
groupA <- groupA %>%  tidyr::gather("key", "value",  CD14_Macrophages1_A.CD14_Macrophages1_A:Tumor_A.Tumor_A) %>% 
  dplyr::filter(value != "NA") 
groupA <- groupA %>% mutate(group = "A")


exA <- microenv %>% dplyr::select(contains("_A.")) %>% colnames()
groupB <- microenv %>% dplyr::select(!contains(exA)) 
groupB <- groupB %>%  tidyr::gather("key", "value",  CD14_Macrophages1_B.CD14_Macrophages1_B:Tumor_B.Tumor_B) %>% 
  dplyr::filter(value != "NA") 
groupB <- groupB %>% mutate(group = "B")

total <- rbind(groupA, groupB)
total <- total %>% mutate(key = str_replace_all(key, "_A", "")) 
total <- total %>% mutate(key = str_replace_all(key, "_B", "")) 


## Repeat the above for the pvalue object

pexB <- pmicroenv %>% dplyr::select(contains("_B.")) %>% colnames()
pgroupA <- pmicroenv %>% dplyr::select(!contains(pexB)) 
pgroupA <- pgroupA %>%  tidyr::gather("key", "value",  CD14_Macrophages1_A.CD14_Macrophages1_A:Tumor_A.Tumor_A) %>% 
  dplyr::filter(value != "NA") 
pgroupA <- pgroupA %>% mutate(group = "A")



pexA <- pmicroenv %>% dplyr::select(contains("_A.")) %>% colnames()
pgroupB <- pmicroenv %>% dplyr::select(!contains(pexA)) 
pgroupB <- pgroupB %>%  tidyr::gather("key", "value",  CD14_Macrophages1_B.CD14_Macrophages1_B:Tumor_B.Tumor_B) %>% 
  dplyr::filter(value != "NA") 
pgroupB <- pgroupB %>% mutate(group = "B")

ptotal <- rbind(pgroupA, pgroupB)
ptotal <- ptotal %>% mutate(key = str_replace_all(key, "_A", "")) 
ptotal <- ptotal %>% mutate(key = str_replace_all(key, "_B", "")) %>% dplyr::rename(pvalue = value) %>% select(id_cp_interaction, pvalue, key, group)  %>% dplyr::filter(pvalue <= 0.05) 

## Add pvalue to means object

total <- left_join(total, ptotal,by = c("id_cp_interaction"="id_cp_interaction", "key"="key", "group"="group"))

# 3. Addition of annotations

## "cell_comb_id"  describes each existing pair of cell types (e.g. celltype1_celltype2 and celltype2_celltype1 get the same ID)


total <- total %>% mutate(celltypes = key) %>% separate(key, c("celltype_a", "celltype_b"), sep = "\\.") 


combinations <- total %>% dplyr::select(celltype_a, celltype_b) %>% unique()  %>%
  dplyr::filter(celltype_a != celltype_b) %>% arrange(celltype_a)
combinations <- combinations[c(1:10, 12:20, 23:30, 34:40, 45:50, 56:60, 67:70, 78:80, 89:90, 100),] # filter out indidivdual celltype pairs


total <- total %>% mutate(cell_comb_id = "same")

# Make data frame that contains combination ID and cell types

for (i in 1:55){
  combinations$cell_comb_id[i] <- i
}
df <- data.frame(celltype_a = rep(NA, 11),celltype_b = rep(NA, 11), cell_comb_id = 56:66)
un <- unique(total$celltype_a)
for (i in 1:11){
  df$celltype_a[i] <- un[i]
  df$celltype_b[i] <- un[i]
  df$cell_comb_id[i] <- un[i]
}
combinations <- rbind(combinations,df)

# Now add cell_comb_id to total 

for (i in 1:55){
  for (j in 1:length(total$celltype_a)){
    if(total$celltype_a[j] == combinations[i,1] & total$celltype_b[j] == combinations[i,2]){
      total$cell_comb_id[j]<- i }
  }}

for (i in 1:55){
  for (j in 1:length(total$celltype_a)){
    if(total$celltype_b[j] == combinations[i,1] & total$celltype_a[j] == combinations[i,2] & total$cell_comb_id[j] == "same"){
      
      total$cell_comb_id[j]<- i }
  }}

# now add for autocrine interactions

single_cell <- unique(total$celltype_a)

for (i in 1:length(total$celltype_a)){
  for (j in 1:length(single_cell)){
    if(total$cell_comb_id[i] == "same" & total$celltype_a[i] == single_cell[j]) {
      total$cell_comb_id[i] <- single_cell[j]
    }
  }}


## Annotation of TME_patterns that the interacting celltypes belong to



N3MC <- c("NK_cells", "Myofibroblast1", "T_conv1", "CD14_Macrophages1", "Myeloid_Dendritic", "Tumor")
CP2E <- c("CD14_Macrophages2", "Myofibroblast2", "Plasmacytoid_Dendritic", "T_CD8_1", "T_CD8_2", "Tumor")
total <- total %>% 
  mutate( TME = ifelse(celltype_a %in% N3MC & celltype_b %in% N3MC, "N3MC", ifelse(celltype_a %in% CP2E & celltype_b %in% CP2E, "CP2E", "mixed")))


# 4. Count number of Interactions

# Step 1: split by group

groupA <- total %>% dplyr::filter(group == "A")
groupB <- total %>% dplyr::filter(group == "B")

combinations <- total %>% select(cell_comb_id, TME, celltype_a, celltype_b) %>% unique()

groupA_counts <- groupA %>% dplyr::count(cell_comb_id) %>% mutate(group = "A")
groupB_counts <- groupB %>% dplyr::count(cell_comb_id) %>% mutate(group = "B")
group_counts <- left_join(groupA_counts, groupB_counts, by = "cell_comb_id")

group_counts <- group_counts %>% dplyr::rename(n.A = n.x, n.B = n.y, group.A = group.x, group.B = group.y)

group_counts <- left_join(group_counts, combinations, by = "cell_comb_id") # add TME information
tumor_cepp <- group_counts %>% dplyr::filter(cell_comb_id == "Tumor") %>% mutate(TME = "CP2E")
group_counts <- rbind(group_counts, tumor_cepp)


range <- c(group_counts$n.A, group_counts$n.B)

CP2E_no_interactions <- group_counts %>% select(-n.A, -group.A) %>%  
  dplyr::filter(TME == "CP2E") %>% 
  dplyr::rename(n = n.B) %>%  
  ggplot(aes(x = celltype_b, y = celltype_a, fill = n)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlGnBu", limits = c(0, max(range))) +
  ggtitle("Number of Interactions") +
  #geom_text(aes(label = n))+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

ggsave("CP2E_interactions.pdf", plot = CP2E_no_interactions, path = "~/")


N3MC_no_interactions <- group_counts %>% select(-n.B, -group.B) %>%  
  dplyr::filter(TME == "N3MC") %>% 
  dplyr::rename(n = n.A) %>%  
  ggplot(aes(x = celltype_b, y = celltype_a, fill = n)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlGnBu", limits = c(0, max(range))) +
  ggtitle("Number of Interactions") +
  #geom_text(aes(label = n))+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

ggsave("N3MC_no_interactions", plot = N3MC_no_interactions, path = "~/")

# 5. Correlation analysis between number of interactions and sequencing depths per celltype

# Comparing number of interactions to number of cells
metadata <-  read.table('~/cellphonedb_meta.txt', header = T, sep = "\t")
countdata <- read.table('~/cellphonedb_count.txt', header = T, sep = "\t")
raw_counts <- read.csv("~/count_raw.csv")


counts <- colSums(countdata)

N3MC_a <- c("Tumor_A", "Myofibroblast1_A", "NK_cells_A", "CD14_Macrophages1_A", "T_conv1_A", "Myeloid_Dendritic_A")
CP2E_b <- c("Tumor_B", "Myofibroblast2_B", "CD14_Macrophages2_B", "T_CD8_1_B", "T_CD8_2_B", "Plasmacytoid_Dendritic_B")
combinations_all <- total  %>% dplyr::filter((group == "A"  & celltype_a %in% N3MC &celltype_b %in% N3MC) | (group == "B"  & celltype_a %in% CP2E & celltype_b %in% CP2E)) %>% 
  dplyr::select(celltype_a, celltype_b, celltypes, TME, cell_comb_id) %>% unique() 

t <- combinations_all %>% dplyr::filter(cell_comb_id == "Tumor") %>% 
  mutate(TME = "CP2E") 
TME_info <- rbind(combinations_all, t)# make sure to have both Tumor-N3MC and Tumor_CP2E

# Create table with celltype counts


cell_type_counts <-  metadata %>% dplyr::count(cell_type) %>% 
  mutate(group = ifelse(cell_type %in% N3MC_a, "N3MC", ifelse(cell_type %in% CP2E_b, "CP2E", NA))) %>% 
  dplyr::filter(!is.na(group)) %>% 
  mutate(celltype_a = str_replace_all(cell_type, "_A", "")) %>% 
  mutate(celltype_a = str_replace_all(celltype_a, "_B", "")) %>% 
  mutate(TME = group) %>% 
  select(-cell_type, - group)


cell_type_counts <- left_join(cell_type_counts, TME_info, by = c("celltype_a", "TME")) %>%  
  mutate(combined_counts = NA)

for(i in 1:length(cell_type_counts$cell_comb_id)){
  a <- as.numeric(cell_type_counts$n[i])
  b <- cell_type_counts %>% filter(celltype_b[i] == celltype_a) %>% dplyr::filter(TME == cell_type_counts$TME[i]) %>% select(n)
  b <- as.numeric(b[1,1]) 
  cell_type_counts$combined_counts[i] <- sum(a,b)
}


cell_type_counts <- cell_type_counts %>% dplyr::arrange(cell_comb_id) 
# need to get rid of every other row as it is just a duplicate with reversed celltypes
cell_type_counts <- cell_type_counts[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61:72),]

# Create table of raw counts


counts_raw <- colSums(raw_counts[,-1])

# combining counts and celltype assignment
raw_counts_celltype <- cbind(metadata,counts_raw)

raw_counts_celltype <- raw_counts_celltype %>% 
  group_by(cell_type) %>% 
  dplyr::summarise(cummulative_counts = sum(counts_raw)) %>% 
  mutate(group = ifelse(cell_type %in% coi_a, "N3MC", ifelse(cell_type %in% coi_b, "CP2E", NA))) %>% dplyr::filter(!is.na(group))  %>% 
  mutate(celltype_a = str_replace_all(cell_type, "_A", "")) %>% 
  mutate(celltype_a = str_replace_all(celltype_a, "_B", "")) %>% 
  select(- cell_type)

raw_counts_celltype <- raw_counts_celltype %>% mutate(TME = group) %>% select(- group)
raw_counts_celltype <- left_join(raw_counts_celltype, TME_info, by = c("celltype_a", "TME")) %>% 
  mutate(combined_counts = NA)

for(i in 1:length(raw_counts_celltype$cell_comb_id)){
  a <- as.numeric(raw_counts_celltype$cummulative_counts[i])
  b <-raw_counts_celltype%>% filter(celltype_b[i] == celltype_a) %>% dplyr::filter(TME == raw_counts_celltype$TME[i]) %>% 
    select(cummulative_counts)
  b <- as.numeric(b[1,1]) 
  raw_counts_celltype$combined_counts[i] <- sum(a,b)
}



raw_counts_celltype <-raw_counts_celltype %>% dplyr::arrange(cell_comb_id) 
raw_counts_celltype <- raw_counts_celltype[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61:72),]

# Creating matrix of interaction counts


CP2E_interaction_counts <- group_counts %>% select(-n.A, -group.A) %>%  
  dplyr::filter(TME == "CP2E") %>% 
  dplyr::rename(n = n.B) 
# need to get rid of every other row as it is just a duplicate with reversed celltypes
CP2E_interaction_counts <- CP2E_interaction_counts[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,32,33,34,35,36),]


N3MC_interaction_counts <- group_counts %>% select(-n.B, -group.B) %>%  
  dplyr::filter(TME == "N3MC") %>% 
  dplyr::rename(n = n.A)
# need to get rid of every other row as it is just a duplicate with reversed celltypes
N3MC_interaction_counts<- N3MC_interaction_counts[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,32,33,34,35,36),]

interaction_counts <- full_join(CP2E_interaction_counts, N3MC_interaction_counts)
interaction_counts[is.na(interaction_counts)] <- 0


# Comparison of number of interactions vs. number of cells



cell_type_counts <- cell_type_counts %>% select(celltypes,TME, combined_counts, cell_comb_id)
interaction_counts<- interaction_counts %>% unite(col = celltypes, celltype_a, celltype_b, sep = ".", remove = F ) %>% 
  dplyr::rename(interactions = n) %>%  select(interactions, cell_comb_id, TME)

corr_cell_number <- full_join(cell_type_counts, interaction_counts, by = c("cell_comb_id", "TME")) %>% unique() %>% 
  select(celltypes, interactions, combined_counts)



corr_cell_number_plot <- ggscatter(corr_cell_number, x = "interactions", y = "combined_counts", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Number of interactions", ylab = "Number of cells per celltype pair")


ggsave("spearman_correlation_cellnumber_vs_interactions.pdf", plot = corr_cell_number_plot, path = "~/")

# Comparison of number of interactions vs. number of mean counts per celltype pair

raw_counts_celltype <- raw_counts_celltype %>% select(celltypes,TME, combined_counts, cell_comb_id)

corr_mean_counts <- full_join(raw_counts_celltype, interaction_counts, by = c("cell_comb_id", "TME")) %>% unique() %>% 
  select(celltypes, interactions, combined_counts)

# Get mean counts by dividing total counts by number of cells per interacting celltypes

corr_mean_counts <- corr_mean_counts %>% mutate(combined_counts = corr_mean_counts$combined_counts/corr_cell_number$combined_counts)

corr_mean_counts_plot <- ggscatter(corr_mean_counts, x = "interactions", y = "combined_counts", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Number of interactions", ylab = "Number of mean counts per celltype pair")

ggsave("spearman_correlation_cellnumber_vs_interactions.pdf", plot = corr_plot1, path = "~/")








# 6. Compare interactions of specific receptor-families and the corresponding ligand-families 

ligands <- read.csv("~/LIGANDS.csv", header = T, sep = ";")
receptors <- read.csv("~/RECEPTORS.csv", header = T, fileEncoding = "UCS-2LE", sep = "\t")
combinations_all <- total  %>% dplyr::filter((group == "A"  & celltype_a %in% N3MC &celltype_b %in% N3MC) | (group == "B"  & celltype_a %in% CP2E & celltype_b %in% CP2E)) %>% 
  dplyr::select(celltype_a, celltype_b, celltypes, TME, cell_comb_id) %>% unique() 





## Add complex_name to "total" to also cover those receptors and ligands that are complexes
complex <- complex %>% dplyr::filter(is_complex == "True") %>% select(id_cp_interaction, complex_name, gene_name) %>% unique()
total <- left_join(total, complex,by = "id_cp_interaction" )

total <- total %>% 
  mutate(gene_a = ifelse(!is.na(complex_name) & gene_a == "", complex_name, gene_a)) %>% 
  mutate(gene_b = ifelse(!is.na(complex_name) & gene_b == "", complex_name, gene_b))

## Add all complex names to ligands and receptors

ligand_spread <- ligands %>% 
  gather("key", "gene_name", X:X.104) %>% 
  dplyr::filter(gene_name != "") %>% 
  select(-key) 

ligand_complex <- complex %>% 
  dplyr::filter(gene_name %in% ligand_spread$gene_name)

additional_ligands <- full_join(ligand_complex, ligand_spread, by = "gene_name") %>% 
  select(gene_name, LIGAND.FAMILY, complex_name) %>% 
  unique()


ligands_final <- full_join(ligand_spread, additional_ligands) 

## Add all complex names to receptors

receptor_spread <- receptors %>% 
  gather("key", "gene_name", X:X.65) %>% 
  dplyr::filter(gene_name != "") %>% 
  select(-key) 

receptor_complex <- complex %>% 
  dplyr::filter(gene_name %in% receptor_spread$gene_name)

additional_receptors <- full_join(receptor_complex, receptor_spread, by = "gene_name") %>% 
  select(gene_name, RECEPTOR.FAMILY, complex_name) %>% unique()


receptors_final <- full_join(receptor_spread, additional_receptors) 

#split by receptor-/ligand-family
rec <- receptors_final %>% 
  group_by(RECEPTOR.FAMILY) %>% group_split()
lig <- ligands_final %>% 
  group_by(LIGAND.FAMILY) %>% group_split()


# Collect all interactions in one final map
final_map <- tibble(rowname = "test")


# iterate over receptor-families, skip those where there's no interaction present
for (i in c(2:4,8,10,12,1,13,15:16)){
  
  data <- total %>% #filter only those interactions where "Tumor" cells are on the receiving end, i.e. Tumor expresses the receptor
    dplyr::filter((group == "A"  & celltype_a %in% N3MC & celltype_b %in% N3MC) | (group == "B"  & celltype_a %in% CP2E & celltype_b %in% CP2E)) %>% 
    dplyr::filter(celltypes != "Tumor.Tumor") %>%  
    dplyr::filter((celltype_a == "Tumor" & receptor_a == "True") | (celltype_b == "Tumor" & receptor_b == "True"))
  
  # receptors and ligands
  
  recep <- rec[[i]]
  liga <- lig[[i]]
  
  
  # filter for receptor and ligand pairs
  figured <- data %>% 
    dplyr::select(TME, interacting_pair, pvalue, celltypes, cell_comb_id, gene_a, gene_b, receptor_a, receptor_b, complex_name) %>% 
    dplyr::filter((receptor_a == "True" & gene_a %in% recep$gene_name & receptor_b == "False" & gene_b %in% liga$gene_name) | 
                    (receptor_b == "True" & gene_b %in% recep$gene_name & receptor_a == "False" & gene_a %in% liga$gene_name ) | 
                    (complex_name != "" & complex_name != "NA" & complex_name %in% recep$complex_name) | 
                    (complex_name != "" & complex_name != "NA" & complex_name %in% liga$complex_name))  
  
  # take negative log of pvalue 
  figured <- figured %>% 
    dplyr::select(-receptor_a, -receptor_b, -gene_a, -gene_b, -complex_name) %>% 
    unique() %>% 
    mutate(pvalue = ifelse(pvalue == "0", 0.0009, pvalue)) %>%
    mutate(neg_log_pvalue = -log(pvalue)) %>% 
    tidyr::spread(interacting_pair, neg_log_pvalue) %>%
    dplyr::select(-pvalue)
  
  figured[is.na(figured)] <- 0 #set NAs to 0
  figured <- figured %>% select_if(~any(. != 0)) #filter out rows that have only 0s
  
 # summarize interactions based on originating celltype
  fd <- figured  %>% 
    dplyr::arrange(TME)  %>% 
    dplyr::select(-TME, -celltypes) %>% 
    group_by(cell_comb_id) %>% 
    dplyr::summarise_all(funs(sum))
  
  comb_all <- combinations_all %>% 
    dplyr::select(celltype_a, celltype_b, cell_comb_id, TME) %>% 
    dplyr::filter(celltype_a != "Tumor") %>% 
    dplyr::select(-celltype_b)
  
  fd <- inner_join(fd, comb_all, by = "cell_comb_id") %>% 
    arrange(TME) %>% 
    column_to_rownames(var = "celltype_a") %>% 
    select(-cell_comb_id, -TME) %>%  
    t() 
  
  rows <- rownames(fd) 
  k <- fd  %>% as_tibble() %>% 
    mutate(rowname = rows) %>% 
    dplyr::mutate(interaction_family = as.character(rec[[i]][1,1])) 
  final_map <- full_join(final_map, k)
  
}


final_map <- final_map %>% dplyr::filter(rowname != "test") # get rid of first row

annot <- final_map %>% 
  select(rowname, interaction_family) %>% 
  group_by(rowname) %>% # some receptors/ligands are cross listed between families, here they are grouped together
  dplyr::summarise(interaction_families = paste(interaction_family, collapse=", ")) 


TME_id <- data %>% select(TME, interacting_pair, value, celltype_a, celltype_b, cell_comb_id) %>% 
  unique() %>% spread(interacting_pair, value) %>% 
  select(celltype_a, celltype_b, TME ) %>% 
  dplyr::filter(celltype_a != "Tumor") %>%  
  arrange(TME) %>% 
  dplyr::select(-celltype_b) %>% 
  column_to_rownames(var ="celltype_a")

final_map[is.na(final_map)] <- 0
final_map <- left_join(final_map, annot, by = "rowname") %>% 
  select(-interaction_family) %>%
  unique() %>% select(-interaction_families)  %>% 
  column_to_rownames("rowname")

annot <- annot %>% arrange(interaction_families) %>% 
  column_to_rownames("rowname")

final_map <- final_map[,rownames(TME_id)]


final_map <- final_map[rownames(annot),]


heatmap <- pheatmap(final_map, 
               color = brewer.pal(n = 9, name = "OrRd"),
               angle_col = "45", 
               cluster_cols = F, cluster_rows = F, 
               annotation_row = annot, annotation_col = TME_id,
               annotation_names_row = T)


save_pheatmap_pdf <- function(x, filename, width=9, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(heatmap, "heatmap_receptor_families.pdf")


# 7. Print table of all relevant interactions (CellphoneDB results)

cpdb_results <- total %>% dplyr::filter(TME != "mixed") %>%  
  mutate(avg_expression = value) %>% 
  select( - complex_name, - gene_name,  -value, - group, - cell_comb_id)


write.xlsx(cpdb_results, file="CellphoneDB_supplement.xlsx", sheetName="CellphoneDB_interactions", row.names=FALSE)
write.xlsx(receptors, file="CellphoneDB_supplement.xlsx", sheetName="Receptor_families", append=TRUE, row.names=FALSE)
write.xlsx(ligands, file="CellphoneDB_supplement.xlsx", sheetName="Ligand_families", append=TRUE, row.names=FALSE)
