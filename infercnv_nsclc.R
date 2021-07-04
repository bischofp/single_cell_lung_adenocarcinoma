library(tidyverse)
library(Seurat)
library(infercnv)
library(cowplot)
library(ggdendro)
library(dendextend)
setwd("/extra/flo/sc/scnsclc")

## load data ------------------------------------------

## load preprocessed seurat object  
seu_epi_pri <- read_rds("_data/computed/epi_anno.RDS")
seu_epi_pri2 <- read_rds("_data/computed/epi.RDS")

meta_tbl <- as_tibble(FetchData(seu_epi_pri, c("cell_id", "sample_id", "main_cell_type", "tissue_type", "patient_id")))

table(meta_tbl$main_cell_type)
table(meta_tbl$tissue_type)

## global vars
# cell_type_colors <- setNames(unique(seu_epi_pri$cell_type_epi_color), unique(seu_epi_pri$cell_type_epi_custom))

cols.use <- c(
  p018 = "#ee2f00",
  p019 = "#ef3e24",
  p023 = "#f37120",
  p024 = "#feba11",
  p027 = "#b3d234",
  p028 = "#84c440",
  p029 = "#75c000",
  p030 = "#6dbe45",
  p031 = "#71bf44",
  p032 = "#6fc48c",
  p033 = "#6fcbcd",
  p034 = "#35c5f3",
  CMS1 = "#eba83a", CMS2 = "#027eb5", CMS3 = "#d684ae", CMS4 = "#00a881",
  `CMS1,CMS2` = "#779378", `CMS1,CMS3` = "#E19674", `CMS1,CMS4` = "#76A85E", 
  `CMS2,CMS3` = "#6C81B2", `CMS2,CMS4` = "#01939B",  
  G1 = "#94b6d2",
  S = "#dc4040",
  G2M = "#7aa77f",
  `NA` = "grey",
  Normal = "steelblue",
  `Tumor core` = "red",
  `Tumor border` = "pink",
  CNA = "red",
  CNN = "grey"
)

pids <- sort(unique(seu_epi_pri$patient_id))
pids_with_tumor <- FetchData(seu_epi_pri, c("tissue_type", "patient_id")) %>% as_tibble() %>% distinct(tissue_type, patient_id) %>% filter(tissue_type == "Tumor") %>% pull(patient_id) %>% sort

## prepare infercnv inputs ---------------------------------

## 1. count matrix: cells x genes
rcm <- seu_epi_pri@assays$RNA@counts
#cell_subset <- seu_epi_pri$cell_id[seu_epi_pri$tissue_type != "Normal"]
#rcm <- rcm[,cell_subset]
#colnames(rcm) <- str_replace_all(colnames(rcm), ":", "_")

## 2. two column annotation table: cell_id, cell_type 
af <- FetchData(seu_epi_pri, c("main_cell_type", "tissue_type", "patient_id", "cell_id")) %>%
  # filter(cell_id %in% cell_subset) %>%
  mutate(cell_type_infercnv = ifelse(tissue_type != "Normal", paste0("malignant_", patient_id), paste0("normal_", patient_id)),
         cell_id = str_replace_all(cell_id, ":", "_")) %>% 
  select(cell_id, cell_type_infercnv) %>% 
  as_tibble
write.table(af, file = "_data/infercnv_nsclc/annotation.txt", sep = "\t", col.names = F, row.names = F)

## 3. four column gene location table: gene, chr, start, end
## used the one from infercnv website
gof <- read_tsv("_data/infercnv_nsclc/gencode_v21_gen_pos.complete.txt", col_names = c("gene", "chr", "start", "end"))
common_genes <- intersect(gof$gene, rownames(rcm))

## 4. defining cell types for baseline reference
baseline_cell_types <- sort(unique(af$cell_type_infercnv)[!str_detect(unique(af$cell_type_infercnv), "malignant")])

## run infer cnv --------------------------------------------
 
## run time: ~12h
 
### create infercnv object
cnv_obj <- infercnv::CreateInfercnvObject(
 raw_counts_matrix = rcm[common_genes,],
 annotations_file = "_data/infercnv_nsclc/annotation.txt",
 gene_order_file = "_data/infercnv_nsclc/gencode_v21_gen_pos.complete.txt",
 ref_group_names = baseline_cell_types
) 

#cnv_obj <- read_rds("_data/infercnv_nsclc/out/preliminary.infercnv_obj")

## run infercnv
#cnv_obj <- infercnv::run(
#  cnv_obj, cutoff=0.1, #out_dir="_data/infercnv_nsclc/out", 
#  cluster_by_groups=T, denoise=F, HMM=F, 
#  num_threads=parallel::detectCores()
#)
#write_rds(cnv_obj, #"_data/infercnv_nsclc/infercnv_results.rds")

cnv_obj <- read_rds("_data/infercnv_nsclc/infercnv_results.rds")


## custom plotting --------------------------------------------

gof <- filter(gof, chr %in% paste0("chr", 1:22),
              gene %in% common_genes) %>% 
  mutate(chr = ordered(paste0(" ", chr), levels = paste0(" chr", 1:22)))

## infercnv expression results
expr_all <- cnv_obj@expr.data
expr_tumor <- cnv_obj@expr.data[,seu_epi_pri$tissue_type != "Normal"]

expr_all_cut <- expr_all
expr_all_cut[expr_all > 1.15] <- 1.15
expr_all_cut[expr_all < 0.85] <- 0.85

gene_chr_list <- cnv_obj@gene_order %>% 
  mutate(x = rownames(.)) %>% 
  group_by(chr) %>% 
  select(-start, -stop) %>% 
  nest(x) %>% 
  {setNames(.$data, .$chr)} %>% 
  lapply(function(x) x$x)


## infercnv dendrogram
dendro_lines <- read_lines("_data/infercnv_nsclc/out/infercnv.observations_dendrogram.txt")
dendro_list <- lapply(dendro_lines, function(x) phylogram::read.dendrogram(text = x)) %>% setNames(pids_with_tumor)
cell_order <- lapply(dendro_lines, function(x) str_split(x, "\\(|:|,|\\)") %>% unlist %>% .[str_detect(., "p0")]) %>% unlist

               
## find clones
clone_tbl <- lapply(dendro_list, function(x) enframe(cutree(x, k = 4), "cell_id", "clone")) %>% 
  bind_rows(.id = "patient_id") %>% 
  mutate(clone = as.character(clone)) %>% 
  bind_rows(
    FetchData(seu_epi_pri, c("patient_id", "cell_id")) %>%
      mutate(cell_id = str_replace_all(cell_id, ":", "_"),
             clone = "Normal") %>%
      filter(seu_epi_pri$tissue_type == "Normal")
  )

clone_scores <- clone_tbl %>% 
  group_by(patient_id, clone) %>%
  # do(cna_score = tryCatch(max(unlist(lapply(gene_chr_list, function(x) mean(colSums(abs(log2(expr_all[x,.$cell_id])))/nrow(expr_all))))), error = function(e) NA)) %>%
  do(cna_score = tryCatch(quantile(unlist(lapply(gene_chr_list, function(x) mean(abs(colSums(log2(expr_all[x,.$cell_id])))))), 0.75), error = function(e) NA)) %>%
  ungroup() %>%
  unnest(cna_score) %>%
  spread(clone, cna_score) %>%
  mutate(mean_n = mean(Normal, na.rm = T)) %>%
  gather(clone, cna_score, -mean_n, -patient_id) %>%
  mutate(cna_score = cna_score/mean_n) %>%
  mutate(max_n = max(cna_score[.$clone == "Normal"], na.rm = T)) %>%
  mutate(cna_clone = ifelse(cna_score > max_n, "CNA", "CNN")) %>%
  mutate(cna_clone = ifelse(cna_score > 1.15, "CNA", "CNN")) %>%
  mutate(cna_clone = ifelse(patient_id == "p018" & cna_score < 1.45, "CNN", cna_clone)) %>%
  mutate(cna_clone = ifelse(patient_id == "p030" & cna_score > 1, "CNA", cna_clone)) %>%
  mutate(cna_clone = ifelse(patient_id == "p033" & cna_score > 0.9, "CNA", cna_clone)) %>%
  mutate(cna_clone = ifelse(clone == "Normal", "CNN", cna_clone)) %>%
  select(patient_id, cna_score, cna_clone, clone) %>%
  # filter(clone != "Normal") %>% 
  left_join(clone_tbl, by = c("clone", "patient_id")) %>%
  left_join(enframe(seu_epi_pri$tissue_type, "cell_id", "tissue_type")) %>%
  mutate(cell_id = ordered(cell_id, levels = cell_order)) %>%
  mutate(helper_var_c = " CN",
         helper_var_t = " Tissue")

ggplot(clone_scores, aes(clone, cna_score, color = patient_id, shape = cna_clone)) + 
  #geom_boxplot() + 
  geom_point() +
  geom_hline(yintercept = 1.15) +
  scale_color_manual(values=cols.use) + 
  facet_wrap(~patient_id) + 
  coord_flip()
                                            
write_tsv(select(clone_scores, -helper_var_c), "_data/_tab/infercnv_clone_scores_nsclc.tsv")
                    
segments_tbl <- lapply(dendro_list, function(x) segment(dendro_data(x)) %>% mutate(helper_var_d = "Dendrogram")) %>% 
  bind_rows(.id = "patient_id") 

## cell annotation
anno_data <- as_tibble(FetchData(seu_epi_pri, c("patient_id", "tissue_type", "cell_id"))) %>%
  mutate(cell_id = str_replace_all(cell_id, ":", "_"),
         cell_id = ordered(cell_id, levels = cell_order)) %>%
  filter(str_detect(tissue_type, "Tumor")) %>%
  mutate(helper_var_p = " Patient",
         helper_var_c = " Cell type")


## sub plot annotation -----------------------------------

#plot_anno_cell <- ggplot(filter(clone_scores, clone != "Normal")) +
#  geom_tile(aes(helper_var_c, cell_id, fill = tissue_type)) +
#  theme_void() +
#  facet_grid(patient_id~helper_var_t, scales = "free", space = "free") +
#  scale_fill_manual(values = cols.use) +
#  guides(fill = F) +
#  theme(panel.spacing.x = unit(0, "npc"),
#        panel.spacing.y = unit(0.003, "npc"),
#        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
#        strip.text.y = element_blank()) 

plot_anno_patient <- ggplot(anno_data) +
  geom_tile(aes(helper_var_p, cell_id, fill = patient_id)) +
  theme_void() +
  facet_grid(patient_id~helper_var_p, scales = "free", space = "free") +
  scale_fill_manual(values = cols.use) +
  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank())

plot_anno_clone <- ggplot(filter(clone_scores, clone != "Normal", cell_id %in% anno_data$cell_id)) +
  geom_tile(aes(helper_var_c, cell_id, fill = cna_clone)) +
  theme_void() +
  facet_grid(patient_id~helper_var_c, scales = "free", space = "free") +
  scale_fill_manual(values = cols.use) +  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank())
                       
plot_dendro <- ggplot(segments_tbl) +
  geom_segment(aes(x=-y,y=x,xend=-yend,yend=xend),size=0.5) +
  theme_void() +
  facet_grid(patient_id~helper_var_d, scales = "free", space = "free") +
  theme(panel.spacing = unit(0, "npc"),
        panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text = element_blank()) +
  scale_y_continuous(expand=c(0,0))

pg_anno <- plot_grid(plot_dendro, plot_anno_patient, 
                     plot_anno_clone, 
                     #plot_anno_cell, 
                     nrow = 1, align = "h", rel_widths = c(0.7, 0.15, 0.15))

ggsave("_fig/infercnv_heat_anno_nsclc.pdf", pg_anno, width = 1.4, height = 20)

## full plot  -----------------------------------

## joining data 
plot_data <- as.data.frame(expr_tumor) %>% 
  as_tibble(rownames = "gene") %>%
  gather(cell_id, value, -gene) %>%
  left_join(gof, by = "gene") %>% 
  mutate(cell_id = ordered(cell_id, levels = cell_order)) %>%
  #filter(chr == 13) %>% 
  left_join(anno_data, by = "cell_id") %>% 
  mutate(cell_id = ordered(cell_id, levels = cell_order)) %>% 
  mutate(value_cut = ifelse(value > 1.15, 1.15, ifelse(value < 0.85, 0.85, value))) %>% 
  group_by(chr) %>% 
  mutate(rank = rank(start)) %>%
  ungroup


plot_data_sub <- filter(plot_data, chr == " chr21", patient_id == "p019")

plot_cnv_sub <- ggplot(plot_data_sub) +
  ggrastr::geom_tile_rast(aes(rank, cell_id, fill=value_cut)) +
  theme_void() +
  facet_grid(patient_id~chr, scales = "free", space = "free") +
  scale_fill_gradient2(midpoint = 1, low = scales::muted("blue"), high = scales::muted("red")) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank()) +
  labs(fill = "inferCNV\nexpression")

plot_cnv <- ggplot(plot_data) +
  ggrastr::geom_tile_rast(aes(rank, cell_id, fill=value_cut), raster.dpi = 50) +
  theme_void() +
  facet_grid(patient_id~chr, scales = "free", space = "free") +
  scale_fill_gradient2(midpoint = 1, low = scales::muted("blue"), high = scales::muted("red")) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank()) +
  labs(fill = "inferCNV\nexpression")
                       
pg <- plot_grid(plot_dendro, plot_anno_patient, 
                plot_anno_clone, 
                #plot_anno_cell, 
                plot_cnv, nrow = 1, align = "h", 
                rel_widths = c(0.1, 0.02, 0.02, 0.86))

ggsave("_fig/infercnv_heat_nsclc.pdf", pg, width = 15, height = 20)
ggsave("_fig/infercnv_heat_nsclc.png", pg, width = 15, height = 20)


