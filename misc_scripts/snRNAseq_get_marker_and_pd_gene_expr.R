# Description: Script for to obtain for single-nucleus marker gene expression and SNCA expression

# Load packages -----------------------------------------------------------

library(data.table)
library(ggpubr)
library(readxl)
library(tidyverse)
library(Seurat)
library(stringr)

# Source fixed file paths
source(here::here("R", "file_paths.R"))

source(here::here("docs", "figures", "colour_schemes_and_themes.R"))

# Load sample info ---------------------------------------------------------

sample_info <-
  read_excel(
    path = file.path(path_to_raw_data, 
                     "sample_details/20201229_MasterFile_SampleInfo.xlsx"),
    sheet = "SampleInfo",
    skip = 1
  ) %>% 
  dplyr::select(
    sample_id = CaseNo, 
    disease_group = Disease_Group
  ) %>% 
  dplyr::mutate(
    sample_id =
      case_when(sample_id == "C036" ~ "C36",
                sample_id == "C048" ~ "C48",
                TRUE ~ sample_id),
    disease_group =
      as.factor(disease_group) %>% 
      fct_relevel(.,
                  c("Control", "PD", "PDD", "DLB"))
  ) %>% 
  dplyr::distinct()

# Load seurat objects ------------------------------------------------------

file_paths <- 
  list.files(
    path = file.path(path_to_results,
                     "snRNA/seurat_objects"), 
    pattern = "SO.",
    full.names = T
  )

so_list <-
  setNames(
    object = 
      file_paths %>% 
      lapply(., function(file_path){
        
        readRDS(file_path)
        
      }),
    nm = 
      file_paths %>% 
      basename() %>% 
      str_remove(".rds") %>% 
      str_remove("SO.")
  )

# Add cell type and disease group in metadata
for(i in 1:length(so_list)){
  
  cell_type <- names(so_list)[i]
  
  metadata <-
    so_list[[i]]@meta.data %>% 
    tibble::rownames_to_column() %>% 
    dplyr::left_join(sample_info, by = c("orig.ident" = "sample_id")) %>% 
    dplyr::mutate(cell_type = cell_type) %>% 
    column_to_rownames() %>% 
    dplyr::select(cell_type, disease_group)
  
  so_list[[i]] <- AddMetaData(so_list[[i]], metadata = metadata)
  
}

# Merge seurat objects
so_merged <- 
  merge(x = so_list$Astrocytes, 
        y = c(
          so_list$Excitatory,
          so_list$Inhibitory,
          so_list$Microglia,
          so_list$Oligodendrocytes,
          so_list$OPC,
          so_list$Vascular
        )
  )

rm(so_list)

# Re-factor cell types
so_merged@meta.data$cell_type <-
  so_merged@meta.data$cell_type %>% 
  as.factor() %>% 
  fct_relevel(
    ., 
    c("Excitatory",
      "Inhibitory",
      "Astrocytes",
      "Oligodendrocytes",
      "OPC",
      "Microglia",
      "Vascular")
  )

# Normalised merged seurat -------------------------------------------------

so_merged <- 
  Seurat::NormalizeData(
    so_merged, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )

so_merged <- 
  Seurat::FindVariableFeatures(
    so_merged, 
    selection.method = "mean.var.plot", 
    nfeatures = 2000
  )

all_genes <- rownames(so_merged)

so_merged <- 
  Seurat::ScaleData(
    so_merged, 
    features = all_genes
  )

# Fetch marker gene data for plotting -------------------------------------

# marker_genes <- 
#   tibble(
#     gene = c("SLC17A7", "GAD2", "AQP4", "PLP1", "PDGFRA", "P2RY12", "FLT1"),
#     marker = c("Excitatory", "Inhibitory", "Astrocytes", "Oligodendrocytes", "OPC", "Microglia", "Vascular")
#   )

marker_genes <- 
  tibble(
    gene = c("SLC17A7", "CAMK2A", "NRGN",
             "GAD1", "GAD2", 
             "MBP", "MOBP", "PLP1",
             "PDGFRA", "VCAN",
             "CD74", "CSF1R", "C3",
             "AQP4", "GFAP",
             "FLT1", "CLDN5"),
    marker = c(rep("Excitatory", 3),
               rep("Inhibitory", 2),
               rep("Oligodendrocytes", 3),
               rep("OPC", 2),
               rep("Microglia", 3),
               rep("Astrocytes", 2),
               rep("Vascular", 2))
  )


marker_genes_data <-
  Seurat::FetchData(object = so_merged, vars = marker_genes$gene, slot = 'data') %>% 
  as_tibble(rownames = "cell_id") %>% 
  dplyr::inner_join(so_merged@meta.data %>% 
                      as_tibble(rownames = "cell_id") %>% 
                      dplyr::rename(sample_id = orig.ident)) 

# Fetch PD gene data for plotting -----------------------------------------
pd_genes <- 
  read_delim(
    file.path(path_to_raw_data, 
              "misc/PD_risk_genes/Blauwendraat2019_PDgenes.txt"),
    delim = "\t"
  )

pd_genes_data <-
  Seurat::FetchData(object = so_merged, vars = pd_genes$Gene, slot = 'data') %>% 
  as_tibble(rownames = "cell_id") %>% 
  dplyr::inner_join(so_merged@meta.data %>% 
                      as_tibble(rownames = "cell_id") %>% 
                      dplyr::rename(sample_id = orig.ident)) 

# UMAP for excitatory and inhibitory ----------------------------------------

so_ex_inh <- 
  subset(
    x = so_merged, 
    subset = (cell_type %in% c("Excitatory", "Inhibitory"))
  )

so_ex_inh <-
  Seurat::RunPCA(
    so_ex_inh,
    npcs = 30,
    verbose = FALSE
  )

so_ex_inh <-
  Seurat::FindNeighbors(
    so_ex_inh,
    dims = 1:30
  )

so_ex_inh <-
  Seurat::FindClusters(
    so_ex_inh,
    resolution = 0.5
  )

so_ex_inh <-
  Seurat::RunUMAP(
    so_ex_inh,
    dims = 1:30
  )

# Plot excitatory and inhibitory neurons & SNCA expression across neurons
fig_list <- vector(mode = "list", length = 2)

fig_list[[1]] <- 
  Seurat::DimPlot(
    so_ex_inh, 
    reduction = "umap", 
    group.by = "cell_type"
  )

fig_list[[2]] <- 
  Seurat::FeaturePlot(
    so_ex_inh, 
    features = "SNCA", 
    reduction = "umap"
  )

fig_list_data <-
  fig_list %>% 
  lapply(., function(fig){
    
    fig$data
    
  })

# Save data ---------------------------------------------------------------

out_dir <- file.path(path_to_results, "snRNA/seurat_objects")
saveRDS(
  object = so_merged, 
  file = file.path(out_dir, "so_normalised_all.rds"), 
  compress = FALSE
)
saveRDS(
  object = so_ex_inh, 
  file = file.path(out_dir, "so_normalised_excitatory_inhibitory.rds"),
  compress = FALSE
)
saveRDS(
  object = fig_list_data, 
  file = file.path(out_dir, "plot_data_umap_snca.rds")
)
saveRDS( 
  setNames(
    object = list(marker_genes, marker_genes_data),
    nm = c("marker_genes", "expression_data")
  ), 
  file = file.path(out_dir, "norm_expr_marker_genes.rds")
)
data.table::fwrite( 
  pd_genes_data, 
  file = file.path(out_dir, "norm_expr_pd_genes.txt"),
  sep = "\t"
)

