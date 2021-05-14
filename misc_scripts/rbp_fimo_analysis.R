# Description: Script to run FIMO, scanning for RBP motifs

#---Load Libraries and data----------------------------------------------------------------------------
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(doParallel)
library(foreach)
library(tidyverse)

# Source fixed file paths
source(here::here("R", "file_paths.R"))

#---Functions------------------------------------------------------------------------------------------

source(here::here("R", "rbp_functions.R"))
source(here::here("R", "get_regions.R"))

#---Main-----------------------------------------------------------------------------------------------

print(str_c(Sys.time(), " - starting script..."))

# Object first created in cluster_validation_btwn_datasets.Rmd
# Granges object with intron counts across all samples
clu <- 
  readRDS(
  file.path(path_to_results,
            "leafcutter/intron_clustering/tissue_polyA_all_clusters_gtfannotated.Rds")
  )

# Import ensembl version 97 .gtf (GRCh38) to get exon start/end. Keep only chr 1-22, X and Y.
ref <- import("/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")
ref <- ref %>% keepSeqlevels(c(1:22, "X", "Y"), pruning.mode = "coarse") 

# Extract adjoining exons/regions for all introns
clu$index <- c(1:length(clu))

clu_intron_exon <- get_adjoining_exons(junc_metadata = clu[, !str_detect(colnames(mcols(clu)), "_SJ_leafcutter")],
                                       index_colname = "index", 
                                       ref_gtf = ref, 
                                       novel_exon_size = 50)

# Add cluster ids
clu_intron_exon <- 
  clu_intron_exon %>% 
  dplyr::inner_join(mcols(clu)[, c("cluster_id", "index")] %>% as_tibble()) %>% 
  dplyr::select(intron_exon_index = index, cluster_id, everything())

clu_proximal <- get_proximal_intronic_regions(intron_exon_coordinates = clu_intron_exon,
                                              exonic_region = 50,
                                              intronic_region = 500)

# Load leafcutter list and create query lists
leafcutter_list <- 
  readRDS(
    file.path(
      path_to_results,
      "leafcutter/diff_splicing_PCaxes/allcomparisons_leafcutter_ds_PCaxes.Rds"
    )
  )

query_list <- creating_RBP_query_annotations(leafcutter_list = leafcutter_list,
                                             intron_exon_coordinates = clu_intron_exon, 
                                             proximal_intronic_regions = clu_proximal,
                                             additional_bp = 50)

# Save query list
saveRDS(
  query_list, 
  file.path(
    path_to_raw_data,
    "rbp_analysis/rbp_query_list.Rds"
  )
)

# Save clu_proximal
saveRDS(
  clu_proximal, 
  file.path(
    path_to_raw_data,
    "rbp_analysis/rbp_intron_exon_prox_coord.Rds"        
  )
)

print(str_c(Sys.time(), " - creating fasta files"))

# Filter to include only comparisons to control
query_list <- query_list[str_detect(names(query_list), "Control")]

# Export fasta files from query list
query_list %>% 
  lapply(., function(list){
    
    for(i in 1:length(list)){
      
      file_name = str_c(unique(list[[i]]$comparison), ":", names(list)[i])
      
      print(str_c(Sys.time(), " - creating: ", file_name, ".fasta"))
      
      generate_fasta_from_gr(gr = list[[i]],
                             genome = BSgenome.Hsapiens.UCSC.hg38.masked,
                             style = "UCSC",
                             remove_star_strand = T,
                             output_path = file.path(path_to_raw_data, "rbp_analysis/rbp_density_query_fasta"),
                             file_name = file_name,
                             reduce = T)
      
    }
    
  })

print(str_c(Sys.time(), " - running fimo analysis"))

# Run MEME analysis
run_fimo_analysis(path_to_fimo ="/home/ssethi/meme/bin/fimo",
                  output_dir = file.path(path_to_results, "rbp_analysis/rbp_density/"),
                  path_to_query_fasta_dir = file.path(path_to_raw_data, "rbp_analysis/rbp_density_query_fasta"),
                  # fasta_filter = "prox", # Comment in to run only proximal intronic regions 
                  path_to_pwm_dir = file.path(path_to_raw_data, "rbp_analysis/attract_db"), 
                  cores = 9)


print(str_c(Sys.time(), " - done!"))
