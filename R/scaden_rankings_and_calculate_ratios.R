#' Determine cell-type rankings and calculate ratios.
#'
#' Function to determine rankings of cell-type within a sample using their
#' proportions and to calculate ratios from scaden predictions and snRNA-seq
#' cell-type assignments.
#'
#' @param prediction_dir Directory where scaden prediction files are located.
#'   Prediction files should follow the naming pattern:
#'   scaden_predictions_cells.$CELLS:samples:$SAMPLES, where $CELLS refers to
#'   the number of cells used to simulate bulk profiles, and $SAMPLES to the
#'   number of bulk samples simulated per individual.
#' @param unclassified_cell_type Vector containing name used for cells that are
#'   considered unclassified in the snRNA-sequencing experiment. Default = NULL.
#' @param sample_info Dataframe of sample metadata, with any pertinent columns
#'   to be joined to the final proportions dataframe e.g. disease group.
#' @param snRNAseq_proportions Matrix of snRNA_seq proportions, as calculated
#'   using the calculate_snRNA_proportions() function.
#' @param cell_type_class Dataframe with cell types and their overarching class
#'   e.g. endothelial cells and perictyes might both be considered vascular.
#'   Column containing cell types should be named 'Celltype'.
#' @param glia List of column names considered glia. Should be provided in the
#'   format vars(a, b, c), without any quotes around a, b and c.
#' @param nonneuron List of column names considered non-neuronal. Should be
#'   provided in the format vars(a, b, c), without any quotes around a, b and c.
#' @param neuron List of column names considered neuronal. Should be provided in
#'   the format vars(a, b, c), without any quotes around a, b and c.
#' @param astrocyte Column name for any astrocyte related cell types. Should be
#'   provided in the format vars(a, b, c), without any quotes around a, b and c.
#' @param oligos Column name for any oligodendrocyte related cell types. Should
#'   be provided in the format vars(a, b, c), without any quotes around a, b and
#'   c.
#' @param microglia Column name for any microglia related cell types. Should be
#'   provided in the format vars(a, b, c), without any quotes around a, b and c.
#'
#' @return List of 3 dataframes. (1) Ranks: Proportions from scaden and
#'   snRNAseq, with within-sample rankings of each cell type. (2) Ratios:
#'   glia-to-neuron ratios and non-neuron-to-neuron ratios per sample based on
#'   proportions derived from scaden or snRNA-seq. (3) Glial proportions:
#'   proportion of one glial cell-type compared to the total glial proportion
#'   per sample based on proportions derived from scaden or snRNA-seq.
#' @export
#' 


scaden_rankings_and_calculate_ratios <- function(prediction_dir, 
                                                 unclassified_cell_type = NULL, 
                                                 sample_info, 
                                                 snRNAseq_proportions = NULL, 
                                                 cell_type_class, 
                                                 glia, 
                                                 nonneuron, 
                                                 neuron,
                                                 astrocyte,
                                                 oligos,
                                                 microglia){
  
  library(tidyverse)
  library(stringr)
  
  scaden_df <- data.frame(
    file_path = list.files(path = prediction_dir, full.names = T, pattern = "scaden_predictions"),
    cells = list.files(path = prediction_dir,
                       pattern = "scaden_predictions") %>% 
      str_replace("scaden_predictions_", "") %>% 
      str_replace(":.*", "") %>% 
      str_replace(".*\\.", ""),
    samples = list.files(path = prediction_dir, pattern = "scaden_predictions") %>% 
      str_replace("scaden_predictions_", "") %>% 
      str_replace(".*:", "") %>% 
      str_replace(".txt", "") %>% 
      str_replace(".*\\.", "")) 
  
  for(i in 1:nrow(scaden_df)){
    
    if(is.null(unclassified_cell_type)){
      
      scaden <- read_delim(file = scaden_df$file_path[i] %>% 
                             as.character(),
                           delim = "\t") %>% 
        dplyr::rename(sample_id = X1) 
      
    } else{
      
      scaden <- read_delim(file = scaden_df$file_path[i] %>% 
                             as.character(),
                           delim = "\t") %>% 
        dplyr::rename(sample_id = X1) %>% 
        dplyr::mutate(!!unclassified_cell_type := NA)
      
    }
    
    ranks <- scaden %>% 
      tidyr::gather(key = "Celltype", value = "cell_type_proportion", -sample_id) %>% 
      dplyr::group_by(sample_id) %>% 
      dplyr::arrange(sample_id, desc(cell_type_proportion)) %>% 
      # Add column with proportion in %, ranking and source
      dplyr::mutate(cell_type_proportion = round((cell_type_proportion * 100), 6), 
                    rank_within_sample = rank(-cell_type_proportion, ties.method = "average"),
                    source = "Scaden deconvolution",
                    cells = scaden_df$cells[i],
                    samples = scaden_df$samples[i])
    
    ratios <- scaden  %>% 
      # New columns with ratios
      # !!! is the unquote-splice operate, which will unquote each element of a list of variables, as provided as an argument in the function with vars()
      dplyr::mutate(glia_neuron_ratio = rowSums(dplyr::select(.,!!!glia), na.rm = T)/rowSums(dplyr::select(.,!!!neuron), na.rm = T),
                    nonneuron_neuron_ratio = rowSums(dplyr::select(.,!!!nonneuron), na.rm = T)/rowSums(dplyr::select(.,!!!neuron), na.rm = T)) %>% 
      dplyr::select(sample_id, glia_neuron_ratio, nonneuron_neuron_ratio) %>%
      tidyr::gather(key = ratio, value = value, -sample_id) %>%
      dplyr::mutate(source = "Scaden deconvolution",
                    cells = scaden_df$cells[i],
                    samples = scaden_df$samples[i])
    
    glial_proportions <- scaden %>% 
      dplyr::mutate(astro_glial_proportion = rowSums(dplyr::select(.,!!!astrocyte), na.rm = T)/rowSums(dplyr::select(.,!!!glia), na.rm = T),
                    oligo_glial_proportion = rowSums(dplyr::select(.,!!!oligos), na.rm = T)/rowSums(dplyr::select(.,!!!glia), na.rm = T),
                    micro_glial_proportion = rowSums(dplyr::select(.,!!!microglia), na.rm = T)/rowSums(dplyr::select(.,!!!glia), na.rm = T)) %>% 
      dplyr::select(sample_id, astro_glial_proportion, oligo_glial_proportion, micro_glial_proportion) %>%
      tidyr::gather(key = proportion, value = value, -sample_id) %>%
      dplyr::mutate(source = "Scaden deconvolution",
                    cells = scaden_df$cells[i],
                    samples = scaden_df$samples[i])
    
    if(i == 1){
      
      master_ranks <- ranks
      master_ratios <- ratios
      master_glia <- glial_proportions
      
    } else {
      
      master_ranks <- master_ranks %>% 
        dplyr::bind_rows(ranks)
      master_ratios <- master_ratios %>% 
        dplyr::bind_rows(ratios)
      master_glia <- master_glia %>% 
        dplyr::bind_rows(glial_proportions)
      
    }
    
  }
  
  # Join single cell data and sample_info
  master_ranks <- master_ranks %>% 
    purrr::when(!is.null(snRNAseq_proportions) ~ dplyr::bind_rows(., snRNAseq_proportions %>% 
                                                                   tidyr::gather(key = "Celltype", value = "cell_type_proportion", -sample_id) %>% 
                                                                   dplyr::group_by(sample_id) %>% 
                                                                   dplyr::arrange(sample_id, desc(cell_type_proportion)) %>% 
                                                                   # Add column with proportion in %, ranking and source
                                                                   dplyr::mutate(cell_type_proportion = cell_type_proportion * 100,
                                                                                 rank_within_sample = order(order(cell_type_proportion, decreasing=TRUE)),
                                                                                 source = "snRNA-seq")),
                ~ .) %>%
    dplyr::inner_join(cell_type_class) %>% 
    dplyr::inner_join(sample_info)
  
  master_ratios <- master_ratios %>%
    purrr::when(!is.null(snRNAseq_proportions) ~ dplyr::bind_rows(., snRNAseq_proportions %>% 
                                                                    dplyr::mutate(glia_neuron_ratio = rowSums(dplyr::select(.,!!!glia), na.rm = T)/rowSums(dplyr::select(.,!!!neuron), na.rm = T),
                                                                                  nonneuron_neuron_ratio = rowSums(dplyr::select(.,!!!nonneuron), na.rm = T)/rowSums(dplyr::select(.,!!!neuron), na.rm = T)) %>% 
                                                                    dplyr::select(sample_id, glia_neuron_ratio, nonneuron_neuron_ratio) %>%
                                                                    tidyr::gather(key = ratio, value = value, -sample_id) %>% 
                                                                    dplyr::mutate(source = "snRNA-seq")),
                ~ .) %>%
    dplyr::inner_join(sample_info)
  
  master_glia <- master_glia %>% 
    purrr::when(!is.null(snRNAseq_proportions) ~ dplyr::bind_rows(., snRNAseq_proportions %>% 
                                                                   dplyr::mutate(astro_glial_proportion = rowSums(dplyr::select(.,!!!astrocyte), na.rm = T)/rowSums(dplyr::select(.,!!!glia), na.rm = T),
                                                                                 oligo_glial_proportion = rowSums(dplyr::select(.,!!!oligos), na.rm = T)/rowSums(dplyr::select(.,!!!glia), na.rm = T),
                                                                                 micro_glial_proportion = rowSums(dplyr::select(.,!!!microglia), na.rm = T)/rowSums(dplyr::select(.,!!!glia), na.rm = T)) %>% 
                                                                   dplyr::select(sample_id, astro_glial_proportion, oligo_glial_proportion, micro_glial_proportion) %>%
                                                                   tidyr::gather(key = proportion, value = value, -sample_id) %>%
                                                                   dplyr::mutate(source = "snRNA-seq")),
                ~ .) %>%
    dplyr::inner_join(sample_info)
  
  master_list <- setNames(list(master_ranks, master_ratios, master_glia), c("ranks", "ratios", "glial_proportions"))
  
  return(master_list)
  
}
