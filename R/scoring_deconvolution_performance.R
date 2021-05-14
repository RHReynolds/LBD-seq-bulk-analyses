#' Score performance of deconvolution instances
#'
#' Function to score the performance of various instances of one or several
#' deconvolution algorithms. Function is based on the assumptions that the
#' truths derived from von Bartheld et al.
#' (https://www.ncbi.nlm.nih.gov/pubmed/27187682) represent reality. These
#' ground truths include: (1) the glia-to-neuron ratio in cerebral cortex is
#' ~1.5 (2) within the glial cell population, oligodendrocytes account for
#' 37-76% of glial cells, astrocytes for 17-47% of glial cells and microglia
#' 5-17% and (3) within glial cell types, oligodendrocytes tend to make up the
#' largest proportion, followed by astrocytes and finally microglia.
#'
#' @param results_list List of 3 dataframes, containing: (1) Ranks: Proportions
#'   from deconvolution approach, with within-sample rankings of each cell type.
#'   (2) Ratios: glia-to-neuron ratios and non-neuron-to-neuron ratios per
#'   sample based on proportions derived from deconvolution. (3) Glial
#'   proportions: proportion of one glial cell-type compared to the total glial
#'   proportion per sample based on proportions derived from deconvolution.
#' @param true_GNR "True" glia-to-neuron ratio as determined from literature.
#'   According to von Bartheld et al.
#'   (https://www.ncbi.nlm.nih.gov/pubmed/27187682) in cerebral cortex grey
#'   matter this is ~1.5.
#' @param deconvolution_variables Column names that denote differentiate between
#'   deconvolution instances that user wishes to judge performance of. E.g. for
#'   deconvolution performed by scaden, this could be the number of cells used
#'   to generate one bulk simulated sample and the number of bulk simulated
#'   samples generated. Should be provided in the format vars(a, b, c), without
#'   any quotes around a, b and c.
#' @param glia Vector of cell types considered to be glia i.e. oligodendrocyte-,
#'   astrocyte- and microglia-related cell types.
#' @param oligo Vector with cell-type name used to refer to mature
#'   oligodendrocyte population.
#' @param astrocyte Vector with cell-type name used to refer to mature astrocyte
#'   population.
#' @param microglia Vector with cell-type name used to refer to microglial
#'   population.
#' @param mean_glial_proportions Named vectors (astro, oligo, microglia) with
#'   the mean proportion of each glial cell type compared to the total glial
#'   population. Default is set and calculated using values from 3 papers cited
#'   in Table 3 of https://www.ncbi.nlm.nih.gov/pubmed/27187682.Papers include
#'   Pope 1958, Pope 1959 and Pelvig 2008.
#' @param use_ranges Logical. If TRUE, an additional score
#'   "*_glial_proportion_within_range" will be computed.
#' @param ranges If use_ranges = TRUE, use this argument to provide the expected
#'   glial proportions for astrocytes, oligodendrocytes and microglia.
#'
#' @return Dataframe with 10 scores. (1) GNR_score: ranking of instances from
#'   best (equivalent t0 rank of 1) to worst by the deviation of its median
#'   glia-to-neuron ratio from the true ratio. (2-4)
#'   *_glial_proportion_within_range: ranking of each instance from best
#'   (equivalent to rank of 1) to worst based on the proportion of samples
#'   within a deconvolution instance that have a proportion of the glial
#'   cell-type in question that falls within expected ranges. (5-7)
#'   *_glial_proportion_deviation_from_mean: ranking of each instance from best
#'   (equivalent to rank 1) to worst by the deviation of its median
#'   astro/oligo/microglial glial proportion from the expected mean (8)
#'   glia_ranking_score: ranking of deconvolution instance form best (equivalent
#'   to rank of 1) to worst by the proportion of samples wherein
#'   oligodendrocytes make up the biggest glial cell-type population, followed
#'   by astrocytes and then microglia. (9) overall_score: based on score of
#'   instance in scores 1-5. For each score, award score of 1 to the
#'   deconvolution instances with the highest rank (where 1 = best). (10)
#'   overall_rank: final ranking of deconvolution instances based on
#'   overall_score i.e. highest score receives rank of 1. Note: if values tie
#'   upon ranking, the average of their combined ranks is computed.
#' @export

scoring_deconvolution_performance <- function(results_list, true_GNR, deconvolution_variables, 
                                              glia, oligo, astrocyte, microglia,
                                              mean_glial_proportions = c( "astro" = 0.3675, "oligo" = 0.522, "microglia" = 0.1088),
                                              use_ranges = FALSE, 
                                              ranges = data.frame(range = c("min", "max"),
                                                                  astro = c(0.17, 0.47),
                                                                  oligo = c(0.37, 0.76),
                                                                  microglia = c(0.05, 0.17))){

  # Rank each instance of deconvolution by its deviation from the "true" glia-to-neuron ratio
  GNR_rank <- results_list$ratios %>% 
    dplyr::filter(ratio == "glia_neuron_ratio") %>% 
    # Calculate absolute deviation from "true" GNR per row
    dplyr::mutate(deviation_from_truth = abs(true_GNR - value)) %>%
    dplyr::group_by(!!!deconvolution_variables) %>% 
    # Summarise across deconvolution variables by median. 
    dplyr::summarise(median_deviation_from_truth = median(deviation_from_truth, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    # Rank performance. Lower median = lower deviation from true GNR = higher rank
    dplyr::mutate(GNR_rank = rank(median_deviation_from_truth, ties.method = "average"))
  
  # Score each instance of deconvolution by deviation of the astro/oligo/microglial to total glial proportion from the expected mean
  glia_proportion_mean_score <- results_list$glial_proportions %>% 
    dplyr::mutate(deviation_from_truth = ifelse(proportion == "astro_glial_proportion", abs(mean_glial_proportions["astro"] - value),
                                                ifelse(proportion == "oligo_glial_proportion", abs(mean_glial_proportions["oligo"] - value),
                                                       ifelse(proportion == "micro_glial_proportion", abs(mean_glial_proportions["microglia"] - value), NA)))) %>% 
    dplyr::mutate(proportion = str_replace(proportion, "_proportion", "_prop_deviation_from_mean_rank")) %>% 
    dplyr::group_by(!!!deconvolution_variables, proportion) %>% 
    dplyr::summarise(median_deviation_from_truth = median(deviation_from_truth, na.rm = T)) %>% 
    dplyr::group_by(proportion) %>% 
    dplyr::mutate(proportion_score = rank(median_deviation_from_truth, ties.method = "average")) %>%
    dplyr::select(-median_deviation_from_truth) %>% 
    tidyr::spread(key = proportion, value = proportion_score)

  # Score each instance of deconvolution by the proportion of samples wherein glial cell-type rankings
  # follow the expected ranking of largest --> smallest i.e. oligo --> astro --> microglia
  glia_rank_score <- results_list$ranks %>% 
    # Select main glial cell-types i.e. oligodendrocytes, astrocytes and microglia
    dplyr::filter(Celltype %in% glia) %>% 
    dplyr::select(-rank_within_sample) %>% 
    dplyr::group_by(!!!deconvolution_variables, sample_id) %>% 
    # Add two new columns.
    # 1) Rank glial cell types by their cell-type proportions
    # 2) Check whether glial ranking follows expected ranking of 1) oligodendrocyte, 2) astrocyte, 3) microglia
    dplyr::mutate(glia_rank = rank(-cell_type_proportion, ties.method = "average"),
                  reflects_expected_rank = ifelse(Celltype %in% oligo & glia_rank == 1, TRUE,
                                                  ifelse(Celltype %in% astrocyte & glia_rank == 2, TRUE,
                                                         ifelse(Celltype %in% microglia & glia_rank == 3, TRUE, FALSE)))) %>% 
    # Summarise by checking whether all logical values in reflects_expected_rank for a single sample in a single deconvolution run = TRUE 
    # i.e. all glial cell types have the expected rank
    dplyr::summarise(all_true = all(reflects_expected_rank)) %>% 
    dplyr::group_by(!!!deconvolution_variables) %>% 
    # Across all samples in a single deconvolution run, check proportion where glial cell types follow expected ranking
    dplyr::summarise(proportion_true = mean(all_true)) %>% 
    dplyr::ungroup() %>% 
    # Rank proportions i.e. higher proportion = higher rank.
    dplyr::mutate(glia_ranking_score = rank(-proportion_true, ties.method = "average"))
  
  if(use_ranges == TRUE){
    
    # Score each instance of deconvolution by the proportion of samples wherein the proportions of
    # one glial cell-type compared to the total glial proportion falls within expected range
    glia_proportion_score <- results_list$glial_proportions %>%
      # Glial cell type proportion within range?
      dplyr::mutate(within_range = ifelse(proportion == "astro_glial_proportion" & value >= ranges$astro[1] & value <= ranges$astro[2], TRUE,
                                          ifelse(proportion == "oligo_glial_proportion" & value >= ranges$oligo[1] & value <= ranges$oligo[2], TRUE,
                                                 ifelse(proportion == "micro_glial_proportion" & value >= ranges$microglia[1] & value <= ranges$microglia[2], TRUE, FALSE)))) %>%
      dplyr::mutate(proportion = str_replace(proportion, "_proportion", "_prop_within_range_rank")) %>% 
      dplyr::group_by(!!!deconvolution_variables, proportion) %>%
      # Divide number of samples where proportions fell within range by total number of samples
      dplyr::summarise(proportion_within_range = sum(within_range)/n()) %>%
      dplyr::group_by(proportion) %>%
      # Rank performance
      dplyr::mutate(proportion_score = rank(-proportion_within_range, ties.method = "average")) %>%
      dplyr::select(-proportion_within_range) %>%
      tidyr::spread(key = proportion, value = proportion_score)
    
    # Join all scores
    all_scores <- GNR_rank %>% 
      dplyr::select(-median_deviation_from_truth) %>% 
      dplyr::inner_join(glia_proportion_score) %>% 
      dplyr::inner_join(glia_proportion_mean_score) %>% 
      dplyr::inner_join(glia_rank_score %>% 
                          dplyr::select(-proportion_true))
    
  } else{
    
    # Join all scores
    all_scores <- GNR_rank %>% 
      dplyr::select(-median_deviation_from_truth) %>% 
      dplyr::inner_join(glia_proportion_mean_score) %>% 
      dplyr::inner_join(glia_rank_score %>% 
                          dplyr::select(-proportion_true))
    
  }
  
  # Calculate an overall score based on performance in 5 other scores
  # Thereafter rank deconvolution instances by highest score
  overall_score <- all_scores %>% 
    dplyr::inner_join(all_scores %>% 
                        tidyr::gather(key = "scoring_system", value = "rank", -cells, -samples) %>% 
                        dplyr::group_by(scoring_system) %>% 
                        # For each score, award score of 1 to the deconvolution instances with the highest rank (where 1 = best)
                        dplyr::mutate(score = ifelse(rank == min(rank), 1, 0)) %>% 
                        dplyr::group_by(!!!deconvolution_variables) %>% 
                        # Sum score across all scoring systems
                        dplyr::summarise(overall_score = sum(score)) %>% 
                        dplyr::ungroup() %>% 
                        dplyr::mutate(overall_rank = rank(-overall_score, ties.method = "average")))
  
  return(overall_score)
  
}