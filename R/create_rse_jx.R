#' Convert output of STAR junction alignment to a RangedSummarisedExperiment
#' class
#'
#' @param sj_df dataframe. Output of \code{load_sj_df}, which loads splice
#'   junctions from STAR alignment into a dataframe.
#' @param sample_info dataframe. Sample metadata for all samples that were
#'   aligned with STAR. Ensure that sample names in \code{sample_info} match the
#'   sample names in the input \code{sj_df}.
#' @param sample_id_col chr. Name of column in \code{sample_info} that provides
#'   sample ids.
#'
#' @return RangedSummarisedExperiment. Output of STAR junction alignment in
#'   RangedSummarisedExperiment class.
#' @export
#' 
create_rse_jx <- function(sj_df, sample_info, sample_id_col){
  
  library(tidyverse)
  library(stringr)
  devtools::load_all("/home/rreynolds/packages/dasper/")
  
  sample_id <- enquo(sample_id_col)
  
  #----Create row_ranges object----
  row_ranges <-
    sj_df %>% 
    # only require distinct junctions
    dplyr::distinct(chr, intron_start, intron_end, strand) %>% 
    # filter only for chromosomes 1:22, X, Y, MT
    dplyr::filter(chr %in% c(str_c(1:22), "X", "Y", "M", "MT")) %>% 
    # change strand format to granges strand format
    # add junction ids for easier pairing with assay data
    dplyr::mutate(strand = ifelse(strand == 0, "*", 
                                  ifelse(strand == 1, "+", "-")),
                  junction_id = str_c(chr, ":", intron_start, "-", intron_end, ":", strand)) %>% 
    GenomicRanges::makeGRangesFromDataFrame(.,
                                            start.field = "intron_start",
                                            end.field = "intron_end",
                                            seqnames.field = "chr",
                                            strand.field = "strand",
                                            keep.extra.columns = TRUE) %>% 
    dasper::annotate_junc_ref(junc_metadata = .,
                                gtf = "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")
  
  #----Create counts table--------
  counts <- 
    sj_df %>% 
    dplyr::select(-intron_motif, -in_annotation, -multi_map_reads_junction, -max_splice_alignment_overhang) %>% 
    dplyr::filter(chr %in% c(str_c(1:22), "X", "Y", "M", "MT")) %>% 
    dplyr::mutate(strand = ifelse(strand == 0, "*", 
                                  ifelse(strand == 1, "+", "-")),
                  junction_id = str_c(chr, ":", intron_start, "-", intron_end, ":", strand)) %>% 
    dplyr::select(Sample, junction_id, unique_reads_junction) %>% 
    tidyr::spread(key = Sample, value = unique_reads_junction) %>% 
    replace(is.na(.), 0)
  
  # All rows in row_ranges also in counts table?
  junc_match <- all(c(row_ranges[, "junction_id"] %>% 
                        as.data.frame() %>% 
                        .[["junction_id"]]) %in% c(counts %>% 
                                                     .[["junction_id"]]))
  
  # Junctions in row_ranges in same order as junctions in counts table?
  junc_order <- all(c(row_ranges[, "junction_id"] %>% 
                        as.data.frame() %>% 
                        .[["junction_id"]]) == c(counts %>% 
                                                   .[["junction_id"]]))
  
  if(junc_match == "TRUE"){
    if(junc_order == "FALSE"){
      
      ids <- match(c(row_ranges[, "junction_id"] %>% 
                       as.data.frame() %>% 
                       .[["junction_id"]]), c(counts %>% 
                                                .[["junction_id"]]))
      
      counts <- counts[ids,]
      
    } else{
      print("Order of junctions in row_ranges matches order of junctions in counts table.")
    }
    
  } else{
    stop("Not all junctions in row_ranges are present in the counts table.")
  }
  
  # Transform counts tibble to counts matrix
  counts_matrix <- counts %>% 
    dplyr::select(-junction_id) %>% 
    as.matrix()
  
  # Rename counts_matrix rows
  rownames(counts_matrix) <- counts %>% 
    .[["junction_id"]]
  
  #----Create colData-------------
  
  # All samples in counts matrix also in sample info?
  sample_match <- all(colnames(counts_matrix) %in% c(sample_info %>% 
                                                       dplyr::filter(!!rlang::sym(sample_id_col) %in% colnames(counts_matrix)) %>% 
                                                       .[[sample_id_col]]))
  
  # Samples in counts matrix in same order as in sample info?
  sample_order <- all(colnames(counts_matrix) %in% c(sample_info %>% 
                                                       dplyr::filter(!!rlang::sym(sample_id_col) == colnames(counts_matrix)) %>% 
                                                       .[[sample_id_col]]))
  
  # Check sample names match in metadata and counts file and are also in the correct order
  if(sample_match == "TRUE"){
    if(sample_order == "FALSE"){
      
      ids <- match(c(colnames(counts_matrix)),
                   c(sample_info %>% 
                       .[[sample_id_col]]))
      
      sample_info <- sample_info[ids,]
      
    } else{
      print("Order of samples in counts matrix matches order of samples in sample info.")
    }
    
  } else{
    stop("Not all samples in counts matrix are present in the sample info.")
  }
  
  colData <- sample_info %>% 
    as.data.frame() %>% 
    dplyr::filter(!!rlang::sym(sample_id_col) %in% colnames(counts_matrix)) 
  
  # Add row names to colData dataframe
  row.names(colData) <- colData %>% .[[sample_id_col]]
  
  rse_jx <- SummarizedExperiment(assays = list(counts = counts_matrix),
                                 rowRanges = row_ranges,
                                 colData = colData %>% 
                                   dplyr::select(-!!sample_id))
  
}