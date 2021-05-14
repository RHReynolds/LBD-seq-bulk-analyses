#' Convert leafcutter intron output to input that can be used with \code{annotate_junc_ref}
#'
#' \code{leafcutter_results} takes as input the output of a leafcutter
#' differential splicing analysis (i.e. data from outputted
#' "*_effect_sizes.txt") and converts it to a granges object that can be used
#' either in the \code{dasper::annotate_junc_ref} if users wish to annotate their
#' junctions, or \code{dasper::plot_sashimi} to visualise differential splicing events.
#'
#' @param leafcutter_results dataframe. Output of leafcutter differential
#'   splicing analysis from "*_effect_sizes.txt" file.
#' @param use_strand logical. Whether to use strand information. Default = TRUE.
#'
#' @return gr. Leafcutter output in granges format.
#' @export

convert_leafcutter <- function(leafcutter_results, use_strand = FALSE){

  # Convert leafcutter results to grange object
  leafcutter_gr <- leafcutter_results %>%
    tidyr::separate(col = "intron",
                    into = c("chr", "intron_start", "intron_end", "cluster_id"),
                    sep = ":") %>%
    tidyr::separate(., col = "cluster_id",
                    into = c("prefix", "cluster_id", "strand"),
                    sep = "_") %>%
    # Conditional pipe that removes strand column if use_strand == FALSE
    purrr::when(
      use_strand == FALSE ~ dplyr::select(., -strand),
      ~ .
    ) %>%
    purrr::when(
      use_strand == FALSE ~ tidyr::unite(., col = "cluster_id", c("prefix", "cluster_id"), sep = "_", remove = FALSE),
      ~ tidyr::unite(., col = "cluster_id", c("prefix", "cluster_id", "strand"), sep = "_", remove = FALSE)
    ) %>%
    dplyr::select(-prefix) %>%
    dplyr::mutate(cluster_id = str_c(chr, ":", cluster_id),
                  chr = str_replace(chr, "chr", "")) %>%
    # Conditional pipe to include/exclude strand field in makeGRangesFromDataFrame
    purrr::when(
      use_strand == FALSE ~ GenomicRanges::makeGRangesFromDataFrame(.,
                                                                    start.field = "intron_start",
                                                                    end.field = "intron_end",
                                                                    seqnames.field = "chr",
                                                                    keep.extra.columns = TRUE),
      ~ GenomicRanges::makeGRangesFromDataFrame(.,
                                                start.field = "intron_start",
                                                end.field = "intron_end",
                                                seqnames.field = "chr",
                                                strand.field = "strand",
                                                keep.extra.columns = TRUE)
    )

  # Leafcutter intron definition adds an extra +1 to intron ends.
  # For best matching to ref. annotation must remove 1 bp
  GenomicRanges::end(leafcutter_gr) <- GenomicRanges::end(leafcutter_gr) - 1

  return(leafcutter_gr)

}
