#' Extract intron and adjoining exon co-ordinates from junction co-ordinates.
#'
#' Function to extract intron and adjoining exon co-ordinates from
#' junction/intron(s) annotated using \code{dasper::annotate_junc_ref}. For
#' those junctions/introns with start/end positions that overlap more than one
#' exon boundary, the longest exon will be used. For those junctions/introns
#' with start/end positions that do not overlap an annotated exon, the user can
#' specify the region that should be used to define the "adjoining exon".
#'
#' @author David Zhang
#' Original author: David Zhang Original code adapted to choose largest
#' associated exon (previous smallest). Function name changed (originally
#' \code{get_regions_to_obtain_cov}), in addition to names of some arguments.
#'
#' @param junc_metadata gr. Junction metadata output of
#'   \code{dasper::annotate_junc_ref}. Must also include a metadata column that
#'   uniquely indexes each row in the gr.
#' @param index_colname chr. Name of metadata column in \code{junc_metadata}
#'   that uniquely indexes each row in the gr.
#' @param ref_gtf gtf with exons loaded using \code{rtracklayer}.
#' @param novel_exon_size int. Size of "exon" to add when an end of a junction
#'   does not overlap an annotated exon. Default is 20.
#'
#' @return df. Dataframe with columns:
#' @export

get_adjoining_exons <- function(junc_metadata, index_colname, ref_gtf, novel_exon_size = 20){
  
  ref_gtf_df <- ref_gtf %>% 
    as.data.frame() %>% 
    tibble::as_tibble() 
  
  # get all junc indexes and their exons/genes in a long table format
  junc_exon_intron_gene_coords <- 
    tibble(index = junc_metadata %>% mcols() %>% .[[index_colname]],
           chr = seqnames(junc_metadata) %>% as.character(), 
           start_intron = start(junc_metadata), 
           end_intron = end(junc_metadata)) %>% 
    dplyr::full_join(.numlist_to_df(gr = junc_metadata, colname = "exon_id_start")) %>% 
    dplyr::full_join(.numlist_to_df(gr = junc_metadata, colname = "exon_id_end")) %>% 
    dplyr::full_join(.numlist_to_df(gr = junc_metadata, colname = "gene_id_junc"))
  
  # keep only the largest associated exons
  exon_widths_uniq <- 
    ref_gtf_df %>% 
    dplyr::filter(type == "exon", !duplicated(exon_id)) %>% 
    dplyr::select(exon_id, width, start, end)
  
  junc_exon_intron_gene_coords <- 
    junc_exon_intron_gene_coords %>% 
    .keep_largest_exon(exon_widths_uniq, col_to_join_by = "exon_id_start") %>% 
    .keep_largest_exon(exon_widths_uniq, col_to_join_by = "exon_id_end") 
  
  # fill in the NA exon coords with a putative novel exon coords
  junc_exon_intron_gene_coords <- 
    junc_exon_intron_gene_coords %>% 
    dplyr::left_join(exon_widths_uniq %>% dplyr::select(exon_id, start, end), by = c("exon_id_start" = "exon_id")) %>% 
    dplyr::left_join(exon_widths_uniq %>% dplyr::select(exon_id, start, end), by = c("exon_id_end" = "exon_id"), suffix = c("_exon_up", "_exon_down")) %>% 
    dplyr::mutate(end_exon_up = ifelse(is.na(exon_id_start), start_intron - 1, end_exon_up), 
                  start_exon_up = ifelse(is.na(exon_id_start), end_exon_up - (novel_exon_size - 1), start_exon_up), 
                  start_exon_down = ifelse(is.na(exon_id_end), end_intron + 1, start_exon_down), 
                  end_exon_down = ifelse(is.na(exon_id_end), start_exon_down + (novel_exon_size - 1), end_exon_down)) %>% 
    dplyr::select(index, chr, start_exon_up, end_exon_up, start_intron, end_intron, start_exon_down, end_exon_down, gene_id_junc)
  
  return(junc_exon_intron_gene_coords)
  
}



#' Unlist list-like column in granges object
#'
#' @author David Zhang
#'
#' @param gr granges object
#' @param colname chr. Column to unlist.
#'
#' @return tibble with unlisted items that can be joined back to original
#'   granges object using the column "index"
#' @export
#' 

.numlist_to_df <- function(gr, colname){
  
  numlist <- gr %>% mcols() %>% .[[colname]]
  
  numlist_unlisted <- numlist %>% unlist()
  
  df <- 
    tibble(index = names(numlist_unlisted) %>% as.integer(), 
           !!colname := numlist_unlisted)
  
  return(df)
  
}

#' Function to keep largest exon where more than one exon id is associated with
#' an intron start/end.
#'
#' @author David Zhang
#' Original author: David Zhang Original code adapted to choose largest
#' associated exon (previous smallest).
#'
#' @param junc_exon_intron_gene_coords tibble with intron co-ordinates, exons
#'   ids at start/end of intron and associated gene id
#' @param exon_widths_uniq tibble with unique exon ids, their width, and their
#'   start/end co-ordinates
#' @param col_to_join_by chr. Name of column in join
#'   junc_exon_intron_gene_coords to join by.
#'
#' @return \code{junc_exon_intron_gene_coords} returned with largest exon in
#'   those cases where more than one exon id assoicated with start/end

.keep_largest_exon <- function(junc_exon_intron_gene_coords, exon_widths_uniq, col_to_join_by){
  
  col_to_join <- "exon_id"
  names(col_to_join) <- col_to_join_by
  
  junc_exon_intron_gene_coords <- 
    junc_exon_intron_gene_coords %>% 
    dplyr::left_join(exon_widths_uniq %>% dplyr::select(exon_id, width), by = col_to_join) %>% 
    dplyr::mutate(width = ifelse(is.na(width), 0, width)) %>% # changing NAs to 0 so NA juncs are not removed by the filter 
    dplyr::group_by(index) %>% 
    dplyr::filter(width == max(width), # filter for max width i.e. largest exon
           !duplicated(width)) %>% # sometimes exons have the same width, so arbitrarily keeping one
    dplyr::select(-width) %>% # remove width column in the end
    ungroup()
  
  return(junc_exon_intron_gene_coords)
  
}



#' Get proximal intronic regions
#'
#' Function extract proximal intronic regions from a set of introns. For those
#' introns with length <= 1000, the entire intronic region will be used +/- the
#' user-defined exonic region. For those introns with length > 1000, proximal
#' intronic regions will be generated for the start and end site of the intron.
#'
#' @param intron_exon_coordinates df. Output of \code{get_adjoining_exons()},
#'   with co-ordinates for intron start/end, and start/end co-ordinates of
#'   upstream and downstream exons.
#' @param exonic_region int. Number of base pairs to conisder as "exonic"
#'   region. Default is 50, as based on definition in
#'   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-01982-9.
#' @param intronic_region int. Number of base pairs to conisder as "intronic"
#'   region. Default is 500, as based on definition in
#'   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-01982-9.
#'
#'
#' @return Dataframe with co-ordinates for proximal intronic regions.
#' @export
#'
#' 

get_proximal_intronic_regions <- function(intron_exon_coordinates, exonic_region = 50, intronic_region = 500){
  
  # Generate dataframe with intron lengths
  intron_length <- intron_exon_coordinates %>% 
    dplyr::mutate(intron_length = end_intron - start_intron,
                  cluster_id_start_end = str_c(cluster_id, start_intron, end_intron, sep = ":")) %>% 
    dplyr::select(-contains("exon_up"), -contains("exon_down"))
  
  # For introns <= 1000 bp, use entire intronic region +/- exonic regions to define proximal intronic region
  intron_under_1000bp <- intron_length %>% 
    dplyr::filter(intron_length <= 1000) %>% 
    dplyr::mutate(start_prox = start_intron - exonic_region,
                  end_prox = end_intron + exonic_region,
                  prox_prime_position = NA) %>% 
    dplyr::select(intron_exon_index, cluster_id, chr, start_intron, end_intron, start_prox, end_prox, gene_id_junc, prox_prime_position)
  
  # For introns > 1000 bp, calculate a proximal intronic region for the start and end of the intron
  # Add whether these proximal intronic regions are 5' or 3' with respect to the intron (as dependent on the strand on which the intron is located)
  intron_over_1000bp <- intron_length %>% 
    dplyr::filter(intron_length > 1000) %>% 
    tidyr::gather(key = "intron_location", value = "position", start_intron, end_intron) %>% 
    dplyr::arrange(cluster_id_start_end) %>% 
    dplyr::mutate(start_prox = case_when(intron_location == "start_intron" ~ position - exonic_region,
                                         intron_location == "end_intron" ~ position - intronic_region),
                  end_prox = case_when(intron_location == "start_intron" ~ position + intronic_region,
                                       intron_location == "end_intron" ~ position + exonic_region),
                  prox_prime_position = case_when(str_detect(cluster_id, "\\+") & intron_location == "start_intron" ~ "five_prime",
                                                  str_detect(cluster_id, "\\+") & intron_location == "end_intron" ~ "three_prime",
                                                  str_detect(cluster_id, "-") & intron_location == "start_intron" ~ "three_prime",
                                                  str_detect(cluster_id, "-") & intron_location == "end_intron" ~ "five_prime")) %>% 
    tidyr::separate(cluster_id_start_end, into = c(NA, NA, "start_intron", "end_intron"), sep = ":", remove = T) %>% 
    dplyr::mutate(start_intron = as.integer(start_intron),
                  end_intron = as.integer(end_intron)) %>% 
    dplyr::select(intron_exon_index, cluster_id, chr, start_intron, end_intron, start_prox, end_prox, gene_id_junc, prox_prime_position)
  
  # Bind both dataframes together
  proximal_regions <- intron_under_1000bp %>% 
    dplyr::bind_rows(intron_over_1000bp) %>% 
    dplyr::arrange(intron_exon_index)

  return(proximal_regions)
  
  
}
  

