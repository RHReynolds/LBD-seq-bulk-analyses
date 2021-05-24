#' Function to filter ATtRACT database for RBPs of interest.
#'
#' @param path_to_attract_db Chr. Path to "ATtRACT_db.txt" file downloaded from:
#'   https://attract.cnic.es/download
#' @param organism_filter Chr. Organism to filter for. Default is "Homo_sapiens"
#' @param length_filter Int. Length of motif to filter for. Default is length >=
#'   5.
#' @param score_filter Chr. Quality score to filter for. Default is "1.000000**"
#' @param GO_filter Chr. GO term
#' @param GO_type
#'
#' @return Returns a filtered version of the ATtract_db.txt file in dataframe
#'   format.
#' @export
#' 

attract_db_filter <- function(path_to_attract_db, 
                              organism_filter = "Homo_sapiens", 
                              length_filter = 7,
                              score_filter = "1.000000**",
                              GO_filter = NULL,
                              GO_type = c("bp", "cc", "mf")){
  
  
  
  attract_db <- read.table(path_to_attract_db, header=TRUE, sep="\t")
  
  # Filters: human, motif length >6
  attract_db_filtered <-  attract_db %>% 
    dplyr::filter(Organism %in% organism_filter, Len >= length_filter, Score %in% score_filter) %>%
    dplyr::group_by(Gene_id) %>% 
    dplyr::slice(which.max(Len))
  
  if(!is.null(GO_filter)){
    
    library(GO.db)
    
    # Get associated offspring terms
    if(GO_type == "bp"){
      
      GO_offspring <- GOBPOFFSPRING[[GO_filter]]
      
    } else if(GO_type == "cc"){
      
      GO_offspring <- GOCCOFFSPRING[[GO_filter]]
      
    } else{
      
      GO_offspring <- GOMFOFFSPRING[[GO_filter]]
      
    }
    
    if(length(GO_offspring) < 1){
      
      print(str_c("No OFFSPRING associated with: ", GO_filter, ". Using this GO term for biomart query."))
      GO_offspring <- GO_filter
      
    } 
    print(str_c("Using OFFSPRING GO terms from: ", GO_filter, ". Number of offspring: ", length(unique(GO_offspring))))
    
    GO_genes <- query_biomart(mart = 38,
                              ensembl_version = "http://jul2019.archive.ensembl.org", #v97
                              attributes=c("hgnc_symbol", "ensembl_gene_id", "go_id"),
                              filter = c("go"), 
                              values = GO_offspring)
    
    print(str_c("Number of unique hgnc_symbols found: ", length(unique(GO_genes$hgnc_symbol))))
    
    attract_db_filtered <- attract_db_filtered %>% 
      dplyr::inner_join(GO_genes %>% 
                          dplyr::distinct(hgnc_symbol),
                        by = c("Gene_name" = "hgnc_symbol"))
    
  }
  
  return(attract_db_filtered)
  
}

# Load query_biomart function
source(here::here("R", "biomart_df.R"))

#' Filter PWMs and return commands that when run in command line will convert
#' PWM.txt to PWM.meme.
#'
#' @param db_list List. Output of \code{attract_db_filter} in a named list
#'   format.
#' @param path_to_pl_script Chr. Path to the \code{attract_pwm_human.pl} script.
#' @param path_to_pwm Chr. Path to postional weight matrices to be filtered.
#' @param output_path Chr. Path where files should be outputted.
#' @param path_to_chen2meme Chr. Path to the MEME application \code{chen2meme}
#'
#' @return Will output two types of file in the output path: (1) the filtered
#'   ATtRACT_db.txt file, with an additional suffix (as determined by the name
#'   of the element in db_list) and (2) the filter pwm.txt file, with an
#'   additional suffix (as determined by the name of the element in db_list).
#'   Will also output a list of commands that can be used to convert these
#'   pwm.txt files into the correct .meme format.
#' @export
#' 

filter_pwm_and_convert <- function(db_list,
                                   path_to_pl_script,
                                   path_to_pwm,
                                   output_path,
                                   path_to_chen2meme){
  
  commands_list <- vector(mode = "list", length = length(db_list))
  names(commands_list) <- names(db_list)
  
  for(i in 1:length(db_list)){
    
    file_suffix <- names(db_list[i])
    
    write.table(db_list[[i]], 
                str_c(output_path, "ATtRACT_db_", file_suffix, ".txt"), 
                quote = FALSE, row.names=FALSE, sep="\t")
    
    perl_args <- str_c(path_to_pl_script, 
                       " -i ", path_to_pwm,
                       " -h ", output_path, "ATtRACT_db_", file_suffix, ".txt",
                       " > ", output_path,
                       "/pwm_", file_suffix, ".txt")
    
    system2(command = "perl", args = perl_args)
    
    meme_conversion_args <- str_c(path_to_chen2meme,
                                  " ", output_path, "pwm_", file_suffix, ".txt",
                                  " > ", output_path, "pwm_", file_suffix, ".meme")
    
    commands_list[[i]] <- meme_conversion_args
    
    
  }
  
  return(commands_list)
  
}


#' Generate granges objects for RBP analysis.
#'
#' For each pairwise comparison, this function will generate 3 or 6 granges
#' objects for each comparison, including: (1) intron_only: granges object using
#' start/end co-ordinates of introns; (2) intron_additional_bp: granges object
#' using start/end co-ordinates of introns +/- x bp ; (3) intron_exon: granges
#' object using start co-ordinate for upstream exon and end co-ordinate for
#' downstream exon; (4) prox_intron: granges object using start co-ordinate for
#' upstream proximal region and end co-orgindate for downstream proximal region;
#' (5) prox_intron_five_prime: granges object derived from prox_intron filtered
#' only for 5' proximal intronic regions; (6) prox_intron_three_prime: granges
#' object derived from prox_intron filtered only for 5' proximal intronic
#' regions.
#'
#' All granges object will retain comparison, cluster_id, p.adjust, logef, and
#' deltapsi as metadata columns.
#'
#' @param leafcutter_list list. Named list of leafcutter outputs, including
#'   "_cluster_significance" and "_effect_sizes" outputs (names in named list:
#'   cluster_significance = "_cluster_significance.txt"; intron_usage =
#'   "_effect_sizes.txt").
#' @param intron_exon_coordinates df. Output of \code{get_adjoining_exons()},
#'   with co-ordinates for intron start/end, and start/end co-ordinates of
#'   upstream and downstream exons.
#' @param proximal_intronic_regions df. Output of
#'   \code{get_proximal_intronic_regions()}, with co-ordinates for proximal
#'   intronic regions. Default is NULL.
#' @param additional_bp int. Number of base pairs to add on to intron start/end
#'   co-ordinates. Default is 50.
#'
#' @return List of pairwise comparisons with 3 or 6 granges objects each.
#' @export
#' 

creating_RBP_query_annotations <- function(leafcutter_list, intron_exon_coordinates, proximal_intronic_regions = NULL, additional_bp = 50){
  
  # Create merge of cluster_significance and intron usage
  clu_significance_intron_usage <-
    leafcutter_list$cluster_significance %>% 
    dplyr::filter(status == "Success",
                  !str_detect(genes, ",")) %>% 
    dplyr::mutate(cluster_id = str_replace(cluster, "chr", "")) %>% 
    dplyr::select(-cluster) %>% 
    dplyr::inner_join(leafcutter_list$intron_usage %>% 
                        tidyr::separate(intron, into = c("chr", "start", "end", "cluster") , sep = ":") %>% 
                        dplyr::mutate(chr = str_replace(chr, "chr", ""),
                                      cluster_id = str_c(chr, ":", cluster)) %>% 
                        dplyr::select(-cluster), by = c("comparison", "cluster_id"))
  
  if(!is.null(proximal_intronic_regions)){
    
    # Join intron/exon with proximal intronic regions
    intron_exon_coordinates <- intron_exon_coordinates %>% 
      dplyr::inner_join(proximal_intronic_regions)
    
  }
  
  # Recall that leafcutter adds +1 bp to intron ends, which interferes with annotation.
  # Thus, annotated clu will need to subtract 1 from intron ends in clu_significance_intron_usage to join two dataframes.
  clu_significance_intron_usage_with_coord <- 
    clu_significance_intron_usage %>% 
    dplyr::mutate(start_intron = as.integer(start),
                  end_intron = as.integer(end) - 1,
                  strand = cluster_id %>% str_replace(".*_", "")) %>% 
    dplyr::inner_join(intron_exon_coordinates,
                      by = c("cluster_id", "chr", "start_intron", "end_intron"))
  
  if(!is.null(proximal_intronic_regions)){
    clu_significance_intron_usage_with_coord <- 
      clu_significance_intron_usage_with_coord %>% 
      dplyr::select(comparison, cluster_id, p.adjust, chr, strand, start_exon_up, end_exon_up, start_intron, end_intron, 
                    start_exon_down, end_exon_down, start_prox, end_prox, prox_prime_position, logef, deltapsi, contains("index"))
    
  } else{
    
    clu_significance_intron_usage_with_coord <- 
      clu_significance_intron_usage_with_coord %>% 
      dplyr::select(comparison, cluster_id, p.adjust, chr, strand, start_exon_up, end_exon_up, start_intron, end_intron, 
                    start_exon_down, end_exon_down, logef, deltapsi, contains("index"))
    
  }
  
  comparisons <- clu_significance_intron_usage_with_coord$comparison %>% unique()
  RBP_lists <- vector(mode = "list", length = length(comparisons))
  names(RBP_lists) <- comparisons
  
  for(i in 1:length(RBP_lists)){
    
    # Intron only gr
    intron_only <- GenomicRanges::makeGRangesFromDataFrame(df = clu_significance_intron_usage_with_coord %>% 
                                                             dplyr::filter(comparison == names(RBP_lists)[i]) %>% 
                                                             dplyr::select(-start_exon_up, -end_exon_up, -start_exon_down, -end_exon_down, -contains("prox")) %>% 
                                                             dplyr::distinct() %>% 
                                                             dplyr::arrange(p.adjust, -abs(deltapsi)), 
                                                           seqnames.field = "chr", 
                                                           start.field = "start_intron", 
                                                           end.field = "end_intron", 
                                                           strand.field = "strand", 
                                                           keep.extra.columns = T)
    
    # Intron +/- x bp gr
    intron_additional_bp <- intron_only
    start(intron_additional_bp) <- start(intron_additional_bp) - additional_bp
    end(intron_additional_bp) <- end(intron_additional_bp) + additional_bp
    
    # Introns w. exons
    intron_exon <- GenomicRanges::makeGRangesFromDataFrame(df = clu_significance_intron_usage_with_coord %>% 
                                                             dplyr::filter(comparison == names(RBP_lists)[i]) %>% 
                                                             dplyr::select(-end_exon_up, -start_intron, -end_intron, -start_exon_down, -contains("prox")) %>%
                                                             dplyr::distinct() %>% 
                                                             dplyr::arrange(p.adjust, -abs(deltapsi)), 
                                                           seqnames.field = "chr", 
                                                           start.field = "start_exon_up", 
                                                           end.field = "end_exon_down", 
                                                           strand.field = "strand", 
                                                           keep.extra.columns = T)
    
    if(!is.null(proximal_intronic_regions)){
      
      # Proximal intronic regions
      prox_intron <- GenomicRanges::makeGRangesFromDataFrame(df = clu_significance_intron_usage_with_coord %>% 
                                                               dplyr::filter(comparison == names(RBP_lists)[i]) %>% 
                                                               dplyr::select(-start_exon_up, -end_exon_up, -start_intron, 
                                                                             -end_intron, -start_intron, -start_exon_down, -end_exon_down) %>% 
                                                               dplyr::group_by(comparison, cluster_id, start_prox, end_prox, prox_prime_position) %>% 
                                                               dplyr::top_n(., n = 1, abs(deltapsi)) %>% 
                                                               dplyr::arrange(p.adjust, -abs(deltapsi)), 
                                                             seqnames.field = "chr", 
                                                             start.field = "start_prox", 
                                                             end.field = "end_prox", 
                                                             strand.field = "strand", 
                                                             keep.extra.columns = T)
      
      # 5' proximal intronic regions
      prox_intron_five_prime <- prox_intron[!is.na(prox_intron$prox_prime_position)]
      prox_intron_five_prime <- prox_intron_five_prime[prox_intron_five_prime$prox_prime_position == "five_prime"]
      mcols(prox_intron_five_prime) <- 
        mcols(prox_intron_five_prime)[, c(!colnames(mcols(prox_intron_five_prime)) %in% c("prox_prime_position"))]
      
      # 3' proximal intronic regions
      prox_intron_three_prime <- prox_intron[!is.na(prox_intron$prox_prime_position)]
      prox_intron_three_prime <- prox_intron_three_prime[prox_intron_three_prime$prox_prime_position == "three_prime"]
      mcols(prox_intron_three_prime) <- 
        mcols(prox_intron_three_prime)[, c(!colnames(mcols(prox_intron_three_prime)) %in% c("prox_prime_position"))]
      
      
      RBP_lists[[i]] <- setNames(list(intron_only, intron_additional_bp, intron_exon, 
                                      prox_intron, prox_intron_five_prime, prox_intron_three_prime),
                                 c("intron_only", "intron_additional_bp", "intron_exon",
                                   "prox_intron", "prox_intron_five_prime", "prox_intron_three_prime"))
      
    } else{
      
      RBP_lists[[i]] <- setNames(list(intron_only, intron_additional_bp, intron_exon),
                                 c("intron_only", "intron_additional_bp", "intron_exon"))
      
    }
    
  }
  
  return(RBP_lists)
  
}



#' Generate fasta from granges object.
#'
#' @param gr gr. Granges object
#' @param genome BSgenome object for full masked genome sequence e.g.
#'   BSgenome.Hsapiens.UCSC.hg38.masked.
#' @param style chr. Vector with the style of chromosome to be used. Default is
#'   "UCSC". See GenomeInfoDb::mapSeqlevels for more details.
#' @param remove_star_strand logical. Function to extract sequence cannot handle
#'   star (*) strand. Thus, if set to TRUE, any entires with star strand will be
#'   removed. If set to FALSE, start strand will be converted to '+'. Default =
#'   TRUE.
#' @param output_path chr. Path where file should be saved.
#' @param file_name chr. File name.
#' @param reduce logical. This specifies whether granges object should reduce
#'   any overlapping ranges. Default is TRUE.
#'
#' @return
#' @export
#' 

generate_fasta_from_gr <- function(gr, genome, style = "UCSC", remove_star_strand = T, output_path, file_name, reduce = T){
  
  if(reduce == TRUE){
    
    comparison <- gr$comparison %>% unique()
    
    gr <- gr %>% 
      GenomicRanges::reduce()
    
    gr$comparison <- comparison
    
  }
  
  # adding "chr" in front of seqnames
  newStyle <- GenomeInfoDb::mapSeqlevels(seqlevels(gr), style = style)
  gr <- GenomeInfoDb::renameSeqlevels(gr, newStyle)
  
  # BSgenome::getSeq does not run if gr strand includes "*"
  # Either remove these entries or convert to "+"
  if(remove_star_strand == T){
    
    gr <- gr[!strand(gr) == "*"]
    
  } else{
    
    strand(gr) <- ifelse(strand(gr) == "*", "+", strand(gr))
    
  }
  
  # Extract sequences (using the masked version)
  seq <- BSgenome::getSeq(genome, gr)
  
  if(reduce == TRUE){
    
    seq@ranges@NAMES <- str_c(gr$comparison, seqnames(gr), start(gr), end(gr), strand(gr), sep = ":")
    
  } else{
    
    seq@ranges@NAMES <- str_c(gr$comparison, gr$cluster_id, seqnames(gr), start(gr), end(gr), strand(gr), sep = ":")
    
  }
  
  writeXStringSet(seq, filepath = str_c(output_path, "/", file_name, ".fasta"), format="fasta")
  
}



#' Function to run fimo analysis.
#'
#' @param path_to_query_fasta_dir chr. Path to directory containing query
#'   fastas.
#' @param fasta_filter chr. Pattern to filter desired fastas. Only necessary if
#'   user only wants to run analysis on certain fastas in fasta_dir. Default =
#'   NULL.
#' @param pwm_filter chr. Pattern to filter desired pwms. Only necessary if user
#'   only wants to run analysis on certain pwms in pwm_dir. Default = NULL.
#' @param path_to_pwm_dir chr. Path to directory containing pwm_*.meme file to
#'   be used.
#' @param path_to_fimo chr. Path to fimo software.
#' @param output_dir chr. Path to output directory to be used for results of
#'   analysis. For each pwm, a directory will be created within the output
#'   directory. Within this pwm directory, for each query fasta, a directory
#'   will be created using the name of the query fasta.
#' @param cores integer. Number of cores to parallelise across. Default = 1.
#'
#' @return Fimo run on combinations of pwms and query fastas.
#' @export
#' 

run_fimo_analysis <- function(path_to_query_fasta_dir, fasta_filter = NULL, pwm_filter = NULL, path_to_pwm_dir, path_to_fimo, output_dir, cores = 1){
  
  query_fastas_df <- 
    tibble(fasta_path = list.files(path = path_to_query_fasta_dir, pattern = ".fasta", full.names = T),
           fasta_name = list.files(path = path_to_query_fasta_dir, pattern = ".fasta") %>% str_replace(., "\\.fasta", "")) %>% 
    tidyr::separate(fasta_name, into = c("comparison", "query"), sep = ":", remove = F)
  pwms_df <- 
    tibble(pwm_path = list.files(path = path_to_pwm_dir, pattern = ".meme", full.names = T),
           pwm_name = list.files(path = path_to_pwm_dir, pattern = ".meme") %>% 
             str_replace(., "\\.meme", ""))
  
  if(!is.null(fasta_filter)){
    
    query_fastas_df <- 
      query_fastas_df %>% 
      dplyr::filter(str_detect(fasta_name, fasta_filter))
    
  }
  
  if(!is.null(pwm_filter)){
    
    pwms_df <- 
      pwms_df %>% 
      dplyr::filter(str_detect(pwm_name, pwm_filter))
    
  }
  
  for(i in 1:nrow(pwms_df)){
    
    pwm <- pwms_df[i, ]
    
    pwm_dir <- str_c(output_dir, "/", pwm$pwm_name)
    
    if (!dir.exists(pwm_dir)){
      dir.create(pwm_dir)
    } else {
      print("Dir already exists!")
    }
    
    # Run in parallel
    cl <- parallel::makeCluster(cores)
    
    # Register clusters
    doParallel::registerDoParallel(cl)
    foreach::getDoParWorkers() 
    
    foreach::foreach(j = 1:nrow(query_fastas_df),
                     .verbose = TRUE,
                     .packages = c("tidyverse", "stringr")) %dopar% {
                       
                       fasta <- query_fastas_df[j, ]
                       
                       print(str_c(Sys.time(), " - running fimo on: ", pwm$pwm_name, "/", fasta$fasta_name))
                       
                       fimo_args <- str_c(" --bfile --uniform--",
                                          " -oc ", output_dir, "/", pwm$pwm_name, "/", fasta$fasta_name,
                                          " ", as.character(pwm$pwm_path),
                                          " ", as.character(fasta$fasta_path))
                       
                       system2(command = path_to_fimo, args = fimo_args)
                       
                     }
    
    # Stop cluster
    stopCluster(cl)
    
  }
  
}

#' Calculate median enrichment score per RBP across all sequences and sum
#' enrichment score per queried sequence across all RBPs.
#'
#' Function will calculate an enrichment score per RBP and sequence. Enrichment
#' score = (count of RBP motif within query sequence)/(query sequence
#' length/100). From this score, a median is calculated for each RBP across all
#' queried sequences. Further, a sum of the score is calculated for each
#' sequence across all RBPs.
#'
#' The median can be used to plot a heatmap for each sequence, while the sum can
#' be used to calculate a cumulative density function.
#'
#' Original function authored by Sid. Original code adapted to allow it to run
#' within an R session.
#'
#' @author Sid Sethi
#' @source https://github.com/sid-sethi/rbp_binding_analysis
#'
#' @param results_dir chr. Path to directory containing directories (each of
#'   which will contain a fimo.tsv file).
#' @param significance_threshold chr. How should sequences be thresholded? Users
#'   have choice of either using FDR-thresholding (e.g. q.value < 0.05), or a
#'   variable p-value threshold that is determined by sequence length.
#' @param p_value int. P-value threshold for query sequences of length < 1000
#'   bp. For sequences of length > 1000 bp, the p-value threshold will be varied
#'   according to their length. That is, the p-value threshold will =
#'   p_value/(sequence_length/1000). Default = 1e-04.
#' @param q_value int. FDR threshold to use. Default = 0.05.
#' @param core integer. Number of cores to parallelise across. Default = 1.
#'
#' @return A list with 3 lists per fimo analysis. These lists include: (1)
#'   significant_results: results filtered by p-value/q-value; (2)
#'   RBP_median_enrich: median enrichment score per RBP; (3)
#'   sequence_sum_enrich: a sum of enrichment scores per analyses sequence.
#' @export
#' 

summarise_fimo_results <- function(results_dir, significance_threshold = c("fdr", "variable"), p_value = 1e-4, q_value = 0.05, cores = 1){
  
  results_df <- tibble(file_path = list.files(results_dir, pattern = "fimo.tsv", full.names = T, recursive = T),
                       file_name = list.files(results_dir, pattern = "fimo.tsv", full.names = T, recursive = T) %>% 
                         str_replace("/fimo.tsv", "") %>% 
                         str_replace(".*/", ""))  
  
  # Run in parallel
  cl <- parallel::makeCluster(cores)
  
  # Register clusters
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers() 
  
  results_list <- foreach::foreach(i = 1:nrow(results_df),
                                   .verbose = TRUE,
                                   .packages = c("tidyverse")) %dopar% {
                                     
                                     # Read in .tsv with results
                                     results <- read.table(results_df$file_path[i] %>% as.character(), header = T) %>% 
                                       as_tibble()
                                     
                                     if(!significance_threshold %in% c("fdr", "variable")) stop("User must specify significance threshold as 'fdr' or 'variable'.")
                                     
                                     if(significance_threshold == "fdr"){
                                       
                                       results <- results %>% 
                                         tidyr::separate(sequence_name, c(NA, NA, "seq_start", "seq_end", "seq_strand"), sep=":", remove = FALSE, convert = TRUE) %>% 
                                         dplyr::mutate(same_strand = case_when(seq_strand == strand ~ TRUE,
                                                                               TRUE ~ FALSE)) %>%
                                         dplyr::filter(same_strand == TRUE) %>%
                                         dplyr::mutate(motif_len = (stop - start) + 1,
                                                       sequence_len = (seq_end - seq_start) + 1) %>% 
                                         dplyr::filter(q.value < q_value) %>% 
                                         dplyr::select(-c("seq_start", "seq_end", "seq_strand", "same_strand"))
                                       
                                     }
                                     
                                     if(significance_threshold == "variable"){
                                       
                                       results <- results %>% 
                                         tidyr::separate(sequence_name, c(NA, NA, "seq_start", "seq_end", "seq_strand"), sep=":", remove = FALSE, convert = TRUE) %>% 
                                         dplyr::mutate(same_strand = case_when(seq_strand == strand ~ TRUE,
                                                                               TRUE ~ FALSE)) %>%
                                         dplyr::filter(same_strand == TRUE) %>%
                                         dplyr::mutate(motif_len = (stop - start) + 1,
                                           sequence_len = (seq_end - seq_start) + 1,
                                           p_value_threshold = case_when(sequence_len < 1000 ~ p_value,
                                                                         sequence_len > 1000 ~ (p_value/(sequence_len/1000))),
                                           significant = case_when(p.value < p_value_threshold ~ TRUE,
                                                                   TRUE ~ FALSE)) %>% 
                                         dplyr::filter(significant == TRUE) %>% 
                                         dplyr::select(-c("p_value_threshold", "significant", "seq_start", "seq_end", "seq_strand", "same_strand"))
                                       
                                     }
                                     
                                     # Count number of times a motif of length unique(motif_len) is observed in a sequence
                                     results_groups <- results %>% 
                                       dplyr::group_by(sequence_name, motif_id) %>% 
                                       dplyr::summarise(count = n(),
                                                        motif_len = unique(motif_len))
                                     
                                     # Calculate enrichment score
                                     enrich_score <- results_groups %>% 
                                       tidyr::separate(sequence_name, c(NA, NA, "start", "end", NA), sep=":", remove = FALSE, convert = TRUE) %>% 
                                       dplyr::mutate(sequence_len = (end-start)+1,
                                                     enrich_score = count/(sequence_len/100)) %>% 
                                       dplyr::select(-c("start", "end"))
                                     
                                     # Calculate median enrichment scores for each RBP
                                     enrich_score_median <- enrich_score %>% 
                                       dplyr::group_by(motif_id) %>% 
                                       dplyr::summarise(median_enrich_score = median(enrich_score))
                                     
                                     # Sum enrichment score across each sequence
                                     enrich_score_sum <- enrich_score %>% 
                                       dplyr::group_by(sequence_name) %>% 
                                       dplyr::summarise(sum_enrich_score = sum(enrich_score)) %>% 
                                       dplyr::inner_join(enrich_score %>% 
                                                           dplyr::distinct(sequence_name, sequence_len))
                                     
                                     setNames(list(results, enrich_score_median, enrich_score_sum),
                                              c("significant_results", "RBP_median_enrich", "sequence_sum_enrich"))
                                     
                                   }
  
  # Stop cluster
  stopCluster(cl)
  
  names(results_list) <- results_df$file_name
  
  return(results_list)
  
}

#' Calculate enrichment score per RBP and sequence and filter for select RBPs.
#'
#' Function will calculate an enrichment score per RBP and sequence. Enrichment
#' score = (count of RBP motif within query sequence)/(query sequence
#' length/100). Thereafter, it will filter for genes supplied by the user.
#'
#' The output can be used to plot ecdf/density for individual RBPs (as opposed
#' to summarising enrichment score across all RBPs).
#'
#' Original function authored by Sid. Original code adapted to allow filtering
#' for individual genes.
#'
#' @author Sid Sethi
#' @source https://github.com/sid-sethi/rbp_binding_analysis
#'
#' @param results_dir chr. Path to directory containing directories (each of
#'   which will contain a fimo.tsv file).
#' @param significance_threshold chr. How should sequences be thresholded? Users
#'   have choice of either using FDR-thresholding (e.g. q.value < 0.05), or a
#'   variable p-value threshold that is determined by sequence length.
#' @param p_value int. P-value threshold for query sequences of length < 1000
#'   bp. For sequences of length > 1000 bp, the p-value threshold will be varied
#'   according to their length. That is, the p-value threshold will =
#'   p_value/(sequence_length/1000). Default = 1e-04.
#' @param q_value int. FDR threshold to use. Default = 0.05.
#' @param genes chr. Vector of genes (either hgnc symbols or ensembl ids) to
#'   filter results for.
#' @param core integer. Number of cores to parallelise across. Default = 1.
#'
#' @return A list with one dataframe per fimo analysis. This dataframe will
#'   contain enrichment scores for each sequence and RBP.
#' @export
#' 

summarise_fimo_results_for_individ_RBPs <- function(results_dir, 
                                                    significance_threshold = c("fdr", "variable"), 
                                                    p_value = 1e-4, 
                                                    q_value = 0.05, 
                                                    genes = NULL,
                                                    cores = 1){
  
  results_df <- tibble(file_path = list.files(results_dir, pattern = "fimo.tsv", full.names = T, recursive = T),
                       file_name = list.files(results_dir, pattern = "fimo.tsv", full.names = T, recursive = T) %>% 
                         str_replace("/fimo.tsv", "") %>% 
                         str_replace(".*/", ""))  
  
  # Run in parallel
  cl <- parallel::makeCluster(cores)
  
  # Register clusters
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers() 
  
  results_list <- foreach::foreach(i = 1:nrow(results_df),
                                   .verbose = TRUE,
                                   .packages = c("tidyverse")) %dopar% {
                                     
                                     # Read in .tsv with results
                                     results <- read.table(results_df$file_path[i] %>% as.character(), header = T) %>% 
                                       as_tibble()
                                     
                                     if(!significance_threshold %in% c("fdr", "variable")) stop("User must specify significance threshold as 'fdr' or 'variable'.")
                                     
                                     if(significance_threshold == "fdr"){
                                       
                                       results <- results %>% 
                                         tidyr::separate(sequence_name, c(NA, NA, "seq_start", "seq_end", "seq_strand"), sep=":", remove = FALSE, convert = TRUE) %>% 
                                         dplyr::mutate(same_strand = case_when(seq_strand == strand ~ TRUE,
                                                                               TRUE ~ FALSE)) %>%
                                         dplyr::filter(same_strand == TRUE) %>% 
                                         dplyr::mutate(motif_len = (stop - start) + 1,
                                                       sequence_len = (seq_end - seq_start) + 1) %>% 
                                         dplyr::filter(q.value < q_value) %>% 
                                         dplyr::select(-c("seq_start", "seq_end", "seq_strand", "same_strand"))
                                       
                                     }
                                     
                                     if(significance_threshold == "variable"){
                                       
                                       results <- results %>% 
                                         tidyr::separate(sequence_name, c(NA, NA, "seq_start", "seq_end", "seq_strand"), sep=":", remove = FALSE, convert = TRUE) %>% 
                                         dplyr::mutate(same_strand = case_when(seq_strand == strand ~ TRUE,
                                                                               TRUE ~ FALSE)) %>%
                                         dplyr::filter(same_strand == TRUE) %>% 
                                         dplyr::mutate(motif_len = (stop - start) + 1,
                                                       sequence_len = (seq_end - seq_start) + 1,
                                                       p_value_threshold = case_when(sequence_len < 1000 ~ p_value,
                                                                                     sequence_len > 1000 ~ (p_value/(sequence_len/1000))),
                                                       significant = case_when(p.value < p_value_threshold ~ TRUE,
                                                                               TRUE ~ FALSE)) %>% 
                                         dplyr::filter(significant == TRUE) %>% 
                                         dplyr::select(-c("p_value_threshold", "significant", "seq_start", "seq_end", "seq_strand", "same_strand"))
                                       
                                     }
                                     
                                     # Count number of times a motif of length unique(motif_len) is observed in a sequence
                                     results_groups <- results %>% 
                                       dplyr::group_by(sequence_name, motif_id) %>% 
                                       dplyr::summarise(count = n(),
                                                        motif_len = unique(motif_len))
                                     
                                     # Calculate enrichment score
                                     enrich_score <- results_groups %>% 
                                       tidyr::separate(sequence_name, c(NA, NA, "start", "end", NA), sep=":", remove = FALSE, convert = TRUE) %>% 
                                       tidyr::separate(motif_id, into = c("hgnc_symbol", "ensembl_id", "matrix_id"), sep = ":", remove = T) %>% 
                                       dplyr::mutate(sequence_len = (end-start)+1,
                                                     enrich_score = count/(sequence_len/100)) %>% 
                                       dplyr::select(-c("start", "end")) 
                                     
                                     if(any(str_detect(genes, "ENSG"))){
                                       
                                       enrich_score %>% 
                                         dplyr::filter(ensembl_id %in% genes)
                                       
                                     } else{
                                       
                                       enrich_score %>% 
                                         dplyr::filter(hgnc_symbol %in% genes)
                                       
                                     }
                                     
                                   }
  
  # Stop cluster
  stopCluster(cl)
  
  names(results_list) <- results_df$file_name
  
  return(results_list)
  
}

#' ECDF plot of enrichment scores across sequences.
#'
#' @param results_df df. Dataframe version of sequence_sum_enrich df outputted
#'   by \code{summarise_fimo_results}.
#' @param facet_var var. Variable(s) to facet by. Provide in the format vars(x,
#'   y).
#' @param colour_var chr. Variable to colour by.
#' @param colour_palette chr. Colours to use for colour palette
#' @param ... var. Variable(s) to group by, in order to calculate median per
#'   group. List column names to be used, without quotes, as is typical of
#'   tidyverse format.
#'
#' @return ECDF plot.
#' @export
#' 

fimo_ecdf_plot <- function(results_df, facet_var, colour_var, colour_palette, ...){
  
  ecdf_medians <- 
    results_df %>% 
    dplyr::group_by(...) %>% 
    dplyr::summarise(median = median(sum_enrich_score))
  
  p <- results_df %>%
    ggplot(., aes(x = sum_enrich_score, colour = .data[[colour_var]])) +
    stat_ecdf(geom = "step", lwd = 0.8, alpha =0.8, na.rm = T) +
    geom_vline(data = ecdf_medians, aes(xintercept= median, colour = .data[[colour_var]]), linetype = "dashed") +
    scale_x_continuous(name = "Number of RBP motifs per 100 nt", trans="log2") +
    labs(x = "Number of RBP motifs binding per 100 nt", y = "CDF") +
    scale_colour_manual(values = colour_palette) 
  
  if(length(facet_var) == 1){
    
    p +
      facet_wrap({{ facet_var }}) +
      theme_rhr
    
  } else if(length(facet_var) == 2){
    
    p +
      facet_grid(rows = {{ facet_var[1] }}, col = {{ facet_var[2] }}) +
      theme_rhr
    
  } else {
    
    stop("facet_var should be length 1 or 2")
    
  }
  
}

#' Density plot of enrichment scores across sequences.
#'
#' @param results_df df. Dataframe version of fimo results outputted by
#'   \code{summarise_fimo_results} or
#'   \code{summarise_fimo_results_for_individ_RBPs}.
#' @param x_var chr. Variable to facet by.
#' @param facet_var var. Variable(s) to facet by. Provide in the format vars(x,
#'   y).
#' @param colour_var chr. Variable to colour by.
#' @param colour_palette chr. Colours to use for colour palette.
#' @param scales chr. Supplied to \code{scales} argument in \code{facet_wrap}. Should scales
#'   be fixed (\code{"fixed"}, the default), free (\code{"free"}), or free in
#'   one dimension (\code{"free_x"}, \code{"free_y"})?
#' @param ... var. Variable(s) to group by, in order to calculate median per
#'   group. List column names to be used, without quotes, as is typical of
#'   tidyverse format.
#'
#' @return Density plot.
#' @export
#' 

fimo_density_plot <- function(results_df, x_var, facet_var, colour_var, colour_palette, scales = NULL, ...){
  
  medians <- 
    results_df %>% 
    dplyr::group_by(...) %>% 
    dplyr::summarise(median = median(.data[[x_var]]))
  
  results_df %>%
    ggplot(., aes(x = .data[[x_var]], colour = .data[[colour_var]], fill = .data[[colour_var]])) +
    geom_density(alpha = 0.5) +
    geom_vline(data = medians, aes(xintercept= median, colour = .data[[colour_var]]), linetype = "dashed") +
    scale_x_continuous(name = "Number of RBP motifs per 100 nt", trans="log2") +
    labs(x = "Number of RBP motifs binding per 100 nt", y = "Desnity") +
    scale_colour_manual(values = colour_palette) +
    scale_fill_manual(values = colour_palette) +
    facet_wrap({{ facet_var }}, scales = scales, ncol = 3) +
    theme_rhr
  
}

theme_rhr <-  theme_bw(base_family = "Helvetica") + 
  theme(panel.grid.major.x = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(vjust = 0.6),
        axis.title = element_text(size = 10),
        panel.spacing = unit(0.1, "lines"))

#' Function to run ame analysis.
#'
#' @param path_to_query_fasta_dir chr. Path to directory containing query
#'   fastas.
#' @param path_to_pwm_dir chr. Path to directory containing pwm_*.meme file to
#'   be used.
#' @param path_to_ame chr. Path to ame software.
#' @param output_dir chr. Path to output directory to be used for results of
#'   analysis. For each pwm, a directory will be created within the output
#'   directory. Within this pwm directory, for each query fasta, a directory
#'   will be created using the name of the query fasta.
#' @param query_name chr. Name of sequence to be used as query in AME analysis.
#' @param control_name chr. Name of sequence to be used as control in AME
#'   analysis.
#' @param scoring chr. One of c("avg", "totalhits"). Average is default. See AME
#'   for details.
#' @param cores integer. Number of cores to parallelise across. Default = 1.
#'
#' @return Fimo run on combinations of pwms and query fastas.
#' @export
#' 

run_ame_analysis <- function(path_to_query_fasta_dir, path_to_pwm_dir, path_to_ame, output_dir, query_name = "query", control_name = "control", scoring = NULL, cores = 1){
  
  query_fastas_df <- 
    tibble(fasta_path = list.files(path = path_to_query_fasta_dir, pattern = ".fasta", full.names = T, recursive = T),
           fasta_name = list.files(path = path_to_query_fasta_dir, pattern = ".fasta", recursive = T) %>% 
             str_replace(., "\\.fasta", "") %>% 
             str_replace(., ".*/", "")) %>% 
    tidyr::separate(fasta_name, into = c("comparison", "query", "list_type"), sep = ":", remove = F)
  pwms_df <- 
    tibble(pwm_path = list.files(path = path_to_pwm_dir, pattern = ".meme", full.names = T),
           pwm_name = list.files(path = path_to_pwm_dir, pattern = ".meme") %>% 
             str_replace(., "\\.meme", ""))
  
  if(is.null(scoring)){
    
    score_param <- "avg"
    
  } else{
    
    if(scoring %in% c("avg", "totalhits")){
      
      score_param <- scoring
      
    }else{
      
      print("Incorrect scoring param provided. Should be 'avg' or 'totalhits'. Defaulting to 'avg'.")
      
      score_param <- "avg"
      
    }
    
  }
  
  for(i in 1:nrow(pwms_df)){
    
    pwm <- pwms_df[i, ]
    
    pwm_dir <- str_c(output_dir, "/", pwm$pwm_name)
    
    if (!dir.exists(pwm_dir)){
      dir.create(pwm_dir)
    } else {
      print("Dir already exists!")
    }
    
    group_df <- query_fastas_df %>% 
      dplyr::distinct(comparison, query)
    
    # Run in parallel
    cl <- parallel::makeCluster(cores)
    
    # Register clusters
    doParallel::registerDoParallel(cl)
    foreach::getDoParWorkers() 
    
    foreach::foreach(j = 1:nrow(group_df),
                     .verbose = TRUE,
                     .packages = c("tidyverse", "stringr")) %dopar% {
                       
                       group <- group_df[j, ]
                       
                       fasta <- group %>% 
                         dplyr::inner_join(query_fastas_df)
                       
                       ame_args <- str_c(" --scoring ", score_param,
                                         " --oc ", output_dir, "/", pwm$pwm_name, "/", group$comparison, ":", group$query, ":", query_name, "_vs_", control_name, 
                                         " --control ", as.character(fasta %>% 
                                                                       dplyr::filter(list_type == control_name) %>% 
                                                                       .[["fasta_path"]]),
                                         " ", as.character(fasta %>% 
                                                             dplyr::filter(list_type == query_name) %>% 
                                                             .[["fasta_path"]]),
                                         " ", as.character(pwm$pwm_path))
                       
                       system2(command = path_to_ame, args = ame_args)
                       
                     }
    
    # Stop cluster
    stopCluster(cl)
    
  }
  
}

#' Summarise AME analysis results
#'
#' Function will import ame.tsv results across a results directory.
#'
#' @param results_dir chr. Path to directory containing directories (each of
#'   which will contain a ame.tsv file).
#' @param core integer. Number of cores to parallelise across. Default = 1.
#'
#' @return A list with 1 dataframe per AME analysis.
#' @export
#' 

summarise_ame_results <- function(results_dir, cores = 1){
  
  results_df <- tibble(file_path = list.files(results_dir, pattern = "ame.tsv", full.names = T, recursive = T),
                       file_name = list.files(results_dir, pattern = "ame.tsv", full.names = T, recursive = T) %>% 
                         str_replace("/ame.tsv", "") %>% 
                         str_replace(".*/", ""))  
  
  # Run in parallel
  cl <- parallel::makeCluster(cores)
  
  # Register clusters
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers() 
  
  results_list <- foreach::foreach(i = 1:nrow(results_df),
                                   .verbose = TRUE,
                                   .packages = c("tidyverse")) %dopar% {
                                     
                                     # Read in .tsv with results
                                     results <- tryCatch(read.table(results_df$file_path[i] %>% as.character(), 
                                                                    header = T), error=function(e) NULL)
                                     
                                     if(is.null(results) == TRUE){
                                       
                                       results <- tibble()
                                       
                                     } else {
                                       
                                       results %>% 
                                         as_tibble() %>% 
                                         dplyr::mutate(motif_DB = motif_DB %>% 
                                                         str_replace(".*/", "") %>% 
                                                         str_replace(".meme", "")) %>% 
                                         dplyr::rename(p_val = p.value,
                                                       bonferroni_adj_p_val = adj_p.value,
                                                       expected_n_enriched_motifs = E.value,
                                                       n_tests = tests,
                                                       percent_TP = X.TP,
                                                       percent_FP = X.FP)
                                       
                                     }
                                     
                                     
                                     
                                   }
  
  # Stop cluster
  stopCluster(cl)
  
  names(results_list) <- results_df$file_name
  
  return(results_list)
  
}

#' Summarise AME sequence analysis results
#'
#' Function will import sequences.tsv results across a results directory.
#'
#' @param results_dir chr. Path to directory containing directories (each of
#'   which will contain a sequences.tsv file).
#' @param core integer. Number of cores to parallelise across. Default = 1.
#'
#' @return A list with 1 dataframe per AME analysis.
#' @export
#' 

summarise_ame_seq_results <- function(results_dir, cores = 1){
  
  results_df <- tibble(file_path = list.files(results_dir, pattern = "sequences.tsv", full.names = T, recursive = T),
                       file_name = list.files(results_dir, pattern = "sequences.tsv", full.names = T, recursive = T) %>% 
                         str_replace("/sequences.tsv", "") %>% 
                         str_replace(".*/", ""))  
  
  # Run in parallel
  cl <- parallel::makeCluster(cores)
  
  # Register clusters
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers() 
  
  results_list <- foreach::foreach(i = 1:nrow(results_df),
                                   .verbose = TRUE,
                                   .packages = c("tidyverse")) %dopar% {
                                     
                                     # Read in .tsv with results
                                     results <- tryCatch(read.table(results_df$file_path[i] %>% as.character(), 
                                                                    header = T), error=function(e) NULL)
                                     
                                     if(is.null(results) == TRUE){
                                       
                                       results <- tibble()
                                       
                                     } else {
                                       
                                       results %>% 
                                         as_tibble() %>% 
                                         dplyr::mutate(motif_DB = motif_DB %>% 
                                                         str_replace(".*/", "") %>% 
                                                         str_replace(".meme", ""))
                                       
                                     }
                                     
                                     
                                     
                                   }
  
  # Stop cluster
  stopCluster(cl)
  
  names(results_list) <- results_df$file_name
  
  return(results_list)
  
}

