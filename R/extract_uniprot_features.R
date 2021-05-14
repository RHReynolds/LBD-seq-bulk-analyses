# Description: Functions to extract uniprot features as granges object and to overlap with granges

# Code modified from maser::mapProteinFeaturesToEvents() - https://github.com/DiogoVeiga/maser/blob/master/R/mappingEvents.R
create_granges_uniprot <- function(tracks, by = c("feature", "category"), ncores = 1){
  
  df <-  maser::availableFeaturesUniprotKB()
  
  if (by == "feature") {
    if (!any(tracks %in% as.vector(df$Name))) {
      stop(cat("\"tracks\" arg is invalid."))
    }
    features <- tracks
  }
  Category <- NULL
  if (by == "category") {
    if (!any(tracks %in% as.vector(df$Category))) {
      stop(cat("\"tracks\" arg is invalid."))
    }
    df_filt <- dplyr::filter(df, Category %in% tracks)
    features <- as.vector(df_filt$Name)
  }
  
  features_Gr <- mclapply(features,  createGRangesUniprotKBtrack, 
                          mc.cores = ncores)
  
  # Remove class "try-error"
  class <- lapply(features_Gr, class) %>% unlist()
  
  features_Gr <- GRangesList(features_Gr[which(class == "GRanges")])
  names(features_Gr) <- features[which(class == "GRanges")]
  
  return(features_Gr)
  
}

# From maser - https://github.com/DiogoVeiga/maser/blob/master/R/queryUniprotKB.R

createGRangesUniprotKBtrack <- function(track_name){
  
  Name <- NULL
  track_df <- urlTracksUniprotKB()
  
  if(!track_name %in% track_df$Name){
    stop(cat("Unknown track name."))
  }
  
  track <- dplyr::filter(track_df, Name %in% track_name)
  # bed <- read.table(as.character(track$URL), header = FALSE, sep = "\t", 
  #                   quote = NULL,
  #                   stringsAsFactors = FALSE)
  
  bed <- data.table::fread(as.character(track$URL), header = FALSE, sep = "\t", 
                           quote = "", stringsAsFactors = FALSE,
                           data.table = FALSE, showProgress = FALSE)
  
  colnames(bed)[1:6] <- c("chr", "start", "end", "Uniprot_ID", "V5", "strand") 
  res <- strsplit(bed$V14, ";")
  
  name <- vapply(seq_along(res), function(i){
    
    if (length(unlist(res[i])) > 1){
      #return(paste0(bed$Uniprot_ID[i], ":", res[[i]][[2]]))
      return(paste0(bed$Uniprot_ID[i], ":", bed$V14[i]))
    }else{
      return(paste0(bed$Uniprot_ID[i], ":", "NA"))
    }
    
  }, character(1)
  )
  
  bed <- cbind(bed, Name = name)
  bed.gr <- methods::as(bed, "GRanges")
  
  GenomeInfoDb::genome(bed.gr) <- "hg38"
  
  return(bed.gr)
  
}


# From maser - https://github.com/DiogoVeiga/maser/blob/master/R/queryUniprotKB.R
urlTracksUniprotKB <- function(){
  
  trackMetadata <- paste0("ftp://ftp.uniprot.org/pub/databases/uniprot/",
                          "current_release/knowledgebase/genome_annotation_tracks/",
                          "UP000005640_9606_tracks.txt")
  
  data <- readLines(trackMetadata)
  
  track_df <- data.frame()
  
  track_meta <- lapply(seq_along(data), function(i){
    
    #read 1st line metadata
    aux <- gsub("\"", "", data[(i*2)-1])
    tokens <- strsplit(aux, split = "=", fixed = FALSE)
    values <- tokens[[1]]
    
    #read 2nd line metadata
    aux <- gsub("\"", "", data[i*2])
    tokens <- strsplit(aux, split = "=", fixed = FALSE)
    values2 <- tokens[[1]]
    
    trackName <- gsub(" description", "", values[2])
    trackName <- gsub("UniProtKB ", "", trackName)
    
    trackDesc <- gsub(" type", "", values[3])
    
    trackFolder <- gsub(" url", "", values[8])
    trackFolder <- gsub("_hub", "_beds", trackFolder)
    trackFolder <- gsub("/hg38", "", trackFolder)
    
    trackFile <- gsub(" description", "", values2[2])
    trackFile <- gsub(".bb", ".bed", trackFile)
    
    ftp <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/"
    trackUrl <- paste0(ftp, trackFolder, "/", trackFile)
    
    return(data.frame(Name = trackName, Description = trackDesc, 
                      URL = as.character(trackUrl)))
    
  })
  track_df <- do.call(rbind, track_meta)
  track_df <- dplyr::filter(track_df, !is.na(Name))
  
  return(track_df)
}

overlap_granges_with_uniprot <- function(junc_metadata, uniprot){
  
  # Ensure same seqnames styles
  if(!any(GenomeInfoDb::seqlevelsStyle(junc_metadata) %in% GenomeInfoDb::seqlevelsStyle(uniprot))){
    
    GenomeInfoDb::seqlevelsStyle(uniprot) <- GenomeInfoDb::seqlevelsStyle(junc_metadata)[1]
    
  }
  
  junc_metadata_df <- junc_metadata %>% 
    as_tibble() %>% 
    dplyr::select(intron_chr = seqnames, intron_start = start, intron_end = end, intron_width = width, intron_strand = strand, 
                  cluster_id, contains("gene"), contains("transcript"), contains("exon"), junc_cat)
  
  overlapping_features <- uniprot %>% 
    lapply(., function(gr){
      
      overlaps <-
        GenomicRanges::findOverlaps(junc_metadata, gr, type = "within")
      
      
      gr[subjectHits(overlaps)] %>% 
        as_tibble() %>% 
        dplyr::select(feature_chr = seqnames, feature_start = start, feature_end = end, feature_width = width,
                      feature_strand = strand, uniprot_id = Uniprot_ID, feature_name = Name) %>% 
        dplyr::bind_cols(junc_metadata_df[queryHits(overlaps),])
      
      
      
    }) %>% 
    qdapTools::list_df2df(col1 = "feature_category")
  
  return(overlapping_features)
  
}

