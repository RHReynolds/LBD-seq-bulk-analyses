#' Dataframe-based biomart query.
#'
#' @param dataframe Dataframe with gene names.
#' @param columnToFilter Name of column in dataframe, which contains gene names.
#' @param mart Specify genome build.
#' @param ensembl_version character scalar. If using GRCh38 genome build, can
#'   specify ensembl version (as given by
#'   \code{biomaRt::listEnsemblArchives()}). If no argument provided, will
#'   default to ensembl version 97.
#' @param attributes Vector of attributes to extract from BioMart.
#' @param filter Vector of filter to be used for BioMart query.
#'
#' @return Original dataframe together with extracted biomart query.
#' @export
#' 

biomart_df <- function(dataframe, columnToFilter, mart = 38, ensembl_version = NULL, attributes, filter){
  
  require(biomaRt)
  library(tidyverse)
  
  if(mart != 38 && mart != 37) stop("Mart must be 38 or 37...")
  
  if(mart == 38){
    
    if(is.null(ensembl_version)){
      
      ensembl_mart <- useMart(host = "jul2019.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
      
    } else{

      ensembl_mart <- useMart(host = ensembl_version, biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")      
      
    }
    
  }else if(mart == 37){
    
    ensembl_mart <-
      useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    
  }
  
  # Query genes as a vector
  genes <- dataframe %>% .[[columnToFilter]] %>% unique()
  print(str_c("Number of unique genes to search: ", length(genes)))
  
  # BioMart search
  biomart_query <- getBM(attributes = attributes, filters = filter, values = genes , mart = ensembl_mart)
  print(str_c("Number of matches found:", nrow(biomart_query)))
  
  # Create new data frame with positional information + remainder of the original dataframe
  # First requires creating join vector for the by argument in inner_join
  join_vector <- filter
  names(join_vector) <- columnToFilter
  Merged <- inner_join(dataframe, biomart_query, by = join_vector)
  
  return(Merged)
  
}

#' Vector-based biomart query.
#'
#' @param mart Specify genome build. 
#' @param ensembl_version character scalar. If using GRCh38 genome build, can
#'   specify ensembl version (as given by
#'   \code{biomaRt::listEnsemblArchives()}). If no argument provided, will
#'   default to ensembl version 97.
#' @param attributes Vector of attributes to extract from BioMart.
#' @param filter Vector of filter to be used for BioMart query.
#' @param values
#'
#' @return Original dataframe together with extracted biomart query.
#' @export
#'

query_biomart <- function(mart = 38, ensembl_version = NULL, attributes, filter, values){
  
  require(biomaRt)
  library(tidyverse)
  
  if(mart != 38 && mart != 37) stop("Mart must be 38 or 37...")
  
  if(mart == 38){
    
    if(is.null(ensembl_version)){
      
      ensembl_mart <- useMart(host = "jul2019.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
      
    } else{
      
      ensembl_mart <- useMart(host = ensembl_version, biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")      
      
    }
    
  }else if(mart == 37){
    
    ensembl_mart <-
      useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    
  }
  
  # BioMart search
  biomart_query <- getBM(attributes = attributes, filters = filter, values = values , mart = ensembl_mart)
  
  return(biomart_query)

}

