#' Function for extracting duplicated values in a defined column.
#' 
#' @param df Dataframe.
#' @param col_w_duplicates Column name in quotation marks.
#' 
#' @return Dataframe with duplicated rows.
#'

keep_duplicates <- function(df, col_w_duplicates){
  df[duplicated(df[[col_w_duplicates]], fromLast = TRUE) | duplicated(df[[col_w_duplicates]], fromLast = FALSE),]
}
