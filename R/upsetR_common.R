# functions by docmanny from https://github.com/hms-dbmi/UpSetR/issues/85

# get_intersect_members() takes as arguments a dataframe that has been formatted as a binary 
#     table, such as movies from the UpSetR vignette; as well as a series of strings with the names of 
#     columns you wish to test membership for.
get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}



#modified version of the fromList function in the package that conserves the item names as row names
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}



#test
#listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), 
#                  two = c(1, 2, 4, 5,10), 
#                  three = c(1, 5, 6, 7, 8, 9,1, 12, 13))
#x <- fromList(listInput)
#get_intersect_members(x, 'one', 'three') # points that are in one and three
#get_intersect_members(x, 'one') # points that are only in one 


#upset(fromList(listInput), sets = names(listInput), keep.order = TRUE, nintersects = 25, order.by = "freq")
