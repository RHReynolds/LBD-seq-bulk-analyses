#' Principal components (cor)relation with experimental covariates
#'
#' Modfied from the pcaExplorer package
#' (https://github.com/federicomarini/pcaExplorer/blob/master/R/correlatePCs.R)
#'
#' Computes the significance of (cor)relations between PCA scores and the sample
#' experimental covariates, using Kruskal-Wallis test for categorial variables
#' and the \code{cor.test} based on Spearman's correlation for continuous
#' variables
#'
#' @param pcaobj A \code{prcomp} object
#' @param coldata A \code{data.frame} object containing the experimental
#'   covariates
#' @param pcs A numeric vector, containing the corresponding PC number
#'
#' @return A \code{list} object with computed statistics (either Kruskal-Wallis
#'   rank sum statistic or the estimated value of association (rho)) and p
#'   values for each covariate and for each principal component
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3,betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' pcaobj <- prcomp(t(assay(rlt)))
#' correlatePCs(pcaobj,colData(dds))
#'
#'
#' @export

correlatePCs <- function (pcaobj, coldata, pcs = 1:4){
  coldataTypes <- sapply(coldata, class)
  x <- pcaobj$x
  res_pval <- matrix(NA, nrow = length(pcs), ncol = ncol(coldata))
  colnames(res_pval) <- colnames(coldata)
  rownames(res_pval) <- paste0("PC_", pcs)
  
  res_stat <- matrix(NA, nrow = length(pcs), ncol = ncol(coldata))
  colnames(res_stat) <- colnames(coldata)
  rownames(res_stat) <- paste0("PC_", pcs)
  for (i in 1:ncol(coldata)) {
    for (j in pcs) {
      if (coldataTypes[i] %in% c("factor", "character")) {
        if (length(levels(coldata[, i])) > 1) {
          res_stat[j, i] <- kruskal.test(x[, j], coldata[,i])$statistic
          res_pval[j, i] <- kruskal.test(x[, j], coldata[,i])$p.value
        }
      }
      else {
        res_stat[j, i] <- cor.test(x[, j], coldata[, i], method = "spearman")$estimate
        res_pval[j, i] <- cor.test(x[, j], coldata[, i], method = "spearman")$p.value
        
      }
    }
  }
  
  res <- setNames(list(res_stat, res_pval), c("statistic", "pvalue"))
  return(res)
}