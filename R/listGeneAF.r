#' list of gene attributes and filters.
#'
#' @param data_set, choose one of '1,2,3,4,5'.'1' = 'btaurus_gene_ensembl','2' = 'ggalluse_gene_ensembl','3' = 'ecaballus_gene_ensembl','4' = 'sscrofa_gene_ensembl','5' = 'oaries_gene_ensembl'.
#' @return result
#' @export
#' @importFrom biomaRt listMarts
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listFilters
#' @export
#' @examples
#' listGeneAF(1)
listGeneAF <- function(data_set) {
    for (pkg in c("biomaRt")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("the ", pkg, " package needed for this function to work. Please install it.", 
                sep = ""), call. = FALSE)
        }
    }
    genedataset <- switch(data_set, `1` = "btaurus_gene_ensembl", `2` = "ggallus_gene_ensembl", 
        `3` = "ecaballus_gene_ensembl", `4` = "sscrofa_gene_ensembl", `5` = "oaries_gene_ensembl")
    if (data_set == 2) {
        ensembl <- useEnsembl("ensembl", version = 85)
        message("The version of chicken QTL is 4.0,and the version of gene is 4.0!")
    } else {
        ensembl <- useEnsembl("ensembl")
    }
    martlist <- listMarts(ensembl)
    mart <- useMart(martlist[1, 1], dataset = genedataset)
    result <- listFilters(mart)
    return(result)
}
