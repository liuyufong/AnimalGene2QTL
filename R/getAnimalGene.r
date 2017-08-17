#' Retrieve gene data by QTL data.
#'
#' @param gene_attributes, Attributes you want to retrieve. A possible list of gene_attributes can be retrieved using the function listAttributes {biomaRt}.
#' @param qtl_filters, qtl_filters (one or more) that should be used in the query. A possible list of filters can be retrieved using the function listFilters.
#' @param qtl_values, Values of the qtl_filters, e.g. vector of IDs. If multiple qtl_filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the qtl_filters in the qtl_filters argument.
#' @param data_set, choose one of '1,2,3,4,5'.'1' = 'btaurus_gene_ensembl',2' = 'ggalluse_gene_ensembl',3' = 'ecaballus_gene_ensembl','4' = 'sscrofa_gene_ensembl',‘5' = 'oaries_gene_ensembl'.
#' @return result
#' @export
#' @importFrom RSQLite dbConnect
#' @importFrom RSQLite dbGetQuery
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite SQLite
#' @importFrom biomaRt listMarts
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom methods as
#' @import knitr
#' @import AnimalQTLDB
#' @examples
#' gene_attributes <- c('ensembl_gene_id');
#' qtl_filters <- c('QTL_ID');
#' qtl_values <- c('64577', '2199', '2354');
#' getAnimalGene(gene_attributes, qtl_filters, qtl_values, data_set = 2);
getAnimalGene <- function(gene_attributes, qtl_filters, qtl_values, data_set) {
    if (missing(gene_attributes)) {
        stop("Argument 'gene_attributes' must be specified.")
    }
    if (is.list(qtl_filters) && !missing(qtl_values)) {
        warning("Argument 'qtl_values' should not be used when argument 'qtl_filters' is a list and will be ignored.")
    }
    if (is.list(qtl_filters) && is.null(names(qtl_filters))) {
        stop("Argument 'qtl_filters' must be a named list when sent as a list.")
    }
    if (!is.list(qtl_filters) && qtl_filters != "" && missing(qtl_values)) {
        stop("Argument 'qtl_values' must be specified.")
    }
    if (length(qtl_filters) > 0 && length(qtl_values) == 0) {
        stop("qtl_values argument contains no data.")
    }
    if (data_set < 1 || data_set > 5) {
        stop("Argument 'data_set' must be a selected in '1','2','3','4','5'.'1'='btaurus_gene_ensembl','2'='ggalluse_gene_ensembl','3'='ecaballus_gene_ensembl','4'='sscrofa_gene_ensembl','5'='oaries_gene_ensembl' ")
    }
    if (is.list(gene_attributes)) {
        gene_attributes <- t(gene_attributes)
    }
    if (is.list(qtl_filters)) {
        qtl_values = qtl_filters
        qtlrow <- nrow(qtl_values)
        qtl_filters = names(qtl_filters)
    }
    if (!is.list(qtl_filters)) {
        qtlrow <- nrow(qtl_values)
    }
    if (!is.list(qtl_values)) {
        qtlrow <- length(qtl_values)
    }
    flength <- length(qtl_filters)
    attlength <- length(gene_attributes)
    for (pkg in c("RSQLite", "biomaRt", "AnimalQTLDB")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("the ", pkg, " package needed for this function to work. Please install it.",
                sep = ""), call. = FALSE)
        }
    }
    con <- dbConnect(SQLite(), system.file("extdata", "animalqtldb.db", package = "AnimalQTLDB"))
    dataset <- switch(data_set, `1` = "btaurus_gene_ensembl", `2` = "ggallus_gene_ensembl",
        `3` = "ecaballus_gene_ensembl", `4` = "sscrofa_gene_ensembl", `5` = "oaries_gene_ensembl")
    datatable <- switch(data_set, `1` = "QTL_Btau_gff", `2` = "QTL_GG_gff", `3` = "QTL_EquCab_gff",
        `4` = "QTL_SS_gff", `5` = "QTL_OAR_gff")
    a1 <- "select Chromosome,Chstart,Chend from "
    result <- data.frame()
    if (data_set == 2) {
        ensembl <- useEnsembl("ensembl", version = 85)
        message("The version of chicken QTL is 4.0,and the version of gene is 4.0!")
    } else {
        ensembl <- useEnsembl("ensembl")
    }
    martlist <- listMarts(ensembl)
    mart <- useMart(martlist[1, 1], dataset = dataset)
    for (i in seq_len(qtlrow)) {
        for (n in seq_len(flength)) {
            if (n == 1) {
                if (is.list(qtl_values)) {
                  a2 <- paste("where", qtl_filters[n], " = '", qtl_values[i, n], "'")
                } else if (!is.list(qtl_values)) {
                  a2 <- paste("where", qtl_filters[n], " = '", qtl_values[i], "'")
                }
            }
            if (n > 1) {
                a2 <- paste(a2, "and", qtl_filters[n], " = '", qtl_values[i, n], "'")
            }
        }
        query <- paste(a1, datatable, a2)
        query <- gsub(pattern = "' ", replacement = "'", query)
        query <- gsub(pattern = " '", replacement = "'", query)
        singleRowQuery <- dbGetQuery(con, query)
        geneNA <- character()
        if (nrow(singleRowQuery) >= 1) {
            singleRowgene <- getBM(attributes = gene_attributes, filters = c("chromosome_name",
                "start", "end"), values = list(singleRowQuery[, 1], singleRowQuery[, 2], singleRowQuery[,
                3]), mart = mart)
            geneRow <- nrow(singleRowgene)
            if (geneRow < 1) {
                geneNA[seq_len(attlength)] <- "NA"
                geneNA <- as(geneNA, "list")
                if (is.list(qtl_values)) {
                  qtlgene <- data.frame(geneNA, qtl_values[i, ])
                  colnames(qtlgene)[seq_len(attlength + flength)] <- c(gene_attributes, qtl_filters)
                  result <- rbind(result, qtlgene)
                } else if (!is.list(qtl_values)) {
                  qtlgene <- data.frame(geneNA, qtl_values[i])
                  colnames(qtlgene)[seq_len(attlength + flength)] <- c(gene_attributes, qtl_filters)
                  result <- rbind(result, qtlgene)
                }
                geneNA <- NULL
            }
            if (geneRow >= 1) {
                for (j in seq_len(geneRow)) {
                  if (is.list(qtl_values)) {
                    qtlgene <- data.frame(singleRowgene[j, seq_len(attlength)], qtl_values[i,
                      ])
                    colnames(qtlgene)[seq_len(attlength + flength)] <- c(gene_attributes, qtl_filters)
                    result <- rbind(result, qtlgene)
                  } else if (!is.list(qtl_values)) {
                    qtlgene <- data.frame(singleRowgene[j, seq_len(attlength)], qtl_values[i])
                    colnames(qtlgene)[seq_len(attlength + flength)] <- c(gene_attributes, qtl_filters)
                    result <- rbind(result, qtlgene)
                  }
                }
            }
        }
        if (nrow(singleRowQuery) < 1) {
            geneNA[seq_len(attlength)] <- "NO this 'value' in QTL database"
            geneNA <- as(geneNA, "list")
            if (is.list(qtl_values)) {
                qtlgene <- data.frame(geneNA, qtl_values[i, ])
                colnames(qtlgene)[seq_len(attlength + flength)] <- c(gene_attributes, qtl_filters)
                result <- rbind(result, qtlgene)
            } else if (!is.list(qtl_values)) {
                qtlgene <- data.frame(geneNA, qtl_values[i])
                colnames(qtlgene)[seq_len(attlength + flength)] <- c(gene_attributes, qtl_filters)
                result <- rbind(result, qtlgene)
            }
            geneNA <- NULL
        }
    }
    dbDisconnect(con)
    result <- unique(result)
    return(result)
}
