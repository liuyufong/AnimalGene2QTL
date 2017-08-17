#' Retrieve gene data by QTL data.
#'
#' @param snp_attributes, Attributes you want to retrieve. A possible list of snp_attributes can be retrieved using the function listAttributes {biomaRt}.
#' @param qtl_filters, qtl_filters (one or more) that should be used in the query. A possible list of filters can be retrieved using the function listFilters.
#' @param qtl_values, Values of the qtl_filters, e.g. vector of IDs. If multiple qtl_filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the qtl_filters in the qtl_filters argument.
#' @param data_set, choose one of '1,2,3,4,5'.'1' = 'btaurus_gene_ensembl', '2' = 'ggalluse_gene_ensembl', '3' = 'ecaballus_gene_ensembl','4' = 'sscrofa_gene_ensembl', '5' = 'oaries_gene_ensembl'.
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
#' @import knitr
#' @import AnimalQTLDB
#' @examples
#' snp_attributes <- c('refsnp_id');
#' qtl_filters <- c('QTL_ID');
#' qtl_values <- c('64577', '2199', '2354');
#' getSNPbyQTL(snp_attributes, qtl_filters, qtl_values, data_set = 2);
getSNPbyQTL <- function(snp_attributes, qtl_filters, qtl_values, data_set) {
    if (missing(snp_attributes)) {
        stop("Argument 'snp_attributes' must be specified.")
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
        stop("Argument 'data_set' must be selected in '1','2','3','4','5'.'1' = 'btaurus_snp','2' = 'ggallus_snp','3' = 'ecaballus_snp','4' = 'sscrofa_snp','5' = 'oaries_snp' ")
    }
    if (is.list(snp_attributes)) {
        snp_attributes <- t(snp_attributes)
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
    qtlflength <- length(qtl_filters)
    snpattlength <- length(snp_attributes)
    for (pkg in c("RSQLite", "biomaRt", "AnimalQTLDB")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("the ", pkg, " package needed for this function to work. Please install it.", 
                sep = ""), call. = FALSE)
        }
    }
    con <- dbConnect(SQLite(), system.file("extdata", "animalqtldb.db", package = "AnimalQTLDB"))
    snpdataset <- switch(data_set, `1` = "btaurus_snp", `2` = "ggallus_snp", `3` = "ecaballus_snp", 
        `4` = "sscrofa_snp", `5` = "oaries_snp")
    qtldatatable <- switch(data_set, `1` = "QTL_Btau_gff", `2` = "QTL_GG_gff", `3` = "QTL_EquCab_gff", 
        `4` = "QTL_SS_gff", `5` = "QTL_OAR_gff")
    a1 <- paste("select", toString(qtl_filters), ",Chromosome,Chstart,Chend from ")
    a2 <- "where"
    result <- data.frame()
    if (data_set == 2) {
        ensembl <- useEnsembl("ensembl", version = 85)
        message("The version of chicken QTL is 4.0,and the version of SNP is 4.0!")
        host <- "grch37.ensembl.org"
    } else {
        ensembl <- useEnsembl("ensembl")
        host <- "www.ensembl.org"
    }
    martlist <- listMarts(ensembl)
    if (data_set == 2) {
        mart <- useMart(martlist[2, 1], dataset = snpdataset, host = host)
    } else {
        mart <- useMart(martlist[3, 1], dataset = snpdataset, host = host)
    }
    for (i in seq_len(qtlflength)) {
        if (i == 1) {
            a3 <- paste(qtl_filters[i], " in (")
        } else {
            a3 <- paste(a3, " and ", qtl_filters[i], " in (")
        }
        for (j in seq(qtlrow)) {
            if (j == qtlrow) {
                if (is.list(qtl_values)) {
                  a3 <- paste(a3, "'", qtl_values[j, i], "')")
                }
                if (!is.list(qtl_values)) {
                  a3 <- paste(a3, "'", qtl_values[j], "')")
                }
            } else {
                if (is.list(qtl_values)) {
                  a3 <- paste(a3, "'", qtl_values[j, i], "',")
                }
                if (!is.list(qtl_values)) {
                  a3 <- paste(a3, "'", qtl_values[j], "',")
                }
            }
        }
    }
    query <- paste(a1, qtldatatable, a2, a3)
    query <- gsub(pattern = "' ", replacement = "'", query)
    query <- gsub(pattern = " '", replacement = "'", query)
    qtlquery <- dbGetQuery(con, query)
    qtlrow <- NROW(qtlquery)
    geneNA <- character()
    if (qtlrow >= 1) {
        snp_all <- getBM(attributes = c(snp_attributes, "chr_name", "chrom_start", "chrom_end"), 
            filters = c("chr_name", "start", "end"), values = list(qtlquery[, qtlflength + 1], 
                qtlquery[, qtlflength + 2], qtlquery[, qtlflength + 3]), mart = mart)
        snprow <- NROW(snp_all)
        snp_result <- data.frame()
        for (k in seq_len(qtlrow)) {
            sg_snp <- snp_all[which(snp_all[snpattlength + 1] == qtlquery[k, qtlflength + 1]), 
                ]
            sg_snp <- sg_snp[which(sg_snp[snpattlength + 2] >= qtlquery[k, qtlflength + 2]), 
                ]
            sg_snp <- sg_snp[which(sg_snp[snpattlength + 3] <= qtlquery[k, qtlflength + 3]), 
                ]
            sg_snp <- sg_snp[-(snpattlength + 1):-(snpattlength + 3)]
            sg_snp[(snpattlength + 1):(snpattlength + 2)] <- qtlquery[k, (qtlflength + 2):(qtlflength + 
                3)]
            snp_result <- rbind(snp_result, sg_snp)
            rm(sg_snp)
            gc()
        }
        qtlquery <- qtlquery[-(qtlflength + 1)]
        result <- merge(qtlquery, snp_result, all = TRUE)
        rm(snp_result)
        rm(qtlquery)
        result <- result[-1:-2]
        gc()
    } else {
        result <- "There is NULL QTL database by searching for input QTL data. Please check your data."
    }
    dbDisconnect(con)
    result <- unique(result)
    gc()
    return(result)
}
