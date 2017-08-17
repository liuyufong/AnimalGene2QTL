#'Retrieve QTL data by gene.
#'
#' @param qtl_attributes, Attributes you want to retrieve. A possible list of attributes can be retrieved using the function listQTLAF().
#' @param gene_filters, Filters (one or more) that should be used in the query. A possible list of filters can be retrieved using the function listGeneAF().
#' @param gene_values, Values of the filter, e.g. vector of IDs. If multiple filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the filters in the filters argument.
#' @param data_set, choose one of '1,2,3,4,5'.'1' = 'btaurus_gene_ensembl', '2' = 'ggalluse_gene_ensembl','3' = 'ecaballus_gene_ensembl', '4' = 'sscrofa_gene_ensembl','5' = 'oaries_gene_ensembl'.
#' @param snp, Default 'FALSE'.Detemine whether retrieve SNP data by TRUE or FALSE.
#' @param snp_attributes, SNP attributes you want to retrieve.A possible list of attributes can be retrieved using the function listSNPAttributes().
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
#' gene_filters <- c('ensembl_gene_id');
#' gene_values <- c('ENSBTAG00000009851');
#' getAnimalQTL(qtl_attributes=c('QTL_ID'), gene_filters, gene_values, data_set=1, snp=TRUE, 'refsnp_id');
getAnimalQTL <- function(qtl_attributes, gene_filters, gene_values, data_set, snp = "", snp_attributes = "") {
    if (missing(qtl_attributes)) {
        stop("Argument 'qtl_attributes' must be specified.")
    }
    if (missing(gene_filters)) {
        stop("Argument 'gene_filters' must be specified.")
    }
    if (missing(data_set)) {
        stop("Argument 'data_set' must be specified.")
    }
    if (snp == TRUE) {
        if (missing(snp_attributes)) {
            stop("Argument 'snp_attributes' must be specified when snp is TRUE.")
        }
        if (is.list(snp_attributes)) {
            snp_attributes <- t(snp_attributes)
        }
        if (NROW(gene_values) > 11) {
            stop("gene value must be less than 10 rows if you retrieve SNP because the SNP data is huge.")
        }
    }
    if (snp == FALSE) {
        if (!missing(snp_attributes)) {
            stop("Argument 'snp_attributes' must not be specified when snp is FALSE.")
        }
    }
    if (data_set < 1 || data_set > 5) {
        stop("Argument 'data_set' must be selected in '1','2','3','4','5'.
         '1' = 'btaurus_gene_ensembl',
         '2' = 'ggalluse_gene_ensembl',
         '3' = 'ecaballus_gene_ensembl',
         '4' = 'sscrofa_gene_ensembl',
         '5' = 'oaries_gene_ensembl' ")
    }
    if (is.list(qtl_attributes)) {
        qtl_attributes <- t(qtl_attributes)
    }
    if (is.list(gene_filters)) {
        gene_values <- gene_filters
        grow <- nrow(gene_values)
        gene_filters <- names(gene_filters)
    }
    if (!is.list(gene_filters)) {
        grow <- nrow(gene_values)
    }
    if (!is.list(gene_values)) {
        grow <- length(gene_values)
    }
    gflength <- length(gene_filters)
    qtlattlength <- length(qtl_attributes)
    snpattlength <- length(snp_attributes)
    for (pkg in c("RSQLite", "biomaRt", "AnimalQTLDB")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("the ", pkg, " package needed for this function to work.Please install it.", 
                sep = ""), call. = FALSE)
        }
    }
    con <- dbConnect(SQLite(), system.file("extdata", "animalqtldb.db", package = "AnimalQTLDB"))
    genedataset <- switch(data_set, `1` = "btaurus_gene_ensembl", `2` = "ggallus_gene_ensembl", 
        `3` = "ecaballus_gene_ensembl", `4` = "sscrofa_gene_ensembl", `5` = "oaries_gene_ensembl")
    qtldatabase <- switch(data_set, `1` = "QTL_Btau_gff", `2` = "QTL_GG_gff", `3` = "QTL_EquCab_gff", 
        `4` = "QTL_SS_gff", `5` = "QTL_OAR_gff")
    snpdataset <- switch(data_set, `1` = "btaurus_snp", `2` = "ggallus_snp", `3` = "ecaballus_snp", 
        `4` = "sscrofa_snp", `5` = "oaries_snp")
    a1 <- "select Chstart , Chend ,"
    a2 <- "from "
    a3 <- " where Chromosome = '"
    a4 <- "' and ((Chstart >= '"
    a5 <- "' and Chend <= '"
    a6 <- "')"
    a7 <- "or (Chstart <= '"
    a8 <- "' and Chend >= '"
    a9 <- "'))"
    b1 <- "select Chstart , Chend ,"
    b2 <- "from "
    b3 <- " where Chromosome = '"
    b4 <- "' and ((Chstart < '"
    b5 <- "' and Chend < '"
    b6 <- "' and Chend > '"
    b7 <- "')"
    b8 <- "or (Chstart > '"
    b9 <- "' and Chend > '"
    b10 <- "' and Chstart < '"
    b11 <- "'))"
    for (n in seq_len(qtlattlength)) {
        if (n == 1) {
            att <- qtl_attributes[n]
        }
        if (n > 1) {
            att <- paste(att, ",")
            att <- paste(att, qtl_attributes[n])
        }
    }
    if (data_set == 2) {
        ensembl <- useEnsembl("ensembl", version = 85)
        message("The version of chicken QTL is 4.0,and the version of gene is 4.0!")
        host <- "grch37.ensembl.org"
    } else {
        ensembl <- useEnsembl("ensembl")
        host <- "www.ensembl.org"
    }
    martlist <- listMarts(ensembl)
    mart <- useMart(martlist[1, 1], dataset = genedataset)
    chromchr <- getBM(attributes = c(gene_filters, "chromosome_name", "start_position", "end_position"), 
        filters = gene_filters, values = gene_values, mart = mart)
    chromrow <- NROW(chromchr)
    result <- data.frame()
    if (chromrow < 1) {
        if (!is.list(gene_values)) {
            gene_values <- t(t(gene_values))
            gene_values <- as.data.frame(gene_values)
        }
        NULL_QTL <- gene_values
        NULL_QTL[(NCOL(NULL_QTL) + 1):(NCOL(NULL_QTL) + qtlattlength + 1)] <- NA
        NULL_QTL[(qtlattlength + 1):(qtlattlength + gflength)] <- NULL_QTL[seq_len(gflength)]
        NULL_QTL[seq_len(qtlattlength)] <- NA
        colnames(NULL_QTL) <- c(qtl_attributes, gene_filters, "length in the region (bp)")
        result <- NULL_QTL
        if (snp == TRUE) {
            result[(NCOL(result) + 1):(NCOL(result) + snpattlength + 1)] <- NA
            colnames(result) <- c(qtl_attributes, gene_filters, "length in the region (bp)", 
                snp_attributes, "snp score(%)")
        }
    }
    if (chromrow >= 1) {
        if (snp == TRUE) {
            if (data_set == 2) {
                snp_mart <- useMart(martlist[2, 1], dataset = snpdataset, host = host)
            } else {
                snp_mart <- useMart(martlist[3, 1], dataset = snpdataset, host = host)
            }
            single_SNP <- getBM(attributes = c(snp_attributes, "chr_name", "chrom_start", "chrom_end"), 
                filters = c("chr_name", "start", "end"), values = list(chromchr[, gflength + 
                  1], chromchr[, gflength + 2], chromchr[, gflength + 3]), mart = snp_mart)
            snprow <- nrow(single_SNP)
        }
        for (i in seq_len(chromrow)) {
            query1 <- paste(a1, att, a2, qtldatabase, a3, chromchr[i, gflength + 1], a4, chromchr[i, 
                gflength + 2], a5, chromchr[i, gflength + 3], a6, a7, chromchr[i, gflength + 
                2], a8, chromchr[i, gflength + 3], a9)
            query2 <- paste(b1, att, b2, qtldatabase, b3, chromchr[i, gflength + 1], b4, chromchr[i, 
                gflength + 2], b5, chromchr[i, gflength + 3], b6, chromchr[i, gflength + 2], 
                b7, b8, chromchr[i, gflength + 2], b9, chromchr[i, gflength + 3], b10, chromchr[i, 
                  gflength + 3], b11)
            query1 <- gsub(pattern = "' ", replacement = "'", query1)
            query1 <- gsub(pattern = " '", replacement = "'", query1)
            query2 <- gsub(pattern = "' ", replacement = "'", query2)
            query2 <- gsub(pattern = " '", replacement = "'", query2)
            single_QTL1 <- dbGetQuery(con, query1)
            single_QTL2 <- dbGetQuery(con, query2)
            gc()
            Nrow <- nrow(single_QTL1)
            NSErow <- nrow(single_QTL2)
            snpresult <- data.frame()
            snqresult <- data.frame()
            if (Nrow >= 1) {
                single_QTL1[(NCOL(single_QTL1) + 1):(NCOL(single_QTL1) + gflength)] <- chromchr[i, 
                  seq_len(gflength)]
                single_QTL1[NCOL(single_QTL1) + 1] <- "ALL"
                colnames(single_QTL1) <- c("Chstart", "Chend", qtl_attributes, gene_filters, 
                  "length in the region (bp)")
                if (snp == TRUE) {
                  if (snprow >= 1) {
                    sg_SNP <- single_SNP[which(single_SNP[snpattlength + 1] == chromchr[i, gflength + 
                      1]), ]
                    for (k in seq_len(Nrow)) {
                      sg_SNP1 <- sg_SNP[which(sg_SNP[snpattlength + 2] >= single_QTL1[k, 1]), 
                        ]
                      sg_SNP1 <- sg_SNP1[which(sg_SNP1[snpattlength + 3] <= single_QTL1[k, 2]), 
                        ]
                      sg_SNP2 <- sg_SNP[which(sg_SNP[snpattlength + 2] < single_QTL1[k, 1]), 
                        ]
                      sg_SNP3 <- sg_SNP[which(sg_SNP[snpattlength + 3] > single_QTL1[k, 2]), 
                        ]
                      sg_SNP1 <- sg_SNP1[seq_len(snpattlength)]
                      sg_SNP2 <- sg_SNP2[seq_len(snpattlength)]
                      sg_SNP3 <- sg_SNP3[seq_len(snpattlength)]
                      sg_SNP1[(NCOL(sg_SNP1) + 1):(NCOL(sg_SNP1) + 2)] <- single_QTL1[k, seq_len(2)]
                      sg_SNP2[(NCOL(sg_SNP2) + 1):(NCOL(sg_SNP2) + 2)] <- single_QTL1[k, seq_len(2)]
                      sg_SNP3[(NCOL(sg_SNP3) + 1):(NCOL(sg_SNP3) + 2)] <- single_QTL1[k, seq_len(2)]
                      colnames(sg_SNP1) <- c(snp_attributes, "Chstart", "Chend")
                      colnames(sg_SNP2) <- c(snp_attributes, "Chstart", "Chend")
                      colnames(sg_SNP3) <- c(snp_attributes, "Chstart", "Chend")
                      snpresult <- rbind(snpresult, sg_SNP1)
                      snqresult <- rbind(snqresult, sg_SNP2)
                      snqresult <- rbind(snqresult, sg_SNP3)
                      rm(sg_SNP1)
                      rm(sg_SNP2)
                      rm(sg_SNP3)
                      gc()
                    }
                    snpresult <- as.data.frame(lapply(snpresult[1:NCOL(snpresult)], unlist)) 
                    single_QTL1 <- as.data.frame(lapply(single_QTL1[1:NCOL(single_QTL1)], unlist))
                    single_QTL1 <- merge(single_QTL1, snpresult, all = TRUE)
                    rm(snpresult)
                    gc()
                    if (NROW(snqresult) > 0) {
                      snqresult[(NCOL(snqresult) + 1):(NCOL(snqresult) + qtlattlength + gflength + 
                        1)] <- "NA"
                      snqresult[(qtlattlength + gflength + 2):NCOL(snqresult)] <- snqresult[seq_len(snpattlength + 
                        2)]
                      snqresult[seq_len(qtlattlength + gflength + 1)] <- NA
                      colnames(snqresult) <- c(qtl_attributes, gene_filters, "length in the region (bp)", 
                        snp_attributes, "Chstart", "Chend")
                      snqresult <- as.data.frame(lapply(snqresult[seq_len(NCOL(snqresult))], 
                        unlist))
                      single_QTL1 <- merge(single_QTL1, snqresult, all = TRUE)
                      rm(snqresult)
                      gc()
                      for (l in (3 + qtlattlength):(3 + qtlattlength)) {
                        single_QTL1[is.na(single_QTL1[l]), l] <- single_QTL1[1, l]
                      }
                    }
                    single_QTL1[NCOL(single_QTL1) + 1] <- NA
                    colnames(single_QTL1) <- c("Chstart", "Chend", qtl_attributes, gene_filters, 
                      "length in the region (bp)", snp_attributes, "snp score(%)")
                    rm(sg_SNP)
                    gc()
                  }
                }
                result <- rbind(result, single_QTL1[, 3:NCOL(single_QTL1)])
                rm(single_QTL1)
                gc()
            }
            if (NSErow >= 1) {
                if (snp == TRUE) {
                  snpresult <- data.frame()
                  snqresult <- data.frame()
                  sg_SNP <- single_SNP[which(single_SNP[snpattlength + 1] == chromchr[i, gflength + 
                    1]), ]
                }
                for (k in seq_len(NSErow)) {
                  if (single_QTL2[k, 2] > chromchr[i, gflength + 3]) {
                    singlebp <- chromchr[i, gflength + 3] - single_QTL2[k, 1]
                    singlebp <- as.character(singlebp)
                  }
                  if (single_QTL2[k, 2] < chromchr[i, gflength + 3]) {
                    singlebp <- single_QTL2[k, 2] - chromchr[i, gflength + 2]
                    singlebp <- as.character(singlebp)
                  }
                  if (k == 1) {
                    bp <- singlebp
                  }
                  if (k > 1) {
                    bp <- rbind(bp, singlebp)
                  }
                  if (snp == TRUE) {
                    if (snprow >= 1) {
                      sg_SNP1 <- sg_SNP[which(sg_SNP[snpattlength + 2] >= single_QTL2[k, 1]), 
                        ]
                      sg_SNP1 <- sg_SNP1[which(sg_SNP1[snpattlength + 3] <= single_QTL2[k, 2]), 
                        ]
                      sg_SNP2 <- sg_SNP[which(sg_SNP[snpattlength + 2] < single_QTL2[k, 1]), 
                        ]
                      sg_SNP3 <- sg_SNP[which(sg_SNP[snpattlength + 3] > single_QTL2[k, 2]), 
                        ]
                      sg_SNP1 <- sg_SNP1[seq_len(snpattlength)]
                      sg_SNP2 <- sg_SNP2[seq_len(snpattlength)]
                      sg_SNP3 <- sg_SNP3[seq_len(snpattlength)]
                      sg_SNP1[(NCOL(sg_SNP1) + 1):(NCOL(sg_SNP1) + 2)] <- single_QTL2[k, seq_len(2)]
                      sg_SNP2[(NCOL(sg_SNP2) + 1):(NCOL(sg_SNP2) + 2)] <- single_QTL2[k, seq_len(2)]
                      sg_SNP3[(NCOL(sg_SNP3) + 1):(NCOL(sg_SNP3) + 2)] <- single_QTL2[k, seq_len(2)]
                      colnames(sg_SNP1) <- c(snp_attributes, "Chstart", "Chend")
                      colnames(sg_SNP2) <- c(snp_attributes, "Chstart", "Chend")
                      colnames(sg_SNP3) <- c(snp_attributes, "Chstart", "Chend")
                      snpresult <- rbind(snpresult, sg_SNP1)
                      snqresult <- rbind(snqresult, sg_SNP2)
                      snqresult <- rbind(snqresult, sg_SNP3)
                      rm(sg_SNP1)
                      rm(sg_SNP2)
                      rm(sg_SNP3)
                      gc()
                    }
                  }
                }
                single_QTL2[(NCOL(single_QTL2) + 1):(NCOL(single_QTL2) + gflength)] <- chromchr[i, 
                  seq_len(gflength)]
                single_QTL2[NCOL(single_QTL2) + 1] <- bp
                colnames(single_QTL2) <- c("Chstart", "Chend", qtl_attributes, gene_filters, 
                  "length in the region (bp)")
                if (snp == TRUE) {
                  snpresult <- as.data.frame(lapply(snpresult[seq_len(NCOL(snpresult))], unlist))
                  single_QTL2 <- as.data.frame(lapply(single_QTL2[1:NCOL(single_QTL2)], unlist))
                  single_QTL2 <- merge(single_QTL2, snpresult, all = TRUE)
                  rm(snpresult)
                  gc()
                  if (NROW(snqresult) > 0) {
                    snqresult[(NCOL(snqresult) + 1):(NCOL(snqresult) + qtlattlength + gflength + 
                      1)] <- "NA"
                    snqresult[(qtlattlength + gflength + 2):NCOL(snqresult)] <- snqresult[seq_len(snpattlength + 
                      2)]
                    snqresult[seq_len(qtlattlength + gflength + 1)] <- NA
                    colnames(snqresult) <- c(qtl_attributes, gene_filters, "length in the region (bp)", 
                      snp_attributes, "Chstart", "Chend")
                    snqresult <- as.data.frame(lapply(snqresult[seq_len(NCOL(snqresult))], unlist))
                    single_QTL2 <- merge(single_QTL2, snqresult, all = TRUE)
                    rm(snqresult)
                    gc()
                    for (l in (3 + qtlattlength):(3 + qtlattlength)) {
                      single_QTL2[is.na(single_QTL2[l]), l] <- single_QTL2[1, l]
                    }
                  }
                  single_QTL2[NCOL(single_QTL2) + 1] <- NA
                  colnames(single_QTL2) <- c("Chstart", "Chend", qtl_attributes, gene_filters, 
                    "length in the region (bp)", snp_attributes, "snp score(%)")
                }
                result <- rbind(result, single_QTL2[3:NCOL(single_QTL2)])
                rm(single_QTL2)
                if (snp == TRUE) {
                  rm(sg_SNP)
                }
                gc()
            }
            if (Nrow < 1 && NSErow < 1) {
                single_QTL1[1, 3:NCOL(single_QTL1)] <- NA
                single_QTL1[(NCOL(single_QTL1) + 1):(NCOL(single_QTL1) + gflength)] <- chromchr[i, 
                  seq_len(gflength)]
                single_QTL1[NCOL(single_QTL1) + 1] <- NA
                single_QTL1[1, 1] <- chromchr[i, gflength + 2]
                single_QTL1[1, 2] <- chromchr[i, gflength + 3]
                colnames(single_QTL1) <- c("Chstart", "Chend", qtl_attributes, gene_filters, 
                  "length in the region (bp)")
                if (snp == TRUE) {
                  if (snprow >= 1) {
                    sg_SNP <- single_SNP[which(single_SNP[snpattlength + 1] == chromchr[i, gflength + 
                      1]), ]
                    sg_SNP <- sg_SNP[which(sg_SNP[snpattlength + 2] >= chromchr[i, gflength + 
                      2]), ]
                    sg_SNP <- sg_SNP[which(sg_SNP[snpattlength + 3] <= chromchr[i, gflength + 
                      3]), ]
                    sg_SNP <- sg_SNP[-(snpattlength + 1)]
                    sg_SNP[snpattlength + 1] <- chromchr[i, gflength + 2]
                    sg_SNP[snpattlength + 2] <- chromchr[i, gflength + 3]
                  }
                  sg_SNP <- as.data.frame(lapply(sg_SNP[1:NCOL(sg_SNP)], unlist)) 
                  single_QTL1 <- as.data.frame(lapply(single_QTL1[1:NCOL(single_QTL1)], unlist))
                  single_QTL1 <- merge(single_QTL1, sg_SNP, all = TRUE)
                  rm(sg_SNP)
                  single_QTL1[NCOL(single_QTL1) + 1] <- NA
                  colnames(single_QTL1) <- c("Chstart", "Chend", qtl_attributes, gene_filters, 
                    "length in the region (bp)", snp_attributes, "snp score(%)")
                }
                result <- rbind(result, single_QTL1[3:NCOL(single_QTL1)])
                rm(single_QTL1)
                gc()
            }
        }
    }
    dbDisconnect(con)
    if (!is.list(gene_values)) {
        gene_values <- t(t(gene_values))
        gene_values <- as.data.frame(gene_values)
    }
    colnames(gene_values) <- gene_filters
    if (snp == TRUE) {
        result$`snp score(%)` <- as.character(result$`snp score(%)`)
        if (chromrow > 0) {
            result[is.na(result[1]), NCOL(result)] <- 0.667
            result[!is.na(result[1]), NCOL(result)] <- 1
            result[is.na(result[NCOL(result)]), NCOL(result)] <- 0.333
        } else {
            result[is.na(result[NCOL(result)]), NCOL(result)] <- 0.333
        }
    }
    result <- merge(result, gene_values, all = TRUE)
    if (snp == TRUE) {
        result[is.na(result[NCOL(result)]), NCOL(result)] <- 0.333
    }
    result <- unique(result)
    return(result)
    gc()
}
