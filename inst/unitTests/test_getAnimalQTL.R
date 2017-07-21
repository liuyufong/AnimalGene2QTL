test_getAnimalQTL <- function() {
checkEquals(getAnimalQTL('QTL_ID', 'ensembl_gene_id', 
'ENSGALG00000003446', 2)[1,2], "24936")
checkTrue(is.na(getAnimalQTL('QTL_ID', 'ensembl_gene_id', 
'ENSGALG00000028410', 2)[1,2]))
}