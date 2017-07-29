#' list of snp attributes and filters.
#'
#' @param data_set, choose one of "1,2,3,4,5".'1' = 'btaurus_snp',
#' '2' = 'ggallus_snp','3' = 'ecaballus_snp',
#' '4' = 'sscrofa_snp','5' = 'oaries_s
#' @return result
#' @export
#' @importFrom biomaRt listMarts
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listAttributes
#' @export
#' @examples
#' listSNPAttributes(1)
listSNPAttributes<-function(data_set){
for(pkg in c("biomaRt")){
if(!requireNamespace(pkg,quietly = TRUE)){
stop(paste("the ",pkg," package needed for this function to work.
Please install it.",sep = ""),
call. = FALSE)
}
}
snpdataset <- switch(data_set,
"1" = "btaurus_snp",
"2" = "ggallus_snp",
"3" = "ecaballus_snp",
"4" = "sscrofa_snp",
"5" = "oaries_snp");
if(data_set == 2){
ensembl<-useEnsembl("ensembl",version = 85);
message("The version of chicken QTL is 4.0,and the
version of SNP is 4.0!");
host<-"grch37.ensembl.org";
}else{
ensembl<-useEnsembl("ensembl");
host<-"www.ensembl.org";
}
martlist<-listMarts(ensembl);
mart <- useMart(martlist[2,1], dataset = snpdataset, host=host);
result <- listAttributes(mart);
return(result);
}
