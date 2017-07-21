#' list of snp attributes and filters.
#' 
#' @param data_set, choose one of "1,2,3,4,5".'1' = 'btaurus_snp',
#' '2' = 'ggallus_snp','3' = 'ecaballus_snp',
#' '4' = 'sscrofa_snp','5' = 'oaries_s
#' @return result
#' @export
#' @importFrom biomaRt useMart
#' @importFrom  biomaRt listFilters
#' @export
#' @examples 
#' listSNPFilters(1)
listSNPFilters<-function(data_set){
for(pkg in c("biomaRt")){
if(!requireNamespace(pkg,quietly = TRUE)){
stop(paste("the ",pkg," package needed for this function to work. 
Please install it.",sep = ""),call. = FALSE)}}
snpdataset <- switch (data_set,
"1" = "btaurus_snp",
"2" = "ggallus_snp",
"3" = "ecaballus_snp",
"4" = "sscrofa_snp",
"5" = "oaries_snp");
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = snpdataset);
result <- listFilters(ensembl);
return(result);
}