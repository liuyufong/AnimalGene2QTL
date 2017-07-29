#' Retrieve QTL data by SNP.
#'
#' @param qtl_attributes, Attributes you want to retrieve. A possible
#' list of attributes can be retrieved using the function listQTLAF().
#' @param snp_filters, Filters (one or more) that should be used in
#' the query. A possible list of filters can be retrieved using the
#' function listSNPFilters().
#' @param snp_values, Values of the filter, e.g. vector of IDs. If
#' multiple filters are specified then the argument should be a list
#' of vectors of which the position of each vector corresponds to the
#' position of the filters in the filters argument.
#' @param data_set, choose one of "1,2,3,4,5".'1' = 'btaurus_snp',
#' '2' = 'ggallus_snp','3' = 'ecaballus_snp',
#' '4' = 'sscrofa_snp','5' = 'oaries_snp'.
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
#' snp_filters <- c('snp_filter');
#' snp_values <- c('rs3136845');
#' getQTLbySNP(qtl_attributes=c('QTL_ID'),snp_filters,
#' snp_values,data_set=2);
getQTLbySNP <- function(qtl_attributes, snp_filters, snp_values, data_set){
if(missing(qtl_attributes)){
stop("Argument 'qtl_attributes' must be specified.")
}
if(missing(snp_filters)){
stop("Argument 'snp_filters' must be specified.")
}
if(missing(data_set)){
stop("Argument 'data_set' must be specified.")
}
if(data_set < 1 || data_set > 5){
stop("Argument 'data_set' must be selected in '1','2','3','4','5'.
'1' = 'btaurus_snp','2' = 'ggallus_snp','3' = 'ecaballus_snp',
'4' = 'sscrofa_snp','5' = 'oaries_snp' ")
}
if(is.list(qtl_attributes)){
qtl_attributes <- t(qtl_attributes);
}
if(is.list(snp_filters)){
snp_values = snp_filters;
srow <- NROW(snp_values);
snp_filters = names(snp_filters);
}
if(!is.list(snp_filters)){
srow <- NROW(snp_values);
}
if(!is.list(snp_values)){
srow <- length(snp_values);
}
sflength <- length(snp_filters);
qtlattlength <- length(qtl_attributes);
for(pkg in c("RSQLite","biomaRt","AnimalQTLDB")){
if(!requireNamespace(pkg,quietly = TRUE)){
stop(paste("the ",pkg," package needed for this function to work.
Please install it.",sep = ""),
call. = FALSE)
}
}
con <- dbConnect(SQLite(), system.file("extdata","animalqtldb.db",
package = "AnimalQTLDB"));
snpdataset <- switch(data_set,"1" = "btaurus_snp",
"2" = "ggallus_snp",
"3" = "ecaballus_snp",
"4" = "sscrofa_snp",
"5" = "oaries_snp");
qtldatabase <- switch(data_set,
"1" = "QTL_Btau_gff",
"2" = "QTL_GG_gff",
"3" = "QTL_EquCab_gff",
"4" = "QTL_SS_gff",
"5" = "QTL_OAR_gff");
a1 <- "select";
a2 <- "from ";
a3 <- " where Chromosome = '";
a4 <- "' and Chstart <= '";
a5 <- "' and Chend >= '";
a6 <- "'";
for(n in 1:qtlattlength){
if(n == 1){
att <- qtl_attributes[n];
}
if(n > 1){
att <- paste(att, ",");
att <- paste(att, qtl_attributes[n]);
}
}
result <- list();
if(data_set == 2){
ensembl<-useEnsembl("ensembl",version = 85);
message("The version of chicken QTL is 4.0,and the
version of gene is 4.0!");
host<-"grch37.ensembl.org";
}else{
ensembl<-useEnsembl("ensembl");
host<-"www.ensembl.org";
}
martlist<-listMarts(ensembl);
if(data_set == 2){
mart <- useMart(martlist[2,1], dataset = snpdataset, host=host);
}else{
mart <- useMart(martlist[3,1], dataset = snpdataset, host=host);
}
for(r in 1:srow){
if(is.list(snp_values)){
chromchr <- getBM(attributes=c('chr_name','chrom_start',
'chrom_end'),filters=snp_filters,values=snp_values[r,],
mart=mart);
}
if(!is.list(snp_values)){
chromchr <- getBM(attributes=c('chr_name','chrom_start',
'chrom_end'), filters=snp_filters, values=snp_values[r],
mart=mart);
}
chromrow <- NROW(chromchr);
for(i in 1:chromrow){
query <- paste(a1, att, a2, qtldatabase, a3, chromchr[i,1],
a4, chromchr[i,2], a5, chromchr[i,3],a6);
query <- gsub(pattern = "' ", replacement = "'", query);
query <- gsub(pattern = " '", replacement = "'", query);
single_QTL <- dbGetQuery(con, query);
NROW <- NROW(single_QTL);
if(NROW >= 1){
for(j in 1:NROW){
if(is.list(snp_values)){
Nlist <- data.frame(snp_values[r,1:sflength],
single_QTL[j,1:length(qtl_attributes)]);
}
if(!is.list(snp_values)){
Nlist <- data.frame(snp_values[r],
single_QTL[j,1:length(qtl_attributes)]);
}
if(names(Nlist)[1:sflength] != snp_filters){
colnames(Nlist) <- c(snp_filters, qtl_attributes);
}
result <- rbind(result, Nlist);
}
}
if(NROW < 1){
na <- list();
for(m in 1:length(qtl_attributes)){
na[m] <- "NA";
}
if(is.list(snp_values)){
Nlist <- data.frame(snp_values[r,1:sflength], na);
}
if(!is.list(snp_values)){
Nlist <- data.frame(snp_values[r], na);
}
if(names(Nlist)[1:sflength] != snp_filters){
colnames(Nlist) <- c(snp_filters ,qtl_attributes);
}
result <- rbind(result, Nlist);
}
}
}
dbDisconnect(con);
result <- unique(result);
colnames(result) <- c(snp_filters ,qtl_attributes);
return (result);
}
