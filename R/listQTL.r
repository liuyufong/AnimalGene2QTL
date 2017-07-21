#' list of QTL database
#' 
#' @return result
#' @export
#' @examples 
#' listQTL()
listQTL<-function(){
i<-c('QTL_Btau_gff', 
'QTL_GG_gff', 'QTL_EquCab_gff', 'QTL_SS_gff', 'QTL_OAR_gff');
j<-c('QTL_Btau_4.6', 
'QTL_GG_4.0', 'QTL_EquCab2.0', 'QTL_SS_10.2', 'QTL_OAR_3.1');
result <- data.frame(i,j);
colnames(result) <- c('QTL','version');
return(result);
}