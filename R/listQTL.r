#' list of QTL database
#' 
#' @return result
#' @export
#' @examples 
#' listQTL()
listQTL<-function(){
i<-c('QTL_Btau_gff', 
'QTL_GG_gff', 'QTL_EquCab_gff', 'QTL_SS_gff', 'QTL_OAR_gff');
j<-c('UMD3.1', 
'GG_4.0', 'EquCab2.0', 'SS_10.2', 'OAR_3.1');
result <- data.frame(i,j);
colnames(result) <- c('QTL','version');
return(result);
}