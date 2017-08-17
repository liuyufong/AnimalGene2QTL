#' list of QTL database
#' 
#' @return result
#' @export
#' @importFrom RSQLite dbConnect
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite SQLite
#' @import knitr
#' @examples 
#' listQTL()
listQTL <- function() {
    con <- dbConnect(SQLite(), system.file("extdata", "animalqtldb.db", package = "AnimalQTLDB"))
    if (isS4(con)) {
        result <- dbGetQuery(con, "SELECT * from QTL_Version")
    }
    dbDisconnect(con)
    colnames(result) <- c("QTL", "version")
    return(result)
}
