#' list of QTL attributes and filters.
#' 
#' @return result
#' @export
#' @importFrom RSQLite dbConnect
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite dbGetQuery
#' @importFrom RSQLite SQLite
#' @import knitr
#' @examples 
#' listQTLAF()
listQTLAF <- function() {
    con <- dbConnect(SQLite(), system.file("extdata", "animalqtldb.db", package = "AnimalQTLDB"))
    if (isS4(con)) {
        result <- dbGetQuery(con, "SELECT * from QTL_Field_Name")
    }
    dbDisconnect(con)
    colnames(result) <- c("name", "description")
    return(result)
}
