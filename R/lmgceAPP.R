#' \code{\link{lmgce}} Shiny application
#'
#' A Shiny application to execute \code{\link{lmgce}}
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @return \code{NULL}. This function is called for its side effect
#' (launching the app).
#'
#' @export
lmgceAPP <- function() {
  appDir <- system.file("shiny-examples", "lmgceAPP", package = "GCEstim")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `GCEstim`.",
         call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
