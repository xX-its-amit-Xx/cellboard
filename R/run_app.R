#' Launch the cellboard dashboard
#'
#' Builds and runs the Shiny application. With no arguments it opens the app in
#' your browser; on Posit Connect the returned `shiny.appobj` is served directly.
#'
#' @param max_upload_mb Maximum upload size in megabytes (default 1024).
#' @param ... Passed to [shiny::shinyApp()].
#'
#' @return A Shiny app object (invisibly launched when interactive).
#' @export
#'
#' @examples
#' \dontrun{
#'   run_app()
#' }
run_app <- function(max_upload_mb = 1024, ...) {
  options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  shiny::shinyApp(ui = app_ui(), server = app_server, ...)
}
