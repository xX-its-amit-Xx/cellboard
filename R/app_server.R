# The dashboard server. Flow:
#   load (example / upload) -> rv$raw
#   "Compute QC" button     -> compute_qc() -> rv$qc
#   threshold sliders        -> thresholds() (reactive)
#   apply_qc_filter()        -> filt() (reactive) -> plots, tables, downloads
# Everything downstream of filt() updates live as the sliders move.

#' cellboard Shiny server
#'
#' @param input,output,session Standard Shiny server arguments.
#' @return Used for its side effects.
#' @export
app_server <- function(input, output, session) {

  rv <- shiny::reactiveValues(raw = NULL, qc = NULL, source = NULL)

  # ---- Load data ----------------------------------------------------------
  shiny::observeEvent(input$load_example, {
    sce <- tryCatch(example_sce(), error = function(e) {
      shiny::showNotification(paste("Example unavailable:", conditionMessage(e)),
                              type = "error", duration = NULL)
      NULL
    })
    if (!is.null(sce)) {
      rv$raw <- sce; rv$qc <- NULL
      rv$source <- "Example PBMC (simulated)"
      shiny::showNotification(
        sprintf("Loaded example: %d genes x %d cells.", nrow(sce), ncol(sce)),
        type = "message")
      bslib::nav_select("nav", "Data")
    }
  })

  shiny::observeEvent(input$file, {
    fdf <- input$file
    shiny::req(fdf)
    sce <- tryCatch({
      if (nrow(fdf) == 1L) {
        ext <- tolower(tools::file_ext(fdf$name[1]))
        if (ext %in% c("rds", "rdata")) {
          read_input(fdf$datapath[1], type = "sce")
        } else {
          read_10x(fdf$datapath[1])         # assume a 10x .h5
        }
      } else {
        # Several files selected -> reconstruct a 10x directory by name.
        tmp <- tempfile("upload10x_")
        dir.create(tmp)
        file.copy(fdf$datapath, file.path(tmp, fdf$name), overwrite = TRUE)
        read_10x(tmp)
      }
    }, error = function(e) {
      shiny::showNotification(paste("Could not read input:", conditionMessage(e)),
                              type = "error", duration = NULL)
      NULL
    })
    if (!is.null(sce)) {
      rv$raw <- sce; rv$qc <- NULL
      rv$source <- if (nrow(fdf) == 1L) fdf$name[1] else
        sprintf("10x upload (%d files)", nrow(fdf))
      shiny::showNotification(
        sprintf("Loaded %d genes x %d cells.", nrow(sce), ncol(sce)),
        type = "message")
      bslib::nav_select("nav", "Data")
    }
  })

  # ---- Compute QC ---------------------------------------------------------
  shiny::observeEvent(input$compute, {
    if (is.null(rv$raw)) {
      shiny::showNotification("Load a dataset first.", type = "warning")
      return(invisible())
    }
    shiny::withProgress(message = "Computing QC metrics", value = 0.1, {
      shiny::incProgress(0.2, detail = "mito %, genes, UMIs")
      mp <- input$mito_pattern %||% ""
      if (!nzchar(mp)) mp <- "^[Mm][Tt]-"
      qc <- tryCatch({
        compute_qc(rv$raw, mito_pattern = mp,
                   run_doublets = isTRUE(input$run_doublets))
      }, error = function(e) {
        shiny::showNotification(paste("QC failed:", conditionMessage(e)),
                                type = "error", duration = NULL)
        NULL
      })
      shiny::incProgress(0.6, detail = "doublets")
      rv$qc <- qc
    })
    if (!is.null(rv$qc)) {
      shiny::showNotification("QC metrics computed.", type = "message")
      bslib::nav_select("nav", "QC overview")
    }
  })

  # ---- Threshold controls (data-driven) -----------------------------------
  output$threshold_ui <- shiny::renderUI({
    qc <- rv$qc
    if (is.null(qc)) {
      return(shiny::helpText("Compute QC to enable threshold sliders."))
    }
    cd <- coldata_df(qc)
    th <- default_thresholds(qc)
    max_g <- max(cd[[.cb_cols[["nGene"]]]], na.rm = TRUE)
    max_u <- ceiling(max(cd[[.cb_cols[["nUMI"]]]], na.rm = TRUE))
    has_db <- isTRUE(S4Vectors::metadata(qc)$cellboard$qc$ran_doublets) &&
      !all(is.na(cd[[.cb_cols[["doublet_score"]]]]))

    shiny::tagList(
      shiny::sliderInput("genes_range", "Genes per cell",
        min = 0, max = max_g,
        value = c(min(th$min_genes %||% 0, max_g), max_g),
        step = max(1, round(max_g / 100))),
      shiny::sliderInput("umi_range", "UMIs per cell",
        min = 0, max = max_u,
        value = c(min(th$min_umi %||% 0, max_u), max_u),
        step = max(1, round(max_u / 100))),
      shiny::sliderInput("max_mito", "Max mitochondrial %",
        min = 0, max = 100, value = th$max_mito %||% 100, step = 0.5),
      if (has_db) shiny::tagList(
        shiny::checkboxInput("remove_doublet_class",
          "Discard scDblFinder doublets", TRUE),
        shiny::sliderInput("max_doublet", "Max doublet score",
          min = 0, max = 1, value = 1, step = 0.01)
      )
    )
  })

  thresholds <- shiny::reactive({
    shiny::req(rv$qc)
    cd <- coldata_df(rv$qc)
    max_g <- max(cd[[.cb_cols[["nGene"]]]], na.rm = TRUE)
    max_u <- max(cd[[.cb_cols[["nUMI"]]]], na.rm = TRUE)
    gr <- input$genes_range %||% c(0, max_g)
    ur <- input$umi_range %||% c(0, max_u)
    qc_thresholds(
      min_genes   = gr[1],
      max_genes   = if (gr[2] < max_g) gr[2] else NULL,
      min_umi     = ur[1],
      max_umi     = if (ur[2] < max_u) ur[2] else NULL,
      max_mito    = input$max_mito %||% 100,
      max_doublet = if (!is.null(input$max_doublet) && input$max_doublet < 1)
        input$max_doublet else NULL,
      remove_doublet_class = isTRUE(input$remove_doublet_class)
    )
  })

  filt <- shiny::reactive({
    shiny::req(rv$qc)
    apply_qc_filter(rv$qc, thresholds())
  })

  # ---- Data tab -----------------------------------------------------------
  output$n_genes <- shiny::renderText(if (is.null(rv$raw)) "-" else fmt_int(nrow(rv$raw)))
  output$n_cells <- shiny::renderText(if (is.null(rv$raw)) "-" else fmt_int(ncol(rv$raw)))
  output$data_source <- shiny::renderText(rv$source %||% "-")

  output$data_status <- shiny::renderUI({
    if (is.null(rv$raw)) {
      return(bslib::card_body(shiny::p(
        "No data loaded. Use ", shiny::strong("Load example PBMC"),
        " in the sidebar, or upload your own, then click ",
        shiny::strong("Compute QC metrics"), ".")))
    }
    bslib::card_body(
      shiny::p(shiny::strong("Source: "), rv$source),
      shiny::p(shiny::strong("Assays: "),
               paste(SummarizedExperiment::assayNames(rv$raw), collapse = ", ")),
      shiny::p(if (is.null(rv$qc))
        "QC not yet computed - click \"Compute QC metrics\"." else
        "QC computed. See the QC overview and Filtering tabs.")
    )
  })

  output$rowdata_preview <- DT::renderDT({
    shiny::req(rv$raw)
    rd <- as.data.frame(SummarizedExperiment::rowData(rv$raw))
    if (!ncol(rd)) rd <- data.frame(feature = utils::head(rownames(rv$raw), 50))
    DT::datatable(utils::head(rd, 50), rownames = TRUE,
                  options = list(pageLength = 5, scrollX = TRUE, dom = "tip"))
  })

  # ---- QC overview value boxes -------------------------------------------
  med <- function(col, fmt) shiny::renderText({
    if (is.null(rv$qc)) return("-")
    v <- stats::median(coldata_df(rv$qc)[[col]], na.rm = TRUE)
    fmt(v)
  })
  output$med_umi   <- med(.cb_cols[["nUMI"]],  fmt_int)
  output$med_genes <- med(.cb_cols[["nGene"]], fmt_int)
  output$med_mito  <- med(.cb_cols[["pct_mito"]], function(v) sprintf("%.1f%%", v))
  output$pct_doublet <- shiny::renderText({
    if (is.null(rv$qc)) return("-")
    cl <- coldata_df(rv$qc)[[.cb_cols[["doublet_class"]]]]
    if (all(is.na(cl))) return("n/a")
    sprintf("%.1f%%", 100 * mean(cl == "doublet", na.rm = TRUE))
  })

  # ---- QC overview plots --------------------------------------------------
  violin <- function(metric) plotly::renderPlotly({
    shiny::req(rv$qc)
    plotly::ggplotly(qc_violin(filt()$sce, metric, thresholds())) |>
      plotly::config(displaylogo = FALSE)
  })
  output$violin_mito    <- violin("pct_mito")
  output$violin_genes   <- violin("nGene")
  output$violin_umi     <- violin("nUMI")
  output$violin_doublet <- violin("doublet_score")

  # ---- Filtering value boxes ---------------------------------------------
  output$f_total   <- shiny::renderText(if (is.null(rv$qc)) "-" else fmt_int(filt()$n_total))
  output$f_keep    <- shiny::renderText(if (is.null(rv$qc)) "-" else fmt_int(filt()$n_keep))
  output$f_discard <- shiny::renderText(if (is.null(rv$qc)) "-" else fmt_int(filt()$n_discard))
  output$f_pct     <- shiny::renderText({
    if (is.null(rv$qc)) return("-")
    f <- filt()
    if (!f$n_total) return("-")
    sprintf("%.1f%%", 100 * f$n_keep / f$n_total)
  })

  # ---- Filtering plots & tables ------------------------------------------
  output$scatter <- plotly::renderPlotly({
    shiny::req(rv$qc)
    plotly::ggplotly(
      qc_scatter(filt()$sce, "nUMI", "nGene", colour = "status",
                 thresholds = thresholds())
    ) |> plotly::config(displaylogo = FALSE)
  })

  output$reason_table <- DT::renderDT({
    shiny::req(rv$qc)
    br <- summarise_qc(filt())$by_reason
    if (!nrow(br)) {
      br <- data.frame(Filter = "(no active filters)", `Cells failing` = 0,
                       check.names = FALSE)
    } else {
      br <- data.frame(Filter = br$description,
                       `Cells failing` = br$cells_failing,
                       check.names = FALSE)
    }
    DT::datatable(br, rownames = FALSE,
                  options = list(dom = "t", ordering = FALSE))
  })

  output$cell_table <- DT::renderDT({
    shiny::req(rv$qc)
    cd <- coldata_df(filt()$sce)
    df <- data.frame(
      Barcode = rownames(cd),
      nUMI = cd[[.cb_cols[["nUMI"]]]],
      nGene = cd[[.cb_cols[["nGene"]]]],
      pct_mito = round(cd[[.cb_cols[["pct_mito"]]]], 2),
      doublet_score = round(cd[[.cb_cols[["doublet_score"]]]], 3),
      doublet_class = cd[[.cb_cols[["doublet_class"]]]],
      status = ifelse(cd$discard, "Discarded", "Retained"),
      stringsAsFactors = FALSE
    )
    DT::datatable(df, rownames = FALSE, filter = "top",
                  options = list(pageLength = 10, scrollX = TRUE)) |>
      DT::formatStyle("status",
        backgroundColor = DT::styleEqual(
          c("Retained", "Discarded"), c("#E3F0D4", "#FBE3D6")))
  })

  # ---- Export -------------------------------------------------------------
  output$dl_sce <- shiny::downloadHandler(
    filename = function() "cellboard_filtered_sce.rds",
    content = function(file) {
      shiny::req(rv$qc)
      saveRDS(filt()$filtered, file)
    }
  )

  output$dl_csv <- shiny::downloadHandler(
    filename = function() "cellboard_cell_metrics.csv",
    content = function(file) {
      shiny::req(rv$qc)
      cd <- coldata_df(filt()$sce)
      keep_cols <- intersect(
        c(.cb_cols, "discard"),
        names(cd))
      out <- data.frame(barcode = rownames(cd), cd[, keep_cols, drop = FALSE])
      utils::write.csv(out, file, row.names = FALSE)
    }
  )

  output$dl_report <- shiny::downloadHandler(
    filename = function() "cellboard_qc_report.html",
    content = function(file) {
      shiny::req(rv$qc)
      shiny::withProgress(message = "Rendering report", value = 0.3, {
        path <- render_qc_report(
          filt()$sce, thresholds(),
          output_file = basename(tempfile(fileext = ".html")),
          title = "cellboard QC report",
          sample_name = rv$source %||% "dataset"
        )
        shiny::incProgress(0.6)
        file.copy(path, file, overwrite = TRUE)
      })
    }
  )

  output$settings_summary <- shiny::renderPrint({
    if (is.null(rv$qc)) { cat("No QC computed yet.\n"); return(invisible()) }
    th <- thresholds(); f <- filt()
    cat("Data source:", rv$source %||% "-", "\n")
    cat(sprintf("Cells: %d input -> %d retained (%.1f%%), %d discarded\n",
                f$n_total, f$n_keep, 100 * f$n_keep / f$n_total, f$n_discard))
    cat("\nActive thresholds:\n")
    utils::str(Filter(Negate(is.null), unclass(th)))
  })

  # ---- About --------------------------------------------------------------
  output$about <- shiny::renderUI({
    ver <- tryCatch(as.character(utils::packageVersion("cellboard")),
                    error = function(e) "dev")
    shiny::HTML(sprintf('
      <p><strong>cellboard</strong> (v%s) is an interactive dashboard for
      quality control of single-cell RNA-seq data.</p>
      <ol>
        <li><strong>Load</strong> a 10x, SingleCellExperiment or Seurat dataset
            (or the bundled example).</li>
        <li><strong>Compute QC</strong>: mitochondrial %%, genes, UMIs and
            scDblFinder doublet scores.</li>
        <li><strong>Filter</strong> with live sliders and watch retained vs
            discarded cells update.</li>
        <li><strong>Export</strong> the filtered object and a reproducible
            report.</li>
      </ol>
      <p>Built on Bioconductor (SingleCellExperiment, scDblFinder) and deployable
      to Posit Connect. Licensed under GPL-3.</p>', ver))
  })
}
