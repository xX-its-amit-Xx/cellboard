# The dashboard UI. Built with bslib so it themes cleanly and lays out
# responsively. Dynamic, data-dependent controls (the threshold sliders) and
# all metric values are rendered server-side via uiOutput()/textOutput().

#' cellboard Shiny UI
#'
#' @return A bslib page definition.
#' @export
app_ui <- function() {
  theme <- bslib::bs_theme(
    version = 5, preset = "cosmo",
    primary = "#2C7FB8", "navbar-bg" = "#1f4f6b"
  )

  sidebar <- bslib::sidebar(
    width = 330, title = "Data & thresholds",
    bslib::accordion(
      open = c("Load data", "Thresholds"),
      bslib::accordion_panel(
        "Load data", icon = shiny::icon("upload"),
        shiny::actionButton(
          "load_example", "Load example PBMC",
          icon = shiny::icon("vial"),
          class = "btn-outline-primary btn-sm w-100 mb-2"
        ),
        shiny::fileInput(
          "file", "...or upload your own",
          multiple = TRUE,
          accept = c(".rds", ".h5", ".hdf5", ".mtx", ".gz", ".tsv"),
          width = "100%"
        ),
        shiny::helpText(
          "Accepts a .rds (SingleCellExperiment or Seurat), a 10x .h5, ",
          "or the three 10x files (matrix.mtx.gz, barcodes.tsv.gz, ",
          "features.tsv.gz) selected together."
        )
      ),
      bslib::accordion_panel(
        "Compute QC", icon = shiny::icon("calculator"),
        shiny::textInput("mito_pattern", "Mitochondrial gene pattern",
                         value = "^[Mm][Tt]-"),
        shiny::checkboxInput("run_doublets",
                             "Run doublet detection (scDblFinder)", TRUE),
        shiny::actionButton(
          "compute", "Compute QC metrics",
          icon = shiny::icon("gears"),
          class = "btn-primary btn-sm w-100"
        )
      ),
      bslib::accordion_panel(
        "Thresholds", icon = shiny::icon("sliders"),
        shiny::uiOutput("threshold_ui")
      )
    )
  )

  bslib::page_navbar(
    title = bslib::tooltip(
      shiny::span(shiny::icon("dna"), "cellboard"),
      "Interactive single-cell RNA-seq QC"
    ),
    id = "nav", theme = theme, sidebar = sidebar,
    fillable = "QC overview",

    # ---- Data tab ----------------------------------------------------------
    bslib::nav_panel(
      "Data", icon = shiny::icon("table"),
      bslib::layout_columns(
        fill = FALSE,
        cb_value_box("Genes", "n_genes", "dna", "primary"),
        cb_value_box("Cells", "n_cells", "circle-nodes", "primary"),
        cb_value_box("Source", "data_source", "file-import", "secondary")
      ),
      bslib::card(
        bslib::card_header("Dataset"),
        shiny::uiOutput("data_status"),
        bslib::card_body(DT::DTOutput("rowdata_preview"), max_height = "320px")
      )
    ),

    # ---- QC overview tab ---------------------------------------------------
    bslib::nav_panel(
      "QC overview", icon = shiny::icon("chart-area"),
      bslib::layout_columns(
        fill = FALSE,
        cb_value_box("Median UMIs", "med_umi", "layer-group", "info"),
        cb_value_box("Median genes", "med_genes", "dna", "info"),
        cb_value_box("Median mito %", "med_mito", "fire", "warning"),
        cb_value_box("Doublets", "pct_doublet", "clone", "danger")
      ),
      bslib::layout_columns(
        col_widths = c(6, 6, 6, 6),
        cb_plot_card("Mitochondrial %", "violin_mito"),
        cb_plot_card("Genes per cell", "violin_genes"),
        cb_plot_card("UMIs per cell", "violin_umi"),
        cb_plot_card("Doublet score", "violin_doublet")
      )
    ),

    # ---- Filtering tab -----------------------------------------------------
    bslib::nav_panel(
      "Filtering", icon = shiny::icon("filter"),
      bslib::layout_columns(
        fill = FALSE,
        cb_value_box("Input cells", "f_total", "circle-nodes", "secondary"),
        cb_value_box("Retained", "f_keep", "circle-check", "success"),
        cb_value_box("Discarded", "f_discard", "circle-xmark", "danger"),
        cb_value_box("% retained", "f_pct", "percent", "primary")
      ),
      bslib::layout_columns(
        col_widths = c(7, 5),
        cb_plot_card("UMIs vs genes (retained vs discarded)", "scatter"),
        bslib::card(
          bslib::card_header("Cells failing each filter"),
          DT::DTOutput("reason_table")
        )
      ),
      bslib::card(
        bslib::card_header("Per-cell QC metrics"),
        bslib::card_body(DT::DTOutput("cell_table"))
      )
    ),

    # ---- Export tab --------------------------------------------------------
    bslib::nav_panel(
      "Export", icon = shiny::icon("download"),
      bslib::layout_columns(
        col_widths = c(6, 6),
        bslib::card(
          bslib::card_header("Download filtered data"),
          bslib::card_body(
            shiny::p("Export the QC'd object and per-cell metrics."),
            shiny::downloadButton("dl_sce",
              "Filtered SingleCellExperiment (.rds)",
              class = "btn-outline-primary w-100 mb-2"),
            shiny::downloadButton("dl_csv",
              "Per-cell metrics (.csv)",
              class = "btn-outline-primary w-100")
          )
        ),
        bslib::card(
          bslib::card_header("Reproducible report"),
          bslib::card_body(
            shiny::p("Render a self-contained HTML report of the current ",
                     "thresholds and filtering outcome."),
            shiny::downloadButton("dl_report",
              "Generate QC report (.html)",
              class = "btn-primary w-100")
          )
        )
      ),
      bslib::card(
        bslib::card_header("Current settings"),
        shiny::verbatimTextOutput("settings_summary")
      )
    ),

    # ---- About tab ---------------------------------------------------------
    bslib::nav_panel(
      "About", icon = shiny::icon("circle-info"),
      bslib::card(
        bslib::card_header("About cellboard"),
        bslib::card_body(shiny::uiOutput("about"))
      )
    ),
    bslib::nav_spacer(),
    bslib::nav_item(
      shiny::tags$a(
        shiny::icon("github"), "Source",
        href = "https://github.com/ashenoy/cellboard", target = "_blank"
      )
    )
  )
}

# ---- small UI builders ------------------------------------------------------

# A value box whose value is filled by a server-side textOutput.
cb_value_box <- function(title, output_id, icon, theme) {
  bslib::value_box(
    title = title,
    value = shiny::textOutput(output_id, inline = TRUE),
    showcase = shiny::icon(icon),
    theme = theme
  )
}

# A card wrapping a full-bleed plotly output.
cb_plot_card <- function(title, output_id) {
  bslib::card(
    full_screen = TRUE,
    bslib::card_header(title),
    bslib::card_body(plotly::plotlyOutput(output_id, height = "100%"),
                     padding = 0)
  )
}
