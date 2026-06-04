# Shiny server logic, exercised headlessly with shiny::testServer().

test_that("server wires load -> compute -> filter", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("shiny")

  shiny::testServer(app_server, {
    session$setInputs(load_example = 1)
    session$setInputs(run_doublets = FALSE, mito_pattern = "^[Mm][Tt]-")
    session$setInputs(compute = 1)
    session$setInputs(max_mito = 100)        # don't discard anything by mito

    expect_gt(filt()$n_total, 0)
    expect_lte(filt()$n_keep, filt()$n_total)
    expect_match(output$f_total, "[0-9]")
    expect_match(output$med_mito, "%")
  })
})

test_that("server filters more cells with a stricter mito threshold", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("shiny")

  shiny::testServer(app_server, {
    session$setInputs(load_example = 1, run_doublets = FALSE)
    session$setInputs(compute = 1)
    session$setInputs(max_mito = 100)
    keep_loose <- filt()$n_keep
    session$setInputs(max_mito = 10)
    expect_lte(filt()$n_keep, keep_loose)
  })
})

test_that("app_ui builds a page object", {
  skip_if_not_installed("bslib")
  ui <- app_ui()
  expect_true(inherits(ui, c("bslib_page", "shiny.tag", "shiny.tag.list")))
})
