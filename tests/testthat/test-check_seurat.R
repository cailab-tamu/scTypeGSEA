dta <- check_seurat(pbmc_raw, min.cells = 1, min.features = 10)
info <- try(Seurat::Project(dta), silent = TRUE)
test_that("check_seurat works", {
  # check it becomes a Seurat object
  expect_equal(grepl("Error", info), FALSE)
})
