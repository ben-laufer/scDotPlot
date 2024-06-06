test_that("base plot works", {
    data("pbmc_small", package = "SeuratObject")
    pbmc_small |>
        scDotPlot(cluster = FALSE,
                  features = Seurat::VariableFeatures(pbmc_small),
                  group = "RNA_snn_res.1") |>
        expect_s3_class("ggplot")
})

test_that("plot annotation works", {
    data("pbmc_small", package = "SeuratObject")
    pbmc_small |>
        scDotPlot(features = Seurat::VariableFeatures(pbmc_small),
                  group = "RNA_snn_res.1") |>
        expect_s3_class("aplot")
})
