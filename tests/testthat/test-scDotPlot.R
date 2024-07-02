test_that("base plot works", {
    data("pbmc_small", package = "SeuratObject")
    pbmc_small %>%
        scDotPlot(cluster = FALSE,
                  features = VariableFeatures(.),
                  group = "RNA_snn_res.1") %T>%
        expect_s3_class("ggplot") %>%
        expect_doppelganger("base plot", .)
})

test_that("plot annotation works", {
    data("pbmc_small", package = "SeuratObject")
    pbmc_small %>%
        scDotPlot(features = VariableFeatures(.),
                  group = "RNA_snn_res.1") %T>%
        expect_s3_class("aplot") %>%
        expect_doppelganger("plot annotation", .)
})
