#' @title scDotPlot
#' @description Create dot plot of gene expression profiles that can be
#'  annotated with hierarchical clustering from \link[ggtree]{ggtree} using
#'  \link[aplot:plot-insertion]{aplot}
#' @param object An object with normalized data
#' \itemize{
#'  \item \link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}
#'  \item \link[SeuratObject:Seurat-class]{Seurat}
#'  \item A data.frame with the following columns: "NumDetected", "Feature", "Group", "Average"
#' }
#' @param features Character vector with genes of interest
#' @param group Column name from colData/metadata of the object to group cells by
#' @param block Column name from colData of a SingleCellExperiment object
#'  to use as a blocking factor (e.g. batch or sample)
#' @param swap_rownames Column name from rowData of a SingleCellExperiment
#'  object to match to features
#' @param scale Logical indicating whether the data should be scaled and centered
#' @param AverageThreshold Numeric specifying threshold for average expression,
#'  where values below AverageThreshold and NumDetectedThreshold are transparent
#' @param NumDetectedThreshold Numeric specifying threshold for fraction of cells,
#'  where values below AverageThreshold and NumDetectedThreshold are transparent
#' @param cluster Logical specifying whether to perform hierarchical clustering
#'  analysis
#' @param groupAnno Cell annotations that are stored as names of columns
#'  in colData of sce with annotations
#' @param featureAnno Feature annotations that are stored as names of rows
#'  in rowData of sce with annotations
#' @param treeWidth Numeric specifying width of the row tree relative to
#'  the dotPlot
#' @param treeHeight Numeric specifying height of the column tree relative to
#'  the dotPlot
#' @param annoWidth Numeric specifying the width of the row annotation relative
#'  to the dotPlot
#' @param annoHeight Numeric specifying height of the column annotation relative
#'  to the dotPlot
#' @param annoColors A list with a name for each annotation that contains
#'  a named vector of colors, where the name is the pairing of values to colors
#' @param dotColors A character vector specifying the colors to be used
#'  in the gradient for the dots. If scale is set to TRUE, the first 3 colors
#'  will be used for the negative, zero, and positive values, respectively.
#' @param groupLegends Logical indicating whether to show legends for group
#'  annotations
#' @param featureLegends Logical indicating whether to show legends for feature
#'  annotations
#' @param fontSize Numeric specifying the base font size in pts
#' @param fontFamily Character specifying the base font family
#' @param flipPlot Logical indicating whether to flip the x and y coordinates
#' @param ... Additional unused arguments
#' @details The data for the dot plot is generated from different sources
#'  depending on the object:
#' \itemize{
#'  \item \code{SingleCellExperiment:} \link[scater]{plotDots}
#'  \item \code{Seurat:} \link[Seurat]{DotPlot}
#' }
#' @return
#'  \itemize{
#'      \item When \code{cluster = TRUE}, a \code{aplot} object
#'      \item When \code{cluster = FALSE}, a \code{ggplot2} object
#' }
#' @references \url{https://yulab-smu.top/pkgdocs/aplot.html#a-single-cell-example}
#' @export
#'
scDotPlot <- function(object,
                      ...){
    UseMethod("scDotPlot")
}

#' @rdname scDotPlot
#' @name scDotPlot
#' @order 1
#' @importFrom SingleCellExperiment rowData colData
#' @importFrom rlang syms sym
#' @importFrom purrr set_names
#' @importFrom scater plotDots
#' @importFrom tibble as_tibble
#' @importFrom BiocGenerics as.data.frame
#' @importFrom dplyr pull rename mutate left_join select distinct
#' @importFrom stringr str_sort
#' @import ggplot2
#' @importFrom cli cli_abort
#' @importFrom magrittr %>% %$%
#' @examples
#' data("pbmc_small", package = "SeuratObject")
#' pbmc_small |>
#'     scDotPlot(features = Seurat::VariableFeatures(pbmc_small),
#'               group = "RNA_snn_res.1")
#' @export
#'
scDotPlot.SingleCellExperiment <- function(object,
                                           features = features,
                                           group = NULL,
                                           block = NULL,
                                           swap_rownames = NULL,
                                           scale = FALSE,
                                           cluster = TRUE,
                                           AverageThreshold = ifelse(scale == FALSE, 0, -Inf),
                                           NumDetectedThreshold = 0.01,
                                           groupAnno = FALSE,
                                           featureAnno = FALSE,
                                           treeWidth = 0.1,
                                           treeHeight = 0.1,
                                           annoWidth = 0.05,
                                           annoHeight = 0.02,
                                           annoColors = NULL,
                                           dotColors = NULL,
                                           groupLegends = TRUE,
                                           featureLegends = TRUE,
                                           fontSize = 11,
                                           fontFamily = "",
                                           flipPlot = FALSE,
                                           ...){

    featureCheck <- if(is.null(swap_rownames)){
        rownames(object)
    }else if(swap_rownames %in% colnames(SingleCellExperiment::rowData(object))){
        SingleCellExperiment::rowData(object)[ ,swap_rownames]
    }

    if(!all(features %in% featureCheck)){
        missingFeatures <- features[!features %in% featureCheck]
        cli::cli_abort(c("The following features are not in the rowData of the object:",
                         missingFeatures %>%
                             purrr::set_names(rep("x", length(missingFeatures)))))
    }

    object %>%
        scater::plotDots(features = features,
                         group = group,
                         block = block,
                         swap_rownames = swap_rownames,
                         scale = scale,
                         center = scale) %$%
        data %>%
        tibble::as_tibble() %>%
        {if(!isFALSE(groupAnno)){
            dplyr::left_join(.,
                             object %>%
                                 SingleCellExperiment::colData() %>%
                                 BiocGenerics::as.data.frame() %>%
                                 dplyr::select(!!!rlang::syms(groupAnno),
                                               Group = !!rlang::sym(group)) %>%
                                 dplyr::distinct(),
                             by = "Group")
        }else{
            .
        }} %>%
        {if(group %in% groupAnno){
            dplyr::mutate(., !!group := Group)
        }else{
            .
        }} %>%
        {if(!isFALSE(featureAnno)){
            dplyr::left_join(.,
                             object %>%
                                 SingleCellExperiment::rowData() %>%
                                 BiocGenerics::as.data.frame() %>%
                                 tibble::rownames_to_column("rownames") %>%
                                 dplyr::select(!!!rlang::syms(featureAnno),
                                               Feature = !!rlang::sym(ifelse(!is.null(swap_rownames),
                                                                             swap_rownames,
                                                                             "rownames"))) %>%
                                 dplyr::distinct(),
                             by = "Feature")
        }else{
            .
        }} %>%
        scDotPlot.default(group = NULL,
                          features = features,
                          scale = scale,
                          cluster = cluster,
                          AverageThreshold = AverageThreshold,
                          NumDetectedThreshold = NumDetectedThreshold,
                          groupAnno = groupAnno,
                          featureAnno = featureAnno,
                          treeWidth = treeWidth,
                          treeHeight = treeHeight,
                          annoWidth = annoWidth,
                          annoHeight = annoHeight,
                          annoColors = annoColors,
                          dotColors = dotColors,
                          groupLegends = groupLegends,
                          featureLegends = featureLegends,
                          fontSize = fontSize,
                          fontFamily = fontFamily,
                          flipPlot = flipPlot,
                          ...)
}

#' @rdname scDotPlot
#' @name scDotPlot
#' @order 2
#' @importFrom Seurat DotPlot FetchData DefaultAssay
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom dplyr rename mutate left_join select distinct
#' @importFrom rlang sym syms
#' @importFrom cli cli_abort
#' @importFrom magrittr %>% %$%
#' @export
#'
scDotPlot.Seurat <- function(object,
                             features = features,
                             group = NULL,
                             block = NULL,
                             swap_rownames = NULL,
                             scale = FALSE,
                             cluster = TRUE,
                             AverageThreshold = ifelse(scale == FALSE, 0, -Inf),
                             NumDetectedThreshold = 0.01,
                             groupAnno = FALSE,
                             featureAnno = FALSE,
                             treeWidth = 0.1,
                             treeHeight = 0.1,
                             annoWidth = 0.05,
                             annoHeight = 0.02,
                             annoColors = NULL,
                             dotColors = NULL,
                             groupLegends = TRUE,
                             featureLegends = TRUE,
                             fontSize = 11,
                             fontFamily = "",
                             flipPlot = FALSE,
                             ...){
    if(!is.null(block) | !is.null(swap_rownames)){
        cli::cli_abort("The {.field block} and {.field swap_rownames} \
               arguments are not supported for Seurat objects.\
               See {.fn Seurat::as.SingleCellExperiment} \
               for information about how to convert the object.")
    }

    object %>%
        Seurat::DotPlot(object = .,
                        features = features,
                        group.by = group,
                        scale = scale) %$%
        data %>%
        tibble::as_tibble() %>%
        dplyr::rename(Average = avg.exp.scaled,
                      NumDetected = pct.exp,
                      Feature = features.plot,
                      Group = id) %>%
        dplyr::mutate(NumDetected = NumDetected/100) %>%
        {if(!isFALSE(groupAnno)){
            dplyr::left_join(.,
                             object %>%
                                 Seurat::FetchData(c(groupAnno, group)) %>%
                                 tibble::as_tibble() %>%
                                 dplyr::rename("Group" = !!rlang::sym(group)) %>%
                                 dplyr::distinct(),
                             by = "Group",
                             relationship = "many-to-many")
        }else{
            .
        }} %>%
        {if(group %in% groupAnno){
            dplyr::mutate(., !!group := Group)
        }else{
            .
        }} %>%
        {if(!isFALSE(featureAnno)){
            dplyr::left_join(.,
                             object[[Seurat::DefaultAssay(object)]][[]] %>%
                                 tibble::rownames_to_column("rownames") %>%
                                 dplyr::select(!!!rlang::syms(featureAnno),
                                               Feature = "rownames") %>%
                                 dplyr::distinct(),
                             by = "Feature")
        }else{
            .
        }} %>%
        scDotPlot.default(group = NULL,
                          features = features,
                          scale = scale,
                          cluster = cluster,
                          AverageThreshold = AverageThreshold,
                          NumDetectedThreshold = NumDetectedThreshold,
                          groupAnno = groupAnno,
                          featureAnno = featureAnno,
                          treeWidth = treeWidth,
                          treeHeight = treeHeight,
                          annoWidth = annoWidth,
                          annoHeight = annoHeight,
                          annoColors = annoColors,
                          dotColors = dotColors,
                          groupLegends = groupLegends,
                          featureLegends = featureLegends,
                          fontSize = fontSize,
                          fontFamily = fontFamily,
                          flipPlot = flipPlot,
                          ...)
}

#' @rdname scDotPlot
#' @name scDotPlot
#' @order 3
#' @export
#'
scDotPlot.default <- function(object,
                              features = NULL,
                              group = NULL,
                              block = NULL,
                              swap_rownames = NULL,
                              scale = FALSE,
                              cluster = TRUE,
                              AverageThreshold = ifelse(scale == FALSE, 0, -Inf),
                              NumDetectedThreshold = 0.01,
                              groupAnno = FALSE,
                              featureAnno = FALSE,
                              treeWidth = 0.1,
                              treeHeight = 0.1,
                              annoWidth = 0.05,
                              annoHeight = 0.02,
                              annoColors = NULL,
                              dotColors = NULL,
                              groupLegends = TRUE,
                              featureLegends = TRUE,
                              fontSize = 11,
                              fontFamily = "",
                              flipPlot = FALSE,
                              ...){
    dotPlot <- object %>%
        {if(!is.null(features)){
            dplyr::mutate(., Feature = factor(Feature, levels = features))
        }else{
            .
        }} %>%
        .baseDotPlot(group = group,
                     scale = scale,
                     AverageThreshold = AverageThreshold,
                     NumDetectedThreshold = NumDetectedThreshold,
                     dotColors = dotColors,
                     fontSize = fontSize,
                     fontFamily = fontFamily,
                     flipPlot = flipPlot)

    if(any(cluster == TRUE, !is.null(groupAnno), !is.null(featureAnno))){
        dotPlot %>%
            .annotateDotPlot(cluster = cluster,
                             groupAnno = groupAnno,
                             featureAnno = featureAnno,
                             treeWidth = treeWidth,
                             treeHeight = treeHeight,
                             annoWidth = annoWidth,
                             annoHeight = annoHeight,
                             annoColors = annoColors,
                             groupLegends = groupLegends,
                             featureLegends = featureLegends,
                             fontSize = fontSize,
                             fontFamily = fontFamily,
                             flipPlot = flipPlot)
    }else{
        dotPlot
    }
}
