#' @title Cluster Dot Plot
#' @inheritParams scDotPlot
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom stats dist hclust
#' @importFrom rlang sym
#' @importFrom purrr set_names map reduce
#' @import ggplot2
#' @import ggtree
#' @import aplot
#' @importFrom magrittr %$% %>%
#' @return A aplot
#' @keywords internal
#'
.annotateDotPlot <- function(dotPlot,
                             cluster = TRUE,
                             groupAnno = FALSE,
                             featureAnno = FALSE,
                             treeWidth = 0.1,
                             treeHeight = 0.1,
                             annoWidth = 0.05,
                             annoHeight = 0.02,
                             annoColors = NULL,
                             groupLegends = TRUE,
                             featureLegends = TRUE,
                             flipPlot = FALSE){

    plotMatrix <- dotPlot %$%
        data %>%
        tibble::as_tibble() %>%
        dplyr::select(Feature, Group, Average) %>%
        tidyr::pivot_wider(names_from = Group, values_from = Average) %>%
        tibble::column_to_rownames("Feature") %>%
        as.matrix()

    if(!all(groupAnno == FALSE)){
        colLabels <- groupAnno %>%
            purrr::set_names() %>%
            purrr::map(.createAnno,
                       annoType = "column",
                       dotPlot = dotPlot,
                       annoColors = annoColors,
                       groupLegends = groupLegends,
                       flipPlot = flipPlot)
    }

    if(!all(featureAnno == FALSE)){
        rowLabels <- featureAnno %>%
            purrr::set_names() %>%
            purrr::map(.createAnno,
                       annoType = "row",
                       dotPlot = dotPlot,
                       annoColors = annoColors,
                       featureLegends = featureLegends,
                       flipPlot = flipPlot)
    }

    treeRow <- plotMatrix %>%
        {if(flipPlot == FALSE){
            .
        }else if(flipPlot == TRUE){
            t(.)
        }} %>%
        stats::dist() %>%
        stats::hclust() %>%
        ggtree::ggtree(branch.length = "none")

    treeCol <- plotMatrix %>%
        {if(flipPlot == FALSE){
            t(.)
        }else if(flipPlot == TRUE){
            .
        }} %>%
        stats::dist() %>%
        stats::hclust() %>%
        ggtree::ggtree(branch.length = "none") +
        ggtree::layout_dendrogram()

    (dotPlot +
            ggplot2::scale_y_discrete(position = "right")) %>%
        {if(!all(groupAnno == FALSE)){
            purrr::reduce(ifelse(flipPlot == FALSE, colLabels, rowLabels),
                          ~ aplot::insert_top(.x, .y, height = annoHeight),
                          .init = .)
        }else{
            .
        }} %>%
        {if(!all(featureAnno == FALSE)){
            purrr::reduce(ifelse(flipPlot == FALSE, rowLabels, colLabels),
                          ~ aplot::insert_left(.x, .y, width = annoWidth),
                          .init = .)
        }else{
            .
        }} %>%
        {if(cluster == TRUE){
            aplot::insert_left(., treeRow, width = treeWidth) %>%
                aplot::insert_top(treeCol, height = treeHeight)
        }else{
            .
        }}
}

#' @title Dot Plot Base
#' @inheritParams scDotPlot
#' @importFrom dplyr mutate case_when
#' @importFrom rlang sym
#' @importFrom scales muted
#' @import ggplot2
#' @importFrom magrittr %>%
#' @return A ggplot2
#' @keywords internal
#'
.baseDotPlot <- function(plotData,
                         group = NULL,
                         scale = NULL,
                         AverageThreshold = NULL,
                         NumDetectedThreshold = NULL,
                         dotColors = NULL,
                         flipPlot = FALSE){

    stopifnot(c("NumDetected", "Feature", "Group", "Average") %in% names(plotData))
    if(is.null(group)){
        group <- "Group"
    }

    if(is.null(dotColors)){
        if(isTRUE(scale)){
            dotColors <- c(scales::muted("blue"), "white", scales::muted("red"))
        }else{
            dotColors <- c("#BFBFBF", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF")
        }
    }

    (plotData %>%
            dplyr::mutate(alpha = dplyr::case_when(Average > !!AverageThreshold &
                                                       NumDetected > NumDetectedThreshold ~ 1,
                                                   .default = 0.1)) %>%
            {if(flipPlot == FALSE){
                ggplot2::ggplot(., ggplot2::aes(x = !!rlang::sym(group),
                                                y = Feature))
            }else if(flipPlot == TRUE){
                ggplot2::ggplot(., ggplot2::aes(x = Feature,
                                                y = !!rlang::sym(group)))
            }} +
            ggplot2::aes(size = NumDetected,
                         fill = Average,
                         alpha = I(alpha)) +
            ggplot2::geom_point(shape = 21,
                                color = "black") +
            ggplot2::scale_size_continuous(labels = scales::label_percent()) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.title = ggplot2::element_blank(),
                           axis.line = ggplot2::element_blank(),
                           axis.ticks = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(angle = 90,
                                                               vjust = 0.5,
                                                               hjust = 1)) +
            ggplot2::guides(fill = ggplot2::guide_colorbar(title = ifelse(scale == FALSE,
                                                                          "Log Counts",
                                                                          "Scaled Expression"),
                                                           frame.colour = "black",
                                                           ticks.colour = "black",
                                                           ticks.linewidth = 0.5,
                                                           order = 1),
                            size = ggplot2::guide_legend(title = "Percent of Cells",
                                                         order = 2)) +
            ggplot2::ylab(NULL)) %>%
        {if(scale == TRUE){
            . + ggplot2::scale_fill_gradient2(low = dotColors[1],
                                              mid = dotColors[2],
                                              high = dotColors[3])
        }else if(scale == FALSE){
            . + ggplot2::scale_fill_gradientn(colors = dotColors)
        }}
}

#' @title Create column annotations
#' @inheritParams scDotPlot
#' @importFrom ggsci scale_fill_d3 scale_fill_cosmic
#' @import ggplot2
#' @importFrom purrr pluck
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %$% %>%
#' @return A ggplot2
#' @keywords internal
#'
.createAnno <- function(annoLabel,
                        annoType = c("column", "row"),
                        dotPlot = dotPlot,
                        annoColors = NULL,
                        groupLegends = TRUE,
                        featureLegends = TRUE,
                        flipPlot = FALSE){

    annoType <- match.arg(annoType)

    plotData <- dotPlot %$%
        data

    p1 <- plotData %>%
        {if(annoType == "column" & flipPlot == FALSE){
            ggplot2::ggplot(., ggplot2::aes(x = Group,
                                            y = 1,
                                            fill = !!rlang::sym(annoLabel)))
        }else if(annoType == "column" & flipPlot == TRUE){
            ggplot2::ggplot(., ggplot2::aes(x = 1,
                                            y = Group,
                                            fill = !!rlang::sym(annoLabel)))
        }else if(annoType == "row" & flipPlot == FALSE){
            ggplot2::ggplot(., ggplot2::aes(x = 1,
                                            y = Feature,
                                            fill = !!rlang::sym(annoLabel)))
        }else if(annoType == "row" & flipPlot == TRUE){
            ggplot2::ggplot(., ggplot2::aes(x = Feature,
                                            y = 1,
                                            fill = !!rlang::sym(annoLabel)))
        }} +
        ggplot2::geom_tile(color = "black",
                           linewidth = 0.25) +
        ggplot2::theme_void()

    labelLevels <- plotData %>%
        purrr::pluck(annoLabel) %>%
        unique() %>%
        length()

    if(annoLabel %in% names(annoColors)){
        p1 <- p1 + ggplot2::scale_fill_manual(values = annoColors[[annoLabel]])
    }else if(annoType == "column"){
        if(labelLevels > 3 & labelLevels < 20){
            p1 <- p1 + ggsci::scale_fill_d3("category20")
        }else if(labelLevels <= 3){
            p1 <- p1 + ggsci::scale_fill_cosmic("signature_substitutions")
        }else if(labelLevels > 20){
            pal <- grDevices::colorRampPalette(ggsci::pal_d3("category20b")(20))(labelLevels)
            p1 <- p1 + ggplot2::scale_fill_manual(values = pal)
        }
    }else if(annoType == "row"){
        p1 <- p1 + ggplot2::scale_fill_brewer(palette = "Dark2")
    }

    if((annoType == "column" & groupLegends == FALSE) |
       (annoType == "row" & featureLegends == FALSE)){
        p1 <- p1 + ggplot2::guides(fill = ggplot2::guide_none())
    }

    p1
}
