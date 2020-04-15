#' Launch geneSetVis.
#' @title Launch geneSetVis.
#' @description Launch geneSetVis.
#' @keywords gene-set enrichment
#' @export
#' @return Shiny application.
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import shinytest
#' @import testthat
#' @importFrom dplyr arrange group_by rename select summarise
#' @import ggplot2
#' @import ggupset
#' @import tidyr
#' @import stringr
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom DT datatable renderDataTable
#' @import formattable
#' @import stats
#' @import knitr
#' @importFrom rmdformats html_clean
#' @import org.Hs.eg.db
#' @import org.Mmu.eg.db
#' @import org.Mm.eg.db
#' @import AnnotationHub
#' @import ReactomePA
#' @import clusterProfiler
#' @import DOSE
#' @import enrichplot
#' @import fgsea
#' @import enrichR
#' @importFrom utils stack
#'
#setwd("/Users/onwuzu/Documents/BimberLab/geneSetVis/")
launchGeneSetVis <- function(...)
{
  shiny::runApp(appDir = system.file("app", package = "geneSetVis"),
                ...)
}
