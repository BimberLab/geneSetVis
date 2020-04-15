#' Launch geneSetVis.
#' @title Launch geneSetVis.
#' @description Launch geneSetVis.
#' @keywords gene-set enrichment
#' @export
#' @return Shiny application.
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import flexdashboard
#' @import shinytest
#' @import testthat
#' @import dplyr
#' @import ggplot2
#' @import ggupset
#' @import tidyr
#' @import stringr
#' @import plotly
#' @import DT
#' @import formattable
#' @import stats
#' @import knitr
#' @import rmdformats
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
#'
#setwd("/Users/onwuzu/Documents/BimberLab/geneSetVis/")
launchGeneSetVis <- function(...)
{
  shiny::runApp(appDir = system.file(".", package = "geneSetVis"),
                ...)
}
