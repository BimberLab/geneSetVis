#' Launch geneSetVis.
#' @title Launch geneSetVis.
#' @description This will launch the geneSetVis app.
#' @keywords gene-set enrichment
#' @export
#' @return Shiny application.
#' @rawNamespace import(shiny, except = c('dataTableOutput', 'renderDataTable'))
#' @import shinydashboard
#' @import shinyWidgets
#' @import shinytest
#' @importFrom dplyr %>% arrange group_by rename select summarise
#' @import ggplot2
#' @import ggupset
#' @import tidyr
#' @import stringr
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom DT datatable renderDataTable
#' @import formattable
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
#' @import shinybusy
#'
launchGeneSetVis <- function(...) {
  gsvis_package <<- TRUE

  #Load libraries via global, to ensure we stay in sync:
  globalFile <- system.file("app/global.R", package = "geneSetVis")
  source(globalFile)

  shinyAppDir(appDir = system.file("app", package = "geneSetVis"))
}


