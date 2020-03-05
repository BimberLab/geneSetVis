library(shiny)
library(shinydashboard)
library(flexdashboard)

library(dplyr)
library(ggplot2)
library(ggupset)
library(tidyr)
library(reshape2)
library(stringr)
library(readxl) 
library(datapasta)
library(plotly)
library(cowplot)
library(kableExtra)
library(DT)
library(formattable)
library(gridExtra)
library(grid)
library(Matrix)

library(testthat)

library(Seurat)
library(OOSAP)

library(org.Hs.eg.db)
library(AnnotationHub)
library(ReactomePA)
library(clusterProfiler)
library(fgsea) 

library(biomaRt)
library(STRINGdb) 
#library(RDAVIDWebService) 
library(msigdbr)


#########################################################################
wd <- "."
setwd(wd)

source('fxs.R', local = TRUE)
source('ui.R')
source('server.R')

#########################################################################
shinyApp(ui = ui, server = server)

