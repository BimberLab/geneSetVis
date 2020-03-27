library(shiny)
library(shinydashboard)
library(flexdashboard)

library(shinyMatrix)
library(shinyBS)
library(rsconnect)

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
library(rhandsontable)

library(testthat)
library(uuid)

#library(Seurat)
#library(OOSAP)

library(org.Hs.eg.db)
library(AnnotationHub)
library(ReactomePA)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(fgsea) 

library(biomaRt)
library(STRINGdb) 
#library(RDAVIDWebService)
library(msigdbr)


options(repos = BiocManager::repositories())


wd <- '.'
setwd(wd)


source('fxs.R', local = TRUE)
source('modules/stringdb.R')
source('modules/msigdb.R')
source('modules/reactome.R')
source('tabs.R')
source('info.R')
source('ui.R')
source('server.R')

