
options(repos=structure(BiocManager::repositories()))
#packrat::init(project = '.')
#warnings()
#packrat::restore(project = '.')
#packrat::snapshot()

packrat::on()

#recordTest("/Users/onwuzu/Documents/BimberLab/geneSetVis/geneSetVis/", loadTimeout=100000)

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(flexdashboard)
# library(shinytest)
# shinytest::installDependencies()

library(dplyr)
library(ggplot2)
library(ggupset)
library(tidyr)
library(stringr)
library(plotly)
library(DT)
library(formattable)
library(stats)

library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(org.Mm.eg.db)
library(AnnotationHub)
library(ReactomePA)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(fgsea) 
library(enrichR) 

library(biomaRt)
library(STRINGdb) 
library(msigdbr)


wd <- '.'
setwd(wd)


source('fxs.R', local = TRUE)
source('GeneAliasing.R', local = TRUE)
source('modules/stringdb.R', local = TRUE)
source('modules/msigdb.R', local = TRUE)
source('modules/reactome.R', local = TRUE)
source('modules/david.R', local = TRUE)
source('modules/dose.R', local = TRUE)
source('modules/ncg.R', local = TRUE)
source('modules/dgn.R', local = TRUE)
source('modules/enrichr.R', local = TRUE)
source('tabs.R', local = TRUE)
source('ui.R', local = TRUE)
source('server.R', local = TRUE)

