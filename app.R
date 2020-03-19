library(shiny)
library(shinydashboard)
library(flexdashboard)

library(dplyr)
library(ggplot2)
library(ggupset)



options(repos = BiocManager::repositories())

#########################################################################
wd <- '.'
setwd(wd)

source('fxs.R', local = TRUE)
source('ui.R')
source('server.R')

#########################################################################
shinyApp(ui = ui, server = server)

