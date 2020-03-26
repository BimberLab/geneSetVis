source('fxs.R', local = TRUE)
source('info.R', local = TRUE)

source('modules/stringdb.R')
source('modules/msigdb.R')
source('modules/reactome.R')

server = function(input, output, session) {
  
  options(shiny.maxRequestSize=50*1024^2) 
  preferences <- reactiveValues(use_webgl = TRUE)
  
  envir <- reactiveValues(
    gene_list = NULL
  )

  ##-------------------##
  ## Sample data
  ##-------------------##

  observeEvent(input$submit, {
    gene_list <- read.table(text = gsub("(?<=[a-z])\\s+", "\n", perl = TRUE, x = input$areaInput),
                            header = FALSE,
                            col.names = c("gene", "avg_logFC"),
                            quote = "",
                            allowEscapes = T)

    envir$gene_list <- gene_list
  })

  output$inputTable <- renderTable({
    validate(need(envir$gene_list, "Please enter the gene list above and hit submit"))

    req(input$submit)
    envir$gene_list
  })

  stringDbModule(session, input, output, envir)
  msigdbModule(session, input, output, envir)
  reactomeModule(session, input, output, envir)
}


