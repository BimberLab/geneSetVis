
server = function(input, output, session) {
  
  options(shiny.maxRequestSize=50*1024^2) 
  
  
  #shinyOptions(cache = diskCache("./cache"))
  appDiskCache <- diskCache("./cache")
  
  
  envir <- reactiveValues(
    gene_list = NULL
  )


  #---------------------------
  observeEvent(input$submit, {
    gene_list <- read.table(text = gsub("(?<=[a-z])\\s+", "\n", perl = TRUE, x = input$areaInput),
                            header = FALSE,
                            col.names = c("gene", "avg_logFC"),
                            quote = "",
                            allowEscapes = T)
    
    if (!is.null(input$select_gene_conversion)) {
      withProgress(message = 'Translating genes..', {
        print(paste0('print gene conversion options:', input$select_gene_conversion))
        gene_list_tr <- TranslateGeneNames(ensemblIds = gene_list$gene, geneSymbols = gene_list$gene, davidEmail = 'oosap.ohsu.edu', 
                                            useEnsembl = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F), 
                                            useSTRINGdb = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F), 
                                            useDAVID = ifelse('useDAVID' %in% input$select_gene_conversion, T, F))
        
        gene_list <- dplyr::bind_cols(gene_list[1:nrow(gene_list),], gene_list_tr[1:nrow(gene_list_tr),])
      })
      
    }
    
    envir$gene_list <- gene_list
  })

  output$inputTable <- DT::renderDataTable({
    validate(need(envir$gene_list, "Please enter the gene list and hit submit"))

    req(input$submit)
    envir$gene_list %>% 
    DT::datatable(
      extensions = c('Buttons'),
      options = list(
        #autoWidth = TRUE, autoHeight = TRUE, scrollX = TRUE, scrollY = TRUE,
        dom = 'Bfrtip',
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 10,
        scrollX = TRUE,
        buttons = list(
          list(extend = "collection", text = 'Show All', action = DT::JS("function ( e, dt, node, config ) { dt.page.len(-1); dt.ajax.reload(); }")),
          list(extend = 'collection', text = 'Download/Copy', buttons = c('copy', 'csv', 'excel') )
        )
      )
    ) 
  })

  stringDbModule(session, input, output, envir, appDiskCache)
  msigdbModule(session, input, output, envir, appDiskCache)
  reactomeModule(session, input, output, envir, appDiskCache)
}


