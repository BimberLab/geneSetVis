
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
        cacheKey <- makeDiskCacheKey(list(gene_list, input$select_gene_conversion), 'genelist')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')
          gene_list_tr <- TranslateGeneNames(ensemblIds = gene_list$gene, geneSymbols = gene_list$gene, davidEmail = 'oosap.ohsu.edu',
                                              useEnsembl = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
                                              useSTRINGdb = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
                                              useDAVID = ifelse('useDAVID' %in% input$select_gene_conversion, T, F))
        
          gene_list <- dplyr::bind_cols(gene_list[1:nrow(gene_list),], gene_list_tr[1:nrow(gene_list_tr),])
          appDiskCache$set(key = cacheKey, value = gene_list)
        } else {
          gene_list <- cacheVal
        }
      })
    }
    
    envir$gene_list <- gene_list
  })

  output$inputTable <- DT::renderDataTable(server = FALSE, {
    validate(need(envir$gene_list, "Please enter the gene list and hit submit"))

    req(input$submit)
    envir$gene_list %>% 
    DT::datatable(
      extensions = c('Buttons'),
      options = list(
        #autoWidth = TRUE, autoHeight = TRUE, scrollX = TRUE, scrollY = TRUE,
        dom = 'lBfrtip',
        lengthMenu = list(c(15, 30, 50, -1), c('15', '30', '50', 'All')),
        pageLength = 10,
        scrollX = TRUE,
        buttons = list(
          #list(extend = "collection", text = 'Show All', action = DT::JS("function ( e, dt, node, config ) { dt.page.len(-1); dt.ajax.reload(); }")),
          list(extend = 'collection', text = 'Download/Copy', buttons = c('copy', 'csv', 'excel') )
        )
      )
    ) 
  })

  stringDbModule(session, input, output, envir, appDiskCache)
  msigdbModule(session, input, output, envir, appDiskCache)
  reactomeModule(session, input, output, envir, appDiskCache)
  davidModule(session, input, output, envir, appDiskCache)
  doseModule(session, input, output, envir, appDiskCache)
  ncgModule(session, input, output, envir, appDiskCache)
  dgnModule(session, input, output, envir, appDiskCache)
}


