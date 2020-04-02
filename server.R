
server = function(input, output, session) {
  
  options(shiny.maxRequestSize=50*1024^2) 
  
  
  #shinyOptions(cache = diskCache("./cache"))
  appDiskCache <- diskCache("./cache")
  
  
  envir <- reactiveValues(
    gene_list = NULL
  )


  #---------------------------
  observeEvent(input$submit, {
    if (input$inputType == "Gene only") {
      print('true')
      gene_list <- read.table(text = gsub(",", "\n", perl = TRUE, x = input$areaInput),
                                header = FALSE,
                                col.names = c("gene"),
                                quote = "",
                                allowEscapes = T)
      gene_list <- data.frame(gene = gene_list, avg_logFC = NA)
    } else {
      gene_list <- read.table(text = gsub("(?<=[a-z])\\s+", "\n", perl = TRUE, x = input$areaInput),
                              header = FALSE,
                              col.names = c("gene", "avg_logFC"),
                              quote = "",
                              allowEscapes = T)
    }
    
    
    if (input$checkGeneIdTranslate == T) {
      withProgress(message = 'Translating genes..', {
        print(paste0('gene translate: ', input$checkGeneIdTranslate))
        print(paste0('gene id type: ', input$geneIdType))
        cacheKey <- makeDiskCacheKey(list(gene_list, input$checkGeneIdTranslate, input$geneIdType), 'genelist')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')
          
          if (input$geneIdType == 'Symbol') {
            ensemblIds <- NULL
            geneSymbols <- gene_list$gene
          } else {
            ensemblIds <- gene_list$gene
            geneSymbols <- NULL
          }
          
          # gene_list_tr <- TranslateGeneNames(ensemblIds = ensemblIds, geneSymbols = geneSymbols, davidEmail = 'oosap@ohsu.edu', 
          #                                     useEnsembl = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
          #                                     useSTRINGdb = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
          #                                     useDAVID = F)
          
          gene_list_tr <- TranslateToEnsembl(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
          gene_list_tr <- gene_list_tr[, !(colnames(gene_list_tr) %in% c('EnsemblId', 'GeneSymbol'))]
        
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
  
  observe({
    geneColnames <- envir$gene_list
    geneColnames['avg_logFC'] <- NULL
    updateSelectInput(session, "stringdb_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "msigdb_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "reactome_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "david_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "dose_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "ncg_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "dgn_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "enrichr_selectGeneCol", choices = colnames(geneColnames))
    })

  stringDbModule(session, input, output, envir, appDiskCache)
  msigdbModule(session, input, output, envir, appDiskCache)
  reactomeModule(session, input, output, envir, appDiskCache)
  davidModule(session, input, output, envir, appDiskCache)
  doseModule(session, input, output, envir, appDiskCache)
  ncgModule(session, input, output, envir, appDiskCache)
  dgnModule(session, input, output, envir, appDiskCache)
  enrichrModule(session, input, output, envir, appDiskCache)
}


