
enrichrModule <- function(session, input, output, envir, appDiskCache) {
  enrichrResults <- reactiveValues(
    results = NULL,
    asenrichResult = NULL
  )
  
  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(list(envir$gene_list), ignoreInit = T, {
    print('resetting enrichr')
    envir$enrichrRes <- NULL
    envir$enrichRes_selected <- NULL
    errEl <- NULL
    if (!is.null(errEl)) {shinyjs::hide(errEl)}
  })
  
  observeEvent(input$runenrichr_button, {
    withBusyIndicatorServer("runenrichr_button", {
      Sys.sleep(1)
      #TODO: validate input present?
      validate(need(input$enrichr_db != '', "Please select database(s)..."))
      
      print('making enrichr query')
      withProgress(message = 'making enrichr query...', {
        cacheKey <- makeDiskCacheKey(list(envir$gene_list[[input$enrichr_selectGeneCol]], input$enrichr_db), 'enrichr')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')
          
          #envir$enrichrRes <- NULL
          fromType <- ifelse(grepl('id', input$enrichr_selectGeneCol), 'ENSEMBL', 'SYMBOL')
          if (fromType != 'SYMBOL') {
            envir$enrichrRes <- NULL
            stop('MsigDB only accepts gene columns of type SYMBOL')
          }
          enrichrRes <- enrichr(genes = as.vector(envir$gene_list[[input$enrichr_selectGeneCol]]), databases = input$enrichr_db)
          #enrichrRes <- enrichr(genes = as.vector(y$ensembl_gene_id), databases = input$enrichr_db)
          
          appDiskCache$set(key = cacheKey, value = enrichrRes)
        } else {
          print('loading from cache...')
          enrichrRes <- cacheVal
        }
        
        envir$enrichrRes <- enrichrRes
        if ( is.null(envir$enrichrRes) || length(envir$enrichrRes) == 0 ) {stop('No significant enrichment found.')}
        
        #enrichrRes@result$ID <- gsub(pattern = 'umls:', replacement = '', enrichrRes@result$ID)
        #rownames(enrichrRes@result) <- enrichrRes@result$ID
        
        })
        
      })
    })
  
  observe({
    updateSelectInput(session, "enrichr_selectQuery", choices = names(envir$enrichrRes), selected = tail(names(envir$enrichrRes), 1))
  })
  
  #TODO: wait to pass to as as.enrich when searching 
  observeEvent(input$enrichr_selectQuery, ignoreInit = T, {
    # withBusyIndicatorServer("enrichr_selectQuery", {
    #   Sys.sleep(1)
      #print(envir$enrichrRes[[input$enrichr_selectQuery]])
      envir$enrichRes_selected <-
        as.enrichResult(
          pvalueCutoff = 1,
          gseResult = envir$enrichrRes[[input$enrichr_selectQuery]],
          gseGenes = envir$gene_list[[input$enrichr_selectGeneCol]],
          idCol = envir$enrichrRes[[input$enrichr_selectQuery]]$Term,
          padjCol = envir$enrichrRes[[input$enrichr_selectQuery]]$Adjusted.P.value,
          pvalCol = envir$enrichrRes[[input$enrichr_selectQuery]]$P.value,
          geneIDCol = gsub(
            pattern = ';',
            replacement = '/',
            x = envir$enrichrRes[[input$enrichr_selectQuery]]$Genes
          ),
          countCol = noquote(
            str_split_fixed(envir$enrichrRes[[input$enrichr_selectQuery]]$Overlap, "/", 2)[, 1]
          ),
          geneRatioCol = paste(
            noquote(
              str_split_fixed(envir$enrichrRes[[input$enrichr_selectQuery]]$Overlap, "/", 2)[, 1]
            ),
            '/',
            length(envir$gene_list[[input$enrichr_selectGeneCol]]),
            sep = ''
          )
        )
    #})
  })
  
  # output$enrichr_selectQuery_table <- renderDataTable(server = FALSE, {
  #   validate(need(!is.null(envir$enrichrRes) & length(envir$enrichrRes) != 0, ""))
  #   table <- envir$enrichrRes[[input$enrichr_selectQuery]] %>%
  #     dplyr::rename(
  #       'Term Description' = Term,
  #       'p-Value' = P.value,
  #       'p-Value (adj.)' = Adjusted.P.value,
  #       'Genes in Term' = Genes
  #     )
  #   
  #   makeTermsTable(
  #     table = table,
  #     genesDelim = ';',
  #     datasetURL = NULL,
  #     caption = NULL,
  #     includeColumns = c('Term Description', 'p-Value (adj.)', 'p-Value', 'Genes in Term', 'Overlap', 'Odds.Ratio', 'Combined.Score')
  #   )
  
  renderPlotSet(
    output = output,
    key = 'enrichr',
    enrichTypeResult = reactive(envir$enrichRes_selected),
    datasetURL = NULL,
    datasetName = 'enrichr'
  )
  
  output$enrichr_map_stats <- renderText({
    validate(need(!is.null(envir$enrichrRes) & length(envir$enrichrRes) != 0, "No mapped genes."))
    num_genes_mapped <- str_split(noquote(envir$enrichrRes@result$GeneRatio[1]), '/')[[1]][2]
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$gene_list[[input$enrichr_selectGeneCol]]), ' genes were mapped.')
    )
  })
  
  observeEvent(input$enrichr_resource_info, {
    enrichr_resource_info <- list(
      title = "enrichR Resource info",
      text = HTML(
        '<b>enrichR</b><br>
				\"Enrichr contains 164 gene-set libraries...some libraries are borrowed from other tools while many other libraries are newly created and only available in Enrichr. The gene-set libraries provided by Enrichr are divided into six categories: transcription, pathways, ontologies, diseases/drugs, cell types and miscellaneous.\"
				<a href=https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-128#Sec2 target="_blank"><b>Chen et al. (2013)</b></a>
				<p>
				<li><a href=https://amp.pharm.mssm.edu/Enrichr/
				title="enrichR website"
				target="_blank"><b>enrichR website</b></a></li>
				'
      )
    )
    
    showModal(modalDialog(
      enrichr_resource_info[["text"]],
      title = enrichr_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
}