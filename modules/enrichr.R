
enrichrModule <- function(session, input, output, envir, appDiskCache) {
  enrichrResults <- reactiveValues(
    results = NULL
  )
  
  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(list(envir$gene_list), {
    print('resetting enrichr')
    enrichrResults$results <- NULL
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
          
          #if (!require(input$enrichr_OrgDB_input)) install.packages(input$enrichr_OrgDB_input)
          enrichrResults$results <- NULL
          fromType <- ifelse(grepl('id', input$enrichr_selectGeneCol), 'ENSEMBL', 'SYMBOL')
          enrichrRes <- enrichr(genes = as.vector(envir$gene_list[[input$enrichr_selectGeneCol]]), databases = input$enrichr_db)
          
          appDiskCache$set(key = cacheKey, value = enrichrRes)
        } else {
          print('loading from cache...')
          enrichrRes <- cacheVal
        }
        
        enrichrResults$results <- enrichrRes
        if ( is.null(enrichrResults$results) || length(enrichrResults$results) == 0 ) {stop('No significant enrichment found.')}
        
        #enrichrRes@result$ID <- gsub(pattern = 'umls:', replacement = '', enrichrRes@result$ID)
        #rownames(enrichrRes@result) <- enrichrRes@result$ID
        
        output$enrichrResults_selected_ui<- renderUI({
          req(enrichrResults$results)
          selectInput(
            inputId = 'enrichrResults_selected',
            label = 'Select query result to view:',
            choices = names(enrichrResults$results)
          )
        })
        
        output$enrichrResults_selected_table <- renderDataTable(server = FALSE, {
          validate(need(!is.null(enrichrResults$results[[input$enrichrResults_selected]]) & length(enrichrResults$results[[input$enrichrResults_selected]]) != 0, ""))
          table <- enrichrResults$results[[input$enrichrResults_selected]] %>%
            dplyr::rename(
              'Term Description' = Term,
              'p-Value' = P.value,
              'p-Value (adj.)' = Adjusted.P.value,
              'Genes in Term' = Genes
            )
          
          makeTermsTable(
            table = table,
            genesDelim = ';',
            datasetURL = NULL,
            caption = NULL,
            includeColumns = c('Term Description', 'p-Value (adj.)', 'p-Value', 'Genes in Term', 'Overlap', 'Odds.Ratio', 'Combined.Score')
          )
          
        })
        
      })
    })
  })
  
  renderPlotSet(
    output = output,
    key = 'enrichr',
    enrichTypeResult = reactive(enrichrResults$results),
    datasetURL = "https://www.disgenet.org/browser/0/1/0/",
    datasetName = 'enrichr'
  )
  
  output$enrichr_map_stats <- renderText({
    validate(need(!is.null(enrichrResults$results) & length(enrichrResults$results) != 0, "No mapped genes."))
    num_genes_mapped <- str_split(noquote(enrichrResults$results@result$GeneRatio[1]), '/')[[1]][2]
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
				\"Enrichr contains 164 gene-set libraries where some libraries are borrowed from other tools while many other libraries are newly created and only available in Enrichr. The gene-set libraries provided by Enrichr are divided into six categories: transcription, pathways, ontologies, diseases/drugs, cell types and miscellaneous.\"
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