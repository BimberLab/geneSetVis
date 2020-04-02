
# y <- enrichr(genes = as.vector(sample_data$gene), databases = dbs$libraryName[15])
# y[['GO_Molecular_Function_2015']]


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
              'Proteins' = Overlap,
              'p-Value' = P.value,
              'p-Value (adj.)' = Adjusted.P.value,
              'Genes in Term' = Genes
            )
          
          makeTermsTable(
            table = table,
            genesDelim = ',',
            termURL = NULL,
            caption = NULL,
            includeColumns = c('Term Description', 'Proteins', 'p-Value (adj.)', 'p-Value', 'Genes in Term')
          )
          
        })
        
      })
    })
  })
  
  
  
  renderPlotSet(
    output = output,
    key = 'enrichr',
    enrichTypeResult = reactive(enrichrResults$results),
    termURL = "https://www.disgenet.org/browser/0/1/0/",
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
      title = "enrichr Resource info",
      text = HTML(
        '<b>DisGeNET Browser</b><br>
				DisGeNET integrates data from expert curated repositories, GWAS catalogues, animal models and the scientific literature. DisGeNET data are homogeneously annotated with controlled vocabularies and community-driven ontologies. Additionally, several original metrics are provided to assist the prioritization of genotype–phenotype relationships.
				<p>
				<li><a href=https://www.disgenet.org/home/
				title="DisGeNET Browser website"
				target="_blank"><b>DisGeNET Browser website</b></a></li>
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