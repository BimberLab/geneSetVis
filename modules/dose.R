runDOSE <- function(DEtable, geneCol, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable$gene), DEtable[[geneCol]]), ]
  }
  
  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable[[geneCol]]), 'ENTREZID', 'SYMBOL')
    
    do <- enrichDO(entrezIDs, OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

doseModule <- function(session, input, output, envir, appDiskCache) {
  doseResults <- reactiveValues(
    results = NULL
  )
  
  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(list(envir$gene_list), {
    print('resetting dose')
    doseResults$results <- NULL
  })
  
  observeEvent(input$rundose_button, {
    withBusyIndicatorServer("rundose_button", {
      Sys.sleep(1)
      #TODO: validate input present?
      validate(need(input$dose_OrgDB_input != '', "Please select OrgDB..."))
      
      print('making dose query')
      withProgress(message = 'making DOSE query...', {
        cacheKey <- makeDiskCacheKey(list(envir$gene_list[[input$dose_selectGeneCol]], input$dose_OrgDB_input), 'dose')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')
          
          #if (!require(input$dose_OrgDB_input)) install.packages(input$dose_OrgDB_input)
          doseResults$results <- NULL
          fromType <- ifelse(grepl('id', input$dose_selectGeneCol), 'ENSEMBL', 'SYMBOL')
          entrezIDs <- bitr(geneID = envir$gene_list[[input$dose_selectGeneCol]], fromType=fromType, toType="ENTREZID", OrgDb=input$dose_OrgDB_input)
          doseRes <- DOSE::enrichDO(entrezIDs$ENTREZID, readable = T)
          appDiskCache$set(key = cacheKey, value = doseRes)
        } else {
          print('loading from cache...')
          doseRes <- cacheVal
        }
        
        doseResults$results <- doseRes
        if ( is.null(doseResults$results) || nrow(doseResults$results) == 0 ) {stop('No significant enrichment found.')}
        doseRes@result$ID <- gsub(pattern = 'DOID:', replacement = '', doseRes@result$ID)
        rownames(doseRes@result) <- doseRes@result$ID
        
        
      })
    })
  })
  
  renderPlotSet(
    output = output,
    key = 'dose',
    enrichTypeResult = reactive(doseResults$results),
    termURL = "https://www.ebi.ac.uk/ols/ontologies/doid/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FDOID_",
    datasetName = 'dose'
  )
  
  output$dose_map_stats <- renderText({
    validate(need(!is.null(doseResults$results) & length(doseResults$results) != 0, "No mapped genes."))
    if (nrow(doseResults$results@result) > 0) {
      num_genes_mapped <- str_split(noquote(doseResults$results@result$GeneRatio[1]), '/')[[1]][2]
    } else {
      num_genes_mapped <- 0
    }
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$gene_list[[input$dose_selectGeneCol]]), ' genes were mapped.')
    )
  })
  
  observeEvent(input$dose_resource_info, {
    dose_resource_info <- list(
      title = "DOSE Resource info",
      text = HTML(
        '<b>Disease Ontology Browser</b><br>
				The Disease Ontology has been developed as a standardized ontology for human disease with the purpose of providing the biomedical community with consistent, reusable and sustainable descriptions of human disease terms, phenotype characteristics and related medical vocabulary disease concepts through collaborative efforts of researchers at Northwestern University, Center for Genetic Medicine and the University of Maryland School of Medicine, Institute for Genome Sciences.
				<p>
				<li><a href=https://disease-ontology.org/
				title="Disease Ontology Browser website"
				target="_blank"><b>Disease Ontology Browser website</b></a></li>
				'
      )
    )
    
    showModal(modalDialog(
      dose_resource_info[["text"]],
      title = dose_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })
    
  
}