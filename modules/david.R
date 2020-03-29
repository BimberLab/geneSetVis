runDAVID <- function(DEtable, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable$gene), DEtable$gene), ]
  }
  
  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable$gene), 'ENTREZID', 'SYMBOL')
    
    david <- clusterProfiler::enrichDAVID(entrezIDs, david.user = 'oosap@ohsu.edu', OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

davidModule <- function(session, input, output, envir, appDiskCache) {
  davidResults <- reactiveValues(
    results = NULL
  )
  
  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(envir$gene_list, {
    print('resetting david')
    davidResults$results <- NULL
  })
  
  observeEvent(input$rundavid_button, {
    #TODO: validate input present?
    validate(need(input$david_OrgDB_input != '', "Please select OrgDB..."))
    
    print('making david query')
    withProgress(message = 'making DAVID query...', {
      cacheKey <- makeDiskCacheKey(c(envir$gene_list, input$david_OrgDB_input, 'david'))
      cacheVal <- appDiskCache$get(cacheKey)
      if (class(cacheVal) == 'key_missing') {
        print('missing cache key...')
        
        entrezIDs <- bitr(geneID = envir$gene_list$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=input$david_OrgDB_input)
        davidRes <- clusterProfiler::enrichDAVID(entrezIDs$ENTREZID, david.user = 'oosap@ohsu.edu', readable = T)
        appDiskCache$set(key = cacheKey, value = davidRes)
      } else {
        print('loading from cache...')
        davidRes <- cacheVal
      }
      davidResults$results <- davidRes
    })
  })
  
  renderPlotSet(
    output = output,
    key = 'david',
    enrichTypeResult = reactive(davidResults$results),
    termURL = "",
    datasetName = 'david'
  )
  
  output$david_map_stats <- renderText({
    validate(need(!is.null(davidResults$results), "No mapped genes."))
    num_genes_mapped <- str_split(noquote(davidResults$results@result$GeneRatio[1]), '/')[[1]][2]
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$gene_list$gene), ' genes were mapped.')
    )
  })
  
  observeEvent(input$david_resource_info, {
    david_resource_info <- list(
      title = "DAVID Resource info",
      text = HTML(
        '<b>DAVID Browser</b><br>
				The DAVID has been developed as a standardized ontology for human disease with the purpose of providing the biomedical community with consistent, reusable and sustainable descriptions of human disease terms, phenotype characteristics and related medical vocabulary disease concepts through collaborative efforts of researchers at Northwestern University, Center for Genetic Medicine and the University of Maryland School of Medicine, Institute for Genome Sciences.
				<p>
				<li><a href=https://david.ncifcrf.gov/home.jsp
				title="DAVID Browser website"
				target="_blank"><b>DAVID Browser website</b></a></li>
				'
      )
    )
    
    showModal(modalDialog(
      david_resource_info[["text"]],
      title = david_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
}