runDAVID <- function(DEtable, geneCol, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable[[geneCol]]), DEtable[[geneCol]]), ]
  }
  
  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable[[geneCol]]), 'ENTREZID', 'SYMBOL')
    
    david <- clusterProfiler::enrichDAVID(entrezIDs, david.user = 'oosap@ohsu.edu', OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

davidModule <- function(session, input, output, envir, appDiskCache) {
  
  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(list(envir$gene_list), ignoreInit = F, {
    print('resetting david')
    envir$davidRes <- NULL
    errEl <- NULL
    if (!is.null(errEl)) {shinyjs::hide(errEl)}
  })
  
  observeEvent(input$rundavid_button, {
    withBusyIndicatorServer("rundavid_button", {
      Sys.sleep(1)
      #TODO: validate input present?
      validate(need(input$david_OrgDB_input != '', "Please select OrgDB..."))
      
      print('making david query')
      withProgress(message = 'making DAVID query...', {
        cacheKey <- makeDiskCacheKey(list(envir$gene_list[[input$david_selectGeneCol]], input$david_OrgDB_input), 'david')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')
          
          #if (!require(input$david_OrgDB_input)) install.packages(input$david_OrgDB_input)
          envir$davidRes <- NULL
          fromType <- ifelse(grepl('id', input$david_selectGeneCol), 'ENSEMBL', 'SYMBOL')
          entrezIDs <- bitr(geneID = envir$gene_list[[input$david_selectGeneCol]], fromType=fromType, toType="ENTREZID", OrgDb=input$david_OrgDB_input)
          davidRes <- clusterProfiler::enrichDAVID(entrezIDs$ENTREZID, david.user = input$davidUserEmail)
          
          appDiskCache$set(key = cacheKey, value = davidRes)
        } else {
          print('loading from cache...')
          davidRes <- cacheVal
        }
        
        envir$davidRes <- davidRes
        if ( is.null(envir$davidRes)|| nrow(envir$davidRes) == 0 ) {stop('No significant enrichment found.')}
        ##change result replace entrezIDs with SYMBOLS; replace with & to replace
        davidRes@result$geneID <- lapply(davidRes@result$geneID, function(x) sapply(strsplit(as.character(x), "/"), 
                                                                                    function(y) paste(entrezIDs$SYMBOL[match(y, entrezIDs$ENTREZID)], collapse='/')))
        davidRes@gene <- entrezIDs$ENTREZID
        
      })
    })
  })
  
  renderPlotSet(
    output = output,
    key = 'david',
    enrichTypeResult = reactive(envir$davidRes),
    datasetURL = "",
    datasetName = 'david'
  )
  
  output$david_map_stats <- renderText({
    validate(need(!is.null(envir$davidRes) & length(envir$davidRes) != 0, "No mapped genes."))
    if (nrow(envir$davidRes@result) > 0) {
      num_genes_mapped <- str_split(noquote(envir$davidRes@result$GeneRatio[1]), '/')[[1]][2]
    } else {
      num_genes_mapped <- 0
    }
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$gene_list[[input$david_selectGeneCol]]), ' genes were mapped.')
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