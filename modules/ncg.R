runNCG <- function(DEtable, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable$gene), DEtable$gene), ]
  }
  
  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable$gene), 'ENTREZID', 'SYMBOL')
    
    ncg <- DOSE::enrichNCG(entrezIDs, OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ncgModule <- function(session, input, output, envir, appDiskCache) {
  ncgResults <- reactiveValues(
    results = NULL
  )
  
  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(envir$gene_list, {
    print('resetting ncg')
    ncgResults$results <- NULL
  })
  
  observeEvent(input$runncg_button, {
    #TODO: validate input present?
    validate(need(input$ncg_OrgDB_input != '', "Please select OrgDB..."))
    
    print('making ncg query')
    withProgress(message = 'making NCG query...', {
      cacheKey <- makeDiskCacheKey(c(envir$gene_list, input$ncg_OrgDB_input, 'ncg'))
      cacheVal <- appDiskCache$get(cacheKey)
      if (class(cacheVal) == 'key_missing') {
        print('missing cache key...')
        
        entrezIDs <- bitr(geneID = envir$gene_list$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=input$ncg_OrgDB_input)
        ncgRes <- DOSE::enrichNCG(entrezIDs$ENTREZID, readable = T)
        appDiskCache$set(key = cacheKey, value = ncgRes)
      } else {
        print('loading from cache...')
        ncgRes <- cacheVal
      }
      ncgResults$results <- ncgRes
    })
  })
  
  renderPlotSet(
    output = output,
    key = 'ncg',
    enrichTypeResult = reactive(ncgResults$results),
    termURL = "",
    datasetName = 'ncg'
  )
  
  output$ncg_map_stats <- renderText({
    validate(need(!is.null(ncgResults$results), "No mapped genes."))
    if (nrow(ncgResults$results@result) > 0) {
      num_genes_mapped <- str_split(noquote(ncgResults$results@result$GeneRatio[1]), '/')[[1]][2]
    } else {
      num_genes_mapped <- 0
    }
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$gene_list$gene), ' genes were mapped.')
    )
  })
  
  observeEvent(input$ncg_resource_info, {
    ncg_resource_info <- list(
      title = "NCG Resource info",
      text = HTML(
        '<b>Network of Cancer Genes</b><br>
				NCG contains information on duplicability, evolution, protein-protein and microRNA-gene interaction, function, expression and essentiality of 2,372 cancer genes from 273 manually curated publications 
				<p>
				<li><a href=http://ncg.kcl.ac.uk/index.php
				title="Network of Cancer Genes Browser website"
				target="_blank"><b>Network of Cancer Genes Browser website</b></a></li>
				'
      )
    )
    
    showModal(modalDialog(
      ncg_resource_info[["text"]],
      title = ncg_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
}