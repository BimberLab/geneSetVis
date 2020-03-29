runDGN <- function(DEtable, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable$gene), DEtable$gene), ]
  }
  
  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable$gene), 'ENTREZID', 'SYMBOL')
    
    dgn <- DOSE::enrichDGN(entrezIDs, OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

dgnModule <- function(session, input, output, envir, appDiskCache) {
  dgnResults <- reactiveValues(
    results = NULL
  )
  
  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(envir$gene_list, {
    print('resetting dgn')
    dgnResults$results <- NULL
  })
  
  observeEvent(input$rundgn_button, {
    #TODO: validate input present?
    validate(need(input$dgn_OrgDB_input != '', "Please select OrgDB..."))
    
    print('making dgn query')
    withProgress(message = 'making DGN query...', {
      cacheKey <- makeDiskCacheKey(c(envir$gene_list, input$dgn_OrgDB_input, 'dgn'))
      cacheVal <- appDiskCache$get(cacheKey)
      if (class(cacheVal) == 'key_missing') {
        print('missing cache key...')
        
        entrezIDs <- bitr(geneID = envir$gene_list$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=input$dgn_OrgDB_input)
        dgnRes <- DOSE::enrichDGN(entrezIDs$ENTREZID, readable = T)
        appDiskCache$set(key = cacheKey, value = dgnRes)
      } else {
        print('loading from cache...')
        dgnRes <- cacheVal
      }
      dgnRes@result$ID <- gsub(pattern = 'umls:', replacement = '', dgnRes@result$ID)
      dgnResults$results <- dgnRes
    })
  })
  
  renderPlotSet(
    output = output,
    key = 'dgn',
    enrichTypeResult = reactive(dgnResults$results),
    termURL = "https://www.disgenet.org/browser/0/1/0/",
    datasetName = 'dgn'
  )
  
  output$dgn_map_stats <- renderText({
    validate(need(!is.null(dgnResults$results), "No mapped genes."))
    num_genes_mapped <- str_split(noquote(dgnResults$results@result$GeneRatio[1]), '/')[[1]][2]
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$gene_list$gene), ' genes were mapped.')
    )
  })
  
  observeEvent(input$dgn_resource_info, {
    dgn_resource_info <- list(
      title = "DGN Resource info",
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
      dgn_resource_info[["text"]],
      title = dgn_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
}