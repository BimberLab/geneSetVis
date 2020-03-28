runReactomePA <- function(DEtable, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable$gene), DEtable$gene), ]
  }
  
  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable$gene), 'ENTREZID', 'SYMBOL')
    
    go <- enrichGO(entrezIDs, OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)
    pa <- ReactomePA::enrichPathway(entrezIDs)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

reactomeModule <- function(session, input, output, envir, appDiskCache) {
	reactomeResults <- reactiveValues(
		results = NULL
	)

	#NOTE: this should reset our tab whenever the input genes change
	observeEvent(envir$gene_list, {
		print('resetting reactome')
		reactomeResults$results <- NULL
	})

	observeEvent(input$runreactome_button, {
		#TODO: validate input present?
		validate(need(input$reactome_OrgDB_input != '', "Please select OrgDB..."))

		print('making Reactome query')
		withProgress(message = 'making reactomePA query...', {
			# saveFile <- paste0('SavedRuns/', 'running', '_reactomePA_result', '.rds', sep = '')
			# if (file.exists(saveFile)) {
			# 	reactomeResults$results <- readRDS(saveFile)
			# } else {
		  cacheKey <- makeDiskCacheKey(c(envir$gene_list, input$reactome_OrgDB_input))
		  cacheVal <- appDiskCache$get(cacheKey)
		  if (class(cacheVal) == 'key_missing') {
		    print('missing cache key...')
		    
				entrezIDs <- bitr(geneID = envir$gene_list$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=input$reactome_OrgDB_input)
				reactomePAres <- ReactomePA::enrichPathway(entrezIDs$ENTREZID, readable = T)
				#saveRDS(reactomePAres, file = saveFile)
				#reactomeResults$results <- reactomePAres
				appDiskCache$set(key = cacheKey, value = reactomePAres)
		  } else {
		    print('loading from cache...')
		    reactomePAres <- cacheVal
		  }
		  reactomeResults$results <- reactomePAres
		})
	})

	renderPlotSet(output = output,
		key = 'reactome',
		enrichTypeResult = reactive(reactomeResults$results),
		termURL = "https://reactome.org/PathwayBrowser/#/",
		datasetName = 'Reactome'
	)
	
	output$reactome_map_stats <- renderText({
	  validate(need(!is.null(reactomeResults$results), "No mapped genes."))
	  num_genes_mapped <- str_split(noquote(reactomeResults$results@result$GeneRatio[1]), '/')[[1]][2]
	  HTML(
	    '<b>Mapped genes</b><br>',
	    paste0(num_genes_mapped, ' out of ', length(envir$gene_list$gene), ' genes were mapped.')
	  )
	})
	
	observeEvent(input$reactome_resource_info, {
	  showModal(
	    modalDialog(
	      reactome_resource_info[["text"]],
	      title = reactome_resource_info[["title"]],
	      easyClose = TRUE,
	      footer = NULL
	    )
	  )
	})
}