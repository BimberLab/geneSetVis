runReactomePA <- function(DEtable, geneCol, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable[[geneCol]]), DEtable[[genecol]]), ]
  }

  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable[[geneCol]]), 'ENTREZID', 'SYMBOL')
    
    pa <- ReactomePA::enrichPathway(entrezIDs)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

reactomeModule <- function(session, input, output, envir, appDiskCache) {
	reactomeResults <- reactiveValues(
		results = NULL
	)

	#NOTE: this should reset our tab whenever the input genes change
	observeEvent(list(envir$gene_list), {
		print('resetting reactome')
		reactomeResults$results <- NULL
	})

	observeEvent(input$runreactome_button, {
	  withBusyIndicatorServer("runreactome_button", {
	    Sys.sleep(1)
	    #TODO: validate input present?
	    validate(need(input$reactome_OrgDB_input != '', "Please select OrgDB..."))
	    
	    print('making Reactome query')
	    withProgress(message = 'making reactomePA query...', {
	      cacheKey <- makeDiskCacheKey(list(envir$gene_list[[input$reactome_selectGeneCol]], input$reactome_OrgDB_input), 'reactome')
	      cacheVal <- appDiskCache$get(cacheKey)
	      if (class(cacheVal) == 'key_missing') {
	        print('missing cache key...')
	        
	        #if (!require(input$reactome_OrgDB_input)) install.packages(input$reactome_OrgDB_input)
	        reactomeResults$results <- NULL
	        fromType <- ifelse(grepl('id', input$reactome_selectGeneCol), 'ENSEMBL', 'SYMBOL')
	        entrezIDs <- bitr(geneID = envir$gene_list[[input$reactome_selectGeneCol]], fromType=fromType, toType="ENTREZID", OrgDb=input$reactome_OrgDB_input)
	        reactomeRes <- ReactomePA::enrichPathway(entrezIDs$ENTREZID, readable = T)
	        appDiskCache$set(key = cacheKey, value = reactomeRes)
	      } else {
	        print('loading from cache...')
	        reactomeRes <- cacheVal
	      }
	      
	      reactomeResults$results <- reactomeRes
	      if ( is.null(reactomeResults$results) | nrow(reactomeResults$results) == 0 ) {stop('No significant enrichment found.')}
	      
	    })
	  })
	})

	renderPlotSet(
	  output = output,
		key = 'reactome',
		enrichTypeResult = reactive(reactomeResults$results),
		termURL = "https://reactome.org/PathwayBrowser/#/",
		datasetName = 'Reactome'
	)

	output$reactome_map_stats <- renderText({
	  validate(need(!is.null(reactomeResults$results) & length(reactomeResults$results) != 0, "No mapped genes."))
	  if (nrow(reactomeResults$results@result) > 0) {
	    num_genes_mapped <- str_split(noquote(reactomeResults$results@result$GeneRatio[1]), '/')[[1]][2]
	  } else {
	    num_genes_mapped <- 0
	  }
	  HTML(
	    '<b>Mapped genes</b><br>',
	    paste0(num_genes_mapped, ' out of ', length(envir$gene_list[[input$reactome_selectGeneCol]]), ' genes were mapped.')
	  )
	})

	observeEvent(input$reactome_resource_info, {
		reactome_resource_info <- list(
			title = "Reactome Resource info",
			text = HTML(
				'<b>Reactome Pathway Browser</b><br>
				Reactome is a free, open-source, curated and peer-reviewed pathway database. It provides bioinformatics tools for the visualization, interpretation and analysis of pathway knowledge to support basic research, genome analysis, modeling, systems biology and education.
				<p>
				<li><a href=https://reactome.org/PathwayBrowser/
				title="Reactome Pathway Browser website"
				target="_blank"><b>Reactome Pathway Browser website</b></a></li>
				'
			)
		)

		showModal(modalDialog(
			reactome_resource_info[["text"]],
			title = reactome_resource_info[["title"]],
			easyClose = TRUE,
			footer = NULL
		))
	})
}