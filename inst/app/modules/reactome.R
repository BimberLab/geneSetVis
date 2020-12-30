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

	#NOTE: this should reset our tab whenever the input genes change
	observeEvent(list(envir$geneList), ignoreInit = F, {
		envir$reactomeRes <- NULL
		errEl <- NULL
		if (!is.null(errEl)) {shinyjs::hide(errEl)}
	})

	observeEvent(input$runreactome_button, {
	  withBusyIndicatorServer("runreactome_button", {
			validate(need(!is.null(envir$geneList) && nrow(envir$geneList) > 0, "Please enter genes into Load Data tab"))
			validate(need(input$reactome_selectGeneCol != '', "Please select gene column to use"))
	    validate(need(input$reactome_OrgDB_input != '', "Please select OrgDB..."))

			shinybusy::show_modal_spinner(text = 'Querying Reactome. This might take some time.')

			cacheKey <- makeDiskCacheKey(list(envir$geneList[[input$reactome_selectGeneCol]], input$reactome_OrgDB_input), 'reactome')
			cacheVal <- appDiskCache$get(cacheKey)
			if (class(cacheVal) == 'key_missing') {
				print('missing cache key...')

				#if (!require(input$reactome_OrgDB_input)) install.packages(input$reactome_OrgDB_input)
				envir$reactomeRes <- NULL
				fromType <- ifelse(grepl('id', input$reactome_selectGeneCol), 'ENSEMBL', 'SYMBOL')
				entrezIDs <- clusterProfiler::bitr(geneID = envir$geneList[[input$reactome_selectGeneCol]], fromType=fromType, toType="ENTREZID", OrgDb=input$reactome_OrgDB_input)
				reactomeRes <- ReactomePA::enrichPathway(entrezIDs$ENTREZID, readable = T)
				appDiskCache$set(key = cacheKey, value = reactomeRes)
			} else {
				print('loading from cache...')
				reactomeRes <- cacheVal
			}

			envir$reactomeRes <- reactomeRes
			remove_modal_progress()

			if ( is.null(envir$reactomeRes) || nrow(envir$reactomeRes) == 0 ) {
				stop('No significant enrichment found.')
			}
	  })
	})

	renderPlotSet(
	  output = output,
		key = 'reactome',
		enrichTypeResult = reactive(envir$reactomeRes),
		datasetURL = "https://reactome.org/PathwayBrowser/#/",
		datasetName = 'Reactome',
		namedGeneList = envir$namedGeneList
	)

	output$reactome_map_stats <- renderText({
	  validate(need(!is.null(envir$reactomeRes) & length(envir$reactomeRes) != 0, "No mapped genes."))
	  if (nrow(envir$reactomeRes@result) > 0) {
	    num_genes_mapped <- stringr::str_split(noquote(envir$reactomeRes@result$GeneRatio[1]), '/')[[1]][2]
	  } else {
	    num_genes_mapped <- 0
	  }
	  HTML(
	    '<b>Mapped genes</b><br>',
	    paste0(num_genes_mapped, ' out of ', length(envir$geneList[[input$reactome_selectGeneCol]]), ' genes were mapped.')
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
				<p>
				<p>
				<b>enrichplot</b><br>
				The plots are produced using <a href=https://bioconductor.org/packages/release/bioc/html/enrichplot.html target="_blank"><b>enrichplot</b></a> by
				<a href=https://github.com/YuLab-SMU/enrichplot target="_blank"><b>Yu G (2019) </b></a>
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
