runSTRINGdb <- function(DEtable, geneCol, maxHitsToPlot = 200, refSpeciesNum = 9606, scoreThreshold = 0, cachedir = "") {
  #DEtable <- geneList
  #geneCol <- 'gene'
  string_db <- STRINGdb::STRINGdb$new(
      version = '10',
      species = refSpeciesNum,
      score_threshold = scoreThreshold,
      input_directory = cachedir
    )

	##dedup table to remove multiple tests
	if (!is.null(DEtable$test)){
		DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
		DEtable <- DEtable[match(unique(DEtable$gene), DEtable$gene), ]
	}

	return_list = list()
	tryCatch({
		clusterTable <- DEtable

		if (nrow(DEtable) > 0) {
			cluster.map <- string_db$map(as.data.frame(clusterTable), geneCol, removeUnmappedRows = FALSE)
			hits <- cluster.map$STRING_id
			# STRING does not support lists with more than 400 genes
			if (length(hits) > 400) {
			  hits <- hits[1:400]
			  print('STRING will only map the first 400 of your genes.')
			}

			if ( sum(!is.na(hits)) == 0 ) {stop('No mapped genes.')}

			max_hits_to_plot <- cluster.map$STRING_id[1:maxHitsToPlot]

			enrichmentGO <- string_db$get_enrichment(hits,
			                           category = 'Process',
			                           methodMT = 'fdr',
			                           iea = TRUE)

			enrichmentKEGG <- string_db$get_enrichment(hits,
			                           category = 'KEGG',
			                           methodMT = 'fdr',
			                           iea = TRUE)


			hit_term_proteins <- string_db$get_term_proteins(enrichmentGO$term_id, hits)
			hit_term_genes <- hit_term_proteins %>%
				dplyr::select(term_id, preferred_name) %>%
				dplyr::group_by(term_id) %>%
				dplyr::summarize('hit_term_genes' = paste0(preferred_name, collapse = ','))

			enrichmentGO <- merge(hit_term_genes, enrichmentGO)


			hit_term_proteins <- string_db$get_term_proteins(enrichmentKEGG$term_id, hits)
			hit_term_genes <- hit_term_proteins %>%
				dplyr::select(term_id, preferred_name) %>%
				dplyr::group_by(term_id) %>%
				dplyr::summarize('hit_term_genes' = paste0(preferred_name, collapse = ','))

			enrichmentKEGG <- merge(hit_term_genes, enrichmentKEGG)


			#string_db$get_png(max_hits_to_plot, file = paste('network.png', sep = ''))


			#network <- string_db$plot_network(max_hits_to_plot)

			#______
			##payload mechanism for upregulated vs downregulated genes:
			##adds a color column for up vs downregulated genes
			cluster.color <- string_db$add_diff_exp_color(cluster.map, logFcColStr = 'avg_logFC')
			# post payload information to the STRING server
			payload_id <- string_db$post_payload(cluster.color$STRING_id, colors = cluster.color$color)
			#string_db$plot_network(hits, payload_id = payload_id)

			##clustering/community algorithms: ”fastgreedy”, ”walktrap”, ”spinglass”, ”edge.betweenness”.
			# networkClustersList <-
			#   string_db$get_clusters(max_hits_to_plot, algorithm = 'fastgreedy')
			# par(mfrow = c(2, 2))
			# for (j in seq(1:length(networkClustersList))) {
			#   string_db$plot_network(networkClustersList[[j]], payload_id = payload_id)
			# }

			link <- string_db$get_link(hits[!is.na(hits)])

			graph <- string_db$get_graph()

			addSubset = paste('hits', sep = '')
			return_list[[addSubset]] <- hits

			#addSubset = paste('network', sep = '')
			#return_list[[addSubset]] <- network

			addSubset = paste('GO', sep = '')
			return_list[[addSubset]] <- enrichmentGO

			addSubset = paste('KEGG', sep = '')
			return_list[[addSubset]] <- enrichmentKEGG

			addSubset = paste('link', sep = '')
			return_list[[addSubset]] <- link


		}

	}, error = function(e) {
		cat('\nERROR :', conditionMessage(e), '\n')
	})

	return(return_list)
}

stringdbModule <- function(session, input, output, envir, appDiskCache) {

	#NOTE: this should reset our tab whenever the input genes change
	observeEvent(list(envir$geneList), ignoreInit = F, {
		print('resetting stringdb')
	  envir$stringdbRes <- NULL
		errEl <- NULL
		if (!is.null(errEl)) {shinyjs::hide(errEl)}
	})

	stringdbSpecies <- STRINGdb::get_STRING_species(version = '10')
	observeEvent(input$runstringdb_button, {
	  withBusyIndicatorServer("runstringdb_button", {
	    Sys.sleep(1)
	    validate(need(input$stringdb_maxHitsToPlot_input != '', "Please type in maxHitsToPlot..."))
	    validate(need(input$stringdb_scoreThreshold_input != '', "Please type in scoreThreshold..."))
	    validate(need(input$stringdb_refSpecies_input != '', "Please type in refSpecies..."))
	    refSpeciesNum = stringdbSpecies$species_id[stringdbSpecies$compact_name == input$stringdb_refSpecies_input]

	    print('making StringDB query')
	    withProgress(message = 'making STRING query..', {
	      cacheKey <- makeDiskCacheKey(list(envir$geneList, input$stringdb_selectGeneCol, input$stringdb_maxHitsToPlot_input, input$stringdb_refSpecies_input, input$stringdb_scoreThreshold_input), 'stringdb')
	      cacheVal <- appDiskCache$get(cacheKey)
	      if (class(cacheVal) == 'key_missing') {
	        print('missing cache key...')

	        stringdbRes <- runSTRINGdb(
	          DEtable = envir$geneList,
	          geneCol = input$stringdb_selectGeneCol,
	          maxHitsToPlot = input$stringdb_maxHitsToPlot_input,
	          refSpeciesNum = refSpeciesNum,
	          scoreThreshold = input$stringdb_scoreThreshold_input,
	          cachedir = envir$cachedir
	        )

	        appDiskCache$set(key = cacheKey, value = stringdbRes)

	      } else {
	        print('loading from cache...')
	        stringdbRes <- cacheVal
	      }

	      envir$stringdbRes <- stringdbRes
	      if (is.null(envir$stringdbRes) | length(envir$stringdbRes) == 0) {stop('No significant enrichment found.')}

	    })
	  })
	})

	output$stringdb_map_stats <- renderText({
	  validate(need(!is.null(envir$stringdbRes) & length(envir$stringdbRes) != 0, "No mapped genes."))
	  num_genes_mapped <- sum(!is.na(envir$stringdbRes[['hits']]))
	  HTML(
	    '<b>Mapped genes</b><br>',
	    paste0(num_genes_mapped, ' out of ', length(envir$geneList[[input$stringdb_selectGeneCol]]), ' genes were mapped.'),
	    '<p>',
	    hyperlink_text(href_base = envir$stringdbRes[['link']], link_text = 'View mapped genes on string-db website', href_cont = NULL)
	  )
	})


	output$stringdb_network <- renderPlot({
	  validate(need(!is.null(envir$stringdbRes), "Please Run STRINGdb on input..."))
	  toSubset <- paste('network', sep = '')
	  envir$stringdbRes[[toSubset]]
	})

	output$stringdb_network_png <- renderImage(deleteFile = F, {
		validate(need(!is.null(envir$stringdbRes), "Please Run STRINGdb on input..."))
	  cacheKey <- makeDiskCacheKey(list(envir$geneList, input$stringdb_selectGeneCol, input$stringdb_maxHitsToPlot_input, input$stringdb_refSpecies_input, input$stringdb_scoreThreshold_input), 'stringdbpng')
	  cacheVal <- appDiskCache$get(cacheKey)
	  if (class(cacheVal) == 'key_missing') {
	    print('missing cache key...')

	    png_file <- paste('network', '.png', sep = '')
	    plot_png <- list(src = paste('', png_file, sep = ''), height='100%', width='100%')
	    #plot_png
	    appDiskCache$set(key = cacheKey, value = plot_png)

	  } else {
	    print('loading from cache...')
	    plot_png <- cacheVal
	  }

	  plot_png

	})

	# TODO: download entire dataset
	# server = FALSE
	output$stringdb_GO <- DT::renderDataTable({
	  validate(need(!is.null(envir$stringdbRes$GO) & length(envir$stringdbRes$GO) != 0, ""))
	  toSubset <- paste('GO', sep = '')
	  table <- envir$stringdbRes[[toSubset]] %>%
	    dplyr::rename(
	      'Term Description' = term_description,
	      'Term ID' = term_id,
	      'Proteins' = proteins,
	      'Hits' = hits,
	      'p-Value' = pvalue,
	      'p-Value (adj.)' = pvalue_fdr,
	      'Genes in Term' = hit_term_genes
	    )

	  makeTermsTable(
	    table = table,
	    genesDelim = ',',
	    datasetURL = "https://www.ebi.ac.uk/QuickGO/term/",
	    caption = NULL,
	    includeColumns = c('Term Description', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term')
	  )

	})

	# TODO: download entire dataset
	# server = FALSE
	output$stringdb_KEGG <- DT::renderDataTable({
		validate(need(!is.null(envir$stringdbRes$KEGG) & length(envir$stringdbRes$KEGG) != 0, ""))
	  toSubset <- paste('KEGG', sep = '')
	  table <- envir$stringdbRes[[toSubset]] %>%
	    dplyr::rename(
	      'Term Description' = term_description,
	      'Term ID' = term_id,
	      'Proteins' = proteins,
	      'Hits' = hits,
	      'p-Value' = pvalue,
	      'p-Value (adj.)' = pvalue_fdr,
	      'Genes in Term' = hit_term_genes
	    )

	  makeTermsTable(
	    table = table,
	    genesDelim = ',',
	    datasetURL = "https://www.genome.jp/dbget-bin/www_bget?map",
	    caption = NULL,
	    includeColumns = c('Term Description', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term')
	  )
	})

	observeEvent(input$stringdb_resource_info, {
		stringdb_resource_info <- list(
			title = "STRING Resource info",
			text = HTML(
				'<b>STRING</b><br>
					STRING is a database of known and predicted protein-protein interactions. \n
					The interactions include direct (physical) and indirect (functional) associations; they stem from computational prediction, \n
					from knowledge transfer between organisms, and from interactions aggregated from other (primary) databases.
					<p>
				<li><a href=https://string-db.org/ title="Official string-db website" target="_blank"><b>Official string-db website</b></a></li>
				'
			)
		)

		showModal(modalDialog(
			stringdb_resource_info[["text"]],
			title = stringdb_resource_info[["title"]],
			easyClose = TRUE,
			footer = NULL
		))
	})
}
