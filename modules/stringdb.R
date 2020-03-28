runSTRINGdb <- function(DEtable, maxHitsToPlot = 200, refSpeciesNum = 9606, scoreThreshold = 0) {
	string_db <-
	STRINGdb::STRINGdb$new(
	version = '10',
	species = refSpeciesNum,
	score_threshold = scoreThreshold,
	input_directory = ''
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
			cluster.map <-
			string_db$map(clusterTable, 'gene', removeUnmappedRows = FALSE)
			hits <- cluster.map$STRING_id
			if ( sum(!is.na(hits)) == 0 ) {stop('No mapped genes.')}

			max_hits_to_plot <- cluster.map$STRING_id[1:maxHitsToPlot]

			enrichmentGO <-
			string_db$get_enrichment(hits,
			category = 'Process',
			methodMT = 'fdr',
			iea = TRUE)

			enrichmentKEGG <-
			string_db$get_enrichment(hits,
			category = 'KEGG',
			methodMT = 'fdr',
			iea = TRUE)


			hit_term_proteins <-
			string_db$get_term_proteins(enrichmentGO$term_id, hits)
			hit_term_genes <- hit_term_proteins %>%
				dplyr::select(term_id, preferred_name) %>%
				dplyr::group_by(term_id) %>%
				dplyr::summarize('hit_term_genes' = paste0(preferred_name, collapse = ','))

			enrichmentGO <- merge(hit_term_genes, enrichmentGO)


			hit_term_proteins <-
			string_db$get_term_proteins(enrichmentKEGG$term_id, hits)
			hit_term_genes <- hit_term_proteins %>%
				dplyr::select(term_id, preferred_name) %>%
				dplyr::group_by(term_id) %>%
				dplyr::summarize('hit_term_genes' = paste0(preferred_name, collapse = ','))

			enrichmentKEGG <- merge(hit_term_genes, enrichmentKEGG)


			string_db$get_png(max_hits_to_plot, file = paste('network.png', sep = ''))


			#network <- string_db$plot_network(max_hits_to_plot)

			#______
			##payload mechanism for upregulated vs downregulated genes:
			##adds a color column for up vs downregulated genes
			cluster.color <-
			string_db$add_diff_exp_color(cluster.map, logFcColStr = 'avg_logFC')
			# post payload information to the STRING server
			payload_id <-
			string_db$post_payload(cluster.color$STRING_id, colors = cluster.color$color)
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

stringDbModule <- function(session, input, output, envir, appDiskCache) {
	stringResults <- reactiveValues(
		results = NULL
	)

	#NOTE: this should reset our tab whenever the input genes change
	observeEvent(envir$gene_list, {
		print('resetting stringdb')
		stringResults$results <- NULL
	})

	stringdbSpecies <- STRINGdb::get_STRING_species(version = '10')
	observeEvent(input$runstringdb_button, {
		validate(need(input$stringdb_maxHitsToPlot_input != '', "Please type in maxHitsToPlot..."))
		validate(need(input$stringdb_scoreThreshold_input != '', "Please type in scoreThreshold..."))
		validate(need(input$stringdb_refSpecies_input != '', "Please type in refSpecies..."))
		refSpeciesNum = stringdbSpecies$species_id[stringdbSpecies$compact_name == input$stringdb_refSpecies_input]

		print('making StringDB query')
		withProgress(message = 'making STRING query..', {
			# saveFile <- paste0('SavedRuns/', 'running', '_string_result', '.rds', sep = '')
			# if (file.exists(saveFile)) {
			# 	stringRes <- readRDS(saveFile)
		  
		  cacheKey <- makeDiskCacheKey(c(envir$gene_list, input$stringdb_maxHitsToPlot_input, refSpeciesNum, input$stringdb_scoreThreshold_input))
		  cacheVal <- appDiskCache$get(cacheKey)
		  if (class(cacheVal) == 'key_missing') {
		    print('missing cache key...')
		    
		    stringRes <- runSTRINGdb(
		      DEtable = envir$gene_list, 
		      maxHitsToPlot = input$stringdb_maxHitsToPlot_input, 
		      refSpeciesNum = refSpeciesNum, 
		      scoreThreshold = input$stringdb_scoreThreshold_input
		    )
		    
		    appDiskCache$set(key = cacheKey, value = stringRes)
		  
			} else {
				# stringRes <- runSTRINGdb(DEtable = envir$gene_list,
				# maxHitsToPlot = input$stringdb_maxHitsToPlot_input,
				# refSpeciesNum = refSpeciesNum,
				# scoreThreshold = input$stringdb_scoreThreshold_input)
				# saveRDS(stringRes, file = saveFile)
			  print('loading from cache...')
			  stringRes <- cacheVal
			}
			stringResults$results <- stringRes
		})
	})
	
	output$string_map_stats <- renderText({
	  validate(need(!is.null(stringResults$results), "No mapped genes."))
	  num_genes_mapped <- sum(!is.na(stringResults$results[['hits']]))
	  HTML(
	    '<b>Mapped genes</b><br>',
	    paste0(num_genes_mapped, ' out of ', length(envir$gene_list$gene), ' genes were mapped.'),
	    '<p>',
	    hyperlink_text(url = stringResults$results[['link']], text = 'View mapped genes on string-db website', hide = NULL)
	  )
	})

	output$stringdb_network <- renderPlot({
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
	  toSubset <- paste('network', sep = '')
	  stringResults$results[[toSubset]]
	})

	output$stringdb_network_png <- renderImage(deleteFile = T, {
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
		png_file <- paste('network', '.png', sep = '')
		list(src = paste('', png_file, sep = ''), height='100%', width='100%')
	})

	# TODO: download entire dataset
	output$stringdb_GO <- renderDataTable({
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
	  validate(need(!is.null(stringResults$results), "No mapped genes."))
	  toSubset <- paste('GO', sep = '')
	  table <- stringResults$results[[toSubset]] %>% 
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
	    termURL = "https://www.ebi.ac.uk/QuickGO/term/",
	    caption = NULL,
	    includeColumns = c('Term Description', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term')
	  )
	  
	})

	# TODO: download entire dataset
	output$stringdb_KEGG <- renderDataTable({
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
	  toSubset <- paste('KEGG', sep = '')
	  table <- stringResults$results[[toSubset]] %>% 
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
	    termURL = "https://www.genome.jp/dbget-bin/www_bget?map",
	    caption = NULL,
	    includeColumns = c('Term Description', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term')
	  )
	})

	observeEvent(input$stringdb_resource_info, {
		showModal(
		modalDialog(
		stringdb_resource_info[["text"]],
		title = stringdb_resource_info[["title"]],
		easyClose = TRUE,
		footer = NULL
		)
		)
	})
}