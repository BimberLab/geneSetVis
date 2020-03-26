source('fxs.R', local = TRUE)

#' @title plot STRINGdb networks
#'
#' @description Takes Differential Expression table and plots STRINGdb networks
#' @param DEtable A DE table
#' @param numHits The num of mapped hits to plot
#' @param refSpeciesNum The dataset (see STRINGdb docs) to use as a reference; 9606=Human, ?=Rhesus Macaque
#' @return The PNGs of network plots
#' @keywords STRINGdb
#' @import STRINGdb
#' @export
#' @importFrom
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

			link <- string_db$get_link(hits)

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

stringDbModule <- function(session, input, output, envir) {
	stringResults <- reactiveValues(
		results = NULL
	)
	
	stringdbSpecies <- STRINGdb::get_STRING_species(version = '10')
	observeEvent(input$runstringdb_button, {
		validate(need(input$stringdb_maxHitsToPlot_input != '', "Please type in maxHitsToPlot..."))
		validate(need(input$stringdb_scoreThreshold_input != '', "Please type in scoreThreshold..."))
		validate(need(input$stringdb_refSpecies_input != '', "Please type in refSpecies..."))
		refSpeciesNum = stringdbSpecies$species_id[stringdbSpecies$compact_name == input$stringdb_refSpecies_input]

		print('making StringDB query')
		withProgress(message = 'making STRING query..', {
			saveFile <- paste0('SavedRuns/', 'running', '_string_result', '.rds', sep = '')
			if (file.exists(saveFile)) {
				stringRes <- readRDS(saveFile)
			} else {
				stringRes <- runSTRINGdb(DEtable = envir$gene_list,
				maxHitsToPlot = input$stringdb_maxHitsToPlot_input,
				refSpeciesNum = refSpeciesNum,
				scoreThreshold = input$stringdb_scoreThreshold_input)
				saveRDS(stringRes, file = saveFile)
			}
			stringResults$results <- stringRes
		})
	})

	output$string_map_stats <- renderText({
		req(stringResults$results)
		paste0(sum(!is.na(stringResults$results[['hits']])), ' out of ', length(stringResults$results[['hits']]), 'genes were mapped.')
	})

	output$num_of_mapped <- flexdashboard::renderValueBox({
		req(stringResults$results)
		shinydashboard::box(
		title = 'Number of genes mapped',
		width = 6,
		background = 'light-blue',
		sum(!is.na(stringResults$results[['hits']]))
		)
	})

	output$num_of_total_genes <- flexdashboard::renderValueBox({
		req(stringResults$results)
		extract <- paste('hits', sep = '')
		shinydashboard::box(
		title = 'Number of genes total',
		width = 6,
		background = 'light-blue',
		length(stringResults$results[[extract]])
		)
	})

	output$stringdb_network <- renderPlot({
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
		extract <- paste('network', sep = '')
		stringResults$results[[extract]]
	})

	output$stringdb_network_png <- renderImage(deleteFile = F, {
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
		png_file <- paste('network', '.png', sep = '')
		list(src = paste('', png_file, sep = ''), height='90%', width='90%')
	})

	# TODO: download entire dataset
	output$stringdb_GO <- renderDataTable({
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
		extract <- paste('GO', sep = '')
		table <- stringResults$results[[extract]] %>% dplyr::rename(
		'Term Description' = term_description,
		'Term ID' = term_id,
		'Proteins' = proteins,
		'Hits' = hits,
		'p-Value' = pvalue,
		'p-Value (adj.)' = pvalue_fdr,
		'Genes in Term' = hit_term_genes
		)

		table$'Term Description' <- hyperlink_text(text = table$'Term Description', url = "https://www.ebi.ac.uk/QuickGO/term/", hide = table$'Term ID')
		table$'Genes in Term' <- multi_hyperlink_text(labels = table$'Genes in Term', links = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=")
		table <- table %>% dplyr::select(c('Term Description', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value'))

		DT::datatable(
		table,
		filter = 'bottom',
		selection = 'multiple',
		escape = FALSE,
		autoHideNavigation = TRUE,
		rownames = FALSE,
		extensions = c('Buttons'),
		class = 'cell-border stripe',
		options = list(
		dom = 'Bfrtip',
		lengthMenu = c(15, 30, 50, 100),
		pageLength = 10,
		buttons = list(
		'colvis',
		list(
		extend = 'collection',
		text = 'Download/Copy',
		buttons = c('copy', 'csv', 'excel')
		)
		)
		)
		)
	})

	# TODO: download entire dataset
	output$stringdb_KEGG <- renderDataTable({
		validate(need(!is.null(stringResults$results), "Please Run STRINGdb on input..."))
		extract <- paste('KEGG', sep = '')
		table <- stringResults$results[[extract]] %>% dplyr::rename(
		'Term Description' = term_description,
		'Term ID' = term_id,
		'Proteins' = proteins,
		'Hits' = hits,
		'p-Value' = pvalue,
		'p-Value (adj.)' = pvalue_fdr,
		'Genes in Term' = hit_term_genes
		)

		table$'Term Description' <- hyperlink_text(text = table$'Term Description', url = "https://www.genome.jp/dbget-bin/www_bget?map", hide = table$'Term ID')

		table$'Genes in Term' <- multi_hyperlink_text(labels = table$'Genes in Term', links = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=")

		table <- table %>%
		dplyr::select(c('Term Description', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value'))
		#table$'geneID' <- gsub('/', ',', x = table$'geneID')

		DT::datatable(
		table,
		filter = 'bottom',
		selection = 'multiple',
		escape = FALSE,
		autoHideNavigation = TRUE,
		rownames = FALSE,
		extensions = c('Buttons'),
		class = 'cell-border stripe',
		options = list(
		dom = 'Bfrtip',
		lengthMenu = c(15, 30, 50, 100),
		pageLength = 10,
		buttons = list(
		'colvis',
		list(
		extend = 'collection',
		text = 'Download/Copy',
		buttons = c('copy', 'csv', 'excel')
		)
		)
		)
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