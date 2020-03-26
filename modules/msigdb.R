source('fxs.R', local = TRUE)


runMSigDB <- function(DEtable, species) {
	##get species datdset
	human.msig = msigdbr::msigdbr(species = species)

	##subset columms of interest: gene-set name (gs_name) and gene symbols or enterez id
	msigTerm = human.msig %>% dplyr::select(gs_name, gene_symbol, gs_cat, gs_subcat) %>% as.data.frame()

	##dedup table to remove multiple tests
	if (!is.null(DEtable$test)){
		DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
		DEtable <- DEtable[match(unique(DEtable$gene), DEtable$gene), ]
	}

	return_list = list()
	tryCatch({
		clusterTable <- DEtable

		if (nrow(clusterTable) > 0) {
			##Use the gene sets data frame for clusterProfiler (for genes as gene symbols)
			msig_enricher <-
			clusterProfiler::enricher(gene = clusterTable$gene, TERM2GENE = msigTerm)
			#msig_enricher_plot <- dotplot(msig_enricher)

			#clusterProfiler::geneInCategory()
			#geneInCategory(msig_enricher)[as.data.frame(msig_enricher)$ID == 'WINTER_HYPOXIA_METAGENE'][1]

			# enricher_KEGG <- enrichKEGG(
			#   clusterTable$gene,
			#   organism = 'hsa',
			#   keyType = 'kegg',
			#   pAdjustMethod = 'BH'
			# )

			#msig_enricher <- as.data.frame(msig_enricher)
			#msig_enricher$geneID <- gsub(x = msig_enricher$geneID, pattern = '/', replacement = ',')
			addSubset = paste('enricher_result', sep = '')
			return_list[[addSubset]] <- msig_enricher

			#.......................................
			##Use the gene sets data frame for fgsea.
			msig_geneSet = human.msig %>% split(x = .$gene_symbol, f = .$gs_name)

			##name the marker genes with their avgLogFC
			ranks <- clusterTable$avg_logFC
			ranks <- setNames(ranks, clusterTable$gene)

			set.seed(1234)
			fgsea_results <- fgsea(
			pathways = msig_geneSet,
			stats = ranks,
			minSize = 5,
			maxSize = 600,
			nperm = 10000
			)

			threshold <- 0.001
			sigPathways.sum <- sum(fgsea_results[, padj < threshold])

			topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n = 10), pathway]
			topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n = 10), pathway]
			topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

			fgsea_gtable <- fgsea::plotGseaTable(
				pathways = msig_geneSet[topPathways],
				stats = ranks,
				fgseaRes = fgsea_results,
				gseaParam = 0.5,
				render = F,
				colwidths = c(5, 3, 0.8, 1.2, 1.2)
			)

			plot(fgsea_gtable)

			addSubset = paste('fgsea_result', sep = '')
			return_list[[addSubset]] <- fgsea_results

			addSubset = paste('fgsea_gtable', sep = '')
			return_list[[addSubset]] <- fgsea_gtable

			addSubset = paste('fgsea_ranks', sep = '')
			return_list[[addSubset]] <- ranks

			addSubset = 'msig_geneSet'
			return_list[[addSubset]] <- msig_geneSet
		}
	})

	return(return_list)
}

msigdbModule <- function(session, input, output, envir) {
	msigdbResults <- reactiveValues(
		results = NULL
	)

	observeEvent(input$msigdbr_category_input, {
		req(input$msigdbr_category_input)

		species.msig <- msigdbr::msigdbr(species = input$msigdbr_species_input)
		subcat <- unique(filter(species.msig, species.msig$gs_cat == input$msigdbr_category_input)$gs_subcat)

		#See: https://shiny.rstudio.com/reference/shiny/1.2.0/updateSelectInput.html
		updateSelectInput(session, "msigdbr_subcategory_input",
			choices = c('', subcat),
			selected = ''
		)
	})

	observeEvent(input$runmsigdbr_button, {
		print('making MSigDB query')

		req(input$msigdbr_species_input)

		#NOTE: the only purpose of this save is for rapid debugging.
		#Should ultimately review Shiny's existing caching mechanism and use this
		withProgress(message = 'making MSigDB query..', {
			saveFile <- paste0('SavedRuns/', 'running', '_msig_result', '.rds', sep = '')
			if (file.exists(saveFile)) {
				msigdbrRes <- readRDS(saveFile)
			} else {
				msigdbrRes <- runMSigDB(DEtable = envir$gene_list, species = input$msigdbr_species_input)
				saveRDS(msigdbrRes, file = saveFile)
			}
		})
		msigdbResults$result <- msigdbrRes

		msigdbResults$enricher_result <- msigdbrRes[['enricher_result']]
		msigdbResults$fgsea_result <- as.enrichResult(
			result = msigdbrRes[['fgsea_result']],
			inputIds = msigdbrRes[['enricher_result']]@gene,
			geneSet = msigdbrRes[['enricher_result']]@geneSets
		)
	})

	renderPlotSet(output = output,
		key = 'fgsea',
		enrichTypeResult = reactive(msigdbResults$fgsea_result),
		termURL = 'https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=',
		datasetName = 'MSigDB'
	)

	renderPlotSet(output = output,
		key = 'enricher',
		enrichTypeResult = reactive(msigdbResults$enricher_result),
		termURL = 'https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=',
		datasetName = 'MSigDB'
	)

	output$num_of_mapped_enricher <- flexdashboard::renderValueBox({
		num_genes_mapped <- str_split(noquote(msigdbResults$enricher_result@result$GeneRatio[1]), '/')[[1]][2]
		shinydashboard::box(
		title = 'Number of genes mapped',
		width = 6,
		background = 'light-blue',
		num_genes_mapped
		)
	})

	output$num_of_total_genes_enricher <- flexdashboard::renderValueBox({
		num_genes_total <- length(msigdbResults$enricher_result@gene)
		shinydashboard::box(
		title = 'Number of genes total',
		width = 6,
		background = 'light-blue',
		num_genes_total
		)
	})

	output[["fgsea_table_PPI"]] <- renderPlot({
		req(length(input$fgsea_table_cell_clicked) > 0)
		info_list <- input$fgsea_table_cell_clicked

		cellVal <- strsplit(info_list[["value"]], perl = T, split = 'target=\"_blank\">')
		cellVal <- strsplit(cellVal[[1]][2], perl = T, split = "</a>")
		cellVal <- cellVal[[1]]

		fgsea::plotEnrichment(msigdbResults$result[['msig_geneSet']][[cellVal]], msigdbResults$result[['fgsea_ranks']]) +
		ggplot2::labs(title = cellVal)
	})

	observeEvent(input$fgsea_table_cell_clicked, {
		showModal(modalDialog(
			plotOutput("fgsea_table_PPI")
		))
	})

	observeEvent(input$msigdbr_resource_info, {
		showModal(modalDialog(
			msigdbr_resource_info[["text"]],
			title = msigdbr_resource_info[["title"]],
			easyClose = TRUE,
			footer = NULL
		))
	})
}