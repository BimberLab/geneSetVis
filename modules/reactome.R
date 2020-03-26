source('fxs.R', local = TRUE)

reactomeModule <- function(session, input, output, envir) {
	reactomeResults <- reactiveValues(
		results = NULL
	)

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

	observeEvent(input$runreactome_button, {
		#TODO: validate input present?
		validate(need(input$reactome_OrgDB_input != '', "Please select OrgDB..."))

		print('making Reactome query')
		withProgress(message = 'making reactomePA query...', {
			saveFile <- paste0('SavedRuns/', 'running', '_reactomePA_result', '.rds', sep = '')
			if (file.exists(saveFile)) {
				reactomeResults$results <- readRDS(saveFile)
			} else {
				entrezIDs <- bitr(geneID = envir$gene_list$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=input$reactome_OrgDB_input)
				reactomePAres <- ReactomePA::enrichPathway(entrezIDs$ENTREZID, readable = T)
				saveRDS(reactomePAres, file = saveFile)

				reactomeResults$results <- reactomePAres
			}
		})
	})

	renderPlotSet(output = output,
		key = 'reactome',
		enrichTypeResult = reactive(reactomeResults$results),
		termURL = "https://reactome.org/PathwayBrowser/#/",
		datasetName = 'Reactome'
	)
}