
runMSigDB <- function(DEtable, geneCol, species, category = NULL, subcategory = NULL) {
	##get species datdset
  species.msig = msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)

	##subset columms of interest: gene-set name (gs_name) and gene symbols or enterez id
  msigTerm = species.msig %>% dplyr::select(gs_name, gene_symbol, entrez_gene, gs_cat, gs_subcat) %>% as.data.frame()

	##dedup table to remove multiple tests
	if (!is.null(DEtable$test)){
		DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
		DEtable <- DEtable[match(unique(DEtable[[geneCol]]), DEtable[[geneCol]]), ]
	}

	return_list = list()
	tryCatch({
		clusterTable <- DEtable

		if (nrow(clusterTable) > 0) {

			##Use the gene sets data frame for clusterProfiler (for genes as gene symbols)
			msig_enricher <- clusterProfiler::enricher(gene = clusterTable[[geneCol]], TERM2GENE = msigTerm)
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
			msig_geneSet = species.msig %>% split(x = .$gene_symbol, f = .$gs_name)

			##name the marker genes with their avgLogFC
			ranks <- clusterTable$avg_logFC
			ranks <- setNames(ranks, clusterTable[[geneCol]])

			set.seed(1234)
			fgsea_results <- fgsea::fgsea(
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

msigdbModule <- function(session, input, output, envir, appDiskCache) {

	#NOTE: this should reset our tab whenever the input genes change
	observeEvent(list(envir$geneList), ignoreInit = F, {
		print('resetting msigdb')
		envir$msigdbRes <- NULL
		envir$msigdbRes_fgsea <- NULL
		errEl <- NULL
		if (!is.null(errEl)) {shinyjs::hide(errEl)}
	})

	print(getwd())
  msigdbSubcategories <- read.table(file = '../msigdb_categories.txt', sep = '\t', header = T, stringsAsFactors = FALSE)
	#msigdbSubcategories <- read.table(file = system.file('msigdb_categories.txt', package = 'geneSetVis'), sep = '\t', header = T, stringsAsFactors = FALSE)
	msigdbSubcategories <- msigdbSubcategories[!is.na(msigdbSubcategories$Subcategory) & msigdbSubcategories$Subcategory != '',]
	msigdbSubcategories$SubcategoryLabel <- paste0(msigdbSubcategories$Subcategory, ': ', msigdbSubcategories$SubcategoryLabel)

	observeEvent(input$msigdb_category_input, {
		req(input$msigdb_category_input)

		#NOTE: these do not appear to be specieis-specific, at least on the website.  This call introduces a lot of lag time
		#species.msig <- msigdbr::msigdbr(species = input$msigdbr_species_input)

		subcat <- msigdbSubcategories[msigdbSubcategories$Category == input$msigdb_category_input,]
		subcat <- subcat[c('Subcategory', 'SubcategoryLabel')]
		subcatValues <- subcat$Subcategory
		names(subcatValues) <- subcat$SubcategoryLabel
		subcatValues <- sort(subcatValues)

		updateSelectInput(session, "msigdb_subcategory_input",
			choices = c('', subcatValues),
			selected = ''
		)
	})

	observeEvent(input$runmsigdb_button, {
	  withBusyIndicatorServer("runmsigdb_button", {
	    Sys.sleep(1)

	    print('making MSigDB query')

	    req(input$msigdb_species_input)

	    withProgress(message = 'making MSigDB query..', {
	      category <- input$msigdb_category_input
	      subcategory <- input$msigdb_subcategory_input

	      if (category == '') {category <- NULL}
	      if (subcategory == '') {subcategory <- NULL}

	      cacheKey <- makeDiskCacheKey(list(envir$geneList, input$msigdb_selectGeneCol, input$msigdb_species_input, category, subcategory), 'msigdb')
	      cacheVal <- appDiskCache$get(cacheKey)
	      if (class(cacheVal) == 'key_missing') {
	        print('missing cache key...')

	        type <- ifelse(grepl('id', input$msigdb_selectGeneCol), 'ENSEMBL', 'SYMBOL')
	        if (type != 'SYMBOL') {
	          envir$msigdbRes$enricher_result <- NULL
	          envir$msigdbRes$fgsea_result <- NULL
	          stop('MsigDB only accepts gene columns of type SYMBOL')
	          }

	        msigdbRes <- runMSigDB(
	          DEtable = envir$geneList,
	          geneCol = input$msigdb_selectGeneCol,
	          species = input$msigdb_species_input,
	          category = category,
	          subcategory = subcategory
	        )
	        appDiskCache$set(key = cacheKey, value = msigdbRes)
	      } else {
	        print('loading from cache...')
	        msigdbRes <- cacheVal
	      }
	    })

	    envir$msigdbRes <- msigdbRes

	    if (is.null(envir$msigdbRes$enricher_result) | nrow(envir$msigdbRes$enricher_result) == 0) {
	      envir$msigdbRes$enricher_result <- NULL
	      envir$msigdbRes$fgsea_result <- NULL
	      stop('No significant enrichment found.')
	      }

	    if ( !is.null(envir$msigdbRes$fgsea_result) & nrow(envir$msigdbRes$fgsea_result) != 0 ) {

	      fgsea_geneIDCol <-
	        getEnrichResGeneID(
	          gseResult = msigdbRes$fgsea_result,
	          idCol = msigdbRes$fgsea_result$pathway,
	          gseGenes = msigdbRes$enricher_result@gene,
	          geneSet = msigdbRes$msig_geneSet,
	          idColName = 'pathway'
	        )

	      envir$msigdbRes_fgsea <- as.enrichResult(
	        gseType = 'FGSEA',
	        gseResult = msigdbRes$fgsea_result,
	        gseGenes = msigdbRes$enricher_result@gene,
	        idCol = msigdbRes$fgsea_result$pathway,
	        padjCol = msigdbRes$fgsea_result$padj,
	        pvalCol = msigdbRes$fgsea_result$pval,
	        geneIDCol = fgsea_geneIDCol,
	        pvalueCutoff = 0.5,
	        #?size is not Count
	        #countCol = msigdbRes$fgsea_result$size,
	        countCol = lapply(stringr::str_split(fgsea_geneIDCol, pattern = '/'), length),
	        geneRatioCol = paste(
	          lapply(stringr::str_split(fgsea_geneIDCol, pattern = '/'), length),
	          '/',
	          length(msigdbRes$enricher_result@gene),
	          sep = ''
	        )
	      )
	    }
	  })
	})

	renderPlotSet(
	  output = output,
	  key = 'enricher',
	  enrichTypeResult = reactive(envir$msigdbRes$enricher_result),
	  datasetURL = 'https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=',
	  datasetName = 'MSigDB',
	  namedGeneList = envir$namedGeneList
	)

	renderPlotSet(
	  output = output,
	  key = 'fgsea',
	  enrichTypeResult = reactive(envir$msigdbRes_fgsea),
	  datasetURL = 'https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=',
	  datasetName = 'MSigDB',
	  caption = 'Click on Term Description cell to view enrichment plot.',
	  namedGeneList = envir$namedGeneList
	)


	output$msigdb_map_stats <- renderText({
	  validate(need(!is.null(envir$msigdbRes$enricher_result), "No mapped genes."))
	  num_genes_mapped <- stringr::str_split(noquote(envir$msigdbRes$enricher_result@result$GeneRatio[1]), '/')[[1]][2]
	  HTML(
	    '<b>Mapped genes</b><br>',
	    paste0(num_genes_mapped, ' out of ', length(envir$geneList$gene), ' genes were mapped.')
	  )
	})



	output[["fgsea_table_PPI"]] <- renderPlot({
		info_list <- input$fgsea_table_cell_clicked

		cellVal <- strsplit(info_list[["value"]], perl = T, split = 'target=\"_blank\">')
		cellVal <- strsplit(cellVal[[1]][2], perl = T, split = "</a>")
		cellVal <- cellVal[[1]]

		fgsea::plotEnrichment(envir$msigdbRes$msig_geneSet[[cellVal]], envir$msigdbRes$fgsea_ranks) +
		ggplot2::labs(title = cellVal)
	})

	observeEvent(input$fgsea_table_cell_clicked, {
	  req(length(input$fgsea_table_cell_clicked) > 0)
		showModal(modalDialog(
			plotOutput("fgsea_table_PPI")
		))
	})

	observeEvent(input$msigdb_resource_info, {
		msigdb_resource_info <- list(
			title = "MSigDB Resource info",
			text = HTML(
			'<b>MSigDB</b><br>
				The Molecular Signatures Database (MSigDB) is a collection of gene sets
				originally created for use with the Gene Set Enrichment Analysis (GSEA) software.
				Gene homologs are provided by HUGO Gene Nomenclature Committee at the European Bioinformatics Institute
				which integrates the orthology assertions predicted for human genes by
				eggNOG, Ensembl Compara, HGNC, HomoloGene, Inparanoid, NCBI Gene Orthology, OMA, OrthoDB, OrthoMCL, Panther, PhylomeDB, TreeFam and ZFIN.
				For each human equivalent within each species, only the ortholog supported by the largest number of databases is used.
				<p>
				<li><a href="https://www.gsea-msigdb.org/gsea/msigdb/index.jsp"
				title="link to Official MSigDB"
				target="_blank"><b>Official MSigDB website</b></a></li>
				<p>
				<p>
				<b>enrichplot</b><br>
				The plots are produced using <a href=https://bioconductor.org/packages/release/bioc/html/enrichplot.html target="_blank"><b>enrichplot</b></a> by
				<a href=https://github.com/YuLab-SMU/enrichplot target="_blank"><b>Yu G (2019) </b></a>
				'
			)
		)

		showModal(modalDialog(
			msigdb_resource_info[["text"]],
			title = msigdb_resource_info[["title"]],
			easyClose = TRUE,
			footer = NULL
		))
	})
}
