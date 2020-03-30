
runMSigDB <- function(DEtable, species, category = NULL, subcategory = NULL) {
	##get species datdset
  species.msig = msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)

	##subset columms of interest: gene-set name (gs_name) and gene symbols or enterez id
  msigTerm = species.msig %>% dplyr::select(gs_name, gene_symbol, gs_cat, gs_subcat) %>% as.data.frame()

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
			msig_enricher <- clusterProfiler::enricher(gene = clusterTable$gene, TERM2GENE = msigTerm)
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

msigdbModule <- function(session, input, output, envir, appDiskCache) {
	msigdbResults <- reactiveValues(
		results = NULL,
		enricher_result = NULL,
		fgsea_result = NULL
	)

	#NOTE: this should reset our tab whenever the input genes change
	observeEvent(envir$gene_list, {
		print('resetting msigdb')
		msigdbResults$results <- NULL
		msigdbResults$enricher_result <- NULL
		msigdbResults$fgsea_result <- NULL
	})

	msigSubcategories <- read.table(file = './data/msigdb_categories.txt', sep = '\t', header = T, stringsAsFactors = FALSE)
	msigSubcategories <- msigSubcategories[!is.na(msigSubcategories$Subcategory) & msigSubcategories$Subcategory != '',]
	msigSubcategories$SubcategoryLabel <- paste0(msigSubcategories$Subcategory, ': ', msigSubcategories$SubcategoryLabel)

	observeEvent(input$msigdbr_category_input, {
		req(input$msigdbr_category_input)

		#NOTE: these do not appear to be specieis-specific, at least on the website.  This call introduces a lot of lag time
		#species.msig <- msigdbr::msigdbr(species = input$msigdbr_species_input)

		subcat <- msigSubcategories[msigSubcategories$Category == input$msigdbr_category_input,]
		subcat <- subcat[c('Subcategory', 'SubcategoryLabel')]
		subcatValues <- subcat$Subcategory
		names(subcatValues) <- subcat$SubcategoryLabel
		subcatValues <- sort(subcatValues)

		updateSelectInput(session, "msigdbr_subcategory_input",
			choices = c('', subcatValues),
			selected = ''
		)
	})

	observeEvent(input$runmsigdbr_button, {
		print('making MSigDB query')

		req(input$msigdbr_species_input)

		withProgress(message = 'making MSigDB query..', {
		  category <- input$msigdbr_category_input
		  subcategory <- input$msigdbr_subcategory_input
		  
		  if (category == '') {category <- NULL}     
		  if (subcategory == '') {subcategory <- NULL} 
		  
		  cacheKey <- makeDiskCacheKey(list(envir$gene_list, input$msigdbr_species_input, category, subcategory), 'msigdb')
		  cacheVal <- appDiskCache$get(cacheKey)
		  if (class(cacheVal) == 'key_missing') {
		    print('missing cache key...')
		    msigdbrRes <- runMSigDB(
		      DEtable = envir$gene_list,
		      species = input$msigdbr_species_input,
		      category = category,
		      subcategory = subcategory
		    )
		    appDiskCache$set(key = cacheKey, value = msigdbrRes)
			} else {
			  print('loading from cache...')
			  msigdbrRes <- cacheVal
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
	
	renderPlotSet(
	  output = output,
	  key = 'enricher',
	  enrichTypeResult = reactive(msigdbResults$enricher_result),
	  termURL = 'https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=',
	  datasetName = 'MSigDB'
	)
	
	renderPlotSet(
	  output = output,
	  key = 'fgsea',
	  enrichTypeResult = reactive(msigdbResults$fgsea_result),
	  termURL = 'https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=',
	  datasetName = 'MSigDB',
	  caption = 'Click on Term Description cell to view enrichment plot.'
	)
	
	
	output$msig_map_stats <- renderText({
	  validate(need(!is.null(msigdbResults$enricher_result), "No mapped genes."))
	  num_genes_mapped <- str_split(noquote(msigdbResults$enricher_result@result$GeneRatio[1]), '/')[[1]][2]
	  HTML(
	    '<b>Mapped genes</b><br>',
	    paste0(num_genes_mapped, ' out of ', length(envir$gene_list$gene), ' genes were mapped.')
	  )
	})
	


	output[["fgsea_table_PPI"]] <- renderPlot({
		info_list <- input$fgsea_table_cell_clicked

		cellVal <- strsplit(info_list[["value"]], perl = T, split = 'target=\"_blank\">')
		cellVal <- strsplit(cellVal[[1]][2], perl = T, split = "</a>")
		cellVal <- cellVal[[1]]

		fgsea::plotEnrichment(msigdbResults$result[['msig_geneSet']][[cellVal]], msigdbResults$result[['fgsea_ranks']]) +
		ggplot2::labs(title = cellVal)
	})

	observeEvent(input$fgsea_table_cell_clicked, {
	  req(length(input$fgsea_table_cell_clicked) > 0)
		showModal(modalDialog(
			plotOutput("fgsea_table_PPI")
		))
	})

	observeEvent(input$msigdbr_resource_info, {
		msigdbr_resource_info <- list(
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
				'
			)
		)

		showModal(modalDialog(
			msigdbr_resource_info[["text"]],
			title = msigdbr_resource_info[["title"]],
			easyClose = TRUE,
			footer = NULL
		))
	})
}