runDGN <- function(DEtable, geneCol, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable[[geneCol]]), DEtable[[geneCol]]), ]
  }

  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable[[geneCol]]), 'ENTREZID', 'SYMBOL')

    dgn <- DOSE::enrichDGN(entrezIDs, OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

dgnModule <- function(session, input, output, envir, appDiskCache) {

  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(list(envir$geneList), ignoreInit = F, {
    envir$dgnRes <- NULL
    errEl <- NULL
    if (!is.null(errEl)) {shinyjs::hide(errEl)}
  })

  observeEvent(input$rundgn_button, {
    withBusyIndicatorServer("rundgn_button", {
      validate(need(!is.null(envir$geneList) && nrow(envir$geneList) > 0, "Please enter genes into Load Data tab"))
      validate(need(input$dgn_selectGeneCol != '', "Please select gene column to use"))
      validate(need(input$dgn_OrgDB_input != '', "Please select OrgDB..."))

      shinybusy::show_modal_spinner(text = 'Querying DisGeNET. This might take some time.')

      cacheKey <- makeDiskCacheKey(list(envir$geneList[[input$dgn_selectGeneCol]], input$dgn_OrgDB_input), 'dgn')
      cacheVal <- appDiskCache$get(cacheKey)
      if (class(cacheVal) == 'key_missing') {
        print('missing cache key...')

        #if (!require(input$dgn_OrgDB_input)) install.packages(input$dgn_OrgDB_input)
        envir$dgnRes <- NULL
        fromType <- ifelse(grepl('id', input$dgn_selectGeneCol), 'ENSEMBL', 'SYMBOL')
        entrezIDs <- clusterProfiler::bitr(geneID = envir$geneList[[input$dgn_selectGeneCol]], fromType=fromType, toType="ENTREZID", OrgDb=input$dgn_OrgDB_input)
        dgnRes <- DOSE::enrichDGN(entrezIDs$ENTREZID, readable = T)
        appDiskCache$set(key = cacheKey, value = dgnRes)
      } else {
        print('loading from cache...')
        dgnRes <- cacheVal
      }

      remove_modal_progress()

      if ( is.null(dgnRes) || nrow(dgnRes) == 0 ) {
        stop('No significant enrichment found.')
      }

      dgnRes@result$ID <- gsub(pattern = 'umls:', replacement = '', dgnRes@result$ID)
      rownames(dgnRes@result) <- dgnRes@result$ID

      envir$dgnRes <- dgnRes
    })
  })

  renderPlotSet(
    output = output,
    key = 'dgn',
    enrichTypeResult = reactive(envir$dgnRes),
    datasetURL = "https://www.disgenet.org",
    datasetName = 'DisGeNET',
    namedGeneList = envir$namedGeneList
  )

  output$dgn_map_stats <- renderText({
    validate(need(!is.null(envir$dgnRes) & length(envir$dgnRes) != 0, "No mapped genes."))
    num_genes_mapped <- stringr::str_split(noquote(envir$dgnRes@result$GeneRatio[1]), '/')[[1]][2]
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$geneList[[input$dgn_selectGeneCol]]), ' genes were mapped.')
    )
  })

  observeEvent(input$dgn_resource_info, {
    dgn_resource_info <- list(
      title = "DGN Resource info",
      text = HTML(
        '<b>DisGeNET Browser</b><br>
				DisGeNET integrates data from expert curated repositories, GWAS catalogues, animal models and the scientific literature. DisGeNET data are homogeneously annotated with controlled vocabularies and community-driven ontologies. Additionally, several original metrics are provided to assist the prioritization of genotype–phenotype relationships.
				<p>
				<li><a href=https://www.disgenet.org/home/
				title="DisGeNET Browser website"
				target="_blank"><b>DisGeNET Browser website</b></a></li>
				<p>
				<p>
				<b>enrichplot</b><br>
				The plots are produced using <a href=https://bioconductor.org/packages/release/bioc/html/enrichplot.html target="_blank"><b>enrichplot</b></a> by
				<a href=https://github.com/YuLab-SMU/enrichplot target="_blank"><b>Yu G (2019) </b></a>
				'
      )
    )

    showModal(modalDialog(
      dgn_resource_info[["text"]],
      title = dgn_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })


}
