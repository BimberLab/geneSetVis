runDAVID <- function(DEtable, geneCol, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable[[geneCol]]), DEtable[[geneCol]]), ]
  }

  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable[[geneCol]]), 'ENTREZID', 'SYMBOL')

    david <- clusterProfiler::enrichDAVID(entrezIDs, david.user = 'oosap@ohsu.edu', OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

davidModule <- function(session, input, output, envir, appDiskCache) {

  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(list(envir$geneList), ignoreInit = F, {
    envir$davidRes <- NULL
    errEl <- NULL
    if (!is.null(errEl)) {shinyjs::hide(errEl)}
  })

  observeEvent(input$rundavid_button, {
    withBusyIndicatorServer("rundavid_button", {
      validate(need(!is.null(envir$geneList) && nrow(envir$geneList) > 0, "Please enter genes into Load Data tab"))
      validate(need(input$david_selectGeneCol != '', "Please select gene column to use"))
      validate(need(input$david_OrgDB_input != '', "Please select OrgDB"))
      validate(need(input$davidUserEmail != '', "Please enter email"))

      shinybusy::show_modal_spinner(text = 'Querying DAVID. This might take some time.')

      cacheKey <- makeDiskCacheKey(list(envir$geneList[[input$david_selectGeneCol]], input$david_OrgDB_input), 'david')
      cacheVal <- appDiskCache$get(cacheKey)
      if (class(cacheVal) == 'key_missing') {
        print('missing cache key...')

        #if (!require(input$david_OrgDB_input)) install.packages(input$david_OrgDB_input)
        envir$davidRes <- NULL
        fromType <- ifelse(grepl('id', input$david_selectGeneCol), 'ENSEMBL', 'SYMBOL')
        entrezIDs <- clusterProfiler::bitr(geneID = envir$geneList[[input$david_selectGeneCol]], fromType=fromType, toType="ENTREZID", OrgDb=input$david_OrgDB_input)
        davidRes <- clusterProfiler::enrichDAVID(entrezIDs$ENTREZID, david.user = input$davidUserEmail)

        appDiskCache$set(key = cacheKey, value = davidRes)
      } else {
        print('loading from cache...')
        davidRes <- cacheVal
      }

      ##change result replace entrezIDs with SYMBOLS; replace with & to replace
      davidRes <- clusterProfiler::setReadable(davidRes, input$david_OrgDB_input, 'ENTREZID')

      envir$davidRes <- davidRes
      remove_modal_progress()

      if ( is.null(envir$davidRes)|| nrow(envir$davidRes) == 0 ) {
        stop('No significant enrichment found.')
      }
    })
  })

  renderPlotSet(
    output = output,
    key = 'david',
    enrichTypeResult = reactive(envir$davidRes),
    datasetURL = NULL,
    datasetName = 'david',
    namedGeneList = envir$namedGeneList
  )

  output$david_map_stats <- renderText({
    validate(need(!is.null(envir$davidRes) & length(envir$davidRes) != 0, "No mapped genes."))
    if (nrow(envir$davidRes@result) > 0) {
      num_genes_mapped <- stringr::str_split(noquote(envir$davidRes@result$GeneRatio[1]), '/')[[1]][2]
    } else {
      num_genes_mapped <- 0
    }
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$geneList[[input$david_selectGeneCol]]), ' genes were mapped.')
    )
  })

  observeEvent(input$david_resource_info, {
    david_resource_info <- list(
      title = "DAVID Resource info",
      text = HTML(
        '<b>DAVID Browser</b><br>
				The DAVID has been developed as a standardized ontology for human disease with the purpose of providing the biomedical community with consistent, reusable and sustainable descriptions of human disease terms, phenotype characteristics and related medical vocabulary disease concepts through collaborative efforts of researchers at Northwestern University, Center for Genetic Medicine and the University of Maryland School of Medicine, Institute for Genome Sciences.
				<p>
				<li><a href=https://david.ncifcrf.gov/home.jsp
				title="DAVID Browser website"
				target="_blank"><b>DAVID Browser website</b></a></li>
				<p>
				<p>
				<b>enrichplot</b><br>
				The plots are produced using <a href=https://bioconductor.org/packages/release/bioc/html/enrichplot.html target="_blank"><b>enrichplot</b></a> by
				<a href=https://github.com/YuLab-SMU/enrichplot target="_blank"><b>Yu G (2019) </b></a>
				'
      )
    )

    showModal(modalDialog(
      david_resource_info[["text"]],
      title = david_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })


}
