runNCG <- function(DEtable, geneCol, species) {
  ##dedup table to remove multiple tests
  if (!is.null(DEtable$test)){
    DEtable <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
    DEtable <- DEtable[match(unique(DEtable[[geneCol]]), DEtable[[geneCol]]), ]
  }

  return_list = list()
  tryCatch({
    entrezIDs <- mapIds(org.Hs.eg.db, as.character(DEtable[[geneCol]]), 'ENTREZID', 'SYMBOL')

    ncg <- DOSE::enrichNCG(entrezIDs, OrgDb = human, pvalueCutoff = 1, qvalueCutoff = 1)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ncgModule <- function(session, input, output, envir, appDiskCache) {

  #NOTE: this should reset our tab whenever the input genes change
  observeEvent(list(envir$geneList), ignoreInit = F, {
    print('resetting ncg')
    envir$ncgRes <- NULL
    errEl <- NULL
    if (!is.null(errEl)) {shinyjs::hide(errEl)}
  })

  observeEvent(input$runncg_button, {
    withBusyIndicatorServer("runncg_button", {
      Sys.sleep(1)
      #TODO: validate input present?
      validate(need(input$ncg_OrgDB_input != '', "Please select OrgDB..."))

      print('making ncg query')
      withProgress(message = 'making NCG query...', {
        cacheKey <- makeDiskCacheKey(list(envir$geneList[[input$ncg_selectGeneCol]], input$ncg_OrgDB_input), 'ncg')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')

          #if (!require(input$ncg_OrgDB_input)) install.packages(input$ncg_OrgDB_input)
          envir$ncgRes <- NULL
          fromType <- ifelse(grepl('id', input$ncg_selectGeneCol), 'ENSEMBL', 'SYMBOL')
          entrezIDs <- clusterProfiler::bitr(geneID = envir$geneList[[input$ncg_selectGeneCol]], fromType=fromType, toType="ENTREZID", OrgDb=input$ncg_OrgDB_input)
          ncgRes <- DOSE::enrichNCG(entrezIDs$ENTREZID, readable = T)
          appDiskCache$set(key = cacheKey, value = ncgRes)
        } else {
          print('loading from cache...')
          ncgRes <- cacheVal
        }

        envir$ncgRes <- ncgRes
        if ( is.null(envir$ncgRes) || nrow(envir$ncgRes) == 0 ) {stop('No significant enrichment found.')}

      })
    })
  })

  renderPlotSet(
    output = output,
    key = 'ncg',
    enrichTypeResult = reactive(envir$ncgRes),
    datasetURL = "",
    datasetName = 'ncg',
    namedGeneList = envir$namedGeneList
  )

  output$ncg_map_stats <- renderText({
    validate(need(!is.null(envir$ncgRes) & length(envir$ncgRes) != 0, "No mapped genes."))
    if (nrow(envir$ncgRes@result) > 0) {
      num_genes_mapped <- stringr::str_split(noquote(envir$ncgRes@result$GeneRatio[1]), '/')[[1]][2]
    } else {
      num_genes_mapped <- 0
    }
    HTML(
      '<b>Mapped genes</b><br>',
      paste0(num_genes_mapped, ' out of ', length(envir$geneList[[input$ncg_selectGeneCol]]), ' genes were mapped.')
    )
  })

  observeEvent(input$ncg_resource_info, {
    ncg_resource_info <- list(
      title = "NCG Resource info",
      text = HTML(
        '<b>Network of Cancer Genes</b><br>
				NCG contains information on duplicability, evolution, protein-protein and microRNA-gene interaction, function, expression and essentiality of 2,372 cancer genes from 273 manually curated publications
				<p>
				<li><a href=http://ncg.kcl.ac.uk/index.php
				title="Network of Cancer Genes Browser website"
				target="_blank"><b>Network of Cancer Genes Browser website</b></a></li>
				<p>
				<p>
				<b>enrichplot</b><br>
				The plots are produced using <a href=https://bioconductor.org/packages/release/bioc/html/enrichplot.html target="_blank"><b>enrichplot</b></a> by
				<a href=https://github.com/YuLab-SMU/enrichplot target="_blank"><b>Yu G (2019) </b></a>
				'
      )
    )

    showModal(modalDialog(
      ncg_resource_info[["text"]],
      title = ncg_resource_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    ))
  })


}
