
as.enrichResult <- function(result, inputIds, geneSet) {
  if (nrow(result) == 0) {stop('No terms in input.')}

  gene <- inputIds
  gene.length <- length(gene)
  
  result <- result %>% 
    dplyr::rename('Count' = size, 'p.adjust' = padj, 'pvalue' = pval, Description = 'pathway') %>% 
    dplyr::arrange(p.adjust)

  result$GeneRatio <- paste(result$Count, '/', gene.length, sep = '')
  result$size <- result$Count
  result$ID <- result$Description
  
  geneSetsOI <- geneSet[c(result$Description)]
  genesInGeneSet <- lapply(geneSetsOI, intersect, y=gene)
  genesInGeneSet.stack <- stack(genesInGeneSet) %>% 
    rename(ind = 'Description') %>% group_by(Description) %>%
    summarise(geneID = paste(values, collapse = '/'))
  
  result <- merge(result, genesInGeneSet.stack, by = 'Description')
  rownames(result) <- result$Description

  new(
    'enrichResult',
    result         = result,
    pvalueCutoff   = 0.05,
    pAdjustMethod  = 'UNKNOWN',
    #qvalueCutoff   = 1,
    gene           = as.character(gene),
    #universe       = extID,
    geneSets       = geneSet,
    organism       = 'UNKNOWN',
    keytype        = 'UNKNOWN',
    ontology       = 'UNKNOWN',
    readable       = T
  )
}


hyperlink_text <- function(url, text, hide=text) {
  for (t in text) {
   s <-  paste0('<a href="',url,hide,'" target="_blank">',text,'</a>')
   return(s)
    }
}


multi_hyperlink_text <- function(labels, links){
  out <- mapply(
    function(text, url){
      dat <- hyperlink_text(text, url = url)
      dat <- split(dat, seq_along(text))
      }, 
    text = strsplit(labels, split = ","),
    url = strsplit(links, split = ","), SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  
  out <- sapply(out, paste, collapse=",")
  return(as.list(out))
}


makeDiskCacheKey <- function(parameters) {
  hashConcat <- ''
  for (param in parameters) {
   hashSub <-  substr(digest::digest(param), 1, 10)
      
   hashConcat <- paste0(hashConcat, hashSub)
  }
  return(hashConcat)
}


makeTermsTable <- function(table, genesDelim,
                           termURL, 
                           caption = NULL, 
                           includeColumns = c('Term Description', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term')) {
  
  table$'Genes in Term' <- gsub(pattern = genesDelim, replacement = ',', x = table$'Genes in Term')
  
  table$'Term Description' <- hyperlink_text(url = termURL, text = table$'Term Description', hide = table$'Term ID')
  table$'Genes in Term' <- multi_hyperlink_text(labels = table$'Genes in Term', links = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=")
  
  table <- table %>%
    dplyr::select(tidyselect::all_of(includeColumns)) %>%
    dplyr::arrange(`p-Value (adj.)`) %>%
    DT::datatable(
      caption = caption,
      filter = 'bottom',
      selection = 'single',
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c('Buttons'),
      class = 'cell-border stripe',
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(15, 30, 50, 100, -1), c('15', '30', '50', '100', 'All')),
        pageLength = 10,
        scrollX = TRUE,
        buttons = list(
          'colvis',
          list(
            extend = 'collection',
            text = 'Download/Copy',
            buttons = c('copy', 'csv', 'excel')
          )
        )
      )
    ) # %>% formatStyle( 0, target= 'row',color = 'black', backgroundColor = NULL, fontWeight = NULL, lineHeight='50%')
}


renderPlotSet <- function(output, key, enrichTypeResult, termURL, datasetName = NULL, caption = NULL) {
  output[[paste(key, 'table', sep = '_')]] <- renderDataTable(server = FALSE, {
    er <- enrichTypeResult()
    validate(need(!is.null(er), paste0('')))
    validate(need(class(er) == 'enrichResult', 'Input should be of enrichResult type...'))
    validate(need(nrow(er) != 0, 'No enriched terms found.'))
    table <- er %>% as.data.frame() %>% 
      dplyr::rename(
      'Term Description' = Description,
      'Term ID' = ID,
      'geneID' = geneID,
      'Hits' = Count,
      'p-Value (adj.)' = pvalue,
      'p-Value' = p.adjust,
      'Genes in Term' = geneID
    )  
    
    makeTermsTable(table = table, genesDelim = '/',
                   termURL = termURL, 
                   caption = caption,
                   includeColumns = c('Term Description', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term'))  })
  
  output[[paste(key, 'dotplot', sep = '_')]] <- renderPlotly({
    er <- enrichTypeResult()
    validate(need(!is.null(er), paste0('Please Run ', datasetName,' on input...')))
    validate(need(nrow(er) != 0, 'No enriched terms found.'))
    enrichplot::dotplot(er)
  })
  
  output[[paste(key, 'emapplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er), paste0('Please Run ', datasetName,' on input...')))
    validate(need(nrow(er) != 0, 'No enriched terms found.'))
    enrichplot::emapplot(er)
  })
  
  output[[paste(key, 'cnetplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er), paste0('Please Run ', datasetName,' on input...')))
    validate(need(nrow(er) != 0, 'No enriched terms found.'))
    enrichplot::cnetplot(er)
  })
  
  output[[paste(key, 'upsetplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er), paste0('Please Run ', datasetName,' on input...')))
    validate(need(nrow(er) != 0, 'No enriched terms found.'))
    enrichplot::upsetplot(er)
  })
  
  output[[paste(key, 'heatplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er), paste0('Please Run ', datasetName,' on input...')))
    validate(need(nrow(er) != 0, 'No enriched terms found.'))
    enrichplot::heatplot(er)
  })
}


makeTabBox <- function(title, key) {
  shinydashboard::tabBox(
    title = title,   
    side = 'right',
    height = NULL,
    selected = 'Table',
    width = 16,
    tabPanel('Table', dataTableOutput(paste(key, 'table', sep = '_'))),
    tabPanel('Dot Plot', plotly::plotlyOutput(paste(key, 'dotplot', sep = '_'))),
    tabPanel('Emap Plot', plotOutput(paste(key, 'emapplot', sep = '_'))),
    tabPanel('Cnet Plot', plotOutput(paste(key, 'cnetplot', sep = '_'))),
    tabPanel('Upset Plot', plotOutput(paste(key, 'upsetplot', sep = '_'))),
    tabPanel('Heat Plot', plotOutput(paste(key, 'heatplot', sep = '_')))
  )
}


