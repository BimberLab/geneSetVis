
getEnrichResGeneID <- function(gseResult, idCol, idColName, gseGenes, geneSet) {
  
  geneSetsOI <- geneSet[idCol]
  
  genesInGeneSet <- lapply(geneSetsOI, intersect, y=gseGenes)
  genesInGeneSet.stack <- stack(genesInGeneSet) %>% 
    rename(ind = idColName) %>% group_by(.dots = idColName) %>%
    summarise(geneID = paste(values, collapse = '/'))
  
  gseResult <- merge(gseResult, genesInGeneSet.stack, by = idColName)
  
  return(gseResult$geneID)
}

as.enrichResult <- function( gseType = 'GSE', gseResult, gseGenes, idCol, descCol = idCol, geneIDCol, countCol, pvalCol, padjCol, geneRatioCol, 
                            bgRatioCol = NULL,  qvalCol = NULL, pvalueCutoff = 0.05, pAdjustMethod = '', qvalueCutoff = 0, 
                            universe = '', geneSets = list(), organism = '', keytype = '', ontology = '', readable = T) {
  
  #if (nrow(gseResult) == 0) {stop(paste0('No terms in ', gseType, ' result.'))}
  if (nrow(gseResult) != 0) {
    result <- NULL
    result$ID <- idCol
    result$Description <- descCol
    result$GeneRatio <- geneRatioCol
    result$BgRatio <- bgRatioCol
    result$pvalue <- pvalCol
    result$p.adjust <- padjCol
    #result$qvalue <- qvalCol
    result$geneID <- geneIDCol
    result$Count <- as.integer(countCol)
    
    result <- data.frame(result)
    
    result <- result %>% dplyr::arrange(p.adjust)
    
    rownames(result) <- result$ID
  } else {
    result <- data.frame(NULL)
  }
  
  
  
  new(
    Class = 'enrichResult',
    result         = result,
    pvalueCutoff   = pvalueCutoff,
    pAdjustMethod  = pAdjustMethod,
    #qvalueCutoff   = qvalueCutoff,
    gene           = as.character(gseGenes),
    universe       = universe,
    geneSets       = geneSets,
    organism       = organism,
    keytype        = keytype,
    ontology       = ontology,
    readable       = readable
  )
}


hyperlink_text <- function(href_base, href_cont, link_text=href_cont) {
  h <-  paste0('<a href="',href_base,href_cont,'" target="_blank">',link_text,'</a>')
  return(h)
}


multi_hyperlink_text <- function(labels, links){
  out <- mapply(
    function(hrefs, link_texts){
      ret <- hyperlink_text(href_base = hrefs, href_cont = link_texts)
      ret <- split(ret, seq_along(link_texts))
      }, 
    link_texts = strsplit(labels, split = ","),
    hrefs = strsplit(links, split = ","), SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  
  out <- sapply(out, paste, collapse=",")
  return(as.list(out))
}


makeDiskCacheKey <- function(inputList, prefix) {
  return(paste0(str_to_lower(prefix), digest::digest(inputList)))
}


makeTermsTable <- function(table, genesDelim,
                           datasetURL, 
                           caption = NULL, 
                           includeColumns = c('Term Description', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term')) {
  
  table$'Genes in Term' <- gsub(pattern = genesDelim, replacement = ',', x = table$'Genes in Term')
  table$'Genes in Term' <- multi_hyperlink_text(labels = table$'Genes in Term', links = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=")
  
  if(!is.null(datasetURL)) {
    table$'Term Description' <- hyperlink_text(href_base = datasetURL, href_cont = table$'Term ID', link_text = table$'Term Description')
  }
  
  
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


renderPlotSet <- function(output, key, enrichTypeResult, datasetURL, datasetName = NULL, caption = NULL) {
  output[[paste(key, 'table', sep = '_')]] <- renderDataTable(server = FALSE, {
    er <- enrichTypeResult()
    validate(need(!is.null(er) & nrow(er) != 0, ''))
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
                   datasetURL = datasetURL, 
                   caption = caption,
                   includeColumns = c('Term Description', 'Hits', 'p-Value (adj.)', 'p-Value', 'Genes in Term'))  })
  
  output[[paste(key, 'dotplot', sep = '_')]] <- renderPlotly({
    er <- enrichTypeResult()
    validate(need(!is.null(er) & nrow(er) != 0, 'No enriched terms.'))
    enrichplot::dotplot(er) 
  })
  
  output[[paste(key, 'emapplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er) & nrow(er) != 0, 'No enriched terms.'))
    enrichplot::emapplot(er)
  })
  
  output[[paste(key, 'cnetplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er) & nrow(er) != 0, 'No enriched terms.'))
    enrichplot::cnetplot(er)
  })
  
  output[[paste(key, 'upsetplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er) & nrow(er) != 0, 'No enriched terms.'))
    enrichplot::upsetplot(er)
  })
  
  output[[paste(key, 'heatplot', sep = '_')]] <- renderPlot({
    er <- enrichTypeResult()
    validate(need(!is.null(er) & nrow(er) != 0, 'No enriched terms.'))
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

#https://github.com/daattali/advanced-shiny/tree/master/busy-indicator
withBusyIndicatorCSS <- "
.btn-loading-container {
margin-left: 10px;
font-size: 1.2em;
}
.btn-done-indicator {
color: green;
}
.btn-err {
margin-top: 10px;
color: red;
}
"

withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    shinyjs::useShinyjs(),
    singleton(tags$head(
      tags$style(withBusyIndicatorCSS)
    )),
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      shinyjs::hidden(
        icon("spinner", class = "btn-loading-indicator fa-spin"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    shinyjs::hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

withBusyIndicatorServer <- function(buttonId, expr) {
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })
  
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                                       time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}
