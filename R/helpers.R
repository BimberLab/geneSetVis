
getEnrichResGeneID <- function(gseResult, idCol, idColName, gseGenes, geneSet) {

  geneSetsOI <- geneSet[idCol]

  genesInGeneSet <- lapply(geneSetsOI, intersect, y=gseGenes)
  genesInGeneSet.stack <- stack(genesInGeneSet) %>%
    dplyr::rename(ind = idColName) %>% dplyr::group_by(.dots = idColName) %>%
    dplyr::summarise(geneID = paste(values, collapse = '/'))

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

