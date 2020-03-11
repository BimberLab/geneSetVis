#' @title plot STRINGdb networks
#'
#' @description Takes Differential Expression table and plots STRINGdb networks 
#' @param DEtable A DE table
#' @param numHits The num of mapped hits to plot
#' @param refSpeciesNum The dataset (see STRINGdb docs) to use as a reference; 9606=Human, ?=Rhesus Macaque
#' @return The PNGs of network plots
#' @keywords STRINGdb
#' @import STRINGdb
#' @export
#' @importFrom

##TODO: save pngs
##change db version


runSTRINGdb <- function(DEtable, maxHitsToPlot, refSpeciesNum, scoreThreshold) {
  string_db <-
    STRINGdb$new(
      version = '10',
      species = refSpeciesNum,
      score_threshold = scoreThreshold,
      input_directory = ''
    )
  
  ##dedup table to remove multiple tests
  DEtable.dedup <-
    DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
  DEtable.dedup <-
    DEtable.dedup[match(unique(DEtable.dedup$gene), DEtable.dedup$gene), ]
  
  DEtable.split <- split(DEtable.dedup, DEtable.dedup$cluster)
  
  return_list = list()
  for (i in as.vector(names(DEtable.split))) {
    tryCatch({
      clusterTable <- DEtable.split[[i]]
      
      if (nrow(DEtable.split[[i]]) > 0) {
        cluster.map <-
          string_db$map(clusterTable, 'gene', removeUnmappedRows = FALSE)
        hits <- cluster.map$STRING_id
        
        max_hits_to_plot <- cluster.map$STRING_id[1:maxHitsToPlot]
        
        enrichmentGO <-
          string_db$get_enrichment(hits,
                                   category = 'Process',
                                   methodMT = 'fdr',
                                   iea = TRUE)
        
        enrichmentKEGG <-
          string_db$get_enrichment(hits,
                                   category = 'KEGG',
                                   methodMT = 'fdr',
                                   iea = TRUE)
        
        
        hit_term_proteins <- string_db$get_term_proteins(enrichmentGO$term_id, hits)
        hit_term_genes <- hit_term_proteins %>% 
          dplyr::select(term_id, preferred_name) %>% 
          dplyr::group_by(term_id) %>% 
          dplyr::summarize('hit_term_genes' = paste0(preferred_name, collapse = ', '))
        
        enrichmentGO <- merge(hit_term_genes, enrichmentGO)
        
        
        hit_term_proteins <- string_db$get_term_proteins(enrichmentKEGG$term_id, hits)
        hit_term_genes <- hit_term_proteins %>% 
          dplyr::select(term_id, preferred_name) %>% 
          dplyr::group_by(term_id) %>% 
          dplyr::summarize('hit_term_genes' = paste0(preferred_name, collapse = ', '))
        
        enrichmentKEGG <- merge(hit_term_genes, enrichmentKEGG)
        
        
        string_db$get_png(max_hits_to_plot, file = paste(i, 'network.png', sep = '_'))
        
        
        network <- string_db$plot_network(max_hits_to_plot)
        
        #______
        ##payload mechanism for upregulated vs downregulated genes:
        ##adds a color column for up vs downregulated genes
        cluster.color <-
          string_db$add_diff_exp_color(cluster.map, logFcColStr = 'avg_logFC')
        # post payload information to the STRING server
        payload_id <-
          string_db$post_payload(cluster.color$STRING_id, colors = cluster.color$color)
        string_db$plot_network(hits, payload_id = payload_id)
        
        ##clustering/community algorithms: ”fastgreedy”, ”walktrap”, ”spinglass”, ”edge.betweenness”.
        networkClustersList <-
          string_db$get_clusters(max_hits_to_plot, algorithm = 'fastgreedy')
        par(mfrow = c(2, 2))
        for (j in seq(1:length(networkClustersList))) {
          string_db$plot_network(networkClustersList[[j]], payload_id = payload_id)
        }
        
        link <- string_db$get_link()
        
        addSubset = paste(i, 'hits', sep = '_')
        return_list[[addSubset]] <- hits
        
        addSubset = paste(i, 'network', sep = '_')
        return_list[[addSubset]] <- network
        
        addSubset = paste(i, 'GO', sep = '_')
        return_list[[addSubset]] <- enrichmentGO
        
        addSubset = paste(i, 'KEGG', sep = '_')
        return_list[[addSubset]] <- enrichmentKEGG
        
        addSubset = paste(i, 'link', sep = '_')
        return_list[[addSubset]] <- link
        
      }
      
    }, error = function(e) {
      cat('\nERROR :', conditionMessage(e), '\n')
    })
    
  }
  return(return_list)
}




runMSigDB <- function(DEtable, species) {
  ##get species datdset
  human.msig = msigdbr(species = species)
  
  ##subset columms of interest: gene-set name (gs_name) and gene symbols or enterez id
  msigTerm = human.msig %>% dplyr::select(gs_name, gene_symbol, gs_cat, gs_subcat) %>% as.data.frame()
  
  ##dedup table to remove multiple tests
  DEtable.dedup <- DEtable[with(DEtable, order(p_val_adj, decreasing = F)), ]
  DEtable.dedup <- DEtable.dedup[match(unique(DEtable.dedup$gene), DEtable.dedup$gene),]
  
  DEtable.split <- split(DEtable.dedup, DEtable.dedup$cluster)
  
  return_list = list()
  for (i in as.vector(names(DEtable.split))) {
    tryCatch({
      clusterTable <- DEtable.split[[i]]
      
      if (nrow(clusterTable) > 0) {
        ##Use the gene sets data frame for clusterProfiler (for genes as gene symbols)
        msig_enricher <- enricher(gene = clusterTable$gene, TERM2GENE = msigTerm)
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
        addSubset = paste(i, 'enricher_result', sep = '_')
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
        print(paste0(
          sigPathways.sum,
          ' significant pathways. pval < ',
          threshold
        ))
        
        topPathwaysUp <-
          fgsea_results[ES > 0][head(order(pval), n = 10), pathway]
        topPathwaysDown <-
          fgsea_results[ES < 0][head(order(pval), n = 10), pathway]
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
        
        
        fgsea_gtable <-
          plotGseaTable(
            pathways = msig_geneSet[topPathways],
            stats = ranks,
            fgseaRes = fgsea_results,
            gseaParam = 0.5,
            render = F,
            colwidths = c(5, 3, 0.8, 1.2, 1.2)
          )
        
        
        plot(fgsea_gtable)
        
        addSubset = paste(i, 'fgsea_results', sep = '_')
        return_list[[addSubset]] <- fgsea_results
        
        addSubset = paste(i, 'fgsea_gtable', sep = '_')
        return_list[[addSubset]] <- fgsea_gtable
        
        addSubset = paste(i, 'fgsea_ranks', sep = '_')
        return_list[[addSubset]] <- ranks
        
        addSubset = 'msig_geneSet'
        return_list[[addSubset]] <- msig_geneSet
        
      }
    })
  }
  return(return_list)
}




as.enrichResult_internal <- function(result, inputIds, geneSet) {
  if (is.null(inputIds)) {
    print('null geneIds')
  }
  gene <- inputIds
  gene.length <- length(gene)
  
  result <- result %>% 
    dplyr::rename('Count' = size, 'p.adjust' = padj, 'pvalue' = pval, Description = 'pathway') 
    # %>% dplyr::arrange(dplyr::desc(p.adjust))
  
  result <- result[order(pvalue),]
  
  result$GeneRatio <- paste(result$Count, '/', gene.length, sep = '')
  result$size <- result$Count
  result$ID <- result$Description
  
  rownames(result) <- result$Description
  
  #result$qvalue <- result$pvalue
  geneSetsOI <- geneSet[c(result$Description)]
  genesInGeneSet <- lapply(geneSetsOI, intersect, y=gene)
  genesInGeneSet.stack <- stack(genesInGeneSet) %>% 
    rename(ind = 'Description') %>% group_by(Description) %>% 
    summarise(geneID = paste(values, collapse = '/'))
  result <- merge(result, genesInGeneSet.stack, by = 'Description')
  
  new('enrichResult',
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

as.enrichResult <- function(result, inputIds, geneSet) {
  e <- as.enrichResult_internal(result = result, inputIds = inputIds, geneSet = geneSet)
  rownames(e@result) <- e@result$Description
  
  return(e)
}


renderPlotSet <- function(output, key, enrichTypeResult) {
  #print(paste0('called render plot: ', is.null(enrichTypeResult)))
  
  ##add if is not null exception
  output[[paste(key, 'table', sep = '_')]] <- renderDataTable({
    req(!is.null(enrichTypeResult()))
    enrichTypeResult() %>% as.data.frame() %>% dplyr::rename(
      'Term Description' = Description,
      'geneID' = geneID,
      'Hits' = Count,
      'p-Value (adj.)' = pvalue,
      'p-Value' = p.adjust,
      #'q-Value' = qvalue
    ) %>% 
      dplyr::select(c('Term Description', 'Hits', 'p-Value (adj.)', 'p-Value', 'geneID', dplyr::everything())) %>% 
      DT::datatable(
        #table,
        filter = 'bottom',
        selection = 'single',
        escape = FALSE,
        autoHideNavigation = TRUE,
        rownames = FALSE,
        extensions = c('Buttons'),
        class = 'cell-border stripe',
        options = list(
          dom = 'Bfrtip',
          lengthMenu = c(15, 30, 50, 100),
          pageLength = 10,
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
  })
  
  output[[paste(key, 'dotplot', sep = '_')]] <- renderPlotly({
    req(!is.null(enrichTypeResult()))
    enrichplot::dotplot(enrichTypeResult())
  })
  
  output[[paste(key, 'emapplot', sep = '_')]] <- renderPlot({
    req(!is.null(enrichTypeResult()))
    enrichplot::emapplot(enrichTypeResult())
  })
  
  output[[paste(key, 'cnetplot', sep = '_')]] <- renderPlot({
    req(!is.null(enrichTypeResult()))
    enrichplot::cnetplot(enrichTypeResult())
  })
  
  output[[paste(key, 'upsetplot', sep = '_')]] <- renderPlot({
    req(!is.null(enrichTypeResult()))
    enrichplot::upsetplot(enrichTypeResult())
  })
  
  output[[paste(key, 'heatplot', sep = '_')]] <- renderPlot({
    req(!is.null(enrichTypeResult()))
    enrichplot::heatplot(enrichTypeResult())
  })
}

makeTabBox <- function(title, key) {
  tabBox(
    title = title,   
    side = 'right',
    height = NULL,
    selected = 'Dot Plot',
    width = 16,
    tabPanel('Table', dataTableOutput(paste(key, 'table', sep = '_'))),
    tabPanel('Dot Plot', plotlyOutput(paste(key, 'dotplot', sep = '_'))),
    tabPanel('Emap Plot', plotOutput(paste(key, 'emapplot', sep = '_'))),
    tabPanel('Cnet Plot', plotOutput(paste(key, 'cnetplot', sep = '_'))),
    tabPanel('Upset Plot', plotOutput(paste(key, 'upsetplot', sep = '_'))),
    tabPanel('Heat Plot', plotOutput(paste(key, 'heatplot', sep = '_')))
  )
}


