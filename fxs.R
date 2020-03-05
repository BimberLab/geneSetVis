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
      version = "10",
      species = refSpeciesNum,
      score_threshold = scoreThreshold,
      input_directory = ""
    )
  
  ##dedup table to remove multiple tests
  DEtable.dedup <-
    DEtable[with(DEtable, order(p_val_adj, decreasing = F)),]
  DEtable.dedup <-
    DEtable.dedup[match(unique(DEtable.dedup$gene), DEtable.dedup$gene), ]
  
  DEtable.split <- split(DEtable.dedup, DEtable.dedup$cluster)
  
  return_list = c()
  for (i in as.vector(names(DEtable.split))) {
    tryCatch({
      clusterTable <- DEtable.split[[i]]
      
      if (nrow(DEtable.split[[i]]) > 0) {
        print(string_db$get_link(hits))
        cluster.map <-
          string_db$map(clusterTable, "gene", removeUnmappedRows = FALSE)
        hits <- cluster.map$STRING_id
        if (sum(is.na(hits)) > 0){
          stop()
        }
        
        max_hits_to_plot <- cluster.map$STRING_id[1:maxHitsToPlot]
        print(string_db$get_link(hits))
        
        enrichmentGO <-
          string_db$get_enrichment(hits,
                                   category = "Process",
                                   methodMT = "fdr",
                                   iea = TRUE)
        print(string_db$get_link(hits))
        enrichmentKEGG <-
          string_db$get_enrichment(hits,
                                   category = "KEGG",
                                   methodMT = "fdr",
                                   iea = TRUE)
        print(string_db$get_link(hits))
        
        
        hit_term_proteins <- string_db$get_term_proteins(enrichmentGO$term_id, hits)
        hit_term_genes <- hit_term_proteins %>% 
          dplyr::select(term_id, preferred_name) %>% 
          dplyr::group_by(term_id) %>% 
          dplyr::summarize("hit_term_genes" = paste0(preferred_name, collapse = ", "))
        
        enrichmentGO <- merge(hit_term_genes, enrichmentGO)
        
        
        hit_term_proteins <- string_db$get_term_proteins(enrichmentKEGG$term_id, hits)
        hit_term_genes <- hit_term_proteins %>% 
          dplyr::select(term_id, preferred_name) %>% 
          dplyr::group_by(term_id) %>% 
          dplyr::summarize("hit_term_genes" = paste0(preferred_name, collapse = ", "))
        
        enrichmentKEGG <- merge(hit_term_genes, enrichmentKEGG)
        
        
        link_interactions <- string_db$get_link(hits)
        
        string_db$get_png(max_hits_to_plot, file = paste(i, "network.png", sep = "_"))
        
        
        string_db$plot_network(max_hits_to_plot)
        
        #______
        ##payload mechanism for upregulated vs downregulated genes:
        ##adds a color column for up vs downregulated genes
        cluster.color <-
          string_db$add_diff_exp_color(cluster.map, logFcColStr = "avg_logFC")
        # post payload information to the STRING server
        payload_id <-
          string_db$post_payload(cluster.color$STRING_id, colors = cluster.color$color)
        string_db$plot_network(hits, payload_id = payload_id)
        
        ##clustering/community algorithms: ”fastgreedy”, ”walktrap”, ”spinglass”, ”edge.betweenness”.
        networkClustersList <-
          string_db$get_clusters(max_hits_to_plot, algorithm = "fastgreedy")
        par(mfrow = c(2, 2))
        for (j in seq(1:length(networkClustersList))) {
          string_db$plot_network(networkClustersList[[j]], payload_id = payload_id)
        }
        
        hits_name = paste(i, "hits", sep = "_")
        return_list[[hits_name]] <- hits
        
        network_name = paste(i, "network", sep = "_")
        return_list[[network_name]] <- network
        
        go_name = paste(i, "GO", sep = "_")
        return_list[[go_name]] <- enrichmentGO
        
        kegg_name = paste(i, "KEGG", sep = "_")
        return_list[[kegg_name]] <- enrichmentKEGG
        
      }
      
    }, error = function(e) {
      cat("\nERROR :", conditionMessage(e), "\n")
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
  
  return_list = c()
  for (i in as.vector(names(DEtable.split))) {
    tryCatch({
      clusterTable <- DEtable.split[[i]]
      
      if (nrow(clusterTable) > 0) {
        ##Use the gene sets data frame for clusterProfiler (for genes as gene symbols)
        msig_enricher <- enricher(gene = clusterTable$gene, TERM2GENE = msigTerm)
        msig_enricher_plot <- dotplot(msig_enricher)
        
        #clusterProfiler::geneInCategory()
        #geneInCategory(msig_enricher)[as.data.frame(msig_enricher)$ID == "WINTER_HYPOXIA_METAGENE"][1]
        
        # enricher_KEGG <- enrichKEGG(
        #   clusterTable$gene,
        #   organism = "hsa",
        #   keyType = "kegg",
        #   pAdjustMethod = "BH"
        # )
        
        msig_enricher <- as.data.frame(msig_enricher)
        msig_enricher$geneID <- gsub(x = msig_enricher$geneID, pattern = "/", replacement = ",")
        enricher_name = paste(i, "enricher", sep = "_")
        return_list[[enricher_name]] <- as.data.frame(msig_enricher)
        
        enricher_plot_name = paste(i, "enricher_plot", sep = "_")
        return_list[[enricher_plot_name]] <- msig_enricher_plot
        
        #.......................................
        ##Use the gene sets data frame for fgsea.
        msg_list = human.msig %>% split(x = .$gene_symbol, f = .$gs_name)
        
        ##name the marker genes with their avgLogFC
        ranks <- clusterTable$avg_logFC
        ranks <- setNames(ranks, clusterTable$gene)
        
        set.seed(1234)
        fgsea_results <- fgsea(
          pathways = msg_list,
          stats = ranks,
          minSize = 5,
          maxSize = 600,
          nperm = 10000
        )
        
        threshold <- 0.001
        sigPathways.sum <- sum(fgsea_results[, padj < threshold])
        print(paste0(
          sigPathways.sum,
          " significant pathways. pval < ",
          threshold
        ))
        
        #topPathways <- fgsea_results[head(order(pval), n=15)][order(NES), pathway]
        topPathwaysUp <-
          fgsea_results[ES > 0][head(order(pval), n = 10), pathway]
        topPathwaysDown <-
          fgsea_results[ES < 0][head(order(pval), n = 10), pathway]
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
        
        
        fgsea_gtable <-
          plotGseaTable(
            pathways = msg_list[topPathways],
            stats = ranks,
            fgseaRes = fgsea_results,
            gseaParam = 0.5,
            render = F,
            colwidths = c(5, 3, 0.8, 1.2, 1.2)
          )
        
        
        plot(fgsea_gtable)
        
        fgsea_results_name = paste(i, "fgsea_results", sep = "_")
        return_list[[fgsea_results_name]] <- fgsea_results
        
        fgsea_gtable_name = paste(i, "fgsea_gtable", sep = "_")
        return_list[[fgsea_gtable_name]] <- fgsea_gtable
        
        num_genes_name = paste(i, "num_genes", sep = "_")
        return_list[[num_genes_name]] <- length(clusterTable$gene)
        
        input_genes_name = paste(i, "input_genes", sep = "_")
        return_list[[input_genes_name]] <- clusterTable$gene
        
        # plot the most significantly enriched pathway
        #plotEnrichment(msg_list[[head(fgsea_results[order(pval), ], 1)$pathway]], ranks)
        #+ labs(title=head(fgsea_results[order(pval), ], 1)$pathway)
      }
    })
  }
  return(return_list)
}




as.enrichResult_internal <- function(result, inputIds) {
  
  gene <- inputIds
  gene.length <- length(gene)
  
  result <- result %>% dplyr::rename("Count" = size, "p.adjust" = padj, "pvalue" = pval, Description = "pathway") # %>% dplyr::arrange(dplyr::desc(p.adjust))
  
  result <- result[order(pvalue),]
  
  result$GeneRatio <- paste(result$Count, "/", gene.length, sep = "")
  result$size <- result$Count
  result$ID <- result$Description
  
  rownames(result) <- result$Description
  
  #result$qvalue <- result$pvalue
  geneSetsOI <- msg_list[c(result$Description)]
  genesInGeneSet <- lapply(geneSetsOI, intersect, y=gene)
  genesInGeneSet.stack <- stack(genesInGeneSet) %>% 
    rename(ind = "Description") %>% group_by(Description) %>% 
    summarise(geneID = paste(values, collapse = "/"))
  result <- merge(result, genesInGeneSet.stack, by = "Description")
  
  new("enrichResult",
      result         = result, 
      pvalueCutoff   = 0.05,
      pAdjustMethod  = "UNKNOWN",
      #qvalueCutoff   = 1,
      gene           = as.character(gene),
      #universe       = extID,
      geneSets       = msg_list,
      organism       = "UNKNOWN",
      keytype        = "UNKNOWN",
      ontology       = "UNKNOWN",
      readable       = T
  )
  
}

as.enrichResult <- function(result, inputIds) {
  e <- as.enrichResult_internal(result = result, inputIds = inputIds)
  rownames(e@result) <- e@result$Description
  
  return(e)
}

makePlotSet <- function(input, output, session, outputKey, enrichTypeResList) {
  #output <- list()
  reactive({
  selectedCluster <- input$msigdbr_select_cluster_input
  selectedCluster <- selectedCluster()
  enrichTypeResList <- reactive({enrichTypeResList()})
  extract <- paste(selectedCluster, "fgsea_results", sep = "_")
  extract2 <- paste(selectedCluster, "input_genes", sep = "_")
  
  e <- as.enrichResult(result = enrichTypeResList[[extract]], inputIds = enrichTypeResList[[extract2]])

  #output[[paste0(outputKey, '_table')]] <- codeToRenderTable(...)
  
  output[[paste(outputKey, 'dotplot', sep = "_")]] <- renderPlotly({plotly::ggplotly(enrichplot::dotplot(e))})
  
  output[[paste(outputKey, 'emapplot', sep = "_")]] <- renderPlotly({plotly::ggplotly(enrichplot::emapplot(e))})
  
  output[[paste(outputKey, 'cnetplot', sep = "_")]] <- renderPlotly({plotly::ggplotly(enrichplot::cnetplot(e))})
  
  output[[paste(outputKey, 'upsetplot', sep = "_")]] <- renderPlotly({plotly::ggplotly(enrichplot::upsetplot(e))})
  
  output[[paste(outputKey, 'heatplot', sep = "_")]] <- renderPlotly({plotly::ggplotly(enrichplot::heatplot(e))})
})
}




#makePlotSet('enrichr', ....)

#makePlotSet('fgsea', ..)

