
source('fxs.R', local = TRUE)

server = function(input, output, session) {
  
  preferences <- reactiveValues(use_webgl = TRUE)
  
  ##---------------------------------------------##
  ## Sidebar menu
  ##---------------------------------------------##
  output$sidebar_menu <- renderMenu({
    sidebarMenu(
      id = 'sidebar',
      menuItem(
        text = 'Load data',
        tabName = 'loadData',
        icon = icon(NULL),
        selected = TRUE
      ),
      menuItem(
        text = 'STRINGdb',
        tabName = 'stringdb',
        icon = icon(NULL),
        selected = TRUE
      ),
      menuItem(
        text = 'MsigDB',
        tabName = 'msigdbr',
        icon = icon(NULL)
      ), 
      menuItem(
        text = 'clusterProfiler',
        tabName = 'clusterprofiler',
        icon = icon(NULL)
      )
    )
  })
  
  
  output$runstringdb_select_parameters <- renderUI({
    stringdbSpecies <- STRINGdb::get_STRING_species(version = '10')
    fluidRow(
      selectInput(
        inputId = 'stringdb_refSpecies_input',
        label = 'Reference species:',
        selected = 'Homo sapiens',
        choices = stringdbSpecies$compact_name
      ),
      numericInput(
        inputId = 'stringdb_maxHitsToPlot_input',
        label = 'Max. number of hits to plot:',
        value = 200
      ), 
      numericInput(
        inputId = 'stringdb_scoreThreshold_input',
        label = 'Score threshold:',
        value = 0
      ), 
      textInput(
        inputId = 'stringdb_runname',
        label = 'Run name:',
        placeholder = 'Run1'
        
      )
    )
  })
  
  eventReactive(input$runstringdb_button, {
    req(input$stringdb_select_cluster_input)
    req(input$stringdb_maxHitsToPlot_input)
    req(input$stringdb_scoreThreshold_input)
    req(input$stringdb_refSpecies_input)
    refSpeciesNum = stringdbSpecies$species_id[stringdbSpecies$compact_name == input$stringdb_refSpecies_input]
    
    stringRes <- runSTRINGdb(DEtable = sample_data, 
                             maxHitsToPlot = input$stringdb_maxHitsToPlot_input, 
                             refSpeciesNum = refSpeciesNum, 
                             scoreThreshold = input$stringdb_scoreThreshold_input)
    saveRDS(stringRes, paste0('SavedRuns/', input$stringdb_runname, '_string', '.rds', sep = ''))
  })
  
  
  output$runmsigdbr_select_parameters <- renderUI({
    msigdbrSpecies <- data.frame(Species = msigdbr::msigdbr_show_species())
    species.msig <- msigdbr::msigdbr(species = 'Homo sapiens')
    
    fluidRow(
      selectInput(
        inputId = 'msigdbr_species_input',
        label = 'Reference species',
        selected = 'Homo sapiens',
        choices = msigdbrSpecies$Species
      ),
      selectInput(
        inputId = 'msigdbr_category_input',
        label = 'Select category (optional):',
        selected = '',
        choices = c('', unique(species.msig$gs_cat))
      )
    )
  })
  
  observeEvent(input$msigdbr_category_input, { 
    subcat <- filter(species.msig, species.msig$gs_cat == input$msigdbr_category_input)
    output$runmsigdbr_select_parameters_sub <- renderUI({
      selectInput(
        inputId = 'msigdbr_subcategory_input',
        label = 'Select subcategory (optional):',
        selected = '',
        choices = c('', unique(subcat$gs_subcat))
      )
    })
  })
  
  output$runmsigdbr_select_parameters_cont <- renderUI({
    textInput(inputId = 'msigdbr_runname',
              label = 'Run name',
              placeholder = 'Run1')
  })
  
  eventReactive(input$runmsigdbr_button, {
    req(input$msigdbr_species_input)
    req(input$msigdbr_runname)
    
    msigdbrRes <- runMSigDB(DEtable = sample_data, species = input$msigdbr_species_input)
    saveRDS(msigdbrRes, paste0('SavedRuns/', input$msigdbr_runname, '_msigdbr', '.rds', sep = ''))
  })
  
  
  ##-------------------##
  ## Sample data
  ##-------------------##
  data_rds_path <- reactive({
    if (is.null(input$input_file) || is.na(input$input_file)) {
      data_rds_path <- 'testdata/cbndObj02_x7.res02.markers.rds'
    } else {
      req(input$input_file)
      data_rds_path <- input$input_file
    }
    data_rds_path
  })
  
  sample_data <- reactive({
    sample_data <- readRDS(data_rds_path())
  })
  
  data_basename <- reactive({
    data_basename <- basename(data_rds_path()) %>% tools::file_path_sans_ext()
  })
  
  string_results <- reactive({
    string_result_rds_name <- paste(data_basename(), 'string', 'result', sep = '_')
    string_result_rds_path <- paste('SavedRuns/', string_result_rds_name, '.rds', sep = '')
    if (file.exists(string_result_rds_path)) {
      string_results <- readRDS(string_result_rds_path)
    } else {
      print('RunSTRINGdb not yet ran on input...')
    }
    string_results
  })
  
  msig_results <- reactive({
    msig_result_rds_name <- paste(data_basename(), 'msig', 'result', sep = '_')
    msig_result_rds_path <- paste('SavedRuns/', msig_result_rds_name, '.rds', sep = '')
    if (file.exists(msig_result_rds_path)) {
      msig_results <- readRDS(msig_result_rds_path)
      print("loading msig_result")
    } else {
      print('RunMSigDB not yet ran on input...')
    }
    msig_results
  })
  
  msig_results_enricher <- reactive({
    toSubset <- paste(input$msigdbr_select_cluster_input, 'enricher_result', sep = '_')
    msig_results_enricher <- msig_results()[[toSubset]]
    msig_results_enricher
  })

  msig_results_fgsea <- reactive({
    toSubset <- paste(input$msigdbr_select_cluster_input, 'fgsea_results', sep = '_')
    msig_results_fgsea <- msig_results()[[toSubset]]
    msig_results_fgsea
  })
    

  ##--------------------##
  ## Tabs - server.R
  ##--------------------##
  #######################################################################################################################    
  output$stringdb_select_run_UI <- renderUI({
    if ( is.null(sample_data()$cluster) ) {
      #textOutput(NULL)
      #textOutput('enriched_pathways_by_sample_table_missing')
    } else {
      fileInput('stringdb_select_run', 'Load Run...', multiple = F, accept = c('.rds'))
    }
  })
  
  output$stringdb_select_cluster_UI <- renderUI({
    if ( is.null(sample_data()$cluster) ) {
      #textOutput(NULL)
      #textOutput('enriched_pathways_by_sample_table_missing')
    } else {
      selectInput(
        inputId = 'stringdb_select_cluster_input',
        label = 'Select group (or cluster)',
        choices = levels(sample_data()$cluster)
      )
    }
  })
  
  output$num_of_mapped <- renderValueBox({
    string_results <- string_results()
    extract <- paste(input$stringdb_select_cluster_input, 'hits', sep = '_')
    box(
      title = 'Number of genes mapped',
      width = 6,
      background = 'light-blue',
      sum(!is.na(string_results[[extract]]))
    )
  })
  
  output$num_of_total_genes <- renderValueBox({
    string_results <- string_results()
    extract <- paste(input$stringdb_select_cluster_input, 'hits', sep = '_')
    box(
      title = 'Number of genes total',
      width = 6,
      background = 'light-blue',
      length(string_results[[extract]])
    )
  })
  
  #TODO: pval and pct.1 filter input
  output$stringdb_network <- renderPlot({
    # req(input$stringdb_select_cluster_input)
    string_results <- string_results()
    extract <- paste(input$stringdb_select_cluster_input, 'network', sep = '_')
    string_results[[extract]]
  })
  
  output$stringdb_network_png <- renderImage(deleteFile = F, {
    png_file <- paste(input$stringdb_select_cluster_input, '_network', '.png', sep = '')
    list(src = paste('SavedData/', png_file, sep = ''), height='90%', width='90%')
    #img(src = png_file, height='60%', width='60%', align='right')
  })
  
  output$stringdb_GO <- renderDataTable({
    string_results <- string_results()
    extract <- paste(input$stringdb_select_cluster_input, 'GO', sep = '_')
    table <- string_results[[extract]]
    table <- table %>% dplyr::rename(
      'Term Description' = term_description,
      'Term ID' = term_id,
      'Proteins' = proteins,
      'Hits' = hits,
      'p-Value' = pvalue,
      'p-Value (adj.)' = pvalue_fdr
    ) %>% 
      dplyr::select(c('Term Description', 'Term ID', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', dplyr::everything())) %>%
      DT::datatable(
        #table,
        filter = 'bottom',
        selection = 'multiple',
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
      )
  })
  
  output$stringdb_select_GO_ann <- renderUI({
    string_results <- string_results()
    extract <- paste(input$stringdb_select_cluster_input, 'GO', sep = '_')
    selectInput(
      inputId = 'stringdb_select_GO_input',
      label = NULL,
      choices = c('', levels(as.factor(string_results[[extract]][['term_id']]))), 
      selected = NULL, 
      multiple = F
    )
  })
  
  output$stringdb_select_GO_ann_output <- renderText({
    req(input$stringdb_select_GO_input)
    string_results <- string_results()
    extract <- paste(input$stringdb_select_cluster_input, 'hits', sep = '_')
    GO_hits(GO_term = input$stringdb_select_GO_input, hits = string_results[[extract]])
  })
  
  # extract <- paste(input$stringdb_select_cluster_input, 'hits', sep = '_')
  # reactGO_hits <- reactive(GO_hits(GO_term = input$stringdb_select_GO_input, hits = string_results[[extract]]))
  # output$stringdb_select_GO_ann_output <- renderText({
  #   observeEvent(input$stringdb_select_GO_input, {
  #     reactGO_hits()
  #   })}
  # )
  
  stringdb_GO_resource_info <- list(
    title = "STRING Resource info",
    text = HTML(
      '<b>STRING</b><br>
      STRING is a database of known and predicted protein-protein interactions. \n
      The interactions include direct (physical) and indirect (functional) associations; they stem from computational prediction, \n
      from knowledge transfer between organisms, and from interactions aggregated from other (primary) databases.
      <p>
    <li><a href=https://string-db.org/ title="Official string-db website" target="_blank"><b>Official string-db website</b></a></li>
    '
    )
  )
  
  # info box
  observeEvent(input[["stringdb_GO_resource_info"]], {
    showModal(
      modalDialog(
        stringdb_GO_resource_info[["text"]],
        title = stringdb_GO_resource_info[["title"]],
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  
  ##--------------------##
  # TODO: download entire dataset
  output$stringdb_KEGG <- renderDataTable({
    string_results <- string_results()
    extract <- paste(input$stringdb_select_cluster_input, 'KEGG', sep = '_')
    table <- string_results[[extract]]
    table <- table %>% dplyr::rename(
      'Term Description' = term_description,
      'Term ID' = term_id,
      'Proteins' = proteins,
      'Hits' = hits,
      'p-Value' = pvalue,
      'p-Value (adj.)' = pvalue_fdr
    ) %>% 
      dplyr::select(c('Term Description', 'Term ID', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', dplyr::everything())) %>%
      DT::datatable(
        #table,
        filter = 'bottom',
        selection = 'multiple',
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
      )
  })
  
  
  output$msigdbr_select_run_UI <- renderUI({
    if ( is.null(sample_data()$cluster) ) {
      #textOutput(NULL)
      #textOutput('enriched_pathways_by_sample_table_missing')
    } else {
      fileInput(inputId = 'msigdbr_select_run', label = 'Load Run...', multiple = F, accept = c('.rds'))
    }
  })
  
  
  output$msigdbr_select_cluster_UI <- renderUI({
    if ( is.null(sample_data()$cluster) ) {
      #textOutput(NULL)
      #textOutput('enriched_pathways_by_sample_table_missing')
    } else {
      selectInput(
        inputId = 'msigdbr_select_cluster_input',
        label = 'Select group (or cluster)',
        choices = levels(sample_data()$cluster)
      )
    }
  })
  
  
  #TODO: pval and pct.1 filter input
  output$msigdbr_select_cluster_enricher_table <- renderDataTable({
    table <- msig_results_enricher() %>% as.data.frame() %>% dplyr::rename(
      'Term Description' = Description,
      'geneID' = geneID,
      'Hits' = Count,
      'p-Value (adj.)' = pvalue,
      'p-Value' = p.adjust,
      'q-Value' = qvalue
    ) %>% 
      dplyr::select(c('Term Description', 'Hits', 'p-Value (adj.)', 'p-Value', 'geneID', dplyr::everything())) %>% 
      DT::datatable(
        #table,
        filter = 'bottom',
        selection = 'multiple',
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
    table
  })
  
  
  output$msigdbr_select_cluster_enricher_plot <- renderPlot({
    e_plot <- enrichplot::dotplot(msig_results_enricher())
    e_plot
  })
  
  
  output$num_of_mapped_enricher <- renderValueBox({
    #msig_results_enricher <- msig_results_enricher()
    num_genes_mapped <- str_split(noquote(msig_results_enricher()@result$GeneRatio[1]), '/')[[1]][2]
    box(
      title = 'Number of genes mapped',
      width = 6,
      background = 'light-blue',
      num_genes_mapped
    )
  })
  
  
  output$num_of_total_genes_enricher <- renderValueBox({
    num_genes_total <- length(msig_results_enricher()@gene)
    box(
      title = 'Number of genes total',
      width = 6,
      background = 'light-blue',
      num_genes_total
    )
  })
  
  
  output$msigdbr_select_cluster_fgsea_table <- renderDataTable({
    msig_results_fgsea <- msig_results_fgsea() %>% as.data.frame() %>% dplyr::rename(
      'Term Description' = pathway,
      'Hits' = size,
      'p-Value' = pval,
      'p-Value (adj.)' = padj
    ) %>% 
      dplyr::select(c('Term Description', 'Hits', 'p-Value (adj.)', 'p-Value', dplyr::everything())) %>% 
      DT::datatable(
        #table,
        filter = 'bottom',
        selection = 'multiple',
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
      )
  })
  
  
  output$msigdbr_select_cluster_fgsea_gtable <- renderPlot({
    toSubset <- paste(input$msigdbr_select_cluster_input, 'fgsea_gtable', sep = '_')
    msig_results <- msig_results()
    plot(msig_results[[toSubset]])
  })
  
  
  output$msigdbr_resource_info <- renderText({
    'The Molecular Signatures Database (MSigDB) is a collection of gene sets
  originally created for use with the Gene Set Enrichment Analysis (GSEA) software.
  Gene homologs are provided by HUGO Gene Nomenclature Committee at the European Bioinformatics Institute
  which integrates the orthology assertions predicted for human genes by
  eggNOG, Ensembl Compara, HGNC, HomoloGene, Inparanoid, NCBI Gene Orthology, OMA, OrthoDB, OrthoMCL, Panther, PhylomeDB, TreeFam and ZFIN.
  For each human equivalent within each species, only the ortholog supported by the largest number of databases is used.
  '
  })
  
  output$david_select_cluster_UI <- renderUI({
    if ( is.null(sample_data()$cluster) ) {
      #textOutput(NULL)
      #textOutput('enriched_pathways_by_sample_table_missing')
    } else {
      selectInput(
        inputId = 'david_select_cluster_input',
        label = NULL,
        choices = levels(sample_data()$cluster)
      )
    }
  })
  
  
  msig_results_fgsea_asEnricher <- reactive({
    msig_results_fgsea_asEnricher <- as.enrichResult(result = msig_results_fgsea(),
                                                  inputIds = msig_results_enricher()@gene,
                                                  geneSet = msig_results_enricher()@geneSets)
    msig_results_fgsea_asEnricher
  })
  
  
  renderPlotSet(output = output, key = 'fgsea', enrichTypeResult = msig_results_fgsea_asEnricher())
  
  renderPlotSet(output = output, key = 'enricher', enrichTypeResult = msig_results_enricher())
  
  output[["fgsea_table_PPI"]] <- renderPlot({
    info_list <- input$fgsea_table_cell_clicked
    print(info_list[["value"]])
    
    toSubset <- paste(input$msigdbr_select_cluster_input, 'fgsea_ranks', sep = '_')
    
    plotEnrichment(msig_results()[['msig_geneSet']][[info_list[["value"]]]], msig_results()[[toSubset]]) + 
      labs(title=info_list[["value"]])
  })
  
  
  observeEvent(input$fgsea_table_cell_clicked, {
    showModal(
      modalDialog(
        # This plot does display
        plotOutput("fgsea_table_PPI")
      )
    )
  })
  
  
  #################################################################################################################
  #################################################################################################################
  output$clusterprofiler_select_run_UI <- renderUI({
    if ( is.null(sample_data()$cluster) ) {
      #textOutput(NULL)
      #textOutput('enriched_pathways_by_sample_table_missing')
    } else {
      fileInput(inputId = 'clusterprofiler_select_run', label = 'Load Run...', multiple = F, accept = c('.rds'))
    }
  })
  
  
  output$clusterprofiler_select_cluster_UI <- renderUI({
    if ( is.null(sample_data()$cluster) ) {
      #textOutput(NULL)
      #textOutput('enriched_pathways_by_sample_table_missing')
    } else {
      selectInput(
        inputId = 'clusterprofiler_select_cluster_input',
        label = 'Select group (or cluster)',
        choices = levels(sample_data()$cluster)
      )
    }
  })
  
  
  #******************************************************************************************************************
}

