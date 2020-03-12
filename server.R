
source('fxs.R', local = TRUE)
source('info.R', local = TRUE)

server = function(input, output, session) {
  
  options(shiny.maxRequestSize=50*1024^2) 
  preferences <- reactiveValues(use_webgl = TRUE)
  
  
  ##-------------------##
  ## Sample data
  ##-------------------##
  areaTable <- NULL
  areaTable <- eventReactive(input$submit, {
    validate(need(input$areaInput_runname != '', "Please type in a run name..."))
    areaTable <- read.table(text = gsub("(?<=[a-z])\\s+", "\n", perl = TRUE, x = input$areaInput),
        header = FALSE,
        col.names = c("gene", "avg_logFC"),
        quote = "",
        allowEscapes = T)
    saveRDS(areaTable, file = paste0('SavedData/', input$areaInput_runname, '.rds', sep = ''))
    #updateTextAreaInput(session, "areaInput", value = paste("HMOX 2.00 \nABCA4 -1.50", x))
    areaTable
  })
  
  output$inputTable <- renderTable({
    req(input$submit)
    areaTable()
  })
  
  observeEvent(input$input_file, {
    output$inputTable <- renderTable({
      req(input$input_file)
      data.table::data.table(readRDS(input$input_file$datapath))
    })
  })
  
  sample_data_path <- reactive({
    areaTableFile <- paste0('SavedData/', input$areaInput_runname, '.rds')
    if ( !file.exists(areaTableFile) ) {
      if ( is.null(input$input_file) || is.na(input$input_file) ) {
        sample_data_path <- 'testdata/cbndObj02_x7.res02.markers.rds'
      } else {
        req(input$input_file)
        sample_data_path <- input$input_file$datapath
      }
    } else {
      sample_data_path <- areaTableFile
    }
    sample_data_path
  })
  
  sample_data <- reactive({readRDS(sample_data_path())})
  
  data_basename <- reactive({
    data_basename <- basename(sample_data_path()) %>% tools::file_path_sans_ext()
  })
  
  ##load results if ran previously
  string_results <- reactive({
    string_result_rds_name <- paste(data_basename(), 'string', 'result', sep = '_')
    string_result_rds_path <- paste('SavedRuns/', string_result_rds_name, '.rds', sep = '')
    if (file.exists(string_result_rds_path)) {
      string_results <- readRDS(string_result_rds_path)
    } else {
      string_results <- NULL
      print('RunSTRINGdb not yet run on input...')
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
      msig_results <- NULL
      print('RunMSigDB not yet ran on input...')
    }
    msig_results
  })
  
  msig_results_enricher <- reactive({
    toSubset <- paste(input$msigdbr_select_cluster_input, 'enricher_result', sep = '_')
    msig_results()[[toSubset]]
  })
  
  msig_results_fgsea <- reactive({
    toSubset <- paste(input$msigdbr_select_cluster_input, 'fgsea_results', sep = '_')
    msig_results_fgsea <- msig_results()[[toSubset]]
    
    req(!is.null(msig_results_enricher()))
    as.enrichResult(result = msig_results_fgsea,
                    inputIds = msig_results_enricher()@gene,
                    geneSet = msig_results_enricher()@geneSets)
  })
  
  
  ##--------------------##
  ## stringdb
  ##--------------------##
  stringdbSpecies <- STRINGdb::get_STRING_species(version = '10')
  output$runstringdb_select_parameters <- renderUI({
    req(stringdbSpecies)
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
      )
    )
  })
  
  string_results <- eventReactive(input$runstringdb_button, {
    stringdbSpecies <- STRINGdb::get_STRING_species(version = '10')
    validate(need(input$stringdb_maxHitsToPlot_input != '', "Please type in maxHitsToPlot..."))
    validate(need(input$stringdb_scoreThreshold_input != '', "Please type in scoreThreshold..."))
    validate(need(input$stringdb_refSpecies_input != '', "Please type in refSpecies..."))
    refSpeciesNum = stringdbSpecies$species_id[stringdbSpecies$compact_name == input$stringdb_refSpecies_input]
    
    withProgress(message = 'making STRING query...', {
    stringRes <- runSTRINGdb(DEtable = sample_data(), 
                             maxHitsToPlot = input$stringdb_maxHitsToPlot_input, 
                             refSpeciesNum = refSpeciesNum, 
                             scoreThreshold = input$stringdb_scoreThreshold_input)
    })
    saveRDS(stringRes, paste0('SavedRuns/', data_basename(), '_string_result', '.rds', sep = ''))
    
    stringRes
  })
  
  output$stringdb_select_run_UI <- renderUI({
    fileInput(
      'stringdb_select_run',
      'Load Run...',
      multiple = F,
      accept = c('.rds')
    )
  })
  
  output$stringdb_select_cluster_UI <- renderUI({
    if (is.null(sample_data()$cluster)){
      choices = "NA"
    } else {
      choices = levels(sample_data()$cluster)
    }
    selectInput(
      inputId = 'stringdb_select_cluster_input',
      label = 'Select group (or cluster)',
      choices = choices,
      selected = choices[1]
    )
  })
  
  output$num_of_mapped <- renderValueBox({
    req(string_results())
    extract <- paste(input$stringdb_select_cluster_input, 'hits', sep = '_')
    box(
      title = 'Number of genes mapped',
      width = 6,
      background = 'light-blue',
      sum(!is.na(string_results()[[extract]]))
    )
  })
  
  output$num_of_total_genes <- renderValueBox({
    req(string_results())
    extract <- paste(input$stringdb_select_cluster_input, 'hits', sep = '_')
    box(
      title = 'Number of genes total',
      width = 6,
      background = 'light-blue',
      length(string_results()[[extract]])
    )
  })
  
  output$stringdb_network <- renderPlot({
    validate(need(!is.null(string_results()), "Please Run STRINGdb on input..."))
    extract <- paste(input$stringdb_select_cluster_input, 'network', sep = '_')
    string_results()[[extract]]
  })
  
  output$stringdb_network_png <- renderImage(deleteFile = F, {
    validate(need(!is.null(string_results()), "Please Run STRINGdb on input..."))
    png_file <- paste(input$stringdb_select_cluster_input, '_network', '.png', sep = '')
    list(src = paste('SavedData/', png_file, sep = ''), height='90%', width='90%')
  })
  
  # TODO: download entire dataset
  output$stringdb_GO <- renderDataTable({
    validate(need(!is.null(string_results()), "Please Run STRINGdb on input..."))
    extract <- paste(input$stringdb_select_cluster_input, 'GO', sep = '_')
    string_results()[[extract]] %>% dplyr::rename(
      'Term Description' = term_description,
      'Term ID' = term_id,
      'Proteins' = proteins,
      'Hits' = hits,
      'p-Value' = pvalue,
      'p-Value (adj.)' = pvalue_fdr
    ) %>% 
      dplyr::select(c('Term Description', 'Term ID', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', dplyr::everything())) %>%
      DT::datatable(
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
  
  # TODO: download entire dataset
  output$stringdb_KEGG <- renderDataTable({
    validate(need(!is.null(string_results()), "Please Run STRINGdb on input..."))
    extract <- paste(input$stringdb_select_cluster_input, 'KEGG', sep = '_')
    string_results()[[extract]] %>% dplyr::rename(
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
  
  observeEvent(input[["stringdb_resource_info"]], {
    showModal(
      modalDialog(
        stringdb_resource_info[["text"]],
        title = stringdb_resource_info[["title"]],
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  ##--------------------##
  ## msigdbr
  ##--------------------##
  output$runmsigdbr_select_parameters <- renderUI({
    msigdbrSpecies <- data.frame(Species = msigdbr::msigdbr_show_species())
    species.msig <- msigdbr::msigdbr(species = 'Homo sapiens')  ##similary category in all species so hardcoded for now
    req(species.msig)
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
    species.msig <- msigdbr::msigdbr(species = 'Homo sapiens')
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
  
  msig_results <- eventReactive(input$runmsigdbr_button, {
    req(input$msigdbr_species_input)
    
    withProgress(message = 'making MSigDB query..', {
    msigdbrRes <- runMSigDB(DEtable = sample_data(), species = input$msigdbr_species_input)
    saveRDS(msigdbrRes, paste0('SavedRuns/', data_basename(), '_msig_result', '.rds', sep = ''))
    })
    msigdbrRes
  })
  
  output$msigdbr_select_run_UI <- renderUI({
    fileInput(
      inputId = 'msigdbr_select_run',
      label = 'Load Run...',
      multiple = F,
      accept = c('.rds')
    )
  })
  
  output$msigdbr_select_cluster_UI <- renderUI({
    if (is.null(sample_data()$cluster)){
      choices = "NA"
    } else {
      choices = levels(sample_data()$cluster)
    }
    selectInput(
      inputId = 'msigdbr_select_cluster_input',
      label = 'Select group (or cluster)',
      choices = choices,
      selected = choices[1]
    )
  })
  
  
  output$num_of_mapped_enricher <- renderValueBox({
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
  
  
  output$msigdbr_select_cluster_fgsea_gtable <- renderPlot({
    toSubset <- paste(input$msigdbr_select_cluster_input, 'fgsea_gtable', sep = '_')
    msig_results <- msig_results()
    plot(msig_results[[toSubset]])
  })
  
  
  renderPlotSet(output = output, key = 'fgsea', enrichTypeResult = msig_results_fgsea)
  renderPlotSet(output = output, key = 'enricher', enrichTypeResult = msig_results_enricher)
  
  
  output[["fgsea_table_PPI"]] <- renderPlot({
    req(input$fgsea_table_cell_clicked)
    info_list <- input$fgsea_table_cell_clicked
    print(info_list[["value"]])
    
    toSubset <- paste(input$msigdbr_select_cluster_input, 'fgsea_ranks', sep = '_')
    
    plotEnrichment(msig_results()[['msig_geneSet']][[info_list[["value"]]]], msig_results()[[toSubset]]) + 
      labs(title=info_list[["value"]])
  })
  
  observeEvent(input$fgsea_table_cell_clicked, {
    showModal(
      modalDialog(
        plotOutput("fgsea_table_PPI")
      )
    )
  })
  
  observeEvent(input[["msigdbr_resource_info"]], {
    showModal(
      modalDialog(
        msigdbr_resource_info[["text"]],
        title = msigdbr_resource_info[["title"]],
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  
  ##--------------------##
  ## clusterprofiler
  ##--------------------##
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
  
  observeEvent(input[["clusterprofiler_resource_info"]], {
    showModal(
      modalDialog(
        clusterprofiler_resource_info[["text"]],
        title = clusterprofiler_resource_info[["title"]],
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  
  #*********************************************************
}

