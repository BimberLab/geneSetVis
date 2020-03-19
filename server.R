source('fxs.R', local = TRUE)
source('info.R', local = TRUE)

server = function(input, output, session) {
  
  options(shiny.maxRequestSize=50*1024^2) 
  preferences <- reactiveValues(use_webgl = TRUE)
  
  envir <- reactiveValues(sample_data = NULL)
  ##-------------------##
  ## Sample data
  ##-------------------##
  #sample_data <- NULL
  observeEvent(input$submit, {
    sample_data <- read.table(text = gsub("(?<=[a-z])\\s+", "\n", perl = TRUE, x = input$areaInput),
                            header = FALSE,
                            col.names = c("gene", "avg_logFC"),
                            quote = "",
                            allowEscapes = T)
    saveRDS(sample_data, file = paste0('SavedData/', 'running', '.rds', sep = ''))
    string_result <- NULL
    msig_result <- NULL
    # if (file.exists('SavedRuns/running_string_result.rds')) 
    #   file.remove('SavedRuns/running_string_result.rds')
    # if (file.exists('SavedRuns/running_msig_result.rds')) 
    #   file.remove('SavedRuns/running_msig_result.rds')
    #updateTextAreaInput(session, "areaInput", value = paste("HMOX 2.00 \nABCA4 -1.50", x))
    envir$sample_data <- sample_data
    sample_data
  })
  #*
  output$inputTable <- renderTable({
    req(input$submit)
    envir$sample_data
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
  
  observeEvent(c(input$runstringdb_button, input$submit), {
    #if (input$submit == 0)    return(NULL)
    
    stringdbSpecies <- STRINGdb::get_STRING_species(version = '10')
    validate(need(input$stringdb_maxHitsToPlot_input != '', "Please type in maxHitsToPlot..."))
    validate(need(input$stringdb_scoreThreshold_input != '', "Please type in scoreThreshold..."))
    validate(need(input$stringdb_refSpecies_input != '', "Please type in refSpecies..."))
    refSpeciesNum = stringdbSpecies$species_id[stringdbSpecies$compact_name == input$stringdb_refSpecies_input]
    
    withProgress(message = 'making STRING query...', {
      stringRes <- runSTRINGdb(DEtable = envir$sample_data, 
                               maxHitsToPlot = input$stringdb_maxHitsToPlot_input, 
                               refSpeciesNum = refSpeciesNum, 
                               scoreThreshold = input$stringdb_scoreThreshold_input)
    })
    saveRDS(stringRes, paste0('SavedRuns/', 'running', '_string_result', '.rds', sep = ''))
    envir$string_results <- stringRes
  })
  
  output$string_map_stats <- renderText({
    req(envir$string_results)
    paste0(sum(!is.na(envir$string_results[['hits']])), ' out of ', length(envir$string_results[['hits']]), 'genes were mapped.')
  })
  
  output$num_of_mapped <- flexdashboard::renderValueBox({
    req(envir$string_results)
    shinydashboard::box(
      title = 'Number of genes mapped',
      width = 6,
      background = 'light-blue',
      sum(!is.na(envir$string_results[['hits']]))
    )
  })
  
  output$num_of_total_genes <- flexdashboard::renderValueBox({
    req(envir$string_results)
    extract <- paste('hits', sep = '')
    shinydashboard::box(
      title = 'Number of genes total',
      width = 6,
      background = 'light-blue',
      length(envir$string_results[[extract]])
    )
  })
  
  output$stringdb_network <- renderPlot({
    validate(need(!is.null(envir$string_results), "Please Run STRINGdb on input..."))
    extract <- paste('network', sep = '')
    envir$string_results[[extract]]
  })
  
  output$stringdb_network_png <- renderImage(deleteFile = F, {
    validate(need(!is.null(envir$string_results), "Please Run STRINGdb on input..."))
    png_file <- paste('network', '.png', sep = '')
    list(src = paste('', png_file, sep = ''), height='90%', width='90%')
  })
  
  # TODO: download entire dataset
  output$stringdb_GO <- renderDataTable({
    validate(need(!is.null(envir$string_results), "Please Run STRINGdb on input..."))
    extract <- paste('GO', sep = '')
    table <- envir$string_results[[extract]] %>% dplyr::rename(
      'Term Description' = term_description,
      'Term ID' = term_id,
      'Proteins' = proteins,
      'Hits' = hits,
      'p-Value' = pvalue,
      'p-Value (adj.)' = pvalue_fdr, 
      'Genes in Term' = hit_term_genes
    ) %>% 
      dplyr::select(c('Term Description', 'Term ID', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', dplyr::everything())) 
    
    table$'Term ID' <- hyperlink_text(text = table$'Term ID', url = "https://www.ebi.ac.uk/QuickGO/term/")
    
    #table$'geneID' <- gsub('/', ',', x = table$'geneID')
    table$'Genes in Term' <- multi_hyperlink_text(labels = table$'Genes in Term', links = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=")
    
      DT::datatable(
        table,
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
    validate(need(!is.null(envir$string_results), "Please Run STRINGdb on input..."))
    extract <- paste('KEGG', sep = '')
    table <- envir$string_results[[extract]] %>% dplyr::rename(
      'Term Description' = term_description,
      'Term ID' = term_id,
      'Proteins' = proteins,
      'Hits' = hits,
      'p-Value' = pvalue,
      'p-Value (adj.)' = pvalue_fdr,
      'Genes in Term' = hit_term_genes
    ) %>% 
      dplyr::select(c('Term Description', 'Term ID', 'Proteins', 'Hits', 'p-Value (adj.)', 'p-Value', dplyr::everything())) 
    
    table$'Term ID' <- hyperlink_text(text = table$'Term ID', url = "https://www.genome.jp/dbget-bin/www_bget?map")
    
    #table$'geneID' <- gsub('/', ',', x = table$'geneID')
    table$'Genes in Term' <- multi_hyperlink_text(labels = table$'Genes in Term', links = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=")
    
      DT::datatable(
        table,
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
  
  observeEvent(c(input$runmsigdbr_button, input$submit), {
    #if (is.null(msig_result)) {}
    # if (input$submit == 0)    return(NULL)
    # msig_result <- NULL
    req(input$msigdbr_species_input)
    
    withProgress(message = 'making MSigDB query..', {
      msigdbrRes <- runMSigDB(DEtable = envir$sample_data, species = input$msigdbr_species_input)
      #saveRDS(msigdbrRes, paste0('SavedRuns/', 'running', '_msig_result', '.rds', sep = ''))
    })
    envir$msig_result <- msigdbrRes
    
    
    envir$msig_result_enricher <- msigdbrRes[['enricher_result']]
    
    
    envir$msig_result_fgsea <- as.enrichResult(result = msigdbrRes[['fgsea_result']],
                                                inputIds = msigdbrRes[['enricher_result']]@gene,
                                                geneSet = msigdbrRes[['enricher_result']]@geneSets)
    
    renderPlotSet(output = output,
                  key = 'fgsea',
                  enrichTypeResult = envir$msig_result_fgsea)
    
    renderPlotSet(output = output,
                  key = 'enricher',
                  enrichTypeResult = envir$msig_result_enricher)
    #}
    
    # observe({
    #   #if (req(input$submit) == 0) {
    #     #updateSelectInput(...)
    #     renderPlotSet(output = output,
    #                   key = 'fgsea',
    #                   enrichTypeResult = envir$msig_result_fgsea)
    #     
    #     renderPlotSet(output = output,
    #                   key = 'enricher',
    #                   enrichTypeResult = envir$msig_result_enricher)
    #   #}
    # })

    # renderPlotSet(output = output,
    #               key = 'fgsea',
    #               enrichTypeResult = envir$msig_result_fgsea)
    # 
    # renderPlotSet(output = output,
    #               key = 'enricher',
    #               enrichTypeResult = envir$msig_result_enricher)
    # 
    


    ##update msig_result()
  })
  
  # print(paste0('envir$msig_result: \n', head(envir$msig_result)))
  # print(paste0('envir$msig_result_enricher: \n', head(envir$msig_result_enricher)))
  # print(paste0('envir$msig_result_fgsea: \n', head(envir$msig_result_fgsea)))
 
  
  output$num_of_mapped_enricher <- flexdashboard::renderValueBox({
    num_genes_mapped <- str_split(noquote(envir$msig_result_enricher@result$GeneRatio[1]), '/')[[1]][2]
    shinydashboard::box(
      title = 'Number of genes mapped',
      width = 6,
      background = 'light-blue',
      num_genes_mapped
    )
  })
  
  output$num_of_total_genes_enricher <- flexdashboard::renderValueBox({
    num_genes_total <- length(envir$msig_result_enricher@gene)
    shinydashboard::box(
      title = 'Number of genes total',
      width = 6,
      background = 'light-blue',
      num_genes_total
    )
  })
  
  
  
  
  
  
  output[["fgsea_table_PPI"]] <- renderPlot({
    req(length(input$fgsea_table_cell_clicked) > 0)
    info_list <- input$fgsea_table_cell_clicked
    
    cellVal <- strsplit(info_list[["value"]], perl = T, split = 'target=\"_blank\">')
    cellVal <- strsplit(cellVal[[1]][2], perl = T, split = "</a>")
    cellVal <- cellVal[[1]]
    print(cellVal)
    
    fgsea::plotEnrichment(envir$msig_result[['msig_geneSet']][[cellVal]], envir$msig_result[['fgsea_ranks']]) + 
      ggplot2::labs(title = cellVal)
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
}