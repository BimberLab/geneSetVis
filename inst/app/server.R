
server = function(input, output, session) {

  options(shiny.maxRequestSize=50*1024^2)

  appDiskCache <- diskCache(system.file("cache", package = 'geneSetVis'), max_size = 75*1024^2, evict = 'lru', logfile = stdout())


  envir <- reactiveValues(
    gene_list = NULL,
    stringRes = NULL,
    msigRes = NULL
  )


  #---------------------------
  observeEvent(input$demo1, {
    updateTextInput(session, 'areaInput', value = demo1)
    updateRadioButtons(session, 'inputType', selected = "Gene & avg. LogFC")
    updateRadioButtons(session, 'geneIdType', selected = "Symbol")
  })

  observeEvent(input$demo2, {
    updateTextInput(session, 'areaInput', value = demo2)
    updateRadioButtons(session, 'inputType', selected = "Gene only")
    updateRadioButtons(session, 'geneIdType', selected = "Symbol")
  })

  observeEvent(input$submit, {
    if (is.null(input$fileInput)){
      if (input$inputType == "Gene only") {
        gene_list <- read.table(text = gsub(",", "\n", perl = TRUE, x = input$areaInput),
                                header = FALSE,
                                col.names = c("gene"),
                                quote = "",
                                allowEscapes = T)
        gene_list <- data.frame(gene = gene_list, avg_logFC = NA)
      } else {
        gene_list <- read.table(text = gsub("(?<=[a-z])\\s+", "\n", perl = TRUE, x = input$areaInput),
                                header = FALSE,
                                col.names = c("gene", "avg_logFC"),
                                quote = "",
                                allowEscapes = T)
      }
    } else {
      #TODO: rm excel skip lines
      fileType <- tools::file_ext(input$fileInput)
      if (fileType == 'xlsx') {
        gene_list <- readxl::read_excel(path = input$fileInput$datapath, sheet = 1, skip = 1, col_names = T)
      }
      if (fileType == 'csv') {
        gene_list <- read.csv(file = input$fileInput$datapath, header = T, sep = ',')
      }

      gene_list <- gene_list %>% dplyr::select(gene, avg_logFC, p_val_adj) %>% dplyr::filter(p_val_adj <= 0.5)
    }

    if (input$checkGeneIdTranslate == T) {
      withProgress(message = 'Translating genes..', {
        print(paste0('gene translate: ', input$checkGeneIdTranslate))
        print(paste0('gene id type: ', input$geneIdType))
        cacheKey <- makeDiskCacheKey(list(gene_list, input$checkGeneIdTranslate, input$geneIdType), 'genelist')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')

          if (input$geneIdType == 'Symbol') {
            ensemblIds <- NULL
            geneSymbols <- gene_list$gene
          } else {
            ensemblIds <- gene_list$gene
            geneSymbols <- NULL
          }

          # gene_list_tr <- TranslateGeneNames(ensemblIds = ensemblIds, geneSymbols = geneSymbols, davidEmail = 'oosap@ohsu.edu',
          #                                     useEnsembl = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
          #                                     useSTRINGdb = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
          #                                     useDAVID = F)

          gene_list_tr <- TranslateToEnsembl(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
          gene_list_tr <- gene_list_tr[, !(colnames(gene_list_tr) %in% c('EnsemblId', 'GeneSymbol'))]

          gene_list <- dplyr::bind_cols(gene_list[1:nrow(gene_list),], gene_list_tr[1:nrow(gene_list_tr),])
          appDiskCache$set(key = cacheKey, value = gene_list)
        } else {
          gene_list <- cacheVal
        }
      })
    }

    envir$gene_list <- gene_list
  })

  output$inputTable <- DT::renderDataTable(server = FALSE, {
    validate(need(envir$gene_list, "Please enter the gene list and hit submit"))

    req(input$submit)
    envir$gene_list %>%
    DT::datatable(
      extensions = c('Buttons'),
      options = list(
        #autoWidth = TRUE, autoHeight = TRUE, scrollX = TRUE, scrollY = TRUE,
        dom = 'lBfrtip',
        lengthMenu = list(c(15, 30, 50, -1), c('15', '30', '50', 'All')),
        pageLength = 10,
        scrollX = TRUE,
        buttons = list(
          #list(extend = "collection", text = 'Show All', action = DT::JS("function ( e, dt, node, config ) { dt.page.len(-1); dt.ajax.reload(); }")),
          list(extend = 'collection', text = 'Download/Copy', buttons = c('copy', 'csv', 'excel') )
        )
      )
    )
  })

  observe({
    geneColnames <- envir$gene_list
    geneColnames['avg_logFC'] <- NULL
    updateSelectInput(session, "stringdb_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "msigdb_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "reactome_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "david_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "dose_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "ncg_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "dgn_selectGeneCol", choices = colnames(geneColnames))
    updateSelectInput(session, "enrichr_selectGeneCol", choices = colnames(geneColnames))
    })

  stringdbModule(session, input, output, envir, appDiskCache)
  msigdbModule(session, input, output, envir, appDiskCache)
  reactomeModule(session, input, output, envir, appDiskCache)
  davidModule(session, input, output, envir, appDiskCache)
  doseModule(session, input, output, envir, appDiskCache)
  ncgModule(session, input, output, envir, appDiskCache)
  dgnModule(session, input, output, envir, appDiskCache)
  enrichrModule(session, input, output, envir, appDiskCache)


  observeEvent(input$make_report, {
    if (isTruthy(input$fileInput$name)) {
      runname <- tools::file_path_sans_ext(basename(input$fileInput$name))
    } else {
      runname <- Sys.time()
    }
    stringdbRes <- envir$stringdbRes
    msigdbRes <- envir$msigdbRes
    reactomeRes <- envir$reactomeRes
    davidRes <- envir$davidRes
    doseRes <- envir$doseRes
    dgnRes <- envir$dgnRes
    ncgRes <- envir$ncgRes
    enrichrRes <- envir$enrichrRes
    rmarkdown::render(input = system.file('report/report.Rmd', package = 'geneSetVis'), output_format = 'html_clean', output_file = paste0(runname,'_Report.html'), output_dir = system.file('report/exports/', package = 'geneSetVis'))
  })
}

