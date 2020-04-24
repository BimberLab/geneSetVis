
server = function(input, output, session) {

  options(shiny.maxRequestSize=50*1024^2)

  envir <- reactiveValues(
    cachedir = NULL,
    geneList = NULL,
    namedGeneList = NULL,
    stringRes = NULL,
    msigRes = NULL
  )

  if (Sys.getenv("GENESETVIS_CACHE") == "") {
    cachedir <- paste0(Sys.getenv("TMPDIR"), "geneSetVis-cache")
  } else {
    cachedir <-  Sys.getenv("GENESETVIS_CACHE")
  }
  print(paste0('cache directory in ', cachedir))
  envir$cachedir <- cachedir
  appDiskCache <- diskCache(cachedir, max_size = 75*1024^2, evict = 'lru', logfile = stdout())


  demo1 <- "HMOX1	1.0596510 \nRNF167	0.9790608 \nHSPA5	0.7293491 \nCDKN1A	0.7265868 \nFCGR2B	0.6369659 \nPFN1	0.5453499 \nLAPTM5	0.5164539 \nAHNAK	0.5045917 \nFN1	0.4090008 \nS100A10	0.3566574 \nVIM	0.3409602 \nYWHAZ	0.2911121 \nFTH1.1	0.2733286 \nPDIA3	0.2555106 \nATP5MPL	-0.2565952 \nLAMTOR4	-0.2574608 \nSMDT1	-0.2589715 \nCOX5A	-0.2610802 \nMTDH	-0.2619066 \nNDUFA2	-0.2638782 \nCOX6C	-0.2679750 \nCOX8A	-0.2756591 \nNDUFA1	-0.2781574 \nH2AFJ	-0.2827520 \nTOMM7	-0.2955068 \nRPL23	-0.3009606 \nCOX7C	-0.3324625 \nCASP1	-0.3531754 \nRPS21	-0.3921719 \nRPL38	-0.3928734 \nFOS	-0.8496947 \nIGFBP1	-2.2179911 \nPPT1	0.2956121 \nHEXB	0.2665466 \nNINJ1	0.3056079 \nFGL2	0.2589270 \nLDHA	0.2736736 \nCD59	-0.3042252 \nGSN	0.2728750 \nANXA2	0.2990603 \nLGALS3	0.2911058 \nSLC2A3	0.4835044 \nMT-CO2	-0.3797473 \nPLIN2	0.2974303 \nPLAUR	0.2979632 \nPPP1R15A	0.3040476"
  demo2 <- "HMOX1, RNF167, HSPA5, CDKN1A, FCGR2B, PFN1, LAPTM5, AHNAK, FN1, S100A10, VIM, YWHAZ, FTH1.1, PDIA3, ATP5MPL, LAMTOR4, SMDT1, COX5A, MTDH, NDUFA2, COX6C, COX8A, NDUFA1, H2AFJ, TOMM7, RPL23, COX7C, CASP1, RPS21, RPL38, FOS, IGFBP1, PPT1, HEXB, NINJ1, FGL2, LDHA, CD59, GSN, ANXA2, LGALS3	, SLC2A3, MT-CO2, PLIN2, PLAUR, PPP1R15A"

  if (exists('gsvis_package')) {
    source(system.file('app/GeneAliasing.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/helpers.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/uiElements.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/stringdb.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/msigdb.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/reactome.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/david.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/dose.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/ncg.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/dgn.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
    source(system.file('app/modules/enrichr.R', package = 'geneSetVis', mustWork = TRUE), local = TRUE)
  } else {
    source('GeneAliasing.R', local = TRUE)
    source('helpers.R', local = TRUE)
    source('uiElements.R', local = TRUE)
    source('modules/stringdb.R', local = TRUE)
    source('modules/msigdb.R', local = TRUE)
    source('modules/reactome.R', local = TRUE)
    source('modules/david.R', local = TRUE)
    source('modules/dose.R', local = TRUE)
    source('modules/ncg.R', local = TRUE)
    source('modules/dgn.R', local = TRUE)
    source('modules/enrichr.R', local = TRUE)
  }



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
        geneList <- read.table(text = gsub(",", "\n", perl = TRUE, x = input$areaInput),
                                header = FALSE,
                                col.names = c("gene"),
                                quote = "",
                                allowEscapes = T)
        geneList <- data.frame(gene = geneList, avg_logFC = NA)
      } else {
        geneList <- read.table(text = gsub("(?<=[a-z])\\s+", "\n", perl = TRUE, x = input$areaInput),
                                header = FALSE,
                                col.names = c("gene", "avg_logFC"),
                                quote = "",
                                allowEscapes = T)
      }
    } else {
      #TODO: rm excel skip lines
      fileType <- tools::file_ext(input$fileInput)
      if (fileType == 'xlsx') {
        geneList <- readxl::read_excel(path = input$fileInput$datapath, sheet = 1, skip = 1, col_names = T)
      }
      if (fileType == 'csv') {
        geneList <- read.csv(file = input$fileInput$datapath, header = T, sep = ',')
      }

      geneList <- geneList %>% dplyr::select(gene, avg_logFC, p_val_adj) %>% dplyr::filter(p_val_adj <= 0.5)
    }

    if (input$checkGeneIdTranslate == T) {
      withProgress(message = 'Translating genes..', {
        print(paste0('gene translate: ', input$checkGeneIdTranslate))
        print(paste0('gene id type: ', input$geneIdType))
        cacheKey <- makeDiskCacheKey(list(geneList, input$checkGeneIdTranslate, input$geneIdType), 'genelist')
        cacheVal <- appDiskCache$get(cacheKey)
        if (class(cacheVal) == 'key_missing') {
          print('missing cache key...')

          if (input$geneIdType == 'Symbol') {
            ensemblIds <- NULL
            geneSymbols <- geneList$gene
          } else {
            ensemblIds <- geneList$gene
            geneSymbols <- NULL
          }

          # geneList_tr <- TranslateGeneNames(ensemblIds = ensemblIds, geneSymbols = geneSymbols, davidEmail = 'oosap@ohsu.edu',
          #                                     useEnsembl = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
          #                                     useSTRINGdb = ifelse('useSTRINGdb' %in% input$select_gene_conversion, T, F),
          #                                     useDAVID = F)

          geneList_tr <- TranslateToEnsembl(ensemblIds = ensemblIds, geneSymbols = geneSymbols)
          geneList_tr <- geneList_tr[, !(colnames(geneList_tr) %in% c('EnsemblId', 'GeneSymbol'))]

          geneList <- dplyr::bind_cols(geneList[1:nrow(geneList),], geneList_tr[1:nrow(geneList_tr),])
          appDiskCache$set(key = cacheKey, value = geneList)
        } else {
          geneList <- cacheVal
        }
      })
    }

    envir$geneList <- geneList
    envir$namedGeneList <- setNames(geneList$avg_logFC, geneList$gene)
  })

  output$inputTable <- DT::renderDataTable({
    validate(need(envir$geneList, "Please enter the gene list and hit submit"))

    req(input$submit)
    envir$geneList %>%
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
    geneColnames <- envir$geneList
    geneColnames['avg_logFC'] <- NULL
    updateSelectInput(session, "stringdb_selectGeneCol", choices = colnames(geneColnames)[1])
    updateSelectInput(session, "msigdb_selectGeneCol", choices = colnames(geneColnames)[1])
    updateSelectInput(session, "reactome_selectGeneCol", choices = colnames(geneColnames)[1])
    updateSelectInput(session, "david_selectGeneCol", choices = colnames(geneColnames)[1])
    updateSelectInput(session, "dose_selectGeneCol", choices = colnames(geneColnames)[1])
    updateSelectInput(session, "ncg_selectGeneCol", choices = colnames(geneColnames)[1])
    updateSelectInput(session, "dgn_selectGeneCol", choices = colnames(geneColnames)[1])
    updateSelectInput(session, "enrichr_selectGeneCol", choices = colnames(geneColnames)[1])
    })

  stringdbModule(session, input, output, envir, appDiskCache)
  msigdbModule(session, input, output, envir, appDiskCache)
  reactomeModule(session, input, output, envir, appDiskCache)
  davidModule(session, input, output, envir, appDiskCache)
  doseModule(session, input, output, envir, appDiskCache)
  ncgModule(session, input, output, envir, appDiskCache)
  dgnModule(session, input, output, envir, appDiskCache)
  enrichrModule(session, input, output, envir, appDiskCache)


  runname <- reactive({
    if (isTruthy(input$fileInput$name)) {
      tools::file_path_sans_ext(basename(input$fileInput$name))
    } else {
      paste0("geneSetVis", ".html")
      #paste0("geneSetVis-",gsub(":","-",Sys.Date()), ".html")
    }
  })

  output$downloadReport <- downloadHandler(
    filename = function() {
      runname()
    },
    content = function(file) {
      runname <- tools::file_path_sans_ext(runname())
      stringdbRes <- envir$stringdbRes
      msigdbRes <- envir$msigdbRes
      reactomeRes <- envir$reactomeRes
      davidRes <- envir$davidRes
      doseRes <- envir$doseRes
      dgnRes <- envir$dgnRes
      ncgRes <- envir$ncgRes
      enrichrRes <- envir$enrichrRes
      namedGeneList <- envir$namedGeneList


      if (exists('gsvis_package')) {
        file.copy(system.file('app/intdata/template_report.Rmd', package = 'geneSetVis'), paste0(envir$cachedir, "/geneSetVis-exports", "/template_report.Rmd"))
      } else {
        file.copy('intdata/template_report.Rmd', paste0(envir$cachedir, "/geneSetVis-exports", "/template_report.Rmd"))
      }
      output_path <-rmarkdown::render(
        input = paste0(envir$cachedir, "/geneSetVis-exports", "/template_report.Rmd"),
        output_format = 'html_clean',
        intermediates_dir = paste0(envir$cachedir, "/geneSetVis-exports"),
        output_dir = paste0(envir$cachedir, "/geneSetVis-exports")
        #output_file = paste0(runname, '_Report.html')
      )
      file.copy(output_path, file)
    }
  )

}

