
tab_info <- shinydashboard::tabItem(
  tabName = 'appInfo',
  fluidRow(
    shinydashboard::box(
      title = tagList(p('App Information', style = "padding-right: 5px; display: inline")),
      side = 'right',
      status = 'primary',
      solidHeader = TRUE,
      height = NULL,
      width = 16,
      uiOutput('app_info')
    )
  )
)

tab_load_data <- shinydashboard::tabItem(
  tabName = 'loadData',
  fluidRow(
    #column = 6,
    shinydashboard::box(
      title = 'Input Gene List',
      status = 'primary',
      solidHeader = TRUE,
      height = '800px',
      textAreaInput(
        inputId = 'areaInput',
        label = tagList(
          'Enter Gene-sets:',
          div(style = "display:inline-block", em('Example')),
          div(style = "display:inline-block", actionLink('demo1', '#1')),
          div(style = "display:inline-block", actionLink('demo2', '#2'))
        ),
        height = '150px',
        width = NULL,
        value = NULL,
        cols = NULL,
        rows = NULL,
        placeholder = NULL,
        resize = NULL
      ),
      fileInput(
        'fileInput',
        label = tagList(
          'or Upload Gene-sets:',
          div(style = "font-weight:10px; font-size:12px;", em("Enter column headers in file. If p_val column present, will filter out >0.5")),
          div(style = "display:inline-block; font-weight:normal; font-size:12px;", textInput("file_geneCol", label = NULL, value = "gene")),
          div(style = "display:inline-block; font-weight:normal; font-size:12px;", textInput("file_avgLogFCcol", label = NULL, value = "avg_logFC")),
          div(style = "display:inline-block; font-weight:normal; font-size:12px;", textInput("file_pvaladjCol", label = NULL, value = "p_val_adj")),
        ),
        multiple = F,
        accept = c(".xlsx", ".xls", "csv")
      ),
      radioButtons("inputType",
                  label = "Input Type:",
                  choices = list("Gene only", "Gene & avg. LogFC"),
                  selected = 'Gene & avg. LogFC'),
      radioButtons("geneIdType",
                   label = "Select Gene ID Type:",
                   choices = list("Ensembl" = "Ensembl",
                                  "Symbol" = "Symbol"),
                   selected = 'Symbol'),
      withBusyIndicatorUI(actionButton('submit', 'Submit')),
    ),
    shinydashboard::box(
      title = 'Parsed Genes',
      status = 'primary',
      solidHeader = TRUE,
      height = '800px',
      DT::dataTableOutput('inputTable')
    )
  )
)

if (exists('gsvis_package')) {
  stringdbSpecies <- read.csv(file = system.file('app/intdata/stringdb_species.v10.txt', package = 'geneSetVis'), header = T, stringsAsFactors = F, check.names = F, quote = "", sep = "\t")
} else {
  stringdbSpecies <- read.csv(file = 'intdata/stringdb_species.v10.txt', header = T, stringsAsFactors = F, check.names = F, quote = "", sep = "\t")
}
tab_stringdb <- shinydashboard::tabItem(
  tabName = 'stringdb',
  shinydashboard::box(
    title = tagList(p('Run STRINGdb', style = "padding-right: 5px; display: inline"),
                    actionButton(
                      inputId = "stringdb_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information.."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    flowLayout(
      selectInput(
        inputId = 'stringdb_selectGeneCol',
        label = 'Select gene column to use:',
        selected = 'gene',
        choices = ''
      ),
      selectInput(
        inputId = 'stringdb_refSpecies_input',
        label = 'Reference species:',
        selected = 'Homo sapiens',
        choices = stringdbSpecies$official_name_NCBI
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
    ),
    withBusyIndicatorUI(actionButton('runstringdb_button', 'Run'))
  ),
  #browser(),
  shinydashboard::tabBox(
    title = NULL,
    side = 'right',
    height = NULL,
    selected = 'Mapped',
    width = 16,
    tabPanel('Mapped', uiOutput('stringdb_map_stats'))
  ),
  # shinydashboard::tabBox(
  #   title = NULL,
  #   side = 'right',
  #   height = NULL,
  #   selected = 'Network (PNG)',
  #   width = 16,
  #   tabPanel( 'Network (PNG)', imageOutput('stringdb_network_png'))
  # ),
  shinydashboard::tabBox(
    title = NULL,
    side = 'right',
    height = NULL,
    selected = 'GO',
    width = 16,
    tabPanel('GO', DT::dataTableOutput('stringdb_GO')),
    tabPanel('KEGG', DT::dataTableOutput('stringdb_KEGG'))
  )
)

if (exists('gsvis_package')) {
  msigdbCategories <- read.table(file = system.file('app/intdata/msigdb_categories.txt', package = 'geneSetVis'), sep = '\t', header = T, stringsAsFactors = FALSE)
} else {
  msigdbCategories <- read.table(file = 'intdata/msigdb_categories.txt', sep = '\t', header = T, stringsAsFactors = FALSE)
}
msigdbCategories <- unique(msigdbCategories[c('Category', 'CategoryLabel')])
msigdbCategories$CategoryLabel <- paste0(msigdbCategories$Category, ': ', msigdbCategories$CategoryLabel)
msigdb_categories <- msigdbCategories$Category
names(msigdb_categories) <- msigdbCategories$CategoryLabel
rm(msigdbCategories)

tab_msigdb <- shinydashboard::tabItem(
  tabName = 'msigdb',
  shinydashboard::box(
    title = tagList(p('Run MSigDB', style = "padding-right: 5px; display: inline"),
                    actionButton(
                      inputId = "msigdb_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information.."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    HTML('<a href="https://www.gsea-msigdb.org/gsea/msigdb/index.jsp" target="_blank" style="font-weight: bold;">Click here for more information on available MSigDB collections</a><br><br>'),
  	flowLayout(
  	  selectInput(
  	    inputId = 'msigdb_selectGeneCol',
  	    label = 'Select gene column to use:',
  	    selected = 'gene',
  	    choices = ''
  	  ),
      selectInput(
        inputId = 'msigdb_species_input',
        label = 'Reference species',
        selected = 'Homo sapiens',
        choices = sort(msigdbr::msigdbr_species()$species_name)
      ),
      selectInput(
        inputId = 'msigdb_category_input',
        label = 'Select category (optional):',
        selected = '',
        choices = c('', msigdb_categories)
      ),
      selectInput(
        inputId = 'msigdb_subcategory_input',
        label = 'Select subcategory (optional):',
        selected = '',
        choices = NULL
      )
    ),
    withBusyIndicatorUI(actionButton('runmsigdb_button', 'Run'))
  ),
  shinydashboard::tabBox(
    title = NULL,
    side = 'right',
    height = NULL,
    selected = 'Mapped',
    width = 16,
    tabPanel('Mapped', uiOutput('msigdb_map_stats'))
  ),
  makeTabBox(title = 'Enricher', key = 'enricher'),
  makeTabBox(title = 'FGSEA', key = 'fgsea')
)


tab_reactome <- shinydashboard::tabItem(
  tabName = 'reactome',
  shinydashboard::box(
    title = tagList(p('Run Reactome', style = "padding-right: 5px; display: inline"),
                    actionButton(
                      inputId = "reactome_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information.."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    flowLayout(selectInput(
      inputId = 'reactome_selectGeneCol',
      label = 'Select gene column to use:',
      selected = 'gene',
      choices = ''
    ),
      selectInput(
        inputId = 'reactome_OrgDB_input',
        label = 'OrgDB:',
        selected = 'org.Hs.eg.db',
        #choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Rn.eg.db', 'org.Mm.eg.db')
        choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Mm.eg.db')
      )
    ),
    withBusyIndicatorUI(actionButton('runreactome_button', 'Run')),
  ),
  shinydashboard::tabBox(
    title = NULL,
    side = 'right',
    height = NULL,
    selected = 'Mapped',
    width = 16,
    tabPanel('Mapped', uiOutput('reactome_map_stats'))
  ),
  makeTabBox(title = 'Reactome', key = 'reactome')
)


tab_dose <- shinydashboard::tabItem(
  tabName = 'dose',
  shinydashboard::box(
    title = tagList(p('Run DOSE', style = "padding-right: 5px; display: inline"),
                    actionButton(
                      inputId = "dose_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information.."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    flowLayout(
      selectInput(
        inputId = 'dose_selectGeneCol',
        label = 'Select gene column to use:',
        selected = 'gene',
        choices = ''
      ),
      selectInput(
        inputId = 'dose_OrgDB_input',
        label = 'OrgDB:',
        selected = 'org.Hs.eg.db',
        #choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Rn.eg.db', 'org.Mm.eg.db')
        choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Mm.eg.db')
      )
    ),
    withBusyIndicatorUI(actionButton('rundose_button', 'Run')),
  ),
  shinydashboard::tabBox(
    title = NULL,
    side = 'right',
    height = NULL,
    selected = 'Mapped',
    width = 16,
    tabPanel('Mapped', uiOutput('dose_map_stats'))
  ),
  makeTabBox(title = 'DOSE', key = 'dose')
  #sliderInput('pvalueCutoff', label = 'p Value cutoff:', min = 0, max = 1, value = 0.5)
)


tab_ncg <- shinydashboard::tabItem(
  tabName = 'ncg',
  shinydashboard::box(
    title = tagList(p('Run NCG', style = "padding-right: 5px; display: inline"),
                    actionButton(
                      inputId = "ncg_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information.."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    flowLayout(
      selectInput(
        inputId = 'ncg_selectGeneCol',
        label = 'Select gene column to use:',
        selected = 'gene',
        choices = ''
      ),
      selectInput(
        inputId = 'ncg_OrgDB_input',
        label = 'OrgDB:',
        selected = 'org.Hs.eg.db',
        #choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Rn.eg.db', 'org.Mm.eg.db')
        choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Mm.eg.db')
      )
    ),
    withBusyIndicatorUI(actionButton('runncg_button', 'Run')),
  ),
  shinydashboard::tabBox(
    title = NULL,
    side = 'right',
    height = NULL,
    selected = 'Mapped',
    width = 16,
    tabPanel('Mapped', uiOutput('ncg_map_stats'))
  ),
  makeTabBox(title = 'NCG', key = 'ncg')
)


tab_dgn <- shinydashboard::tabItem(
  tabName = 'dgn',
  shinydashboard::box(
    title = tagList(p('Run DisGeNET', style = "padding-right: 5px; display: inline"),
                    actionButton(
                      inputId = "dgn_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information.."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    flowLayout(
      selectInput(
        inputId = 'dgn_selectGeneCol',
        label = 'Select gene column to use:',
        selected = 'gene',
        choices = ''
      ),
      selectInput(
        inputId = 'dgn_OrgDB_input',
        label = 'OrgDB:',
        selected = 'org.Hs.eg.db',
        #choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Rn.eg.db', 'org.Mm.eg.db')
        choices = c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Mm.eg.db')
      )
    ),
    withBusyIndicatorUI(actionButton('rundgn_button', 'Run')),
  ),
  shinydashboard::tabBox(
    title = NULL,
    side = 'right',
    height = NULL,
    selected = 'Mapped',
    width = 16,
    tabPanel('Mapped', uiOutput('dgn_map_stats'))
  ),
  makeTabBox(title = 'DisGeNET', key = 'dgn')
)


tab_enrichr <- shinydashboard::tabItem(
  tabName = 'enrichr',
  shinydashboard::box(
    title = tagList(p('Run enrichR', style = "padding-right: 5px; display: inline"),
                    actionButton(
                      inputId = "enrichr_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information.."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    HTML('<a href="https://amp.pharm.mssm.edu/Enrichr/#stats" target="_blank" style="font-weight: bold;">Click here for more information on available enrichR databases/libraries</a><br><br>'),
    flowLayout(
      selectInput(
        inputId = 'enrichr_selectGeneCol',
        label = 'Select gene column to use:',
        selected = 'gene',
        choices = ''
      ),
      shinyWidgets::pickerInput(inputId = "enrichr_db",
                                label = 'Select database(s) to query:',
                                choices = c('', enrichR::listEnrichrDbs()$libraryName),
                                options = list(`actions-box` = TRUE),
                                multiple = T)
      ),
    withBusyIndicatorUI(actionButton('runenrichr_button', 'Run'))
  ),
    selectInput(
      inputId = 'enrichr_selectQuery',
      label = 'Select query result to view:',
      choices = ''
  ),
  #uiOutput('enrichrResults_selected_ui'),
  #dataTableOutput('enrichrResults_selected_table'),
  # ),
  # shinydashboard::tabBox(
  #   title = NULL,
  #   side = 'right',
  #   height = NULL,
  #   selected = 'Mapped',
  #   width = 16,
  #   tabPanel('Mapped', uiOutput('enrichr_map_stats'))
  # ),
  makeTabBox(title = 'enrichR', key = 'enrichr')
)

