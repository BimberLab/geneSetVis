##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 1 - loadData
##----------------------------------------------------------------------------##

tab_load_data <- tabItem(
  ## more like tabId, attached to menuItem(tabName = ...)
  tabName = "loadData",
  fluidRow(
    column(12,
           titlePanel("Load data"),
           fileInput(
             inputId = "input_file",
             label = "Select input data (.rds file)",
             multiple = FALSE,
             accept = c(".rds"),
             width = '350px',
             buttonLabel = "Browse...",
             placeholder = "No file selected"
           )
    )
  ), 
  textAreaInput(
    inputId = "data_input",
    label = "Enter Gene and/or Average Log Fold Change (avg. LogFC)",
    value = "",
    width = NULL,
    height = NULL,
    cols = NULL,
    rows = NULL,
    placeholder = NULL,
    resize = NULL
  ),
  actionButton("submit", "Submit"), 
  tags$hr(style="border-color: blue;"), 
  fluidRow(box(
    title = "Run STRINGdb",
    uiOutput("runstringdb_select_parameters"),
    actionButton("runstringdb_button", "Run")
  )
  ), 
  fluidRow(box(
    title = "Run MSigDB",
    uiOutput("runmsigdbr_select_parameters"),
    uiOutput("runmsigdbr_select_parameters_sub"),
    uiOutput("runmsigdbr_select_parameters_cont"),
    actionButton("runmsigdbr_button", "Run")))
)

##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 2 - STRINGdb
##----------------------------------------------------------------------------##
tab_stringdb <- tabItem(
  tabName = "stringdb",
  box(
    title = "Select input",
    status = "primary",
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput("stringdb_select_run_UI"),
      uiOutput("stringdb_select_cluster_UI")
    )
  ),
  fluidRow(
    valueBoxOutput("num_of_mapped"),
    valueBoxOutput("num_of_total_genes")
  ),
  box(
    title = "Network",
    status = "primary",
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      plotOutput("stringdb_network")
    )
  ),
  box(
    title = "Network (PNG)",
    status = "primary",
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      imageOutput("stringdb_network_png")
    )
  ),
  tabBox(
    title = "GO",   
    side = "right",
    height = NULL,
    selected = "Table",
    width = 16,
    tabPanel("Resource",
             textOutput("stringdb_GO_resource_info"), 
             style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
    tabPanel("Explore", 
             uiOutput("stringdb_select_GO_ann"), 
             textOutput("stringdb_select_GO_ann_output")),
    tabPanel("Table", dataTableOutput("stringdb_GO"), 
             #style = "height:500px; overflow-y: scroll;overflow-x: scroll;", 
             collapsible = TRUE
    )
  ),
  tabBox(
    title = "KEGG",   
    side = "right",
    height = NULL,
    selected = "Table",
    width = 16,
    tabPanel("Resource",
             textOutput("stringdb_KEGG_resource_info"), 
             style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
    tabPanel("Explore", 
             uiOutput("stringdb_select_KEGG_ann")),
    tabPanel("Table", dataTableOutput("stringdb_KEGG"), 
             #style = "height:500px; overflow-y: scroll;overflow-x: scroll;", 
             collapsible = TRUE
    )
  )
)

##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 3 - MSigDBr
##----------------------------------------------------------------------------##

tab_msigdbr <- tabItem(
  tabName = "msigdbr",
  box(
    title = "Select input",
    status = "primary",
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput("msigdbr_select_run_UI"),
      uiOutput("msigdbr_select_cluster_UI")
    )
  ),
  tabBox(
    title = "clusterProfiler (Enricher)",   
    side = "right",
    height = NULL,
    selected = "Table",
    width = 16,
    tabPanel("Resource",
             textOutput("msigdbr_enricher_resource_info"), 
             style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
    tabPanel("Table", dataTableOutput("msigdbr_select_cluster_enricher_table"), 
             #style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
    ),
    tabPanel("Plot", plotOutput("msigdbr_select_cluster_enricher_plot"), 
             #style = "height:500px; overflow-y: scroll;overflow-x: scroll;", 
             collapsible = TRUE),
    tabPanel("Info", 
             fluidRow(
               valueBoxOutput("num_of_mapped_enricher"),
               valueBoxOutput("num_of_total_genes_enricher")), 
             uiOutput("msigdbr_enricher_select_PA_ann"), 
             textOutput("msigdbr_enricher_select_PA_ann_output")
    ), 
    tabBox(
      title = "FGSEA",   
      side = "right",
      height = NULL,
      selected = "Table",
      width = 16,
      tabPanel("Resource",
               textOutput("msigdbr_fgsea_resource_info"), 
               #style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
      ),
      tabPanel("Table", dataTableOutput("msigdbr_select_cluster_fgsea_table"), 
               #style = "height:500px; overflow-y: scroll;overflow-x: scroll;", 
               collapsible = TRUE
      ),
      tabPanel("GTable", plotOutput("msigdbr_select_cluster_fgsea_gtable"))
    )
  ),
  tabBox(
    title = "FGSEA-clusterProfiler Plots",   
    side = "right",
    height = NULL,
    selected = "dotplot",
    width = 16,
    # tabPanel("dotplot", plotOutput("fgsea_cpplots_select_cluster_dotplot_plot")),
    # tabPanel("emapplot", plotOutput("fgsea_cpplots_select_cluster_emapplot_plot")),
    # tabPanel("cnetplot", plotOutput("fgsea_cpplots_select_cluster_cnetplot_plot")),
    # tabPanel("upsetplot", plotOutput("fgsea_cpplots_select_cluster_upsetplot_plot")),
    # tabPanel("heatplot", plotOutput("fgsea_cpplots_select_cluster_heatplot_plot")), 
    tabPanel("dotplot", plotlyOutput("fgsea_dotplot")),
    tabPanel("emapplot", plotlyOutput("fgsea_emapplot")),
    tabPanel("cnetplot", plotlyOutput("fgsea_cnetplot")),
    tabPanel("upsetplot", plotlyOutput("fgsea_upsetplot")),
    tabPanel("heatplot", plotlyOutput("fgsea_heatplot"))
  ),
)

tab_clusterprofiler <- tabItem(
  tabName = "clusterprofiler",
  box(
    title = "Select Cluster",
    status = "primary",
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput("clusterprofiler_select_run_UI"),
      uiOutput("clusterprofiler_select_cluster_UI")
    )
  )
)

