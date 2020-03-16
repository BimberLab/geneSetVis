
source('fxs.R', local = TRUE)


##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 1 - loadData
##----------------------------------------------------------------------------##

tab_load_data <- tabItem(
  tabName = 'loadData',
  fluidRow(
  textAreaInput(
    inputId = 'areaInput',
    label = 'Enter Gene and Average Log Fold Change (avg. LogFC)',
    value = 'CDKN1A	0.7265868 \nHMOX1	1.0596510 \nENSMMUG00000062894	0.9929236 \nRNF167	0.9790608  \nFTH1.1	0.2733286 \nVIM	0.3409602 \nFN1	0.4090008 \nENSMMUG00000049833	0.9146529 \nFYB1	-0.4637513 \nRPL37A	-0.6478795 \nTXNIP	-0.6576339',
    width = NULL,
    height = '300px',
    cols = NULL,
    rows = NULL,
    placeholder = NULL,
    resize = NULL
  ),
  actionButton('submit', 'Submit'), 
  tags$hr(style='border-color: blue;'),
  tableOutput(('inputTable'))
  )
)

##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 2 - STRINGdb
##----------------------------------------------------------------------------##
tab_stringdb <- tabItem(
  tabName = 'stringdb',
  box(
    title = tagList(p('Run STRINGdb', style = "padding-right: 5px; display: inline"), 
                    actionButton(
                      inputId = "stringdb_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information for this panel."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput('runstringdb_select_parameters'),
      actionButton('runstringdb_button', 'Run')
    )),
  fluidRow(
    valueBoxOutput('num_of_mapped'),
    valueBoxOutput('num_of_total_genes')
  ),
  box(
    title = 'Network (PNG)',
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      imageOutput('stringdb_network_png')
    )
  ),
  tabBox(
    title = 'GO',
    side = 'right',
    height = NULL,
    selected = 'Table',
    width = 16,
    tabPanel('Table', dataTableOutput('stringdb_GO'), 
             #style = 'height:500px; overflow-y: scroll;overflow-x: scroll;', 
             collapsible = TRUE)
  ),
  tabBox(
    title = 'KEGG',   
    side = 'right',
    height = NULL,
    selected = 'Table',
    width = 16,
    tabPanel('Table', dataTableOutput('stringdb_KEGG'), 
             #style = 'height:500px; overflow-y: scroll;overflow-x: scroll;', 
             collapsible = TRUE)
  )
)


##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 3 - MSigDBr
##----------------------------------------------------------------------------##
tab_msigdbr <- tabItem(
  tabName = 'msigdbr',
  box(
    title = tagList(p('Run MSigDB', style = "padding-right: 5px; display: inline"), 
                    actionButton(
                      inputId = "msigdbr_resource_info",
                      label = "info",
                      icon = NULL,
                      class = "btn-xs",
                      title = "Show additional information for this panel."
                    )),
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput('runmsigdbr_select_parameters'),
      uiOutput('runmsigdbr_select_parameters_sub'),
      actionButton('runmsigdbr_button', 'Run')
    )
  ),
  
  makeTabBox(title = 'Enricher', key = 'enricher'),
  makeTabBox(title = 'FGSEA', key = 'fgsea')
)

