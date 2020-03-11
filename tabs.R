
source('fxs.R', local = TRUE)


##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 1 - loadData
##----------------------------------------------------------------------------##

tab_load_data <- tabItem(
  tabName = 'loadData',
  fluidRow(
    column(12,
           titlePanel('Load data'),
           fileInput(
             inputId = 'input_file',
             label = 'Select input data (.rds file)',
             multiple = FALSE,
             accept = c('.rds'),
             width = '350px',
             buttonLabel = 'Browse...',
             placeholder = 'No file selected'
           )
    )
  ), 
  textInput(
    inputId = 'areaInput_runname',
    label = 'Run name:',
    placeholder = 'Run1'
  ),
  textAreaInput(
    inputId = 'areaInput',
    label = 'Enter Gene and/or Average Log Fold Change (avg. LogFC)',
    value = '',
    width = NULL,
    height = NULL,
    cols = 3,
    rows = NULL,
    placeholder = "HMOX 2.00 \nABCA4 -1.50",
    resize = NULL
  ),
  actionButton('submit', 'Submit'), 
  tags$hr(style='border-color: blue;'),
  tableOutput(('inputTable'))
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
  box(
    title = 'Select input',
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput('stringdb_select_run_UI'),
      uiOutput('stringdb_select_cluster_UI')
    )
  ),
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
    tabPanel('Explore', 
             uiOutput('stringdb_select_GO_ann'), 
             textOutput('stringdb_select_GO_ann_output')),
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
    tabPanel('Explore', 
             uiOutput('stringdb_select_KEGG_ann')),
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
  box(
    title = 'Select input',
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput('msigdbr_select_run_UI'),
      uiOutput('msigdbr_select_cluster_UI')
    )
  ),
  makeTabBox(title = 'Enricher', key = 'enricher'),
  makeTabBox(title = 'FGSEA', key = 'fgsea')
)

##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 4 - clusterProfiler
##----------------------------------------------------------------------------##
tab_clusterprofiler <- tabItem(
  tabName = 'clusterprofiler',
  box(
    title = tagList(p('Select Enricher-type', style = "padding-right: 5px; display: inline"), 
                    actionButton(
                      inputId = "clusterprofiler_resource_info",
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
      
    )
  ),
  box(
    title = 'Select Cluster',
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      uiOutput('clusterprofiler_select_run_UI'),
      uiOutput('clusterprofiler_select_cluster_UI')
    )
  )
)


