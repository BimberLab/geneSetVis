
source('fxs.R', local = TRUE)


##----------------------------------------------------------------------------##
## Tabs -- UI.R
## 1 - loadData
##----------------------------------------------------------------------------##

tab_load_data <- shinydashboard::tabItem(
  tabName = 'loadData',
  fluidRow(
  textAreaInput(
    inputId = 'areaInput',
    label = 'Enter Gene and Average Log Fold Change (avg. LogFC)',
    value = "ENSMMUG00000051300	1.0497763 \nENSMMUG00000002075	1.0150748 \nENSMMUG00000028701	0.9422637 \nENSMMUG00000053137	0.9164004 \nENSMMUG00000049855	0.5698137 \nHMOX1	1.0596510 \nRNF167	0.9790608 \nHSPA5	0.7293491 \nCDKN1A	0.7265868 \nFCGR2B	0.6369659 \nPFN1	0.5453499 \nLAPTM5	0.5164539 \nAHNAK	0.5045917 \nFN1	0.4090008 \nS100A10	0.3566574 \nVIM	0.3409602 \nYWHAZ	0.2911121 \nFTH1.1	0.2733286 \nPDIA3	0.2555106 \nATP5MPL	-0.2565952 \nLAMTOR4	-0.2574608 \nSMDT1	-0.2589715 \nCOX5A	-0.2610802 \nMTDH	-0.2619066 \nNDUFA2	-0.2638782 \nCOX6C	-0.2679750 \nCOX8A	-0.2756591 \nNDUFA1	-0.2781574 \nH2AFJ	-0.2827520 \nTOMM7	-0.2955068 \nRPL23	-0.3009606 \nCOX7C	-0.3324625 \nCASP1	-0.3531754 \nRPS21	-0.3921719 \nRPL38	-0.3928734 \nFOS	-0.8496947 \nIGFBP1	-2.2179911 \nPPT1	0.2956121 \nHEXB	0.2665466 \nNINJ1	0.3056079 \nFGL2	0.2589270 \nLDHA	0.2736736 \nCD59	-0.3042252 \nGSN	0.2728750 \nENSMMUG00000051444	0.4035713 \nANXA2	0.2990603 \nLGALS3	0.2911058 \nSLC2A3	0.4835044 \nMT-CO2	-0.3797473 \nPLIN2	0.2974303 \nPLAUR	0.2979632 \nENSMMUG00000062238	0.4565557 \nPPP1R15A	0.3040476 \nMamu-DPA1	0.3566418 \nENSMMUG00000008111	0.274841",
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
tab_stringdb <- shinydashboard::tabItem(
  tabName = 'stringdb',
  shinydashboard::box(
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
    #textOutput('string_map_stats')
    flexdashboard::valueBoxOutput('num_of_mapped'),
    flexdashboard::valueBoxOutput('num_of_total_genes')
  ),
  shinydashboard::box(
    title = 'Network (PNG)',
    status = 'primary',
    solidHeader = TRUE,
    width = 16,
    collapsible = TRUE,
    tagList(
      imageOutput('stringdb_network_png')
    )
  ),
  shinydashboard::tabBox(
    title = 'GO',
    side = 'right',
    height = NULL,
    selected = 'Table',
    width = 16,
    tabPanel('Table', dataTableOutput('stringdb_GO'), 
             #style = 'height:500px; overflow-y: scroll;overflow-x: scroll;', 
             collapsible = TRUE)
  ),
  shinydashboard::tabBox(
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
tab_msigdbr <- shinydashboard::tabItem(
  tabName = 'msigdbr',
  shinydashboard::box(
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

