
source('fxs.R', local = TRUE)
source('tabs.R', local = TRUE)


ui = shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(
    title = span('geneSetVis', style = 'color: white; font-size: 28px; font-weight: bold')
  ),
  shinydashboard::dashboardSidebar(tags$head(tags$style(
    HTML('.content-wrapper {overflow-x: scroll;}')
  )), 
  shinydashboard::sidebarMenu(
    id = 'sidebar',
    shinydashboard::menuItem(
      text = 'Load data',
      tabName = 'loadData',
      icon = icon(NULL),
      selected = TRUE
    ),
    shinydashboard::menuItem(
      text = 'STRINGdb',
      tabName = 'stringdb',
      icon = icon(NULL),
      selected = FALSE
    ),
    shinydashboard::menuItem(
      text = 'MsigDB',
      tabName = 'msigdbr',
      icon = icon(NULL)
    ),
    shinydashboard::menuItem(
      text = 'Reactome',
      tabName = 'reactome',
      icon = icon(NULL)
    )
  )
  ),
  shinydashboard::dashboardBody(
    tags$script(HTML('$("body").addClass("fixed");')),
    shinydashboard::tabItems(tab_load_data,
             tab_stringdb,
             tab_msigdbr,
             tab_reactome)
  )
)
