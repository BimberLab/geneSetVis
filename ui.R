
source('fxs.R', local = TRUE)
source('tabs.R', local = TRUE)


ui = dashboardPage(
  dashboardHeader(
    title = span('geneSetVis', style = 'color: white; font-size: 28px; font-weight: bold')
  ),
  dashboardSidebar(tags$head(tags$style(
    HTML('.content-wrapper {overflow-x: scroll;}')
  )), 
  sidebarMenu(
    id = 'sidebar',
    menuItem(
      text = 'Load data',
      tabName = 'loadData',
      icon = icon(NULL),
      selected = TRUE
    ),
    menuItem(
      text = 'STRINGdb',
      tabName = 'stringdb',
      icon = icon(NULL),
      selected = FALSE
    ),
    menuItem(
      text = 'MsigDB',
      tabName = 'msigdbr',
      icon = icon(NULL)
    )
  )
  ),
  dashboardBody(
    tags$script(HTML('$("body").addClass("fixed");')),
    tabItems(tab_load_data,
             tab_stringdb,
             tab_msigdbr)
  )
)
