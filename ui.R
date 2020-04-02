
ui = shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = 'geneSetVis'),
  shinydashboard::dashboardSidebar(
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
      ),
      shinydashboard::menuItem(
        text = 'DAVID',
        tabName = 'david',
        icon = icon(NULL)
      ),
      shinydashboard::menuItem(
        text = 'DOSE',
        tabName = 'dose',
        icon = icon(NULL)
      ),
      shinydashboard::menuItem(
        text = 'NCG',
        tabName = 'ncg',
        icon = icon(NULL)
      ),
      shinydashboard::menuItem(
        text = 'DGN',
        tabName = 'dgn',
        icon = icon(NULL)
      ),
      shinydashboard::menuItem(
        text = 'enrichR',
        tabName = 'enrichr',
        icon = icon(NULL)
      )
    )
  ),
  shinydashboard::dashboardBody(
    shinydashboard::tabItems(tab_load_data,
                             tab_stringdb,
                             tab_msigdbr,
                             tab_reactome,
                             tab_david,
                             tab_dose,
                             tab_ncg,
                             tab_dgn,
                             tab_enrichr)
  )
)
