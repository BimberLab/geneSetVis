
source('fxs.R', local = TRUE)
source('tabs.R', local = TRUE)


ui = dashboardPage(
  dashboardHeader(
    title = span('GeneSetVis', style = 'color: white; font-size: 28px; font-weight: bold')
  ),
  dashboardSidebar(tags$head(tags$style(
    HTML('.content-wrapper {overflow-x: scroll;}')
  )),
  sidebarMenu(sidebarMenuOutput('sidebar_menu'))),
  dashboardBody(
    tags$script(HTML('$("body").addClass("fixed");')),
    tabItems(tab_load_data,
             tab_stringdb,
             tab_msigdbr, 
             tab_clusterprofiler)
  )
)

