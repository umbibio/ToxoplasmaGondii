
ui <- fluidPage(
  
  title = 'Pseudo Time Expression Profiler',
  
  
  h1('Expression Curve'),
  
  fluidRow(
    column(4, plotOutput('x1', height = 300)),
    column(4, plotOutput('x2', height = 300)),
    column(4, DT::dataTableOutput('x3'))
    
    #column(4, DT::dataTableOutput('x3'), hight = 300)
    #column(4, plotOutput('x3', height = 300))
  ),
  
  hr(),
  
  fluidRow(
    column(12, DT::dataTableOutput('x4'))
  )
  
)