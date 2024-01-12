#libs
library(tidyverse)
library(shiny)
library(plotly)
library(Matrix)

#main
load("demo_data.RData")

ui <- fluidPage (
  h1("samovaR: Microbiome generator v.0"),
  tags$hr(),
  mainPanel(
  fluidRow(column(5, selectInput("init", label = "Species", choices = species, selected = "Bifidobacterium bifidum")),
           column(5, sliderInput("init_level", label = "Initial level for generate", min = 0, max = 1, value = 0.02))),
  plotlyOutput("plot"))
)

server <- function (input, output) {
  
  output$plot <- renderPlotly({
    
    req(input$init)
    
    init <- input$init
    init_level <- input$init_level
    
    res <- prediction(init = init,
                      init_level = init_level, 
                      Rsq_cl = Rsq_cl, Rsq_il = Rsq_il, Rsq_pr = Rsq_pr,
                      res_sp = res_sp,
                      res_scale = res_scale,
                      res_sp_scale = res_sp_scale,
                      kmeans = kmeans,
                      trsh = trsh,
                      species = species,
                      min_res = min_res_sp_scale)
    
    res <- res[res$res_scale>trsh,]
    res$res_scale <- res$res_scale / sum(res$res_scale)
    
    fig <- res %>% 
      plot_ly(labels = ~sp, values = ~res_scale,
              textposition = 'inside',
              textinfo = 'percent',
              insidetextfont = list(color = '#FFFFFF'),
              hoverinfo = ~sp,
              marker = list(colors = colors, line = list(color = '#FFFFFF', width = 1)),
              showlegend = FALSE) %>% 
      add_pie(hole = 0.4) %>%
      layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    
    fig
    })
}

shinyApp(ui, server)
