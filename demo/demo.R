# SPECIFY AND RUN BEFORE EXCECUTING THE SCRIPT ----
wd <- "PATH/TO/DEMO"
#load(paste0(sd, "/demo_data.RData"))
#####

# MAIN ----

#libs
library(tidyverse)
library(shiny)
library(plotly)

ui <- fluidPage (
  h1("samovaR: Microbiome generator v.0.2"),
  tags$hr(),
  mainPanel(
    fluidRow(column(4, actionButton("generate", label = "Generate data", width = 400)),
             column(4, actionButton("add", label = "Add data", width = 400)),
             column(2, downloadButton("save", label = "Save data",))),
    tags$br(),
    fluidRow(column(5, selectInput("init", label = "Species", choices = species, selected = "Bifidobacterium bifidum", width = 500)),
             column(5, sliderInput("init_level", label = "Initial level for generate", min = 0, max = 1, value = 0.02, width = 500))),
    tags$br(),
    h2("Generated data"),
    plotlyOutput("plot", width = 1000),
    tags$hr(),
    h2("Preview table before download"),
    tableOutput("preview"),
    textOutput("specifics")
  )
)

server <- function (input, output) {
  #load("~/bioinformatics/R/samovaR/not_to_pull/results/demo/demo_data.RData")
  results <- as.data.frame(results[,1])
  colnames(results) <- "sp"
  counter <- 0
  
  output$plot <- renderPlotly({
    req(input$generate)
    
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
  
  output$preview <- renderTable({
    req(input$add)
    counter <<- counter + 1
    colnames(res)[2] <- paste(counter, input$init, input$init_level, sep = "_")
    results <<- full_join(results, res, by = "sp")
    head(results)
  })
  
  output$save <- downloadHandler(filename = "generated.csv",
                                 content = function(file) {
                                   write.csv(results, file, row.names = F)
                                   })
  
  output$specifics <- renderText({
    req(input$add)
    if(ncol(results) - 1 == 1){
      paste0("Data frame of ", nrow(results), " species, generated ", ncol(results) - 1, " sample")
    } else {
      paste0("Data frame of ", nrow(results), " species, generated ", ncol(results) - 1, " sample")
    }
  })
}

shinyApp(ui, server)
