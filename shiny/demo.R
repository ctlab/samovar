#prepare
library(tidyverse)
library(shiny)
library(plotly)

#devtools::install_github("https://github.com/ctlab/samovar")
library(samovaR)

# UI ----
ui <- fluidPage (
  # head
  h1("samovaR: Microbiome generator v.0.9.0"),
  tags$hr(),
  img(src = "logo_compressed.png", align = "right"),
  textOutput("description"),
  tags$hr(),

  # main panels
  uiOutput("main_panel"),
  tags$br(),
  tags$hr(),
  tags$br(),

  # generation panels
  uiOutput("generation_panel")
)


# UX ----

server <- function (input, output, session) {

  # misc ----
  selectInput_glm <- function(...) {
    selectInput(choices = c(
      "binomial", "gaussian", "Gamma",
      "inverse.gaussian", "poisson","quasi",
      "quasibinomial", "quasipoisson"
    ), selected = "gaussian", ...
    )
  }

  # imports ----
  load("demo_data/demo.RData")

  if(!exists("samovar_bases")) samovar_bases <<- list()
  if(!exists("teatree_bases")) teatree_bases <<- list()
  if(!exists("iter")) iter <<- 1

  # UX: stable ----
  output$description <-
    renderText(
      "For generation, select properties of your database and than select it, or load one of pre-built samova.R databases"
    )

  output$probability_description <-
    renderUI(
      tags$ul(
        tags$li('If "simple": P(A|B) = sum(A&B)/sum(A|B).'),
        tags$li('If "oriented": P(A|B) = P(A&B|B)'),
        tags$li('If "compositional": P(A|B) = P(A&B|B), and than sampled one of represented occurence conditions')
      )
    )

  # UX: render UI for database selection ----
  output$main_panel <- renderUI({

    ## import ----
    ### meshID
    meshID_table <<-
      read.csv("demo_data/meshID_phenotypes.csv", row.names = 1)
    meshID_choices <<-
      paste0(meshID_table$term, " (", meshID_table$disease, ")")


    ## out ----
    fluidPage(
      # left
      column(3,
             h3("Build samova.R database"),
             tags$br(),
             textInput("name_db", label = NULL, placeholder = "Name of DB"),
             actionButton("generate_db", label = "Generate DB"),
             tags$br(),
             tags$br(),
             tags$br(),
             tabsetPanel(
               tabPanel(
                 "Import",
                 tags$br(),
                 selectInput("meshID",
                             label = "Data group from GMRepo",
                             choices = meshID_choices,
                             selected = "Health (D006262)"),
                 selectInput(
                   "at_level",
                   label = "Summarize at level",
                   choices = c("species", "genus")
                 ),
                 numericInput("number_to_process",
                              label = "Number of runs",
                              value = 1500,
                              min = 0, max = 10000, step = 100)
               ),
               tabPanel(
                 "Treshholds",
                 tags$br(),
                 sliderInput("treshhold_species", label = "Minimum number of species per sample", min = 0, max = 100, value = 3, step = 1),
                 sliderInput("treshhold_samples", label = "Minimum number of samples per species", min = 0, max = 100, value = 3, step = 1),
                 sliderInput("treshhold_amount", label = "Minimum amount above 0", min = 1, max = 9, value = 3, step = 0.1, pre = "10^(-", post = ")")
               ),
               tabPanel(
                 "Cluster sizes",
                 tags$br(),
                 sliderInput("min_cluster_size", label = "Minimal cluster size", min = 4, max = 80, value = 40, step = 1),
                 sliderInput("max_cluster_size", label = "Minimal cluster size", min = 6, max = 100, value = 60, step = 1)
               ),
               tabPanel(
                 "Connections",
                 tags$br(),
                 selectInput("probability_calculation",
                             label = "Probability calculation",
                             choices =  c("oriented", "simple", "compositional"),
                             selected = "simple"),
                 uiOutput("probability_description"),
                 selectInput_glm("inner_model", label = "GLM method for species connections"),
                 selectInput_glm("inter_model", label = "GLM method for cluster connections")
               )
             )
             ),

      column(1),
      # right
      column(2,
             h3("Use samova.R database"),
             tags$br(),
             uiOutput("db_built"),
             checkboxInput("init_show_top", label = "Show only top species", value = T),
             sliderInput("init_top_count", width = "100%",
                         label = "Show top species",
                         min = 1, max = 50, value = 10),
             fluidRow(selectInput("init_type", label = NULL,
                                  choices = c("tile", "column"), selected = "tile")),
             tags$br(),
             uiOutput("download_button_init"),
             textOutput("error_text")),
      column(4,
             plotlyOutput("composition_init"))
    )
  })

  # create a database ----
  observeEvent(input$generate_db,
               withProgress(message = "Creating DB", value = 0, {
                 name <- paste0(iter, ". ", input$name_db)

    N <- 7
    incProgress(1/N, detail = "Download data")

    samovar_catch <- try({
    mesh_ids <- meshID_table$disease[which(input$meshID == meshID_choices)]

    teatree <- GMrepo_type2data(
      mesh_ids = mesh_ids,
      number_to_process = input$number_to_process,
      at_level = input$at_level
    )

    incProgress(1/N, detail = "Filter")
    tealeaves <- teatree %>%
      teatree_trim(treshhold_species = input$treshhold_species,
                   treshhold_samples = input$treshhold_samples,
                   treshhold_amount = 10^(-input$treshhold_amount))

    rownames(tealeaves$data) <- tealeaves$species <- teatree$species[as.numeric(tealeaves$species)]
    teatree_bases[[name]] <<- tealeaves$data %>% data.frame

    incProgress(1/N, detail = "Rescale")
    teabag <- tealeaves %>%
      tealeaves_pack()

    incProgress(1/N, detail = "Clusterize")
    concotion <- teabag %>%
      teabag_brew(min_cluster_size = input$min_cluster_size,
                  max_cluster_size = input$max_cluster_size)

    incProgress(1/N, detail = "Build")
    samovar <- concotion %>%
      concotion_pour(probability_calculation = input$probability_calculation,
                     inner_model = input$inner_model,
                     inter_model = input$inter_model,
                     cluster_connection = "mean")

    incProgress(1/N, detail = "Copy to environment")

    iter <<- iter + 1
    samovar_bases[[name]] <<- samovar
    updateSelectInput(session = session,
                      inputId = "db_use",
                      choices = names(samovar_bases),
                      selected = names(samovar_bases)[length(names(samovar_bases))])

    incProgress(1/N, detail = "Done")
    })
    if (is.character(samovar_catch)) {
      output$error_text <- renderText(paste0("Warning: ", samovar_catch))
    }
  }))

  # UX: load databases ----
  output$db_built <- renderUI({
    selectInput("db_use",
                label = NULL,
                choices = names(samovar_bases),
                selected = "1. Health_1500"
                )
    }
  )

  # UX: explore database ----
  output$composition_init <- renderPlotly({
    req(input$db_use)

    data <- teatree_bases[[input$db_use]]

    if(input$db_use %in% names(samovar_bases)) {
      if(!input$init_show_top) {
        top <- F
      } else {
        top <- input$init_top_count
      }
      data %>%
        viz_composition(type = input$init_type,
                        top = top, interactive = T)
    }
  })


  # UX: generate ----
  observeEvent(input$db_use, {

    idu <- input$db_use
    data <- teatree_bases[[idu]]

    output$downloadData1 <- downloadHandler(
      filename = function() {
        paste0("new-", idu, "-", Sys.Date(), ".csv")
      },
      content = function(con) {
        write.csv(data, con)})

    output$download_button_init <- renderUI(downloadButton("downloadData1", "Download generated data"))

    output$generation_panel <- renderUI({
      mainPanel(
        h3("Generate"),
        tags$br(),
        fluidPage(
          column(3,
                 numericInput("number_of_samples", "Number of samples", value = 100,
                                 min = 1, max = 10000, step = 100),
                 checkboxInput("avoid_zero_generations", "Avoid zero-based predictions", value = T),
                 checkboxInput("gen_show_top", label = "Show only top species", value = T),
                 sliderInput("gen_top_count", width = "100%",
                             label = "Show top species",
                             min = 1, max = 50, value = 10),
                 fluidRow(selectInput("gen_type", label = NULL,
                                      choices = c("tile", "column"), selected = "tile")),
                 tags$br(),
                 actionButton("generate", "Generate samples"),
                 tags$br(),
                 uiOutput("download_button"),
                 textOutput("error_text_generation")),
          column(8, plotlyOutput("composition_gen")),
        ),
      )
    })
  })

  observeEvent(input$generate, {
    generate_error <- try({
      idu <- input$db_use
      inos <- input$number_of_samples

      withProgress(message = "Generate data", value = 0, {
        N <- 2
        incProgress(1/N, detail = "Generate data")
      new_data <- samovar_bases[[idu]] %>%
        samovar_boil(N = inos, avoid_zero_generations = F)
      incProgress(1/N, detail = "Done")
      })
      sp <- rownames(new_data$data)

      new_data <- new_data$data %>%
        data.frame()%>%
        apply(2, as.numeric) %>%
        data.frame()

      rownames(new_data) <- sp

      new_data[is.na(new_data)] <- 0
      new_data[new_data < 0] <- 0

      s1 <- apply(new_data, 1, function(x) sum(x) > 0)
      s2 <- apply(new_data, 2, function(x) sum(x) > 0)
      new_data <<- new_data[s1,s2]

      output$downloadData <- downloadHandler(
        filename = function() {
          paste0("new-", idu, "-", inos, "-", Sys.Date(), ".csv")
        },
        content = function(con) {
          write.csv(new_data, con)})

      output$download_button <- renderUI(downloadButton("downloadData", "Download generated data"))
      })

    if(is.character(generate_error)) {
      output$error_text_generation <- generate_error
    }
  })

  generated_check <- reactive({
    list(input$gen_show_top,input$gen_top_count,input$generate, input$gen_type)
  })

  output$composition_gen <- renderPlotly({
      req(input$generate)
      if(!input$gen_show_top) {
        top <- F
      } else {
        top <- input$gen_top_count
      }
      viz_composition(new_data,
                      type = input$gen_type,
                      top = top, interactive = T)
    })
}

shinyApp(ui, server)
