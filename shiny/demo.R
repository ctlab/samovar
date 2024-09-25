#prepare
library(tidyverse)
library(shiny)
library(plotly)
library(samovaR)

# UI ----
ui <- fluidPage (
  h1("samovaR: Microbiome generator v.1.0"),
  tags$hr(),
  textOutput("description"),
  tags$hr(),
)

# UX ----

server <- function (input, output) {

  # UX: stable ----
  output$description <- renderText("For generation, select properties of your database and than select it, or load one of pre-built samova.R databases")

  # UX: render UI for database selection ----


  # create a database ----


  # UX: load databases ----


  # UX: explore database ----


  # UX: generate ----

}

shinyApp(ui, server)
