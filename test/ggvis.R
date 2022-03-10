## ggvis

library(ggvis)
library(dplyr)

# 入门
mtcars %>% 
  ggvis(~wt, ~mpg) %>%
  layer_points(fill = ~factor(cyl), 
               size := 25, shape := "diamond", 
               stroke := "red") %>% 
  group_by(cyl) %>%
  layer_model_predictions(model = "lm", se = TRUE)

# 动态输入
mtcars %>% 
  ggvis(~wt, ~mpg, 
        size := input_slider(10, 100),
        opacity := input_slider(0, 1)
  ) %>% 
  layer_points()

# 鼠标提示
mtcars %>% ggvis(~wt, ~mpg) %>% 
  layer_points() %>% 
  add_tooltip(function(df) df$wt)

# use in Shiny
# ui.R
library(shiny)

# Define UI for miles per gallon application
shinyUI(sidebarLayout(
  sidebarPanel(
    sliderInput("size", "Area", 10, 1000, 500)
  ),
  mainPanel(
    uiOutput("ggvis_ui"),
    ggvisOutput("ggvis")
  )
))
# server.R
library(shiny)
library(ggvis)
library(dplyr)
mpgData <- mtcars
mpgData$am <- factor(mpgData$am, labels = c("Automatic", "Manual"))
# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  input_size <- reactive(input$size)
  
  mtcars %>% 
    ggvis(~disp, ~mpg, size := input_size) %>%
    layer_points() %>%
    bind_shiny("ggvis", "ggvis_ui")
})