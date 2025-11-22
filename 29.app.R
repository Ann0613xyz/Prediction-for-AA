library(shiny)
library(xgboost)
library(Matrix)
library(tidyverse)

setwd="C:\\Users\\anxin\\Desktop\\integrated_model"
xgb_model_final <- readRDS("xgb_model_final.rds")
feature_genes <- readRDS("feature_genes.rds")
ui <- fluidPage(
  titlePanel("Prediction for alopecia areata (AA)-based on XGBoost"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Please input the expression levels of the following genes (normalized):"),
      numericInput("input_KRT83", 
                   paste0(feature_genes[1], " expression level"), 
                   value = 0, 
                   min = 0, max = 50), 
      
      numericInput("input_PPP1R1C", 
                   paste0(feature_genes[2], " expression level"), 
                   value = 0,
                   min = 0, max = 50),
      
      numericInput("input_PIRT", 
                   paste0(feature_genes[3], " expression level"), 
                   value = 0,
                   min = 0, max = 50),
      
      actionButton("predict_button", "make a prediction", class = "btn-primary")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("prediction result", 
                 h3("prediction result:"),
                 
                 verbatimTextOutput("prediction_prob"),

                 h3("diagnostic classification (threshold > 0.5):"),
                 verbatimTextOutput("prediction_class"),
                 
                 br(),

                 p("Note: 0 represents 'con' (control group), 1 represents 'AA' (disease group)."),
                 p("This model predicts the probability that this sample belongs to 'AA' (disease group).")
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
  prediction_output <- eventReactive(input$predict_button, {

    input_data <- data.frame(
      KRT83 = input$input_KRT83,
      PPP1R1C = input$input_PPP1R1C,
      PIRT = input$input_PIRT
    )

    pred_matrix <- Matrix(as.matrix(input_data), sparse = TRUE)
    dpred <- xgb.DMatrix(data = pred_matrix)

    prediction_prob <- predict(xgb_model_final, newdata = dpred)
    

    prediction_class <- ifelse(prediction_prob > 0.5, 1, 0)

    list(prob = prediction_prob, class = prediction_class)
  })
  

  output$prediction_prob <- renderPrint({
    req(prediction_output()) 
    prob <- prediction_output()$prob
    cat(sprintf("AA (disease group) predicted probability: %.4f", prob))
  })

  output$prediction_class <- renderPrint({
    req(prediction_output())
    class <- prediction_output()$class
    class_label <- ifelse(class == 1, "AA (disease group)", "con (control group)")
    cat(class_label)
  })
  
}


shinyApp(ui = ui, server = server)

# # ##https://yuai.shinyapps.io/XGboost_for_alopecia_areata/

