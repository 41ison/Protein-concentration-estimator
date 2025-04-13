# This Shiny app allows users to upload BSA standard and sample absorbance data,
# perform blank corrections, fit a regression model, and calculate protein concentrations.
# The results can be visualized in plots and downloaded as a CSV file.

# requireNamespace is used to check if the package is installed
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!requireNamespace("ggtext", quietly = TRUE)) {
  install.packages("ggtext")
}

# load required packages
library(shiny)
library(tidyverse)
library(ggtext)

# Define UI for the application
ui <- fluidPage(
  titlePanel("Protein Concentration Estimator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("stdFile", "Upload BSA Standard Data (CSV)", accept = ".csv"),
      fileInput("sampleFile", "Upload Sample Absorbance Data (CSV)", accept = ".csv"),
      numericInput("dilution", "Sample dilution factor", value = 1, min = 1),
      textInput("blankValues", "Blank absorbance values (comma-separated)", value = "0.069, 0.073, 0.089"),
      selectInput("fitModel", "Regression Model", choices = c("Linear", "2nd-order Polynomial")),
      downloadButton("downloadData", "Download Results")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Standard Curve", plotOutput("stdPlot")),
        tabPanel("Sample Concentrations", 
                 plotOutput("samplePlot"),
                 tableOutput("sampleTable"))
      )
    )
  )
)

server <- function(input, output) {
  # Blank correction
  blank <- reactive({
    blanks <- as.numeric(unlist(strsplit(input$blankValues, ",")))
    mean(blanks, na.rm = TRUE)
  })
  
  observe({
    if (!is.null(input$stdFile)) {
      std_preview <- read_csv(input$stdFile$datapath, show_col_types = FALSE)
      print("Standard File Preview:")
      print(head(std_preview))
    }
    if (!is.null(input$sampleFile)) {
      sample_preview <- read_csv(input$sampleFile$datapath, show_col_types = FALSE)
      print("Sample File Preview:")
      print(head(sample_preview))
    }
  })
  
  # Standard data input and long format transformation
  bsa_data <- reactive({
    req(input$stdFile)
    std_df <- read_csv(input$stdFile$datapath, show_col_types = FALSE)
    std_df_long <- std_df %>%
      mutate(across(where(is.numeric), ~ . - blank())) %>%
      pivot_longer(-BSA_ng, names_to = "Std", values_to = "Abs")
    std_df_long
  })
  
  # Reactive model fitting
  linear_model <- reactive({
    req(bsa_data())
    if (input$fitModel == "Linear") {
      lm(Abs ~ BSA_ng, data = bsa_data())
    } else {
      lm(Abs ~ poly(BSA_ng, 2, raw = TRUE), data = bsa_data())
    }
  })
  
  # Plot of standard curve
  output$stdPlot <- renderPlot({
    req(bsa_data(), linear_model())
    
    model <- linear_model()
    r2 <- round(summary(model)$r.squared, 3)
    coeffs <- round(coef(model), 5)
    
    formula_text <- if (input$fitModel == "Linear") {
      paste0("y = ", coeffs[1], " + ", coeffs[2], "x")
    } else {
      paste0("y = ", coeffs[1], " + ", coeffs[2], "x + ", coeffs[3], "x²")
    }
    
    # Defining the formula outside the plot call
    formula_to_use <- if (input$fitModel == "Linear") {
      y ~ x
    } else {
      y ~ poly(x, 2, raw = TRUE)
    }
    
    bsa_data() %>%
      ggplot(aes(x = BSA_ng, y = Abs)) +
      geom_point(size = 3, alpha = 0.5) +
      stat_smooth(method = "lm",
                  formula = formula_to_use,
                  se = TRUE, linewidth = 1, color = "firebrick") +
      scale_x_continuous(breaks = c(125, 250, 500, 750, 1000, 1500, 2000)) +
      annotate("text",
               x = min(bsa_data()$BSA_ng),
               y = max(bsa_data()$Abs),
               label = paste0("R² = ", r2, "\n", formula_text),
               size = 5, fontface = "bold", hjust = 0) +
      labs(x = "BSA (ng/µL)", y = "Absorbance (562 nm)") +
      theme_bw() +
      theme(text = element_text(size = 18),
            axis.text = element_text(face = "bold", color = "black"),
            axis.ticks = element_line(color = "black"))
  })
  
  # Sample processing
  sample_results <- reactive({
    req(input$sampleFile, linear_model())
    sample_df <- read_csv(input$sampleFile$datapath, show_col_types = FALSE)
    
    sample_df_adj <- sample_df %>%
      mutate(across(everything(), ~ . - blank())) %>%
      pivot_longer(cols = everything(), names_to = "Sample", values_to = "Abs") %>%
      group_by(Sample) %>%
      mutate(
        meanAbs = mean(Abs),
        sd = sd(Abs)
      ) %>%
      distinct(Sample, meanAbs, sd) %>%
      mutate(
        concentration = if (input$fitModel == "Linear") {
          (meanAbs - coef(linear_model())[1]) / coef(linear_model())[2]
        } else {
          # Solve quadratic: a*x^2 + b*x + c = 0
          a <- coef(linear_model())[3]
          b <- coef(linear_model())[2]
          c <- coef(linear_model())[1] - meanAbs
          (-b + sqrt(b^2 - 4*a*c)) / (2*a)
        },
        corrected_ug_uL = round(concentration * input$dilution / 1000, 2)
      )
    sample_df_adj
  })
  
  # Table output
  output$sampleTable <- renderTable({
    req(sample_results())
    sample_results() %>%
      select(Sample, corrected_ug_uL, sd) %>%
      rename("Protein Conc. (µg/µL)" = corrected_ug_uL, "SD" = sd) %>%
      mutate(across(everything(), ~ round(., 2)))
  })
  
  # Plot of sample concentrations with standard deviation bars
  output$samplePlot <- renderPlot({
    req(sample_results())
    sample_results() %>%
      ggplot(aes(x = Sample, y = corrected_ug_uL)) +
      geom_col(fill = "skyblue", width = 0.6) +
      geom_errorbar(aes(ymin = corrected_ug_uL - sd, ymax = corrected_ug_uL + sd),
                    width = 0.2, color = "black") +
      labs(x = "Sample", y = "Protein Conc. (µg/µL)") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(color = "black"))
  })
  
  # Download handler for results
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("protein_concentrations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(sample_results(), file)
    }
  )
}

# Run the application
shinyApp(ui, server)