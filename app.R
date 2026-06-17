library(shiny)

source("metrics_2x2.R")

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      body {
        background-color: #f7f7f7;
        font-family: Arial, sans-serif;
        color: #222222;
      }

      .title {
        font-size: 28px;
        font-weight: 700;
        margin-top: 20px;
        margin-bottom: 5px;
      }

      .subtitle {
        font-size: 14px;
        color: #555555;
        margin-bottom: 20px;
      }

      .card {
        background-color: white;
        border-radius: 14px;
        padding: 22px;
        margin-bottom: 22px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.10);
        border: 1px solid #e5e5e5;
      }

      .card-title {
        font-size: 18px;
        font-weight: 700;
        margin-bottom: 12px;
      }

      .small-text {
        font-size: 13px;
        color: #555555;
        margin-top: 10px;
      }

      .metric-card {
        background-color: white;
        border-radius: 12px;
        padding: 18px;
        margin-bottom: 14px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.08);
        border: 1px solid #e5e5e5;
        min-height: 100px;
      }

      .metric-title {
        font-size: 12px;
        font-weight: 700;
        text-transform: uppercase;
        color: #777777;
        margin-bottom: 8px;
      }

      .metric-value {
        font-size: 28px;
        font-weight: 700;
        color: #006d6f;
      }

      table {
        font-size: 13px;
      }

      .form-group {
        margin-bottom: 12px;
      }

      .control-label {
        font-size: 13px;
        font-weight: 600;
      }
    "))
  ),
  
  div(class = "title", "Binary classifier evaluation using 2×2 contingency tables"),
  div(
    class = "subtitle",
    "Interactive tool for calculating evaluation metrics, standard errors and visual summaries from a 2×2 contingency table."
  ),
  
  fluidRow(
    
    column(
      width = 4,
      
      div(
        class = "card",
        div(class = "card-title", "Contingency table"),
        
        fluidRow(
          column(4, ""),
          column(4, strong("Reference +")),
          column(4, strong("Reference −"))
        ),
        
        fluidRow(
          column(4, strong("Prediction +")),
          column(4, numericInput("a", NULL, value = NULL, min = 0)),
          column(4, numericInput("b", NULL, value = NULL, min = 0))
        ),
        
        fluidRow(
          column(4, strong("Prediction −")),
          column(4, numericInput("c", NULL, value = NULL, min = 0)),
          column(4, numericInput("d", NULL, value = NULL, min = 0))
        ),
        
        hr(),
        
        checkboxInput("se", "Compute standard error (SE)", TRUE),
        
        radioButtons(
          "se_method",
          "SE estimation method",
          choices = c(
            "Automatic: analytical when available, bootstrap otherwise" = "auto",
            "Bootstrap for all metrics" = "bootstrap"
          ),
          selected = "auto"
        ),
        
        numericInput(
          "B",
          "Bootstrap replications (B)",
          value = 1000,
          min = 100,
          step = 100
        ),
        
        div(
          class = "small-text",
          "Automatic mode uses analytical SE for proportion-based metrics when available and bootstrap for the remaining metrics. Bootstrap mode applies resampling to all metrics. Metric values are displayed with 3 decimals and SE values with 5 decimals."
        )
      ),
      
      div(
        class = "card",
        div(class = "card-title", "Summary indicators"),
        
        fluidRow(
          column(6, uiOutput("card_accuracy")),
          column(6, uiOutput("card_f1"))
        ),
        
        fluidRow(
          column(6, uiOutput("card_mcc")),
          column(6, uiOutput("card_kappa"))
        )
      )
    ),
    
    column(
      width = 8,
      
      div(
        class = "card",
        div(class = "card-title", "2×2 table"),
        tableOutput("tabla"),
        plotOutput("confusion_plot", height = 620)
      ),
      
      div(
        class = "card",
        tabsetPanel(
          
          tabPanel(
            "Basic measures",
            br(),
            tableOutput("basic"),
            plotOutput("basic_plot", height = 500)
          ),
          
          tabPanel(
            "Composite metrics",
            br(),
            tableOutput("composite"),
            plotOutput("composite_plot", height = 500)
          ),
          
          tabPanel(
            "Association",
            br(),
            tableOutput("association"),
            plotOutput("association_plot", height = 500)
          ),
          
          tabPanel(
            "Chance-corrected agreement",
            br(),
            tableOutput("agreement"),
            plotOutput("agreement_plot", height = 500)
          ),
          
          tabPanel(
            "Full summary",
            br(),
            tableOutput("summary")
          )
        )
      )
    )
  )
)

server <- function(input, output){
  
  tab <- reactive({
    vals <- c(input$a, input$b, input$c, input$d)
    validate(need(all(!is.na(vals)), "Please fill all cells."))
    
    matrix(
      c(input$a, input$b,
        input$c, input$d),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        Prediction = c("Positive", "Negative"),
        Reference  = c("Positive", "Negative")
      )
    )
  })
  
  output$tabla <- renderTable(tab(), rownames = TRUE)
  
  format_table <- function(x){
    x <- as.data.frame(x)
    
    if ("statistic" %in% names(x)) {
      x$statistic <- formatC(
        as.numeric(x$statistic),
        format = "f",
        digits = 3
      )
    }
    
    if ("se" %in% names(x)) {
      x$se <- formatC(
        as.numeric(x$se),
        format = "f",
        digits = 5
      )
    }
    
    x
  }
  
  method_basic <- reactive({
    if (input$se_method == "auto") {
      "analytic"
    } else {
      "bootstrap"
    }
  })
  
  method_other <- reactive({
    "bootstrap"
  })
  
  method_basic_label <- reactive({
    if (!input$se) {
      "Not computed"
    } else if (input$se_method == "auto") {
      "Analytical"
    } else {
      "Bootstrap"
    }
  })
  
  method_other_label <- reactive({
    if (!input$se) {
      "Not computed"
    } else {
      "Bootstrap"
    }
  })
  
  add_labels_and_error_bars <- function(bp, values, se_vals = NULL, offset = 0.04){
    
    if (!is.null(se_vals) &&
        input$se &&
        all(is.finite(se_vals)) &&
        any(se_vals > 0)) {
      
      upper <- values + 1.96 * se_vals
      
      arrows(
        x0 = bp,
        y0 = values - 1.96 * se_vals,
        x1 = bp,
        y1 = upper,
        angle = 90,
        code = 3,
        length = 0.04
      )
      
      label_y <- upper + offset
      
    } else {
      
      label_y <- values + offset
    }
    
    text(
      x = bp,
      y = label_y,
      labels = round(values, 3),
      cex = 0.85
    )
  }
  metric_card <- function(title, value){
    div(
      class = "metric-card",
      div(class = "metric-title", title),
      div(class = "metric-value", value)
    )
  }
  
  basic_vals <- reactive({
    t <- tab()
    
    x <- rbind(
      "Sensitivity" = unlist(ind_sens(t, input$se, method_basic(), "haldane", input$B)),
      "Specificity" = unlist(ind_spec(t, input$se, method_basic(), "haldane", input$B)),
      "PPV" = unlist(ind_ppv(t, input$se, method_basic(), "haldane", input$B)),
      "NPV" = unlist(ind_npv(t, input$se, method_basic(), "haldane", input$B)),
      "FPR" = unlist(ind_fpr(t, input$se, method_basic(), "haldane", input$B)),
      "FNR" = unlist(ind_fnr(t, input$se, method_basic(), "haldane", input$B)),
      "Accuracy" = unlist(ind_acc(t, input$se, method_basic(), "haldane", input$B))
    )
    
    x <- as.data.frame(x)
    x$`SE method` <- method_basic_label()
    x
  })
  
  composite_vals <- reactive({
    t <- tab()
    
    x <- rbind(
      "F1-score" = unlist(ind_f1(t, input$se, method_other(), "haldane", input$B)),
      "Balanced accuracy" = unlist(ind_bacc(t, input$se, method_other(), "haldane", input$B)),
      "Youden index" = unlist(ind_youden(t, input$se, method_other(), "haldane", input$B)),
      "Markedness" = unlist(ind_markedness(t, input$se, method_other(), "haldane", input$B))
    )
    
    x <- as.data.frame(x)
    x$`SE method` <- method_other_label()
    x
  })
  
  association_vals <- reactive({
    t <- tab()
    
    x <- rbind(
      "Phi" = unlist(ind_phi(t, input$se, method_other(), "haldane", input$B)),
      "MCC" = unlist(ind_mcc(t, input$se, method_other(), "haldane", input$B)),
      "Yule Q" = unlist(ind_yule_q(t, input$se, method_other(), "haldane", input$B)),
      "Yule Y" = unlist(ind_yule_y(t, input$se, method_other(), "haldane", input$B)),
      "Diagnostic OR" = unlist(ind_dor(t, input$se, method_other(), "haldane", input$B))
    )
    
    x <- as.data.frame(x)
    x$`SE method` <- method_other_label()
    x
  })
  
  agreement_vals <- reactive({
    t <- tab()
    
    x <- rbind(
      "Cohen's kappa" = unlist(ind_kappa(t, input$se, method_other(), "haldane", input$B)),
      "Gwet AC1" = unlist(ind_ac1(t, input$se, method_other(), "haldane", input$B)),
      "Delta" = unlist(ind_delta(t, input$se, method_other(), "haldane", input$B))
    )
    
    x <- as.data.frame(x)
    x$`SE method` <- method_other_label()
    x
  })
  
  output$card_accuracy <- renderUI({
    m <- as.data.frame(basic_vals())
    metric_card("Accuracy", round(as.numeric(m["Accuracy", "statistic"]), 3))
  })
  
  output$card_f1 <- renderUI({
    m <- as.data.frame(composite_vals())
    metric_card("F1-score", round(as.numeric(m["F1-score", "statistic"]), 3))
  })
  
  output$card_mcc <- renderUI({
    m <- as.data.frame(association_vals())
    metric_card("MCC", round(as.numeric(m["MCC", "statistic"]), 3))
  })
  
  output$card_kappa <- renderUI({
    m <- as.data.frame(agreement_vals())
    metric_card("Kappa", round(as.numeric(m["Cohen's kappa", "statistic"]), 3))
  })
  
  output$basic <- renderTable(format_table(basic_vals()), rownames = TRUE)
  
  output$basic_plot <- renderPlot({
    par(mar = c(11, 5, 4, 2) + 0.1)
    
    m <- as.data.frame(basic_vals())
    
    bp <- barplot(
      as.numeric(m$statistic),
      ylim = c(0, 1.20),
      las = 2,
      ylab = "Value",
      names.arg = rownames(m),
      main = "Basic evaluation measures",
      col = "#43958f",
      border = "#2f6f6a",
      cex.names = 0.85
    )
    
    add_labels_and_error_bars(
      bp = bp,
      values = as.numeric(m$statistic),
      se_vals = as.numeric(m$se),
      offset = 0.04
    )
  })
  
  output$composite <- renderTable(format_table(composite_vals()), rownames = TRUE)
  
  output$composite_plot <- renderPlot({
    par(mar = c(11, 5, 4, 2) + 0.1)
    
    m <- as.data.frame(composite_vals())
    
    bp <- barplot(
      as.numeric(m$statistic),
      ylim = c(0, 1.20),
      las = 2,
      ylab = "Value",
      names.arg = rownames(m),
      main = "Composite metrics",
      col = "#43958f",
      border = "#2f6f6a",
      cex.names = 0.85
    )
    
    add_labels_and_error_bars(
      bp = bp,
      values = as.numeric(m$statistic),
      se_vals = as.numeric(m$se),
      offset = 0.04
    )
  })
  
  output$association <- renderTable(format_table(association_vals()), rownames = TRUE)
  
  output$association_plot <- renderPlot({
    par(mar = c(11, 5, 4, 2) + 0.1)
    
    m <- as.data.frame(association_vals())
    
    values <- as.numeric(m$statistic)
    se_vals <- as.numeric(m$se)
    
    upper <- values + ifelse(is.finite(se_vals), 1.96 * se_vals, 0)
    lower <- values - ifelse(is.finite(se_vals), 1.96 * se_vals, 0)
    
    ymax <- max(upper, values, na.rm = TRUE)
    ymin <- min(lower, values, na.rm = TRUE)
    
    bp <- barplot(
      values,
      ylim = c(min(0, ymin) * 1.20, ymax * 1.25),
      las = 2,
      ylab = "Value",
      names.arg = rownames(m),
      main = "Association measures",
      col = "#43958f",
      border = "#2f6f6a",
      cex.names = 0.85
    )
    
    abline(h = 0, lty = 2)
    
    text(
      x = bp,
      y = values + 0.04 * ymax,
      labels = round(values, 3),
      cex = 0.85
    )
    
    if (input$se &&
        all(is.finite(se_vals)) &&
        any(se_vals > 0)) {
      
      arrows(
        x0 = bp,
        y0 = values - 1.96 * se_vals,
        x1 = bp,
        y1 = values + 1.96 * se_vals,
        angle = 90,
        code = 3,
        length = 0.04
      )
    }
  })
  
  output$agreement <- renderTable(format_table(agreement_vals()), rownames = TRUE)
  
  output$agreement_plot <- renderPlot({
    par(mar = c(11, 5, 4, 2) + 0.1)
    
    m <- as.data.frame(agreement_vals())
    
    bp <- barplot(
      as.numeric(m$statistic),
      ylim = c(-1, 1.20),
      las = 2,
      ylab = "Value",
      names.arg = rownames(m),
      main = "Chance-corrected agreement measures",
      col = "#43958f",
      border = "#2f6f6a",
      cex.names = 0.85
    )
    
    abline(h = 0, lty = 2)
    
    add_labels_and_error_bars(
      bp = bp,
      values = as.numeric(m$statistic),
      se_vals = as.numeric(m$se),
      offset = 0.05
    )
  })
  
  output$summary <- renderTable({
    x <- rbind(
      basic_vals(),
      composite_vals(),
      association_vals(),
      agreement_vals()
    )
    
    format_table(x)
  }, rownames = TRUE)
  
  output$confusion_plot <- renderPlot({
    
    t <- tab()
    
    TP <- t["Positive", "Positive"]
    FP <- t["Positive", "Negative"]
    FN <- t["Negative", "Positive"]
    TN <- t["Negative", "Negative"]
    
    par(mar = c(4, 7, 3, 2), xpd = NA)
    
    plot(
      x = c(-0.9, 2.05),
      y = c(-0.35, 2.05),
      type = "n",
      axes = FALSE,
      xlab = "",
      ylab = "",
      asp = 1
    )
    
    rect(0, 1, 1, 2, col = "#d9eeeb", border = "black", lwd = 2)
    rect(1, 1, 2, 2, col = "#f4f4f4", border = "black", lwd = 2)
    rect(0, 0, 1, 1, col = "#f4f4f4", border = "black", lwd = 2)
    rect(1, 0, 2, 1, col = "#d9eeeb", border = "black", lwd = 2)
    
    text(0.5, 1.5, TP, cex = 2.0, font = 2, col = "#004f52")
    text(1.5, 1.5, FP, cex = 2.0, font = 2, col = "#004f52")
    text(0.5, 0.5, FN, cex = 2.0, font = 2, col = "#004f52")
    text(1.5, 0.5, TN, cex = 2.0, font = 2, col = "#004f52")
    
    text(0.5, -0.12, "Reference +", cex = 1.2, col = "gray35")
    text(1.5, -0.12, "Reference −", cex = 1.2, col = "gray35")
    text(1, -0.28, "Reference", cex = 1.35, font = 2)
    
    text(-0.28, 1.5, "Prediction +", cex = 1.2, col = "gray35")
    text(-0.28, 0.5, "Prediction −", cex = 1.2, col = "gray35")
    text(-0.55, 1, "Prediction", cex = 1.35, font = 2, srt = 90)
    
    title("Confusion matrix", cex.main = 1.6, font.main = 2)
  })
}

shinyApp(ui, server)

