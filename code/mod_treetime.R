# =====================================================================
# Module: TreeTime TDRP Explorer
# =====================================================================

json_safe_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  unname(as.numeric(x))
}
json_safe_character <- function(x) {
  if (is.null(x)) return(character(0))
  unname(as.character(x))
}
json_safe_scalar_numeric <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  y <- suppressWarnings(as.numeric(x[[1]]))
  if (length(y) == 0 || is.na(y)) default else unname(y)
}

mod_treetime_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h5("TreeTime ML Analyzer"),
      p("Explore macro-evolutionary trends from maximum-likelihood phylogenetic trees.", style = "font-size: 0.85em; color: gray;"),
      fileInput(ns("tt_tree"), "Upload TreeTime Tree (.nexus/.nwk)", accept = c(".nexus", ".nwk")),
      numericInput(ns("tt_min_time"), "Min Time-scale (t):", value = 0.0, step = 0.1),
      sliderInput(ns("tt_max_rate_q"), "Max Rate Quantile:", min = 0.9, max = 1.0, value = 0.99),
      selectInput(ns("tt_rate_col"), "Rate Column:", choices = c("Auto-detect", "mut_rate", "rate", "clock.rate")),
      
      div(
        style = "display: flex; gap: 8px; margin-bottom: 15px;",
        actionButton(ns("run_tt"), "Run TDRP Explorer", class = "btn-primary", style = "flex: 1;"),
        actionButton(ns("reset_tt"), "Reset", class = "btn-secondary", style = "flex: 0 0 auto;")
      ),
      
      hr(),
      h5("Robustness Diagnostic"),
      sliderInput(ns("tt_trim_q"), "Trim Extremes (Top/Bottom %):", min = 1, max = 15, value = 5, step = 1, post = "%"),
      actionButton(ns("run_tt_robustness"), "Run Sensitivity Check", class = "btn-primary", style = "width: 100%; margin-top: 5px;")
    ),
    mainPanel(
      uiOutput(ns("tt_summary")) %>% withSpinner(color = "#2C3E50"),
      plotlyOutput(ns("tt_plot"), height = "480px") %>% withSpinner(color = "#2C3E50"),
      hr(),
      div(
        h6(strong("Robustness Diagnostic Results"), style = "margin-bottom: 8px;"),
        uiOutput(ns("tt_robustness_ui")) %>% withSpinner(color = "#2C3E50")
      ),
      hr(),
      div(
        style = "display: flex; gap: 10px;",
        downloadButton(ns("dl_tt_csv"), "Download Parameters (.csv)", class = "btn btn-outline-primary btn-sm"),
        downloadButton(ns("dl_tt_pdf"), "Export Vector PDF", class = "btn btn-outline-danger btn-sm")
      )
    )
  )
}

mod_treetime_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    tt_state <- reactiveValues(
      df = NULL,
      model = NULL,
      current_plot = NULL,
      params_df = NULL
    )
    
    observeEvent(input$reset_tt, {
      tt_state$df <- NULL
      tt_state$model <- NULL
      tt_state$current_plot <- NULL
      tt_state$params_df <- NULL
      output$tt_summary <- renderUI({ NULL })
      output$tt_plot <- renderPlotly({ NULL })
      output$tt_robustness_ui <- renderUI({ NULL })
      
      updateNumericInput(session, "tt_min_time", value = 0.0)
      updateSliderInput(session, "tt_max_rate_q", value = 0.99)
      updateSelectInput(session, "tt_rate_col", selected = "Auto-detect")
      updateSliderInput(session, "tt_trim_q", value = 5)
      
      showNotification("TreeTime Explorer workspace reset.", type = "message")
    })
    
    observeEvent(input$run_tt, {
      req(input$tt_tree)
      
      output$tt_robustness_ui <- renderUI({ NULL })
      
      tryCatch({
        parsed <- suppressWarnings(parse_beast_tree(input$tt_tree$datapath, origin = "treetime"))
        
        col <- if(input$tt_rate_col == "Auto-detect") NULL else input$tt_rate_col
        raw_df <- extract_node_features(parsed$data, col)
        df <- filter_extremes(raw_df, input$tt_min_time, input$tt_max_rate_q)
        
        model <- fit_pow_model(df)
        if (is.null(model)) stop("Model fitting failed to converge.")
        
        coefs_raw <- coef(model)
        mu0 <- json_safe_scalar_numeric(coefs_raw["mu_0"])
        muinf <- json_safe_scalar_numeric(coefs_raw["mu_infty"])
        lam <- json_safe_scalar_numeric(coefs_raw["lambda"])
        
        tt_state$df <- df
        tt_state$model <- model
        tt_state$params_df <- data.frame(
          Parameter = c("\u03bc\u2080 (Initial rate parameter)", "\u03bc\u221e (Long-term substitution rate)", "\u03bb (Decay Coefficient)"), 
          Value = c(mu0, muinf, lam)
        )
        
        output$tt_summary <- renderUI({
          div(class="alert alert-success", style="padding:10px; margin-bottom:15px;", 
              HTML(sprintf("<b>Parsed successfully!</b><br>Decay Coefficient (\u03bb) = %.3f<br>Initial rate parameter (\u03bc<sub>0</sub>) = %.3e<br>Long-term substitution rate (\u03bc<sub>\u221e</sub>) = %.3e", lam, mu0, muinf))
          )
        })
        
        t_min <- max(min(df$branch.length, na.rm = TRUE), 1e-5)
        t_max <- max(df$branch.length, na.rm = TRUE)
        t_seq <- exp(seq(log(t_min), log(t_max), length.out = 100))
        
        curve_df <- predict_pow_with_ci(model, t_seq)
        
        plot_df <- as.data.frame(df, stringsAsFactors = FALSE)
        plot_df$log_time <- log10(pmax(as.numeric(plot_df$branch.length), 1e-12))
        plot_df$log_rate <- log10(pmax(as.numeric(plot_df$rate), 1e-12))
        
        curve_df <- as.data.frame(curve_df, stringsAsFactors = FALSE)
        curve_df$log_time <- log10(pmax(as.numeric(curve_df$branch.length), 1e-12))
        curve_df$log_rate_pred <- log10(pmax(as.numeric(curve_df$rate_pred), 1e-12))
        curve_df$log_upper_ci <- log10(pmax(as.numeric(curve_df$upper_ci), 1e-12))
        curve_df$log_lower_ci <- log10(pmax(as.numeric(curve_df$lower_ci), 1e-12))
        
        output$tt_plot <- renderPlotly({
          plot_df_json <- plot_df
          curve_df_json <- curve_df
          
          plot_df_json$log_time <- json_safe_numeric(plot_df_json$log_time)
          plot_df_json$log_rate <- json_safe_numeric(plot_df_json$log_rate)
          curve_df_json$log_time <- json_safe_numeric(curve_df_json$log_time)
          curve_df_json$log_upper_ci <- json_safe_numeric(curve_df_json$log_upper_ci)
          curve_df_json$log_lower_ci <- json_safe_numeric(curve_df_json$log_lower_ci)
          curve_df_json$log_rate_pred <- json_safe_numeric(curve_df_json$log_rate_pred)
          
          node_values <- if("node" %in% names(plot_df_json)) json_safe_character(plot_df_json$node) else rep("", nrow(plot_df_json))
          
          hover_text <- vapply(seq_len(nrow(plot_df_json)), function(i) {
            paste0("Node: ", node_values[i], "<br>Log10(Time): ", signif(plot_df_json$log_time[i], 4), "<br>Log10(Rate): ", signif(plot_df_json$log_rate[i], 4))
          }, FUN.VALUE = character(1))
          
          param_annotation_text <- paste0(
            "<b>DeepPoW Parameters:</b><br>",
            "\u03bc<sub>0</sub> = ", format(mu0, scientific = TRUE, digits = 3), "<br>",
            "\u03bc<sub>\u221e</sub> = ", format(muinf, scientific = TRUE, digits = 3), "<br>",
            "\u03bb = ", signif(lam, 3)
          )
          
          p <- plotly::plot_ly() %>%
            plotly::add_trace(x = json_safe_numeric(curve_df_json$log_time), y = json_safe_numeric(curve_df_json$log_upper_ci), type = "scatter", mode = "lines", line = list(color = "rgba(0,0,0,0)"), showlegend = FALSE, hoverinfo = "none") %>%
            plotly::add_trace(x = json_safe_numeric(curve_df_json$log_time), y = json_safe_numeric(curve_df_json$log_lower_ci), type = "scatter", mode = "lines", fill = "tonexty", fillcolor = "rgba(217,83,79,0.15)", line = list(color = "rgba(0,0,0,0)"), name = "95% CI", hoverinfo = "none") %>%
            plotly::add_trace(x = json_safe_numeric(plot_df_json$log_time), y = json_safe_numeric(plot_df_json$log_rate), type = "scatter", mode = "markers", text = hover_text, hoverinfo = "text", marker = list(size = 6, color = "rgba(44,62,80,0.6)", line = list(color = "rgba(44,62,80,0.9)", width = 0.5)), name = "Branch Rates") %>%
            plotly::add_trace(x = json_safe_numeric(curve_df_json$log_time), y = json_safe_numeric(curve_df_json$log_rate_pred), type = "scatter", mode = "lines", line = list(color = "#D9534F", width = 2.5), name = "PoW Fit", hoverinfo = "none") %>%
            plotly::layout(
              title = list(text = "<b>Time-Dependent Rate Phenomenon (TreeTime Explorer)</b>", font = list(size = 12)),
              xaxis = list(title = "<b>Log<sub>10</sub>(Time-scale, years)</b>"),
              yaxis = list(title = "<b>Log<sub>10</sub>(Evolutionary Rate, μ)</b>"),
              annotations = list(list(x = 0.05, y = 0.95, xref = "paper", yref = "paper", text = as.character(param_annotation_text), showarrow = FALSE, bgcolor = "rgba(255,255,255,0.85)", bordercolor = "#CCCCCC", borderwidth = 1)),
              plot_bgcolor = "#FFFFFF", paper_bgcolor = "#FFFFFF", showlegend = TRUE,
              legend = list(x = 0.05, y = 0.12, xanchor = "left", yanchor = "bottom", bgcolor = "rgba(255,255,255,0.85)", bordercolor = "#CCCCCC", borderwidth = 1)
            ) %>%
            plotly::config(displaylogo = FALSE)
          
          tt_state$current_plot <- p
          p
        })
        
      }, error = function(e) {
        showNotification(paste("TreeTime Explorer Error:", e$message), type = "error")
      })
    })
    
    observeEvent(input$run_tt_robustness, {
      req(tt_state$df)
      
      tryCatch({
        trim_val <- input$tt_trim_q / 100
        robust_res <- robustness_pow_model(tt_state$df, trim_quantile = trim_val)
        
        if (is.null(robust_res) || is.null(robust_res$comparison_table)) {
          stop("Robustness check failed. The trimmed dataset might be too small or failed to converge.")
        }
        
        output$tt_robustness_ui <- renderUI({
          tbl <- robust_res$comparison_table
          
          tbl[] <- lapply(tbl, function(col) {
            if (is.character(col) || is.factor(col)) {
              col <- as.character(col)
              col <- gsub("mu_0", "\u03bc\u2080", col)
              col <- gsub("mu_infty", "\u03bc\u221e", col)
              col <- gsub("lambda", "\u03bb", col)
            }
            col
          })
          
          shift_val <- robust_res$lambda_shift
          full_mu_0 <- as.numeric(tbl$Full_Data[1])
          full_mu_infty <- as.numeric(tbl$Full_Data[2])
          full_lambda <- as.numeric(tbl$Full_Data[3])
          
          flat_reason <- ""
          if (is.na(full_lambda) || full_lambda <= 1.01e-4) {
            flat_reason <- "\u03bb is at the lower fitting bound (10⁻⁴; no detectable decay)"
          } else if (is.na(full_mu_0) || is.na(full_mu_infty) || full_mu_0 <= full_mu_infty) {
            flat_reason <- "\u03bc\u2080 \u2264 \u03bc\u221e (inverse or flat trajectory)"
          }
          is_paradoxical_or_flat <- (nchar(flat_reason) > 0)
          
          if (is.na(shift_val)) {
            alert_div <- div(class = "alert alert-danger", style = "padding: 10px; font-size: 0.85em; margin-top: 10px;", HTML(sprintf("<b>Convergence Failed:</b> The model could not converge after trimming %d%% outliers. The TDRP signal might be entirely driven by extreme noise.", input$tt_trim_q)))
          } else if (!is.na(full_lambda) && full_lambda > 150) {
            alert_div <- div(class = "alert alert-danger", style = "padding: 10px; font-size: 0.85em; margin-top: 10px;", HTML(sprintf("<b>Ultra-Rapid Decay Warning:</b> The decay coefficient (\u03bb = %.1f) exceeds the heuristic artifact-screening threshold (>150). This may indicate a fitting artifact or sensitivity to sampling structure and warrants further investigation.", full_lambda)))
          } else if (is_paradoxical_or_flat && shift_val < 0.20) {
            alert_div <- div(class = "alert alert-secondary", style = "padding: 10px; font-size: 0.85em; margin-top: 10px;", HTML(sprintf("<b>Stable but Flat/Inverse:</b> The parameter space is stable (shift %.1f%%), but <b>%s</b>. This indicates a lack of genuine TDRP signal on this ML tree.", shift_val * 100, flat_reason)))
          } else if (shift_val < 0.20) {
            alert_div <- div(class = "alert alert-success", style = "padding: 10px; font-size: 0.85em; margin-top: 10px;", HTML(sprintf("<b>Globally Robust:</b> The decay coefficient (\u03bb) shifted by only %.1f%% after removing outliers.", shift_val * 100)))
          } else if (shift_val < 1.00) {
            alert_div <- div(class = "alert alert-warning", style = "padding: 10px; font-size: 0.85em; margin-top: 10px;", HTML(sprintf("<b>Moderately Sensitive:</b> The decay coefficient (\u03bb) shifted by %.1f%%. Please manually review terminal branches.", shift_val * 100)))
          } else {
            alert_div <- div(class = "alert alert-danger", style = "padding: 10px; font-size: 0.85em; margin-top: 10px;", HTML(sprintf("<b>Severe Artifact Detected:</b> The decay coefficient (\u03bb) shifted by %.1f%%. The fitted TDRP pattern may be strongly influenced by extreme branch-rate observations.", shift_val * 100)))
          }
          
          tagList(
            tags$table(
              class = "table table-sm table-bordered",
              style = "background-color: #ffffff; font-size: 0.9em;",
              tags$thead(tags$tr(lapply(colnames(tbl), function(col) tags$th(as.character(col))))),
              tags$tbody(lapply(seq_len(nrow(tbl)), function(i) {
                tags$tr(lapply(tbl[i, , drop = FALSE], function(val) tags$td(as.character(val))))
              }))
            ),
            alert_div
          )
        })
      }, error = function(e) {
        showNotification(paste("Robustness Check Error:", e$message), type = "error")
      })
    })
    
    output$dl_tt_csv <- downloadHandler(
      filename = function() { paste0("TreeTime_PoW_Params_", Sys.Date(), ".csv") },
      content = function(file) {
        req(tt_state$params_df)
        write.csv(tt_state$params_df, file, row.names = FALSE)
      }
    )
    
    output$dl_tt_pdf <- downloadHandler(
      filename = function() { paste0("TreeTime_TDRP_Plot_", Sys.Date(), ".pdf") },
      content = function(file) {
        req(tt_state$df, tt_state$model)
        export_pow_pdf(file, tt_state$df, tt_state$model)
      }
    )
    
  })
}
