# =====================================================================
# DeepPoW: Interactive Platform for Deep Evolutionary Rate Modeling
# Main Application Script (JSON-safe / Plotly-safe Version)
# =====================================================================

library(shiny)
library(bslib)
library(plotly)
library(DT)
library(shinycssloaders)
library(xml2)
library(ape)
library(nleqslv)
library(ggplot2)
library(shinyFiles) 

options(shiny.maxRequestSize = 500 * 1024^2)


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

json_safe_object <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.atomic(x) && !is.null(names(x))) return(as.list(x))
  x
}

# ---------------------------------------------------------------------
# About (UI & Server)
# ---------------------------------------------------------------------
mod_about_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    style = "max-width: 1000px; padding: 20px;",
    
    div(style = "text-align: center; margin-bottom: 30px;",
        h2(style = "font-weight: bold; color: #0F172A;", "DeepPoW v1.0"),
        p(style = "font-size: 16px; color: #475569;", "Presented By DeepPoW Development Team")
    ),
    
    hr(style = "margin-bottom: 30px;"),
    
    fluidRow(
      column(6,
             div(style = "background-color: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px; padding: 20px; height: 100%;",
                 h4(style = "color: #1E293B; font-weight: bold; margin-top: 0;", "About & Contact"),
                 p(style = "font-size: 14px; color: #334155;", 
                   "DeepPoW is a rigorous and interactive platform designed to diagnose and calibrate time-dependent evolutionary rate patterns in viral evolution."),
                 tags$ul(style = "padding-left: 20px; font-size: 14px; color: #334155;",
                   tags$li(tags$strong("Project Director: "), "Prof. Wu"),
                   tags$li(tags$strong("Lead Developer: "), "Dr. Gao"),
                   tags$li(tags$strong("Contact Email: "), "raindy@fafu.edu.cn (F. Gao)")
                 )
             )
      ),
      
      column(6,
             div(style = "background-color: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px; padding: 20px; height: 100%;",
                 h4(style = "color: #1E293B; font-weight: bold; margin-top: 0;", "Docs & Code"),
                 p(style = "font-size: 14px; color: #334155;", "Browse the full user manual, or view the source code, dependencies, and local deployment guide on GitHub."),
                 div(style = "margin-top: 15px; display: flex; gap: 10px;",
                     tags$a(href = "guide.pdf", target = "_blank", class = "btn btn-outline-primary btn-sm", " Quick Guide (PDF)"),
                     tags$a(href = "manual.pdf", target = "_blank", class = "btn btn-outline-primary btn-sm", " User Manual (PDF)"),
                     tags$a(href = "https://github.com/BioEasy/DeepPoW", target = "_blank", class = "btn btn-outline-secondary btn-sm", "GitHub Repository")
                 ),
                 p(style = "font-size: 12px; color: #64748B; margin-top: 15px;", "* Please refer to README.md for full tech stack and dependency details.")
             )
      )
    ),
    
    br(),
    
    div(style = "background-color: #EFF6FF; border: 1px solid #BFDBFE; border-radius: 8px; padding: 20px;",
        h4(style = "color: #1E293B; font-weight: bold; margin-top: 0;", "How to Cite"),
        p(style = "font-size: 14px; color: #1E293B; margin-bottom: 15px;", 
          "Please cite our article if you use DeepPoW in your research."),
        div(style = "background-color: #FFFFFF; border-left: 4px solid #3B82F6; padding: 10px 15px; font-family: monospace; font-size: 13px; color: #334155; margin-bottom: 15px;",
            "Zhang, S., Shen, J., Jin, S., Du, Z., Wu, Z. and Gao, F. (2026). DeepPoW: A rigorous and interactive computational framework for calibrating the time-dependent rate phenomenon in viral evolution. Manuscript in preparation."
        ),
        tags$a(href = "citation.ris", download = "DeepPoW_Citation.ris", class = "btn btn-primary btn-sm", "Download DeepPoW Citation (.RIS)"),
    hr(style = "margin: 20px 0; border-top: 1px solid #93C5FD;"),
    p(style = "font-size: 14px; color: #1E293B; margin-bottom: 10px;", 
          "DeepPoW implements the mechanistic PoW model for time-dependent rate calibration. If your analysis utilizes this framework, please also cite the foundational theoretical work:"),
        div(style = "background-color: #F8FAFC; border-left: 4px solid #94A3B8; padding: 10px 15px; font-family: monospace; font-size: 13px; color: #475569; margin-bottom: 0;",
            "Ghafari, M., Simmonds, P., Pybus, O. G. and Katzourakis, A. (2021). A mechanistic evolutionary model explains the time-dependent pattern of substitution rates in viruses. Current Biology, 31, 4689–4696.e5."
        )
    ),
    div(style = "text-align: center; margin-top: 40px; font-size: 12px; color: #94A3B8;",
        "Copyright (C)2026 DeepPoW Development Team. Released under the MIT License."
    )
  )
}

mod_about_server <- function(id) {
  moduleServer(id, function(input, output, session) {
  })
}


code_dir <- "code"
if (!dir.exists(code_dir)) {
  stop(paste("Critical Error: The directory '", code_dir, "' does not exist!"))
}

module_files <- list.files(
  code_dir,
  pattern = "\\.R$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(module_files) == 0) {
  stop("Critical Error: No module files matching '.R' were found in 'code/'!")
}

for (f in module_files) {
  source(f)
  message("Successfully loaded module: ", basename(f))
}

# ==========================================
# 2. UI
# ==========================================

ui <- page_navbar(
  title = "DeepPoW",
  bg = "#2C3E50",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    secondary = "#34495E"
  ),

  header = tags$head(
    tags$style(HTML("
      .navbar {
        --bs-navbar-color: #FFFFFF !important;
        --bs-navbar-hover-color: #FFFFFF !important;
        --bs-navbar-active-color: #FFFFFF !important;
        --bs-navbar-brand-color: #FFFFFF !important;
        --bs-navbar-brand-hover-color: #FFFFFF !important;
      }
      .btn-primary, .btn-primary:hover, .btn-primary:active,
      .btn-primary:focus {
        background-color: #2C3E50 !important;
        border-color: #2C3E50 !important;
        color: #FFFFFF !important;
      }
      .btn-secondary, .btn-secondary:hover, .btn-secondary:active,
      .btn-secondary:focus {
        background-color: #6C757D !important;
        border-color: #6C757D !important;
        color: #FFFFFF !important;
      }
      .btn-file, .btn-file:hover, .input-group .btn {
        background-color: #34495E !important;
        border-color: #34495E !important;
        color: #FFFFFF !important;
      }
      .alert-success, .card.border-success {
        background-color: #EBF5FB !important;
        border-color: #2C3E50 !important;
        color: #1A252F !important;
      }
      .card.border-success > .card-header.bg-success {
        background-color: #2C3E50 !important;
        color: #FFFFFF !important;
        border-color: #2C3E50 !important;
      }
      .alert-info {
        background-color: #F4F6F7 !important;
        border-color: #34495E !important;
        color: #2C3E50 !important;
      }
      .nav-tabs .nav-link.active, .nav-underline .nav-link.active,
      .nav-pills .nav-link.active {
        color: #2C3E50 !important;
        border-bottom: 3px solid #2C3E50 !important;
        font-weight: bold !important;
        background-color: transparent !important;
      }
      .nav-link { color: #555555 !important; }
      .nav-link:hover { color: #2C3E50 !important; }
      .navbar-nav .nav-link {
        color: #FFFFFF !important;
        opacity: 1.0 !important;
      }
      .navbar-nav .nav-link:hover {
        color: #FFFFFF !important;
        opacity: 1.0 !important;
        background-color: rgba(255, 255, 255, 0.15) !important;
        border-radius: 4px;
      }
      .navbar-nav .nav-link.active {
        color: #FFFFFF !important;
        opacity: 1.0 !important;
        font-weight: bold !important;
        background-color: rgba(255, 255, 255, 0.25) !important;
        border-radius: 4px;
      }
      .navbar-brand {
        color: #FFFFFF !important;
        font-weight: bold !important;
      }
      .dataTables_wrapper .dataTables_paginate .paginate_button {
        background: #F8F9FA !important;
        border: 1px solid #DEE2E6 !important;
        color: #2C3E50 !important;
        border-radius: 4px !important;
        padding: 5px 10px !important;
        margin: 0 2px !important;
      }
      .dataTables_wrapper .dataTables_paginate .paginate_button.current,
      .dataTables_wrapper .dataTables_paginate .paginate_button.current:hover {
        background: #2C3E50 !important;
        border-color: #2C3E50 !important;
        color: #FFFFFF !important;
      }
      .dataTables_wrapper .dataTables_paginate .paginate_button:hover {
        background: #34495E !important;
        border-color: #34495E !important;
        color: #FFFFFF !important;
      }
      .pagination .page-item .page-link {
        color: #2C3E50 !important;
        background-color: #FFFFFF !important;
        border: 1px solid #DEE2E6 !important;
      }
      .pagination .page-item.active .page-link {
        background-color: #2C3E50 !important;
        border-color: #2C3E50 !important;
        color: #FFFFFF !important;
      }
      .pagination .page-item .page-link:hover {
        background-color: #34495E !important;
        border-color: #34495E !important;
        color: #FFFFFF !important;
      }
    "))
  ),


  nav_panel(
    title = "1. PoW Dynamics & Calibration",
    sidebarLayout(
      sidebarPanel(
        div(
          style = "background-color: #f8f9fa; padding: 12px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 15px;",
          radioButtons(
            "pow_mode",
            strong("Select Analysis Pipeline:"),
            choices = c(
              "Mode A: Single MCC Tree (Fast Diagnostics)" = "mode_a",
              "Mode B: Full MCMC Ensemble (HPD Calibration)" = "mode_b"
            ),
            selected = "mode_a"
          )
        ),

        conditionalPanel(
          condition = "input.pow_mode == 'mode_a'",
          h5("Step 1: Ingest Tree"),
          div(
            class = "alert alert-info",
            style = "padding: 10px; font-size: 0.85em; margin-bottom: 10px;",
            HTML("<b>Prerequisite:</b> Please upload an Annotated MCC tree (e.g., generated via BEAST/TreeAnnotator) containing branch-rate metadata. This mode is designed for rapid TDRP diagnostics prior to full MCMC calibration.")
          ),
          fileInput(
            "beast_tree",
            "Upload Time-calibrated Tree (.tree/.tre/.trees/.nexus)",
            accept = c(".tree", ".tre", ".trees", ".nexus")
          ),
          h5("Data Filtering"),
          numericInput(
            "min_time",
            "Min Time-scale (t):",
            value = 0.0,
            step = 0.1
          ),
          sliderInput(
            "max_rate_q",
            "Max Rate Quantile:",
            min = 0.9,
            max = 1.0,
            value = 0.99
          ),
          div(
            style = "display: flex; gap: 8px; margin-top: 10px; margin-bottom: 15px;",
            actionButton(
              "process_tree",
              "Parse Tree",
              class = "btn-primary",
              style = "flex: 1;"
            ),
            actionButton(
              "reset_tree",
              "Reset",
              class = "btn-secondary",
              style = "flex: 0 0 auto;"
            )
          ),
          hr(),
          h5("Step 2: Model & Fit"),
          selectInput(
            "rate_col",
            "Select Rate Column:",
            choices = c(
              "Auto-detect",
              "rate",
              "rate_median",
              "clock.rate",
              "branch_rate"
            )
          ),

          checkboxInput(
            "expert_mode",
            tags$span(style = "color: #D9534F; font-weight: bold; font-size: 0.9em;",
                      "Expert Mode: Disable k > 150 filter (intra-host data only)"),
            value = FALSE
          ),
          div(
            style = "display: flex; gap: 8px; margin-top: 10px;",
            actionButton(
              "run_pow",
              "Fit TDRP Model",
              class = "btn-primary",
              style = "flex: 1;"
            ),
            actionButton(
              "reset_pow_a",
              "Reset Fit",
              class = "btn-secondary",
              style = "flex: 0 0 auto;"
            )
          ),
          hr(),
          h5("Step 3: Robustness Diagnostic"),
          sliderInput("trim_quantile_ui", "Trim Extremes (Top/Bottom %):", min = 1, max = 15, value = 5, step = 1, post = "%"),
          actionButton("run_robustness", "Run Robustness Check", class = "btn-primary", style = "width: 100%; margin-top: 5px;")
        ),

        conditionalPanel(
          condition = "input.pow_mode == 'mode_b'",
          h5("Step 1: Ingest MCMC Ensemble"),
          div(
            class = "alert alert-info",
            style = "padding: 10px; font-size: 0.85em; margin-bottom: 15px;",
            HTML("<b>Prerequisite:</b> This calibration engine strictly processes Bayesian MCMC posterior ensembles (e.g., BEAST log and trees) to account for phylogenetic uncertainty across posterior trees.")
          ),
          fileInput(
            "rate_log",
            "1. Short-term substitution rate log file (*.log)",
            accept = c(".log", ".txt")
          ),
          selectInput(
            "sub_model",
            "Standard HKY or JC69 substitution model:",
            choices = c(
              "HKY85" = "HKY",
              "JC69 (Jukes-Cantor)" = "JC69"
            ),
            selected = "HKY"
          ),
          fileInput(
            "trees_file",
            "2. Ultrametric distance trees file (*.trees)",
            accept = c(".trees", ".tre", ".nex")
          ),
          conditionalPanel(
            condition = "input.sub_model == 'HKY'",
            fileInput(
              "hky_log",
              "3. Standard HKY substitution model log file (*.log)",
              accept = c(".log", ".txt")
            )
          ),
          hr(),
          h5("Step 2: Constraints & Run"),
          selectInput(
            "virus_type",
            HTML("Virus Genome Type (&mu;<sub>Max</sub>)"),
            choices = c(
              "RNA Viruses (3.65e-2)" = "RNA",
              "ssDNA Viruses (2.0e-2)" = "ssDNA",
              "dsDNA Viruses (3.0e-3)" = "dsDNA",
              "Custom" = "Other"
            ),
            selected = "RNA"
          ),
          conditionalPanel(
            condition = "input.virus_type == 'Other'",
            numericInput(
              "mu_max_b",
              "Custom muMax value:",
              value = 0.05,
              min = 1e-5,
              max = 1,
              step = 1e-3
            )
          ),
          numericInput(
            "sample_size",
            "Trees to Sample (N)",
            value = 100,
            min = 10,
            max = 1000,
            step = 10
          ),
          div(
            style = "display: flex; gap: 8px; margin-top: 15px;",
            actionButton(
              "validate_mcmc",
              "Validate & Run Ensemble",
              class = "btn-primary",
              style = "flex: 1;"
            ),
            actionButton(
              "reset_mcmc",
              "Reset",
              class = "btn-secondary",
              style = "flex: 0 0 auto;"
            )
          )
        )
      ),

      mainPanel(
        conditionalPanel(
          condition = "input.pow_mode == 'mode_a'",
          uiOutput("tree_summary") %>%
            withSpinner(color = "#2C3E50"),
          hr(),
          DTOutput("tree_data_preview"),
          hr(),
          plotlyOutput("tdrp_plot", height = "480px") %>%
            withSpinner(color = "#2C3E50"),
          div(
            style = "margin-top: 10px; margin-bottom: 20px; display: flex; gap: 10px;",
            downloadButton(
              "download_plot_pdf",
              "Export Vector PDF",
              class = "btn btn-outline-danger btn-sm"
            )
          ),
          hr(),
          div(
            h6(
              strong("Model Comparison (AIC/BIC)"),
              style = "margin-bottom: 8px;"
            ),
            uiOutput("model_comparison_ui"),
            tags$div(
              style = "font-size: 0.9em; color: #555555; line-height: 1.4; margin-top: 10px;",
              tags$strong("Note: "),
              HTML("AIC = Akaike Information Criterion; BIC = Bayesian Information Criterion; &Delta;AIC / &Delta;BIC = Differences relative to the best model; RSS = Residual Sum of Squares; &mu;<sub>0</sub> = Initial Rate Parameter; &mu;<sub>&infin;</sub> = Long-term Substitution Rate; <i>k</i>= Decay Coefficient.")
            )
          ),
          hr(),
          div(
            h6(strong("Robustness Diagnostic Results"), style = "margin-bottom: 8px;"),
            uiOutput("robustness_ui") %>% withSpinner(color = "#2C3E50")
          ),
          hr(),
          uiOutput("tree_calibration_ui")
        ),

        conditionalPanel(
          condition = "input.pow_mode == 'mode_b'",
          uiOutput("mcmc_summary_ui") %>%
            withSpinner(color = "#2C3E50"),
          hr(),
          h4("PoW-Transformed Ensemble (showing Tree 1)"),
          plotOutput("tree_plot_multi", height = "500px") %>%
            withSpinner(color = "#2C3E50"),
          br(),
          uiOutput("multi_tree_success_ui")
        )
      )
    )
  ),


  nav_panel(
    title = "2. Structural Validation",
    sidebarLayout(
      sidebarPanel(
        h5("1. Tree Input Source"),
        fileInput(
          "custom_pow_tree",
          "Upload Annotated MCC Tree (.tree/.tre/.trees/.nexus/.nwk)",
          accept = c(".tree", ".tre", ".trees", ".nexus", ".nwk")
        ),
        hr(),
        h5("2. Structural Matrix"),
        fileInput(
          "rmsd_matrix",
          "Upload RMSD Matrix (.csv)",
          accept = ".csv"
        ),
        textInput(
          "virus_label",
          "Virus Name for Plot Title:",
          value = "DeepPoW"
        ),
        numericInput(
          "permutations",
          "Mantel Permutations:",
          value = 999,
          min = 99,
          step = 100
        ),
        hr(),
        div(
          style = "display: flex; gap: 8px; margin-top: 10px;",
          actionButton(
            "run_mantel",
            "Run Mantel Test",
            class = "btn-primary",
            style = "flex: 1;"
          ),
          actionButton(
            "reset_mantel",
            "Reset",
            class = "btn-secondary",
            style = "flex: 0 0 auto;"
          )
        )
      ),

      mainPanel(
        plotOutput("struct_vs_phylo_plot", height = "450px") %>%
          withSpinner(color = "#2C3E50"),
        div(
          style = "margin-top: 10px; margin-bottom: 20px;",
          downloadButton(
            "download_mantel_pdf",
            "Export Vector PDF",
            class = "btn btn-outline-secondary btn-sm"
          )
        ),
        hr(),
        uiOutput("mantel_result"),
        hr(),
        h5(strong("Taxon Matching Diagnostic Table (Tree vs. RMSD Matrix)")),
        p(
          "Review this table to ensure all taxon IDs in your tree successfully match the row/column names of your uploaded RMSD matrix.",
          style = "font-size:0.85em; color:gray;"
        ),
        DTOutput("taxon_diag_table") %>%
          withSpinner(color = "#2C3E50")
      )
    )
  ),


  nav_panel(
    title = "3. Advanced Utilities",
    navset_card_underline(
      
      nav_panel(
        "TreeTime TDRP Explorer",
        mod_treetime_ui("treetime_mod")
      ),

      nav_panel(
        "RMSD Matrix Generator",
        sidebarLayout(
          sidebarPanel(
            h5("Local Execution Config"),
            selectInput(
              "rmsd_tool",
              "Alignment Engine:",
              choices = c("Foldseek (Fast)", "TM-align (Classic)")
            ),
            
            div(
              style = "display: flex; gap: 8px; margin-bottom: 10px;",
              div(
                style = "flex: 1;",
                textInput(
                  "custom_exe",
                  "Executable Path (Optional):",
                  value = "",
                  placeholder = "Leave blank to auto-detect",
                  width = "100%"
                )
              ),
              div(
                style = "margin-top: 32px;",
                shinyFilesButton(
                  "custom_exe_btn",
                  "Browse",
                  "Select Executable",
                  multiple = FALSE,
                  class = "btn btn-secondary"
                )
              )
            ),
            
            div(
              style = "display: flex; gap: 8px; margin-bottom: 15px;",
              div(
                style = "flex: 1;",
                textInput(
                  "pdb_dir",
                  "PDB/CIF Directory Path:",
                  value = "./alphafold_pdbs",
                  width = "100%"
                )
              ),
              div(
                style = "margin-top: 32px;",
                shinyDirButton(
                  "pdb_dir_btn",
                  "Browse",
                  "Select PDB/CIF Directory",
                  class = "btn btn-secondary"
                )
              )
            ),
            
            actionButton(
              "generate_rmsd",
              "Generate RMSD Matrix",
              class = "btn-primary",
              style = "width: 100%; padding: 10px; font-weight: bold;"
            )
          ),
          mainPanel(
            h5("RMSD Matrix Results"),
            p(
              "Structural alignment results and execution trace.",
              style = "font-size: 0.85em; color: gray;"
            ),
            navset_card_tab(
              nav_panel(
                "Status", 
                uiOutput("rmsd_status_ui") %>% withSpinner(color = "#2C3E50")
              ),
              nav_panel(
                "Heatmap", 
                plotlyOutput("rmsd_heatmap", height = "500px") %>% withSpinner(color = "#2C3E50"),
                div(style = "margin-top: 15px;", downloadButton("download_heatmap_pdf", "Export Heatmap PDF", class = "btn btn-outline-secondary btn-sm"))
              ),
              nav_panel(
                "Data Table", 
                DTOutput("rmsd_table") %>% withSpinner(color = "#2C3E50"),
                div(style = "margin-top: 15px;", downloadButton("download_rmsd_csv", "Download CSV", class = "btn btn-outline-secondary btn-sm"))
              )
            )
          )
        )
      ),
      
      nav_panel(
        "BEAST 2 XML Patcher",
        sidebarLayout(
          sidebarPanel(
            h5("State-Resume XML Auto-Patcher"),
            p(
              "Injects branch-rate loggers into your BEAST 2 XML so you can resume your state and capture rate data for DeepPoW.",
              style = "font-size: 0.85em; color: gray;"
            ),
            fileInput(
              "beast2_xml",
              "Upload Original XML (.xml)",
              accept = ".xml"
            ),
            actionButton(
              "patch_xml",
              "Patch XML",
              class = "btn-primary",
              style = "width: 100%;"
            ),
            hr(),
            downloadButton(
              "download_patched_xml",
              "Download Patched XML",
              class = "btn-primary",
              style = "width: 100%;"
            )
          ),
          mainPanel(
            h5("Modification Log"),
            uiOutput("xml_patch_log") %>%
              withSpinner(color = "#2C3E50"),
            hr(),
            h5("How to Resume"),
            uiOutput("resume_command_hint")
          )
        )
      )
    )
  ),


  nav_panel(
    title = "4. About",
    mod_about_ui("about_tab")
  )
)

# ==========================================
# 3. Server
# ==========================================

server <- function(input, output, session) {


  volumes <- tryCatch({
    c(WorkingDir = getwd(), shinyFiles::getVolumes()())
  }, error = function(e) {
    fallback <- c(WorkingDir = getwd(), Home = Sys.getenv("HOME"))
    if (.Platform$OS.type == "windows") {
      fallback <- c(fallback, "C:/" = "C:/", "D:/" = "D:/")
    } else {
      fallback <- c(fallback, "Root" = "/")
    }
    fallback
  })
  
  shinyFileChoose(
    input, 
    "custom_exe_btn", 
    roots = volumes, 
    session = session
  )
  
  observeEvent(input$custom_exe_btn, {
    if (!is.integer(input$custom_exe_btn)) {
      file_info <- parseFilePaths(volumes, input$custom_exe_btn)
      if (nrow(file_info) > 0) {
        updateTextInput(session, "custom_exe", value = as.character(file_info$datapath))
      }
    }
  })

  shinyDirChoose(
    input, 
    "pdb_dir_btn", 
    roots = volumes, 
    session = session
  )
  
  observeEvent(input$pdb_dir_btn, {
    if (!is.integer(input$pdb_dir_btn)) {
      dir_path <- parseDirPath(volumes, input$pdb_dir_btn)
      updateTextInput(session, "pdb_dir", value = dir_path)
    }
  })

  # --- State Init ---
  app_state <- reactiveValues(
    tree_obj = NULL,
    raw_df = NULL,
    clean_df = NULL,
    pow_model = NULL,
    current_plot = NULL,
    pow_trees_multi = NULL,
    mcmc_valid = FALSE,
    mcmc_log_text = "Please upload the required BEAST files and click 'Validate & Run Ensemble'.",
    mcmc_results = NULL, 
    patched_xml_doc = NULL,
    patched_xml_log = "Waiting for XML upload...",
    rmsd_status_ui = "<span style='color:gray;'>Waiting for execution...</span>",
    rmsd_matrix_data = NULL,
    robustness_shift = NULL  
  )

  get_mu_max_value <- reactive({
    switch(
      input$virus_type,
      "RNA" = 3.65e-2,
      "ssDNA" = 2.0e-2,
      "dsDNA" = 3.0e-3,
      "Other" = input$mu_max_b
    )
  })

  clear_mode_a_state <- function() {
    app_state$pow_model <- NULL
    app_state$current_plot <- NULL
    app_state$robustness_shift <- NULL

    output$model_comparison_ui <- renderUI({ NULL })
    output$robustness_ui <- renderUI({ NULL }) 
    output$tree_calibration_ui <- renderUI({ NULL })
    output$tdrp_plot <- renderPlotly({ NULL })
  }

  clear_mode_b_state <- function() {
    app_state$pow_trees_multi <- NULL

    output$multi_tree_success_ui <- renderUI({
      NULL
    })
  }

  get_active_tree <- reactive({
    req(input$custom_pow_tree)
    load_pow_tree_robust(input$custom_pow_tree$datapath)
  })

  output$tree_summary <- renderUI({
    p("Please upload a tree file and click 'Parse Tree'.", style = "color: gray; font-style: italic;")
  })

  observeEvent(input$reset_tree, {

    app_state$tree_obj <- NULL
    app_state$raw_df <- NULL
    app_state$clean_df <- NULL

    clear_mode_a_state()

    updateSelectInput(
      session,
      "rate_col",
      choices = c(
        "Auto-detect",
        "rate",
        "rate_median",
        "clock.rate",
        "branch_rate"
      ),
      selected = "Auto-detect"
    )

    output$tree_summary <- renderUI({
      p("Workspace reset. Please upload a new single tree file.", style = "color: gray; font-style: italic;")
    })

    output$tree_data_preview <- renderDT({
      NULL
    })
  })

  observeEvent(input$reset_pow_a, {
    clear_mode_a_state()
    showNotification(
      "Mode A model fit reset.",
      type = "message"
    )
  })

  observeEvent(input$process_tree, {

    req(input$beast_tree)

    tryCatch({

      parsed <- suppressWarnings(
        parse_beast_tree(
          input$beast_tree$datapath,
          origin = "beast"
        )
      )

      updateSelectInput(
        session,
        "rate_col",
        choices = c(
          "Auto-detect",
          colnames(parsed$data)
        ),
        selected = input$rate_col
      )

      rate_col_arg <- if (
        input$rate_col == "Auto-detect"
      ) {
        NULL
      } else {
        input$rate_col
      }

      raw_df <- extract_node_features(
        parsed$data,
        rate_col_arg
      )

      app_state$tree_obj <- parsed$tree
      app_state$raw_df <- raw_df
      app_state$clean_df <- filter_extremes(
        raw_df,
        input$min_time,
        input$max_rate_q
      )

      output$tree_summary <- renderUI({
        div(
          class = "alert alert-success",
          style = "padding: 10px; margin-bottom: 10px;",
          HTML(paste0(
            "<b>Successfully parsed:</b> ", input$beast_tree$name, "<br>",
            "<b>Total nodes extracted:</b> ", nrow(app_state$raw_df), "<br>",
            "<b>Nodes retained after filtering:</b> ", nrow(app_state$clean_df)
          ))
        )
      })

      output$tree_data_preview <- renderDT({

        preview_df <- as.data.frame(
          app_state$clean_df %>% head(50),
          stringsAsFactors = FALSE
        )

        datatable(
          preview_df,
          options = list(
            pageLength = 5L,
            scrollX = TRUE
          ),
          rownames = FALSE
        )
      })

    }, error = function(e) {

      showNotification(
        paste("Error processing tree:", e$message),
        type = "error"
      )
    })
  })
  
  observeEvent(input$run_robustness, {
    req(app_state$clean_df)
    
    tryCatch({
      trim_val <- input$trim_quantile_ui / 100
      robust_res <- robustness_pow_model(app_state$clean_df, trim_quantile = trim_val)
      
      if (is.null(robust_res) || is.null(robust_res$comparison_table)) {
        stop("Robustness check failed. The trimmed dataset might be too small or failed to converge.")
      }
      
      app_state$robustness_shift <- robust_res$lambda_shift
      
      output$robustness_ui <- renderUI({
        tbl <- robust_res$comparison_table
        tbl[] <- lapply(tbl, function(col) {
          if (is.character(col) || is.factor(col)) {
            col <- as.character(col)
            col <- gsub("mu_0", "\u03bc\u2080", col)
            col <- gsub("mu_infty", "\u03bc\u221e", col)
            col <- gsub("lambda", "k", col)
          }
          col
        })
        
        shift_val <- robust_res$lambda_shift
        
        full_mu_0 <- as.numeric(tbl$Full_Data[1])
        full_mu_infty <- as.numeric(tbl$Full_Data[2])
        full_lambda <- as.numeric(tbl$Full_Data[3])
        
        flat_reason <- ""
        if (is.na(full_lambda) || full_lambda <= 1.01e-4) {
          flat_reason <- "k is at the lower fitting bound (10⁻⁴; no detectable decay)"
        } else if (is.na(full_mu_0) || is.na(full_mu_infty) || full_mu_0 <= full_mu_infty) {
          flat_reason <- "\u03bc\u2080 \u2264 \u03bc\u221e (inverse or flat trajectory)"
        }
        is_paradoxical_or_flat <- (nchar(flat_reason) > 0)
        
        if (is.na(shift_val)) {
          alert_div <- div(
            class = "alert alert-danger",
            style = "padding: 10px; font-size: 0.85em; margin-top: 10px;",
            HTML(sprintf("<b>Convergence Failed:</b> The model could not converge after trimming %d%% outliers. The TDRP signal might be entirely driven by extreme noise.", input$trim_quantile_ui))
          )
        } else if (is_paradoxical_or_flat && shift_val < 0.20) {
          alert_div <- div(
            class = "alert alert-secondary",
            style = "padding: 10px; font-size: 0.85em; margin-top: 10px;",
            HTML(sprintf("<b>Stable but Flat/Inverse:</b> The decay coefficient (k) shifted by %.1f%%. The parameter space is stable, but <b>%s</b>, indicating a constant rate rather than a genuine TDRP decay.", shift_val * 100, flat_reason))
          )
        } else if (shift_val < 0.20) {
          alert_div <- div(
            class = "alert alert-success",
            style = "padding: 10px; font-size: 0.85em; margin-top: 10px;",
            HTML(sprintf("<b>Globally Robust:</b> The decay coefficient (k) shifted by %.1f%% after removing outliers. The TDRP signal is highly stable and shows limited sensitivity to the trimming procedure.", shift_val * 100))
          )
        } else if (shift_val < 1.00) {
          alert_div <- div(
            class = "alert alert-warning",
            style = "padding: 10px; font-size: 0.85em; margin-top: 10px;",
            HTML(sprintf("<b>Moderately Sensitive:</b> The decay coefficient (k) shifted by %.1f%%. Please manually review the terminal branches for potential sequencing errors or severe bottlenecks.", shift_val * 100))
          )
        } else {
          alert_div <- div(
            class = "alert alert-danger",
            style = "padding: 10px; font-size: 0.85em; margin-top: 10px;",
            HTML(sprintf("<b>Severe Parameter Instability Detected:</b> The decay coefficient (k) shifted by %.1f%% after outlier trimming. The fitted decay pattern is highly sensitive to extreme branch-rate observations and requires further investigation.", shift_val * 100))
          )
        }
        
        tagList(
          tags$table(
            class = "table table-sm table-bordered",
            style = "background-color: #ffffff; font-size: 0.9em;",
            tags$thead(
              tags$tr(lapply(colnames(tbl), function(col) tags$th(as.character(col))))
            ),
            tags$tbody(
              lapply(seq_len(nrow(tbl)), function(i) {
                tags$tr(lapply(tbl[i, , drop = FALSE], function(val) tags$td(as.character(val))))
              })
            )
          ),
          alert_div
        )
      })
    }, error = function(e) {
      showNotification(paste("Robustness Check Error:", e$message), type = "error")
    })
  })

  observeEvent(input$run_pow, {

    req(app_state$clean_df)
    
    output$robustness_ui <- renderUI({ NULL })
    app_state$robustness_shift <- NULL

    df <- app_state$clean_df

    tryCatch({

      app_state$pow_model <- fit_pow_model(df)

      if (is.null(app_state$pow_model)) {
        stop("Model fitting failed to converge.")
      }

      compare_tbl <- compare_models(
        df,
        app_state$pow_model
      )

      coefs_raw <- coef(app_state$pow_model)

      mu_0_val <- json_safe_scalar_numeric(
        coefs_raw["mu_0"]
      )

      mu_infty_val <- json_safe_scalar_numeric(
        coefs_raw["mu_infty"]
      )

      lambda_val <- json_safe_scalar_numeric(
        coefs_raw["lambda"]
      )

      output$model_comparison_ui <- renderUI({

        pow_delta_aic <- as.numeric(
          compare_tbl$Delta_AIC[
            compare_tbl$Model == "PoW Model (TDRP)"
          ]
        )

        pow_delta_bic <- as.numeric(
          compare_tbl$Delta_BIC[
            compare_tbl$Model == "PoW Model (TDRP)"
          ]
        )

        is_conflict <- (
          !is.na(pow_delta_aic) &&
          !is.na(pow_delta_bic) &&
          pow_delta_aic <= 0 &&
          pow_delta_bic > 0
        )

        tagList(

          tags$table(
            class = "table table-sm table-bordered",
            style = "background-color: #ffffff;",

            tags$thead(
              tags$tr(
                lapply(
                  colnames(compare_tbl),
                  function(col) {
                    tags$th(
                      as.character(col)
                    )
                  }
                )
              )
            ),

            tags$tbody(
              lapply(
                seq_len(nrow(compare_tbl)),
                function(i) {

                  tags$tr(
                    lapply(
                      compare_tbl[i, , drop = FALSE],
                      function(val) {
                        tags$td(
                          as.character(val)
                        )
                      }
                    )
                  )
                }
              )
            )
          ),

          if (is_conflict) {

            div(
              class = "alert alert-info",
              style = "margin-top: 10px; font-size: 0.88em; padding: 10px;",
              HTML(
                "<b> Model Divergence Note (&Delta;AIC vs &Delta;BIC):</b><br>
                 AIC favors the PoW model, whereas BIC favors the Constant Rate model.
                 This may reflect the stronger complexity penalty imposed by BIC.
                 <b>DeepPoW strictly prioritizes BIC.</b>"
              )
            )
          }
        )
      })

      pow_delta_bic <- as.numeric(
        compare_tbl$Delta_BIC[
          compare_tbl$Model == "PoW Model (TDRP)"
        ]
      )

      output$tree_calibration_ui <- renderUI({
        shift_val <- app_state$robustness_shift
        if (!is.na(lambda_val) && lambda_val > 150&& !isTRUE(input$expert_mode)) {
          div(
            class = "card border-danger mb-3",
            div(class = "card-header bg-danger text-white", strong("Ultra-Rapid Decay Detected (k > 150)")),
            div(class = "card-body",
              p(HTML(paste0("<b>Scale Warning:</b> The estimated decay coefficient is extremely high (k = ", signif(lambda_val, 3), "), implying a characteristic rate-decay timescale of only a few days or less, depending on the value of k. This may indicate a fitting artifact (e.g., overfitting to identical sequences from ultra-dense contact tracing) in standard inter-host epidemiological trees.")), style = "font-size: 0.9em; color: #721C24;"),
              p(HTML("<b>&#8594; Calibration is not recommended.</b>"), style = "font-size: 0.85em; font-weight: bold; color: #721C24; margin-bottom: 0;")
            )
          )
        } else if (!is.null(shift_val) && (is.na(shift_val) || shift_val >= 1.0)) {
          div(
            class = "card border-danger mb-3",
            div(class = "card-header bg-danger text-white", strong("Methodological Calibration Vetoed (Potential Fitting Artifact)")),
            div(class = "card-body",
                p(HTML("<b>Robustness Failure:</b> The diagnostic revealed that the TDRP signal is mathematically unstable, driven by outliers, or failed to converge upon trimming. Standard calibration is strongly discouraged."), style = "font-size: 0.9em; color: #721C24;"),
                p(HTML("<b>&#8594; Methodological Calibration Veto:</b> Mode B calibration is strongly discouraged based on the current fit. Exploratory analyses remain available at the user's discretion."), style = "font-size: 0.85em; font-weight: bold; color: #721C24; margin-bottom: 0;")
            )
          )
        } else if (is.na(lambda_val) || lambda_val <= 1.01e-4 || is.na(mu_0_val) || is.na(mu_infty_val) || mu_0_val <= mu_infty_val) {
          div(
            class = "card border-warning mb-3",
            div(class = "card-header bg-warning text-dark", strong("Signal Muffled on Summary Tree (Flat/Inverse Trajectory)")),
            div(class = "card-body",
              p(HTML("<b>Diagnostic Note:</b> The phenomenological fit on this specific MCC tree did not detect a clear decay pattern (\u03bc\u2080 &le; \u03bc\u221e or k at the lower fitting bound). For deep-time evolutionary datasets (e.g., ancient viral origins), standard relaxed-clock summary trees often smooth out genuine time-dependent signals. If your dataset belongs to this category, the mechanistic TDRP signal may still be recovered across the full posterior."), style = "font-size: 0.9em; color: #856404;"),
              p(HTML("<b>&#8594; You may cautiously proceed to Mode B</b> for rigorous MCMC ensemble calibration."), style = "font-size: 0.85em; font-weight: bold; color: #856404; margin-bottom: 0;")
            )
          )
        } else if (is.na(pow_delta_bic) || pow_delta_bic > 0) {
          div(
            class = "card border-warning mb-3",
            div(class = "card-header bg-warning text-dark", strong("Constant Rate Favored (No TDRP Signal Supported)")),
            div(class = "card-body",
              p(HTML("<b>Statistical Result:</b> The Constant Rate model fits the data better (&Delta;BIC > 0). The data on this summary tree does not statistically support time-dependent rate calibration."), style = "font-size: 0.9em; color: #856404;"),
              p(HTML("<b>&#8594; Calibration Not Justified:</b> Using your original BEAST trees is recommended. Mode B calibration is not required unless you are explicitly exploring hidden signals in the full MCMC posterior."), style = "font-size: 0.85em; font-weight: bold; color: #856404; margin-bottom: 0;")
            )
          )
        } else if (!is.null(shift_val) && !is.na(shift_val) && shift_val > 0.2) {
          div(
            class = "card border-warning mb-3",
            div(class = "card-header bg-warning text-dark", strong("Calibration Permitted (Proceed with Caution)")),
            div(class = "card-body",
                p(HTML("<b>Proceed with Caution:</b> The model is moderately sensitive to outliers. Please review your tree branches."), style = "font-size: 0.9em; color: #856404;"),
                p(HTML("<b>&#8594; Action Required:</b> Review your tree branches. If proceeding to Mode B, interpret results with caution."), style = "font-size: 0.85em; font-weight: bold; color: #856404; margin-bottom: 0;")
            )
          )
        } else {
          title_text <- if (!is.null(shift_val) && !is.na(shift_val)) "Tree Calibration (Robustness Verified)" else "Tree Calibration (PoW Supported)"
          div(
            class = "card border-success mb-3",
            div(class = "card-header bg-success text-white", strong(title_text)),
            div(class = "card-body",
              p(HTML("A robust TDRP pattern consistent with time-dependent rate decay was detected."), style = "font-size: 0.9em; color: #155724;"),
              p(HTML("<b>&#8594; Please proceed to Mode B</b> to perform rigorous MCMC ensemble calibration. Single tree extraction is disabled for MCC trees to preserve HPD intervals."), style = "font-size: 0.85em; font-weight: bold; color: #155724; margin-bottom: 0;")
            )
          )
        }
      })

      t_min <- max(
        min(df$branch.length, na.rm = TRUE),
        1e-5
      )

      t_max <- max(
        df$branch.length,
        na.rm = TRUE
      )

      t_seq <- exp(
        seq(
          log(t_min),
          log(t_max),
          length.out = 100
        )
      )

      curve_df <- predict_pow_with_ci(
        app_state$pow_model,
        t_seq
      )

      plot_df <- as.data.frame(
        df,
        stringsAsFactors = FALSE
      )

      plot_df$log_time <- log10(
        pmax(
          as.numeric(plot_df$branch.length),
          1e-12
        )
      )

      plot_df$log_rate <- log10(
        pmax(
          as.numeric(plot_df$rate),
          1e-12
        )
      )

      curve_df <- as.data.frame(
        curve_df,
        stringsAsFactors = FALSE
      )

      curve_df$log_time <- log10(
        pmax(
          as.numeric(curve_df$branch.length),
          1e-12
        )
      )

      curve_df$log_rate_pred <- log10(
        pmax(
          as.numeric(curve_df$rate_pred),
          1e-12
        )
      )

      curve_df$log_upper_ci <- log10(
        pmax(
          as.numeric(curve_df$upper_ci),
          1e-12
        )
      )

      curve_df$log_lower_ci <- log10(
        pmax(
          as.numeric(curve_df$lower_ci),
          1e-12
        )
      )

      output$tdrp_plot <- renderPlotly({

        req(app_state$pow_model)

        plot_df_json <- as.data.frame(
          plot_df,
          stringsAsFactors = FALSE
        )

        curve_df_json <- as.data.frame(
          curve_df,
          stringsAsFactors = FALSE
        )

        plot_df_json$log_time <- json_safe_numeric(
          plot_df_json$log_time
        )

        plot_df_json$log_rate <- json_safe_numeric(
          plot_df_json$log_rate
        )

        curve_df_json$log_time <- json_safe_numeric(
          curve_df_json$log_time
        )

        curve_df_json$log_upper_ci <- json_safe_numeric(
          curve_df_json$log_upper_ci
        )

        curve_df_json$log_lower_ci <- json_safe_numeric(
          curve_df_json$log_lower_ci
        )

        curve_df_json$log_rate_pred <- json_safe_numeric(
          curve_df_json$log_rate_pred
        )

        node_values <- if (
          "node" %in% names(plot_df_json)
        ) {
          json_safe_character(plot_df_json$node)
        } else {
          rep("", nrow(plot_df_json))
        }

        hover_text <- vapply(
          seq_len(nrow(plot_df_json)),
          function(i) {

            paste0(
              "Node: ",
              node_values[i],
              "<br>Log10(Time): ",
              signif(plot_df_json$log_time[i], 4),
              "<br>Log10(Rate): ",
              signif(plot_df_json$log_rate[i], 4)
            )
          },
          FUN.VALUE = character(1)
        )

        hover_text <- json_safe_character(
          hover_text
        )

        param_annotation_text <- paste0(
          "<b>DeepPoW Parameters:</b><br>",
          "\u03bc<sub>0</sub> = ", format(mu_0_val, scientific = TRUE, digits = 3), "<br>",
          "\u03bc<sub>\u221e</sub> = ", format(mu_infty_val, scientific = TRUE, digits = 3), "<br>",
          "k = ", signif(lambda_val, 3)
        )

        p <- plotly::plot_ly()

        p <- p %>%
          plotly::add_trace(
            x = json_safe_numeric(
              curve_df_json$log_time
            ),
            y = json_safe_numeric(
              curve_df_json$log_upper_ci
            ),
            type = "scatter",
            mode = "lines",
            line = list(
              color = "rgba(0,0,0,0)"
            ),
            showlegend = FALSE,
            hoverinfo = "none"
          )

        p <- p %>%
          plotly::add_trace(
            x = json_safe_numeric(
              curve_df_json$log_time
            ),
            y = json_safe_numeric(
              curve_df_json$log_lower_ci
            ),
            type = "scatter",
            mode = "lines",
            fill = "tonexty",
            fillcolor = "rgba(217,83,79,0.15)",
            line = list(
              color = "rgba(0,0,0,0)"
            ),
            name = "95% CI",
            hoverinfo = "none"
          )

        p <- p %>%
          plotly::add_trace(
            x = json_safe_numeric(
              plot_df_json$log_time
            ),
            y = json_safe_numeric(
              plot_df_json$log_rate
            ),
            type = "scatter",
            mode = "markers",
            text = hover_text,
            hoverinfo = "text",
            marker = list(
              size = 6,
              color = "rgba(44,62,80,0.6)",
              line = list(
                color = "rgba(44,62,80,0.9)",
                width = 0.5
              )
            ),
            name = "Branch Rates"
          )

        p <- p %>%
          plotly::add_trace(
            x = json_safe_numeric(
              curve_df_json$log_time
            ),
            y = json_safe_numeric(
              curve_df_json$log_rate_pred
            ),
            type = "scatter",
            mode = "lines",
            line = list(
              color = "#D9534F",
              width = 2.5
            ),
            name = "PoW Fit",
            hoverinfo = "none"
          )

        p <- p %>%
          plotly::layout(

            title = list(
              text = "<b>Time-Dependent Rate Phenomenon (TDRP)</b>",
              font = list(
                size = 12
              )
            ),

            xaxis = list(
              title = "<b>Log<sub>10</sub>(Time-scale, years)</b>"
            ),

            yaxis = list(
              title = "<b>Log<sub>10</sub>(Evolutionary Rate, μ)</b>"
            ),

            annotations = list(
              list(
                x = 0.05,
                y = 0.95,
                xref = "paper",
                yref = "paper",
                text = as.character(
                  param_annotation_text
                ),
                showarrow = FALSE,
                bgcolor = "rgba(255,255,255,0.85)",
                bordercolor = "#CCCCCC",
                borderwidth = 1
              )
            ),

            plot_bgcolor = "#FFFFFF",
            paper_bgcolor = "#FFFFFF",

            showlegend = TRUE,

            legend = list(
              x = 0.05,
              y = 0.12,
              xanchor = "left",
              yanchor = "bottom",
              bgcolor = "rgba(255,255,255,0.85)",
              bordercolor = "#CCCCCC",
              borderwidth = 1
            )
          )

        p <- plotly::config(
          p,
          displaylogo = FALSE
        )

        app_state$current_plot <- p

        p
      })

    }, error = function(e) {

      showNotification(
        paste("PoW Fit Error:", e$message),
        type = "error"
      )
    })
  })

  output$download_plot_pdf <- downloadHandler(

    filename = function() {
      paste0(
        "TDRP_Vector_Plot_",
        Sys.Date(),
        ".pdf"
      )
    },

    content = function(file) {

      req(
        app_state$clean_df,
        app_state$pow_model
      )

      tryCatch({

        export_pow_pdf(
          file,
          app_state$clean_df,
          app_state$pow_model
        )

      }, error = function(e) {

        showNotification(
          paste(
            "PDF Export Failed:",
            e$message
          ),
          type = "error"
        )
      })
    }
  )

  output$mcmc_summary_ui <- renderUI({
    if (!app_state$mcmc_valid) {
      if (grepl("Waiting|Please upload", app_state$mcmc_log_text)) {
        return(p(app_state$mcmc_log_text, style = "color: gray; font-style: italic;"))
      } else {
        return(div(class = "alert alert-danger", HTML(app_state$mcmc_log_text)))
      }
    } else {
      res <- app_state$mcmc_results
      tagList(
        div(
          class = "alert alert-success",
          HTML("<b>Validation & Ensemble Transformation Successful!</b><br>
               <span style='font-size:0.9em; color:#155724;'>Calibrated ensemble trees are ready for download or Structural Validation (Tab 2).</span>")
        ),
        tags$table(class = "table table-sm table-bordered", style = "background-color: #ffffff;",
          tags$thead(tags$tr(tags$th("Parameter"), tags$th("Configuration / Result"))),
          tags$tbody(
            tags$tr(tags$td(strong("Substitution Model")), tags$td(res$sub_model)),
            tags$tr(tags$td(strong("Rate Log")), tags$td(paste0("Matched column '", res$rate_col, "' (", res$rate_valid, " post burn-in)"))),
            tags$tr(tags$td(strong("Trees Processed")), tags$td(paste0(res$target_n, " sampled trees transformed"))),
            tags$tr(tags$td(strong("Active muMax Boundary")), tags$td(paste0(res$active_mumax, " (Type: ", res$virus_type, ")")))
          )
        )
      )
    }
  })

  observeEvent(input$validate_mcmc, {

    req(
      input$rate_log,
      input$trees_file
    )

    if (input$sub_model == "HKY") {
      req(input$hky_log)
    }

    tryCatch({

      rate_df <- read_beast_log(
        input$rate_log$datapath
      )

      rate_col_idx <- grep(
        "meanRate|rate\\.mean|clock\\.rate|clockRate|rate_median",
        colnames(rate_df),
        ignore.case = TRUE
      )[1]

      if (is.na(rate_col_idx)) {
        stop(
          "[Rate Log]: Could not find 'meanRate' column."
        )
      }

      rate_avail <- nrow(rate_df)

      rate_valid <- rate_avail -
        max(
          1,
          floor(rate_avail * 0.1)
        )

      if (input$sub_model == "HKY") {

        hky_df <- read_beast_log(
          input$hky_log$datapath
        )

        kappa_idx <- grep(
          "^kappa$|kappa",
          colnames(hky_df),
          ignore.case = TRUE
        )[1]

        freq_indices <- find_freq_cols(
          colnames(hky_df)
        )

        if (is.na(kappa_idx)) {
          stop(
            "[HKY Log]: Missing 'kappa' column."
          )
        }

        if (any(is.na(freq_indices))) {
          stop(
            "[HKY Log]: Could not locate all 4 nucleotide frequency columns."
          )
        }

        hky_avail <- nrow(hky_df)

        hky_valid <- hky_avail -
          max(
            1,
            floor(hky_avail * 0.1)
          )

      } else {

        hky_valid <- rate_valid
      }

      tree_lines <- readLines(
        input$trees_file$datapath,
        warn = FALSE
      )

      tree_start <- grep(
        "^tree",
        tree_lines,
        ignore.case = TRUE
      )[1]

      if (is.na(tree_start)) {
        stop(
          "[Trees File]: Cannot find lines starting with 'tree'."
        )
      }

      tree_avail <- length(tree_lines) -
        tree_start + 1

      tree_valid <- tree_avail -
        max(
          1,
          floor(tree_avail * 0.1)
        )

      target_n <- input$sample_size

      min_valid <- min(
        rate_valid,
        hky_valid,
        tree_valid
      )

      if (min_valid < target_n) {
        stop(
          sprintf(
            "Insufficient data post burn-in. Requested: %d, Available: %d",
            target_n,
            min_valid
          )
        )
      }

      active_mumax <- get_mu_max_value()
      current_mumax <- get_mu_max_value()

      hky_path <- if (
        input$sub_model == "HKY"
      ) {
        input$hky_log$datapath
      } else {
        NULL
      }

      withProgress(
        message = "Executing MCMC Ensemble PoW transformation...",
        value = 0.5,
        {

          trees <- pow_transform(
            rate_log_path = input$rate_log$datapath,
            trees_path = input$trees_file$datapath,
            hky_log_path = hky_path,
            muMax_val = current_mumax,
            sample_size = input$sample_size,
            sub_model = input$sub_model
          )
        }
      )

      class(trees) <- "multiPhylo"

      app_state$pow_trees_multi <- trees

      app_state$mcmc_results <- list(
        sub_model = input$sub_model,
        rate_col = colnames(rate_df)[rate_col_idx],
        rate_valid = rate_valid,
        target_n = target_n,
        active_mumax = active_mumax,
        virus_type = input$virus_type
      )
      
      app_state$mcmc_valid <- TRUE

      output$multi_tree_success_ui <- renderUI({

        div(

          div(
            class = "alert alert-success",
            "Ensemble Calibration Complete! You can download the transformed multi-tree ensemble below."
          ),

          div(
            style = "display: flex; gap: 10px; margin-top: 10px;",

            downloadButton(
              "download_trees_multi",
              "Download Ensemble (.trees)",
              class = "btn btn-primary"
            ),

            downloadButton(
              "download_plot_multi",
              "Download Plot (PDF)",
              class = "btn btn-outline-secondary"
            )
          )
        )
      })

    }, error = function(e) {

      app_state$mcmc_valid <- FALSE
      app_state$mcmc_log_text <- paste("Validation/Transformation Failed:<br>", e$message)
    })
  })

  observeEvent(input$reset_mcmc, {

    app_state$mcmc_valid <- FALSE
    app_state$mcmc_log_text <- "Please upload the required BEAST files and click 'Validate & Run Ensemble'."

    clear_mode_b_state()

  })

  output$tree_plot_multi <- renderPlot({

    req(
      app_state$pow_trees_multi
    )

    plot(
      app_state$pow_trees_multi[[1]],
      main = "PoW-Transformed Tree (Sample #1 from MCMC Ensemble)",
      cex = 0.8
    )

    add.scale.bar()
  })

  output$download_trees_multi <- downloadHandler(

    filename = function() {
      paste0(
        "PoW_MCMC_Ensemble_",
        Sys.Date(),
        ".trees"
      )
    },

    content = function(file) {

      req(
        app_state$pow_trees_multi
      )

      tryCatch({

        trees_to_write <- app_state$pow_trees_multi

        class(trees_to_write) <- "multiPhylo"

        ape::write.nexus(
          trees_to_write,
          file = file
        )

      }, error = function(e) {

        stop(
          paste(
            "Failed to write trees:",
            e$message
          )
        )
      })
    }
  )

  output$download_plot_multi <- downloadHandler(

    filename = function() {
      paste0(
        "PoW_MultiTree_Plot_",
        Sys.Date(),
        ".pdf"
      )
    },

    content = function(file) {

      req(
        app_state$pow_trees_multi
      )

      pdf(
        file,
        width = 10,
        height = 8
      )

      tryCatch({

        plot(
          app_state$pow_trees_multi[[1]],
          main = "PoW-Transformed Tree",
          cex = 0.8
        )

        add.scale.bar()

      }, finally = {

        dev.off()
      })
    }
  )


  observeEvent(input$reset_mantel, {
    output$struct_vs_phylo_plot <- renderPlot({
      plot.new()
      text(0.5, 0.5, "Click 'Run Mantel Test' to generate congruence plot.", cex = 1.1, col = "#555555")
    })
    output$mantel_result <- renderUI({
      p("Please upload your RMSD matrix and click 'Run Mantel Test'.", style = "color: gray; font-style: italic;")
    })
    output$taxon_diag_table <- renderDT({ NULL })
    updateTextInput(session, "virus_label", value = "Virus name")
    updateNumericInput(session, "permutations", value = 999)
    showNotification("Mantel Test module reset.", type = "message")
  })

  output$taxon_diag_table <- renderDT({

    req(
      input$rmsd_matrix
    )

    tree_obj <- tryCatch(
      {
        get_active_tree()
      },
      error = function(e) {
        NULL
      }
    )

    diag_df <- check_taxon_matching(
      tree_obj,
      input$rmsd_matrix$datapath
    )

    diag_df <- as.data.frame(
      diag_df,
      stringsAsFactors = FALSE
    )

    datatable(
      diag_df,
      options = list(
        pageLength = 8L,
        scrollX = TRUE
      ),
      rownames = FALSE
    )
  })

  mantel_reactive <- eventReactive(
    input$run_mantel,
    {

      req(
        input$rmsd_matrix
      )

      tryCatch({

        active_tree <- get_active_tree()

        if (is.null(active_tree)) {
          stop(
            "Tree object is missing. Please load a valid PoW tree."
          )
        }

        phylo_df <- calculate_phylo_distance(
          active_tree
        )

        matched_df <- match_structure_matrix(
          phylo_df,
          input$rmsd_matrix$datapath
        )

        mantel_res <- run_mantel_test(
          matched_df,
          permutations = input$permutations
        )

        list(
          matched_df = matched_df,
          mantel_res = mantel_res,
          error = NULL
        )

      }, error = function(e) {

        list(
          matched_df = NULL,
          mantel_res = NULL,
          error = e$message
        )
      })
    }
  )


output$mantel_result <- renderUI({

    if (input$run_mantel == 0) {
      return(p("Please upload your RMSD matrix and click 'Run Mantel Test'.", style = "color: gray; font-style: italic;"))
    }

    res <- mantel_reactive()

    if (!is.null(res$error)) {
      return(
        div(
          class = "alert alert-danger",
          style = "padding: 10px;",
          HTML(paste("<b>Mantel Analysis Error:</b><br>", res$error))
        )
      )
    }

    mantel_res <- res$mantel_res
    matched_df <- res$matched_df

    is_sig <- mantel_res$p_value < 0.05
    conclusion_text <- if (is_sig) {
      "A statistically significant correlation was detected between calibrated phylogenetic distances and protein structural distances."
    } else {
      "No statistically significant correlation detected."
    }
    alert_class <- if (is_sig) "alert alert-success" else "alert alert-secondary"

    tagList(
      h6(strong("Mantel Test Results"), style = "margin-bottom: 10px; color: #2C3E50;"),
      tags$table(
        class = "table table-sm table-bordered",
        style = "background-color: #ffffff; margin-bottom: 10px;",
        tags$tbody(
          tags$tr(tags$td(strong("Mantel Statistic (r)")), tags$td(sprintf("%.4f", mantel_res$r_statistic))),
          tags$tr(tags$td(strong("Significance (p-value)")), tags$td(mantel_res$p_string)),
          tags$tr(tags$td(strong("Permutations")), tags$td(mantel_res$permutations)),
          tags$tr(tags$td(strong("Matched Taxa Pairs")), tags$td(nrow(matched_df)))
        )
      ),
      div(
        class = alert_class,
        style = "padding: 10px; margin-bottom: 0;",
        HTML(paste("<b>Conclusion:</b>", conclusion_text))
      )
    )
  })

  output$struct_vs_phylo_plot <- renderPlot({

    if (input$run_mantel == 0) {

      plot.new()

      text(
        0.5,
        0.5,
        "Click 'Run Mantel Test' to generate congruence plot.",
        cex = 1.1,
        col = "#555555"
      )

      return()
    }

    res <- mantel_reactive()

    if (!is.null(res$error)) {

      plot.new()

      text(
        0.5,
        0.5,
        paste0(
          "Analysis Error:\n",
          res$error
        ),
        cex = 0.9,
        col = "#D9534F"
      )

      return()
    }

    plot_mantel_publication(
      res$matched_df,
      res$mantel_res,
      virus_name = input$virus_label
    )
  })

  output$download_mantel_pdf <- downloadHandler(

    filename = function() {
      paste0(
        "Mantel_Validation_Plot_",
        Sys.Date(),
        ".pdf"
      )
    },

    content = function(file) {

      res <- mantel_reactive()

      if (!is.null(res$error)) {
        stop(res$error)
      }

      p <- plot_mantel_publication(
        res$matched_df,
        res$mantel_res,
        virus_name = input$virus_label
      )

      pdf(
        file,
        width = 7,
        height = 5
      )

      tryCatch(
        {
          print(p)
        },
        finally = {
          dev.off()
        }
      )
    }
  )

  mod_treetime_server("treetime_mod")

  observeEvent(input$patch_xml, {

    req(
      input$beast2_xml
    )

    tryCatch({

      patch_result <- patch_beast2_xml(
        input$beast2_xml$datapath
      )

      app_state$patched_xml_log <- patch_result$msg

      if (patch_result$success) {
        app_state$patched_xml_doc <- patch_result$doc
      }

    }, error = function(e) {

      app_state$patched_xml_log <- paste(
        "Error:",
        e$message
      )
    })
  })

  output$xml_patch_log <- renderUI({
    if (app_state$patched_xml_log == "Waiting for XML upload...") {
      p(app_state$patched_xml_log, style = "color: gray; font-style: italic;")
    } else if (grepl("Error:", app_state$patched_xml_log)) {
      div(
        class = "alert alert-danger",
        style = "padding: 10px;",
        HTML(paste("<b>Patch Failed:</b><br>", app_state$patched_xml_log))
      )
    } else {
      div(
        class = "alert alert-success",
        style = "padding: 10px;",
        HTML(gsub("\n", "<br>", app_state$patched_xml_log))
      )
    }
  })

  output$resume_command_hint <- renderUI({
    div(
      style = "background-color: #f8f9fa; padding: 15px; border-left: 4px solid #2C3E50; border-radius: 4px;",
      HTML("1. Download the patched XML file.<br>
            2. Put it in the same folder as your original <code>.state</code> file.<br>
            3. Run the following command in your terminal:<br><br>
            <pre style='background-color: #2C3E50; color: #F8F8F2; padding: 10px; border-radius: 4px; margin-bottom: 0;'>beast -resume original_run.state patched_run.xml</pre>")
    )
  })

  output$download_patched_xml <- downloadHandler(

    filename = function() {
      paste0(
        "patched_resume_",
        Sys.Date(),
        ".xml"
      )
    },

    content = function(file) {

      req(
        app_state$patched_xml_doc
      )

      xml2::write_xml(
        app_state$patched_xml_doc,
        file
      )
    }
  )

  observeEvent(input$generate_rmsd, {
    req(input$pdb_dir)
    
    app_state$rmsd_status_ui <- "<span style='color:gray; font-style: italic;'>Initializing generation process...</span>"
    app_state$rmsd_matrix_data <- NULL 
    
    withProgress(
      message = paste("Running", input$rmsd_tool, "..."),
      detail = "This may take several minutes.",
      value = 0.5,
      {
        tryCatch({
          
          tmp_csv <- tempfile(fileext = ".csv")
          
          log_output <- capture.output({
            run_rmsd_pipeline(
              engine = input$rmsd_tool,
              pdb_dir = input$pdb_dir,
              output_csv = tmp_csv,
              custom_exe_path = trimws(input$custom_exe)
            )
          }, type = "message")
          
          app_state$rmsd_matrix_data <- read.csv(tmp_csv, row.names = 1, check.names = FALSE)
          
          app_state$rmsd_status_ui <- as.character(
            tagList(
              div(
                class = "alert alert-success",
                HTML("<b>Matrix Generation Successfully Completed!</b><br>
                     <span style='font-size:0.9em; color:#2C3E50;'>The structural RMSD matrix has been generated. You can view the Heatmap, Data Table, and download the results from the adjacent tabs. You can also directly upload the downloaded CSV in Tab 2 for Structural Validation.</span>")
              ),
              h6("Execution Trace:", style = "margin-top: 15px; font-weight: bold;"),
              tags$pre(
                style = "background-color: #2C3E50; color: #FFFFFF; padding: 10px; border-radius: 5px; font-size: 0.85em; max-height: 300px; overflow-y: auto;",
                paste(log_output, collapse = "\n")
              )
            )
          )
          
        }, error = function(e) {
          app_state$rmsd_status_ui <- as.character(
            div(
              class = "alert alert-danger",
              HTML(paste0("<b>Process Failed:</b><br>", e$message, 
                          "<br><br><span style='font-size:0.9em;'>Tip: Ensure the directory path is correct and the alignment engine is installed in your system PATH.</span>"))
            )
          )
        })
      }
    )
  })
  
  output$rmsd_status_ui <- renderUI({
    HTML(app_state$rmsd_status_ui)
  })
  
  output$rmsd_heatmap <- renderPlotly({
    req(app_state$rmsd_matrix_data)
    mat <- as.matrix(app_state$rmsd_matrix_data)
    
    plot_ly(
      x = colnames(mat), 
      y = rownames(mat), 
      z = mat, 
      type = "heatmap", 
      colorscale = "Viridis", 
      reversescale = TRUE 
    ) %>%
      layout(
        title = "Structural RMSD Clustering",
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "")
      )
  })
  
  output$rmsd_table <- renderDT({
    req(app_state$rmsd_matrix_data)
    datatable(
      app_state$rmsd_matrix_data,
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })
  
  output$download_rmsd_csv <- downloadHandler(
    filename = function() { paste0("RMSD_Matrix_", Sys.Date(), ".csv") },
    content = function(file) {
      req(app_state$rmsd_matrix_data)
      write.csv(app_state$rmsd_matrix_data, file, row.names = TRUE)
    }
  )
  
  output$download_heatmap_pdf <- downloadHandler(
    filename = function() { paste0("RMSD_Heatmap_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(app_state$rmsd_matrix_data)
      mat <- as.matrix(app_state$rmsd_matrix_data)
      
      pdf(file, width = 8, height = 7)
      df <- as.data.frame(as.table(mat))
      colnames(df) <- c("Y", "X", "RMSD")
      df$Y <- factor(df$Y, levels = rev(rownames(mat)))
      p <- ggplot(df, aes(x = X, y = Y, fill = RMSD)) +
        geom_tile(color = "white") +
        scale_fill_gradientn(colors = rev(hcl.colors(50, "viridis"))) + 
        theme_minimal() +
        labs(
          title = paste("Structural RMSD Clustering (", input$rmsd_tool, ")", sep=""),
          x = "", y = "", fill = "RMSD"
        ) +
        theme(
          axis.text.x = element_text(angle = -45, hjust = 0),
          plot.title = element_text(face = "bold", hjust = 0.5)
        )
      
      print(p)
      dev.off()
    }
  )


  mod_about_server("about_tab")
  
}

shinyApp(ui, server)
