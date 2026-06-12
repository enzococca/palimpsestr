## ============================================================
## palimpsestr Shiny Dashboard
## Interactive interface for Stratigraphic Entanglement Field analysis
## ============================================================

library(shiny)
library(shinydashboard)
library(palimpsestr)
library(DT)

has_plotly  <- requireNamespace("plotly", quietly = TRUE)
has_sf      <- requireNamespace("sf", quietly = TRUE)
has_sqlite  <- requireNamespace("RSQLite", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)
has_pg      <- requireNamespace("RPostgres", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)
has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)

# Helpers for plotly/ggplot fallback
sef_plot_ui <- function(id, height = "500px") {
 if (has_plotly) plotly::plotlyOutput(id, height = height)
 else plotOutput(id, height = height)
}

sef_plot_render <- function(input, output, session, id, plot_fn) {
  if (has_plotly) {
    output[[id]] <- plotly::renderPlotly({
      p <- plot_fn()
      if (is.null(p)) return(NULL)
      tryCatch(as_plotly(p), error = function(e) plotly::ggplotly(p))
    })
  } else {
    output[[id]] <- renderPlot({ plot_fn() })
  }
}

## ============================================================
## UI
## ============================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "palimpsestr"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Data", tabName = "data", icon = icon("upload")),
      menuItem("Analysis", tabName = "analysis", icon = icon("cogs")),
      menuItem("Plots", tabName = "results", icon = icon("chart-bar")),
      menuItem("Tables", tabName = "tables", icon = icon("table")),
      menuItem("Type longevity", tabName = "longevity", icon = icon("clock")),
      menuItem("Maps", tabName = "maps", icon = icon("map")),
      menuItem("Report", tabName = "report", icon = icon("file-alt"))
    )
  ),
  dashboardBody(
    tabItems(

      ## --- Tab 1: Data Input ---
      tabItem(tabName = "data",
        fluidRow(
          box(title = "CSV / Excel", status = "primary", solidHeader = TRUE, width = 4,
              fileInput("csv_file", "Upload CSV / Excel",
                        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
              radioButtons("csv_separator", "Separator (for CSV)",
                           choices = c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t"),
                           selected = ",", inline = TRUE)
          ),
          box(title = "SQLite", status = "info", solidHeader = TRUE, width = 4,
              if (has_sqlite) tagList(
                fileInput("sqlite_file", "SQLite Database", accept = c(".sqlite", ".db")),
                textInput("sqlite_table", "Table", value = "inventario_materiali_table"),
                actionButton("load_sqlite", "Load", icon = icon("database"))
              ) else p("Install RSQLite to use SQLite")
          ),
          box(title = "PostgreSQL", status = "warning", solidHeader = TRUE, width = 4,
              if (has_pg) tagList(
                fluidRow(
                  column(8, textInput("pg_host", "Host", value = "localhost")),
                  column(4, numericInput("pg_port", "Port", value = 5432))
                ),
                textInput("pg_dbname", "Database", value = ""),
                fluidRow(
                  column(6, textInput("pg_user", "User", value = "postgres")),
                  column(6, passwordInput("pg_pass", "Password"))
                ),
                textInput("pg_table", "Table", value = "inventario_materiali_table"),
                textInput("pg_query", "SQL Query (optional)", value = ""),
                actionButton("load_pg", "Connect", icon = icon("plug"))
              ) else p("Install RPostgres to use PostgreSQL")
          )
        ),
        fluidRow(
          box(title = "Column Mapping", status = "success", solidHeader = TRUE,
              width = 12, collapsible = TRUE,
              uiOutput("column_mapping"),
              actionButton("apply_mapping", "Apply Mapping", icon = icon("check"),
                           class = "btn-success")
          )
        ),
        fluidRow(
          box(title = "Merge taf_score from External File", status = "info",
              solidHeader = TRUE, width = 12, collapsible = TRUE, collapsed = TRUE,
              p("Upload a CSV or Excel file with taf_score per context (US). ",
                "The file must have a column matching the context IDs and a taf_score column."),
              fluidRow(
                column(4, fileInput("taf_file", "taf_score File (CSV/Excel)",
                                    accept = c(".csv", ".tsv", ".xlsx", ".xls"))),
                column(3, textInput("taf_join_col", "Context column in taf file", value = "US")),
                column(3, textInput("taf_value_col", "taf_score column", value = "taf_score")),
                column(2, br(), actionButton("merge_taf", "Merge taf_score",
                                             icon = icon("link"), class = "btn-info"))
              ),
              verbatimTextOutput("taf_merge_status")
          )
        ),
        fluidRow(
          box(title = "Calibrated radiocarbon import (rcarbon)",
              status = "info", solidHeader = TRUE, width = 12,
              collapsible = TRUE, collapsed = TRUE,
              p("Optional. Upload a CSV with columns ",
                tags$code("lab_id"), ", ",
                tags$code("cra"), ", ",
                tags$code("error"),
                " (uncalibrated 14C BP and 1-sigma error). The dates will ",
                "be calibrated with rcarbon and merged into the dataset as ",
                tags$code("date_min"), "/", tags$code("date_max"),
                "/", tags$code("date_mid"), "."),
              fluidRow(
                column(4, fileInput("rcarbon_csv",
                                    "CSV with columns lab_id, cra, error",
                                    accept = ".csv")),
                column(4, selectInput("rcarbon_method",
                                      "Reduction method",
                                      choices = c("hpd", "median_iqr",
                                                  "weighted_mean"))),
                column(4, br(), actionButton("rcarbon_run",
                                             "Calibrate and merge",
                                             icon = icon("calendar-check"),
                                             class = "btn-info"))
              ),
              verbatimTextOutput("rcarbon_status")
          )
        ),
        fluidRow(
          box(title = "Data Preview", width = 12, collapsible = TRUE,
              verbatimTextOutput("data_summary"),
              DT::dataTableOutput("data_preview")
          )
        )
      ),

      ## --- Tab 2: Analysis ---
      tabItem(tabName = "analysis",
        fluidRow(
          box(title = "Model Parameters", status = "primary", solidHeader = TRUE, width = 4,
              sliderInput("k", "Number of Phases (K)", min = 2, max = 10, value = 4),
              checkboxInput("use_taf", "Taphonomic weighting (taf_score)", TRUE),
              checkboxInput("use_harris", "Harris constraint (auto from depth)", FALSE),
              numericInput("n_init", "Initialisations", value = 5, min = 1, max = 50),
              numericInput("em_iter", "Max EM iterations", value = 100, min = 10, max = 1000),
              numericInput("seed", "Seed", value = 42),
              hr(),
              h5("Feature Space"),
              checkboxInput("chrono_precision", "Chronological precision (1/tspan)", FALSE),
              checkboxInput("taf_as_feature", "taf_score as feature", FALSE),
              checkboxInput("residuality", "Residuality (date/context mismatch)", FALSE),
              checkboxInput("class_scale", "One-hot class scaling", FALSE),
              textInput("subclass", "Sub-class column (optional)", value = ""),
              actionButton("run_analysis", "Run Analysis",
                           icon = icon("play"), class = "btn-primary btn-lg")
          ),
          box(title = "Phase Count Selection (compare_k)", status = "info", solidHeader = TRUE, width = 8,
              fluidRow(
                column(6, sliderInput("k_range", "K Range", min = 2, max = 10, value = c(2, 7))),
                column(3, numericInput("ck_ninit", "n_init", value = 3, min = 1)),
                column(3, br(), actionButton("run_compare_k", "Compare K", icon = icon("search")))
              ),
              sef_plot_ui("plot_compare_k", "400px"),
              DT::dataTableOutput("tbl_compare_k")
          )
        ),
        fluidRow(
          box(title = "Status", width = 12,
              verbatimTextOutput("analysis_status")
          )
        )
      ),

      ## --- Tab 3: Results (plots) ---
      tabItem(tabName = "results",
        tabBox(width = 12, id = "result_tabs",
          tabPanel("Phase Field", sef_plot_ui("plot_phasefield", "550px")),
          tabPanel("Entropy", sef_plot_ui("plot_entropy", "550px")),
          tabPanel("Energy (ESE)", sef_plot_ui("plot_energy", "550px")),
          tabPanel("Intrusions",
                   numericInput("top_n", "Top N suspects", value = 10, min = 1, max = 50),
                   sef_plot_ui("plot_intrusions", "550px")),
          tabPanel("Vertical Profile", sef_plot_ui("plot_phase_profile", "550px")),
          tabPanel("EM Convergence", sef_plot_ui("plot_convergence", "550px"))
        )
      ),

      ## --- Tab 4: Tables ---
      tabItem(tabName = "tables",
        tabBox(width = 12, id = "table_tabs",
          tabPanel("Phases", DT::dataTableOutput("tbl_phases")),
          tabPanel("Intrusions", DT::dataTableOutput("tbl_intrusions")),
          tabPanel("US Purity", DT::dataTableOutput("tbl_us_summary")),
          tabPanel("Transition Matrix", DT::dataTableOutput("tbl_transition")),
          tabPanel("Model",
                   verbatimTextOutput("txt_print_fit"),
                   verbatimTextOutput("txt_summary_fit"))
        )
      ),

      ## --- Tab: Type longevity ---
      tabItem(tabName = "longevity",
        fluidRow(
          box(title = "Type longevity", status = "primary", solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(4, numericInput("longevity_threshold",
                                       "Posterior threshold",
                                       value = 0.1, min = 0, max = 1, step = 0.05)),
                column(4, br(), downloadButton("longevity_download",
                                               "Download (xlsx)"))
              ),
              plotOutput("longevity_plot", height = "400px"),
              DT::dataTableOutput("longevity_table")
          )
        )
      ),

      ## --- Tab 5: Maps ---
      tabItem(tabName = "maps",
        if (has_sf) tagList(
          fluidRow(
            box(title = "Geometries", status = "primary", solidHeader = TRUE, width = 4,
                fileInput("geom_file", "GeoPackage / GeoJSON / Shapefile",
                          accept = c(".gpkg", ".geojson", ".shp", ".shx", ".dbf", ".prj"),
                          multiple = TRUE),
                textInput("geom_context_col", "US column in geometries", value = "context"),
                selectInput("map_layer", "Layer", choices = c(
                  "Phases" = "phase", "Entropy" = "entropy",
                  "Energy" = "energy", "Intrusions" = "intrusion"
                )),
                actionButton("render_map", "Render Map", icon = icon("map"),
                             class = "btn-primary")
            ),
            box(title = "Excavation Plan Map", width = 8,
                sef_plot_ui("plot_map", "600px")
            )
          )
        ) else h3("Install the sf package for GIS maps")
      ),

      ## --- Tab 6: Report ---
      tabItem(tabName = "report",
        fluidRow(
          box(title = "Interpretive Report", status = "primary", solidHeader = TRUE, width = 4,
              selectInput("report_lang", "Language", choices = c("English" = "en", "Italiano" = "it")),
              actionButton("gen_report", "Generate Report", icon = icon("file-text")),
              hr(),
              h4("Download Results"),
              downloadButton("dl_phases_csv", "Phases (CSV)"),
              br(), br(),
              downloadButton("dl_intrusions_csv", "Intrusions (CSV)"),
              br(), br(),
              downloadButton("dl_us_summary_csv", "US Purity (CSV)"),
              br(), br(),
              downloadButton("dl_all_zip", "All Results (ZIP)")
          ),
          box(title = "Report", width = 8,
              verbatimTextOutput("report_text")
          )
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================

server <- function(input, output, session) {

  rv <- reactiveValues(
    raw_data   = NULL,
    data       = NULL,
    fit        = NULL,
    compare_k  = NULL,
    geometries = NULL
  )

  ## --- Data Loading ---

  # CSV / Excel
  observeEvent(input$csv_file, {
    tryCatch({
      ext <- tolower(tools::file_ext(input$csv_file$name))
      if (ext %in% c("xlsx", "xls")) {
        if (!has_openxlsx)
          stop("Install the openxlsx package to import Excel files")
        rv$raw_data <- openxlsx::read.xlsx(input$csv_file$datapath)
      } else {
        sep <- input$csv_separator
        rv$raw_data <- read.csv(input$csv_file$datapath, sep = sep, stringsAsFactors = FALSE)
      }
      showNotification(paste("Loaded", nrow(rv$raw_data), "records"), type = "message")
    }, error = function(e) showNotification(paste("Import error:", e$message), type = "error"))
  })

  # SQLite
  observeEvent(input$load_sqlite, {
    req(input$sqlite_file)
    tryCatch({
      con <- DBI::dbConnect(RSQLite::SQLite(), input$sqlite_file$datapath)
      on.exit(DBI::dbDisconnect(con))
      rv$raw_data <- DBI::dbReadTable(con, input$sqlite_table)
      showNotification(paste("Loaded", nrow(rv$raw_data), "records from SQLite"), type = "message")
    }, error = function(e) showNotification(paste("SQLite error:", e$message), type = "error"))
  })

  # PostgreSQL
  observeEvent(input$load_pg, {
    tryCatch({
      con <- DBI::dbConnect(RPostgres::Postgres(),
                            host = input$pg_host, port = input$pg_port,
                            dbname = input$pg_dbname, user = input$pg_user,
                            password = input$pg_pass)
      on.exit(DBI::dbDisconnect(con))
      if (nzchar(input$pg_query)) {
        rv$raw_data <- DBI::dbGetQuery(con, input$pg_query)
      } else {
        rv$raw_data <- DBI::dbReadTable(con, input$pg_table)
      }
      showNotification(paste("Loaded", nrow(rv$raw_data), "records from PostgreSQL"), type = "message")
    }, error = function(e) showNotification(paste("PostgreSQL error:", e$message), type = "error"))
  })

  ## --- Column Mapping ---

  output$column_mapping <- renderUI({
    req(rv$raw_data)
    cols <- names(rv$raw_data)
    none <- c("(none)" = "")

    auto <- function(target) {
      m <- match(tolower(target), tolower(cols))
      if (!is.na(m)) cols[m] else ""
    }

    tagList(
      fluidRow(
        column(2, selectInput("col_x", "x (easting)", choices = c(none, cols), selected = auto("x"))),
        column(2, selectInput("col_y", "y (northing)", choices = c(none, cols), selected = auto("y"))),
        column(2, selectInput("col_z", "z (elevation)", choices = c(none, cols), selected = auto("z"))),
        column(2, selectInput("col_dmin", "date_min", choices = c(none, cols), selected = auto("date_min"))),
        column(2, selectInput("col_dmax", "date_max", choices = c(none, cols), selected = auto("date_max"))),
        column(2, selectInput("col_class", "class", choices = c(none, cols), selected = auto("class")))
      ),
      fluidRow(
        column(2, selectInput("col_id", "id (opt.)", choices = c(none, cols), selected = auto("id"))),
        column(2, selectInput("col_context", "context (opt.)", choices = c(none, cols), selected = auto("context"))),
        column(2, selectInput("col_taf", "taf_score (opt.)", choices = c(none, cols), selected = auto("taf_score")))
      )
    )
  })

  observeEvent(input$apply_mapping, {
    req(rv$raw_data)
    tryCatch({
      d <- rv$raw_data
      out <- data.frame(
        x        = as.numeric(d[[input$col_x]]),
        y        = as.numeric(d[[input$col_y]]),
        z        = as.numeric(d[[input$col_z]]),
        date_min = as.numeric(d[[input$col_dmin]]),
        date_max = as.numeric(d[[input$col_dmax]]),
        class    = as.character(d[[input$col_class]]),
        stringsAsFactors = FALSE
      )
      if (nzchar(input$col_id))      out$id <- as.character(d[[input$col_id]])
      else                            out$id <- paste0("R_", seq_len(nrow(out)))
      if (nzchar(input$col_context))  out$context <- as.character(d[[input$col_context]])
      if (nzchar(input$col_taf))      out$taf_score <- as.numeric(d[[input$col_taf]])

      ok <- complete.cases(out[, c("x", "y", "z", "date_min", "date_max", "class")])
      out <- out[ok, ]

      rv$data <- out
      rv$fit <- NULL
      rv$compare_k <- NULL
      showNotification(paste("Dataset ready:", nrow(out), "records,",
                             length(unique(out$class)), "classes"), type = "message")
    }, error = function(e) showNotification(paste("Mapping error:", e$message), type = "error"))
  })

  ## --- Merge taf_score from external file ---

  observeEvent(input$merge_taf, {
    req(rv$data, input$taf_file)
    tryCatch({
      # Read external taf file
      ext <- tolower(tools::file_ext(input$taf_file$name))
      if (ext %in% c("xlsx", "xls")) {
        if (!has_openxlsx) stop("Install openxlsx to import Excel files")
        taf_df <- openxlsx::read.xlsx(input$taf_file$datapath)
      } else {
        # Try comma first, then semicolon, then tab
        taf_df <- tryCatch(
          read.csv(input$taf_file$datapath, stringsAsFactors = FALSE),
          error = function(e) read.csv(input$taf_file$datapath, sep = ";", stringsAsFactors = FALSE)
        )
      }

      join_col <- input$taf_join_col
      val_col  <- input$taf_value_col

      if (!join_col %in% names(taf_df))
        stop(paste0("Column '", join_col, "' not found in taf file. Available: ",
                     paste(names(taf_df), collapse = ", ")))
      if (!val_col %in% names(taf_df))
        stop(paste0("Column '", val_col, "' not found in taf file. Available: ",
                     paste(names(taf_df), collapse = ", ")))

      # Build lookup: context -> taf_score
      taf_lookup <- taf_df[[val_col]]
      names(taf_lookup) <- as.character(taf_df[[join_col]])
      taf_lookup <- as.numeric(taf_lookup)

      # Match contexts in dataset
      if (!"context" %in% names(rv$data))
        stop("Dataset has no 'context' column. Map it first in Column Mapping.")

      # Try matching with and without "US_" prefix
      ctx <- as.character(rv$data$context)
      matched <- taf_lookup[ctx]

      # If no matches, try stripping "US_" prefix from dataset
      if (all(is.na(matched))) {
        ctx_stripped <- sub("^US_", "", ctx)
        matched <- taf_lookup[ctx_stripped]
      }
      # If still no matches, try adding "US_" prefix
      if (all(is.na(matched))) {
        ctx_prefixed <- paste0("US_", ctx)
        matched <- taf_lookup[ctx_prefixed]
      }

      n_matched <- sum(!is.na(matched))
      if (n_matched == 0)
        stop("No matching contexts found. Check context column names.")

      # Merge: update existing taf_score or create new
      rv$data$taf_score <- ifelse(is.na(matched),
                                   if ("taf_score" %in% names(rv$data)) rv$data$taf_score else 0.5,
                                   matched)

      # Also merge any extra columns (sotto_tipologia, residualita, etc.)
      extra_cols <- setdiff(names(taf_df), c(join_col, val_col))
      for (ec in extra_cols) {
        vals <- taf_df[[ec]]
        names(vals) <- as.character(taf_df[[join_col]])
        ec_matched <- vals[ctx]
        if (all(is.na(ec_matched))) ec_matched <- vals[sub("^US_", "", ctx)]
        if (all(is.na(ec_matched))) ec_matched <- vals[paste0("US_", ctx)]
        if (any(!is.na(ec_matched))) {
          rv$data[[ec]] <- as.character(ec_matched)
        }
      }

      showNotification(paste("taf_score merged:", n_matched, "/", nrow(rv$data),
                             "finds matched"), type = "message")
    }, error = function(e) showNotification(paste("Merge error:", e$message), type = "error"))
  })

  output$taf_merge_status <- renderPrint({
    req(rv$data)
    if ("taf_score" %in% names(rv$data)) {
      cat("taf_score present. Mean:", round(mean(rv$data$taf_score, na.rm = TRUE), 3),
          "| Range:", paste(range(rv$data$taf_score, na.rm = TRUE), collapse = "-"),
          "| NA:", sum(is.na(rv$data$taf_score)), "\n")
    } else {
      cat("No taf_score in dataset yet.\n")
    }
  })

  ## --- Data Preview ---

  output$data_summary <- renderPrint({
    req(rv$data)
    cat("Finds:", nrow(rv$data), "\n")
    if ("context" %in% names(rv$data))
      cat("Contexts:", length(unique(rv$data$context)), "\n")
    cat("Classes:", length(unique(rv$data$class)), "\n")
    cat("Date range:", min(rv$data$date_min, na.rm = TRUE), "to",
        max(rv$data$date_max, na.rm = TRUE), "\n")
    cat("Z range:", round(min(rv$data$z, na.rm = TRUE), 2), "-",
        round(max(rv$data$z, na.rm = TRUE), 2), "\n")
    if ("taf_score" %in% names(rv$data))
      cat("Mean taf_score:", round(mean(rv$data$taf_score, na.rm = TRUE), 3), "\n")
  })

  output$data_preview <- DT::renderDataTable({
    req(rv$data)
    DT::datatable(rv$data, options = list(scrollX = TRUE, pageLength = 10))
  })

  ## --- Analysis ---

  observeEvent(input$run_analysis, {
    req(rv$data)
    withProgress(message = "Fitting SEF model...", {
      tryCatch({
        taf_arg <- if (input$use_taf && "taf_score" %in% names(rv$data)) "taf_score" else NULL
        ctx_arg <- if ("context" %in% names(rv$data)) "context" else NULL
        h_arg <- NULL
        if (input$use_harris && !is.null(ctx_arg)) {
          incProgress(0.1, detail = "Building Harris matrix...")
          h_arg <- harris_from_contexts(rv$data, z_col = "z", context_col = "context")
        }
        incProgress(0.2, detail = paste("Fitting K =", input$k, "..."))
        subclass_arg <- if (nzchar(input$subclass)) input$subclass else NULL
        rv$fit <- fit_sef(rv$data, k = input$k,
                          tafonomy = taf_arg, context = ctx_arg, harris = h_arg,
                          n_init = input$n_init, em_iter = input$em_iter,
                          seed = input$seed,
                          chrono_precision = input$chrono_precision,
                          taf_as_feature = input$taf_as_feature,
                          residuality = input$residuality,
                          class_scale = input$class_scale,
                          subclass = subclass_arg)
        showNotification(paste("Fit complete! PDI =", round(pdi(rv$fit), 4)), type = "message")
      }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
    })
  })

  observeEvent(input$run_compare_k, {
    req(rv$data)
    withProgress(message = "Comparing K values...", {
      tryCatch({
        taf_arg <- if (input$use_taf && "taf_score" %in% names(rv$data)) "taf_score" else NULL
        ctx_arg <- if ("context" %in% names(rv$data)) "context" else NULL
        h_arg <- NULL
        if (input$use_harris && !is.null(ctx_arg)) {
          h_arg <- harris_from_contexts(rv$data, z_col = "z", context_col = "context")
        }
        kvals <- seq(input$k_range[1], input$k_range[2])
        subclass_arg <- if (nzchar(input$subclass)) input$subclass else NULL
        rv$compare_k <- compare_k(rv$data, k_values = kvals,
                                   tafonomy = taf_arg, context = ctx_arg, harris = h_arg,
                                   seed = input$seed, n_init = input$ck_ninit,
                                   chrono_precision = input$chrono_precision,
                                   taf_as_feature = input$taf_as_feature,
                                   residuality = input$residuality,
                                   class_scale = input$class_scale,
                                   subclass = subclass_arg)
        showNotification("K comparison complete!", type = "message")
      }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
    })
  })

  output$analysis_status <- renderPrint({
    if (!is.null(rv$fit)) {
      print(rv$fit)
      cat("\nPDI:", pdi(rv$fit), "\n")
      di <- detect_intrusions(rv$fit)
      cat("Intrusions (prob > 0.5):", sum(di$intrusion_prob > 0.5), "/", nrow(rv$data), "\n")
    } else {
      cat("No analysis started. Load data and press 'Run Analysis'.\n")
    }
  })

  ## --- Compare K ---

  sef_plot_render(input, output, session, "plot_compare_k", function() {
    req(rv$compare_k)
    gg_compare_k(rv$compare_k)
  })

  output$tbl_compare_k <- DT::renderDataTable({
    req(rv$compare_k)
    DT::datatable(rv$compare_k, options = list(scrollX = TRUE, pageLength = 10)) |>
      DT::formatRound(columns = c("pdi", "mean_entropy", "bic", "icl"), digits = 4)
  })

  ## --- Result Plots ---

  sef_plot_render(input, output, session, "plot_phasefield", function() {
    req(rv$fit); gg_phasefield(rv$fit)
  })
  sef_plot_render(input, output, session, "plot_entropy", function() {
    req(rv$fit); gg_entropy(rv$fit)
  })
  sef_plot_render(input, output, session, "plot_energy", function() {
    req(rv$fit); gg_energy(rv$fit)
  })
  sef_plot_render(input, output, session, "plot_intrusions", function() {
    req(rv$fit); gg_intrusions(rv$fit, top_n = input$top_n)
  })
  sef_plot_render(input, output, session, "plot_phase_profile", function() {
    req(rv$fit); gg_phase_profile(rv$fit)
  })
  sef_plot_render(input, output, session, "plot_convergence", function() {
    req(rv$fit); gg_convergence(rv$fit)
  })

  ## --- Tables ---

  output$tbl_phases <- DT::renderDataTable({
    req(rv$fit)
    DT::datatable(as_phase_table(rv$fit),
                  options = list(scrollX = TRUE, pageLength = 15)) |>
      DT::formatRound(columns = c("entropy", "energy"), digits = 4)
  })

  output$tbl_intrusions <- DT::renderDataTable({
    req(rv$fit)
    di <- detect_intrusions(rv$fit)
    di <- di[order(di$intrusion_prob, decreasing = TRUE), ]
    if ("direction" %in% names(di)) {
      dir_chr <- as.character(di$direction)
      di$dir_icon <- ifelse(is.na(dir_chr), "",
                     ifelse(dir_chr == "older_than_context",   "DOWN",
                     ifelse(dir_chr == "younger_than_context", "UP",
                     ifelse(dir_chr == "in_context",           "=", ""))))
    }
    round_cols <- intersect(c("intrusion_prob", "chrono_gap"), names(di))
    DT::datatable(di, options = list(scrollX = TRUE, pageLength = 15)) |>
      DT::formatRound(columns = round_cols, digits = 4)
  })

  output$tbl_us_summary <- DT::renderDataTable({
    req(rv$fit)
    tryCatch({
      DT::datatable(us_summary_table(rv$fit),
                    options = list(scrollX = TRUE, pageLength = 20)) |>
        DT::formatPercentage(columns = "purity", digits = 1) |>
        DT::formatRound(columns = c("mean_entropy", "mean_energy", "mean_local_sei"), digits = 4)
    }, error = function(e) DT::datatable(data.frame(note = "Context not available")))
  })

  output$tbl_transition <- DT::renderDataTable({
    req(rv$fit)
    tm <- phase_transition_matrix(rv$fit)
    DT::datatable(as.data.frame.matrix(tm), options = list(scrollX = TRUE))
  })

  output$txt_print_fit <- renderPrint({ req(rv$fit); print(rv$fit) })
  output$txt_summary_fit <- renderPrint({ req(rv$fit); summary(rv$fit) })

  ## --- Type longevity ---

  longevity_data <- reactive({
    req(rv$fit)
    type_longevity(rv$fit, posterior_threshold = input$longevity_threshold)
  })

  output$longevity_table <- DT::renderDataTable({
    tl <- longevity_data()
    show_cols <- intersect(c("class", "longevity_min", "longevity_max",
                             "longevity_span", "dominant_phase", "n_finds"),
                           names(tl))
    DT::datatable(tl[, show_cols, drop = FALSE],
                  options = list(scrollX = TRUE, pageLength = 15)) |>
      DT::formatRound(columns = intersect(c("longevity_min", "longevity_max",
                                            "longevity_span"), show_cols),
                      digits = 2)
  })

  output$longevity_plot <- renderPlot({
    tl <- longevity_data()
    gg_longevity(tl, rv$fit)
  })

  output$longevity_download <- downloadHandler(
    filename = function() paste0("type_longevity_", Sys.Date(), ".xlsx"),
    content = function(file) {
      if (!has_openxlsx) stop("Install the openxlsx package to export xlsx")
      tl <- longevity_data()
      if ("weight_matrix" %in% names(tl) && is.list(tl$weight_matrix)) {
        tl$weight_matrix <- vapply(tl$weight_matrix,
                                   function(w) paste(round(w, 3),
                                                     collapse = ";"), "")
      }
      openxlsx::write.xlsx(tl, file)
    }
  )

  ## --- rcarbon import ---

  observeEvent(input$rcarbon_run, {
    req(input$rcarbon_csv, rv$data)
    if (!requireNamespace("rcarbon", quietly = TRUE)) {
      showNotification("Install package 'rcarbon' first.", type = "error")
      return()
    }
    tryCatch({
      csv <- utils::read.csv(input$rcarbon_csv$datapath, stringsAsFactors = FALSE)
      req_cols <- c("lab_id", "cra", "error")
      missing_cols <- setdiff(req_cols, names(csv))
      if (length(missing_cols))
        stop("CSV is missing column(s): ",
             paste(missing_cols, collapse = ", "))
      cal <- rcarbon::calibrate(x = csv$cra, errors = csv$error,
                                verbose = FALSE)
      ch <- chronology_from_rcarbon(cal, method = input$rcarbon_method,
                                    ids = csv$lab_id)
      # Drop pre-existing date columns to avoid .x/.y suffixes after merge
      drop_cols <- intersect(c("date_min", "date_max", "date_mid"),
                             names(rv$data))
      d <- rv$data[, setdiff(names(rv$data), drop_cols), drop = FALSE]
      rv$data <- merge(d, ch, by.x = "id", by.y = "id", all.x = TRUE)
      showNotification(
        paste("Calibrated dates merged:", nrow(ch), "labels processed."),
        type = "message")
    }, error = function(e)
      showNotification(paste("rcarbon error:", e$message), type = "error"))
  })

  output$rcarbon_status <- renderPrint({
    if (!requireNamespace("rcarbon", quietly = TRUE)) {
      cat("Package 'rcarbon' not installed. ",
          "Run install.packages('rcarbon') to enable this feature.\n",
          sep = "")
    } else {
      cat("rcarbon ready. Upload a CSV with columns lab_id, cra, error.\n")
    }
  })

  ## --- Maps ---

  observeEvent(input$geom_file, {
    req(has_sf)
    tryCatch({
      files <- input$geom_file
      tmpdir <- tempfile()
      dir.create(tmpdir)
      for (i in seq_len(nrow(files))) {
        ext <- tools::file_ext(files$name[i])
        base <- tools::file_path_sans_ext(files$name[1])
        newpath <- file.path(tmpdir, paste0(base, ".", ext))
        file.copy(files$datapath[i], newpath)
      }
      main_file <- list.files(tmpdir, pattern = "\\.(gpkg|geojson|shp)$", full.names = TRUE)[1]
      geom <- sf::st_read(main_file, quiet = TRUE)
      geom <- sf::st_make_valid(geom)

      ctx_col <- input$geom_context_col
      if (ctx_col %in% names(geom) && ctx_col != "context") {
        geom$context <- as.character(geom[[ctx_col]])
      } else if (!"context" %in% names(geom)) {
        for (cn in c("us", "US", "us_s", "context", "nome")) {
          if (cn %in% names(geom)) { geom$context <- as.character(geom[[cn]]); break }
        }
      }

      us_list <- unique(geom$context)
      dissolved <- do.call(rbind, lapply(us_list, function(u) {
        sub <- geom[geom$context == u, ]
        sf::st_sf(context = u, geom = sf::st_union(sub))
      }))
      rv$geometries <- sf::st_cast(dissolved, "MULTIPOLYGON")

      showNotification(paste("Geometries loaded:", nrow(rv$geometries), "units"), type = "message")
    }, error = function(e) showNotification(paste("Geometry error:", e$message), type = "error"))
  })

  observeEvent(input$render_map, {
    req(rv$fit, rv$geometries)
    output$plot_map_trigger <- renderText({ input$render_map })
  })

  sef_plot_render(input, output, session, "plot_map", function() {
    req(rv$fit, rv$geometries)
    input$render_map
    gg_map(rv$fit, geometries = rv$geometries, layer = input$map_layer)
  })

  ## --- Report ---

  report_rv <- reactiveVal(NULL)

  observeEvent(input$gen_report, {
    req(rv$fit)
    tryCatch({
      txt <- capture.output(report_sef(rv$fit, lang = input$report_lang))
      report_rv(paste(txt, collapse = "\n"))
      showNotification("Report generated!", type = "message")
    }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
  })

  output$report_text <- renderPrint({
    if (!is.null(report_rv())) cat(report_rv())
    else cat("Press 'Generate Report' after running the analysis.")
  })

  ## --- Downloads ---

  output$dl_phases_csv <- downloadHandler(
    filename = function() paste0("palimpsestr_phases_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$fit)
      write.csv(as_phase_table(rv$fit), file, row.names = FALSE)
    }
  )

  output$dl_intrusions_csv <- downloadHandler(
    filename = function() paste0("palimpsestr_intrusions_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$fit)
      write.csv(detect_intrusions(rv$fit), file, row.names = FALSE)
    }
  )

  output$dl_us_summary_csv <- downloadHandler(
    filename = function() paste0("palimpsestr_us_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$fit)
      tryCatch(
        write.csv(us_summary_table(rv$fit), file, row.names = FALSE),
        error = function(e) write.csv(data.frame(note = "Context not available"), file, row.names = FALSE)
      )
    }
  )

  output$dl_all_zip <- downloadHandler(
    filename = function() paste0("palimpsestr_results_", Sys.Date(), ".zip"),
    content = function(file) {
      req(rv$fit)
      tmpdir <- tempfile()
      dir.create(tmpdir)
      export_results(rv$fit, dir = tmpdir, prefix = "palimpsestr")
      tryCatch({
        txt_en <- capture.output(report_sef(rv$fit, lang = "en"))
        writeLines(txt_en, file.path(tmpdir, "report_en.md"))
        txt_it <- capture.output(report_sef(rv$fit, lang = "it"))
        writeLines(txt_it, file.path(tmpdir, "report_it.md"))
      }, error = function(e) NULL)
      files <- list.files(tmpdir, full.names = TRUE)
      zip(file, files, flags = "-j")
    }
  )
}

## ============================================================
## Launch
## ============================================================
shinyApp(ui = ui, server = server)
