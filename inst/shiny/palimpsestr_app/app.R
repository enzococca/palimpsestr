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
      menuItem("Dati", tabName = "data", icon = icon("upload")),
      menuItem("Analisi", tabName = "analysis", icon = icon("cogs")),
      menuItem("Grafici", tabName = "results", icon = icon("chart-bar")),
      menuItem("Tabelle", tabName = "tables", icon = icon("table")),
      menuItem("Mappe", tabName = "maps", icon = icon("map")),
      menuItem("Report", tabName = "report", icon = icon("file-alt"))
    )
  ),
  dashboardBody(
    tabItems(

      ## --- Tab 1: Data Input ---
      tabItem(tabName = "data",
        fluidRow(
          box(title = "CSV", status = "primary", solidHeader = TRUE, width = 4,
              fileInput("csv_file", "Carica CSV", accept = ".csv"),
              checkboxInput("csv_sep_semicolon", "Separatore ;", FALSE)
          ),
          box(title = "SQLite", status = "info", solidHeader = TRUE, width = 4,
              if (has_sqlite) tagList(
                fileInput("sqlite_file", "Database SQLite", accept = c(".sqlite", ".db")),
                textInput("sqlite_table", "Tabella", value = "inventario_materiali_table"),
                actionButton("load_sqlite", "Carica", icon = icon("database"))
              ) else p("Installa RSQLite per usare SQLite")
          ),
          box(title = "PostgreSQL", status = "warning", solidHeader = TRUE, width = 4,
              if (has_pg) tagList(
                fluidRow(
                  column(8, textInput("pg_host", "Host", value = "localhost")),
                  column(4, numericInput("pg_port", "Porta", value = 5432))
                ),
                textInput("pg_dbname", "Database", value = ""),
                fluidRow(
                  column(6, textInput("pg_user", "Utente", value = "postgres")),
                  column(6, passwordInput("pg_pass", "Password"))
                ),
                textInput("pg_table", "Tabella", value = "inventario_materiali_table"),
                textInput("pg_query", "Query SQL (opzionale)", value = ""),
                actionButton("load_pg", "Connetti", icon = icon("plug"))
              ) else p("Installa RPostgres per usare PostgreSQL")
          )
        ),
        fluidRow(
          box(title = "Mappatura colonne", status = "success", solidHeader = TRUE,
              width = 12, collapsible = TRUE,
              uiOutput("column_mapping"),
              actionButton("apply_mapping", "Applica mappatura", icon = icon("check"),
                           class = "btn-success")
          )
        ),
        fluidRow(
          box(title = "Anteprima dati", width = 12, collapsible = TRUE,
              verbatimTextOutput("data_summary"),
              DT::dataTableOutput("data_preview")
          )
        )
      ),

      ## --- Tab 2: Analysis ---
      tabItem(tabName = "analysis",
        fluidRow(
          box(title = "Parametri modello", status = "primary", solidHeader = TRUE, width = 4,
              sliderInput("k", "Numero di fasi (K)", min = 2, max = 10, value = 4),
              checkboxInput("use_taf", "Ponderazione tafonomica (taf_score)", TRUE),
              checkboxInput("use_harris", "Vincolo Harris (auto da profondita')", FALSE),
              numericInput("n_init", "Inizializzazioni", value = 5, min = 1, max = 50),
              numericInput("em_iter", "Max iterazioni EM", value = 100, min = 10, max = 1000),
              numericInput("seed", "Seed", value = 42),
              actionButton("run_analysis", "Avvia analisi",
                           icon = icon("play"), class = "btn-primary btn-lg")
          ),
          box(title = "Selezione K (compare_k)", status = "info", solidHeader = TRUE, width = 8,
              fluidRow(
                column(6, sliderInput("k_range", "Range K", min = 2, max = 10, value = c(2, 7))),
                column(3, numericInput("ck_ninit", "n_init", value = 3, min = 1)),
                column(3, br(), actionButton("run_compare_k", "Confronta K", icon = icon("search")))
              ),
              sef_plot_ui("plot_compare_k", "400px"),
              DT::dataTableOutput("tbl_compare_k")
          )
        ),
        fluidRow(
          box(title = "Stato", width = 12,
              verbatimTextOutput("analysis_status")
          )
        )
      ),

      ## --- Tab 3: Results (plots) ---
      tabItem(tabName = "results",
        tabBox(width = 12, id = "result_tabs",
          tabPanel("Phase Field", sef_plot_ui("plot_phasefield", "550px")),
          tabPanel("Entropia", sef_plot_ui("plot_entropy", "550px")),
          tabPanel("Energia (ESE)", sef_plot_ui("plot_energy", "550px")),
          tabPanel("Intrusioni",
                   numericInput("top_n", "Top N sospetti", value = 10, min = 1, max = 50),
                   sef_plot_ui("plot_intrusions", "550px")),
          tabPanel("Profilo verticale", sef_plot_ui("plot_phase_profile", "550px")),
          tabPanel("Convergenza EM", sef_plot_ui("plot_convergence", "550px"))
        )
      ),

      ## --- Tab 4: Tables ---
      tabItem(tabName = "tables",
        tabBox(width = 12, id = "table_tabs",
          tabPanel("Fasi", DT::dataTableOutput("tbl_phases")),
          tabPanel("Intrusioni", DT::dataTableOutput("tbl_intrusions")),
          tabPanel("Purezza US", DT::dataTableOutput("tbl_us_summary")),
          tabPanel("Matrice transizione", DT::dataTableOutput("tbl_transition")),
          tabPanel("Modello",
                   verbatimTextOutput("txt_print_fit"),
                   verbatimTextOutput("txt_summary_fit"))
        )
      ),

      ## --- Tab 5: Maps ---
      tabItem(tabName = "maps",
        if (has_sf) tagList(
          fluidRow(
            box(title = "Geometrie", status = "primary", solidHeader = TRUE, width = 4,
                fileInput("geom_file", "GeoPackage / GeoJSON / Shapefile",
                          accept = c(".gpkg", ".geojson", ".shp", ".shx", ".dbf", ".prj"),
                          multiple = TRUE),
                textInput("geom_context_col", "Colonna US nelle geometrie", value = "context"),
                selectInput("map_layer", "Layer", choices = c(
                  "Fasi" = "phase", "Entropia" = "entropy",
                  "Energia" = "energy", "Intrusioni" = "intrusion"
                )),
                actionButton("render_map", "Genera mappa", icon = icon("map"),
                             class = "btn-primary")
            ),
            box(title = "Mappa su pianta", width = 8,
                sef_plot_ui("plot_map", "600px")
            )
          )
        ) else h3("Installa il pacchetto sf per le mappe GIS")
      ),

      ## --- Tab 6: Report ---
      tabItem(tabName = "report",
        fluidRow(
          box(title = "Report interpretativo", status = "primary", solidHeader = TRUE, width = 4,
              selectInput("report_lang", "Lingua", choices = c("Italiano" = "it", "English" = "en")),
              actionButton("gen_report", "Genera report", icon = icon("file-text")),
              hr(),
              h4("Download risultati"),
              downloadButton("dl_phases_csv", "Fasi (CSV)"),
              br(), br(),
              downloadButton("dl_intrusions_csv", "Intrusioni (CSV)"),
              br(), br(),
              downloadButton("dl_us_summary_csv", "Purezza US (CSV)"),
              br(), br(),
              downloadButton("dl_all_zip", "Tutto (ZIP)")
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

  # CSV
  observeEvent(input$csv_file, {
    tryCatch({
      sep <- if (input$csv_sep_semicolon) ";" else ","
      rv$raw_data <- read.csv(input$csv_file$datapath, sep = sep, stringsAsFactors = FALSE)
      showNotification(paste("Caricati", nrow(rv$raw_data), "record"), type = "message")
    }, error = function(e) showNotification(paste("Errore CSV:", e$message), type = "error"))
  })

  # SQLite
  observeEvent(input$load_sqlite, {
    req(input$sqlite_file)
    tryCatch({
      con <- DBI::dbConnect(RSQLite::SQLite(), input$sqlite_file$datapath)
      on.exit(DBI::dbDisconnect(con))
      rv$raw_data <- DBI::dbReadTable(con, input$sqlite_table)
      showNotification(paste("Caricati", nrow(rv$raw_data), "record da SQLite"), type = "message")
    }, error = function(e) showNotification(paste("Errore SQLite:", e$message), type = "error"))
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
      showNotification(paste("Caricati", nrow(rv$raw_data), "record da PostgreSQL"), type = "message")
    }, error = function(e) showNotification(paste("Errore PostgreSQL:", e$message), type = "error"))
  })

  ## --- Column Mapping ---

  output$column_mapping <- renderUI({
    req(rv$raw_data)
    cols <- names(rv$raw_data)
    none <- c("(nessuna)" = "")

    # Auto-detect matching columns
    auto <- function(target) {
      m <- match(tolower(target), tolower(cols))
      if (!is.na(m)) cols[m] else ""
    }

    tagList(
      fluidRow(
        column(2, selectInput("col_x", "x (est)", choices = c(none, cols), selected = auto("x"))),
        column(2, selectInput("col_y", "y (nord)", choices = c(none, cols), selected = auto("y"))),
        column(2, selectInput("col_z", "z (quota)", choices = c(none, cols), selected = auto("z"))),
        column(2, selectInput("col_dmin", "date_min", choices = c(none, cols), selected = auto("date_min"))),
        column(2, selectInput("col_dmax", "date_max", choices = c(none, cols), selected = auto("date_max"))),
        column(2, selectInput("col_class", "class", choices = c(none, cols), selected = auto("class")))
      ),
      fluidRow(
        column(2, selectInput("col_id", "id (opz.)", choices = c(none, cols), selected = auto("id"))),
        column(2, selectInput("col_context", "context (opz.)", choices = c(none, cols), selected = auto("context"))),
        column(2, selectInput("col_taf", "taf_score (opz.)", choices = c(none, cols), selected = auto("taf_score")))
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

      # Remove incomplete rows
      ok <- complete.cases(out[, c("x", "y", "z", "date_min", "date_max", "class")])
      out <- out[ok, ]

      rv$data <- out
      rv$fit <- NULL
      rv$compare_k <- NULL
      showNotification(paste("Dataset pronto:", nrow(out), "record,",
                             length(unique(out$class)), "classi"), type = "message")
    }, error = function(e) showNotification(paste("Errore mappatura:", e$message), type = "error"))
  })

  ## --- Data Preview ---

  output$data_summary <- renderPrint({
    req(rv$data)
    cat("Reperti:", nrow(rv$data), "\n")
    if ("context" %in% names(rv$data))
      cat("US:", length(unique(rv$data$context)), "\n")
    cat("Classi:", length(unique(rv$data$class)), "\n")
    cat("Date:", min(rv$data$date_min, na.rm = TRUE), "a",
        max(rv$data$date_max, na.rm = TRUE), "\n")
    cat("Z:", round(min(rv$data$z, na.rm = TRUE), 2), "-",
        round(max(rv$data$z, na.rm = TRUE), 2), "\n")
    if ("taf_score" %in% names(rv$data))
      cat("taf_score medio:", round(mean(rv$data$taf_score, na.rm = TRUE), 3), "\n")
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
        rv$fit <- fit_sef(rv$data, k = input$k,
                          tafonomy = taf_arg, context = ctx_arg, harris = h_arg,
                          n_init = input$n_init, em_iter = input$em_iter,
                          seed = input$seed)
        showNotification(paste("Fit completato! PDI =", round(pdi(rv$fit), 4)), type = "message")
      }, error = function(e) showNotification(paste("Errore:", e$message), type = "error"))
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
        rv$compare_k <- compare_k(rv$data, k_values = kvals,
                                   tafonomy = taf_arg, context = ctx_arg, harris = h_arg,
                                   seed = input$seed, n_init = input$ck_ninit)
        showNotification("Confronto K completato!", type = "message")
      }, error = function(e) showNotification(paste("Errore:", e$message), type = "error"))
    })
  })

  output$analysis_status <- renderPrint({
    if (!is.null(rv$fit)) {
      print(rv$fit)
      cat("\nPDI:", pdi(rv$fit), "\n")
      di <- detect_intrusions(rv$fit)
      cat("Intrusioni (prob > 0.5):", sum(di$intrusion_prob > 0.5), "/", nrow(rv$data), "\n")
    } else {
      cat("Nessuna analisi avviata. Carica i dati e premi 'Avvia analisi'.\n")
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
    DT::datatable(di, options = list(scrollX = TRUE, pageLength = 15)) |>
      DT::formatRound(columns = "intrusion_prob", digits = 4)
  })

  output$tbl_us_summary <- DT::renderDataTable({
    req(rv$fit)
    tryCatch({
      DT::datatable(us_summary_table(rv$fit),
                    options = list(scrollX = TRUE, pageLength = 20)) |>
        DT::formatPercentage(columns = "purity", digits = 1) |>
        DT::formatRound(columns = c("mean_entropy", "mean_energy", "mean_local_sei"), digits = 4)
    }, error = function(e) DT::datatable(data.frame(note = "Context non disponibile")))
  })

  output$tbl_transition <- DT::renderDataTable({
    req(rv$fit)
    tm <- phase_transition_matrix(rv$fit)
    DT::datatable(as.data.frame.matrix(tm), options = list(scrollX = TRUE))
  })

  output$txt_print_fit <- renderPrint({ req(rv$fit); print(rv$fit) })
  output$txt_summary_fit <- renderPrint({ req(rv$fit); summary(rv$fit) })

  ## --- Maps ---

  observeEvent(input$geom_file, {
    req(has_sf)
    tryCatch({
      files <- input$geom_file
      # For shapefiles: copy all files to temp dir preserving extensions
      tmpdir <- tempfile()
      dir.create(tmpdir)
      for (i in seq_len(nrow(files))) {
        ext <- tools::file_ext(files$name[i])
        base <- tools::file_path_sans_ext(files$name[1])  # use first file's base name
        newpath <- file.path(tmpdir, paste0(base, ".", ext))
        file.copy(files$datapath[i], newpath)
      }
      # Find the main file (.gpkg, .geojson, or .shp)
      main_file <- list.files(tmpdir, pattern = "\\.(gpkg|geojson|shp)$", full.names = TRUE)[1]
      geom <- sf::st_read(main_file, quiet = TRUE)
      geom <- sf::st_make_valid(geom)

      # Rename context column
      ctx_col <- input$geom_context_col
      if (ctx_col %in% names(geom) && ctx_col != "context") {
        geom$context <- as.character(geom[[ctx_col]])
      } else if (!"context" %in% names(geom)) {
        # Try common names
        for (cn in c("us", "US", "us_s", "context", "nome")) {
          if (cn %in% names(geom)) { geom$context <- as.character(geom[[cn]]); break }
        }
      }

      # Dissolve per US
      us_list <- unique(geom$context)
      dissolved <- do.call(rbind, lapply(us_list, function(u) {
        sub <- geom[geom$context == u, ]
        sf::st_sf(context = u, geom = sf::st_union(sub))
      }))
      rv$geometries <- sf::st_cast(dissolved, "MULTIPOLYGON")

      showNotification(paste("Geometrie caricate:", nrow(rv$geometries), "US"), type = "message")
    }, error = function(e) showNotification(paste("Errore geometrie:", e$message), type = "error"))
  })

  observeEvent(input$render_map, {
    req(rv$fit, rv$geometries)
    # Trigger plot update
    output$plot_map_trigger <- renderText({ input$render_map })
  })

  sef_plot_render(input, output, session, "plot_map", function() {
    req(rv$fit, rv$geometries)
    input$render_map  # dependency
    gg_map(rv$fit, geometries = rv$geometries, layer = input$map_layer)
  })

  ## --- Report ---

  report_rv <- reactiveVal(NULL)

  observeEvent(input$gen_report, {
    req(rv$fit)
    tryCatch({
      txt <- capture.output(report_sef(rv$fit, lang = input$report_lang))
      report_rv(paste(txt, collapse = "\n"))
      showNotification("Report generato!", type = "message")
    }, error = function(e) showNotification(paste("Errore:", e$message), type = "error"))
  })

  output$report_text <- renderPrint({
    if (!is.null(report_rv())) cat(report_rv())
    else cat("Premi 'Genera report' dopo aver eseguito l'analisi.")
  })

  ## --- Downloads ---

  output$dl_phases_csv <- downloadHandler(
    filename = function() paste0("palimpsestr_fasi_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$fit)
      write.csv(as_phase_table(rv$fit), file, row.names = FALSE)
    }
  )

  output$dl_intrusions_csv <- downloadHandler(
    filename = function() paste0("palimpsestr_intrusioni_", Sys.Date(), ".csv"),
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
        error = function(e) write.csv(data.frame(note = "Context non disponibile"), file, row.names = FALSE)
      )
    }
  )

  output$dl_all_zip <- downloadHandler(
    filename = function() paste0("palimpsestr_risultati_", Sys.Date(), ".zip"),
    content = function(file) {
      req(rv$fit)
      tmpdir <- tempfile()
      dir.create(tmpdir)
      export_results(rv$fit, dir = tmpdir, prefix = "palimpsestr")
      # Add report
      tryCatch({
        txt <- capture.output(report_sef(rv$fit, lang = "it"))
        writeLines(txt, file.path(tmpdir, "report_it.md"))
        txt_en <- capture.output(report_sef(rv$fit, lang = "en"))
        writeLines(txt_en, file.path(tmpdir, "report_en.md"))
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
