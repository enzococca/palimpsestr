# Model-output diagnostic plots (v0.18.0). Each function takes a sef_fit and
# returns a plain ggplot, reusing the internal helpers from R/gg_plots.R.

#' Per-phase class composition plot
#'
#' Visualises the estimated per-phase categorical class profile
#' (\code{object$cat_prob}) from a multinomial-class fit.
#'
#' @param object A \code{sef_fit} fitted with \code{class_model = "multinomial"}.
#' @param type Either \code{"heatmap"} (default; phase x class tiles filled by
#'   probability) or \code{"bar"} (stacked bars per phase).
#' @param top_n For \code{type = "bar"}, show only the \code{top_n} classes by
#'   mean probability and collapse the rest into \code{"other"} (default: all).
#' @return A ggplot object.
#' @seealso \code{\link{phase_composition}}, \code{\link{fit_sef}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' gg_phase_composition(fit)
#' @export
gg_phase_composition <- function(object, type = c("heatmap", "bar"), top_n = NULL) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  type <- match.arg(type)
  .check_ggplot()
  if (is.null(object$cat_prob)) {
    stop("gg_phase_composition() requires a fit made with class_model = \"multinomial\" ",
         "(no per-phase class profile exists for a Gaussian-class fit).", call. = FALSE)
  }
  cp <- object$cat_prob
  k <- nrow(cp); classes <- colnames(cp)
  df <- data.frame(
    phase = factor(rep(paste0("phase", seq_len(k)), times = ncol(cp)),
                   levels = paste0("phase", seq_len(k))),
    class = factor(rep(classes, each = k), levels = classes),
    prob  = as.numeric(cp), stringsAsFactors = FALSE)

  if (type == "heatmap") {
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$class, y = .data$phase, fill = .data$prob)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        .viridis_f("P(class | phase)") +
        ggplot2::labs(title = "Per-phase class composition",
                      subtitle = "Estimated categorical class profile of each phase",
                      x = "Material class", y = "Phase", caption = .sef_caption()) +
        .theme_sef() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    )
  }

  if (!is.null(top_n) && top_n < length(classes)) {
    keep <- names(sort(colMeans(cp), decreasing = TRUE))[seq_len(top_n)]
    df$class <- as.character(df$class)
    df$class[!(df$class %in% keep)] <- "other"
    df <- stats::aggregate(prob ~ phase + class, data = df, FUN = sum)
    lev <- c(keep, if (any(df$class == "other")) "other")
    df$class <- factor(df$class, levels = lev)
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$phase, y = .data$prob, fill = .data$class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(title = "Per-phase class composition",
                  subtitle = "Estimated categorical class profile of each phase",
                  x = "Phase", y = "Probability", fill = "Class", caption = .sef_caption()) +
    .theme_sef()
}

#' Directional intrusion plot
#'
#' Diverging lollipop of the signed chronological offset (\code{chrono_gap})
#' for finds classified as residual (\emph{older-than-context}) or latent
#' (\emph{younger-than-context}) by \code{\link{detect_intrusions}}.
#'
#' @param object A \code{sef_fit} object.
#' @param top_n Maximum number of finds to show, by largest absolute gap
#'   (default: 20).
#' @return A ggplot object. When every dated find is in context (e.g. unit-tied
#'   chronology) an explanatory placeholder plot is returned.
#' @seealso \code{\link{detect_intrusions}}, \code{\link{gg_intrusions}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2, context = "context")
#' gg_direction(fit)
#' @export
gg_direction <- function(object, top_n = 20) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  di <- detect_intrusions(object)
  di <- di[!is.na(di$direction) & as.character(di$direction) != "in_context", ]

  if (nrow(di) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, colour = "grey40", size = 4.5,
          label = "All dated finds are in context\n(directional diagnostics require per-find dating)") +
        ggplot2::labs(title = "Directional intrusions",
                      subtitle = "No residual or latent finds detected", caption = .sef_caption()) +
        .theme_sef() +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())
    )
  }

  di <- di[order(-abs(di$chrono_gap)), ]
  if (nrow(di) > top_n) di <- di[seq_len(top_n), ]
  di$id <- factor(di$id, levels = di$id[order(di$chrono_gap)])
  di$direction <- factor(as.character(di$direction),
    levels = c("older_than_context", "younger_than_context"),
    labels = c("older (residual)", "younger (latent feature)"))

  ggplot2::ggplot(di, ggplot2::aes(x = .data$chrono_gap, y = .data$id, colour = .data$direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$chrono_gap, yend = .data$id), linewidth = 0.8) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(name = "Direction", drop = FALSE,
      values = c("older (residual)" = "#0072B2", "younger (latent feature)" = "#D55E00")) +
    ggplot2::labs(title = "Directional intrusions",
                  subtitle = "Signed chronological offset from the stratigraphic-unit envelope",
                  x = "Chronological gap (years)", y = "Find", caption = .sef_caption()) +
    .theme_sef()
}

#' Intrusion ranking plot
#'
#' Non-spatial ranking of the per-find intrusion probability (complements
#' \code{\link{gg_intrusions}}, which is the map). Uses the noise-component
#' posterior when the model was fitted with \code{noise = TRUE}, otherwise the
#' heuristic composite score.
#'
#' @param object A \code{sef_fit} object.
#' @param threshold Reference probability above which finds are highlighted
#'   (default: 0.5).
#' @param top_n Number of highest-scoring finds to show (default: 20).
#' @return A ggplot object.
#' @seealso \code{\link{detect_intrusions}}, \code{\link{gg_intrusions}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2, noise = TRUE)
#' gg_outliers(fit)
#' @export
gg_outliers <- function(object, threshold = 0.5, top_n = 20) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  di <- detect_intrusions(object, intrusion_threshold = threshold)
  has_noise <- !is.null(object$noise_prob)
  sub <- if (has_noise) "Noise-component posterior (model-based outlier probability)"
         else "Heuristic composite score (fit without noise = TRUE)"
  df <- data.frame(id = as.character(di$id), prob = di$intrusion_prob,
                   type = as.character(di$intrusion_type), stringsAsFactors = FALSE)
  df <- df[order(-df$prob), ]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), ]
  df$id <- factor(df$id, levels = rev(df$id))
  df$type <- factor(df$type,
    levels = c("not_flagged", "residual", "latent_feature", "outlier_in_context"))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$prob, y = .data$id, colour = .data$type)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$prob, yend = .data$id), linewidth = 0.7) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(name = "Type", drop = FALSE,
      values = c(not_flagged = "grey70", residual = "#0072B2",
                 latent_feature = "#D55E00", outlier_in_context = "#E69F00")) +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(title = "Intrusion ranking", subtitle = sub,
                  x = "Intrusion probability", y = "Find", caption = .sef_caption()) +
    .theme_sef()
}

#' Within-unit phase coherence plot
#'
#' Shows, per stratigraphic unit, whether the finds were assigned to a single
#' phase (coherent) or split across several (split). A coherent assignment is
#' expected because a unit is one depositional event. An optional second fit
#' (e.g. a \code{class_model = "gaussian"} fit) is shown side by side.
#'
#' @param object A \code{sef_fit} fitted with a \code{context} column.
#' @param compare Optional second \code{sef_fit} over the same data, shown
#'   alongside \code{object}.
#' @return A ggplot object.
#' @seealso \code{\link{fit_sef}}, \code{\link{us_summary_table}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3, context = "context")
#' gg_unit_coherence(fit)
#' @export
gg_unit_coherence <- function(object, compare = NULL) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  if (is.null(object$context)) {
    stop("gg_unit_coherence() requires a fit made with a 'context' column.", call. = FALSE)
  }
  build <- function(fit, label) {
    ctx <- fit$data[[fit$context]]
    up <- tapply(fit$phase, ctx, function(p) length(unique(p)))
    data.frame(context = names(up), model = label,
               coherent = ifelse(up == 1, "coherent", "split"), stringsAsFactors = FALSE)
  }
  df <- build(object, "this fit")
  if (!is.null(compare)) {
    if (!inherits(compare, "sef_fit")) stop("compare must be a sef_fit", call. = FALSE)
    df <- rbind(df, build(compare, "comparison"))
  }
  n_split <- tapply(df$coherent == "split", df$model, sum)
  n_tot   <- tapply(df$context, df$model, length)
  sub <- paste(sprintf("%s: %d of %d units split", names(n_split),
                       as.integer(n_split), as.integer(n_tot)), collapse = "   |   ")
  agg <- as.data.frame(table(model = df$model, coherent = df$coherent))

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$model, y = .data$Freq, fill = .data$coherent)) +
    ggplot2::geom_col(position = "stack", width = 0.6) +
    ggplot2::scale_fill_manual(name = NULL,
      values = c(coherent = "#009E73", split = "#D55E00")) +
    ggplot2::labs(title = "Within-unit phase coherence", subtitle = sub,
                  x = NULL, y = "Stratigraphic units", caption = .sef_caption()) +
    .theme_sef()
}
