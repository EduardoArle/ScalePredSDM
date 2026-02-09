#' Estimate scale-of-effect (SoE) correlations for raster predictors
#'
#' Computes correlations between a response variable measured at point locations
#' and predictor values extracted at the point scale ("point") and at multiple
#' spatial buffer scales.
#'
#' Aggregation functions must be supplied explicitly via `fun_by_pred`,
#' a named character vector mapping each predictor to an aggregation rule
#' (e.g. `"mean"` for continuous predictors, `"sum"` for compositional or binary predictors).
#'
#' At the `"point"` scale, no spatial aggregation is performed; however, the
#' aggregation function is still recorded for consistency.
#'
#' @param predictors SpatRaster. Predictor stack with named layers.
#' @param points_sf sf object with POINT geometry.
#' @param response_col Character. Column in `points_sf` containing the response variable.
#' @param scales_m Numeric vector. Buffer radii (map units, typically metres).
#' @param method Character. Correlation method passed to `stats::cor` (default `"spearman"`).
#' @param include_point Logical. Include the point scale (`"point"`) in the analysis.
#' @param fun_by_pred Named character vector mapping predictors to aggregation
#'   functions (`"mean"` or `"sum"`). Must include all predictors.
#' @param return_plot Logical. If TRUE, returns a list including a ggplot object.
#' @param plot_highlight_best Logical. Highlight the best (max |correlation|) scale per predictor.
#' @param quiet Logical. If FALSE, prints progress messages.
#'
#' @return
#' If `return_plot = FALSE`:
#'   A data.frame with columns `predictor`, `scale`, `agg_fun`,
#'   `correlation`, and `n_complete`.
#'
#' If `return_plot = TRUE`:
#'   A list with elements `soe_df`, `best_scales`, and `plot`.
#'
#' @export
estimate_SoE_correlations <- function(
    predictors,
    points_sf,
    response_col,
    scales_m,
    method = "spearman",
    include_point = TRUE,
    fun_by_pred,
    return_plot = FALSE,
    plot_highlight_best = TRUE,
    quiet = FALSE
) {

  # ---- basic validation ----
  if (!inherits(predictors, "SpatRaster")) {
    stop("predictors must be a terra SpatRaster.")
  }

  pred_names <- names(predictors)
  if (is.null(pred_names) || any(pred_names == "")) {
    stop("predictors must have named layers.")
  }

  if (!inherits(points_sf, "sf")) {
    stop("points_sf must be an sf object.")
  }

  if (!all(sf::st_geometry_type(points_sf) %in% c("POINT", "MULTIPOINT"))) {
    stop("points_sf must contain POINT geometries.")
  }

  if (!(response_col %in% names(points_sf))) {
    stop("response_col not found in points_sf.")
  }

  if (!is.numeric(scales_m) || length(scales_m) < 1) {
    stop("scales_m must be a numeric vector of length >= 1.")
  }

  if (any(!is.finite(scales_m)) || any(scales_m < 0)) {
    stop("scales_m must contain finite, non-negative values.")
  }

  # ---- validate aggregation functions ----
  if (missing(fun_by_pred) || is.null(fun_by_pred)) {
    stop(
      "You must supply 'fun_by_pred', a named character vector mapping ",
      "predictor names to aggregation functions (e.g., c(depth='mean', habitat='sum'))."
    )
  }

  if (!is.character(fun_by_pred) || is.null(names(fun_by_pred))) {
    stop("'fun_by_pred' must be a named character vector.")
  }

  allowed_funs <- c("mean", "sum")
  bad_funs <- setdiff(unique(fun_by_pred), allowed_funs)
  if (length(bad_funs) > 0) {
    stop("Unsupported aggregation functions: ", paste(bad_funs, collapse = ", "))
  }

  missing_fun <- setdiff(pred_names, names(fun_by_pred))
  if (length(missing_fun) > 0) {
    stop(
      "Missing aggregation functions for predictors: ",
      paste(missing_fun, collapse = ", ")
    )
  }

  # ---- clean response ----
  y <- points_sf[[response_col]]
  ok <- is.finite(y)

  if (sum(ok) < 3) {
    stop("Need at least 3 finite response values to compute correlations.")
  }

  pts <- points_sf[ok, , drop = FALSE]
  y   <- y[ok]
  pts_vect <- terra::vect(pts)

  out <- list()

  # ---- 1) point scale ----
  if (isTRUE(include_point)) {
    if (!quiet) message("Scale: point")

    for (pred in pred_names) {
      if (!quiet) message("  Predictor: ", pred)

      r <- predictors[[pred]]
      fun_name <- fun_by_pred[[pred]]

      x <- terra::extract(r, pts_vect)[, 2]
      keep <- is.finite(x) & is.finite(y)

      rho <- suppressWarnings(
        stats::cor(x[keep], y[keep], method = method, use = "complete.obs")
      )

      out[[length(out) + 1]] <- data.frame(
        predictor   = pred,
        scale       = "point",
        agg_fun     = fun_name,
        correlation = rho,
        n_complete  = sum(keep),
        stringsAsFactors = FALSE
      )
    }
  }

  # ---- 2) buffered scales ----
  for (sc in scales_m) {
    if (!quiet) message("Scale: ", sc)

    buf_sf   <- sf::st_buffer(pts, dist = sc)
    buf_vect <- terra::vect(buf_sf)

    for (pred in pred_names) {
      if (!quiet) message("  Predictor: ", pred)

      r <- predictors[[pred]]
      fun_name <- fun_by_pred[[pred]]

      x <- terra::extract(
        r,
        buf_vect,
        fun   = fun_name,
        na.rm = TRUE
      )[, 2]

      keep <- is.finite(x) & is.finite(y)

      rho <- suppressWarnings(
        stats::cor(x[keep], y[keep], method = method, use = "complete.obs")
      )

      out[[length(out) + 1]] <- data.frame(
        predictor   = pred,
        scale       = sc,
        agg_fun     = fun_name,
        correlation = rho,
        n_complete  = sum(keep),
        stringsAsFactors = FALSE
      )
    }
  }

  soe_df <- do.call(rbind, out)

  if (!isTRUE(return_plot)) {
    return(soe_df)
  }

  # ---- identify best scales ----
  best_scales <- NULL
  if (isTRUE(plot_highlight_best)) {
    best_scales <- soe_df |>
      dplyr::group_by(.data$predictor) |>
      dplyr::slice_max(abs(.data$correlation), n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  # ---- plotting ----
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for return_plot = TRUE.")
  }

  numeric_scales <- sort(unique(as.numeric(soe_df$scale[soe_df$scale != "point"])))
  scale_levels <- c("point", as.character(numeric_scales))

  soe_df$scale_ord <- factor(
    soe_df$scale,
    levels = scale_levels,
    ordered = TRUE
  )
  # ensure best_scales has the same ordered scale factor
  if (!is.null(best_scales)) {
    best_scales$scale_ord <- factor(
      best_scales$scale,
      levels = levels(soe_df$scale_ord),
      ordered = TRUE
    )
  }

  p <- ggplot2::ggplot(
    soe_df,
    ggplot2::aes(
      x = scale_ord,
      y = correlation,
      colour = predictor,
      group = predictor
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::labs(
      x = "Scale of effect",
      y = "Correlation",
      colour = "Predictor",
      title = "Scale-of-effect correlations"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (!is.null(best_scales)) {
    p <- p +
      ggplot2::geom_point(
        data = best_scales,
        ggplot2::aes(x = scale_ord, y = correlation),
        inherit.aes = FALSE,
        shape = 8,
        size = 4,
        stroke = 1.2
      )
  }

  list(
    soe_df = soe_df,
    best_scales = best_scales,
    plot = p
  )
}
