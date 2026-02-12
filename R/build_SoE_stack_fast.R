#' Build a scale-of-effect (SoE) raster stack from original predictors
#'
#'
#' Given (i) an original predictor raster stack and (ii) a table of selected
#' scales-of-effect, this function produces a new raster stack in which each
#' predictor is represented at its selected spatial scale.
#'
#' Spatial aggregation is applied using `terra::focal()` for numeric scales.
#' If `scale_m == "point"`, the original predictor raster is returned unchanged.
#'
#' Aggregation functions must be supplied explicitly via `fun_by_pred`,
#' a named character vector mapping each predictor to an aggregation rule
#' (e.g. `"mean"` for continuous variables, `"sum"` for compositional or binary variables).
#'
#' @param predictors SpatRaster. Original predictor layers. Must be named.
#' @param SoE_table data.frame or tibble with columns:
#'   - `predictor`: name of the predictor (must match layer names in `predictors`)
#'   - `scale_m`: numeric scale (map units, typically metres) or `"point"`
#' @param scale_type Character. Interpretation of `scale_m`: `"diameter"` or `"radius"`.
#' @param fun_by_pred Named character vector mapping predictors to aggregation
#'   functions (`"mean"` or `"sum"`). Must include all predictors.
#' @param na.rm Logical. Passed to `terra::focal()`.
#' @param expand Logical. Passed to `terra::focal()`.
#' @param filename Optional filename to write the resulting raster stack.
#' @param overwrite Logical. Overwrite existing file if `filename` is supplied.
#' @param quiet Logical. If FALSE, prints progress messages.
#'
#' @return SpatRaster with one layer per predictor, named by predictor.
#' @export
build_SoE_stack_fast <- function(
    predictors,
    SoE_table,
    scale_type = c("diameter", "radius"),
    fun_by_pred,
    na.rm = TRUE,
    expand = FALSE,
    filename = NULL,
    overwrite = TRUE,
    quiet = FALSE
) {

  # ---- validation ----
  if (!inherits(predictors, "SpatRaster")) {
    stop("predictors must be a terra SpatRaster.")
  }

  pred_names <- names(predictors)
  if (is.null(pred_names) || any(pred_names == "")) {
    stop("predictors must have named layers.")
  }

  if (!is.data.frame(SoE_table)) {
    stop("SoE_table must be a data.frame or tibble.")
  }

  if (!all(c("predictor", "scale_m") %in% names(SoE_table))) {
    stop("SoE_table must contain columns: 'predictor' and 'scale_m'.")
  }

  scale_type <- match.arg(scale_type)

  SoE_table <- as.data.frame(SoE_table)
  SoE_table$predictor <- as.character(SoE_table$predictor)

  missing_preds <- setdiff(SoE_table$predictor, pred_names)
  if (length(missing_preds) > 0) {
    stop(
      "These predictors in SoE_table are missing from predictors: ",
      paste(missing_preds, collapse = ", ")
    )
  }

  # ---- aggregation functions ----
  if (missing(fun_by_pred) || is.null(fun_by_pred)) {
    stop(
      "You must supply 'fun_by_pred', a named character vector mapping ",
      "predictor names to aggregation functions."
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

  # ---- build SoE layers ----
  out_layers <- vector("list", nrow(SoE_table))

  for (i in seq_len(nrow(SoE_table))) {

    pred <- SoE_table$predictor[i]
    sc   <- SoE_table$scale_m[i]
    r    <- predictors[[pred]]

    # --- CASE 1: point scale ---
    if (is.character(sc) && sc == "point") {

      if (!quiet) {
        message("Predictor: ", pred, " | scale: point (native resolution)")
      }

      r_f <- r
      names(r_f) <- pred
      out_layers[[i]] <- r_f
      next
    }

    # --- CASE 2: numeric scale ---
    sc_num <- as.numeric(sc)
    if (!is.finite(sc_num) || sc_num <= 0) {
      stop(
        "scale_m must be a positive number or 'point'. ",
        "Problem with predictor '", pred, "'."
      )
    }

    w <- scale_to_odd_window(r, sc_num, scale_type = scale_type)
    fun_name <- fun_by_pred[[pred]]

    if (!quiet) {
      message(
        "Predictor: ", pred,
        " | scale: ", sc_num,
        " | window: ", w, "x", w,
        " | fun: ", fun_name
      )
    }

    r_f <- terra::focal(
      r,
      w      = matrix(1, nrow = w, ncol = w),
      fun    = fun_name,
      na.rm  = na.rm,
      expand = expand
    )

    names(r_f) <- pred
    out_layers[[i]] <- r_f
  }

  soe_stack <- terra::rast(out_layers)
  names(soe_stack) <- SoE_table$predictor

  if (!is.null(filename)) {
    soe_stack <- terra::writeRaster(
      soe_stack,
      filename,
      overwrite = overwrite
    )
  }

  soe_stack
}
