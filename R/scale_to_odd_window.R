#' Convert a target scale (map units) to an odd focal window size (cells)
#'
#' @param r A single-layer SpatRaster
#' @param scale_m Numeric. Spatial scale in metres.
#' @param scale_type Character. Interpretation of `scale`:
#'   - "diameter": window size ≈ scale / resolution
#'   - "radius":   window size ≈ (2 * scale) / resolution
#'
#' @return Integer odd window size (>= 1)
#' @keywords internal
scale_to_odd_window <- function(r, scale_m, scale_type = c("diameter", "radius")) {

  scale_type <- match.arg(scale_type)

  if (!is.numeric(scale_m) || !is.finite(scale_m) || scale_m <= 0) {
    stop("scale_m must be a positive finite number.")
  }

  res_m <- terra::res(r)[1]

  if (scale_type == "radius") {
    w <- round((2 * scale_m) / res_m)
  } else {
    w <- round(scale_m / res_m)
  }

  # enforce odd
  if (w %% 2 == 0) w <- w + 1

  # IMPORTANT: minimum meaningful focal window is 3x3
  if (w < 3) {
    stop(
      "Selected scale (", scale_m, " m) is smaller than raster resolution (",
      round(res_m, 2), " m). ",
      "Treat this scale as 'point'."
    )
  }

  w
}
