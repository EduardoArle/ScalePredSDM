# =============================================================================
# Anonymous positional-error benchmark for scale-of-effect predictor stacks
# =============================================================================
# Purpose
#   This script demonstrates the complete benchmark used to compare
#   scale-of-effect (SoE) predictor stacks with conventional single-scale
#   predictor stacks under simulated positional error.
#
#   The original high-resolution raster predictors are not included with the
#   anonymous submission. By default, the script therefore runs in DEMO mode,
#   generating a small synthetic raster stack and synthetic abundance data.
#   The benchmark structure is the same as the manuscript workflow:
#     1. select predictor-specific SoE scales from abundance data,
#     2. build SoE raster stacks using ScalePredSDM,
#     3. perturb occurrence coordinates across positional-error distances,
#     4. fit matched single-scale and SoE MaxEnt models, and
#     5. summarise paired differences in AUC and TSS.
#
#   For the real analysis, set CONFIG$data_mode = "real" and provide local
#   rasters and survey points in the anonymous data folders specified below.
#
# Important for anonymous review
#   No local user paths, project names or author-identifying paths are included.
#   The script assumes that the ScalePredSDM package is either already installed
#   or supplied separately as an anonymous source archive.
# =============================================================================

options(stringsAsFactors = FALSE)

# =============================================================================
# 0) CONFIG
# =============================================================================
CONFIG <- list(

  # ---- Data mode ----
  # "demo" generates synthetic input data so the code can run without unpublished
  # raster files. "real" loads local anonymous input files from data_dir.
  data_mode = "demo",

  # ---- Anonymous folders ----
  data_dir = "data",
  out_dir  = "output_anonymous_benchmark_1m",

  # ---- Optional local package archive for anonymous review ----
  # Leave NULL if ScalePredSDM is already installed. If supplying an anonymous
  # package source archive, place it next to this script and set the filename here.
  scalepredsdm_tarball = NULL,  # e.g. "ScalePredSDM_0.1.0.tar.gz"

  # ---- Real-data inputs, used only when data_mode = "real" ----
  abundance_file  = file.path("data", "survey_points.geojson"),
  abundance_layer = NULL,
  single_pred_dir = file.path("data", "single_scale_predictors"),
  single_pred_pat = "\\.(tif|asc|grd)$",

  # ---- Abundance table layout ----
  # Wide format: one species column per survey point. Long format is also accepted
  # if the file already contains species_col and value_col.
  species_cols_wide = c("species_A", "species_B", "species_C", "species_D"),
  id_col      = "site_id",
  species_col = "species",
  value_col   = "value",

  # ---- Species to run ----
  species_list = c("species_A", "species_B", "species_C", "species_D"),

  # ---- Fixed MaxEnt settings ----
  # Use the manuscript settings in real-data mode. Demo labels are generic.
  fc_rm_table = data.frame(
    species = c("species_A", "species_B", "species_C", "species_D"),
    fc      = c("LQ",        "LQ",        "LQ",        "L"),
    rm      = c(2,           2,           3,           1),
    stringsAsFactors = FALSE
  ),

  # ---- Positional-error design ----
  # Demo defaults are intentionally small. For the full manuscript run use:
  # pe_dists_m = c(0, 2, 5, 10, 25, 50, 100, 200, 300), reps = 20
  pe_dists_m = c(0, 10, 25),
  reps       = 2,

  # ---- Coordinate system ----
  # Use a projected CRS in metres. EPSG:32636 was used in the real analysis.
  crs_epsg = 32636,

  # ---- Background and evaluation ----
  # Demo defaults are small. For the full manuscript run use n_bg = 10000.
  n_bg             = 500,
  bg_seed          = 9001,
  kfold_k          = 3,
  pres_threshold_n = 0,   # 0 forces k-fold evaluation in demo mode

  # ---- ScalePredSDM SoE settings ----
  # Demo defaults are small. For the full manuscript run use:
  # soe_scales_m = c(2, 5, 10, 25, 50, 100, 200, 300)
  soe_scales_m         = c(2, 10, 25, 50),
  soe_method           = "spearman",
  include_point        = TRUE,
  scale_type_for_stack = "radius",
  overwrite_soe_stack  = TRUE,

  # Current ScalePredSDM stack builder supports mean/sum aggregation.
  # Predictors not listed here default to "mean".
  agg_fun_by_pred = c(
    depth    = "mean",
    slope    = "mean",
    sos      = "mean",
    eastness = "mean",
    edge     = "sum",
    complex  = "sum"
  ),

  # ---- Predictors ----
  # NULL uses all raster layers.
  keep_predictors = NULL,

  # ---- Misc ----
  seed = 42,
  write_demo_inputs = TRUE
)

# =============================================================================
# 1) Packages and ScalePredSDM checks
# =============================================================================
msg <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))
  flush.console()
}

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(maxnet)
  library(pROC)
})

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || grepl("^/", path) || grepl("^~", path)
}

resolve_path <- function(path, base = NULL) {
  if (is.null(path) || is.na(path)) return(path)
  if (is_absolute_path(path) || is.null(base)) return(path.expand(path))
  file.path(base, path)
}

ensure_scalepredsdm <- function(cfg) {
  if (!requireNamespace("ScalePredSDM", quietly = TRUE)) {
    tarball <- cfg$scalepredsdm_tarball
    if (!is.null(tarball) && file.exists(tarball)) {
      msg("Installing ScalePredSDM from anonymous local archive: %s", tarball)
      install.packages(tarball, repos = NULL, type = "source")
    }
  }

  if (!requireNamespace("ScalePredSDM", quietly = TRUE)) {
    stop(
      "ScalePredSDM is required. For anonymous review, install the package from ",
      "the supplied anonymous source archive, or set CONFIG$scalepredsdm_tarball."
    )
  }

  needed <- c("estimate_SoE_correlations", "build_SoE_stack_fast")
  missing <- setdiff(needed, getNamespaceExports("ScalePredSDM"))
  if (length(missing) > 0) {
    stop("ScalePredSDM is missing exported function(s): ", paste(missing, collapse = ", "))
  }

  desc <- utils::packageDescription("ScalePredSDM")
  msg("Loaded ScalePredSDM %s", desc$Version)
  invisible(TRUE)
}

ensure_scalepredsdm(CONFIG)
set.seed(CONFIG$seed)
dir.create(CONFIG$out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2) Generic helpers
# =============================================================================
call_with_supported_args <- function(fun, args) {
  fml <- names(formals(fun))
  if ("..." %in% fml) return(do.call(fun, args))
  do.call(fun, args[names(args) %in% fml])
}

make_fun_by_pred <- function(pred_names, cfg) {
  funs <- stats::setNames(rep("mean", length(pred_names)), pred_names)
  user_map <- cfg$agg_fun_by_pred
  if (!is.null(user_map)) {
    matched <- intersect(names(user_map), pred_names)
    funs[matched] <- user_map[matched]
  }
  bad <- setdiff(unique(funs), c("mean", "sum"))
  if (length(bad) > 0) stop("Unsupported aggregation function(s): ", paste(bad, collapse = ", "))
  funs
}

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

safe_auc <- function(obs, pred) {
  if (length(unique(obs)) < 2 || length(unique(pred)) < 2) return(NA_real_)
  as.numeric(pROC::roc(obs, pred, quiet = TRUE)$auc)
}

best_tss_threshold <- function(obs, pred, thresholds = seq(0, 1, by = 0.01)) {
  out <- lapply(thresholds, function(th) {
    bin <- as.integer(pred >= th)
    tp <- sum(bin == 1 & obs == 1)
    tn <- sum(bin == 0 & obs == 0)
    fp <- sum(bin == 1 & obs == 0)
    fn <- sum(bin == 0 & obs == 1)
    sens <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    spec <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
    data.frame(threshold = th, sensitivity = sens, specificity = spec, tss = sens + spec - 1)
  }) |> bind_rows() |> filter(!is.na(tss))
  if (nrow(out) == 0) return(data.frame(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_, tss = NA_real_))
  out[which.max(out$tss), , drop = FALSE]
}

normalise_fc <- function(fc) tolower(gsub("[^A-Za-z]", "", fc))

prepare_maxnet_frame <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  x[] <- lapply(x, function(z) as.numeric(z))
  x
}

non_constant_cols <- function(x) {
  vapply(x, function(z) {
    z <- z[is.finite(z)]
    length(unique(z)) > 1 && stats::sd(z) > 1e-10
  }, logical(1))
}

fit_maxnet <- function(x, y, fc, rm) {
  # maxnet can fail during prediction if a training fold contains predictors
  # with zero variance or if formula/data columns are not handled identically.
  # This wrapper keeps only informative numeric predictors for the training set
  # and stores the exact columns needed for prediction.
  x <- prepare_maxnet_frame(x)
  y <- as.integer(y)

  keep_rows <- complete.cases(x) & !is.na(y)
  x <- x[keep_rows, , drop = FALSE]
  y <- y[keep_rows]

  if (length(unique(y)) < 2) stop("maxnet requires both presences and background/absence rows.")

  keep_cols <- non_constant_cols(x)
  if (!any(keep_cols)) stop("No non-constant predictor columns available for maxnet.")
  x <- x[, keep_cols, drop = FALSE]

  f <- maxnet::maxnet.formula(p = y, data = x, classes = normalise_fc(fc))
  mod <- maxnet::maxnet(p = y, data = x, f = f, regmult = rm)
  attr(mod, "predictor_cols") <- names(x)
  mod
}

predict_maxnet_safe <- function(model, x, type = "cloglog") {
  cols <- attr(model, "predictor_cols")
  if (is.null(cols)) cols <- names(x)
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0) {
    stop("Prediction data are missing model predictor column(s): ", paste(missing, collapse = ", "))
  }
  x <- prepare_maxnet_frame(x[, cols, drop = FALSE])
  stats::predict(model, x, type = type)
}

make_folds <- function(n, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

extract_raster_values <- function(rast, xy) {
  vals <- terra::extract(rast, xy)
  vals <- as.data.frame(vals, check.names = FALSE)

  # terra::extract() returns an ID column for vector geometries, but usually
  # not for a coordinate matrix/data.frame. Do not blindly remove the first
  # column, because that can drop the first predictor layer, e.g. depth.
  if ("ID" %in% names(vals)) {
    vals <- vals[, names(vals) != "ID", drop = FALSE]
  }
  if (ncol(vals) == terra::nlyr(rast) + 1 && tolower(names(vals)[1]) %in% c("id", "cell")) {
    vals <- vals[, -1, drop = FALSE]
  }

  # Some terra versions return generic names. Reapply raster layer names when
  # the number of extracted columns matches the number of raster layers.
  if (ncol(vals) == terra::nlyr(rast) && !all(names(rast) %in% names(vals))) {
    names(vals) <- names(rast)
  }

  missing_cols <- setdiff(names(rast), names(vals))
  if (length(missing_cols) > 0) {
    stop("Raster extraction did not return expected predictor columns: ",
         paste(missing_cols, collapse = ", "),
         ". Returned columns were: ", paste(names(vals), collapse = ", "))
  }

  vals[, names(rast), drop = FALSE]
}

extract_xy_vals <- function(rast, sf_pts) {
  xy <- sf::st_coordinates(sf_pts)
  vals <- extract_raster_values(rast, xy)
  data.frame(x = xy[, 1], y = xy[, 2], vals, check.names = FALSE)
}

# =============================================================================
# 3) Demo input generation or real input loading
# =============================================================================
make_demo_inputs <- function(cfg) {
  msg("Generating synthetic demo inputs using vignette-style predictors")

  # The demo follows the ScalePredSDM vignette idea: create patchy virtual
  # predictors from smoothed random fields, simulate abundance from a spatial
  # neighbourhood effect, then run the same SoE and positional-error workflow.
  # Here the raster extent is 1000 x 1000 m and the resolution is 1 m.
  r0 <- terra::rast(
    nrows = 1000, ncols = 1000,
    xmin = 0, xmax = 1000,
    ymin = 0, ymax = 1000,
    crs  = paste0("EPSG:", cfg$crs_epsg)
  )

  set.seed(cfg$seed)

  # 1) Patchy continuous predictor, analogous to the vignette "depth" layer.
  depth_noise <- r0
  terra::values(depth_noise) <- stats::rnorm(terra::ncell(depth_noise))
  depth_s1 <- terra::focal(depth_noise, w = matrix(1, 5, 5),  fun = "mean", na.rm = TRUE)
  depth_s2 <- terra::focal(depth_noise, w = matrix(1, 11, 11), fun = "mean", na.rm = TRUE)
  depth <- 0.6 * depth_s1 + 0.4 * depth_s2
  depth <- terra::scale(depth)
  names(depth) <- "depth"

  # 2) Rare, patchy binary habitat, analogous to the vignette habitat layer.
  hab_noise <- r0
  terra::values(hab_noise) <- stats::runif(terra::ncell(hab_noise))
  hab_smooth <- terra::focal(hab_noise, w = matrix(1, 9, 9), fun = "mean", na.rm = TRUE)
  thr <- stats::quantile(terra::values(hab_smooth), probs = 0.90, na.rm = TRUE)
  complex <- as.numeric(hab_smooth >= thr)
  names(complex) <- "complex"

  # 3) Terrain derivatives derived from the patchy continuous surface.
  slope <- terra::terrain(depth, v = "slope", unit = "degrees", neighbors = 8)
  slope <- terra::ifel(is.na(slope), 0, slope)
  names(slope) <- "slope"

  sos <- terra::terrain(slope, v = "slope", unit = "degrees", neighbors = 8)
  sos <- terra::ifel(is.na(sos), 0, sos)
  names(sos) <- "sos"

  aspect <- terra::terrain(depth, v = "aspect", unit = "radians", neighbors = 8)
  eastness <- cos(aspect)
  eastness <- terra::ifel(is.na(eastness), 0, eastness)
  names(eastness) <- "eastness"

  # 4) Edge is a local boundary around the rare habitat patches.
  # A cell is treated as edge when a 3 x 3 neighbourhood contains both habitat
  # and non-habitat cells. This provides a simple contextual seascape metric.
  complex_local_mean <- terra::focal(complex, w = matrix(1, 3, 3), fun = "mean", na.rm = TRUE)
  edge <- as.numeric(complex_local_mean > 0 & complex_local_mean < 1)
  edge <- terra::ifel(is.na(edge), 0, edge)
  names(edge) <- "edge"

  single_stack <- c(depth, slope, sos, eastness, edge, complex)
  names(single_stack) <- c("depth", "slope", "sos", "eastness", "edge", "complex")

  if (!all(terra::res(single_stack) == c(1, 1))) {
    stop("Demo raster stack was expected to have 1 m resolution, but terra::res() returned: ",
         paste(terra::res(single_stack), collapse = ", "))
  }
  msg("Synthetic demo raster resolution: %.1f x %.1f m",
      terra::res(single_stack)[1], terra::res(single_stack)[2])

  # Survey points. Kept away from the border to reduce NA extraction from
  # neighbourhood operations.
  n <- 90
  pts_df <- data.frame(
    site_id = seq_len(n),
    x = stats::runif(n, 50, 950),
    y = stats::runif(n, 50, 950)
  )
  pts_sf <- sf::st_as_sf(pts_df, coords = c("x", "y"), crs = paste0("EPSG:", cfg$crs_epsg))

  # Simulate abundance from spatial neighbourhood effects, as in the vignette.
  # Species differ in whether they respond most strongly to rare habitat,
  # habitat edge, terrain heterogeneity or the patchy continuous gradient.
  dist_complex <- terra::distance(complex == 1)
  names(dist_complex) <- "dist_complex"
  dist_edge <- terra::distance(edge == 1)
  names(dist_edge) <- "dist_edge"

  focal_complex_50 <- terra::focal(complex, w = matrix(1, 101, 101), fun = "sum", na.rm = TRUE)
  names(focal_complex_50) <- "focal_complex_50"
  focal_edge_25 <- terra::focal(edge, w = matrix(1, 51, 51), fun = "sum", na.rm = TRUE)
  names(focal_edge_25) <- "focal_edge_25"

  sim_stack <- c(single_stack, dist_complex, dist_edge, focal_complex_50, focal_edge_25)
  vals <- extract_raster_values(sim_stack, sf::st_coordinates(pts_sf))

  z_depth <- scale01(vals$depth)
  z_slope <- scale01(vals$slope)
  z_sos <- scale01(vals$sos)
  z_dist_complex <- exp(-vals$dist_complex / 80)
  z_dist_edge <- exp(-vals$dist_edge / 60)
  z_complex_context <- scale01(vals$focal_complex_50)
  z_edge_context <- scale01(vals$focal_edge_25)

  # The low intercepts allow zeros, while the neighbourhood effects ensure a
  # clear SoE signal. These are anonymous synthetic species, not manuscript data.
  mu_A <- 0.3 + 8.0 * z_dist_complex + 1.5 * z_complex_context
  mu_B <- 0.3 + 7.0 * z_dist_edge + 1.5 * z_edge_context
  mu_C <- 0.4 + 5.0 * z_complex_context + 2.0 * z_sos
  mu_D <- 0.4 + 3.0 * scale01(-z_depth) + 2.5 * z_slope

  pts_wide <- pts_sf |>
    mutate(
      species_A = stats::rpois(n, mu_A),
      species_B = stats::rpois(n, mu_B),
      species_C = stats::rpois(n, mu_C),
      species_D = stats::rpois(n, mu_D)
    )

  pts_long <- pts_wide |>
    tidyr::pivot_longer(cols = all_of(cfg$species_cols_wide),
                        names_to = "species", values_to = "value")

  if (isTRUE(cfg$write_demo_inputs)) {
    demo_dir <- file.path(cfg$out_dir, "demo_inputs")
    dir.create(demo_dir, recursive = TRUE, showWarnings = FALSE)
    sf::st_write(pts_wide, file.path(demo_dir, "demo_survey_points.geojson"), delete_dsn = TRUE, quiet = TRUE)
    terra::writeRaster(single_stack, file.path(demo_dir, "demo_single_scale_predictors.tif"), overwrite = TRUE)
  }

  list(pts_long = pts_long, single_stack = single_stack)
}

safe_read_sf <- function(path, layer = NULL) {
  if (is.null(layer)) sf::st_read(path, quiet = TRUE) else sf::st_read(path, layer = layer, quiet = TRUE)
}

load_abundance_long <- function(cfg) {
  path <- resolve_path(cfg$abundance_file, NULL)
  if (!file.exists(path)) stop("Abundance file not found: ", path)
  pts <- safe_read_sf(path, layer = cfg$abundance_layer)
  if (!cfg$id_col %in% names(pts)) pts[[cfg$id_col]] <- seq_len(nrow(pts))

  if (cfg$species_col %in% names(pts) && cfg$value_col %in% names(pts)) {
    pts |>
      select(all_of(c(cfg$id_col, cfg$species_col, cfg$value_col)), geometry) |>
      rename(species = all_of(cfg$species_col), value = all_of(cfg$value_col))
  } else {
    missing_cols <- setdiff(cfg$species_cols_wide, names(pts))
    if (length(missing_cols) > 0) stop("Missing species columns: ", paste(missing_cols, collapse = ", "))
    pts |>
      select(all_of(cfg$id_col), all_of(cfg$species_cols_wide), geometry) |>
      pivot_longer(cols = all_of(cfg$species_cols_wide), names_to = "species", values_to = "value")
  }
}

load_single_scale_predictors <- function(cfg) {
  pred_dir <- resolve_path(cfg$single_pred_dir, NULL)
  if (!dir.exists(pred_dir)) stop("Predictor directory not found: ", pred_dir)
  files <- list.files(pred_dir, pattern = cfg$single_pred_pat, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No raster predictors found in: ", pred_dir)
  r <- terra::rast(files)
  names(r) <- make.names(names(r), unique = TRUE)

  if (!is.null(cfg$keep_predictors)) {
    missing <- setdiff(cfg$keep_predictors, names(r))
    if (length(missing) > 0) stop("Predictors not found: ", paste(missing, collapse = ", "))
    r <- r[[cfg$keep_predictors]]
  }

  if (is.na(terra::crs(r))) terra::crs(r) <- paste0("EPSG:", cfg$crs_epsg)
  r
}

load_inputs <- function(cfg) {
  if (identical(cfg$data_mode, "demo")) return(make_demo_inputs(cfg))
  if (identical(cfg$data_mode, "real")) {
    return(list(
      pts_long = load_abundance_long(cfg),
      single_stack = load_single_scale_predictors(cfg)
    ))
  }
  stop("CONFIG$data_mode must be 'demo' or 'real'.")
}

# =============================================================================
# 4) ScalePredSDM SoE selection and raster generation
# =============================================================================
make_best_scales_for_stack <- function(soe_result) {
  if (is.list(soe_result) && !is.null(soe_result$best_scales)) {
    best <- as.data.frame(soe_result$best_scales)
  } else if (is.data.frame(soe_result)) {
    best <- soe_result |>
      group_by(.data$predictor) |>
      slice_max(abs(.data$correlation), n = 1, with_ties = FALSE) |>
      ungroup() |>
      as.data.frame()
  } else {
    stop("Could not identify the best-scales table returned by ScalePredSDM.")
  }

  if (!"predictor" %in% names(best)) stop("best_scales table lacks a predictor column.")
  if (!"scale" %in% names(best)) {
    if ("scale_m" %in% names(best)) best$scale <- best$scale_m
    if ("radius_m" %in% names(best)) best$scale <- best$radius_m
  }
  if (!"scale" %in% names(best)) stop("best_scales table lacks a scale column.")

  out <- data.frame(
    predictor = as.character(best$predictor),
    scale_m = as.character(best$scale),
    stringsAsFactors = FALSE
  )
  # Avoid ifelse() coercion warnings when ScalePredSDM returns a mix of
  # numeric scales and the character value "point".
  out$scale_m <- vapply(out$scale_m, function(s) {
    if (is.na(s) || tolower(s) == "point") return("point")
    s_num <- suppressWarnings(as.numeric(s))
    if (is.na(s_num)) return(as.character(s))
    as.character(s_num)
  }, character(1))
  out
}

run_scalepredsdm_soe <- function(single_stack, survey_sf, cfg, species) {
  fun_by_pred <- make_fun_by_pred(names(single_stack), cfg)
  survey_sf <- sf::st_transform(survey_sf, cfg$crs_epsg)

  msg("ScalePredSDM SoE selection | species=%s | survey rows=%d", species, nrow(survey_sf))

  soe_fun <- getExportedValue("ScalePredSDM", "estimate_SoE_correlations")
  soe_args <- list(
    predictors = single_stack,
    points_sf = survey_sf,
    response_col = "value",
    scales_m = cfg$soe_scales_m,
    method = cfg$soe_method,
    include_point = cfg$include_point,
    fun_by_pred = fun_by_pred,
    return_plot = TRUE,
    plot_highlight_best = TRUE,
    quiet = FALSE
  )
  soe_result <- call_with_supported_args(soe_fun, soe_args)

  soe_dir <- file.path(cfg$out_dir, "scale_of_effect")
  dir.create(soe_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.list(soe_result) && !is.null(soe_result$soe_df)) {
    write.csv(soe_result$soe_df, file.path(soe_dir, sprintf("soe_correlations_%s.csv", species)), row.names = FALSE)
  }
  if (is.list(soe_result) && !is.null(soe_result$best_scales)) {
    write.csv(soe_result$best_scales, file.path(soe_dir, sprintf("best_scales_raw_%s.csv", species)), row.names = FALSE)
  }

  SoE_table <- make_best_scales_for_stack(soe_result)
  write.csv(SoE_table, file.path(soe_dir, sprintf("best_scales_for_stack_%s.csv", species)), row.names = FALSE)

  stack_file <- file.path(soe_dir, sprintf("soe_stack_%s.tif", species))
  stack_fun <- getExportedValue("ScalePredSDM", "build_SoE_stack_fast")
  stack_args <- list(
    predictors = single_stack,
    SoE_table = SoE_table,
    scale_type = cfg$scale_type_for_stack,
    fun_by_pred = fun_by_pred,
    filename = stack_file,
    overwrite = cfg$overwrite_soe_stack,
    quiet = FALSE
  )
  soe_stack <- call_with_supported_args(stack_fun, stack_args)

  list(soe_result = soe_result, SoE_table = SoE_table, soe_stack = soe_stack)
}

# =============================================================================
# 5) Positional-error benchmark helpers
# =============================================================================
prepare_background_points <- function(cfg, mask_rast, species) {
  out_dir <- file.path(cfg$out_dir, "background")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bg_file <- file.path(out_dir, sprintf("background_%s.csv", species))

  if (file.exists(bg_file)) {
    return(read.csv(bg_file) |>
             transmute(x = as.numeric(x), y = as.numeric(y)) |>
             filter(is.finite(x), is.finite(y)))
  }

  set.seed(cfg$bg_seed)
  pts <- terra::spatSample(mask_rast, size = cfg$n_bg, method = "random",
                           na.rm = TRUE, as.points = TRUE, warn = 0)
  xy <- terra::crds(pts)
  bg <- data.frame(x = xy[, 1], y = xy[, 2])
  write.csv(bg, bg_file, row.names = FALSE)
  bg
}

make_pe_geometry_list <- function(occ_sf, dists_m, crs_epsg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  occ_utm <- sf::st_transform(occ_sf, crs_epsg)
  if (!all(sf::st_geometry_type(occ_utm) %in% "POINT")) stop("Occurrence geometries must be POINT geometries.")

  n <- nrow(occ_utm)
  coords <- sf::st_coordinates(occ_utm)
  angles <- runif(n, 0, 2 * pi)

  out <- list()
  for (d in dists_m) {
    if (d == 0) {
      tmp <- occ_utm
    } else {
      new_pts <- sf::st_as_sf(
        data.frame(x = coords[, 1] + cos(angles) * d,
                   y = coords[, 2] + sin(angles) * d),
        coords = c("x", "y"), crs = sf::st_crs(occ_utm)
      )
      tmp <- occ_utm
      sf::st_geometry(tmp) <- sf::st_geometry(new_pts)
    }
    out[[sprintf("pe_%dm", d)]] <- tmp
  }
  out
}

eval_cv_maxnet <- function(x, y, fc, rm, method = c("auto", "kfold", "jackknife"),
                           k = 5, pres_threshold_n = 50, seed = NULL) {
  method <- match.arg(method)
  pres_idx <- which(y == 1)
  bg_idx <- which(y == 0)

  if (method == "auto") {
    method <- if (length(pres_idx) < pres_threshold_n) "jackknife" else "kfold"
  }

  if (method == "kfold") {
    k <- min(k, length(pres_idx), length(bg_idx))
    if (k < 2) stop("Too few presences or background points for k-fold evaluation.")
    fold_p <- make_folds(length(pres_idx), k, seed)
    fold_b <- make_folds(length(bg_idx), k, seed)
    auc_v <- rep(NA_real_, k)
    tss_v <- rep(NA_real_, k)

    for (i in seq_len(k)) {
      test_idx <- c(pres_idx[fold_p == i], bg_idx[fold_b == i])
      train_idx <- setdiff(seq_along(y), test_idx)
      mod <- fit_maxnet(x[train_idx, , drop = FALSE], y[train_idx], fc, rm)
      pred <- predict_maxnet_safe(mod, x[test_idx, , drop = FALSE], type = "cloglog")
      auc_v[i] <- safe_auc(y[test_idx], pred)
      tss_v[i] <- best_tss_threshold(y[test_idx], pred)$tss
    }
    return(list(auc = mean(auc_v, na.rm = TRUE), tss = mean(tss_v, na.rm = TRUE), method = "kfold"))
  }

  n_pres <- length(pres_idx)
  auc_v <- rep(NA_real_, n_pres)
  tss_v <- rep(NA_real_, n_pres)
  for (i in seq_len(n_pres)) {
    test_p <- pres_idx[i]
    train <- setdiff(seq_along(y), test_p)
    mod <- fit_maxnet(x[train, , drop = FALSE], y[train], fc, rm)
    test_idx <- c(test_p, bg_idx)
    pred <- predict_maxnet_safe(mod, x[test_idx, , drop = FALSE], type = "cloglog")
    auc_v[i] <- safe_auc(y[test_idx], pred)
    tss_v[i] <- best_tss_threshold(y[test_idx], pred)$tss
  }
  list(auc = mean(auc_v, na.rm = TRUE), tss = mean(tss_v, na.rm = TRUE), method = "jackknife")
}

perm_importance_auc <- function(model, x, y) {
  base <- safe_auc(y, predict_maxnet_safe(model, x, type = "cloglog"))
  imp <- vapply(names(x), function(v) {
    x_perm <- x
    x_perm[[v]] <- sample(x_perm[[v]])
    auc_p <- safe_auc(y, predict_maxnet_safe(model, x_perm, type = "cloglog"))
    base - auc_p
  }, numeric(1))
  imp[!is.finite(imp)] <- 0
  perc <- if (sum(imp) == 0) rep(0, length(imp)) else 100 * imp / sum(imp)
  data.frame(variable = names(imp), perm_importance = perc, row.names = NULL)
}

# =============================================================================
# 6) One-species and full benchmark
# =============================================================================
run_species_benchmark <- function(cfg, pts_long, single_stack, species) {
  msg("Benchmark species: %s", species)

  survey <- pts_long |>
    filter(species == !!species, !is.na(value)) |>
    st_transform(cfg$crs_epsg)
  occ <- survey |> filter(value > 0)
  n_pres <- nrow(occ)
  if (n_pres < 3) stop("Too few presences for species: ", species)

  # SoE selection uses all survey rows, including zero values.
  soe_obj <- run_scalepredsdm_soe(single_stack, survey_sf = survey, cfg = cfg, species = species)
  soe_stack <- soe_obj$soe_stack

  bg_xy <- prepare_background_points(cfg, single_stack[[1]], species)

  fr <- cfg$fc_rm_table |> filter(species == !!species)
  if (nrow(fr) != 1) stop("fc/rm settings missing or duplicated for species: ", species)
  fc <- fr$fc
  rm <- fr$rm

  res <- list()
  vi_list <- list()

  for (rep_i in seq_len(cfg$reps)) {
    rep_seed <- cfg$seed + rep_i
    pe_list <- make_pe_geometry_list(occ, cfg$pe_dists_m, cfg$crs_epsg, seed = rep_seed)

    for (pe_name in names(pe_list)) {
      pe_m <- as.integer(gsub("^pe_|m$", "", pe_name))

      for (model_type in c("single_scale", "soe")) {
        pred_stack <- if (model_type == "single_scale") single_stack else soe_stack
        occ_vals <- extract_xy_vals(pred_stack, pe_list[[pe_name]])
        bg_vals <- extract_raster_values(pred_stack, bg_xy[, c("x", "y")])
        bg_df <- cbind(bg_xy, bg_vals)

        x_all <- bind_rows(occ_vals |> select(x, y, everything()), bg_df)
        y_all <- c(rep(1, nrow(occ_vals)), rep(0, nrow(bg_df)))
        x <- x_all |> select(-x, -y)
        keep <- complete.cases(x)
        x <- x[keep, , drop = FALSE]
        y <- y_all[keep]

        mod <- fit_maxnet(x, y, fc = fc, rm = rm)
        ev <- eval_cv_maxnet(x, y, fc, rm, method = "auto", k = cfg$kfold_k,
                             pres_threshold_n = cfg$pres_threshold_n, seed = rep_seed)

        res[[length(res) + 1]] <- data.frame(
          species = species,
          rep = rep_i,
          distance_m = pe_m,
          model_type = model_type,
          fc = fc,
          rm = rm,
          n_pres_original = n_pres,
          n_pres_after_na = sum(y == 1),
          n_bg_after_na = sum(y == 0),
          eval_method = ev$method,
          auc = ev$auc,
          tss = ev$tss,
          stringsAsFactors = FALSE
        )

        vi_list[[length(vi_list) + 1]] <- perm_importance_auc(mod, x, y) |>
          mutate(species = species, rep = rep_i, distance_m = pe_m,
                 model_type = model_type, fc = fc, rm = rm)

        msg("%s | rep=%d | PE=%s | %s | AUC=%.3f | TSS=%.3f",
            species, rep_i, pe_m, model_type, ev$auc, ev$tss)
      }
    }
  }

  list(
    results = bind_rows(res),
    varimp = bind_rows(vi_list),
    SoE_table = soe_obj$SoE_table |> mutate(species = species)
  )
}

summarise_paired_differences <- function(results_all) {
  wide <- results_all |>
    select(species, rep, distance_m, model_type, auc, tss) |>
    tidyr::pivot_wider(names_from = model_type, values_from = c(auc, tss)) |>
    mutate(
      delta_auc = auc_soe - auc_single_scale,
      delta_tss = tss_soe - tss_single_scale
    )

  summary <- wide |>
    group_by(distance_m) |>
    summarise(
      n_pairs = dplyr::n(),
      mean_delta_auc = mean(delta_auc, na.rm = TRUE),
      sd_delta_auc = sd(delta_auc, na.rm = TRUE),
      auc_positive_pct = 100 * mean(delta_auc > 0, na.rm = TRUE),
      mean_delta_tss = mean(delta_tss, na.rm = TRUE),
      sd_delta_tss = sd(delta_tss, na.rm = TRUE),
      tss_positive_pct = 100 * mean(delta_tss > 0, na.rm = TRUE),
      .groups = "drop"
    )

  list(paired = wide, summary = summary)
}

main <- function(cfg) {
  dir.create(file.path(cfg$out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$out_dir, "varimp"), recursive = TRUE, showWarnings = FALSE)

  inputs <- load_inputs(cfg)
  pts_long <- inputs$pts_long
  single_stack <- inputs$single_stack
  terra::crs(single_stack) <- paste0("EPSG:", cfg$crs_epsg)

  sink(file.path(cfg$out_dir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()

  all_res <- list()
  all_vi <- list()
  all_soe <- list()

  for (sp in cfg$species_list) {
    out <- run_species_benchmark(cfg, pts_long, single_stack, sp)
    all_res[[sp]] <- out$results
    all_vi[[sp]] <- out$varimp
    all_soe[[sp]] <- out$SoE_table

    write.csv(out$results, file.path(cfg$out_dir, "results", sprintf("results_%s.csv", sp)), row.names = FALSE)
    write.csv(out$varimp, file.path(cfg$out_dir, "varimp", sprintf("permutation_importance_%s.csv", sp)), row.names = FALSE)
  }

  results_all <- bind_rows(all_res)
  varimp_all <- bind_rows(all_vi)
  soe_all <- bind_rows(all_soe)
  paired <- summarise_paired_differences(results_all)

  write.csv(results_all, file.path(cfg$out_dir, "results", "benchmark_results_all.csv"), row.names = FALSE)
  write.csv(varimp_all, file.path(cfg$out_dir, "varimp", "permutation_importance_all.csv"), row.names = FALSE)
  write.csv(soe_all, file.path(cfg$out_dir, "scale_of_effect", "best_scales_all_species.csv"), row.names = FALSE)
  write.csv(paired$paired, file.path(cfg$out_dir, "results", "paired_differences.csv"), row.names = FALSE)
  write.csv(paired$summary, file.path(cfg$out_dir, "results", "paired_difference_summary.csv"), row.names = FALSE)

  msg("Benchmark complete. Main output: %s", file.path(cfg$out_dir, "results", "benchmark_results_all.csv"))
  invisible(list(results = results_all, varimp = varimp_all, soe = soe_all, paired = paired))
}

# =============================================================================
# 7) Execute
# =============================================================================
if (sys.nframe() == 0) {
  benchmark_out <- main(CONFIG)
}


#Plotting
# =============================================================================
# Plot scale-of-effect correlation curves
# One panel per predictor, coloured lines/dots by species
# Selected SoE scale is marked with a larger open point
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)

# Use the same output folder as the benchmark
out_dir <- CONFIG$out_dir
soe_dir <- file.path(out_dir, "scale_of_effect")

# ---- Helper to standardise scale labels ----
scale_label_fun <- function(x) {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- "point"
  x_chr[x_chr %in% c("0", "native", "point")] <- "point"
  x_chr
}

# ---- Read correlation files ----
corr_files <- list.files(
  soe_dir,
  pattern = "^soe_correlations_.*\\.csv$",
  full.names = TRUE
)

if (length(corr_files) == 0) {
  stop("No SoE correlation files found in: ", soe_dir)
}

read_corr_file <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  species_name <- basename(f) |>
    str_remove("^soe_correlations_") |>
    str_remove("\\.csv$")
  
  names(df) <- tolower(names(df))
  
  if (!"species" %in% names(df)) {
    df$species <- species_name
  }
  
  if (!"predictor" %in% names(df)) {
    pred_alt <- intersect(c("pred", "variable", "layer"), names(df))
    if (length(pred_alt) == 0) {
      stop("No predictor column found in: ", basename(f))
    }
    df <- df |> rename(predictor = all_of(pred_alt[1]))
  }
  
  if (!"scale_m" %in% names(df)) {
    scale_alt <- intersect(c("scale", "radius", "radius_m", "buffer_m"), names(df))
    if (length(scale_alt) == 0) {
      stop("No scale column found in: ", basename(f))
    }
    df <- df |> rename(scale_m = all_of(scale_alt[1]))
  }
  
  cor_alt <- intersect(
    c("correlation", "spearman", "spearman_rho", "rho", "cor", "r"),
    names(df)
  )
  
  if (length(cor_alt) == 0) {
    stop(
      "No correlation column found in: ", basename(f),
      ". Available columns are: ", paste(names(df), collapse = ", ")
    )
  }
  
  df |>
    rename(correlation = all_of(cor_alt[1])) |>
    mutate(
      species = as.character(species),
      predictor = as.character(predictor),
      scale_label = scale_label_fun(scale_m),
      scale_num = suppressWarnings(as.numeric(scale_label))
    ) |>
    select(species, predictor, scale_m, scale_label, scale_num, correlation)
}

corr_df <- map_dfr(corr_files, read_corr_file)

# ---- Read selected best scales ----
best_files <- list.files(
  soe_dir,
  pattern = "^best_scales_for_stack_.*\\.csv$",
  full.names = TRUE
)

read_best_file <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  species_name <- basename(f) |>
    str_remove("^best_scales_for_stack_") |>
    str_remove("\\.csv$")
  
  names(df) <- tolower(names(df))
  
  if (!"species" %in% names(df)) {
    df$species <- species_name
  }
  
  if (!"predictor" %in% names(df)) {
    pred_alt <- intersect(c("pred", "variable", "layer"), names(df))
    if (length(pred_alt) == 0) {
      stop("No predictor column found in: ", basename(f))
    }
    df <- df |> rename(predictor = all_of(pred_alt[1]))
  }
  
  if (!"scale_m" %in% names(df)) {
    scale_alt <- intersect(c("scale", "radius", "radius_m", "buffer_m"), names(df))
    if (length(scale_alt) == 0) {
      stop("No scale column found in: ", basename(f))
    }
    df <- df |> rename(scale_m = all_of(scale_alt[1]))
  }
  
  df |>
    mutate(
      species = as.character(species),
      predictor = as.character(predictor),
      scale_label = scale_label_fun(scale_m)
    ) |>
    select(species, predictor, scale_label)
}

best_df <- map_dfr(best_files, read_best_file)

# ---- Order x-axis: point first, then numeric scales ----
scale_levels <- corr_df |>
  distinct(scale_label) |>
  mutate(scale_num = suppressWarnings(as.numeric(scale_label))) |>
  arrange(is.na(scale_num), scale_num) |>
  pull(scale_label)

scale_levels <- c(
  intersect("point", scale_levels),
  setdiff(scale_levels, "point")
)

corr_df <- corr_df |>
  mutate(scale_label = factor(scale_label, levels = scale_levels))

best_df <- best_df |>
  mutate(scale_label = factor(scale_label, levels = scale_levels))

# ---- Join selected SoE points to their correlation values ----
selected_df <- corr_df |>
  inner_join(
    best_df,
    by = c("species", "predictor", "scale_label")
  )

# ---- Plot ----
p_soe <- ggplot(
  corr_df,
  aes(
    x = scale_label,
    y = correlation,
    colour = species,
    group = species
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 2, alpha = 0.9) +
  geom_point(
    data = selected_df,
    shape = 21,
    fill = "white",
    size = 3.5,
    stroke = 1.2
  ) +
  facet_wrap(~ predictor, scales = "free_y") +
  labs(
    x = "Candidate scale of effect (m)",
    y = "Spearman correlation with abundance",
    colour = "Species"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", colour = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_soe


# =============================================================================
# Plot positional-error effect by species
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)

out_dir <- CONFIG$out_dir
res_file <- file.path(out_dir, "results", "benchmark_results_all.csv")

if (!file.exists(res_file)) {
  stop("Could not find benchmark results file: ", res_file)
}

res <- readr::read_csv(res_file, show_col_types = FALSE)

# ---- Make column names robust across script versions ----
names(res) <- tolower(names(res))

if (!"distance_m" %in% names(res)) {
  if ("pe_m" %in% names(res)) {
    res <- res |> rename(distance_m = pe_m)
  } else {
    stop("Could not find positional-error column. Expected 'distance_m' or 'pe_m'.")
  }
}

if (!"model_type" %in% names(res)) {
  stop("Could not find model_type column.")
}

if (!all(c("auc", "tss") %in% names(res))) {
  stop("Could not find auc and/or tss columns.")
}

# ---- Clean labels ----
res <- res |>
  mutate(
    model_type = case_when(
      model_type %in% c("soe", "SoE", "multi_scale", "multi scale") ~ "Scale-of-effect",
      model_type %in% c("single_scale", "single scale") ~ "Single-scale",
      TRUE ~ as.character(model_type)
    ),
    species = as.character(species),
    distance_m = as.numeric(distance_m)
  )

# ---- Long format for AUC and TSS ----
res_long <- res |>
  select(species, rep, distance_m, model_type, auc, tss) |>
  pivot_longer(
    cols = c(auc, tss),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(metric, auc = "AUC", tss = "TSS")
  )

# ---- Mean and SD per species, distance, model and metric ----
res_sum <- res_long |>
  group_by(species, distance_m, model_type, metric) |>
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se_value = sd_value / sqrt(n),
    .groups = "drop"
  )

# ---- Plot ----
p_pe <- ggplot() +
  # individual replicate values
  geom_point(
    data = res_long,
    aes(
      x = distance_m,
      y = value,
      colour = model_type
    ),
    alpha = 0.25,
    size = 1.4,
    position = position_jitter(width = 2, height = 0)
  ) +
  # mean line
  geom_line(
    data = res_sum,
    aes(
      x = distance_m,
      y = mean_value,
      colour = model_type,
      group = model_type
    ),
    linewidth = 0.9
  ) +
  # mean points
  geom_point(
    data = res_sum,
    aes(
      x = distance_m,
      y = mean_value,
      colour = model_type
    ),
    size = 2.5
  ) +
  # optional uncertainty band as SE error bars
  geom_errorbar(
    data = res_sum,
    aes(
      x = distance_m,
      ymin = mean_value - se_value,
      ymax = mean_value + se_value,
      colour = model_type
    ),
    width = 3,
    alpha = 0.6
  ) +
  facet_grid(metric ~ species, scales = "free_y") +
  scale_x_continuous(
    breaks = sort(unique(res_sum$distance_m))
  ) +
  labs(
    x = "Positional error distance (m)",
    y = "Model performance",
    colour = "Model type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", colour = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_pe


