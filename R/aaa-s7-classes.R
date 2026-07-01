#' S7 classes for utop objects
#'
#' These classes provide the formal object model used by the snake_case utop
#' API. Prefer `utop_object()` and related helpers for user-facing workflows.
#'
#' @param model Variogram model name.
#' @param nugget,unc,cloud,cv Logical parameters.
#' @param r_resol,h_resol,amul,dmul,n_max,n_clus Numeric or integer tuning
#'   parameters.
#' @param cn_areas Number of areas used in cluster-related calculations.
#' @param rs_type,hs_type,clus_type,outfile,wlim_method,uk_trend_support
#'   Character parameters.
#' @param fit_method,g_dist_est,var_clean,max_dist Variogram and distance
#'   parameters.
#' @param partial_overlap Whether partially overlapping areas are supported.
#' @param g_dist_pred Geostatistical distance parameter or computed distance
#'   object, depending on the class.
#' @param wlim,singular_solve,debug_level Additional tuning parameters.
#' @param par_init Initial parameter bounds.
#' @param data Data stored in variogram objects.
#' @param params Numeric variogram model parameters.
#' @param ss_err,criterion Variogram fitting diagnostics.
#' @param observations,prediction_locations Spatial observations and prediction
#'   locations.
#' @param formula Model formula.
#' @param variogram,variogram_cloud,variogram_model Variogram results.
#' @param predictions,simulations Prediction and simulation results.
#' @param d_obs,d_pred,g_dist_obs,g_dist_pred_obs,g_dist_bin,d_bin
#'   Discretisation and geostatistical distance objects.
#' @param overlap_obs,overlap_pred_obs Overlap matrices.
#' @param var_mat_obs,var_mat_pred_obs,var_mat_pred Covariance matrices.
#' @param var_fit,cv_info,weight,removed,uk_residual,check_vario Diagnostics
#'   and intermediate results.
#'
#' @name utop-s7-classes
#' @aliases Utop
#' @aliases UtopParams
#' @aliases UtopVariogram
#' @aliases UtopVariogramCloud
#' @aliases UtopVariogramModel
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
#' @importFrom S7 S7_dispatch S7_inherits class_any class_character
#' @importFrom S7 class_data.frame class_formula class_list class_logical
#' @importFrom S7 class_numeric method methods_register new_class new_generic
#' @importFrom S7 new_property new_S3_class prop
NULL

utop_sf_class <- S7::new_S3_class("sf")
utop_stars_class <- S7::new_S3_class("stars")
utop_matrix_class <- S7::new_S3_class(c("matrix", "array"))

utop_optional <- function(class = S7::class_any) {
  S7::new_property(NULL | class, default = NULL)
}

#' @rdname utop-s7-classes
#' @export
UtopParams <- S7::new_class(
  "UtopParams",
  package = "utop",
  properties = list(
    model = S7::new_property(S7::class_character, default = "Ex1"),
    nugget = S7::new_property(S7::class_logical, default = FALSE),
    unc = S7::new_property(S7::class_logical, default = TRUE),
    r_resol = S7::new_property(S7::class_numeric, default = 100),
    h_resol = S7::new_property(S7::class_numeric, default = 5),
    rs_type = S7::new_property(S7::class_character, default = "rtop"),
    hs_type = S7::new_property(S7::class_character, default = "regular"),
    cloud = S7::new_property(S7::class_logical, default = FALSE),
    amul = S7::new_property(S7::class_numeric, default = 2),
    dmul = S7::new_property(S7::class_numeric, default = 3),
    fit_method = S7::new_property(S7::class_numeric, default = 9),
    g_dist_est = S7::new_property(S7::class_logical, default = FALSE),
    g_dist_pred = S7::new_property(S7::class_logical, default = FALSE),
    var_clean = S7::new_property(S7::class_logical, default = FALSE),
    max_dist = S7::new_property(S7::class_numeric, default = Inf),
    n_max = S7::new_property(S7::class_numeric, default = 10),
    n_clus = S7::new_property(S7::class_numeric, default = 1),
    cn_areas = S7::new_property(S7::class_numeric, default = 100),
    clus_type = utop_optional(S7::class_character),
    outfile = utop_optional(S7::class_character),
    partial_overlap = S7::new_property(S7::class_logical, default = FALSE),
    wlim = S7::new_property(S7::class_numeric, default = 1.5),
    wlim_method = S7::new_property(S7::class_character, default = "all"),
    singular_solve = S7::new_property(S7::class_logical, default = FALSE),
    uk_trend_support = S7::new_property(
      S7::class_character,
      default = "centroid"
    ),
    cv = S7::new_property(S7::class_logical, default = FALSE),
    debug_level = S7::new_property(
      S7::class_numeric,
      default = if (interactive()) 1 else 0
    ),
    par_init = utop_optional()
  ),
  validator = function(self) {
    if (length(self@model) != 1L) {
      return("model must be a character scalar")
    }
    if (!self@model %in% c("Exp", "Sph", "Gau", "Sp1", "Ex1", "Fra")) {
      return(paste("model", self@model, "not implemented"))
    }
    if (length(self@cloud) != 1L) {
      return("cloud must be a logical scalar")
    }
    if (length(self@r_resol) != 1L || self@r_resol <= 0) {
      return("r_resol must be a positive number")
    }
    NULL
  }
)

#' @rdname utop-s7-classes
#' @export
UtopVariogramModel <- S7::new_class(
  "UtopVariogramModel",
  package = "utop",
  properties = list(
    model = S7::new_property(S7::class_character, default = "Ex1"),
    params = S7::new_property(S7::class_numeric, default = c(1, 1, 0, 0, 1)),
    ss_err = utop_optional(S7::class_numeric),
    criterion = utop_optional()
  ),
  validator = function(self) {
    if (length(self@model) != 1L) {
      return("model must be a character scalar")
    }
    if (!is.numeric(self@params) || length(self@params) == 0L) {
      return("params must be a non-empty numeric vector")
    }
    NULL
  }
)

#' @rdname utop-s7-classes
#' @export
UtopVariogram <- S7::new_class(
  "UtopVariogram",
  package = "utop",
  constructor = function(data = data.frame()) {
    S7::new_object(S7::S7_object(), data = data)
  },
  properties = list(data = S7::new_property(S7::class_data.frame)),
  validator = function(self) {
    if (nrow(self@data) == 0L) {
      return(NULL)
    }
    required <- c("np", "dist", "gamma", "a1", "a2")
    missing <- setdiff(required, names(self@data))
    if (length(missing) > 0L) {
      return(paste(
        "missing variogram columns:",
        paste(missing, collapse = ", ")
      ))
    }
    NULL
  }
)

#' @rdname utop-s7-classes
#' @export
UtopVariogramCloud <- S7::new_class(
  "UtopVariogramCloud",
  package = "utop",
  constructor = function(data = data.frame()) {
    S7::new_object(S7::S7_object(), data = data)
  },
  properties = list(data = S7::new_property(S7::class_data.frame))
)

#' @rdname utop-s7-classes
#' @export
Utop <- S7::new_class(
  "Utop",
  package = "utop",
  properties = list(
    observations = utop_optional(),
    prediction_locations = utop_optional(),
    formula = utop_optional(S7::class_formula),
    params = S7::new_property(UtopParams),
    variogram = utop_optional(UtopVariogram),
    variogram_cloud = utop_optional(UtopVariogramCloud),
    variogram_model = utop_optional(UtopVariogramModel),
    predictions = utop_optional(),
    simulations = utop_optional(),
    d_obs = utop_optional(),
    d_pred = utop_optional(),
    g_dist_obs = utop_optional(),
    g_dist_pred = utop_optional(),
    g_dist_pred_obs = utop_optional(),
    g_dist_bin = utop_optional(),
    d_bin = utop_optional(),
    overlap_obs = utop_optional(),
    overlap_pred_obs = utop_optional(),
    var_mat_obs = utop_optional(),
    var_mat_pred_obs = utop_optional(),
    var_mat_pred = utop_optional(),
    var_fit = utop_optional(),
    cv_info = utop_optional(),
    weight = utop_optional(),
    removed = utop_optional(),
    uk_residual = utop_optional(),
    check_vario = utop_optional()
  )
)
