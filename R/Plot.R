#' Plot Force of Infection (FOI) curves from fitted catalytic models
#'
#' Produces a line plot of estimated force of infection (FOI) over age for one
#' or more fitted catalytic models, with optional bootstrap confidence ribbons
#' and optional bars for a known/true FOI by age interval.
#'
#' @param cat_models A single fitted catalytic model object or a named list of such
#'   objects. Each model is expected to provide:
#'   \itemize{
#'     \item \code{foi_t}: function mapping ages and parameters to FOI, e.g. \code{function(t, par) \{...\}}.
#'     \item \code{params_MLE}: named numeric vector (or list) of MLE parameters.
#'     \item \code{bootparams}: (required if \code{confint = TRUE} and non-spline) a matrix of
#'       bootstrap parameter draws, one row per bootstrap.
#'     \item \code{spline_pi_t}: (spline models only) a \code{smooth.spline} fit used by \code{foi_t}.
#'     \item \code{boot_spline_pi_t}: (optional, spline CIs) list of \code{smooth.spline} fits,
#'       one per bootstrap. If absent, ribbons are omitted for spline models even when \code{confint = TRUE}.
#'   }
#'   If \code{cat_models} is unnamed, legend labels are auto-generated as
#'   \code{"Model 1"}, \code{"Model 2"}, etc.
#' @param from,to Numeric scalars. Age range to evaluate and plot (inclusive).
#' @param confint Logical. If \code{TRUE}, draws 95\% ribbons from bootstrap
#'   percentiles. For non-spline models, FOI is recomputed for each age grid point
#'   using rows of \code{bootparams}. For spline models, ribbons are shown only when
#'   \code{boot_spline_pi_t} is supplied.
#' @param true_foi Optional. A list or data frame with components
#'   \code{t} (a two-column matrix of age intervals \[a, b\]) and
#'   \code{foi} (numeric of the same length) to draw semi-transparent rectangles
#'   indicating a known/true FOI per interval.
#' @param line_colors Optional. Either (i) a named character vector of colors where
#'   names match \code{names(cat_models)}, or (ii) an unnamed vector recycled across
#'   models. If omitted, a colorblind-friendly HCL palette (when available) is used.
#' @param xmin,ymin,xmax,ymax Numeric limits passed to \code{coord_cartesian()}.
#'   Use \code{NA} for any bound you want ggplot2 to choose from the data.
#' @param true_fill Fill color for the \code{true_foi} rectangles (default \code{"grey80"}).
#'
#' @return A \code{ggplot} object showing FOI curves (and ribbons/rectangles if requested).
#' @export
#'
#' @details
#' \itemize{
#'   \item FOI is evaluated on an internal grid of 100 ages from \code{from} to \code{to}.
#'   \item Ribbons use the 2.5\% and 97.5\% quantiles at each grid point.
#'   \item For spline models, FOI is computed via \code{foi_t(t, spline_pi_t)}.
#'   \item For non-spline models, FOI is computed via \code{foi_t(t, params_MLE)}.
#' }
#'
#' @examples
#' \dontrun{
#' # Age-intervals and midpoints
#' t <- matrix(
#'   c(
#'     0, 1,
#'     1, 5,
#'     5, 10,
#'     10, 15,
#'     15, 20,
#'     20, 30,
#'     30, 40,
#'     40, 50,
#'     50, 60
#'   ),
#'   ncol = 2, byrow = TRUE
#' )
#' t_mid <- rowMeans(t)
#'
#' # Observed positives and totals
#' y <- c(28, 224, 577, 712, 740, 808, 848, 873, 894)
#' n <- rep(1000, length(y))
#'
#' # True FOI by interval (rounded)
#' true_foi <- round(c(
#'   0.02958701, 0.06343983, 0.18413499, 0.08801705, 0.04843492,
#'   0.03307834, 0.04797940, 0.03155037, 0.01551178
#' ), 4)
#'
#' # Fit Farringtons FOI on exact-age midpoints
#' Farrington_exactt <- FoiFromCatalyticModel(
#'   t = t_mid,
#'   y = y,
#'   n = n,
#'   catalytic_model_type = "SimpleCatalytic",
#'   foi_functional_form = "Farringtons",
#'   lower = c(0, 0, 0),
#'   boot_num = 100
#' )
#'
#' # Plot FOI with bootstrap ribbons and true interval FOI as rectangles
#' plotFOI(
#'   list(Farrington_exactt = Farrington_exactt),
#'   from = 0, to = 60, confint=TRUE,
#'   true_foi = list(foi = true_foi, t = t)
#' )
#' }
#'
#' @seealso \code{\link{FoiFromCatalyticModel}}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_rect
#' @importFrom ggplot2 labs theme_minimal scale_color_manual scale_fill_manual coord_cartesian
#' @importFrom stats setNames
#' @importFrom grDevices hcl.pals hcl.colors rainbow
plotFOI <- function(
    cat_models, from, to, confint = FALSE,
    true_foi = NA,
    line_colors = NULL, # named vector of hex colours for model lines (names must match model names)
    xmin = 0, ymin = 0, xmax = NA, ymax = NA,
    true_fill = "grey80" # background histogram fill
) {

  safe_palette <- function(n) {
    pals <- tryCatch(grDevices::hcl.pals(), error = function(e) character())
    preferred <- c("Okabe-Ito", "Dark 3", "Set 2", "Set 3", "Vibrant", "Warm")
    pick <- intersect(preferred, pals)
    if (length(pick) > 0) grDevices::hcl.colors(n, palette = pick[1]) else grDevices::rainbow(n)
  }

  if (!is.list(cat_models)) cat_models <- list(cat_models)
  if (is.null(names(cat_models))) names(cat_models) <- paste0("Model ", seq_along(cat_models))
  model_names <- names(cat_models)

  plot_data <- lapply(seq_along(cat_models), function(i) {
    model <- cat_models[[i]]
    name  <- model_names[i]
    tau <- model$tau
    t_grid <- seq(from, to, length.out = 100)
    is_spline <- !is.null(model[["spline_pi_t"]])

    foi_grid <- if (is_spline) model$foi_t(t_grid, model$spline_pi_t)
    else  model$foi_t(t_grid, unlist(model$params_MLE))

    if (isTRUE(confint)) {
      if (is_spline) {
        if (!is.null(model$boot_y)) {
          # Build a new spline for each bootstrap y, then evaluate FOI on t_grid
          t_fit <- if (is.matrix(model$t) && ncol(model$t) == 2) rowMeans(model$t) else model$t
          g <- length(t_fit)
          boot_mat <- matrix(model$boot_y, ncol = g, byrow = TRUE)

          boot_foi_mat <- apply(boot_mat, 1, function(yb) {
            bspline <- smooth.spline(t_fit, yb / model$n)
            model$foi_t(t_grid, bspline)
          })

          if (is.vector(boot_foi_mat)) boot_foi_mat <- matrix(boot_foi_mat, nrow = length(t_grid))
          if (nrow(boot_foi_mat) != length(t_grid)) boot_foi_mat <- t(boot_foi_mat)

          foi_CI <- t(apply(boot_foi_mat, 1, function(v)
            quantile(v, probs = c(0.025, 0.975), na.rm = TRUE)
          ))
          colnames(foi_CI) <- c("foi_lower","foi_upper")

          data.frame(age = t_grid, foi = foi_grid, foi_lower = foi_CI[,1], foi_upper = foi_CI[,2],
                     model = name, row.names = NULL)
        } else {
          data.frame(age = t_grid, foi = foi_grid, model = name, row.names = NULL)
        }
      } else {
        stopifnot(!is.null(model$bootparams))
        foi_CI <- matrix(NA_real_, nrow = length(t_grid), ncol = 2)
        for (k in seq_along(t_grid)) {
          boot_fois <- apply(model$bootparams, 1, function(x) model$foi_t(t_grid[k], x))
          foi_CI[k, ] <- as.numeric(quantile(boot_fois, probs = c(0.025, 0.975)))
        }
        colnames(foi_CI) <- c("foi_lower","foi_upper")
        data.frame(age = t_grid, foi = foi_grid, foi_lower = foi_CI[,1], foi_upper = foi_CI[,2],
                   model = name, row.names = NULL)
      }
    } else {
      data.frame(age = t_grid, foi = foi_grid, model = name, row.names = NULL)
    }
  })

  full_df <- do.call(rbind, plot_data)
  full_df$model <- factor(full_df$model, levels = model_names)

  p <- ggplot2::ggplot(full_df, ggplot2::aes(x = age, y = foi, color = model))

  have_true <- !is.na(true_foi)[1]
  if (have_true) {
    stopifnot(is.list(true_foi) || is.data.frame(true_foi))
    stopifnot(all(c("t","foi") %in% names(true_foi)))
    stopifnot(is.matrix(true_foi$t) && ncol(true_foi$t) == 2)
    stopifnot(nrow(true_foi$t) == length(true_foi$foi))

    df_true <- data.frame(
      xmin = true_foi$t[,1],
      xmax = true_foi$t[,2],
      ymin = 0,
      ymax = as.numeric(true_foi$foi)
    )

    p <- p + ggplot2::geom_rect(
      data = df_true,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = true_fill,
      color = "white",
      linewidth = 0.25
    )
  }

  default_cols <- stats::setNames(safe_palette(length(model_names)), model_names)
  if (!is.null(line_colors)) {
    if (!is.null(names(line_colors)) && length(names(line_colors))) {
      final_cols <- default_cols
      overlap <- intersect(names(line_colors), model_names)
      final_cols[overlap] <- line_colors[overlap]
    } else {
      final_cols <- setNames(rep(line_colors, length.out = length(model_names)), model_names)
    }
  } else {
    final_cols <- default_cols
  }

  if (isTRUE(confint) && all(c("foi_lower","foi_upper") %in% names(full_df))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = foi_lower, ymax = foi_upper, fill = model),
      alpha = 0.20, color = NA, show.legend = FALSE
    ) +
      ggplot2::scale_fill_manual(values = final_cols)
  }

  p +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = final_cols, breaks = model_names, name = NULL) +
    ggplot2::labs(x = "Age", y = "FOI") +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
}
