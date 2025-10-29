#' Plot Force of Infection Estimates from Catalytic Models
#'
#' Creates a line plot (optionally with confidence intervals) of the
#' estimated force of infection (FOI) from one or more fitted catalytic
#' models over a specified age range.
#'
#' @param cat_models A fitted catalytic model object, or a named list of such objects.
#'   Each model must contain:
#'   \itemize{
#'     \item `foi_t`: a function taking a vector of ages and a parameter vector, returning FOI estimates.
#'     \item `params_MLE`: a vector of maximum likelihood estimates for model parameters.
#'     \item `params_CI`: (optional) a list of confidence intervals for parameters, where each element is a numeric vector \[lower, upper\].
#'   }
#' @param from Numeric scalar. Lower bound of the age range for plotting.
#' @param to Numeric scalar. Upper bound of the age range for plotting.
#' @param confint Logical. If `TRUE`, plots shaded confidence intervals based on `params_CI`.
#'
#' @return A \code{ggplot2} object showing the FOI curves for the provided model(s).
#' @export
#'
#' @examples
#' # Example using a hypothetical catalytic model object `my_model`
#' # plot_foi_grid(my_model, from = 0, to = 50, confint = TRUE)
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_minimal

# plot_foi_grid <- function(cat_models, from, to, confint = FALSE, true_foi = NA, xmin=0, ymin=0, xmax=NA, ymax=NA) {
#
#   # --- Safe palette helper (handles old R/grDevices) ---
#   safe_palette <- function(n) {
#     pals <- tryCatch(grDevices::hcl.pals(), error = function(e) character())
#     preferred <- c("Okabe-Ito", "Dark 3", "Set 2", "Set 3", "Vibrant", "Warm")
#     pick <- intersect(preferred, pals)
#     if (length(pick) > 0) grDevices::hcl.colors(n, palette = pick[1]) else grDevices::rainbow(n)
#   }
#
#   # Normalize input and preserve order
#   if (!is.list(cat_models)) cat_models <- list(cat_models)
#   if (is.null(names(cat_models))) names(cat_models) <- paste0("Model ", seq_along(cat_models))
#   model_names <- names(cat_models)
#
#   plot_data <- lapply(seq_along(cat_models), function(i) {
#     model <- cat_models[[i]]
#     name  <- model_names[i]
#
#     t_grid <- seq(from, to, length.out = 100)
#     is_spline <- !is.null(model[["spline_pi_t"]])
#
#     # Baseline FOI curve
#     foi_grid <- if (is_spline) {
#       model$foi_t(t_grid, model$spline_pi_t)
#     } else {
#       model$foi_t(t_grid, unlist(model$params_MLE))
#     }
#
#     # Confidence intervals
#     if (isTRUE(confint)) {
#       if (is_spline) {
#         if (!is.null(model$boot_spline_pi_t) && is.list(model$boot_spline_pi_t)) {
#           foi_CI <- matrix(NA_real_, nrow = length(t_grid), ncol = 2)
#           for (k in seq_along(t_grid)) {
#             t <- t_grid[k]
#             boot_fois <- vapply(model$boot_spline_pi_t, function(f) model$foi_t(t, f), numeric(1))
#             foi_CI[k, ] <- as.numeric(quantile(boot_fois, probs = c(0.025, 0.975), na.rm = TRUE))
#           }
#           colnames(foi_CI) <- c("lower", "upper")
#           data.frame(age = t_grid, foi = foi_grid,
#                      foi_lower = foi_CI[, "lower"], foi_upper = foi_CI[, "upper"],
#                      model = name, row.names = NULL)
#         } else {
#           data.frame(age = t_grid, foi = foi_grid, model = name, row.names = NULL)
#         }
#       } else {
#         stopifnot(!is.null(model$bootparams))
#         foi_CI <- matrix(ncol = 2, nrow = length(t_grid))
#         for (k in seq_along(t_grid)) {
#           t <- t_grid[k]
#           boot_fois <- apply(model$bootparams, 1, function(x) model$foi_t(t, x))
#           foi_CI[k, ] <- as.numeric(quantile(boot_fois, probs = c(0.025, 0.975)))
#         }
#         colnames(foi_CI) <- c("lower", "upper")
#         data.frame(age = t_grid, foi = foi_grid,
#                    foi_lower = foi_CI[, "lower"], foi_upper = foi_CI[, "upper"],
#                    model = name, row.names = NULL)
#       }
#     } else {
#       data.frame(age = t_grid, foi = foi_grid, model = name, row.names = NULL)
#     }
#   })
#
#   full_df <- do.call(rbind, plot_data)
#   full_df$model <- factor(full_df$model, levels = model_names)  # preserve input order
#
#   p <- ggplot2::ggplot(full_df, ggplot2::aes(x = age, y = foi, color = model)) +
#     ggplot2::geom_line(linewidth = 0.8)
#
#   if (isTRUE(confint) && all(c("foi_lower", "foi_upper") %in% names(full_df))) {
#     p <- p +
#       ggplot2::geom_ribbon(
#         ggplot2::aes(ymin = foi_lower, ymax = foi_upper, fill = model),
#         alpha = 0.20, color = NA
#       )
#   }
#
#   # "True" FOI overlay with legend entry FIRST
#   have_true <- !is.na(true_foi)[1]
#   if (have_true) {
#     stopifnot(is.list(true_foi) || is.data.frame(true_foi))
#     stopifnot(all(c("t", "foi") %in% names(true_foi)))
#     df_true <- data.frame(age = true_foi$t, foi = true_foi$foi)
#
#     p <- p +
#       ggplot2::geom_line(
#         data = df_true,
#         ggplot2::aes(x = age, y = foi, color = "True Population FOI"),
#         inherit.aes = FALSE, linewidth = 0.9
#       ) +
#       ggplot2::geom_point(
#         data = df_true,
#         ggplot2::aes(x = age, y = foi, color = "True Population FOI"),
#         inherit.aes = FALSE, size = 1.6
#       )
#   }
#
#   # Colors + legend order
#   model_cols <- stats::setNames(safe_palette(length(model_names)), model_names)
#   color_values <- if (have_true) c("True Population FOI" = "black", model_cols) else model_cols
#   legend_breaks <- if (have_true) c("True Population FOI", model_names) else model_names
#
#   p +
#     ggplot2::scale_color_manual(values = color_values, breaks = legend_breaks, name = NULL) +
#     ggplot2::scale_fill_manual(values = model_cols, guide = "none") +
#     ggplot2::labs(
#       x = "Age",
#       y = "FOI"
#     ) +
#     ggplot2::theme_minimal()+
#     ggplot2::coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
# }
#




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
    else            model$foi_t(t_grid, unlist(model$params_MLE))

    if (isTRUE(confint)) {
      if (is_spline) {
        if (!is.null(model$boot_spline_pi_t) && is.list(model$boot_spline_pi_t)) {
          foi_CI <- matrix(NA_real_, nrow = length(t_grid), ncol = 2)
          for (k in seq_along(t_grid)) {
            boot_fois <- vapply(model$boot_spline_pi_t, function(f) model$foi_t(t_grid[k], f), numeric(1))
            foi_CI[k, ] <- as.numeric(quantile(boot_fois, probs = c(0.025, 0.975), na.rm = TRUE))
          }
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
