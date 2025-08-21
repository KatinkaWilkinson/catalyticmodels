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
plot_foi_grid <- function(cat_models, from, to, confint = FALSE) {
  if (!is.list(cat_models) || is.null(names(cat_models))) {
    cat_models <- list(cat_models)
  }

  model_names <- names(cat_models)

  plot_data <- lapply(seq_along(cat_models), function(i) {
    model <- cat_models[[i]]
    name <- model_names[i]

    t_grid <- seq(from, to, length.out = 100)
    foi_grid <- model$foi_t(t_grid, unlist(model$params_MLE))

    if (confint) {
      foi_CI <- matrix(ncol = 2, nrow = length(t_grid))
      for (i in 1:length(t_grid)) {
        t <- t_grid[i]
        boot_fois <- apply(model$bootparams, 1, function(x) {model$foi_t(t, x)})
        foi_CI[i,] <- t(quantile(boot_fois, probs = c(0.025, 0.975)))
      }
      colnames(foi_CI) <- c("lower", "upper")
      print(foi_CI)
      data.frame(
        age = t_grid,
        foi = foi_grid,
        foi_lower = foi_CI[, "lower"],
        foi_upper = foi_CI[, "upper"],
        model = name
      )
    } else {
      data.frame(
        age = t_grid,
        foi = foi_grid,
        model = name
      )
    }
  })

  full_df <- do.call(rbind, plot_data)

  p <- ggplot2::ggplot(full_df, ggplot2::aes(x = age, y = foi, color = model)) +
    ggplot2::geom_line()

  if (confint) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = foi_lower, ymax = foi_upper, fill = model),
        alpha = 0.2,
        color = NA
      )
  }

  p + ggplot2::labs(
    x = "Age",
    y = "Force of Infection (FOI)",
    title = "Estimated Force of Infection"
  ) + ggplot2::theme_minimal()
}

##### BELOW MIGHT BE OUTDATED!!!!!!!!!!!

#
#
# plot_foi_grid <- function(cat_models, from, to, confint = FALSE) {
#   if (!requireNamespace("ggplot2", quietly = TRUE)) {
#     stop("Please install ggplot2 to use this function.")
#   }
#
#   # Coerce to list if a single model is given
#   if (!is.list(cat_models)) {
#     cat_models <- list(cat_models)
#   }
#
#   # Generate model names from list names, or fallback
#   model_names <- names(cat_models)
#   if (is.null(model_names)) {
#     model_names <- paste0("Model_", seq_along(cat_models))
#   } else {
#     missing_names <- which(is.na(model_names) | model_names == "")
#     model_names[missing_names] <- paste0("Model_", missing_names)
#   }
#
#   plot_data <- lapply(seq_along(cat_models), function(i) {
#     model <- cat_models[[i]]
#     name <- model_names[i]
#
#     if (!is.function(model$foi_t)) {
#       stop(paste("Model", name, "has no valid foi_t function."))
#     }
#
#     t_grid <- seq(from, to, length.out = 100)
#     foi_grid <- model$foi_t(t_grid, unlist(model$params_MLE))
#
#     if (confint) {
#       lower_par <- sapply(model$params_CI, function(ci) ci[1])
#       upper_par <- sapply(model$params_CI, function(ci) ci[2])
#       foi_CI <- cbind(
#         lower = model$foi_t(t_grid, lower_par),
#         upper = model$foi_t(t_grid, upper_par)
#       )
#
#       data.frame(
#         age = t_grid,
#         foi = foi_grid,
#         foi_lower = foi_CI[, "lower"],
#         foi_upper = foi_CI[, "upper"],
#         model = name
#       )
#     } else {
#       data.frame(
#         age = t_grid,
#         foi = foi_grid,
#         model = name
#       )
#     }
#   })
#
#   full_df <- do.call(rbind, plot_data)
#
#   p <- ggplot2::ggplot(full_df, ggplot2::aes(x = age, y = foi, color = model)) +
#     ggplot2::geom_line()
#
#   if (confint) {
#     p <- p +
#       ggplot2::geom_ribbon(
#         ggplot2::aes(ymin = foi_lower, ymax = foi_upper, fill = model),
#         alpha = 0.2,
#         color = NA
#       )
#   }
#
#   p + ggplot2::labs(
#     x = "Age",
#     y = "Force of Infection (FOI)",
#     title = "Estimated Force of Infection"
#   ) + ggplot2::theme_minimal()
# }
