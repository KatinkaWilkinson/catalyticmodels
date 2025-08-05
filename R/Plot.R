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
      lower_par <- sapply(model$params_CI, function(ci) ci[1])
      upper_par <- sapply(model$params_CI, function(ci) ci[2])
      foi_CI <- cbind(
        lower = model$foi_t(t_grid, lower_par),
        upper = model$foi_t(t_grid, upper_par)
      )

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



plot_foi_grid <- function(cat_models, from, to, confint = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use this function.")
  }

  # Coerce to list if a single model is given
  if (!is.list(cat_models)) {
    cat_models <- list(cat_models)
  }

  # Generate model names from list names, or fallback
  model_names <- names(cat_models)
  if (is.null(model_names)) {
    model_names <- paste0("Model_", seq_along(cat_models))
  } else {
    missing_names <- which(is.na(model_names) | model_names == "")
    model_names[missing_names] <- paste0("Model_", missing_names)
  }

  plot_data <- lapply(seq_along(cat_models), function(i) {
    model <- cat_models[[i]]
    name <- model_names[i]

    if (!is.function(model$foi_t)) {
      stop(paste("Model", name, "has no valid foi_t function."))
    }

    t_grid <- seq(from, to, length.out = 100)
    foi_grid <- model$foi_t(t_grid, unlist(model$params_MLE))

    if (confint) {
      lower_par <- sapply(model$params_CI, function(ci) ci[1])
      upper_par <- sapply(model$params_CI, function(ci) ci[2])
      foi_CI <- cbind(
        lower = model$foi_t(t_grid, lower_par),
        upper = model$foi_t(t_grid, upper_par)
      )

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
