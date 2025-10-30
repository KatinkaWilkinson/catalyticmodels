#' Compare Models by AIC
#'
#' Computes the Akaike Information Criterion (AIC) for one or more fitted
#' catalytic models and prints a ranked summary (smaller AIC is better).
#'
#' @param cat_models A **named list** of catalytic model objects (e.g.,
#'   the elements returned by \code{FoiFromCatalyticModel()} and collected
#'   into a list). Each model must contain at least:
#'   \itemize{
#'     \item \code{t}, \code{y}, \code{n} (the data used to fit the model),
#'     \item \code{params_MLE} (named vector of MLEs),
#'     \item \code{pi_t}, \code{group_pi}, \code{group_foi} (model functions),
#'     \item \code{rho} and \code{param_names} (as used by the likelihood).
#'   }
#'   Spline-based fits are not supported (i.e., models with \code{spline_pi_t}
#'   produce an error).
#'
#' @details
#' For each model, the function evaluates the negative binomial log-likelihood
#' via \code{neg_total_binom_loglik()} and computes
#' \deqn{\mathrm{AIC} = 2\,(-\log L) + 2k,}
#' where \(k\) is the number of estimated parameters. Results are printed as a
#' ranked list and also returned in a structured object.
#'
#' **Important:** AIC values are only comparable across models fitted to the
#' same data \(\[t, y, n\]\). Ensure you pass models trained on identical
#' datasets.
#'
#' @return
#' If a **single** model is supplied, returns a list with:
#' \itemize{
#'   \item \code{ranking_string}: a formatted line with the model's AIC, and
#'   \item \code{AIC_vals}: a named list containing \code{AIC}, \code{k}, and \code{n}.
#' }
#'
#' If **multiple** models are supplied, a ranked table is printed to the console
#' and the function returns:
#' \itemize{
#'   \item \code{AIC_vals}: a named list (one element per model), each containing
#'         \code{AIC}, \code{k}, and \code{n}.
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose m1 and m2 are fitted catalytic models from FoiFromCatalyticModel()
#' mods <- list(Constant = m1, Farringtons = m2)
#' AIC_comp(mods)  # prints ranking; returns per-model AIC, k, n
#' }
#'
#' @seealso \code{\link{AICc_comp}} for small-sample corrected AIC.
#' @export
AIC_comp <- function(cat_models) {
  if (!is.list(cat_models[[1]][[1]])) {
    stop("cat_models must be a named list of catalytic models outputted by FoiFromCatalyticModels.")
  }

  calc_AIC <- function(cat_model) {
    # not defined for spline-based fits in this package layout
    if (!is.null(cat_model[["spline_pi_t"]])) {
      stop("AIC/AICc not available for spline-based models (no parametric likelihood).")
    }

    # calculate AIC
    num_groups <- length(cat_model$n)
    num_estimated_pars <- length(unlist(cat_model$params_MLE))
    neg_loglik <- neg_total_binom_loglik(par=unlist(cat_model$params_MLE), pi_t=cat_model$pi_t, group_pi=cat_model$group_pi, group_foi=cat_model$group_foi, t=cat_model$t, y=cat_model$y, n=cat_model$n, rho=cat_model$rho, param_names=cat_model$param_names)
    # AIC <- -2*logLik + 2*k
    AIC <- 2*neg_loglik + 2*num_estimated_pars

    return(list(AIC=AIC, k = num_estimated_pars, n = num_groups))
  }

  # Compute AICs
  res_list <- lapply(cat_models, calc_AIC)
  model_names <- names(cat_models)

  # If only one model, just return that model's AIC results
  if (length(res_list) == 1) {
    nm <- model_names[1]
    a  <- res_list[[1]]
    ranking_str <- sprintf("%s: AIC = %.3f",nm, a$AIC)
    AIC_vals <- setNames(list(list(AIC = a$AIC, k = a$k, n = a$n)), nm)
    return(list(ranking_string = ranking_str, AIC_vals = AIC_vals))
  }

  # Multiple models: rank by AIC (smaller is better)
  AIC_vec  <- sapply(res_list, function(x) x$AIC)
  ord <- order(AIC_vec, decreasing = FALSE)

   ranking_str <- (function(names, aic, ord) {
    nm  <- names[ord]
    av  <- aic[ord]
    wnm <- max(nchar(nm))
    wa  <- max(nchar(formatC(av, format = "f", digits = 3, big.mark = ",")))
    lines <- sprintf("%2d) %-*s  AIC = %*s",
                     seq_along(ord), wnm, nm,
                     wa, formatC(av, format = "f", digits = 3, big.mark = ","))
    paste(lines, collapse = "\n")
  })(model_names, AIC_vec, ord)

  cat(ranking_str, "\n")

  AIC_vals <- setNames(
    lapply(seq_along(res_list), function(i) {
      list(AIC = AIC_vec[i], k = res_list[[i]]$k, n = res_list[[i]]$n)}),
    model_names
  )

  return(AIC_vals = AIC_vals)
}

#' Compare Models by AICc
#'
#' Computes the small-sample corrected Akaike Information Criterion (AICc) for
#' one or more fitted catalytic models and prints a ranked summary
#' (smaller AICc is better). Models with undefined AICc are listed but not ranked.
#'
#' @param cat_models A **named list** of catalytic model objects (for example,
#'   elements returned by \code{FoiFromCatalyticModel()} and collected into a
#'   list, or the list produced by \code{FoiFromCatalyticModels()}). Each model
#'   must contain at least:
#'   \itemize{
#'     \item \code{t}, \code{y}, \code{n} (the data used to fit the model),
#'     \item \code{params_MLE} (named vector of MLEs),
#'     \item \code{pi_t}, \code{group_pi}, \code{group_foi} (model functions),
#'     \item \code{rho} and \code{param_names} (as used by the likelihood).
#'   }
#'   Spline-based fits are not supported (models with \code{spline_pi_t} produce
#'   an error).
#'
#' @details
#' For each model, the function evaluates the negative binomial log-likelihood
#' via \code{neg_total_binom_loglik()} and computes
#' \deqn{\mathrm{AIC} = 2\,(-\log L) + 2k,\quad
#'       \mathrm{AICc} = \mathrm{AIC} + \frac{2k(k+1)}{n - k - 1},}
#' where \(k\) is the number of estimated parameters and \(n\) is the number of
#' age groups. If \(n - k - 1 \le 0\), AICc is undefined and reported as \code{NA};
#' such models are listed in an “undefined AICc” note and excluded from the
#' ranking.
#'
#' **Important:** AICc values are only comparable across models fitted to the
#' same data \(\[t, y, n\]\). Ensure the supplied models were trained on
#' identical datasets.
#'
#' @return
#' If a **single** model is supplied, returns a list with:
#' \itemize{
#'   \item \code{ranking_string}: a formatted line with the model's AICc (or \code{NA}),
#'   \item \code{AICc_vals}: a named list containing \code{AICc}, \code{k}, and \code{n}.
#' }
#'
#' If **multiple** models are supplied, a ranked table is printed to the console
#' (undefined-AICc models are listed separately) and the function returns:
#' \itemize{
#'   \item \code{AICc_vals}: a named list (one element per model), each containing
#'         \code{AICc}, \code{k}, and \code{n}.
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose m1 and m2 are fitted catalytic models from FoiFromCatalyticModel()
#' mods <- list(Constant = m1, Farringtons = m2)
#' AICc_comp(mods)  # prints ranking; returns per-model AICc, k, n
#' }
#'
#' @seealso \code{\link{AIC_comp}} for the uncorrected AIC comparison.
#' @export
AICc_comp <- function(cat_models) {
  if (!is.list(cat_models[[1]][[1]])) {
    stop("cat_models must be a named list of catalytic models outputted by FoiFromCatalyticModels.")
  }

  calc_AICc <- function(cat_model) {
    # not defined for spline-based fits in this package layout
    if (!is.null(cat_model[["spline_pi_t"]])) {
      stop("AIC/AICc not available for spline-based models (no parametric likelihood).")
    }

    # calculate AIC
    num_groups <- length(cat_model$n)
    num_estimated_pars <- length(unlist(cat_model$params_MLE))
    neg_loglik <- neg_total_binom_loglik(par=unlist(cat_model$params_MLE), pi_t=cat_model$pi_t, group_pi=cat_model$group_pi, group_foi=cat_model$group_foi, t=cat_model$t, y=cat_model$y, n=cat_model$n, rho=cat_model$rho, param_names=cat_model$param_names)
    # AIC <- -2*logLik + 2*k
    AIC <- 2*neg_loglik + 2*num_estimated_pars
    # AICc <- AIC + (2*k*(k+1)) / (n - k - 1)
    if (num_groups - num_estimated_pars - 1 <= 0) {
      AICc <- NA
    } else {
      AICc <- AIC + (2*num_estimated_pars*(num_estimated_pars+1)) / (num_groups - num_estimated_pars - 1)
    }

    return(list(AICc = AICc, k = num_estimated_pars, n = num_groups))
  }

  res_list    <- lapply(cat_models, calc_AICc)
  model_names <- names(cat_models)

  if (length(res_list) == 1) {
    nm <- model_names[1]
    a  <- res_list[[1]]
    ranking_str <- if (is.finite(a$AICc)) sprintf("%s: AICc = %.3f", nm, a$AICc)
    else sprintf("%s: AICc = NA (undefined)", nm)
    AICc_vals <- setNames(list(list(AICc = a$AICc, k = a$k, n = a$n)), nm)
    return(list(ranking_string = ranking_str, AICc_vals = AICc_vals))
  }

  AICc_vec <- vapply(res_list, function(x) x$AICc, numeric(1))
  finite   <- is.finite(AICc_vec)
  ord      <- order(AICc_vec[finite], decreasing = FALSE)

  ranked_names <- model_names[finite][ord]
  ranked_vals  <- AICc_vec[finite][ord]
  ranking_str <- (function(ranked_names, ranked_vals, model_names, finite, label = "AICc") {
    nm <- ranked_names
    av <- ranked_vals

    # widths for nice alignment
    w_rank <- nchar(as.character(length(nm)))
    w_nm   <- if (length(nm)) max(nchar(nm)) else 0
    w_av   <- if (length(nm)) max(nchar(formatC(av, format = "f", digits = 3, big.mark = ","))) else 0

    # lines like: " 1) ModelName          AICc =   78.184"
    lines <- if (length(nm)) sprintf(
      paste0("%", w_rank, "d) %-", w_nm, "s  ", label, " = %", w_av, "s"),
      seq_along(nm), nm, formatC(av, format = "f", digits = 3, big.mark = ",")
    ) else character(0)

    out <- paste(lines, collapse = "\n")

    # add undefined list (if any), on a new line
    if (isTRUE(any(!finite))) {
      undef <- paste(model_names[!finite], collapse = ", ")
      note  <- paste0("-- undefined ", label, ": ", undef)
      out   <- if (nzchar(out)) paste(out, note, sep = "\n") else note
    }

    out
  })(ranked_names, ranked_vals, model_names, finite)

  cat(ranking_str, "\n")

  AICc_vals <- setNames(
    lapply(seq_along(res_list), function(i) {
      list(AICc = AICc_vec[i], k = res_list[[i]]$k, n = res_list[[i]]$n)
    }),
    model_names
  )

  list(AICc_vals = AICc_vals)
}




