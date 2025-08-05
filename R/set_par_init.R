set_par_init <- function(type, model_fixed_params = NA) {
  if (type == "MuenchGeneral") {
    return(c(k=0.5, l=1, foi=0.1))
  }

  else if (type == "MuenchRestricted") {
    return(c(foi=0.1))
  }

  else if (type == "Griffiths") {
    return(c(gamma0=0.1, gamma1=-5))
  }

  else if (type == "Farringtons") {
    return(c(gamma0=0.1, gamma1=1, gamma2=0.1))
  }

  else if (type == "PiecewiseConstant") {
    num_pieces <- length(model_fixed_params$upper_cutoffs)
    par_vals <- rep(0.1, num_pieces)
    names(par_vals) <- paste0("foi", seq_len(num_pieces))
    return(par_vals)
  }

  else {
    return(NA)
  }
}
