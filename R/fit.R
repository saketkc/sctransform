# Fir NB regression models using different approaches

fit_poisson <- function(umi, model_str, data, theta_estimation_fun) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  dfr <- ncol(umi) - ncol(regressor_data)
  par_mat <- t(apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, TRUE)
    theta <- switch(theta_estimation_fun,
      'theta.ml' = as.numeric(x = theta.ml(y = y, mu = fit$fitted)),
      'theta.mm' = as.numeric(x = theta.mm(y = y, mu = fit$fitted, dfr = dfr)),
      stop('theta_estimation_fun ', theta_estimation_fun, ' unknown - only theta.ml and theta.mm supported at the moment')
    )
    return(c(theta, fit$coefficients))
  }))
  return(par_mat)
}

fit_qpoisson <- function(umi, model_str, data) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  par_mat <- t(apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, FALSE)
    return(c(fit$theta.guesstimate, fit$coefficients))
  }))
  return(par_mat)
}

fit_nb_theta_given <- function(umi, model_str, data, theta_given) {
  par_lst <- lapply(1:nrow(umi), function(j) {
    y <- umi[j, ]
    theta <- theta_given[j]
    fit2 <- 0
    try(fit2 <- glm(as.formula(model_str), data = data, family = negative.binomial(theta=theta)), silent=TRUE)
    if (inherits(x = fit2, what = 'numeric')) {
      return(c(theta, glm(as.formula(model_str), data = data, family = poisson)$coefficients))
    } else {
      return(c(theta, fit2$coefficients))
    }
  })
  return(do.call(rbind, par_lst))
}

fit_nb_fast <- function(umi, model_str, data, theta_estimation_fun) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  dfr <- ncol(umi) - ncol(regressor_data)
  par_mat <- apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, TRUE)
    theta <- switch(theta_estimation_fun,
                    'theta.ml' = as.numeric(x = theta.ml(y = y, mu = fit$fitted)),
                    'theta.mm' = as.numeric(x = theta.mm(y = y, mu = fit$fitted, dfr = dfr)),
                    stop('theta_estimation_fun ', theta_estimation_fun, ' unknown - only theta.ml and theta.mm supported at the moment')
    )
    fit2 <- 0
    try(fit2 <- glm(as.formula(model_str), data = data, family = negative.binomial(theta=theta)), silent=TRUE)
    if (inherits(x = fit2, what = 'numeric')) {
      return(c(theta, fit$coefficients))
    } else {
      return(c(theta, fit2$coefficients))
    }
  })
  return(t(par_mat))
}

fit_nb <- function(umi, model_str, data) {
  par_mat <- apply(umi, 1, function(y) {
    fit <- 0
    try(fit <- glm.nb(as.formula(model_str), data = data), silent=TRUE)
    if (inherits(x = fit, what = 'numeric')) {
      fit <- glm(as.formula(model_str), data = data, family = poisson)
      fit$theta <- as.numeric(x = theta.ml(y = y, mu = fit$fitted))
    }
    return(c(fit$theta, fit$coefficients))
  })
  return(t(par_mat))
}

fit_glmGamPoi <- function(umi, model_str, data) {
  fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(gsub("y", "", model_str)),
                           col_data = data,
                           size_factors = FALSE)
  fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-4)
  colnames(fit$Beta)[match(x = 'Intercept', colnames(fit$Beta))] <- "(Intercept)"
  return(cbind(fit$theta, fit$Beta))
}

fit_glmGamPoi2 <- function(umi, model_str, data) {
  fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(gsub("y", "", model_str)),
                           col_data = data,
                           size_factors = FALSE)
  fit$theta <- 1 / fit$overdispersions
  colnames(fit$Beta)[match(x = 'Intercept', colnames(fit$Beta))] <- "(Intercept)"
  return(cbind(fit$theta, fit$Beta))
}

fit_glmGamPoi6 <- function(umi, model_str, data) {
  # only intercept varies 
  # TODO does away with batch_var which needs to be handled
  new_formula <- "~ 1"
  log10_umi <- data$log_umi 
  stopifnot(!is.null(log10_umi))
  log_umi <- log(10^log10_umi)

  fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(new_formula),#gsub("y", "", model_str)),
                           col_data = data,
                           offset = log_umi,
                           size_factors = FALSE)
  fit$theta <- 1 / fit$overdispersions
  # TODO does away with batch_var which needs to be handled
  model_pars <- cbind(fit$theta,
                       fit$Beta[, "Intercept"],
                        rep(log(10), nrow(umi)))
  dimnames(model_pars) <- list(rownames(umi), c('theta', '(Intercept)', 'log_umi'))
  return(model_pars)
}


fit_glmGamPoi7 <- function(umi, mean){
  
  fit <- glmGamPoi::overdispersion_mle(umi,
                                       mean,
                                       model_matrix = NULL,
                                       # TODO: not sure if it makes a difference?
                                       do_cox_reid_adjustment = TRUE, #!is.null(model_matrix),
                                       global_estimate = FALSE,
                                       subsample = FALSE,
                                       max_iter = 200,
                                       verbose = FALSE
                                       )

  theta <- 1 / fit$estimate
  return (theta)
}

fit_glmGamPoi8 <- function(umi, mean){
  
  fit <- glmGamPoi::overdispersion_mle(umi,
                                       mean,
                                       model_matrix = NULL,
                                       # TODO: not sure if it makes a difference?
                                       do_cox_reid_adjustment = TRUE, #!is.null(model_matrix),
                                       global_estimate = FALSE,
                                       subsample = FALSE,
                                       max_iter = 200,
                                       verbose = FALSE
                                       )

  theta <- 1 / fit$estimate
  return (data.frame(theta=theta))
}
