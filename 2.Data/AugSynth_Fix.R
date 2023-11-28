################################################################################
## Main functions for single-period treatment augmented synthetic controls Method
################################################################################


#' Fit Augmented SCM
#' 
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @param progfunc What function to use to impute control outcomes
#'                 ridge=Ridge regression (allows for standard errors),
#'                 none=No outcome model,
#'                 en=Elastic Net, RF=Random Forest, GSYN=gSynth,
#'                 mcp=MCPanel, 
#'                 cits=Comparitive Interuppted Time Series
#'                 causalimpact=Bayesian structural time series with CausalImpact
#' @param scm Whether the SCM weighting function is used
#' @param fixedeff Whether to include a unit fixed effect, default F 
#' @param cov_agg Covariate aggregation functions, if NULL then use mean with NAs omitted
#' @param ... optional arguments for outcome model
#'
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#'          \item{"data"}{Panel data as matrices}
#'         }
#' @export
single_augsynth <- function(form, unit, time, t_int, data,
                            progfunc = "ridge",
                            scm=T,
                            fixedeff = FALSE,
                            cov_agg=NULL, ...) {
  call_name <- match.call()
  
  form <- Formula::Formula(form)
  unit <- enquo(unit)
  time <- enquo(time)
  
  ## format data
  outcome <- terms(formula(form, rhs=1))[[2]]
  trt <- terms(formula(form, rhs=1))[[3]]
  
  wide <- format_data(outcome, trt, unit, time, t_int, data)
  synth_data <- do.call(format_synth, wide)
  
  treated_units <- data %>% filter(!!trt == 1) %>% distinct(!!unit) %>% pull(!!unit)
  control_units <- data %>% filter(!(!!unit %in% treated_units)) %>% 
    distinct(!!unit) %>% arrange(!!unit) %>% pull(!!unit)
  ## add covariates
  if(length(form)[2] == 2) {
    Z <- extract_covariates(form, unit, time, t_int, data, cov_agg)
  } else {
    Z <- NULL
  }
  
  # fit augmented SCM
  augsynth <- fit_augsynth_internal(wide, synth_data, Z, progfunc, 
                                    scm, fixedeff, ...)
  
  # add some extra data
  augsynth$data$time <- data %>% distinct(!!time) %>%
    arrange(!!time) %>% pull(!!time)
  augsynth$call <- call_name
  augsynth$t_int <- t_int 
  
  augsynth$weights <- matrix(augsynth$weights)
  rownames(augsynth$weights) <- control_units
  
  return(augsynth)
}


#' Internal function to fit augmented SCM
#' @param wide Data formatted from format_data
#' @param synth_data Data formatted from foramt_synth
#' @param Z Matrix of auxiliary covariates
#' @param progfunc outcome model to use
#' @param scm Whether to fit SCM
#' @param fixedeff Whether to de-mean synth
#' @param V V matrix for Synth, default NULL
#' @param ... Extra args for outcome model
#' 
#' @noRd
#' 
fit_augsynth_internal <- function(wide, synth_data, Z, progfunc,
                                  scm, fixedeff, V = NULL, ...) {
  
  n <- nrow(wide$X)
  t0 <- ncol(wide$X)
  ttot <- t0 + ncol(wide$y)
  if(fixedeff) {
    demeaned <- demean_data(wide, synth_data)
    fit_wide <- demeaned$wide
    fit_synth_data <- demeaned$synth_data
    mhat <- demeaned$mhat
  } else {
    fit_wide <- wide
    fit_synth_data <- synth_data
    mhat <- matrix(0, n, ttot)
  }
  if (is.null(progfunc)) {
    progfunc = "none"
  }
  progfunc = tolower(progfunc)
  ## fit augsynth
  if(progfunc == "ridge") {
    # Ridge ASCM
    augsynth <- do.call(fit_ridgeaug_formatted,
                        list(wide_data = fit_wide,
                             synth_data = fit_synth_data,
                             Z = Z, V = V, scm = scm, ...))
  } else if(progfunc == "none") {
    ## Just SCM
    augsynth <- do.call(fit_ridgeaug_formatted,
                        c(list(wide_data = fit_wide, 
                               synth_data = fit_synth_data,
                               Z = Z, ridge = F, scm = T, V = V, ...)))
  } else {
    ## Other outcome models
    progfuncs = c("ridge", "none", "en", "rf", "gsyn", "mcp",
                  "cits", "causalimpact", "seq2seq")
    if (progfunc %in% progfuncs) {
      augsynth <- fit_augsyn(fit_wide, fit_synth_data, 
                             progfunc, scm, ...)
    } else {
      stop("progfunc must be one of 'EN', 'RF', 'GSYN', 'MCP', 'CITS', 'CausalImpact', 'seq2seq', 'None'")
    }
    
  }
  
  augsynth$mhat <- mhat + cbind(matrix(0, nrow = n, ncol = t0), 
                                augsynth$mhat)
  augsynth$data <- wide
  augsynth$data$Z <- Z
  augsynth$data$synth_data <- synth_data
  augsynth$progfunc <- progfunc
  augsynth$scm <- scm
  augsynth$fixedeff <- fixedeff
  augsynth$extra_args <- list(...)
  if(progfunc == "ridge") {
    augsynth$extra_args$lambda <- augsynth$lambda
  } else if(progfunc == "gsyn") {
    augsynth$extra_args$r <- ncol(augsynth$params$factor)
    augsynth$extra_args$CV <- 0
  }
  ##format output
  class(augsynth) <- "augsynth"
  return(augsynth)
}

#' Get prediction of ATT or average outcome under control
#' @param object augsynth object
#' @param att If TRUE, return the ATT, if FALSE, return imputed counterfactual
#' @param ... Optional arguments
#'
#' @return Vector of predicted post-treatment control averages
#' @export
predict.augsynth <- function(object, att = F, ...) {
  # if ("att" %in% names(list(...))) {
  #     att <- list(...)$att
  # } else {
  #     att <- F
  # }
  augsynth <- object
  
  X <- augsynth$data$X
  y <- augsynth$data$y
  comb <- cbind(X, y)
  trt <- augsynth$data$trt
  mhat <- augsynth$mhat
  
  m1 <- colMeans(mhat[trt==1,,drop=F])
  
  resid <- (comb[trt==0,,drop=F] - mhat[trt==0,drop=F])
  
  y0 <- m1 + t(resid) %*% augsynth$weights
  if(att) {
    return(colMeans(comb[trt == 1,, drop = F]) - c(y0))
  } else {
    rnames <- rownames(y0)
    y0_vec <- c(y0)
    names(y0_vec) <- rnames
    return(y0_vec)
  }
}


#' Print function for augsynth
#' @param x augsynth object
#' @param ... Optional arguments
#' @export
print.augsynth <- function(x, ...) {
  augsynth <- x
  
  ## straight from lm
  cat("\nCall:\n", paste(deparse(augsynth$call), sep="\n", collapse="\n"), "\n\n", sep="")
  
  ## print att estimates
  tint <- ncol(augsynth$data$X)
  ttotal <- tint + ncol(augsynth$data$y)
  att_post <- predict(augsynth, att = T)[(tint + 1):ttotal]
  
  cat(paste("Average ATT Estimate: ",
            format(round(mean(att_post),3), nsmall = 3), "\n\n", sep=""))
}


#' Plot function for augsynth
#' @importFrom graphics plot
#' 
#' @param x Augsynth object to be plotted
#' @param inf Boolean, whether to get confidence intervals around the point estimates
#' @param cv If True, plot cross validation MSE against hyper-parameter, otherwise plot effects
#' @param ... Optional arguments
#' @export
plot.augsynth <- function(x, inf = T, cv = F, ...) {
  # if ("se" %in% names(list(...))) {
  #     se <- list(...)$se
  # } else {
  #     se <- T
  # }
  
  augsynth <- x
  
  if (cv == T) {
    errors = data.frame(lambdas = augsynth$lambdas,
                        errors = augsynth$lambda_errors,
                        errors_se = augsynth$lambda_errors_se)
    p <- ggplot2::ggplot(errors, ggplot2::aes(x = lambdas, y = errors)) +
      ggplot2::geom_point(size = 2) + 
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = errors,
                     ymax = errors + errors_se),
        width=0.2, size = 0.5) 
    p <- p + ggplot2::labs(title = bquote("Cross Validation MSE over " ~ lambda),
                           x = expression(lambda), y = "Cross Validation MSE", 
                           parse = TRUE)
    p <- p + ggplot2::scale_x_log10()
    
    # find minimum and min + 1se lambda to plot
    min_lambda <- choose_lambda(augsynth$lambdas,
                                augsynth$lambda_errors,
                                augsynth$lambda_errors_se,
                                F)
    min_1se_lambda <- choose_lambda(augsynth$lambdas,
                                    augsynth$lambda_errors,
                                    augsynth$lambda_errors_se,
                                    T)
    min_lambda_index <- which(augsynth$lambdas == min_lambda)
    min_1se_lambda_index <- which(augsynth$lambdas == min_1se_lambda)
    
    p <- p + ggplot2::geom_point(
      ggplot2::aes(x = min_lambda, 
                   y = augsynth$lambda_errors[min_lambda_index]),
      color = "gold")
    p + ggplot2::geom_point(
      ggplot2::aes(x = min_1se_lambda,
                   y = augsynth$lambda_errors[min_1se_lambda_index]),
      color = "gold") +
      ggplot2::theme_bw()
  } else {
    plot(summary(augsynth, ...), inf = inf)
  }
}


#' Summary function for augsynth
#' @param object augsynth object
#' @param inf Boolean, whether to get confidence intervals around the point estimates
#' @param inf_type Type of inference algorithm. Options are
#'         \itemize{
#'          \item{"conformal"}{Conformal inference (default)}
#'          \item{"jackknife+"}{Jackknife+ algorithm over time periods}
#'          \item{"jackknife"}{Jackknife over units}
#'         }
#' @param ... Optional arguments for inference, for more details for each `inf_type` see
#'         \itemize{
#'          \item{"conformal"}{`conformal_inf`}
#'          \item{"jackknife+"}{`time_jackknife_plus`}
#'          \item{"jackknife"}{`jackknife_se_single`}
#'         }
#' @export
summary.augsynth <- function(object, inf = T, inf_type = "conformal", ...) {
  augsynth <- object
  # if ("inf" %in% names(list(...))) {
  #     inf <- list(...)$inf
  # } else {
  #     inf <- T
  # }
  # if ("inf_type" %in% names(list(...))) {
  #     inf_type <- list(...)$inf_type
  # } else {
  #     inf_type <- "conformal"
  # }
  
  
  summ <- list()
  
  t0 <- ncol(augsynth$data$X)
  t_final <- t0 + ncol(augsynth$data$y)
  
  if(inf) {
    if(inf_type == "jackknife") {
      att_se <- jackknife_se_single(augsynth)
    } else if(inf_type == "jackknife+") {
      att_se <- time_jackknife_plus(augsynth, ...)
    } else if(inf_type == "conformal") {
      att_se <- conformal_inf(augsynth, ...)
    } else {
      stop(paste(inf_type, "is not a valid choice of 'inf_type'"))
    }
    
    att <- data.frame(Time = augsynth$data$time,
                      Estimate = att_se$att[1:t_final])
    if(inf_type == "jackknife") {
      att$Std.Error <- att_se$se[1:t_final]
      att_avg_se <- att_se$se[t_final + 1]
    } else {
      att_avg_se <- NA
    }
    att_avg <- att_se$att[t_final + 1]
    if(inf_type %in% c("jackknife+", "nonpar_bs", "t_dist", "conformal")) {
      att$lower_bound <- att_se$lb[1:t_final]
      att$upper_bound <- att_se$ub[1:t_final]
    }
    if(inf_type == "conformal") {
      att$p_val <- att_se$p_val[1:t_final]
    }
    
  } else {
    t0 <- ncol(augsynth$data$X)
    t_final <- t0 + ncol(augsynth$data$y)
    att_est <- predict(augsynth, att = T)
    att <- data.frame(Time = augsynth$data$time,
                      Estimate = att_est)
    att$Std.Error <- NA
    att_avg <- mean(att_est[(t0 + 1):t_final])
    att_avg_se <- NA
  }
  
  summ$att <- att
  summ$average_att <- data.frame(Estimate = att_avg, Std.Error = att_avg_se)
  if(inf) {
    if(inf_type %in% c("jackknife+", "conformal")) {
      summ$average_att$lower_bound <- att_se$lb[t_final + 1]
      summ$average_att$upper_bound <- att_se$ub[t_final + 1]
      summ$alpha <-  att_se$alpha
    }
    if(inf_type == "conformal") {
      summ$average_att$p_val <- att_se$p_val[t_final + 1]
    }
  }
  summ$t_int <- augsynth$t_int
  summ$call <- augsynth$call
  summ$l2_imbalance <- augsynth$l2_imbalance
  summ$scaled_l2_imbalance <- augsynth$scaled_l2_imbalance
  if(!is.null(augsynth$covariate_l2_imbalance)) {
    summ$covariate_l2_imbalance <- augsynth$covariate_l2_imbalance
    summ$scaled_covariate_l2_imbalance <- augsynth$scaled_covariate_l2_imbalance
  }
  ## get estimated bias
  
  if(tolower(augsynth$progfunc) == "ridge") {
    mhat <- augsynth$ridge_mhat
    w <- augsynth$synw
  } else {
    mhat <- augsynth$mhat
    w <- augsynth$weights
  }
  trt <- augsynth$data$trt
  m1 <- colMeans(mhat[trt==1,,drop=F])
  
  if(tolower(augsynth$progfunc) == "none" | (!augsynth$scm)) {
    summ$bias_est <- NA
  } else {
    summ$bias_est <- m1 - t(mhat[trt==0,,drop=F]) %*% w
  }
  
  
  summ$inf_type <- if(inf) inf_type else "None"
  class(summ) <- "summary.augsynth"
  return(summ)
}

#' Print function for summary function for augsynth
#' @param x summary object
#' @param ... Optional arguments
#' @export
print.summary.augsynth <- function(x, ...) {
  summ <- x
  
  ## straight from lm
  cat("\nCall:\n", paste(deparse(summ$call), sep="\n", collapse="\n"), "\n\n", sep="")
  
  t_final <- nrow(summ$att)
  
  ## distinction between pre and post treatment
  att_est <- summ$att$Estimate
  t_total <- length(att_est)
  t_int <- summ$att %>% filter(Time <= summ$t_int) %>% nrow()
  
  att_pre <- att_est[1:(t_int-1)]
  att_post <- att_est[t_int:t_total]
  
  
  out_msg <- ""
  
  
  # print out average post treatment estimate
  att_post <- summ$average_att$Estimate
  se_est <- summ$att$Std.Error
  if(summ$inf_type == "jackknife") {
    se_avg <- summ$average_att$Std.Error
    
    out_msg <- paste("Average ATT Estimate (Jackknife Std. Error): ",
                     format(round(att_post,3), nsmall=3), 
                     "  (",
                     format(round(se_avg,3)), ")\n")
    inf_type <- "Jackknife over units"
  } else if(summ$inf_type == "conformal") {
    p_val <- summ$average_att$p_val
    out_msg <- paste("Average ATT Estimate (p Value for Joint Null): ",
                     format(round(att_post,3), nsmall=3), 
                     "  (",
                     format(round(p_val,3)), ")\n")
    inf_type <- "Conformal inference"
  } else if(summ$inf_type == "jackknife+") {
    out_msg <- paste("Average ATT Estimate: ",
                     format(round(att_post,3), nsmall=3), "\n")
    inf_type <- "Jackknife+ over time periods"
  } else {
    out_msg <- paste("Average ATT Estimate: ",
                     format(round(att_post,3), nsmall=3), "\n")
    inf_type <- "None"
  }
  
  
  out_msg <- paste(out_msg, 
                   "L2 Imbalance: ",
                   format(round(summ$l2_imbalance,3), nsmall=3), "\n",
                   "Percent improvement from uniform weights: ",
                   format(round(1 - summ$scaled_l2_imbalance,3)*100), "%\n\n",
                   sep="")
  if(!is.null(summ$covariate_l2_imbalance)) {
    
    out_msg <- paste(out_msg,
                     "Covariate L2 Imbalance: ",
                     format(round(summ$covariate_l2_imbalance,3), 
                            nsmall=3),
                     "\n",
                     "Percent improvement from uniform weights: ",
                     format(round(1 - summ$scaled_covariate_l2_imbalance,3)*100), 
                     "%\n\n",
                     sep="")
    
  }
  out_msg <- paste(out_msg, 
                   "Avg Estimated Bias: ",
                   format(round(mean(summ$bias_est), 3),nsmall=3), "\n\n",
                   "Inference type: ",
                   inf_type,
                   "\n\n",
                   sep="")
  cat(out_msg)
  
  if(summ$inf_type == "jackknife") {
    out_att <- summ$att[t_int:t_final,] %>% 
      select(Time, Estimate, Std.Error)
  } else if(summ$inf_type == "conformal") {
    out_att <- summ$att[t_int:t_final,] %>% 
      select(Time, Estimate, lower_bound, upper_bound, p_val)
    names(out_att) <- c("Time", "Estimate", 
                        paste0((1 - summ$alpha) * 100, "% CI Lower Bound"),
                        paste0((1 - summ$alpha) * 100, "% CI Upper Bound"),
                        paste0("p Value"))
  } else if(summ$inf_type == "jackknife+") {
    out_att <- summ$att[t_int:t_final,] %>% 
      select(Time, Estimate, lower_bound, upper_bound)
    names(out_att) <- c("Time", "Estimate", 
                        paste0((1 - summ$alpha) * 100, "% CI Lower Bound"),
                        paste0((1 - summ$alpha) * 100, "% CI Upper Bound"))
  } else {
    out_att <- summ$att[t_int:t_final,] %>% 
      select(Time, Estimate)
  }
  out_att %>%
    mutate_at(vars(-Time), ~ round(., 3)) %>%
    print(row.names = F)
  
  
}

#' Plot function for summary function for augsynth
#' @param x Summary object
#' @param inf Boolean, whether to plot confidence intervals
#' @param ... Optional arguments
#' @export
plot.summary.augsynth <- function(x, inf = T, ...) {
  summ <- x
  # if ("inf" %in% names(list(...))) {
  #     inf <- list(...)$inf
  # } else {
  #     inf <- T
  # }
  
  p <- summ$att %>%
    ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate))
  if(inf) {
    if(all(is.na(summ$att$lower_bound))) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                                                 ymax=Estimate+2*Std.Error),
                                    alpha=0.2)
    } else {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=lower_bound,
                                                 ymax=upper_bound),
                                    alpha=0.2)
    }
    
  }
  p + ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept=summ$t_int, lty=2) +
    ggplot2::geom_hline(yintercept=0, lty=2) + 
    ggplot2::theme_bw()
  
}





################################################################################
## Code for inference
################################################################################

#' Jackknife+ algorithm over time
#' @param ascm Fitted `augsynth` object
#' @param alpha Confidence level
#' @param conservative Whether to use the conservative jackknife+ procedure
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"heldout_att"}{Vector of ATT estimates with the time period held out}
#'          \item{"se"}{Standard error, always NA but returned for compatibility}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
time_jackknife_plus <- function(ascm, alpha = 0.05, conservative = F) {
  wide_data <- ascm$data
  synth_data <- ascm$data$synth_data
  n <- nrow(wide_data$X)
  n_c <- dim(synth_data$Z0)[2]
  Z <- wide_data$Z
  
  t0 <- dim(synth_data$Z0)[1]
  tpost <- ncol(wide_data$y)
  t_final <- dim(synth_data$Y0plot)[1]
  
  jack_ests <- lapply(1:t0, 
                      function(tdrop) {
                        # drop unit i
                        new_data <- drop_time_t(wide_data, Z, tdrop)
                        # refit
                        new_ascm <- do.call(fit_augsynth_internal,
                                            c(list(wide = new_data$wide,
                                                   synth_data = new_data$synth_data,
                                                   Z = new_data$Z,
                                                   progfunc = ascm$progfunc,
                                                   scm = ascm$scm,
                                                   fixedeff = ascm$fixedeff),
                                              ascm$extra_args))
                        # get ATT estimates and held out error for time t
                        # t0 is prediction for held out time
                        est <- predict(new_ascm, att = F)[(t0 +1):t_final]
                        est <- c(est, mean(est))
                        err <- c(colMeans(wide_data$X[wide_data$trt == 1,
                                                      tdrop,
                                                      drop = F]) -
                                   predict(new_ascm, att = F)[t0])
                        list(err, rbind(est + abs(err), est - abs(err), est + err, est))
                      })
  # get errors and jackknife distribution
  held_out_errs <- vapply(jack_ests, `[[`, numeric(1), 1)
  jack_dist <- vapply(jack_ests, `[[`,
                      matrix(0, nrow = 4, ncol = tpost + 1), 2)
  
  out <- list()
  att <- predict(ascm, att = T)
  out$att <- c(att, 
               mean(att[(t0 + 1):t_final]))
  # held out ATT
  out$heldout_att <- c(held_out_errs, 
                       att[(t0 + 1):t_final], 
                       mean(att[(t0 + 1):t_final]))
  
  # out$se <- rep(NA, 10 + tpost)
  if(conservative) {
    qerr <- stats::quantile(abs(held_out_errs), 1 - alpha)
    out$lb <- c(rep(NA, t0), apply(jack_dist[4,,], 1, min) - qerr)
    out$ub <- c(rep(NA, t0), apply(jack_dist[4,,], 1, max) + qerr)
  } else {
    out$lb <- c(rep(NA, t0), apply(jack_dist[2,,], 1, stats::quantile, alpha / 2))
    out$ub <- c(rep(NA, t0), apply(jack_dist[1,,], 1, stats::quantile, 1 - alpha / 2))
  }
  # shift back to ATT scale
  y1 <- predict(ascm, att = F) + att
  y1 <-  c(y1, mean(y1[(t0 + 1):t_final]))
  shifted_lb <- y1 - out$ub
  shifted_ub <- y1 - out$lb
  out$lb <- shifted_lb
  out$ub <- shifted_ub
  out$alpha <- alpha
  
  
  return(out)
}

#' Drop time period from pre-treatment data
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
#' @param t_drop Time to drop
#' @noRd
drop_time_t <- function(wide_data, Z, t_drop) {
  
  new_wide_data <- list()
  new_wide_data$trt <- wide_data$trt
  new_wide_data$X <- wide_data$X[, -t_drop, drop = F]
  new_wide_data$y <- cbind(wide_data$X[, t_drop, drop = F], 
                           wide_data$y)
  
  X0 <- new_wide_data$X[new_wide_data$trt == 0,, drop = F]
  x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1,,
                                        drop = F]),
               ncol=1)
  y0 <- new_wide_data$y[new_wide_data$trt == 0,, drop = F]
  y1 <- colMeans(new_wide_data$y[new_wide_data$trt == 1,, drop = F])
  
  new_synth_data <- list()
  new_synth_data$Z0 <- t(X0)
  new_synth_data$X0 <- t(X0)
  new_synth_data$Z1 <- x1
  new_synth_data$X1 <- x1
  
  return(list(wide_data = new_wide_data,
              synth_data = new_synth_data,
              Z = Z)) 
}

#' Conformal inference procedure to compute p-values and point-wise confidence intervals
#' @param ascm Fitted `augsynth` object
#' @param alpha Confidence level
#' @param stat_func Function to compute test statistic
#' @param type Either "iid" for iid permutations or "block" for moving block permutations; default is "block"
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param grid_size Number of grid points to use when inverting the hypothesis test
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"heldout_att"}{Vector of ATT estimates with the time period held out}
#'          \item{"se"}{Standard error, always NA but returned for compatibility}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"p_val"}{p-value for test of no post-treatment effect}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
conformal_inf <- function(ascm, alpha = 0.05, 
                          stat_func = NULL, type = "block",
                          q = 1, ns = 1000, grid_size = 200) {
  wide_data <- ascm$data
  synth_data <- ascm$data$synth_data
  n <- nrow(wide_data$X)
  n_c <- dim(synth_data$Z0)[2]
  Z <- wide_data$Z
  
  t0 <- dim(synth_data$Z0)[1]
  tpost <- ncol(wide_data$y)
  t_final <- dim(synth_data$Y0plot)[1]
  
  # grid of nulls
  att <- predict(ascm, att = T)
  post_att <- att[(t0 +1):t_final]
  post_sd <- sqrt(mean(post_att ^ 2))
  # iterate over post-treatment periods to get pointwise CIs
  vapply(1:tpost,
         function(j) {
           # fit using t0 + j as a pre-treatment period and get reisduals
           new_wide_data <- wide_data
           new_wide_data$X <- cbind(wide_data$X, wide_data$y[, j, drop = TRUE])
           if(tpost > 1) {
             new_wide_data$y <- wide_data$y[, -j, drop = FALSE]
           } else {
             # set the post period has to be *something*
             new_wide_data$y <- matrix(1, nrow = n, ncol = 1)
           }
           
           
           # make a grid around the estimated ATT
           grid <- seq(att[t0 + j] - 1 * post_sd, att[t0 + j] + 1 * post_sd,
                       length.out = grid_size)
           compute_permute_ci(new_wide_data, ascm, grid, 1, alpha, type,
                              q, ns, stat_func)
         },
         numeric(3)) -> cis
  
  # test a null post-treatment effect
  new_wide_data <- wide_data
  new_wide_data$X <- cbind(wide_data$X, wide_data$y)
  new_wide_data$y <- matrix(1, nrow = n, ncol = 1)
  null_p <- compute_permute_pval(new_wide_data, ascm, 0, ncol(wide_data$y), 
                                 type, q, ns, stat_func)
  
  out <- list()
  att <- predict(ascm, att = T)
  out$att <- c(att, mean(att[(t0 + 1):t_final]))
  # out$se <- rep(NA, t_final)
  # out$sigma <- NA
  out$lb <- c(rep(NA, t0), cis[1, ], NA)
  out$ub <- c(rep(NA, t0), cis[2, ], NA)
  out$p_val <- c(rep(NA, t0), cis[3, ], null_p)
  out$alpha <- alpha
  return(out)
}

#' Compute conformal test statistics
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param h0 Null hypothesis to test
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param stat_func Function to compute test statistic
#' 
#' @return List that contains:
#'         \itemize{
#'          \item{"resids"}{Residuals after enforcing the null}
#'          \item{"test_stats"}{Permutation distribution of test statistics}
#'          \item{"stat_func"}{Test statistic function}
#'         }
#' @noRd
compute_permute_test_stats <- function(wide_data, ascm, h0,
                                       post_length, type,
                                       q, ns, stat_func) {
  # format data
  new_wide_data <- wide_data
  t0 <- ncol(wide_data$X) - post_length
  tpost <- t0 + post_length
  # adjust outcomes for null
  new_wide_data$X[wide_data$trt == 1,(t0 + 1):tpost ] <- new_wide_data$X[wide_data$trt == 1,(t0 + 1):tpost] - h0
  X0 <- new_wide_data$X[new_wide_data$trt == 0,, drop = F]
  x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1,, drop = F]),
               ncol=1)
  
  new_synth_data <- list()
  new_synth_data$Z0 <- t(X0)
  new_synth_data$X0 <- t(X0)
  new_synth_data$Z1 <- x1
  new_synth_data$X1 <- x1
  
  # fit synth with adjusted data and get residuals
  new_ascm <- do.call(fit_augsynth_internal,
                      c(list(wide = new_wide_data,
                             synth_data = new_synth_data,
                             Z = wide_data$Z,
                             progfunc = ascm$progfunc,
                             scm = ascm$scm,
                             fixedeff = ascm$fixedeff),
                        ascm$extra_args))
  resids <- predict(new_ascm, att = T)[1:tpost]
  # permute residuals and compute test statistic
  if(is.null(stat_func)) {
    stat_func <- function(x) (sum(abs(x) ^ q)  / sqrt(length(x))) ^ (1 / q)
  }
  if(type == "iid") {
    test_stats <- sapply(1:ns, 
                         function(x) {
                           reorder <- sample(resids)
                           stat_func(reorder[(t0 + 1):tpost])
                         })
  } else {
    ## increment time by one step and wrap
    test_stats <- sapply(1:tpost,
                         function(j) {
                           reorder <- resids[(0:tpost -1 + j) %% tpost + 1]
                           stat_func(reorder[(t0 + 1):tpost])
                         })
  }
  
  return(list(resids = resids,
              test_stats = test_stats,
              stat_func = stat_func))
}


#' Compute conformal p-value
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param h0 Null hypothesis to test
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param stat_func Function to compute test statistic
#' 
#' @return Computed p-value
#' @noRd
compute_permute_pval <- function(wide_data, ascm, h0,
                                 post_length, type,
                                 q, ns, stat_func) {
  t0 <- ncol(wide_data$X) - post_length
  tpost <- t0 + post_length
  out <- compute_permute_test_stats(wide_data, ascm, h0,
                                    post_length, type, q, ns, stat_func)
  mean(out$stat_func(out$resids[(t0 + 1):tpost]) <= out$test_stats)
}

#' Compute conformal p-value
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param grid Set of null hypothesis to test for inversion
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param stat_func Function to compute test statistic
#' 
#' @return (lower bound of interval, upper bound of interval, p-value for null of 0 effect)
#' @noRd
compute_permute_ci <- function(wide_data, ascm, grid,
                               post_length, alpha, type,
                               q, ns, stat_func) {
  
  grid <- c(grid, 0)
  ps <- sapply(grid, function(x) {
    compute_permute_pval(wide_data, ascm, x, post_length, 
                         type, q, ns, stat_func)
  })
  c(min(grid[ps >= alpha]), max(grid[ps >= alpha]), ps[grid == 0])
}


#' Jackknife+ algorithm over time
#' @param ascm Fitted `augsynth` object
#' @param alpha Confidence level
#' @param conservative Whether to use the conservative jackknife+ procedure
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"heldout_att"}{Vector of ATT estimates with the time period held out}
#'          \item{"se"}{Standard error, always NA but returned for compatibility}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
time_jackknife_plus_multiout <- function(ascm_multi, alpha = 0.05, conservative = F) {
  wide_data <- ascm_multi$data
  data_list <- ascm_multi$data_list
  
  n <- nrow(wide_data$X)
  k <- length(data_list$X)
  
  
  t0 <- min(sapply(data_list$X, ncol))
  tpost <- max(sapply(data_list$y, ncol))
  t_final <- t0 + tpost
  Z <- wide_data$Z
  
  jack_ests <- lapply(1:t0, 
                      function(tdrop) {
                        # drop unit i
                        new_data_list <- drop_time_t_multiout(data_list, Z, tdrop)
                        # refit
                        new_ascm <- do.call(fit_augsynth_multiout_internal,
                                            c(list(wide_list = new_data_list,
                                                   combine_method = ascm_multi$combine_method,
                                                   Z = data_list$Z,
                                                   progfunc = ascm_multi$progfunc,
                                                   scm = ascm_multi$scm,
                                                   fixedeff = ascm_multi$fixedeff,
                                                   outcomes_str = ascm_multi$outcomes),
                                              ascm_multi$extra_args))
                        # get ATT estimates and held out error for time t
                        # t0 is prediction for held out time
                        est <- predict(new_ascm, att = F)[(t0 +1):t_final, , drop = F]
                        est <- rbind(est, colMeans(est))
                        # err <- c(colMeans(wide_data$X[wide_data$trt == 1,
                        #                              tdrop,
                        #                              drop = F]) -
                        #         predict(new_ascm, att = F)[t0])
                        err <- c(predict(new_ascm, att = T)[t0, , drop = F])
                        list(err, t(t(est) + abs(err)), t(t(est) - abs(err)), t(t(est) + err), est)
                      })
  # get errors and jackknife distribution
  held_out_errs <- matrix(vapply(jack_ests, `[[`, numeric(k), 1), nrow = k)
  jack_dist_high <- vapply(jack_ests, `[[`,
                           matrix(0, nrow = tpost + 1, ncol = k), 2)
  jack_dist_low <- vapply(jack_ests, `[[`,
                          matrix(0, nrow = tpost + 1, ncol = k), 3)
  jack_dist_cons <- vapply(jack_ests, `[[`,
                           matrix(0, nrow = tpost + 1, ncol = k), 4)
  
  out <- list()
  att <- predict(ascm_multi, att = T)
  out$att <- rbind(att, 
                   colMeans(att[(t0 + 1):t_final, , drop = F]))
  # held out ATT
  
  out$heldout_att <- rbind(t(held_out_errs), 
                           att[(t0 + 1):t_final, , drop = F], 
                           colMeans(att[(t0 + 1):t_final, , drop = F]))
  if(conservative) {
    qerr <- apply(abs(held_out_errs), 1, 
                  stats::quantile, 1 - alpha, type = 1)
    out$lb <- rbind(matrix(NA, nrow = t0, ncol = k),
                    t(t(apply(jack_dist_cons, 1:2, min)) - qerr))
    out$ub <- rbind(matrix(NA, nrow = t0, ncol = k),
                    t(t(apply(jack_dist_cons, 1:2, max)) + qerr))
    
  } else {
    out$lb <- rbind(matrix(NA, nrow = t0, ncol = k),
                    apply(jack_dist_low, 1:2,
                          stats::quantile, alpha, type = 1))
    out$ub <- rbind(matrix(NA, nrow = t0, ncol = k), 
                    apply(jack_dist_high, 1:2, 
                          stats::quantile, 1 - alpha, type = 1))
  }
  # shift back to ATT scale
  y1 <- predict(ascm_multi, att = F) + att
  y1 <-  rbind(y1, colMeans(y1[(t0 + 1):t_final, , drop = F]))
  shifted_lb <- y1 - out$ub
  shifted_ub <- y1 - out$lb
  out$lb <- shifted_lb
  out$ub <- shifted_ub
  out$alpha <- alpha
  
  
  return(out)
}

#' Drop time period from pre-treatment data
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
#' @param t_drop Time to drop
#' @noRd
drop_time_t_multiout <- function(data_list, Z, t_drop) {
  
  new_data_list <- list()
  new_data_list$trt <- data_list$trt
  new_data_list$X <- lapply(data_list$X,
                            function(x) x[, -t_drop, drop = F])
  new_data_list$y <- lapply(1:length(data_list$y),
                            function(k) {
                              cbind(data_list$X[[k]][, t_drop, drop = F], 
                                    data_list$y[[k]])
                            })
  return(new_data_list)
}


#' Conformal inference procedure to compute p-values and point-wise confidence intervals
#' @param ascm Fitted `augsynth` object
#' @param alpha Confidence level
#' @param stat_func Function to compute test statistic
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param grid_size Number of grid points to use when inverting the hypothesis test (default is 1, so only to test joint null)
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"heldout_att"}{Vector of ATT estimates with the time period held out}
#'          \item{"se"}{Standard error, always NA but returned for compatibility}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"p_val"}{p-value for test of no post-treatment effect}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
conformal_inf_multiout <- function(ascm_multi, alpha = 0.05, 
                                   stat_func = NULL, type = "iid",
                                   q = 1, ns = 1000, grid_size = 1,
                                   lin_h0 = NULL) {
  wide_data <- ascm_multi$data
  data_list <- ascm_multi$data_list
  
  n <- nrow(wide_data$X)
  k <- length(data_list$X)
  
  
  t0 <- min(sapply(data_list$X, ncol))
  tpost <- max(sapply(data_list$y, ncol))
  t_final <- t0 + tpost
  
  # grid of nulls
  att <- predict(ascm_multi, att = T)
  post_att <- att[(t0 +1):t_final,, drop = F]
  post_sd <- apply(post_att, 2, function(x) sqrt(mean(x ^ 2, na.rm = T)))
  # iterate over post-treatment periods to get pointwise CIs
  vapply(1:tpost,
         function(j) {
           # fit using t0 + j as a pre-treatment period and get residuals
           new_data_list <- data_list
           new_data_list$X <- lapply(1:k,
                                     function(i) {
                                       Xi <- cbind(data_list$X[[i]], data_list$y[[i]][, j, drop = TRUE])
                                       colnames(Xi) <- c(colnames(data_list$X[[i]]),
                                                         colnames(data_list$y[[i]])[j])
                                       Xi
                                     })
           
           
           if(tpost > 1) {
             new_data_list$y <- lapply(1:k,
                                       function(i) {
                                         data_list$y[[i]][, -j, drop = FALSE]
                                       })
           } else {
             # set the post period has to be *something*
             new_data_list$y <- lapply(1:k,
                                       function(i) {
                                         x <- matrix(1, nrow = n, ncol = 1)
                                         colnames(x) <- max(as.numeric(colnames(data_list$y[[i]]))) + 1
                                         x
                                       })
           }
           
           
           # make a grid around the estimated ATT
           if(is.null(lin_h0)) {
             grid <- lapply(1:k, 
                            function(i) {
                              seq(att[t0 + j, i] - 2 * post_sd[i], att[t0 + j, i] + 2 * post_sd[i],
                                  length.out = grid_size)
                            })
           } else {
             grid <- seq(min(att[t0 + j, ]) - 2 * max(post_sd),
                         max(att[t0 + j, ]) + 2 * max(post_sd),
                         length.out = grid_size)
           }
           compute_permute_ci_multiout(new_data_list, ascm_multi, grid, 1,
                                       alpha, type, q, ns, lin_h0, stat_func)
         },
         matrix(0, ncol = k, nrow=3)) -> cis
  # # test a null post-treatment effect
  
  new_data_list <- data_list
  new_data_list$X <- lapply(1:k,
                            function(i) {
                              Xi <- cbind(data_list$X[[i]], data_list$y[[i]])
                              colnames(Xi) <- c(colnames(data_list$X[[i]]),
                                                colnames(data_list$y[[i]]))
                              Xi
                            })
  # set post treatment to be *something*
  new_data_list$y <- lapply(1:k,
                            function(i) {
                              data_list$y[[i]][, 1, drop = FALSE]
                            })
  null_p <- compute_permute_pval_multiout(new_data_list, ascm_multi,
                                          numeric(k), 
                                          tpost, type, q, ns, stat_func)
  if(is.null(lin_h0)) {
    grid <- lapply(1:k, 
                   function(i) {
                     seq(min(att[(t0 + 1):tpost, i]) - 4 * post_sd[i],
                         max(att[(t0 + 1):tpost, i]) + 4 * post_sd[i],
                         length.out = grid_size)
                   })
  } else {
    grid <- seq(min(att[t0 + 1, ]) - 3 * max(post_sd),
                max(att[t0 + 1, ]) + 3 * max(post_sd),
                length.out = grid_size)
  }
  null_ci <- compute_permute_ci_multiout(new_data_list, ascm_multi, grid,
                                         tpost, alpha, type, q, ns,
                                         lin_h0, stat_func)
  out <- list()
  att <- predict(ascm_multi, att = T)
  out$att <- rbind(att, apply(att[(t0 + 1):t_final, , drop = F], 2, mean))
  out$lb <- rbind(matrix(NA, nrow = t0, ncol = k),
                  t(matrix(cis[1, ,], nrow = k)),
                  # rep(NA, k)
                  null_ci[1,]
  )
  colnames(out$lb) <- ascm_multi$outcomes
  out$ub <- rbind(matrix(NA, nrow = t0, ncol = k),
                  t(matrix(cis[2, ,], nrow = k)),
                  # rep(NA, k)
                  null_ci[2,]
  )
  colnames(out$ub) <- ascm_multi$outcomes
  out$p_val <- rbind(matrix(NA, nrow = t0, ncol = k),
                     t(matrix(cis[3, ,], nrow = k)),
                     # rep(null_p, k)
                     null_ci[3,])
  colnames(out$p_val) <- ascm_multi$outcomes
  out$alpha <- alpha
  return(out)
}



#' Compute conformal test statistics
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param h0 Null hypothesis to test
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param stat_func Function to compute test statistic
#' 
#' @return List that contains:
#'         \itemize{
#'          \item{"resids"}{Residuals after enforcing the null}
#'          \item{"test_stats"}{Permutation distribution of test statistics}
#'          \item{"stat_func"}{Test statistic function}
#'         }
#' @noRd
compute_permute_test_stats_multiout <- function(data_list, ascm_multi, h0,
                                                post_length, type,
                                                q, ns, stat_func) {
  # format data
  new_data_list <- data_list
  t0 <- ncol(data_list$X[[1]]) - post_length
  tpost <- t0 + post_length
  k <- length(data_list$X)
  # adjust outcomes for null
  for(i in 1:k) {
    new_data_list$X[[k]][data_list$trt == 1,(t0 + 1):tpost ] <- new_data_list$X[[k]][data_list$trt == 1,(t0 + 1):tpost] - h0[i]
  }
  # fit synth with adjusted data and get residuals
  new_ascm <- do.call(fit_augsynth_multiout_internal,
                      c(list(wide_list = new_data_list,
                             combine_method = ascm_multi$combine_method,
                             Z = data_list$Z,
                             progfunc = ascm_multi$progfunc,
                             scm = ascm_multi$scm,
                             fixedeff = ascm_multi$fixedeff,
                             outcomes_str = ascm_multi$outcomes),
                        ascm_multi$extra_args))
  
  resids <- predict(new_ascm, att = T)[1:tpost, , drop = F]
  
  # permute residuals and compute test statistic
  if(is.null(stat_func)) {
    stat_func <- function(x) {
      x <- na.omit(x)
      (sum(abs(x) ^ q)  / sqrt(length(x))) ^ (1 / q)
    }
  }
  if(type == "iid") {
    test_stats <- sapply(1:ns, 
                         function(x) {
                           idxs <- sample(1:nrow(resids))
                           reorder <- resids[idxs, , drop = F]
                           apply(reorder[(t0 + 1):tpost, ,drop = F], 2, stat_func)
                         })
  } else {
    ## increment time by one step and wrap
    test_stats <- sapply(0:(tpost - 1),
                         function(j) {
                           reorder <- resids[(0:(tpost -1) + j) %% tpost + 1, ,drop = F]
                           if(!all(dim(reorder) == dim(resids))) {
                             stop("Error in block resampling")
                           }
                           apply(reorder[(t0 + 1):tpost, , drop = F], 2, stat_func)
                         })
  }
  
  return(list(resids = resids,
              test_stats = matrix(test_stats, nrow = k),
              stat_func = stat_func))
}


#' Compute conformal p-value
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param h0 Null hypothesis to test
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param stat_func Function to compute test statistic
#' 
#' @return Computed p-value
#' @noRd
compute_permute_pval_multiout <- function(data_list, ascm_multi, h0,
                                          post_length, type,
                                          q, ns, stat_func) {
  t0 <- ncol(data_list$X[[1]]) - post_length
  tpost <- t0 + post_length
  
  out <- compute_permute_test_stats_multiout(data_list, ascm_multi, h0,
                                             post_length, type, q, ns, stat_func)
  k <- length(data_list$X)
  
  comb_stat <- mean(apply(out$resids[(t0 + 1):tpost, , drop = F], 2, out$stat_func), na.rm = TRUE)
  comb_test_stats <- apply(out$test_stats, 2, mean, na.rm = TRUE)
  # if(h0 == 0) {
  #   hist(comb_test_stats)
  #   abline(v = comb_stat)
  #   print(mean(comb_stat <= comb_test_stats))
  #   print(1 - mean(comb_stat > comb_test_stats))
  # }
  1 - mean(comb_stat > comb_test_stats)
}

#' Compute conformal p-value
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param grid Set of null hypothesis to test for inversion
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param stat_func Function to compute test statistic
#' 
#' @return (lower bound of interval, upper bound of interval, p-value for null of 0 effect)
#' @noRd
compute_permute_ci_multiout <- function(data_list, ascm_multi, grid,
                                        post_length, alpha, type,
                                        q, ns, lin_h0 = NULL, stat_func) {
  # make sure 0 is in the grid
  if(is.null(lin_h0)) {
    grid <- lapply(grid, function(x) c(x, 0))
    k <- length(grid)
    # get all combinations of grid
    grid <- expand.grid(grid)
    grid_low <- NULL
  } else {
    k <- length(lin_h0)
    # keep track of low dimensional grid
    grid_low <- c(grid, 0)
    # transform into high dimensional grid with linear hypothesis
    grid <- sapply(lin_h0, function(x) x * grid_low)
  }
  ps <- apply(grid, 1,
              function(x) {
                compute_permute_pval_multiout(data_list, ascm_multi, x, 
                                              post_length, type, q, ns, stat_func)
              })
  sapply(1:k, 
         function(i) c(min(grid[ps >= alpha, i]), 
                       max(grid[ps >= alpha, i]), 
                       ps[apply(grid == 0, 1, all)]))
}




#' Drop unit i from data
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
#' @param i Unit to drop
#' @noRd
drop_unit_i <- function(wide_data, Z, i) {
  
  new_wide_data <- list()
  new_wide_data$trt <- wide_data$trt[-i]
  new_wide_data$X <- wide_data$X[-i,, drop = F]
  new_wide_data$y <- wide_data$y[-i,, drop = F]
  
  X0 <- new_wide_data$X[new_wide_data$trt == 0,, drop = F]
  x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1,, drop = F]),
               ncol=1)
  y0 <- new_wide_data$y[new_wide_data$trt == 0,, drop = F]
  y1 <- colMeans(new_wide_data$y[new_wide_data$trt == 1,, drop = F])
  
  new_synth_data <- list()
  new_synth_data$Z0 <- t(X0)
  new_synth_data$X0 <- t(X0)
  new_synth_data$Z1 <- x1
  new_synth_data$X1 <- x1
  new_Z <- if(!is.null(Z)) Z[-i, , drop = F] else NULL
  
  return(list(wide_data = new_wide_data,
              synth_data = new_synth_data,
              Z = new_Z))
}

#' Drop unit i from data
#' @param wide_list (X, y, trt)
#' @param Z Covariates matrix
#' @param i Unit to drop
#' @noRd
drop_unit_i_multiout <- function(wide_list, Z, i) {
  
  new_wide_data <- list()
  new_wide_data$trt <- wide_list$trt[-i]
  new_wide_data$X <- lapply(wide_list$X, function(x) x[-i,, drop = F])
  new_wide_data$y <- lapply(wide_list$y, function(x) x[-i,, drop = F])
  new_Z <- if(!is.null(Z)) Z[-i, , drop = F] else NULL
  
  return(list(wide_list = new_wide_data,
              Z = new_Z))
}


#' Estimate standard errors for single ASCM with the jackknife
#' Do this for ridge-augmented synth
#' @param ascm Fitted augsynth object
#' 
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"se"}{Standard error estimate}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
jackknife_se_single <- function(ascm) {
  
  wide_data <- ascm$data
  synth_data <- ascm$data$synth_data
  n <- nrow(wide_data$X)
  n_c <- dim(synth_data$Z0)[2]
  Z <- wide_data$Z
  
  t0 <- dim(synth_data$Z0)[1]
  tpost <- ncol(wide_data$y)
  t_final <- dim(synth_data$Y0plot)[1]
  errs <- matrix(0, n_c, t_final - t0)
  
  
  # only drop out control units with non-zero weights
  nnz_weights <- numeric(n)
  nnz_weights[wide_data$trt == 0] <- round(ascm$weights, 3) != 0
  # if more than one unit is treated, include them in the jackknife
  if(sum(wide_data$trt) > 1) {
    nnz_weights[wide_data$trt == 1] <- 1
  }
  
  trt_idxs <- (1:n)[as.logical(nnz_weights)]
  
  
  # jackknife estimates
  ests <- vapply(trt_idxs,
                 function(i) {
                   # drop unit i
                   new_data <- drop_unit_i(wide_data, Z, i)
                   # refit
                   new_ascm <- do.call(fit_augsynth_internal,
                                       c(list(wide = new_data$wide,
                                              synth_data = new_data$synth_data,
                                              Z = new_data$Z,
                                              progfunc = ascm$progfunc,
                                              scm = ascm$scm,
                                              fixedeff = ascm$fixedeff),
                                         ascm$extra_args))
                   # get ATT estimates
                   est <- predict(new_ascm, att = T)[(t0 + 1):t_final]
                   c(est, mean(est))
                 },
                 numeric(tpost + 1))
  # convert to matrix
  ests <- matrix(ests, nrow = tpost + 1, ncol = length(trt_idxs))
  ## standard errors
  se2 <- apply(ests, 1,
               function(x) (n - 1) / n * sum((x - mean(x, na.rm = T)) ^ 2))
  se <- sqrt(se2)
  
  out <- list()
  att <- predict(ascm, att = T)
  out$att <- c(att, mean(att[(t0 + 1):t_final]))
  
  out$se <- c(rep(NA, t0), se)
  # out$sigma <- NA
  return(out)
}


#' Compute standard errors using the jackknife
#' @param multisynth fitted multisynth object
#' @param relative Whether to compute effects according to relative time
#' @noRd
jackknife_se_multi <- function(multisynth, relative=NULL, alpha = 0.05, att_weight = NULL) {
  ## get info from the multisynth object
  if(is.null(relative)) {
    relative <- multisynth$relative
  }
  n_leads <- multisynth$n_leads
  n <- nrow(multisynth$data$X)
  att <- predict(multisynth, att=T, att_weight = att_weight)
  outddim <- nrow(att)
  
  J <- length(multisynth$grps)
  
  ## drop each unit and estimate overall treatment effect
  jack_est <- vapply(1:n,
                     function(i) {
                       msyn_i <- drop_unit_i_multi(multisynth, i)
                       pred <- predict(msyn_i[[1]], relative=relative, att=T, att_weight = att_weight)
                       if(length(msyn_i[[2]]) != 0) {
                         out <- matrix(NA, nrow=nrow(pred), ncol=(J+1))
                         out[,-(msyn_i[[2]]+1)] <- pred
                       } else {
                         out <- pred
                       }
                       out
                     },
                     matrix(0, nrow=outddim,ncol=(J+1)))
  
  se2 <- apply(jack_est, c(1,2),
               function(x) (n-1) / n * sum((x - mean(x,na.rm=T))^2, na.rm=T))
  lower_bound <- att - qnorm(1 - alpha / 2) * sqrt(se2)
  upper_bound <- att + qnorm(1 - alpha / 2) * sqrt(se2)
  return(list(att = att, se = sqrt(se2),
              lower_bound = lower_bound, upper_bound = upper_bound))
  
}

#' Helper function to drop unit i and refit
#' @param msyn multisynth_object
#' @param i Unit to drop
#' @noRd
drop_unit_i_multi <- function(msyn, i) {
  
  n <- nrow(msyn$data$X)
  time_cohort <- msyn$time_cohort
  which_t <- (1:n)[is.finite(msyn$data$trt)]
  
  not_miss_j <- which_t %in% setdiff(which_t, i)
  
  # drop unit i from data
  drop_i <- msyn$data
  drop_i$X <- msyn$data$X[-i, , drop = F]
  drop_i$y <- msyn$data$y[-i, , drop = F]
  drop_i$trt <- msyn$data$trt[-i]
  drop_i$mask <- msyn$data$mask[not_miss_j,, drop = F]
  
  if(!is.null(msyn$data$Z)) {
    drop_i$Z <- msyn$data$Z[-i, , drop = F]
  } else {
    drop_i$Z <- NULL
  }
  
  long_df <- msyn$long_df
  unit <- colnames(long_df)[1]
  # make alphabetical, because the ith unit is the index in alphabetical ordering
  long_df <- long_df[order(long_df[, unit, drop = TRUE]),]
  ith_unit <- unique(long_df[,unit, drop = TRUE])[i]
  long_df <- long_df[long_df[,unit, drop = TRUE] != ith_unit,]
  
  # re-fit everything
  args_list <- list(wide = drop_i, relative = msyn$relative,
                    n_leads = msyn$n_leads, n_lags = msyn$n_lags,
                    nu = msyn$nu, lambda = msyn$lambda,
                    V = msyn$V,
                    force = msyn$force, n_factors = msyn$n_factors,
                    scm = msyn$scm, time_w = msyn$time_w,
                    lambda_t = msyn$lambda_t,
                    fit_resids = msyn$fit_resids,
                    time_cohort = msyn$time_cohort, long_df = long_df,
                    how_match = msyn$how_match)
  msyn_i <- do.call(multisynth_formatted, c(args_list, msyn$extra_pars))
  
  # check for dropped treated units/time periods
  if(time_cohort) {
    dropped <- which(!msyn$grps %in% msyn_i$grps)
  } else {
    dropped <- which(!not_miss_j)
  }
  return(list(msyn_i,
              dropped))
}


#' Estimate standard errors for multi outcome ascm with jackknife
#' @param ascm Fitted augsynth object
#' @noRd
jackknife_se_multiout <- function(ascm) {
  
  wide_data <- ascm$data
  wide_list <- ascm$data_list
  n <- nrow(wide_data$X)
  Z <- wide_data$Z
  
  
  # only drop out control units with non-zero weights
  nnz_weights <- numeric(n)
  nnz_weights[wide_data$trt == 0] <- round(ascm$weights, 3) != 0
  
  trt_idxs <- (1:n)[as.logical(nnz_weights)]
  
  # jackknife estimates
  ests <- lapply(trt_idxs,
                 function(i) {
                   # drop unit i
                   new_data <- drop_unit_i_multiout(wide_list, Z, i)
                   # refit
                   new_ascm <- do.call(fit_augsynth_multiout_internal,
                                       c(list(wide = new_data$wide,
                                              combine_method = ascm$combine_method,
                                              Z = new_data$Z,
                                              progfunc = ascm$progfunc,
                                              scm = ascm$scm,
                                              fixedeff = ascm$fixedeff,
                                              outcomes_str = ascm$outcomes),
                                         ascm$extra_args))
                   new_ascm$outcomes <- ascm$outcomes
                   new_ascm$data_list <- ascm$data_list
                   new_ascm$data$time <- ascm$data$time
                   # get ATT estimates
                   est <- predict(new_ascm, att = T)
                   est <- est[as.numeric(rownames(est)) >= ascm$t_int,, drop = F]
                   rbind(est, colMeans(est, na.rm = T))
                 })
  ests <- simplify2array(ests)
  ## standard errors
  se2 <- apply(ests, c(1, 2),
               function(x) (n - 1) / n * sum((x - mean(x, na.rm = T)) ^ 2))
  se <- sqrt(se2)
  out <- list()
  att <- predict(ascm, att = T)
  att_post <- colMeans(att[as.numeric(rownames(att)) >= ascm$t_int,, drop = F],
                       na.rm = T)
  out$att <- rbind(att, att_post)
  t0 <- sum(as.numeric(rownames(att)) < ascm$t_int)
  out$se <- rbind(matrix(NA, t0, ncol(se)), se)
  out$sigma <- NA
  return(out)
}



#' Compute the weighted bootstrap distribution
#' @param multisynth fitted multisynth object
#' @param rweight Function to draw random weights as a function of n (e.g rweight(n))
#' @param relative Whether to compute effects according to relative time
#' @noRd
weighted_bootstrap_multi <- function(multisynth,
                                     rweight = rwild_b,
                                     n_boot = 1000,
                                     alpha = 0.05,
                                     att_weight = NULL,
                                     relative=NULL) {
  ## get info from the multisynth object
  if(is.null(relative)) {
    relative <- multisynth$relative
  }
  
  n <- nrow(multisynth$data$X)
  att <- predict(multisynth, att=T, att_weight = att_weight)
  outddim <- nrow(att)
  n1 <- sum(is.finite(multisynth$data$trt))
  J <- length(multisynth$grps)
  
  
  # draw random weights to get bootstrap distribution
  bs_est <- vapply(1:n_boot,
                   function(i) {
                     Z <- rweight(n)# / sqrt(n1)
                     
                     predict(multisynth, att=T, att_weight = att_weight, bs_weight = Z) - sum(Z) / n1 * att
                   },
                   matrix(0, nrow=outddim,ncol=(J+1)))
  
  se2 <- apply(bs_est, c(1,2),
               function(x) mean((x - mean(x))^2, na.rm=T))
  bias <- apply(bs_est, c(1,2),
                function(x) mean(x, na.rm=T))
  upper_bound <- att - apply(bs_est, c(1,2),
                             function(x) quantile(x, alpha / 2, na.rm = T))
  
  lower_bound <- att - apply(bs_est, c(1,2),
                             function(x) quantile(x, 1 - alpha / 2, na.rm = T))
  
  return(list(att = att,
              bias = bias,
              se = sqrt(se2),
              upper_bound = upper_bound,
              lower_bound = lower_bound))
  
}

#' Bayesian bootstrap
#' @param n Number of units
#' @export
rdirichlet_b <- function(n) {
  Z <- as.numeric(rgamma(n, 1, 1))
  return(Z / sum(Z) * n)
}

#' Non-parametric bootstrap
#' @param n Number of units
#' @export
rmultinom_b <- function(n) as.numeric(rmultinom(1, n, rep(1 / n, n)))

#' Wild bootstrap (Mammen 1993)
#' @param n Number of units
#' @export
rwild_b <- function(n) {
  sample(c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2 ), n,
         replace = TRUE,
         prob = c((sqrt(5) + 1)/ (2 * sqrt(5)), (sqrt(5) - 1) / (2 * sqrt(5))))
}



fit_ridgeaug_formatted <- function (wide_data, synth_data, Z = NULL, lambda = NULL, ridge = T, 
          scm = T, lambda_min_ratio = 1e-08, n_lambda = 20, lambda_max = NULL, 
          holdout_length = 1, min_1se = T, V = NULL, residualize = FALSE, 
          ...) {
  extra_params = list(...)
  if (length(extra_params) > 0) {
    warning("Unused parameters in using ridge augmented weights: ", 
            paste(names(extra_params), collapse = ", "))
  }
  X <- wide_data$X
  y <- wide_data$y
  trt <- wide_data$trt
  lambda_errors <- NULL
  lambda_errors_se <- NULL
  lambdas <- NULL
  X_cent <- apply(X, 2, function(x) x - mean(x[trt == 0]))
  X_c <- X_cent[trt == 0, , drop = FALSE]
  X_1 <- matrix(colMeans(X_cent[trt == 1, , drop = FALSE]), 
                nrow = 1)
  y_cent <- apply(y, 2, function(x) x - mean(x[trt == 0]))
  y_c <- y_cent[trt == 0, , drop = FALSE]
  t0 <- ncol(X_c)
  V <- make_V_matrix(t0, V)
  X_c <- X_c %*% V
  X_1 <- X_1 %*% V
  new_synth_data <- synth_data
  if (!is.null(Z)) {
    Z_cent <- apply(Z, 2, function(x) x - mean(x[trt == 0]))
    Z_c <- Z_cent[trt == 0, , drop = FALSE]
    Z_1 <- matrix(colMeans(Z_cent[trt == 1, , drop = FALSE]), 
                  nrow = 1)
    if (residualize) {
      Xc_hat <- Z_c %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c) %*% 
        X_c
      X1_hat <- Z_1 %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c) %*% 
        X_c
      res_t <- X_1 - X1_hat
      res_c <- X_c - Xc_hat
      X_c <- res_c
      X_1 <- res_t
      X_cent[trt == 0, ] <- res_c
      X_cent[trt == 1, ] <- res_t
      new_synth_data$Z1 <- t(res_t)
      new_synth_data$X1 <- t(res_t)
      new_synth_data$Z0 <- t(res_c)
      new_synth_data$X0 <- t(res_c)
    }
    else {
      sdz <- apply(Z_c, 2, sd)
      sdx <- sd(X_c)
      Z_c <- sdx * t(t(Z_c)/sdz)
      Z_1 <- sdx * Z_1/sdz
      X_c <- cbind(X_c, Z_c)
      X_1 <- cbind(X_1, Z_1)
      new_synth_data$Z1 <- t(X_1)
      new_synth_data$X1 <- t(X_1)
      new_synth_data$Z0 <- t(X_c)
      new_synth_data$X0 <- t(X_c)
      V <- diag(ncol(X_c))
    }
  }
  else {
    new_synth_data$Z1 <- t(X_1)
    new_synth_data$X1 <- t(X_1)
    new_synth_data$Z0 <- t(X_c)
    new_synth_data$X0 <- t(X_c)
  }
  out <- fit_ridgeaug_inner(X_c, X_1, trt, new_synth_data, 
                            lambda, ridge, scm, lambda_min_ratio, n_lambda, lambda_max, 
                            holdout_length, min_1se)
  weights <- out$weights
  synw <- out$synw
  lambda <- out$lambda
  lambdas <- out$lambdas
  lambda_errors <- out$lambda_errors
  lambda_errors_se <- out$lambda_errors_se
  if (!is.null(Z)) {
    if (residualize) {
      no_cov_weights <- weights
      ridge_w <- t(t(Z_1) - t(Z_c) %*% weights) %*% solve(t(Z_c) %*% 
                                                            Z_c) %*% t(Z_c)
      weights <- weights + t(ridge_w)
    }
    else {
      no_cov_weights <- NULL
    }
  }
  l2_imbalance <- sqrt(sum((synth_data$X0 %*% weights - synth_data$X1)^2))
  uni_w <- matrix(1/ncol(synth_data$X0), nrow = ncol(synth_data$X0), 
                  ncol = 1)
  unif_l2_imbalance <- sqrt(sum((synth_data$X0 %*% uni_w - 
                                   synth_data$X1)^2))
  scaled_l2_imabalance <- l2_imbalance/unif_l2_imbalance
  mhat <- matrix(0, nrow = nrow(y), ncol = ncol(y))
  ridge_mhat <- mhat
  if (!is.null(Z)) {
    if (residualize) {
      ridge_mhat <- ridge_mhat + Z_cent %*% solve(t(Z_c) %*% 
                                                    Z_c) %*% t(Z_c) %*% y_c
      yc_hat <- ridge_mhat[trt == 0, , drop = F]
      y_c <- y_c - yc_hat
    }
    else {
      X_cent <- cbind(X_cent, Z_cent)
    }
  }
  if (ridge) {
    ridge_mhat <- ridge_mhat + X_cent %*% solve(t(X_c) %*% 
                                                  X_c + lambda * diag(ncol(X_c))) %*% t(X_c) %*% y_c
  }
  output <- list(weights = weights, l2_imbalance = l2_imbalance, 
                 scaled_l2_imbalance = scaled_l2_imabalance, mhat = mhat, 
                 lambda = lambda, ridge_mhat = ridge_mhat, synw = synw, 
                 lambdas = lambdas, lambda_errors = lambda_errors, lambda_errors_se = lambda_errors_se)
  if (!is.null(Z)) {
    output$no_cov_weights <- no_cov_weights
    z_l2_imbalance <- sqrt(sum((t(Z_c) %*% weights - t(Z_1))^2))
    z_unif_l2_imbalance <- sqrt(sum((t(Z_c) %*% uni_w - t(Z_1))^2))
    z_scaled_l2_imbalance <- z_l2_imbalance/z_unif_l2_imbalance
    output$covariate_l2_imbalance <- z_l2_imbalance
    output$scaled_covariate_l2_imbalance <- z_scaled_l2_imbalance
  }
  return(output)
}

make_V_matrix <- function (t0, V) {
  if (is.null(V)) {
    V <- diag(rep(1, t0))
  }
  else if (is.vector(V)) {
    if (length(V) != t0) {
      stop(paste("`V` must be a vector with", t0, "elements or a", 
                 t0, "x", t0, "matrix"))
    }
    V <- diag(V)
  }
  else if (ncol(V) == 1 & nrow(V) == t0) {
    V <- diag(c(V))
  }
  else if (ncol(V) == t0 & nrow(V) == 1) {
    V <- diag(c(V))
  }
  else if (nrow(V) == t0) {
  }
  else {
    stop(paste("`V` must be a vector with", t0, "elements or a", 
               t0, "x", t0, "matrix"))
  }
  return(V)
}

fit_ridgeaug_inner <- function (X_c, X_1, trt, synth_data, lambda, ridge, scm, lambda_min_ratio, 
                                n_lambda, lambda_max, holdout_length, min_1se) {
  lambda_errors <- NULL
  lambda_errors_se <- NULL
  lambdas <- NULL
  if (scm) {
    syn <- fit_synth_formatted(synth_data)$weights
  }
  else {
    syn <- rep(1/sum(trt == 0), sum(trt == 0))
  }
  if (ridge) {
    if (is.null(lambda)) {
      cv_out <- cv_lambda(X_c, X_1, synth_data, trt, holdout_length, 
                          scm, lambda_max, lambda_min_ratio, n_lambda, 
                          min_1se)
      lambda <- cv_out$lambda
      lambda_errors <- cv_out$lambda_errors
      lambda_errors_se <- cv_out$lambda_errors_se
      lambdas <- cv_out$lambdas
    }
    ridge_w <- t(t(X_1) - t(X_c) %*% syn) %*% solve(t(X_c) %*% 
                                                      X_c + lambda * diag(ncol(X_c))) %*% t(X_c)
  }
  else {
    ridge_w <- matrix(0, ncol = sum(trt == 0), nrow = 1)
  }
  weights <- syn + t(ridge_w)
  return(list(weights = weights, synw = syn, lambda = lambda, 
              lambdas = lambdas, lambda_errors = lambda_errors, lambda_errors_se = lambda_errors_se))
}

fit_synth_formatted <- function (synth_data, V = NULL) {
  t0 <- dim(synth_data$Z0)[1]
  V <- make_V_matrix(t0, V)
  weights <- synth_qp(synth_data$X1, t(synth_data$X0), V)
  l2_imbalance <- sqrt(sum((synth_data$Z0 %*% weights - synth_data$Z1)^2))
  uni_w <- matrix(1/ncol(synth_data$Z0), nrow = ncol(synth_data$Z0), 
                  ncol = 1)
  unif_l2_imbalance <- sqrt(sum((synth_data$Z0 %*% uni_w - 
                                   synth_data$Z1)^2))
  scaled_l2_imbalance <- l2_imbalance/unif_l2_imbalance
  return(list(weights = weights, l2_imbalance = l2_imbalance, 
              scaled_l2_imbalance = scaled_l2_imbalance))
}

synth_qp <- function (X1, X0, V) {
  Pmat <- X0 %*% V %*% t(X0)
  qvec <- -t(X1) %*% V %*% t(X0)
  n0 <- nrow(X0)
  A <- rbind(rep(1, n0), diag(n0))
  l <- c(1, numeric(n0))
  u <- c(1, rep(1, n0))
  settings = osqp::osqpSettings(verbose = FALSE, eps_rel = 1e-08, 
                                eps_abs = 1e-08)
  sol <- osqp::solve_osqp(P = Pmat, q = qvec, A = A, l = l, 
                          u = u, pars = settings)
  return(sol$x)
}
