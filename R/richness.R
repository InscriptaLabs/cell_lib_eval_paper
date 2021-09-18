library(tidyverse)
library(extraDistr)
# Context: consider sampling m times (with replacement)
# from a set of S items where the probabilities of each
# item being selected are given by the vector p.
#
# This function returns the expected number of items
# that will be seen at least n times.
#
# Inputs:
#   p - Vector of probabilities - each must be in the
#       range [0,1] and they must sum to one
#   m - Positive integer specifying the number of items
#       sampled.  Can be scalar or a vector.
#   n - Positive integer specifying the minimum number of
#       observations required.  Must be a scalar.
#   f - Edit efficiency. Can be a scalar or vector.
#       Must be in the range [0,1]
richness_mean <- function(p, m, n=1, f=1) {
  g_sum <- sapply(outer(m,f), function(mf){sum(ppois(n-1,mf*p))})
  if(length(m) > 0 & length(f) > 0)
    g_sum <- matrix(g_sum, nrow=length(m), ncol=length(f))
  out <- length(p) - g_sum

  # Reshape & rename results
  if(length(m)==1) {
    out <- c(out)
    names(out) <- names(f)
  } else if (length(f)==1) {
    out <- c(out)
    names(out) <- names(m)
  } else {
    rownames(out) <- names(m)
    colnames(out) <- names(f)
  }
  return(out)
}

# Returns the variance of the number of unique observations
# among m samples, in the same context as described for the
# function uniq_obs_mean()
#
# Unlike uniq_obs_mean, n must be scalar
richness_var <- function(p, m, n=1, f=1) {
  term1 <- sapply(outer(m,f), function(mf){sum(ppois(n-1,mf*p))})
  term2 <- sapply(outer(m,f), function(mf){sum(ppois(n-1,mf*p)^2)})
  out <- term1 - term2
  out[out < .Machine$double.eps] <- 0
  if(length(m) > 0 & length(f) > 0)
    out <- matrix(out, nrow=length(m), ncol=length(f))

  # Reshape & rename results
  if(length(m)==1) {
    out <- c(out)
    names(out) <- names(f)
  } else if (length(f)==1) {
    out <- c(out)
    names(out) <- names(m)
  } else {
    rownames(out) <- names(m)
    colnames(out) <- names(f)
  }
  return(out)
}


# return the mean and variance of fractional richness when
# the edit frequencies follow a beta distribution with known cv.
# Parameters:
#   F - the effective sampling (edit fraction x cells screened / library size)
#   cv - the CV of the edit frequencies
#   lib_size - library size, needed for variance calculation
#   min_obs - the minimum required number of observations
#   type1_err - type 1 error rate for confidence intervals
# Both F and cv can be vector valued
fractional_richness_beta <- function(F, cv, lib_size, min_obs=1, type1_err=0.05, outer_product=TRUE) {
  k <- 0:(min_obs-1)
  # Pre-populate the results with the optimal value for cases where CV==0
  if (outer_product) {
    matrix_F <- matrix(F, nrow=length(F), ncol=length(cv))
    richness_mean <- 1-exp(-matrix_F)
    array_F <- outer(matrix_F, rep(1, length(k)))
    array_k <- outer(matrix(1, nrow=length(F), ncol=length(cv)), k)
    term <- array_F^array_k / factorial(array_k) * exp(-array_F)
    richness_var <- apply(term * (1-term), 1:2, sum)
  } else if (length(F)==length(cv) || length(F)==1 || length(cv)==1) {
    if(length(F) == 1 & length(cv) > 1)
      F = rep(F, length(cv))
    if(length(F) > 1 & length(cv) == 1)
      cv = rep(cv, length(F))
    richness_mean <- 1-exp(-F)
    t1 <- outer(F, k, "^")
    t2 <- matrix(factorial(k), nrow=length(F), ncol=length(k), byrow=TRUE)
    t3 <- exp(-matrix(F, nrow=length(F), ncol=length(k)))
    term <- t1 / t2 * t3
    richness_var <- apply(term * (1-term), 1, sum)
  } else {
    stop("If outer_product is FALSE then F and cv must be the same length")
  }
  # Fill in results for all cases where CV > 0
  pos <- cv > 0
  if(any(pos)) {
    c2 <- cv[pos]^2
    if (outer_product) {
      matrix_F  <- matrix(F,  nrow=length(F), ncol=length(c2))
      matrix_c2 <- matrix(c2, nrow=length(F), ncol=length(c2), byrow=TRUE)
      array_F  <- outer(matrix_F, rep(1, length(k)))
      array_c2 <- outer(matrix_c2, rep(1, length(k)))
      array_k  <- outer(matrix(1,  nrow=length(F), ncol=length(c2)), k)
    } else {
      matrix_F  <- F
      matrix_c2 <- c2
      array_F  <- outer(F, rep(1, length(k)))
      array_c2 <- outer(c2, rep(1, length(k)))
      array_k  <- outer(rep(1, length(F)), k)
    }
    t1 <- pnbinom(min_obs-1, 1/(1+matrix_F*matrix_c2), size=1/matrix_c2)
    t21 <- (1+2*array_F*array_c2)^(-1/array_c2)
    t22 <- (array_F*array_c2/(1+2*array_F*array_c2))^(2*array_k)
    t23 <- choose(1/array_c2 + array_k - 1, array_k)
    t24 <- choose(1/array_c2 + 2*array_k - 1, array_k)
    t2 <- t21 * t22 * t23 * t24
    if (outer_product) {
      richness_mean[,pos] <- 1 - t1
      richness_var[,pos] <- t1 - apply(t2, 1:2, sum)
    } else {
      richness_mean[pos] <- 1 - t1
      richness_var[pos] <- t1 - apply(t2, 1, sum)
    }
  }
  richness_var <- richness_var/lib_size

  # Use mean and variance to construct confidence intervals
  mult <- qnorm(1-type1_err/2)
  richness_lower <- pmax(0, richness_mean - mult * sqrt(richness_var))
  richness_upper <- pmin(1, richness_mean + mult * sqrt(richness_var))

  if (outer_product) {
    expanded_F  <- matrix(F,  nrow=length(F), ncol=length(cv))
    expanded_cv <- matrix(cv, nrow=length(F), ncol=length(cv), byrow=TRUE)
  } else {
    expanded_F  <- F
    expanded_cv <- cv
  }
  tibble(
    FracSample=c(expanded_F),
    CV=c(expanded_cv),
    MinObs=min_obs,
    Mean=c(richness_mean),
    Var=c(richness_var),
    Lower=c(richness_lower),
    Upper=c(richness_upper)
  )
}


# Return beta parameters given library size and cv
get_beta_param <- function(S, cv) {
  beta_a <- (S-1)/S/cv^2 - 1/S
  beta_b <- (S-1)*beta_a
  if(beta_a <= 0 | beta_b <= 0 | cv < 1e-6)
    return(c(a=NA, b=NA))
  else
    return(c(a=beta_a, b=beta_b))
}

# Negative log likelihood of a beta-binomial parameterized
# by the CV and size of a library
negloglik_beta_binom <- function(cv, S, n, k) {
  beta_a <- (S-1)/S/cv^2 - 1/S
  if(beta_a < 0)
    stop(sprintf("negative beta_a, S=%1.3f, cv=%1.3f", S, cv))
  beta_b <- (S-1)*beta_a
  #-sum(lchoose(n, k) + lbeta(k+beta_a, n-k+beta_b) - lbeta(beta_a, beta_b))
  -sum(dbbinom(k, n, alpha=beta_a, beta=beta_b, log=TRUE))
}

cv_estimate_beta_binom <- function(k, n, lib_size=NULL, alpha=0.05, auto_clean=FALSE, na.rm=FALSE) {
  default_results <- list(
    est = NA,
    sd  = NA,
    lo  = NA,
    hi  = NA,
    nll = NA
  )

  if(length(n)==1)
    n = rep(n, length(k))

  # Handle NA values
  if(any(is.na(k)) | any(is.na(n))) {
    both_good <- !is.na(k) & !is.na(n)
    if((auto_clean | na.rm) & any(both_good)) {
      k <- k[both_good]
      n <- n[both_good]
    } else {
      return(default_results)
    }
  }

  # Ensure k is integer and nonnegative
  if(!is.integer(k)) {
    if(auto_clean)
      k <- as.integer(round(k))
    else
      stop("k must be integer")
  }
  if(any(k<0))
    stop("k must be nonnegative")

  # Ensure n is integer and nonnegative
  if(!is.integer(n)) {
    if(auto_clean)
      n <- as.integer(round(n))
    else
      stop("n must be integer")
  }
  if(any(n<0))
    stop("n must be nonnegative")

  # Ensure all n values are positive
  if(!any(n>0)) {
    return(default_results)
  } else if(any(n==0)) {
    k <- k[n>0]
    n <- n[n>0]
  }

  # Ensure at least one k is positive
  if(!any(k>0))
    return(default_results)

  # Ensure k never exceeds n
  if(any(k>n))
    stop("k must not exceed n")

  if(is.null(lib_size))
    lib_size <- length(k)
  if(length(n)==1)
    n <- rep(n, lib_size)
  res <- optim(1, negloglik_beta_binom, S=lib_size, n=n, k=k, method="Brent", lower=0, upper=sqrt(lib_size-1), hessian=TRUE)
  if(res$hessian > 0)
    cv_beta_sd <- sqrt(c(1/res$hessian))
  else
    cv_beta_sd <- 0
  # Evaluate confidence bounds on the log scale, to ensure conf interval is positive
  conf_bounds <- exp(log(res$par) + c(-1, 1) * qnorm(1-alpha/2) * cv_beta_sd/res$par)
  return(list(
    est = res$par,
    sd  = cv_beta_sd,
    lo  = max(0, conf_bounds[1]),
    hi  = min(sqrt(lib_size-1), conf_bounds[2]),
    nll = res$value
  ))
}

cv <- function(z, na.rm=FALSE) {sd(z, na.rm=na.rm)/mean(z, na.rm=na.rm)}

repr_hist_and_qq_plot <- function(k, n, repr_range=NULL, var_type="", mycolor=c("dodgerblue", "tomato")) {
  if(length(n)==length(k)) {
    to_drop <- is.na(n) | is.na(k) | (n <= 0)
    k <- k[!to_drop]
    n <- n[!to_drop]
  } else if(length(n)==1) {
    if(is.na(n) | n <= 0)
      stop("n must be positive")
    k <- k[!is.na(k)]
  } else {
    stop("n must be of length 1 or the same length as k")
  }
  S <- length(k)
  repr <- S*k/n
  hist_data <- hist(log10(repr), breaks=30, plot=FALSE)
  cv_est <- cv_estimate_beta_binom(k, n, auto_clean=TRUE)
  beta_param <- get_beta_param(S, cv_est$est)
  beta_param_lo <- get_beta_param(S, cv_est$lo)
  beta_param_hi <- get_beta_param(S, cv_est$hi)
  y_vals <- ppoints(S)
  beta_quantiles <- qbeta(y_vals, shape1=beta_param["a"], shape2=beta_param["b"])
  beta_quantiles_lo <- qbeta(y_vals, shape1=beta_param_lo["a"], shape2=beta_param_lo["b"])
  beta_quantiles_hi <- qbeta(y_vals, shape1=beta_param_hi["a"], shape2=beta_param_hi["b"])
  if(length(n)==1) {
    sampling_depth <- n
  } else {
    sampling_depth <- round(mean(n))
  }

  ecdf_data <- tibble(
    Repr = sort(repr),
    Prob = max(hist_data$counts) * y_vals
  )
  beta_cdf_data <- tibble(
    Repr = S*beta_quantiles,
    Repr_Lower = S*beta_quantiles_lo,
    Repr_Upper = S*beta_quantiles_hi,
    Prob = max(hist_data$counts) * y_vals
  )
  if(is.null(repr_range))
    repr_range <- range(repr[repr > 0])

  fig_repr_hist <- tibble(Repr=repr) %>%
    filter(Repr > 0) %>%
    ggplot(aes(x=log10(Repr))) +
      geom_histogram(breaks=hist_data$breaks, alpha=0.7) +
      scale_y_continuous(sec.axis=sec_axis(
        trans = ~./max(hist_data$counts),
        labels=percent
      )) +
      geom_line(
        data=ecdf_data,
        aes(x=log10(Repr), y=Prob),
        lwd=2,
        color=mycolor[1]
      ) +
      geom_line(
        data=beta_cdf_data,
        aes(x=log10(Repr), y=Prob),
        lwd=1,
        color=mycolor[2]
      )
  if(all(!is.na(beta_cdf_data$Repr_Lower)) & all(!is.na(beta_cdf_data$Repr_Upper))) {
    fig_repr_hist <- fig_repr_hist +
      geom_ribbon(
        data=beta_cdf_data,
        aes(xmin=log10(Repr_Lower), xmax=log10(Repr_Upper), y=Prob),
        lwd=1, color=mycolor[2], fill=mycolor[2], alpha=0.25
      )
  }
  fig_repr_hist <- fig_repr_hist +
    coord_cartesian(xlim=log10(repr_range)) +
    theme(legend.position="none") +
    labs(x="Log10(Representation)", y=sprintf("%s Count", var_type))

  fig_qq <- tibble(
    Quantiles=sort(repr)/S,
    BetaQuantiles=beta_quantiles
  ) %>%
    filter(Quantiles > 0) %>%
    ggplot(aes(x=BetaQuantiles, y=Quantiles)) +
      geom_abline(slope=1, intercept=0, alpha=0.5) +
      geom_line(alpha=1) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x="Beta Frequency Quantiles", y="Empirical Quantiles") +
      theme(legend.position = "none")

  return(list(
    fig_repr_hist=fig_repr_hist,
    fig_repr_hist_xlim=log10(repr_range),
    fig_repr_hist_ylim=c(0, max(hist_data$counts)),
    fig_qq_plot=fig_qq,
    frac_obs=mean(repr>0),
    sampling_depth=sampling_depth,
    cv_sample_est=cv(repr),
    cv_bbinom_est=cv_est$est,
    cv_bbinom_sd=cv_est$sd,
    cv_bbinom_lo=cv_est$lo,
    cv_bbinom_hi=cv_est$hi,
    repr_range=repr_range
  ))
}

richness_data_frame <- function(
  p,
  n_test,
  edit_efficiency,
  min_obs=1,
  mean_only=FALSE,
  significance=0.05,
  edit_efficiency_var=0
) {
  d_m <- richness_mean(p, n_test, n=min_obs, f=edit_efficiency)
  if(mean_only) {
    d_l <- d_u <- d_m
  } else {
    d_s <- qnorm(1-significance/2) * sqrt(richness_var(p, n_test, n=min_obs, f=edit_efficiency))
    d_l <- d_m - d_s
    d_u <- d_m + d_s
    if(edit_efficiency_var > 0) {
      d_s_mod <- richness_conf_bound(p, n_test, min_obs, edit_efficiency, edit_efficiency_var, significance=significance) %>%
        lapply(function(z){qnorm(1-significance/2)*sqrt(z)})
      d_l <- d_l - d_s_mod[['lower']]
      d_u <- d_u + d_s_mod[['upper']]
    }
  }
  split_var <- rep(1:length(edit_efficiency),each=length(n_test))
  d_mean  <- c(list(n_test), split(d_m, split_var))
  d_lower <- c(list(n_test), split(d_l, split_var))
  d_upper <- c(list(n_test), split(d_u, split_var))

  f_names <- sprintf("%f",100*edit_efficiency)
  names(d_mean) <- names(d_lower) <- names(d_upper) <- c('Tests',f_names)
  d <-         data.table::melt(as.data.table(d_mean),  id.vars='Tests', variable.name="EditEff", value.name="Mean") %>%
    inner_join(data.table::melt(as.data.table(d_lower), id.vars='Tests', variable.name="EditEff", value.name="Lower"), by=c('Tests','EditEff')) %>%
    inner_join(data.table::melt(as.data.table(d_upper), id.vars='Tests', variable.name="EditEff", value.name="Upper"), by=c('Tests','EditEff'))
  
  d <- d %>% mutate(EditEff = as.numeric(as.character(EditEff)))
  if(mean_only)
    d <- d %>% select(-Lower,-Upper)
  return(d)
}

