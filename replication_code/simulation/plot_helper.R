
guess_rank_hist_bins <- function(max_rank, N) {
  min(max_rank + 1, max(floor(N / 10), 5))
}
# Note that this function is memoized in .onLoad
adjust_gamma <- function(N, L, K=N, conf_level=0.95) {
  if (any(c(K, N, L) < 1)) {
    abort("Parameters 'N', 'L' and 'K' must be positive integers.")
  }
  if (conf_level >= 1 || conf_level <= 0) {
    abort("Value of 'conf_level' must be in (0,1).")
  }
  if (L==1) {
    gamma <- adjust_gamma_optimize(N, K, conf_level)
  }
  else {
    gamma <- adjust_gamma_simulate(N, L, K, conf_level)
  }
  gamma
}

# Adjust coverage parameter to find silmultaneous confidence intervals for the
# ECDF of a sample from the uniform distribution.
# N - length of samples
# K - number of equally spaced evaluation points, i.e. the right ends of the
# partition intervals.
adjust_gamma_optimize <- function(N, K, conf_level=0.95) {
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1) / K
    z1 <- c(0,z)
    z2 <- c(z,1)

    # pre-compute quantiles and use symmetry for increased efficiency.
    x2_lower <- qbinom(gamma / 2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], N)

    # Compute the total probability of trajectories inside the confidence
    # intervals. Initialize the set and corresponding probasbilities known
    # to be 0 and 1 for the starting value z1 = 0.
    x1 <- 0
    p_int <- 1
    for (i in seq_along(z1)) {
      tmp <- p_interior(
        p_int, x1 = x1, x2 = x2_lower[i]: x2_upper[i],
        z1 = z1[i], z2 = z2[i], gamma = gamma, N = N
      )
      x1 <- tmp$x1
      p_int <- tmp$p_int
    }
    abs(conf_level - sum(p_int))
  }
  optimize(target, c(0, 1 - conf_level), conf_level, N = N, K = K)$minimum
}

# Adjust coverage parameter to find silmultaneous confidence intervals for the
# ECDFs of multiple samples (chains) from the uniform distribution.
# N - length of samples (chains).
# L - number of samples (chains).
# K - number of equally spaced evaluation points, i.e. the right ends of the
# partition intervals.
# M - number of simulations used to determine the 'conf_level' middle quantile.
adjust_gamma_simulate <-function(N, L, K, conf_level=0.95, M=5000) {
  gamma <- numeric(M)
  z <- (1:(K - 1)) / K
  if (L > 1){
    n <- N * (L - 1)
    k <- floor(z * N * L)
    for (m in seq_len(M)) {
      u = u_scale(replicate(L, runif(N)))
      scaled_ecdfs <- apply(outer(u, z, "<="), c(2,3), sum)
      gamma[m] <- 2 * min(
        apply(
          scaled_ecdfs, 1, phyper, m = N, n = n, k = k
        ),
        apply(
          scaled_ecdfs - 1, 1, phyper, m = N, n = n, k = k, lower.tail = FALSE
        )
      )
    }

  }
  else {
    for (m in seq_len(M)) {
      u <- runif(N)
      scaled_ecdf <- colSums(outer(u, z, "<="))
      gamma[m] <- 2 * min(
        pbinom(scaled_ecdf, N, z),
        pbinom(scaled_ecdf - 1, N, z, lower.tail = FALSE)
      )
    }
  }
  alpha_quantile(gamma, 1 - conf_level)
}

p_interior <- function(p_int, x1, x2, z1, z2, gamma, N) {
  z_tilde <- (z2 - z1) / (1 - z1)

  N_tilde <- rep(N - x1, each = length(x2))
  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)

  list(p_int = rowSums(p_x2_int), x1 = x2)
}

# 100 * `alpha` percent of the trials are allowed to be rejected.
# In case of ties, return the largest value dominating at most
# 100 * (alpha + tol) percent of the values.
alpha_quantile <- function(gamma, alpha, tol = 0.001) {
  a <- unname(quantile(gamma, probs = alpha))
  a_tol <- unname(quantile(gamma, probs = alpha + tol))
  if (a == a_tol) {
    if (min(gamma) < a) {
      # take the largest value that doesn't exceed the tolerance.
      a <- max(gamma[gamma < a])
    }
  }
  a
}

# Compute simultaneous confidence intervals for one or more samples from the
# standard uniform distribution.
# N - sample length
# L - number of samples
# K - size of uniform partition defining the ECDF evaluation points.
# gamma - coverage parameter for the marginal distribution (binomial for
# one sample and hypergeometric for multiple rank transformed chains).
ecdf_intervals <- function(N, L, K, gamma) {
  lims <- list()
  z <- seq(0,1, length.out = K + 1)
  if (L == 1) {
    lims$lower <- qbinom(gamma / 2, N, z)
    lims$upper <- qbinom(1 - gamma / 2, N, z)
  } else {
    n <- N * (L - 1)
    k <- floor(z * L * N)
    lims$lower <- qhyper(gamma / 2, N, n, k)
    lims$upper <- qhyper(1 - gamma / 2, N, n, k)
  }
  lims$lower <- c(rep(lims$lower[1:K], each=2), lims$lower[K + 1])
  lims$upper <- c(rep(lims$upper[1:K], each=2), lims$upper[K + 1])
  lims
}

# Transform observations in 'x' into their corresponding fractional ranks.
u_scale <- function(x) {
  array(rank(x) / length(x), dim = dim(x), dimnames = dimnames(x))
}

# for each value in 'y', compute the fractional ranks (empirical pit values)
# with respect to 'yrep'.
empirical_pit <- function(y, yrep) {
  (1 +  apply(outer(yrep, y, "<="), 3, sum)) / (1 +length(yrep))
}


ranks_to_empirical_pit <- function(ranks, max_rank){
  (1 + ranks) / (1 + max_rank)
}

#' Compute observed coverage of posterior credible intervals.
#'
#' Uses ranks to compute coverage and surrounding uncertainty of posterior credible intervals.
#' The uncertainty is only approximate (treating coverage for each interval as a set of independent
#' Bernoulli trials, while in fact they are not independent), so for making claims on presence/
#' absence of detectable discrepancies we strongly recommend using [plot_ecdf()] or [plot_ecdf_diff()].
#' The uncertainty about the coverage can however be useful for guiding decisions on whether
#' more SBC steps should be performed (i.e. whether we can rule out that the coverage of
#' the given backend differs too much for our purposes from the optimal value).
#'
#' Note that while coverage of central posterior intervals (with the default `type = "central"`)
#' is often of the biggest practical interest, perfect calibration of central intervals
#' still leaves space for substantial problems with the model (e.g. if the posterior 25% - 50% intervals
#' contain 50% of the true values and the posterior 50% - 75% interval never contains the true value,
#' the central 50% interval still has the ideal 50% coverage), so investigating central
#' intervals should always be accompanied by checks with [plot_ecdf()] or [plot_ecdf_diff()]
#' or by using `type = "leftmost"`, because if all leftmost credible intervals are well calibrated,
#' then all intervals are well calibrated.
#'
#' @param stats a data.frame of rank statistics (e.g. as returned in the `$stats` component of [SBC_results]),
#'   at minimum should have at least `variable`, `rank` and `max_rank` columns)
#' @param width a vector of values between 0 and 1 representing widths of credible intervals for
#'   which we compute coverage.
#' @param prob determines width of the uncertainty interval around the observed coverage
#' @param inteval_type `"central"` to show coverage of central credible intervals
#'   or `"leftmost"` to show coverage of leftmost credible intervals (i.e. the observed CDF).
#' @return A `data.frame` with columns `variable`, `width` (width of the interval as given
#'   in the `width` parameter), `width_represented` the closest width that can be represented by
#'   the ranks in the input (any discrepancy needs to be judged against this rather than `width`),
#'   `estimate` - observed coverage for the interval, `ci_low`, `ci_high` the uncertainty
#'   interval around `estimate` (width of the interval is given by the `prob` argument).
#' @seealso [plot_coverage()]
#' @export
empirical_coverage <- function(stats, width, prob = 0.95, interval_type = "central") {
  stopifnot(is.data.frame(stats))
  # Ensuring backwards compatibility
  if("parameter" %in% names(stats)) {
    if(!("variable" %in% names(stats))) {
      warning("The stats parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      stats$variable <- stats$parameter
    }
  }

  if(!all(c("variable", "rank", "max_rank") %in% names(stats))) {
    stop(SBC_error("SBC_invalid_argument_error",
                   "The stats data.frame needs a 'variable', 'rank' and 'max_rank' columns"))
  }

  stopifnot(is.numeric(width))
  stopifnot(all(width >= 0) && all(width <= 1))

  stopifnot(interval_type %in% c("central", "leftmost"))

  get_low_rank <- function(max_rank, n_ranks_covered) {
    if(interval_type == "central") {
      round(max_rank / 2 - n_ranks_covered / 2)
    } else if(interval_type == "leftmost") {
      rep(0, max(length(n_ranks_covered), length(max_rank)))
    } else {
      stop("Invalid interval_type")
    }
  }

  # Some juggling to reduce memory footprint
  stats_trimmed <- dplyr::select(stats, variable, rank, max_rank)

  variable_was_character <- is.character(stats_trimmed$variable)
  if(variable_was_character) {
    stats_trimmed <- dplyr::mutate(stats_trimmed, variable = factor(variable))
  }


  long <- dplyr::full_join(stats_trimmed, data.frame(width = width), by = character())
  long <- dplyr::mutate(long,
                        n_ranks_covered = round((max_rank + 1) * width),
                        low_rank = get_low_rank(max_rank, n_ranks_covered),
                        high_rank = low_rank + n_ranks_covered - 1,
                        width_represented =  (high_rank - low_rank + 1) / (max_rank + 1),
                        is_covered = rank >= low_rank & rank <= high_rank)

  summ <- dplyr::summarise(
    dplyr::group_by(long, variable, width),
    post_alpha = sum(is_covered) + 1,
    post_beta = dplyr::n() - sum(is_covered) + 1,
    width_represented = unique(width_represented),
    # Special handling if width_represented is either 0 or 1 as in such case,
    # the result can never be different from 0 / 1 and so the CI should collapse to a point
    representable = width_represented > 0 & width_represented < 1,
    ci_low =  dplyr::if_else(representable,
                             qbeta(0.5 - prob / 2, post_alpha, post_beta),
                             width_represented),
    estimate = mean(is_covered),
    ci_high = dplyr::if_else(representable,
                             qbeta(0.5 + prob / 2, post_alpha, post_beta),
                             width_represented),
    .groups = "drop"
  )

  if(variable_was_character) {
    summ <- dplyr::mutate(summ, variable = as.character(variable))
  }


  dplyr::select(summ, -post_alpha, -post_beta, -representable)
}

#' @export
#' @import ggplot2
plot_rank_hist <- function(x, variables = NULL, bins = NULL, prob = 0.95, max_rank = x$max_rank, parameters = NULL) {
  # Ensuring backwards compatibility
  if("parameter" %in% names(x)) {
    if(!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }


  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }



  if(!all(c("variable", "rank") %in% names(x))) {
    stop("The data.frame needs a 'variable' and 'rank' columns")
  }
  n_sims <- dplyr::summarise(dplyr::group_by(x, variable), count = dplyr::n())$count
  if(length(unique(n_sims)) > 1) {
    stop("Differing number of SBC steps per variable not supported.")
  }

  if(is.null(max_rank)) {
    stop("max_rank either has to be supplied explicitly or be a column in the data")
  }
  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across variables is not supported yet.")
  }

  n_sims <- unique(n_sims)

  if(is.null(bins)){
    bins <- guess_rank_hist_bins(max_rank, n_sims)
  } else if(bins > max_rank + 1) {
    stop("Cannot use more bins than max_rank + 1")
  }

  if(!is.null(variables)) {
    x <- dplyr::filter(x, variable %in% variables)
  }

  if(nrow(x) == 0) {
    stop("No data for the selected variables.")
  }

  #CI - taken from https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R


  # Bins can differ by size (at most by 1). Build a CI that is conservative,
  # i.e. includes lower quantile of smalelr bins and higher quantile of larger bins
  larger_bin_size <- ceiling(((max_rank + 1) / bins))
  smaller_bin_size <- floor(((max_rank + 1) / bins))
  ci_lower = qbinom(0.5 * (1 - prob), size=n_sims,prob  =  smaller_bin_size / max_rank)
  ci_mean = qbinom(0.5, size=n_sims,prob  =  1 / bins)
  ci_upper = qbinom(0.5 * (1 + prob), size=n_sims,prob  =  larger_bin_size / max_rank)

  CI_polygon_x <- c(-0.1*max_rank,0,-0.1*max_rank,1.1 * max_rank,max_rank,1.1 * max_rank,-0.1 * max_rank)
  CI_polygon_y <- c(ci_lower,ci_mean,ci_upper,ci_upper,ci_mean,ci_lower,ci_lower)

  #The visualisation style taken as well from   https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
  ggplot(x, aes(x = rank)) +
        annotate("segment", x = 0, xend = max_rank, y = ci_mean, yend = ci_mean, colour = "grey25") +
          geom_polygon(data=data.frame(x= CI_polygon_x,y= CI_polygon_y),aes(x=x,y=y),fill="skyblue",color="skyblue1",alpha=0.33) +
          geom_histogram(breaks =  seq(0, max_rank, length.out = bins + 1), closed = "left" ,fill="#808080",colour="black") +
          labs(y = "count") +
          facet_wrap(~variable, scales = "free_y")

}

#' @export
data_for_ecdf_plots <- function(x,
                                       max_rank,
                                       variables = NULL,
                                       prob = 0.95,
                                       gamma = NULL,
                                       K = NULL,
                                       size = 1,
                                       alpha = 0.33,
                                       combine_variables = NULL,
                                       ecdf_alpha = NULL,
                                       ...,
                                       parameters = NULL) {

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  ranks_matrix <- x
#   if(any(!is.finite(ranks_matrix))) {
#     stop("Ranks may only contain finite values")
#   }

#   if(!is.null(variables)) {
#     ranks_matrix <- ranks_matrix[, variables]
#   }

  N <- nrow(ranks_matrix)
  if (is.null(K)) {
    if(N < 50) {
      K <- max_rank + 1
    } else {
      K <- min(max_rank + 1, N)
    }
  } else if (K == "max") {
    K <- max_rank + 1
  } else if (K == "min") {
    K <- min(max_rank + 1, N, 100)
  } else if (!is.numeric(K) & length(K) == 1) {
    stop("K must be either a single number, \"max\", \"min\" or NULL")
  } else {
    K <- as.integer(K)
  }
  if (is.null(gamma)) {
    gamma <- adjust_gamma(
      N = N,
      L = 1,
      K = K,
      conf_level = prob
    )
  }
  z <- seq(0,1, length.out = K + 1)
  z_twice <- c(0, rep(z[2:(K + 1)], each = 2))

  limits_df <- as.data.frame(ecdf_intervals(
    N = N,
    L = 1,
    K = K,
    gamma = gamma))
  limits_df <- dplyr::mutate(limits_df,
                             x = z_twice,
                             lower = lower / N,
                             upper = upper / N,
                             # The uniform_val needs to be shifted w.r.t z_twice
                             uniform_val =  c(rep(z[1:K], each = 2), 1))

  # Combining pit and ecdf calculations in one function to avoid
  # numerical problems causing issue #79
  base_vals <- floor((0:K) * ((max_rank + 1) / K))
  ecdf_vals <- matrix(nrow = K + 1, ncol = ncol(ranks_matrix))
  colnames(ecdf_vals) <- colnames(ranks_matrix)
  for(i in 1:(K + 1)) {
    # Note: for pit calculations we would use (col + 1) / (max_rank + 1)
    # For ecdf we would use pit <= base_val / (max_rank + 1)
    # So the "+ 1" and "<=" can be subsumed in "<"
    ecdf_vals[i,] <- colMeans(ranks_matrix < base_vals[i])
  }


  ecdf_df <- as.data.frame(ecdf_vals)
  ecdf_df$..z <- z
  ecdf_df <- tidyr::pivot_longer(ecdf_df, -..z, names_to = "variable", values_to = "ecdf")
  ecdf_df <- dplyr::rename(ecdf_df, z = ..z)
  # Allow user-specified grouping of variables + alpha on ecdf line (issue #88)
  if(is.null(ecdf_alpha)) {
    ecdf_alpha <- \(x) sqrt(1/x)
  } else if(is.numeric(ecdf_alpha) & length(ecdf_alpha) == 1) {
    ecdf_alpha_numeric <- ecdf_alpha
    ecdf_alpha <- \(x) ecdf_alpha_numeric
  } else if(!is.function(ecdf_alpha) | nargs(ecdf_alpha) != 1) {
    stop("`ecdf_alpha` must be a function taking a single argument or a single numerical value")
  }
  if (!is.null(combine_variables)) {
    if(is.function(combine_variables)) {
      combine_variables <- combine_variables(unique(ecdf_df$variable))
    }
    if(!is.list(combine_variables) | is.null(names(combine_variables))) {
      stop("`combine_variables` must be a named list or a function returning a named list")
    }

    if(!identical(unique(table(unlist(combine_variables))), 1L)) {
      stop("Duplicated variable names are not allowed in `combine_variables`")
    }
    if(!all(unlist(combine_variables) %in% ecdf_df$variable)) {
      stop("The following variables in `combine_variables` couldn't be found: ",
        paste(unlist(combine_variables)[!unlist(combine_variables) %in% ecdf_df$variable], collapse = ", "))
    }
    display_names <- names(combine_variables)
    for (i in seq_along(combine_variables)) {
      ecdf_df[ecdf_df$variable %in% combine_variables[[i]], "group"] <- display_names[i]
    }
    ecdf_df$group <- factor(ecdf_df$group, levels = display_names, ordered = TRUE)
    ecdf_df <- dplyr::mutate(ecdf_df,
      alpha = ecdf_alpha(length(unique(variable))), .by = group)
  } else {
    ecdf_df$alpha <- ecdf_alpha(1)
    ecdf_df$group <- ecdf_df$variable
  }

  structure(list(limits_df = limits_df, ecdf_df = ecdf_df, K = K, N = N, z = z_twice),
            class = "SBC_ecdf_data")
}

plot_ecdf <- function(x,
                      variables = NULL,
                      K = NULL,
                      gamma = NULL,
                      prob = 0.95,
                      size = 1,
                      alpha = 0.33,
                      combine_variables = NULL,
                      ecdf_alpha = NULL,
                      ...,
                      parameters = NULL) {

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  ecdf_data <-
    data_for_ecdf_plots(x, variables = variables,
                        prob = prob, K = K, gamma = gamma,
                        combine_variables = combine_variables, ecdf_alpha = ecdf_alpha, ...)

  N <- ecdf_data$N
  K <- ecdf_data$K
  z <- ecdf_data$z

  ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, type = "sample ECDF")
  ecdf_df$group <- factor(ecdf_df$group, levels = unique(ecdf_df$group), ordered = TRUE)
  limits_df <- ecdf_data$limits_df
  limits_df$type <- "theoretical CDF"

  # construct figure
  ggplot(ecdf_df, aes(color = type, fill = type)) +
    geom_ribbon(
      data = limits_df,
      aes(x = x, ymax = upper, ymin = lower),
      alpha = alpha,
      size = size) +
    geom_step(
      aes(x = z, y = ecdf, group = variable, alpha = alpha)
    ) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_color_manual(
      name = "",
      values = rlang::set_names(
        c("skyblue1", "black"),
        c("theoretical CDF", "sample ECDF")),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    scale_fill_manual(
      name = "",
      values = c("theoretical CDF" = "skyblue",
                 "sample ECDF" = "transparent"),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    scale_alpha_identity() +
    xlab(NULL) +
    ylab(NULL) +
    facet_wrap(~ group)
}