#' Find local maxima in a numeric vector
#'
#' A vectorized implementation to find local peaks (maxima) in a sequence of
#' numbers. A point is considered a local maximum if it is greater than or
#' equal to its left neighbor and strictly greater than its right neighbor.
#'
#' @param data A numeric vector in which to find local maxima.
#'
#' @return A tibble with two columns: `new_mu_idx` (the 1-based indices of the
#'   maxima) and `mu_prob` (the values of the data at those indices). Returns an
#'   empty tibble if no maxima are found.
#'
#' @importFrom dplyr lag lead
#' @importFrom tibble tibble
#' @noRd
find_local_maxima_vectorized <- function(data) {
  n <- length(data)
  if (n < 3) return(tibble(new_mu_idx = integer(0), mu_prob = numeric(0)))

  # data[i-1] <= data[i] & data[i] > data[i+1]
  # Using lead and lag is safer and more explicit
  max_indices <- which(lag(data, default = -Inf) <= data & lead(data, default = -Inf) < data)

  tibble(new_mu_idx = max_indices, mu_prob = data[max_indices])
}

#' Resolve overlapping intervals by adjusting sigma
#'
#' This helper function takes a set of new peaks (mu and sigma values) and
#' adjusts their sigmas to resolve any direct overlaps. If two intervals
#' (mu +/- sigma) overlap, their sigmas are proportionally resized based on
#' their original widths to eliminate the overlap.
#'
#' @param mu_sig_tibble A tibble containing new peak candidates. Must have
#'   columns `new_mu` and `new_sigma`, and should be sorted by `new_mu`.
#'
#' @return The input tibble with adjusted `new_sigma` values.
#'
#' @importFrom dplyr arrange
#' @noRd
collision_resolver <- function(mu_sig_tibble) {
  if (nrow(mu_sig_tibble) <= 1) {
    return(mu_sig_tibble)
  }

  mu_sig_tibble <- mu_sig_tibble %>% arrange(new_mu)
  mus <- mu_sig_tibble$new_mu
  sigs <- mu_sig_tibble$new_sigma

  overlap_check <- function(a, b) { (a[1] - b[2]) * (a[2] - b[1]) < 0 }

  for (i in seq_along(mus)[-length(mus)]) {
    interval1 <- c(mus[i] - sigs[i], mus[i] + sigs[i])
    interval2 <- c(mus[i+1] - sigs[i+1], mus[i+1] + sigs[i+1])

    if (overlap_check(interval1, interval2)) {
      len_ratio <- sigs[i] / (sigs[i] + sigs[i+1])
      dist <- mus[i+1] - mus[i]
      sigs[i] <- round(dist * len_ratio)
      sigs[i+1] <- round(dist * (1 - len_ratio))
    }
  }

  mu_sig_tibble$new_sigma <- sigs
  return(mu_sig_tibble)
}
