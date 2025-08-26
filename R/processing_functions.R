#' Process a single union of genomic regions to find consensus peaks
#'
#' This is a core processing function that takes a dataframe representing all
#' peaks that fall within a single merged "union" region. It calculates a new
#' set of consensus peaks based on a combined probabilistic landscape derived
#' from the original peaks.
#'
#' @param region_data A tibble containing peak data for a single union region.
#'   Must include columns: `start_union`, `end_union`, `start`, `end`,
#'   and `seqnames`.
#'
#' @return A tibble with the new consensus peaks in BED-like format (seqnames,
#'   start, end), or `NULL` if no peaks are found.
#'
#' @importFrom dplyr filter mutate select group_by summarise arrange desc cross_join
#' @importFrom purrr map reduce
#' @importFrom magrittr %>%
#' @export
process_one_union_region <- function(region_data) {
  # --- Step 1: Get region info and define x-axis ---
  # Using absolute genomic coordinates for the interval
  region_start <- region_data$start_union[1]
  region_end <- region_data$end_union[1]
  union_interval <- seq(from = region_start, to = region_end - 1)

  # Ensure original mu and sigma columns exist. 'mu' is the center of the original peak.
  region_data <- region_data %>%
    mutate(
      mu = round((start + end) / 2),
      sigma = (end - start) / 2
    )

  # --- Step 2: Calculate combined probability landscape ---
  #combined_prob <- process_one_condition(union_interval = union_interval, condition_df = region_data)
  combined_prob <- rowSums(mapply(dnorm,
                          MoreArgs = list(x = union_interval),
                          mean = region_data$mu,
                          sd = region_data$sigma))

  # --- Step 3: Find and rank new peak locations (maxima) ---
  # find_local_maxima_vectorized returns indices relative to the start of the vector
  local_max <- find_local_maxima_vectorized(combined_prob)

  if (nrow(local_max) == 0) return(NULL) # No peaks found, exit early

  # Python-compatible logic for number of peaks to keep
  num_peaks_to_keep <- nrow(region_data)

  # Rank maxima and convert their indices to absolute genomic coordinates
  ranked_maxima <- local_max %>%
    arrange(desc(mu_prob)) %>%
    head(n = num_peaks_to_keep) %>%
    # Convert index (e.g., 10th position) to coordinate (e.g., region_start + 9)
    mutate(new_mu = region_start + new_mu_idx - 1) %>%
    select(new_mu, mu_prob)

  if (nrow(ranked_maxima) == 0) return(NULL)

  new_sigma_df <- ranked_maxima %>%
    cross_join(select(region_data, mu, sigma)) %>%
    mutate(dists = abs(mu - new_mu)) %>%
    group_by(new_mu, mu_prob) %>%
    summarise(
      new_sigma = round(sum(sigma / (dists + 1)) / sum(1 / (dists + 1))),
      .groups = 'drop'
    ) %>%
    arrange(new_mu)

  # --- Step 5: Resolve collisions ---
  final_peaks <- collision_resolver(new_sigma_df)

  # --- Step 6: Format to final BED file structure ---
  final_peaks %>%
    mutate(
      seqnames = region_data$seqnames[1],
      # The new start/end are based on the final, resolved sigmas
      start = round(new_mu - new_sigma),
      end = round(new_mu + new_sigma)
    ) %>%
    select(seqnames, start, end)
}
