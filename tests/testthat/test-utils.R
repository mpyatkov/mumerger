test_that("find_local_maxima_vectorized finds correct peaks", {
  # Simple case with two distinct peaks
  data <- c(1, 2, 5, 2, 3, 8, 4)
  expected <- tibble::tibble(new_mu_idx = c(3L, 6L), mu_prob = c(5, 8))
  expect_equal(find_local_maxima_vectorized(data), expected)

  # Case with a plateau: should identify the start of the plateau
  data_plateau <- c(1, 2, 5, 5 ,5, 3, 1)
  expected_plateau <- tibble::tibble(new_mu_idx = 5L, mu_prob = 5)
  expect_equal(find_local_maxima_vectorized(data_plateau), expected_plateau)
})

test_that("find_local_maxima_vectorized handles edge cases", {
  # Empty vector
  expect_equal(nrow(find_local_maxima_vectorized(numeric(0))), 0)

  # Short vectors (length < 3)
  expect_equal(nrow(find_local_maxima_vectorized(c(1))), 0)
  expect_equal(nrow(find_local_maxima_vectorized(c(1, 2))), 0)

  # No peaks (monotonically increasing)
  expect_equal(nrow(find_local_maxima_vectorized(c(1, 2, 3, 4))), 1)

  # No peaks (monotonically decreasing)
  expect_equal(nrow(find_local_maxima_vectorized(c(4, 3, 2, 1))), 1)
})


test_that("collision_resolver works when there is no overlap", {
  mu_sig_tibble <- tibble::tribble(
    ~new_mu, ~new_sigma,
    100,     10,
    200,     10
  ) # Intervals are [90, 110] and [190, 210] -> No overlap

  # Expect the output to be identical to the input
  expect_equal(collision_resolver(mu_sig_tibble), mu_sig_tibble)
})

test_that("collision_resolver correctly resolves a simple overlap", {
  # mu=100, sigma=20 -> [80, 120]
  # mu=110, sigma=20 -> [90, 130]
  # They clearly overlap.
  mu_sig_tibble <- tibble::tribble(
    ~new_mu, ~new_sigma,
    100,     20,
    110,     20
  )

  resolved <- collision_resolver(mu_sig_tibble)

  # Manually calculate expected result:
  # dist = 10, ratio = 20 / (20+20) = 0.5
  # new_sigma1 = round(10 * 0.5) = 5
  # new_sigma2 = round(10 * (1 - 0.5)) = 5
  expected <- tibble::tribble(
    ~new_mu, ~new_sigma,
    100,     5,
    110,     5
  )

  expect_equal(resolved, expected)

  # Check that the new intervals [95, 105] and [105, 115] touch but don't overlap
  expect_false((resolved$new_mu[1] + resolved$new_sigma[1]) > resolved$new_mu[2] - resolved$new_sigma[2])
})

test_that("collision_resolver handles single-row input", {
  mu_sig_tibble <- tibble::tibble(new_mu = 100, new_sigma = 10)
  expect_equal(collision_resolver(mu_sig_tibble), mu_sig_tibble)
})
