# Script for Marianthi

# Objective: telling whether some CpG islands that vary thoughout the lifetime
# of the bird are under- or over-represented in some functional regions of the
# genome compared to what would be expected by chance (so differently
# distributed than when all CpG islands are taken into consideration).

rm(list = ls())

library(tidyverse)

# Make a dummy table
data <- tibble(
  location = c("Promoter", "Exon", "Intron", "Intergenic"),
  obs = c(100, 95, 300, 700),
  exp = c(4, 2, 6, 6) * 10^6
)

# Compute proportions
data <- data %>%
  mutate(
    prop_obs = obs / sum(obs),
    prop_exp = exp / sum(exp)
  )

# Perform a one-way chi square test
with(data, chisq.test(x = obs, p = prop_exp))

# Permutation version
with(data, chisq.test(x = obs, p = prop_exp, simulate.p.value = TRUE))

# If the test is significant, we now want to perform a post-hoc test to know
# in which categories CpG islands might be under or over-represented.

# What is the total number of trials for the binomial test?
N <- sum(data$obs)

# For each category...
with(data, map2(obs, prop_exp, function(obs, prop_exp) {

  # Perform a binomial test: is the no. of CpG in that category different than expected?
  binom.test(obs, N, p = prop_exp, alternative = "two.sided")

  # Note: the test is two-sided to identify both under and over-representation

}))

# We can do the same and attach the results to our original data
data <- data %>%
  mutate(
    p_value = map2_dbl(
      obs, prop_exp,
      ~ binom.test(.x, N, .y, alternative = "two.sided")$p.value
    )
  )

# Note: What to return in a binomial test?
# (1) the number of successes ("special" CpG falling in that category)
# (2) the number of trials (all "special" CpG)
# (3) the probability of success (computed out of all CpG)
# (4) the P-value

# Correction for multiple testing
data <- data %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni"))

# Check out the results
data

# Note: the Bonferroni correction is quite strict. For many P-values, it may
# substantially lower statisical power to detect anything. For many categories,
# maybe opt for the Benjamini-Hochbeg method ("BH" in p.adjust). In general,
# such corrections are expected in post-hoc tests, because we explore many
# cases one by one, even if in the present case with only four groups, adjusting
# the P-values may not matter much. It would matter, though, with many tests
# being performed, e.g. with chromosomes instead of functional regions.
