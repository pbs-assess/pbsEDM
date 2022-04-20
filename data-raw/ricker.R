# Simulate
ricker <- larkin::sim(
  alpha = 2,
  beta = c(8, 0, 0, 0),
  init = rep(0.1, 8),
  p_bar = c(0.003, 0.917, 0.08),
  sigma = 0.1,
  burn = 100,
  span = 100,
  harvest = 0.2,
  extirp = 1e-06,
  seed = 123
)

# Plot spawner abundance
plot(ricker$spawners, type = "l")

# Write to data/
usethis::use_data(ricker, overwrite = TRUE)
