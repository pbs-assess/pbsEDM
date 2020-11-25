#' Simulate Three Species Dynamics via Hastings & Powell Model
#'
#' @param length_out [integer()] Time series length
#' @param time_step [numeric()] Sample frequency relative to t = R0 * T
#' @param burn_steps [integer()] Number of burn in steps
#' @param x0 [numeric()] Predator state variable.
#' @param y0 [numeric()] Prey state variable.
#' @param z0 [numeric()] Producer state variable.
#' @param a1 [numeric()] Parameter a1
#' @param a2 [numeric()] Parameter a2
#' @param b1 [numeric()] Parameter b1
#' @param b2 [numeric()] Parameter b2
#' @param d1 [numeric()] Parameter d1 
#' @param d2 [numeric()] Parameter d2
#'
#' @return [matrix()] Time series as columns.
#' @export
#'
#' @examples
#' length_out <- 5000
#' time_step <- 1
#' burn_steps <- 5000
#' x0 <- 0.78
#' y0 <- 0.16
#' z0 <- 10.0
#' a1 <- 5.0
#' a2 <- 0.1
#' b1 <- 3.0
#' b2 <- 2.0
#' d1 <- 0.4
#' d2 <- 0.01
#' 
#' s1 <- pbsSim(
#'   length_out = length_out,
#'   time_step = time_step,
#'   burn_steps = burn_steps,
#'   x0 = x0,
#'   y0 = y0,
#'   z0 = z0,
#'   a1 = a1,
#'   a2 = a2,
#'   b1 = b1,
#'   b2 = b2,
#'   d1 = d1,
#'   d2 = d2)
#' 
pbsSim <- function (length_out,
										time_step = 1,
										burn_steps = 10,
	                  x0 = 0.8,
										y0 = 0.2,
										z0 = 9.0,
										a1 = 5.0,
										a2 = 0.1,
										b1 = 2.0,
										b2 = 2.0,
										d1 = 0.4,
										d2 = 0.01) {
	
	# Initialize
	dt <- 0.001
	decimal <- as.integer(-log10(dt))
	length_sim <- as.integer((length_out + burn_steps - 1) * 1/dt)
	time_step <- as.integer(time_step)
	burn_steps <- as.integer(burn_steps)
	# Current value
	x_now <- x0
	y_now <- y0
	z_now <- z0
	# Sample matrix
	m <- matrix(0, nrow = length_out, ncol = 4)
	colnames(m) <- c("time", "predators", "prey", "producers")
	# Count
	row_index <- 1L
	sum_time <- 0
	mod_time <- 0
	# Simulate
	for (i in seq_len(length_sim + 1L)) {
		if (sum_time >= burn_steps) {
			if (mod_time %% time_step == 0) {
				m[row_index, 1] <- sum_time - burn_steps + 1
				m[row_index, 2] <- round(as.double(z_now), decimal)
				m[row_index, 3] <- round(as.double(y_now), decimal)
				m[row_index, 4] <- round(as.double(x_now), decimal)
				row_index <- row_index + 1L
				# cat("time = ", sum_time, "\n")
			}
		}
		# cat(sum_time, mod_time, x_now, y_now, z_now, "\n")
		# Compute rate
		dx_dt <- x_now * (1 - x_now) - y_now * f(x_now, a1, b1)
		dy_dt <- y_now * f(x_now, a1, b1) - z_now * f(y_now, a2, b2) - y_now * d1
		dz_dt <- z_now * f(y_now, a2, b2) - z_now * d2
		# Compute change
		dx <- dx_dt * dt
		dy <- dy_dt * dt
		dz <- dz_dt * dt
		# Update state variables
		x_now <- x_now + dx
		y_now <- y_now + dy
		z_now <- z_now + dz
		# Counts
		sum_time <- round(sum_time + dt, decimal)
		mod_time  <- round(mod_time + dt, decimal) %% 1L
	}
	return(structure(m, class = "pbsSim"))
}

#' Saturating Functional Response
#'
#' @param u [numeric()] State variable value.
#' @param a [numeric()] Parameter a_i.
#' @param b [numeric()] Parameter b_i.
#'
#' @return [numeric()]
#'
f <- function (u, a, b) {
	a * u / (1 + b * u)
}

