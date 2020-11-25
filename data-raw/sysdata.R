# Simulate R/sysdata.rda

# Arguments
length_out <- 5000
time_step <- 1
burn_steps <- 5000
x0 <- 0.78
y0 <- 0.16
z0 <- 10.0
a1 <- 5.0
a2 <- 0.1
b1 <- 3.0
b2 <- 2.0
d1 <- 0.4
d2 <- 0.01

# Simulate
sysdata <- pbsSim(
	length_out = length_out,
	time_step = time_step,
	burn_steps = burn_steps,
	x0 = x0,
	y0 = y0,
	z0 = z0,
	a1 = a1,
	a2 = a2,
	b1 = b1,
	b2 = b2,
	d1 = d1,
	d2 = d2)

# Write to R/sysdata.rda
usethis::use_data(sysdata, internal = TRUE, overwrite = TRUE)
