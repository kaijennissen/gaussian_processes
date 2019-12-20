library(forecast)
library(tidyverse)

# plot 4 - gaussian process -----------------------------------------------

library(MASS)
library(plyr)
library(reshape2)
library(ggplot2)


set.seed(12345)

# Calculates the covariance matrix sigma using a
# simplified version of the squared exponential function.
#
# Although the nested loops are ugly, I've checked and it's about
# 30% faster than a solution using expand.grid() and apply()
#
# Parameters:
#	X1, X2 = vectors
# 	l = the scale length parameter
# Returns:
# 	a covariance matrix
calcSigma <- function(X1, X2, l = 1) {
	Sigma <- matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))
	for (i in 1:nrow(Sigma)) {
		for (j in 1:ncol(Sigma)) {
			Sigma[i, j] <- exp(-0.5 * (abs(X1[i] - X2[j]) / l) ^ 2)
		}
	}
	return(Sigma)
}

# 1. Plot some sample functions from the Gaussian process
# as shown in Figure 2.2(a)

# Define the points at which we want to define the functions
x.star <- seq(-6, 6, len = 200)

# Calculate the covariance matrix
sigma <- calcSigma(x.star, x.star)

# Generate a number of functions from the process
n.samples <- 3
values <-
	matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
	# Each column represents a sample from a multivariate normal distribution
	# with zero mean and covariance sigma
	values[, i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id = "x")

# Plot the result
fig2a <- ggplot(values, aes(x = x, y = value)) +
	geom_rect(
		xmin = -Inf,
		xmax = Inf,
		ymin = -2,
		ymax = 2,
		fill = "grey80"
	) +
	geom_line(aes(group = variable)) +
	theme_bw() +
	scale_y_continuous(lim = c(-2.5, 2.5), name = "output, f(x)") +
	xlab("input, x")
fig2a

# 2. Now let's assume that we have some known data points;
# this is the case of Figure 2.2(b). In the book, the notation 'f'
# is used for f$y below.  I've done this to make the ggplot code
# easier later on.
f <- data.frame(x = c(-5, -3, -1, 0, 3),
				y = c(-2, 0, .5, 2, -1))

# Calculate the covariance matrices
# using the same x.star values as above
x <- f$x
k.xx <- calcSigma(x, x)
k.xxs <- calcSigma(x, x.star)
k.xsx <- calcSigma(x.star, x)
k.xsxs <- calcSigma(x.star, x.star)

# These matrix calculations correspond to equation (2.19)
# in the book.
f.star.bar <- k.xsx %*% solve(k.xx) %*% f$y
cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx) %*% k.xxs

# This time we'll plot more samples.  We could of course
# simply plot a +/- 2 standard deviation confidence interval
# as in the book but I want to show the samples explicitly here.
n.samples <- 25
values <-
	matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
	values[, i] <- mvrnorm(1, f.star.bar, cov.f.star)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id = "x")

# Plot the results including the mean function
# and constraining data points
x_star <-
	tibble(x_star = x.star, f_star = f.star.bar[, , drop = TRUE])
fig2b <- ggplot(values, aes(x = x, y = value)) +
	geom_line(aes(group = variable), colour = "#013848", alpha = .25) +
	geom_line(
		data = x_star,
		aes(x = x_star, y = f_star),
		colour = "#013848",
		size = 1
	) +
	geom_point(data = f, aes(x = x, y = y)) +
	theme_minimal() +
	scale_y_continuous(lim = c(-4, 4), name = "output, f(x)") +
	xlab("input, x")
fig2b

# 3. Now assume that each of the observed data points have some
# normally-distributed noise.

# The standard deviation of the noise
sigma.n <- 0.1

# Recalculate the mean and covariance functions
f.bar.star <-
	k.xsx %*% solve(k.xx + sigma.n ^ 2 * diag(1, ncol(k.xx))) %*% f$y
cov.f.star <-
	k.xsxs - k.xsx %*% solve(k.xx + sigma.n ^ 2 * diag(1, ncol(k.xx))) %*% k.xxs

# Recalulate the sample functions
values <-
	matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
	values[, i] <- mvrnorm(1, f.bar.star, cov.f.star)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id = "x")

# Plot the result, including error bars on the observed points
gg <- ggplot(values, aes(x = x, y = value)) +
	geom_line(aes(group = variable), colour = "grey80") +
	geom_line(
		data = x_star,
		aes(x = x_star, y = f_star),
		colour = "red",
		size = 1
	) +
	geom_errorbar(
		data = f,
		aes(
			x = x,
			y = NULL,
			ymin = y - 2 * sigma.n,
			ymax = y + 2 * sigma.n
		),
		width = 0.2
	) +
	geom_point(data = f, aes(x = x, y = y)) +
	theme_minimal() +
	scale_y_continuous(lim = c(-4, 4), name = "output, f(x)") +
	xlab("input, x")
gg



# plot 5 -  gp for weekly data --------------------------------------------
SS <- periodic_kernel(c(1:5), c(1:5))
SS <- matrix(0, nrow = 5, ncol = 5)
SS[lower.tri(SS)] <- 1:10
#1)
SS <- SS + t(SS)
isSymmetric(SS)
# 2)
SS[upper.tri(SS)] <- t(SS)[upper.tri(SS)]
isSymmetric(SS)
library(RcppZiggurat)
library(xts)
# rm(list = ls())
rw_log <- read_csv2("rw.csv") %>%
	dplyr::select(date, views) %>%
	dplyr::mutate(views = log(views)) %>%
	dplyr::filter(views > 0) %>%
	timetk::tk_xts()

rbf_kernel <- function(X1, X2, l = 1, sig = 1) {
	# RBF-Kernel
	#
	# Parameters:
	#	X1, X2 = vectors
	# 	l = the scale length parameter
	# Returns:
	# 	a covariance matrix
	Sigma <- matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))
	for (i in 1:nrow(Sigma)) {
		for (j in 1:ncol(Sigma)) {
			Sigma[i, j] <- sig ** 2 * exp(-0.5 * (abs(X1[i] - X2[j]) / l) ^ 2)
		}
	}
	return(Sigma)
}


periodic_kernel <- function(X1,
							X2,
							l = 1,
							p = 12,
							sig = 1) {
	# Periodic-Kernel
	#
	# Parameters:
	#	X1, X2 = vectors
	# 	l = the scale length parameter
	# Returns:
	# 	a covariance matrix
	Sigma <- matrix(0, nrow = length(X1), ncol = length(X2))
	for (i in 1:nrow(Sigma)) {
		for (j in 1:ncol(Sigma)) {
			Sigma[i, j] <-
				sig ** 2 * exp(-2 * (sin (pi * abs(
					X1[i] - X2[j]
				) / p) / l) ** 2)
		}
	}
	# TODO
	#Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
	#isSymmetric(SS)
	
	return(Sigma)
	
}



#Matrix::forceSymmetric()

locally_periodoc_kernel <- function(X1,
									X2,
									l = 1,
									p = 12,
									sig = 1) {
	# Locally Periodic-Kernel
	#
	# Parameters:
	#	X1, X2 = vectors
	# 	l = the scale length parameter
	# Returns:
	# 	a covariance matrix
	Sigma <- matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))
	for (i in 1:nrow(Sigma)) {
		for (j in 1:ncol(Sigma)) {
			Sigma[i, j] <- sig ** 2 * exp(-0.5 * (abs(X1[i] - X2[j]) / l) ^ 2)
		}
	}
	return(Sigma)
	
}


l <- 20
sig <- 2.5
sig_white <- 0
kernel <- rbf_kernel

TT <- nrow(rw_log)
xx <-
	(as.numeric(index(rw_log)) - min(as.numeric(index(rw_log))))[sample(1:TT, size =
																			100, replace = FALSE)]
xx <- xx - max(xx) / 2
yy <-
	coredata(rw_log)[, , drop = TRUE][sample(1:TT, size = 100, replace = FALSE)]
yy <- yy - mean(yy)
f <- data.frame(x = xx,
				y = yy)

# Define the points at which we want to define the functions
x.star <- seq(min(xx) * 1.1, max(xx) * 1.1, len = 400)

# Calculate the covariance matrix
sigma <- kernel(x.star, x.star, l = l, sig = sig)


# Calculate the covariance matrices
# using the same x.star values as above
x <- f$x
k.xx <- kernel(x, x, l = l, sig = sig)
k.xxs <- kernel(x, x.star, l = l, sig = sig)
k.xsx <- kernel(x.star, x, l = l, sig = sig)
k.xsxs <- kernel(x.star, x.star, l = l, sig = sig)

# These matrix calculations correspond to equation (2.19)
# in the book.
f.star.bar <-
	k.xsx %*% solve(k.xx + sig_white * diag(nrow(k.xx)), diag(nrow(k.xx))) %*% f$y
#abc <- solve(k.xx %*% f$y, )
cov.f.star <-
	k.xsxs - k.xsx %*% solve(k.xx, diag(nrow(k.xx))) %*% k.xxs

# This time we'll plot more samples.  We could of course
# simply plot a +/- 2 standard deviation confidence interval
# as in the book but I want to show the samples explicitly here.
n.samples <- 25
values <-
	matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
for (i in 1:n.samples) {
	LL <- chol(cov.f.star + 0.01 * diag(nrow(cov.f.star)))
	
	#values[, i] <- mvrnorm(1, f.star.bar, cov.f.star)
	values[, i] <-
		c(f.star.bar + LL %*% zrnorm(nrow(f.star.bar)))#mvrnorm(1, f.star.bar, cov.f.star)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id = "x")

# Plot the results including the mean function
# and constraining data points
x_star <-
	tibble(x_star = x.star, f_star = f.star.bar[, , drop = TRUE])
ggplot(values, aes(x = x, y = value)) +
	geom_line(aes(group = variable), colour = "#013848", alpha = .1) +
	geom_line(
		data = x_star,
		aes(x = x_star, y = f_star),
		colour = "#013848",
		size = 1
	) +
	geom_point(data = f, aes(x = x, y = y)) +
	theme_minimal() +
	scale_y_continuous(lim = c(-5, 5), name = "output, f(x)") +
	xlab("input, x")
