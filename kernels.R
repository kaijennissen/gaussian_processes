# library(RcppZiggurat)
# library(Matrix)
# library(xts)
library(ggplot2)

# SS <- periodic_kernel(c(1:5), c(1:5))
# SS <- matrix(0, nrow = 5, ncol = 5)
# SS[lower.tri(SS)] <- 1:10
# # 1)
# SS <- SS + t(SS)
# isSymmetric(SS)
# # 2)
# SS[upper.tri(SS)] <- t(SS)[upper.tri(SS)]
# isSymmetric(SS)


# Squared Exponential Kernel ----------------------------------------------
squared_exponential_kernel <- function(X1, X2, l = 1, sig = 1) {
	# RBF-Kernel
	# The SE kernel has become the de-facto default kernel for GPs and SVMs.
	# This is probably because it has some nice properties. It is universal,
	# and you can integrate it against most functions that you need to.
	# Every function in its prior has infinitely many derivatives.
	# It also has only two parameters:
	# 	The lengthscale l determines the length of the 'wiggles' in your
	# 		function. In general, you won't be able to extrapolate more than l units
	# 		away from your data.
	# 	The output variance sig^2 determines the average distance of your
	# 		function away from its mean. Every kernel has this parameter out
	#		in front; it's just a scale factor.
	# Parameters:
	# 	X1, X2 = vectors
	# 	l = the scale length parameter
	# Returns:
	# 	a covariance matrix
	Sigma <- matrix(0, nrow = length(X1), ncol = length(X2))
	for (i in 1:nrow(Sigma)) {
		for (j in 1:ncol(Sigma)) {
			Sigma[i, j] <- sig ** 2 * exp(-0.5 * (abs(X1[i] - X2[j]) / l) ^ 2)
		}
	}
	return(Sigma)
}


# Rational Quadratic Kernel -----------------------------------------------
rational_quadratic_kernel <-
	function(X1,
			 X2,
			 l = 1,
			 a = 1,
			 sig = 1) {
		# Rational Quadratic Kernel
		#
		# T his kernel is equivalent to adding together many SE kernels with different
		# lengthscales. So, GP priors with this kernel expect to see functions which
		# vary smoothly across many lengthscales.
		# The parameter alpha
		# 		determines the relative weighting of large-scale and small-scale variations.
		# 		When alpha -> inf, the RQ is identical to the SE.
		#
		# Parameters:
		# 	X1, X2 = vectors
		# 	l = the scale length parameter
		# 	a =
		# 	sig =
		# Returns:
		# 	Sigma = covariance matrix
		
		Sigma <-
			matrix(0, nrow = length(X1), ncol = length(X2))
		for (i in 1:nrow(Sigma)) {
			for (j in 1:ncol(Sigma)) {
				Sigma[i, j] <-
					(sig ** 2 * (1 + (X1[i] - X2[j]) ** 2 / (2 * a * l ** 2)) ** (-a))
			}
		}
		return(Sigma)
	}


# Periodic Kernel ---------------------------------------------------------
periodic_kernel <- function(X1,
							X2,
							l = 1,
							p = pi,
							sig = 1) {
	# Periodic-Kernel
	#
	# The periodic kernel allows one to model functions which repeat themselves
	# exactly. Its parameters are easily interpretable:
	# 	The period p
	# 		simply determines the distance between repititions of the function.
	# The lengthscale l
	# 		determines the lengthscale function in the same way as in the SE kernel.
	#
	# Parameters:
	# 	X1, X2 = vectors
	# 	l = the scale length parameter
	# Returns:
	# 	Sigma = covariance matrix
	Sigma <-
		matrix(0, nrow = length(X1), ncol = length(X2))
	for (i in 1:nrow(Sigma)) {
		for (j in 1:ncol(Sigma)) {
			Sigma[i, j] <-
				(sig ** 2 * exp(-2 * (sin (
					pi * abs(X1[i] - X2[j]) / p
				) / l) ** 2))
		}
	}
	
	# TODO
	# Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
	# isSymmetric(SS)
	
	return(Sigma)
}


# Locally Periodic Kernel -------------------------------------------------
locally_periodoc_kernel <- function(X1,
									X2,
									...) {
	# Locally Periodic-Kernel
	#
	# A SE kernel times a periodic results in functions which are periodic,
	# but which can slowly vary over time.
	# Most periodic functions don't repeat themselves exactly.
	# To add some flexibility to our model, we can consider adding or multiplying
	# a local kernel such as the squared-exp with our periodic kernel.
	# This will allow us to model functions that are only locally periodic -
	# the shape of the repeating part of the function can now change over time.
	#
	# Parameters:
	# 	X1, X2 = vectors
	# 	l = the scale length parameter
	# Returns:
	# 	a covariance matrix
	dots <- list(...)
	Sigma <-
		periodic_kernel(
			X1 = X1,
			X2 = X2,
			l = dots$l1,
			p = dots$p,
			sig = dots$sig1
		) * squared_exponential_kernel(
			X1 = X1,
			X2 = X2,
			l = dots$l2,
			sig = dots$sig2
		)
	return(Sigma)
}


# Linear Kernel -----------------------------------------------------------
linear_kernel <- function(X1,
						  X2,
						  c = 1,
						  sig_b = 1,
						  sig_v = 1) {
	# Linear-Kernel
	#
	#T he linear kernel is not like the others in that it's non-stationary.
	# A stationary covariance function is one that only depends on the relative
	# position of its two inputs, and not on their absolute location.
	# That means that the parameters of the linear kernel are about specifying the origin:
	# The offset c
	# 		determines the x-coordinate of the point that all the lines in the
	# 		posterior go though. At this point, the function will have zero
	#		variance (unless you add noise)
	# The constant variance sig_b^2
	# 		determines how far from 0 the height of the function will be at zero.
	# 		It's a little confusing, becuase it's not specifying that value directly,
	# 		but rather putting a prior on it. It's equivalent to adding an
	# 		uncertain offset to our model.
	#
	# Parameters:
	# 	X1, X2 = vectors
	# 	l = the scale length parameter
	# Returns:
	# 	a covariance matrix
	Sigma <-
		matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))
	for (i in 1:nrow(Sigma)) {
		for (j in 1:ncol(Sigma)) {
			Sigma[i, j] <- sig_b ** 2 + sig_v ** 2 * (X1[i] - c) * (X2[j] - c)
		}
	}
	return(Sigma)
}


# kernel plots ------------------------------------------------------------

fit_plot_kernel <- function(df, x.star, n_samples = 5, kernel, ...) {
	sigma <- kernel(x.star, x.star, ...)
	
	x <- df$x
	k.xx <- kernel(x, x, ...)
	k.xxs <- kernel(x, x.star, ...)
	k.xsx <- kernel(x.star, x, ...)
	k.xsxs <- kernel(x.star, x.star, ...)
	
	f.star.bar <- k.xsx %*% solve(k.xx) %*% df$y
	cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx) %*% k.xxs
	
	n.samples <- n_samples
	values <-
		matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
	for (i in 1:n.samples) {
		values[, i] <- MASS::mvrnorm(1, f.star.bar, cov.f.star)
	}
	
	# values <- cbind(x = x.star,
	# 				#as.data.frame(values),
	# 				ymax=apply(values, MARGIN=1, FUN=max, na.rm=TRUE),
	# 				ymin=apply(values, MARGIN=1, FUN=min, na.rm=TRUE)) %>% as.data.frame()
	# values$abs_max <- abs(values$ymax-f.star.bar)
	# values$abs_min <- abs(values$ymin-f.star.bar)
	# values$y_bar <- f.star.bar
	# 
	# ggplot(values,aes(x = x)) +
	# 	geom_ribbon(aes(ymax=ymax,ymin=ymin), alpha=.25, col="#013848")+
	# 	geom_line(aes(y = y_bar), col="red")+
	# 	geom_point(data = f, aes(x = x, y = y))
		#geom_line(aes(x=x,y=ymax))+
		# geom_line(aes(x=x,y=ymin))+

	values <- cbind(x = x.star,
					as.data.frame(values))
	values <- reshape2::melt(values, id = "x")
	values$col <- abs(values$value -rep(f.star.bar,n_samples))
	#values$col <- (values$col-mean(values$col))/sd(values$col)
	
	x_star <-
		tidyr::tibble(x_star = x.star, f_star = f.star.bar[, , drop = TRUE])
	 fig1 <- 
		ggplot(values, aes(x = x, y = value)) +
		geom_line(aes(group = variable, col=col,alpha=col))+
		scale_color_viridis_c(option = "viridis")+
		#scale_color_gradient(low="#09557f", high="#00a378")+
		scale_alpha_continuous(range=c(.5,.2))+
		geom_line(
			data = x_star,
			aes(x = x_star, y = f_star),
			colour = "#C62F4B", 
			size = .5
		) +
		geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
		theme_minimal() +
		scale_y_continuous(lim = c(1.5 * min(f.star.bar), 1.5 * max(f.star.bar)), name = "output, f(x)") +
		xlab("input, x")
	return(fig1)
}


# examples ----------------------------------------------------------------

# f <- data.frame(x = c(-5, -3, -2, -1, -0.5, 3),
# 				y = c(-2, 0,-1,.5, 2, -1))

x.star <- seq(-20, 20, len = 200)
xx <- seq(from = -5 * pi,
		  to = 5 * pi,
		  by = .5)
yy <- sin(xx)
ff <- data.frame(x = xx,
				 y = yy)

x <- sample(xx, size = 15, replace = FALSE)
y <- sin(x)
f <- data.frame(x = x,
				y = y)

fig1 <-
	fit_plot_kernel(
		df = f,
		x.star = x.star,
		n_sample = 100,
		kernel = rational_quadratic_kernel,
		l = 3,
		sig = 1
	)
fig1 +
	geom_point(data = ff,
			   aes(x = xx, y = yy),
			   color = "red",
			   size = .2)










# rw data ------------------------------------------------------------
# 
# library(tidyverse)
# library(xts)
# # rm(list = ls())
# rw_log <- read_csv2("rw_log_views.csv") %>% timetk::tk_xts()
# 
# l1 <- 10
# l2 <- 10
# sig1 <- 2
# sig2 <- 2
# p <- 1
# sig_white <- 0
# 
# TT <- nrow(rw_log)
# xx <-
# 	(as.numeric(index(rw_log)) - min(as.numeric(index(rw_log))))[sample(1:TT,
# 																		size =
# 																			20, replace = FALSE)]
# xx <- xx - max(xx) / 2
# yy <-
# 	coredata(rw_log)[, , drop = TRUE][sample(1:TT, size = 20, replace = FALSE)]
# yy <- yy - mean(yy)
# f <- data.frame(x = xx,
# 				y = yy)
# 
# x.star <- seq(min(xx) * 1.1, max(xx) * 1.1, len = 400)
# sigma <-
# 	locally_periodoc_kernel(
# 		x.star,
# 		x.star,
# 		l1 = l,
# 		sig1 = sig,
# 		p = p,
# 		l2 = l,
# 		sig2 = sig
# 	)
# 
# x <- f$x
# k.xx <-
# 	locally_periodoc_kernel(
# 		x,
# 		x,
# 		l1 = l,
# 		sig1 = sig,
# 		p = p,
# 		l2 = l,
# 		sig2 = sig
# 	)
# k.xxs <-
# 	locally_periodoc_kernel(
# 		x,
# 		x.star,
# 		l1 = l,
# 		sig1 = sig,
# 		p = p,
# 		l2 = l,
# 		sig2 = sig
# 	)
# k.xsx <-
# 	locally_periodoc_kernel(
# 		x.star,
# 		x,
# 		l1 = l,
# 		sig1 = sig,
# 		p = p,
# 		l2 = l,
# 		sig2 = sig
# 	)
# k.xsxs <-
# 	locally_periodoc_kernel(
# 		x.star,
# 		x.star,
# 		l1 = l,
# 		sig1 = sig,
# 		p = p,
# 		l2 = l,
# 		sig2 = sig
# 	)
# 
# 
# K_inv <-
# 	solve(k.xx + sig_white * diag(nrow(k.xx)), diag(nrow(k.xx)))
# f.star.bar <-
# 	k.xsx %*% K_inv %*% f$y
# cov.f.star <-
# 	k.xsxs - k.xsx %*%  K_inv %*% k.xxs
# 
# n.samples <- 25
# values <-
# 	matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
# for (i in 1:n.samples) {
# 	LL <- chol(cov.f.star + 0.01 * diag(nrow(cov.f.star)))
# 	values[, i] <- MASS::mvrnorm(1, f.star.bar, cov.f.star)
# }
# values <- cbind(x = x.star, as.data.frame(values))
# values %>% as_tibble()
# values <- reshape2::melt(values, id = "x")
# 
# x_star <-
# 	tibble(x_star = x.star, f_star = f.star.bar[, , drop = TRUE])
# ggplot(values, aes(x = x, y = value)) +
# 	geom_line(aes(group = variable), colour = "#013848", alpha = .1) +
# 	geom_line(
# 		data = x_star,
# 		aes(x = x_star, y = f_star),
# 		colour = "#013848",
# 		size = 1
# 	) +
# 	geom_point(data = f, aes(x = x, y = y)) +
# 	theme_minimal() +
# 	scale_y_continuous(lim = c(-5, 5), name = "output, f(x)") +
# 	xlab("input, x")
