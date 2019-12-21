library(tidyverse)
library(ggplot2)
source("./kernels.R")



# gp simulation from squared exponential kernel ---------------------------
sim_function <- function(x){
    return(x*3)
}

x.star <- seq(-20, 20, len = 800)
xx <- seq(from = -5 * pi,
		  to = 5 * pi,
		  by = .5)
yy <- sim_function(xx)
ff <- data.frame(x = xx,
				 y = yy)

x <- sample(xx, size = 10, replace = FALSE)
y <- sim_function(x)
f <- data.frame(x = x,
				y = y)

fig1 <-
	fit_color_plot(
		df = f,
		x.star = x.star,
		n_sample = 25,
		kernel = squared_exponential_kernel,
		l=2,
		sig=1
	)
fig1+
	geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
	geom_point(data = ff,
			   aes(x = xx, y = yy),
			   color = "red",
			   fill="red",
			   size = .2)


fig2 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = squared_exponential_kernel,
        l=2,
        sig=1
    )
fig2+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)


fig3 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = squared_exponential_kernel,
        l=3,
        sig=1
    )
fig3+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)


fig4 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = squared_exponential_kernel,
        l=5,
        sig=1
    )
fig4+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)


fig5 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = squared_exponential_kernel,
        l=1,
        sig=2
    )
fig5+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)

fig6 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = squared_exponential_kernel,
        l=1,
        sig=4
    )
fig6+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)

fig7 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = squared_exponential_kernel,
        l=3,
        sig=4
    )
fig7+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)


# gp simulation from squared exponential kernel ---------------------------
sim_function <- function(x){
    return(sin(x))
}

x.star <- seq(-20, 20, len = 800)
xx <- seq(from = -5 * pi,
          to = 5 * pi,
          by = .5)
yy <- sim_function(xx)
ff <- data.frame(x = xx,
                 y = yy)

x <- sample(xx, size = 10, replace = FALSE)
y <- sim_function(x)
f <- data.frame(x = x,
                y = y)

fig8 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = squared_exponential_kernel,
    )
fig8+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)


fig9 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 25,
        kernel = periodic_kernel,
        sig=.5
    )
fig9+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
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
