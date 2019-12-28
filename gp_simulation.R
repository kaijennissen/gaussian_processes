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


# gp simulation from periodic kernel ---------------------------
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


# source("./kernels.R")
# x.star <- seq(-2*pi, 2*pi, len = 800)
# xx <- seq(from = -2 * pi,
#           to = 2 * pi,
#           by = .25)
# yy <- sim_function(xx)
# ff <- data.frame(x = xx,
#                  y = yy)
# x <- sort(sample(xx, size = 15, replace = FALSE))
# y <- sim_function(x)
# f <- data.frame(x = c(-pi/2, -pi/4,0,pi/4,pi/2,3*4*pi,pi),
#                 y = sim_function(c(-pi/2, -pi/4,0,pi/4,pi/2,3*4*pi,pi)))

fig9 <-
    fit_color_plot(
        df = f,
        x.star = x.star,
        n_sample = 50,
        kernel = periodic_kernel,
        l=1,
        sig=1,
        p=2*pi,
        sigma_wn=0.00001
    )
fig9+
    geom_point(data = f, aes(x = x, y = y),color="#C62F4B" ) +
    geom_point(data = ff,
               aes(x = xx, y = yy),
               color = "red",
               fill="red",
               size = .2)


x.star <- seq(-2*pi, 2*pi, len = 800)
xx <- seq(from = -2 * pi,
          to = 2 * pi,
          by = .25)
yy <- sim_function(xx)
ff <- data.frame(x = xx,
                 y = yy)
x <- sort(sample(xx, size = 25, replace = FALSE))
y <- sim_function(x)
f <- data.frame(x = x,
                y = y)
n_samples <- 50
sig_noise <- 0.5

sigma <- periodic_kernel(x.star, x.star, p=pi, l=pi, sig=1)

k.xx <- periodic_kernel(x, x,  p=pi, l=1, sig=1)
k.xxs <- periodic_kernel(x, x.star,  p=pi, l=1, sig=1)
k.xsx <- periodic_kernel(x.star, x,  p=pi, l=1, sig=1)
k.xsxs <- periodic_kernel(x.star, x.star,  p=pi, l=1, sig=1)

inv_sigma <- solve(k.xx+.01*diag(nrow(k.xx)), diag(nrow(k.xx)))
f.star.bar <- k.xsx %*% inv_sigma  %*% y
cov.f.star <- k.xsxs - k.xsx %*% inv_sigma  %*% k.xxs

values <-
    matrix(rep(0, length(x.star) * n_samples), ncol = n_samples)
for (i in 1:n_samples) {
    L <- chol(cov.f.star+sig_noise*diag(nrow(cov.f.star)))
    values[, i] <- f.star.bar+L%*%rnorm(nrow(f.star.bar))  # MASS::mvrnorm(1, f.star.bar, cov.f.star)
}

values_abs <-
    apply(
        values,
        MARGIN = 2,
        FUN = function(x, y) {
            abs(x - y)
        },
        y = f.star.bar
    )
values_col_max <- apply(values_abs, MARGIN = 1, FUN = max)
values_col_min <- apply(values_abs, MARGIN = 1, FUN = min)
values_col_range <- values_col_max - values_col_min
values_col <-
    (values_abs - matrix(rep(values_col_min, n_samples), ncol = n_samples)) /
    matrix(rep(values_col_range, n_samples), ncol = n_samples)

values <- cbind(x = x.star,
                as.data.frame(values))
values_col <- cbind(x = x.star,
                    as.data.frame(values_col))
values <- reshape2::melt(values, id = "x")
values_col <- reshape2::melt(values_col, id = "x")

values$col <- values_col$value

x_star <-
    tidyr::tibble(x_star = x.star, f_star = f.star.bar[, , drop = TRUE])
fig1 <-
    ggplot(values, aes(x = x, y = value)) +
    geom_line(aes(
        group = variable,
        col = col,
        alpha = col
    )) +
    scale_color_viridis_c(option = "viridis") +
    theme_void() +
    theme(legend.position = "none") +
    scale_y_continuous(lim = c(1.1 * min(f.star.bar), 1.1 * max(f.star.bar)), name = "output, f(x)")
fig1

mape <- function(y_pred, y_true){
    return(mean(abs((y_pred-y_true)/y_true)))
}

h <- 36
y <- log(AirPassengers)
TT <- length(X)
y_train <- window(y, end=c(1957,12))
y_test <- window(y, start=c(1958,1))

y_pred <- auto.arima(y_train) %>% forecast(h=36)
y_pred <- tslm(y_train~trend+season) %>% forecast(h=36)


mape(exp(y_pred$mean), exp(y_test))
cbind(y_pred$mean, y_test)

library(forecast)









