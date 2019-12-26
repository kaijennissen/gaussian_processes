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
    Sigma <- matrix(0, nrow = length(X1), ncol = length(X2))
    for (i in 1:nrow(Sigma)) {
        for (j in 1:ncol(Sigma)) {
            Sigma[i, j] <- sig_b ** 2 + sig_v ** 2 * (X1[i] - c) * (X2[j] - c)
        }
    }
    return(Sigma)
}


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
        
        Sigma <- matrix(0, nrow = length(X1), ncol = length(X2))
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
                            sig=1,
                            p=pi,
                            l=1) {
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

    Sigma <- matrix(0, nrow = length(X1), ncol = length(X2))
    for (i in 1:nrow(Sigma)) {
        for (j in 1:ncol(Sigma)) {
            Sigma[i, j] <-
                (sig ** 2 * exp(-2 * (sin (
                    pi * abs(X1[i] - X2[j]) / p
                ) / l) ** 2))
            # periodic_kernel(abs(X1[i] - X2[j]), l=l,p=p,sig=sig)
        }
    }
   
    #Sigma[Sigma <.Machine$double.eps] <- 0
    # TODO
    # Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
    # isSymmetric(SS)
    
    return(Sigma)
}
# 
# periodic_kernel <- function(x,
#                             sig = 1,
#                             l = 1,
#                             p = pi) {
#     y <- sig ** 2 * exp(-2 * (sin (pi * abs(x) / p) / l) ** 2)
# }


# plot length
# tibble(x = seq(-3 * pi, 3 * pi, pi / 64)) %>%
#     mutate(
#         '.5' = periodic_kernel(x, l = .5, sig = 1, p = pi),
#         '1' = periodic_kernel(x, l = 1, sig = 1, p = pi),
#         '2' = periodic_kernel(x, l = 2, sig = 1, p = pi),
#         '4' = periodic_kernel(x, l = 4, sig = 1, p = pi)
#     ) %>%
#     gather(key = "l", value = "y", -x) %>%
#     ggplot(aes(x = x, y = y, col = l)) +
#     geom_line()+
#     geom_vline(xintercept = c(pi, -pi))
# 
# # plot sigma
# tibble(x = seq(-3 * pi, 3 * pi, pi / 64)) %>%
#     mutate(
#         '.1' = periodic_kernel(x, sig = .1, p = pi),
#         '.5' = periodic_kernel(x, sig = .5, p = pi),
#         '1' = periodic_kernel(x, sig = 1, p = pi),
#         '1.5' = periodic_kernel(x, sig = 1.5, p = pi)
#     ) %>%
#     gather(key = "sig", value = "y", -x) %>%
#     ggplot(aes(x = x, y = y, col = sig)) +
#     geom_line()
# 
# # plot p
# tibble(x = seq(-3 * pi, 3 * pi, pi / 64)) %>%
#     mutate(
#         'pi/2' = periodic_kernel(x, p = pi / 2),
#         'pi*3/2' = periodic_kernel(x, p = 3 / 2 * pi),
#         'pi' = periodic_kernel(x, p = pi)
#     ) %>%
#     gather(key = "pi", value = "y", -x) %>%
#     ggplot(aes(x = x, y = y, col = pi)) +
#     geom_line()

# Locally Periodic Kernel -------------------------------------------------
locally_periodic_kernel <- function(X1,
                                    X2,...) {
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




# kernel plots ------------------------------------------------------------
fit_plot <- function(df, x.star, n_samples = 5, kernel, ...) {
    
    sigma <- kernel(x.star, x.star, ...)
    
    x <- df$x
    k.xx <- kernel(x, x, ...)
    k.xxs <- kernel(x, x.star, ...)
    k.xsx <- kernel(x.star, x, ...)
    k.xsxs <- kernel(x.star, x.star, ...)
    
    
    inv_sigma <- solve(k.xx, diag(nrow(k.xx)))
    f.star.bar <- k.xsx %*% inv_sigma %*% df$y
    cov.f.star <- k.xsxs - k.xsx %*% inv_sigma %*% k.xxs
    
    n.samples <- n_samples
    values <-
        matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
    for (i in 1:n.samples) {
        values[, i] <- MASS::mvrnorm(1, f.star.bar, cov.f.star)
    }
    
    values <- cbind(x = x.star,
                    as.data.frame(values))
    values <- reshape2::melt(values, id = "x")
    
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
        geom_line(
            data = x_star,
            aes(x = x_star, y = f_star),
            colour = "#C62F4B",
            size = .5
        ) +
        geom_point(data = f, aes(x = x, y = y), color = "#C62F4B") +
        theme_void() +
        theme(legend.position = "none") +
        scale_y_continuous(lim = c(1.75 * min(f.star.bar), 1.75 * max(f.star.bar)), name = "output, f(x)") +
        xlab("input, x")
    return(fig1)
}


# colored plots -----------------------------------------------------------
fit_color_plot <- function(df, x.star, n_samples = 5, kernel, ...) {
    sigma <- kernel(x.star, x.star, ...)
    
    x <- df$x
    k.xx <- kernel(x, x, ...)
    k.xxs <- kernel(x, x.star, ...)
    k.xsx <- kernel(x.star, x, ...)
    k.xsxs <- kernel(x.star, x.star, ...)
    
    inv_sigma <- solve(k.xx+.Machine$double.eps*diag(nrow(k.xx)), diag(nrow(k.xx)))
    f.star.bar <- k.xsx %*% inv_sigma  %*% df$y
    cov.f.star <- k.xsxs - k.xsx %*% inv_sigma  %*% k.xxs
    
    n.samples <- n_samples
    values <-
        matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)
    for (i in 1:n.samples) {
        values[, i] <- MASS::mvrnorm(1, f.star.bar, cov.f.star)
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
        scale_y_continuous(lim = c(1.75 * min(f.star.bar), 1.75 * max(f.star.bar)), name = "output, f(x)")
    
    return(fig1)
}
