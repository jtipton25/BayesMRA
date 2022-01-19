library(tidyverse)
library(splines)
set.seed(404)
N <- 100
X <- cbind(1, rnorm(N))
Z <- bs(seq(0, 1, length.out = N), df = 8)
beta <- rnorm(ncol(X))
alpha <- rnorm(ncol(Z), 0, 2)
epsilon <- rnorm(N, 0, 0.5)
y <- X %*% beta + Z %*% alpha + epsilon

layout(matrix(1:2, 2, 1))
plot(X[, 2], y - Z %*% alpha)
abline(beta)
plot(seq(0, 1, length.out = N), y)
lines(seq(0, 1, length.out = N), Z %*% alpha)

U <- cbind(X, Z)
theta <- c(beta, alpha)
all.equal(U %*% theta, X %*% beta + Z %*% alpha)




# gradient
target_fun <- function(y, U, theta) {
    # U = cbind(X, W)
    # theta <- c(beta, alpha)
    return(sum((y - U %*% theta)^2) / (2 * length(y)))
}

# gradient function
gradient_fun <- function(y, U, theta) {
    # U = cbind(X, W)
    # theta <- c(beta, alpha)
    return(t(U) %*% (U %*% theta - y) / length(y))
}

# Defining gradient Descent function

regression_gradient_descent <- function(target_fun, gradient_fun, y, X, W, inits, threshold = 0.01, alpha = 0.001, num_iters = 10, print_every = 1,
                                        minibatch_size = NULL){

    # center and scale the response y and covariates X
    # improves computational stability
    # y_scaled <- scale(y)
    # y_mean <- attr(y_scaled, "scaled:center")
    # y_sd   <- attr(y_scaled, "scaled:scale")

    # remove the intercept term for the scaling
    # X_scaled <- scale(X[, -1])
    # X_mean <- attr(X_scaled, "scaled:center")
    # X_sd   <- attr(X_scaled, "scaled:scale")
    # X_scaled <- cbind(1, X_scaled)
    # U <- cbind(X_scaled, W)
    U <- cbind(X, W)


    # initialize algorithm
    loss <- rep(0, num_iters)
    beta_save <- matrix(0, nrow = num_iters, ncol = ncol(U))
    colnames(beta_save) <- paste0("beta[", 0:(ncol(U)-1),"]")
    beta_old <- inits

    loss[1] <- target_fun(y, U, beta_old)
    beta_save[1, ] <- beta_old
    # first iteration
    if (is.null(minibatch_size)) {
        beta_new <- beta_old - alpha * gradient_fun(y, U, beta_old)
    } else {
        idx <- sample(1:length(y), minibatch_size)
        beta_new <- beta_old - alpha * gradient_fun(y[idx], U[idx, ], beta_old)
    }

    i <- 2
    print_every <- ifelse(print_every > 0, print_every, num_iters)

    while( any(abs(beta_new - beta_old) > threshold) && i <= num_iters){
        beta_old <- beta_new
        if (is.null(minibatch_size)) {
            beta_new <- beta_new - alpha * gradient_fun(y, U, beta_new)
        } else {
            idx <- sample(1:length(y), minibatch_size)
            beta_new <- beta_new - alpha * gradient_fun(y[idx], U[idx, ], beta_new)
        }
        loss[i] <- target_fun(y, U, beta_new)
        # save the estimates on the original data scale (not normalized scale)
        # beta_save[i, ] <- c(beta_new[1] * y_sd + y_mean - sum(beta_new[-1][1:(ncol(X)-1)] * X_mean / X_sd) * y_sd,
        #                     beta_new[-1] * y_sd / X_sd)
        beta_save[i, ] <- beta_new
        if(i %% print_every == 0) {
            message("iteration = ", i, ", loss = ", loss[i])
        }
        i <- i+1
    }

    cbind(as.data.frame(beta_save[1:(i-1), ]), loss = loss[1:(i-1)], iteration = 1:(i-1))
}

# gradient descent function for regression ----

dat <- regression_gradient_descent(target_fun, gradient_fun, y, X, Z, inits = rnorm(ncol(U)), threshold = 0.000001, alpha = 0.01, num_iters = 50000, print_every = 100)


ggplot(dat, aes(x = iteration, y = loss)) +
    geom_point() +
    geom_line() +
    ggtitle("loss as a function of iteration")


# Plot the regression parameters

# convert beta into a data.frame

dat %>%
    pivot_longer(cols = 1:10, names_to = "parameter", values_to = "beta") %>%
    ggplot(aes(x = iteration, y = beta, color = parameter)) +
    geom_line() +
    theme(legend.position = "none") +
    ggtitle("Parameter estimates by iteration")




model$coefficients
tail(dat, n=1)[1:10]

# p <- ggplot(data = NULL, aes(x = X[, 2], y = y - Z %*% alpha)) +
#     geom_point() +
#     geom_abline(data = dat, aes(intercept = `beta[0]`, slope = `beta[1]`, color = iteration), alpha = 0.25) +
#     # geom_abline(intercept = dat$`beta[0]`, slope = dat$`beta[1]`, color = dat$iteration) +
#     scale_color_viridis_c(direction = -1)
# p


# Do this with BayesMRA
library(BayesMRA)
library(spam)
set.seed(11)

N <- 100^2

## setup the spatial process
locs <- as.matrix(
    expand.grid(
        seq(0, 1, length.out = sqrt(N)),
        seq(0, 1, length.out = sqrt(N))
    )
)
# D <- fields::rdist(locs)

## fixed effects include intercept, elevation, and latitude
X <- cbind(1, rnorm(N), locs[, 2])
# X <- cbind(1, as.vector(mvnfast::rmvn(1, rep(0, N), 3 * exp(-D / 20))), locs[, 2])
p <- ncol(X)

beta <- rnorm(ncol(X))

## MRA spatio-temporal random effect
M <- 3
n_coarse <- 10

MRA    <- mra_wendland_2d(locs, M = M, n_coarse = n_coarse, use_spam = TRUE)

# MRA    <- mra_wendland_2d(locs, M = M, n_max_fine_grid = 2^8, use_spam = TRUE)
W <- MRA$W
n_dims <- MRA$n_dims
dims_idx <- MRA$dims_idx

Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(0.999, length(n_dims)), prec_model = "CAR")

tau2         <- 10 * 2^(2 * (1:M - 1))
Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2)

## initialize the random effect
## set up a linear constraint so that each resolution sums to one
A_constraint <- sapply(1:M, function(i){
    tmp = rep(0, sum(n_dims))
    tmp[dims_idx == i] <- 1
    return(tmp)
})

a_constraint <- rep(0, M)
alpha   <- as.vector(rmvnorm.prec.const(n = 1, mu = rep(0, sum(n_dims)), Q = Q_alpha_tau2, A = t(A_constraint), a = a_constraint))

sigma2 <- runif(1, 0.25, 0.5)

y <- X %*% beta + W %*% alpha + rnorm(N, 0, sqrt(sigma2))


# gradient descent function for MRA ----

U <- cbind(X, W)

dat <- regression_gradient_descent(target_fun, gradient_fun, c(y),
                                   X, W, inits = rnorm(ncol(U)),
                                   threshold = 0.000001,
                                   alpha = 0.1, num_iters = 5000, print_every = 100,
                                   minibatch_size = 2^6)


ggplot(dat, aes(x = iteration, y = loss)) +
    geom_point() +
    geom_line() +
    ggtitle("loss as a function of iteration")


# Plot the regression parameters

# convert beta into a data.frame
# results from lm()

dat %>%
    pivot_longer(cols = 1:10, names_to = "parameter", values_to = "beta") %>%
    ggplot(aes(x = iteration, y = beta, color = parameter)) +
    geom_line() +
    ggtitle("Parameter estimates by iteration")


layout(matrix(1:2, 2, 1))
plot(y, cbind(X, W) %*% unlist(dat[nrow(dat), ][1:(ncol(dat) - 2)]))
plot(y, X %*% beta + W %*% alpha)

# fitted MSE and R^2
preds <- cbind(X, W) %*% unlist(dat[nrow(dat), ][1:(ncol(dat) - 2)])
mean((y - preds)^2)
cor(y, preds)^2
# oracle MSE and R^2
mu <- X %*% beta + W %*% alpha
mean((y - mu)^2)
cor(y, mu)^2
