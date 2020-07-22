library(tidyverse)
library(mvnfast)
library(BayesMRA)
library(spam)
library(patchwork)

## test different ways of sampling from Bayesian Gaussian full conditionals

## dense matrices
set.seed(111)
m <- 4
tmp=matrix(rnorm(m^2), m, m)
A=t(tmp)%*% tmp + diag(m)
A_inv = solve(A)
b <- rnorm(m) + (1:m)*4
N <- 100

dat_plot <- data.frame(
    samples_mvnfast = c(rmvn(N, A_inv %*% b, A_inv)),
    samples_spam    = c(rmvnorm.canonical(N, b, A)),
    samples_arma    = c(t(replicate(N, as.vector(rmvn_arma(A, b))))),
    iterations = 1:N,
    par = rep(1:m, each = N)
)

dat_plot %>%
    pivot_longer(cols = c("samples_mvnfast", "samples_spam", "samples_arma"), names_to = "method", values_to = "samples") %>%
    ggplot(aes(x = samples)) +
    geom_histogram() +
    facet_grid(method ~ par)


## sparse matrices
set.seed(111)
m <- 3^2 ## make sure m is a square
A = make_Q_alpha_2d(sqrt(m), 0.9)[[1]]
A_inv = solve(A)
b <- rnorm(m) + (1:m)*4
N <- 5000

## generate samples
method_list     <- c("samples_mvnfast", "samples_spam", "samples_arma", "samples_wrong")
samples_mvnfast <- rmvn(N, A_inv %*% b, A_inv)
samples_spam    <- rmvnorm.canonical(N, b, A)
samples_arma    <- t(replicate(N, as.vector(rmvn_arma(as.matrix(A), b))))
samples_wrong   <- rmvn(N, b, A)

dat_plot <- data.frame(
    samples_mvnfast = c(samples_mvnfast),
    samples_spam    = c(samples_spam),
    samples_arma    = c(samples_arma),
    samples_wrong    = c(samples_wrong),
    iterations = 1:N,
    par = rep(1:m, each = N)
) %>%
    pivot_longer(cols = all_of(method_list), names_to = "method", values_to = "samples")

dat_plot %>%
    ggplot(aes(x = samples)) +
    geom_histogram() +
    facet_grid(method ~ par)



## calculate statistics
p1 <- data.frame(
    truth     = c(A_inv %*% b),
    estimates = c(
        colMeans(samples_mvnfast),
        colMeans(samples_spam),
        colMeans(samples_arma),
        colMeans(samples_wrong)
        ),
    method = rep(method_list, each = m)
) %>%
    ggplot(aes(x = truth, y = estimates)) +
    geom_point() +
    facet_wrap(~ method, scales = "free") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    ggtitle("Estimates of the mean")

p2 <- data.frame(
    truth     = c(A_inv),
    estimates = c(
        c(cov(samples_mvnfast)),
        c(cov(samples_spam)),
        c(cov(samples_arma)),
        c(cov(samples_wrong))
    ),
    method = rep(method_list, each = m^2)
) %>%
    ggplot(aes(x = truth, y = estimates)) +
    geom_point() +
    facet_wrap(~ method, scales = "free") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    ggtitle("Estimates of the covariance")

p1 / p2

