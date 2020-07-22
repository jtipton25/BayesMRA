## investigate prior on tau2
lambda <- rgamma(1000, 0.5, 1/600)
tau2   <- 1 / rgamma(1000, 0.5, lambda)
hist(1 / tau2[1 / tau2 < 100], breaks = 2000, freq = FALSE)
curve(LaplacesDemon::dhalfcauchy(x, scale = sqrt(6)), from = 0, to = 100, add = TRUE, col = "red")
