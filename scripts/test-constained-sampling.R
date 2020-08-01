A_constraint <- matrix(0, nrow(A_alpha), M)
A_constraint <- t(
    sapply(1:M, function(i){
        tmp = rep(0, sum(n_dims))
        tmp[dims_idx == i] <- 1
        return(tmp)
    })
)
a_constraint <- rep(0, M)
alpha   <- as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha, Rstruct = Rstruct,
                                             A = t(A_constraint), a = a_constraint))
sum(alpha)
sum(alpha[dims_idx == 1])
sum(alpha[dims_idx == 2])
sum(alpha[dims_idx == 3])




##
## test constrained alpha simulation (not full conditional form)
##
A_constraint <- sapply(1:M, function(i){
    tmp = rep(0, sum(n_dims))
    tmp[dims_idx == i] <- 1
    return(tmp)
})

a_constraint <- rep(0, M)
alpha   <- as.vector(rmvnorm.prec.const(n = 1, mu = rep(0, sum(n_dims)), Q = Q_alpha_tau2, A = A_constraint, a = a_constraint))
