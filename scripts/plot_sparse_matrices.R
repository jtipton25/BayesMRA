library(BayesMRA)
library(tidyverse)

base_size <- 24
M <- 3
locs <- as.matrix(expand.grid(seq(0, 1, length.out = 40), seq(0, 1, length.out = 40)))
MRA <- mra_wendland_2d(locs, M = 3, max_points = 128800)
Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), rep(1, MRA$M))
W <- MRA$W
pW <- plot_sparse_matrix(W,
                         title="Basis Functions W",
                         base_size = base_size)
# p1 <- plot_sparse_matrix(Q[[1]]) +
#     ggtitle("First resolution") +
#     theme(plot.title = element_text(hjust = 0.5))
# p2 <- plot_sparse_matrix(Q[[2]]) +
#     ggtitle("Second resolution") +
#     theme(plot.title = element_text(hjust = 0.5))
# p3 <- plot_sparse_matrix(Q[[3]]) +
#     ggtitle("Third resolution") +
#     theme(plot.title = element_text(hjust = 0.5))
# (p1 + p2) / p3
# p1 + p2 + p3
# p1 / p2 / p3

pQ <- plot_sparse_matrix(do.call(spam::bdiag.spam, Q),
                         title = "Precision Matrix Q",
                         base_size = base_size)

G <- t(W) %*% W + do.call(spam::bdiag.spam, Q)

pG <- plot_sparse_matrix(G,
                         title = "Kernel Matrix G",
                         base_size = base_size)

G_chol <- chol(G)
p_chol <- plot_sparse_matrix(spam::as.spam(G_chol),
                             title = "Cholesky Factor of G",
                             base_size = base_size)


ggsave(here::here("images", "sparse-matrix.png"),
       (pW + pQ) / (pG + p_chol),
       device = "png",
       width = 16,
       height = 9,
       units = "in")
