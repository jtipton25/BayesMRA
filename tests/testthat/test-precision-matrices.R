context("testing precision matrices")

# test_that("make_Q", {
#
#     expect_error(make_Q(4, "aaa"), "phi must be between -1 and 1")
#     expect_error(make_Q(4, NA), "phi must be between -1 and 1")
#     expect_error(make_Q(4, -2), "phi must be between -1 and 1")
#     expect_error(make_Q(4, 2), "phi must be between -1 and 1")
#     make_Q(4, 0.9)
#     # expect_error(wendland_basis(1:10, 1:10), "radius must be a single positive numeric value")
#     # expect_error(wendland_basis(-10:10, 5), "d must be nonnegative")
#     # expect_error(wendland_basis(c(1:10, NA), 2), "d must not contain missing values")
#     # expect_equal(wendland_basis(1:5, 3), c(0.37717827566936, 0.013971447441955, 0, 0, 0))
# })

test_that("make_Q_alpha_2d", {

    n_dims <- c(4, 16)
    phi <- rep(1, 3)
    expect_error(make_Q_alpha_2d(n_dims, phi), "n_dims and phi must both be vectors of length M.")

    n_dims <- c(4, 16, 32)
    phi <- rep(1, 2)
    expect_error(make_Q_alpha_2d(n_dims, phi), "n_dims and phi must both be vectors of length M.")

    expect_identical(
        make_Q_alpha_2d(4, 0.5),
        list(new("spam",
                 entries = c(2, -0.5, -0.5, -0.5, 3, -0.5, -0.5,
                             -0.5, 3, -0.5, -0.5, -0.5, 2, -0.5, -0.5, 3, -0.5, -0.5, -0.5,
                             -0.5, 4, -0.5, -0.5, -0.5, -0.5, 4, -0.5, -0.5, -0.5, -0.5, 3,
                             -0.5, -0.5, 3, -0.5, -0.5, -0.5, -0.5, 4, -0.5, -0.5, -0.5, -0.5,
                             4, -0.5, -0.5, -0.5, -0.5, 3, -0.5, -0.5, 2, -0.5, -0.5, -0.5,
                             3, -0.5, -0.5, -0.5, 3, -0.5, -0.5, -0.5, 2),
                 colindices = c(1L,
                                2L, 5L, 1L, 2L, 3L, 6L, 2L, 3L, 4L, 7L, 3L, 4L, 8L, 1L, 5L, 6L,
                                9L, 2L, 5L, 6L, 7L, 10L, 3L, 6L, 7L, 8L, 11L, 4L, 7L, 8L, 12L,
                                5L, 9L, 10L, 13L, 6L, 9L, 10L, 11L, 14L, 7L, 10L, 11L, 12L, 15L,
                                8L, 11L, 12L, 16L, 9L, 13L, 14L, 10L, 13L, 14L, 15L, 11L, 14L,
                                15L, 16L, 12L, 15L, 16L),
                 rowpointers = c(1L, 4L, 8L, 12L, 15L,
                                 19L, 24L, 29L, 33L, 37L, 42L, 47L, 51L, 54L, 58L, 62L, 65L),
                 dimension = c(16L, 16L))
        )
    )

    expect_error(make_Q_alpha_2d(4, 1.5), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(4, -1.5), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, 4), c(0.5, 1.5)), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, 4), c(0.5, NA)), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, 4), c(0.5, "aaa")), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, NA), c(0.5, 0.5)), "n_dims must be a vector of integers of length M.")
    expect_error(make_Q_alpha_2d(c(2, "aaa"), c(0.5, 0.5)), "n_dims must be a vector of integers of length M.")
    expect_error(make_Q_alpha_2d(c(2, 3.5), c(0.5, 0.5)), "n_dims must be a vector of integers of length M.")

    n_dims <- c(4, 16)
    phi <- rep(1, 2)

    expect_error(make_Q_alpha_2d(n_dims, phi, use_spam = "TRUE"), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_2d(n_dims, phi, use_spam = 3.5), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_2d(n_dims, phi, use_spam = NA), "use_spam must be either TRUE or FALSE.")

    expect_error(make_Q_alpha_2d(n_dims, phi, prec_model = "AAA"), 'The only valid options for prec_model are \"CAR\" and \"SAR\".')
    expect_error(make_Q_alpha_2d(n_dims, phi, prec_model = NA), 'The only valid options for prec_model are \"CAR\" and \"SAR\".')
    expect_error(make_Q_alpha_2d(n_dims, phi, prec_model = 32), 'The only valid options for prec_model are \"CAR\" and \"SAR\".')

})


test_that("make_Q_alpha_tau2", {

    n_dims <- c(4, 16)
    phi    <- c(1, 1)
    Q_alpha <- make_Q_alpha_2d(n_dims, phi)
    tau2 <- 1:3
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "Q_alpha must be a list of length M and tau2 must be a positive numeric vector of length M.")
    tau2 <- c(1, NA)
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "tau2 must be a positive numeric vector of length M.")
    tau2 <- c(1, "aaa")
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "tau2 must be a positive numeric vector of length M.")
    tau2 <- c(1, -1)
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "tau2 must be a positive numeric vector of length M.")

    tau2 <- 1:2
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2, use_spam = "TRUE"), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2, use_spam = 3.5), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2, use_spam = NA), "use_spam must be either TRUE or FALSE.")

    expect_identical({
        n_dims <- c(2, 4)
        phi    <- c(1, 1)
        make_Q_alpha_2d(n_dims, phi)
    },
    {
        list(
            new("spam",
                entries = c(2, -1, -1, -1, 2, -1, -1, 2, -1,
                            -1, -1, 2),
                colindices = c(1L, 2L, 3L, 1L, 2L, 4L, 1L, 3L, 4L,
                               2L, 3L, 4L),
                rowpointers = c(1L, 4L, 7L, 10L, 13L),
                dimension = c(4L, 4L)
            ),
            new("spam",
                entries = c(2, -1, -1, -1, 3, -1, -1, -1, 3,
                            -1, -1, -1, 2, -1, -1, 3, -1, -1, -1, -1, 4, -1, -1, -1, -1,
                            4, -1, -1, -1, -1, 3, -1, -1, 3, -1, -1, -1, -1, 4, -1, -1, -1,
                            -1, 4, -1, -1, -1, -1, 3, -1, -1, 2, -1, -1, -1, 3, -1, -1, -1,
                            3, -1, -1, -1, 2),
                colindices = c(1L, 2L, 5L, 1L, 2L, 3L, 6L,
                               2L, 3L, 4L, 7L, 3L, 4L, 8L, 1L, 5L, 6L, 9L, 2L, 5L, 6L, 7L, 10L,
                               3L, 6L, 7L, 8L, 11L, 4L, 7L, 8L, 12L, 5L, 9L, 10L, 13L, 6L, 9L,
                               10L, 11L, 14L, 7L, 10L, 11L, 12L, 15L, 8L, 11L, 12L, 16L, 9L,
                               13L, 14L, 10L, 13L, 14L, 15L, 11L, 14L, 15L, 16L, 12L, 15L, 16L),
                rowpointers = c(1L, 4L, 8L, 12L, 15L, 19L, 24L, 29L, 33L,
                                37L, 42L, 47L, 51L, 54L, 58L, 62L, 65L),
                dimension = c(16L, 16L)
            )
        )
    })

    # locs <- matrix(1:20, 10, 2)
    # MRA <- mra_wendland_2d(locs)
    # locs_pred <- matrix(NA, 20, 2)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs_pred <- matrix(1:30, 10, 3)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs_pred <- matrix("11", 10, 2)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs_pred <- 1:10
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs <- matrix(1:30, 10, 3)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs must be a numeric matrix with N rows and 2 columns")
    #
    # locs <- matrix(1:20, 10, 2)
    # locs_pred <- matrix(1:20, 10, 2)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = 3.5), "use_spam must be either TRUE or FALSE")
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = NA), "use_spam must be either TRUE or FALSE")
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = "aaa"), "use_spam must be either TRUE or FALSE")
    #
    # class(MRA) <- NULL
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), 'MRA must be of class "mra_wendland_2d"')
    #
    # class(MRA) <- "XXX"
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), 'MRA must be of class "mra_wendland_2d"')
})
