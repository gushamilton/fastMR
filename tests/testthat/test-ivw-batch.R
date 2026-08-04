test_that("batched IVW variants agree with the scalar mixed-method path", {
  g <- grid_fixture(3, 2)
  direct <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                         methods = c("ivw", "ivw_fe", "ivw_mre"), nboot = 0, threads = 5)
  mixed <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                        methods = c("ivw", "ivw_fe", "ivw_mre", "uwr"),
                        nboot = 0, threads = 1)
  mixed <- mixed[mixed$method_code %in% c("ivw", "ivw_fe", "ivw_mre"), ]
  expect_equal(direct$method_code, mixed$method_code)
  expect_equal(direct$b, mixed$b, tolerance = 1e-14)
  expect_equal(direct$se, mixed$se, tolerance = 1e-14)
  expect_equal(direct$pval, mixed$pval, tolerance = 1e-14)
})

test_that("batched IVW handles zero and signed exposure effects", {
  g <- grid_fixture(2, 2)
  g$exposure_beta[1, 1] <- 0
  g$exposure_beta[1, 2] <- -abs(g$exposure_beta[1, 2])
  g$exposure_beta[2, 3] <- 1e-15
  result <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                          methods = "ivw", nboot = 0, threads = 5)
  expect_equal(nrow(result), 4)
  expect_true(all(is.finite(result$b)))
  expect_true(all(is.finite(result$se)))
})

test_that("batched IVW preserves an exact zero-residual fit", {
  exposure_beta <- matrix(
    c(-0.18857602, 0.06684381, -0.12420417, -0.52960750,
      0.08635273, 0.08257450),
    nrow = 1L
  )
  outcome_beta <- exposure_beta
  exposure_se <- matrix(
    c(0.04247208, 0.01420307, 0.02663656, 0.11156014,
      0.01890862, 0.01566239),
    nrow = 1L
  )
  outcome_se <- exposure_se

  batched <- fast_mr_grid(
    exposure_beta, outcome_beta, exposure_se, outcome_se,
    methods = "ivw", nboot = 0L, threads = 1L
  )
  scalar <- fast_mr_grid(
    exposure_beta, outcome_beta, exposure_se, outcome_se,
    methods = c("ivw", "uwr"), nboot = 0L, threads = 1L
  )
  scalar <- scalar[scalar$method_code == "ivw", , drop = FALSE]

  expect_equal(batched$b, 1, tolerance = 1e-15)
  expect_equal(batched$Q, 0, tolerance = 0)
  expect_equal(batched$se, 0, tolerance = 0)
  expect_equal(batched[c("b", "se", "Q")], scalar[c("b", "se", "Q")],
               tolerance = 0)
})
