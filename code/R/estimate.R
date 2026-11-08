# =============================================================================
#  Empirical functional worst-risk minimization for non-linear models
#
#  Objective (pooled over both environments):
#     L(beta, g) = (1 - gamma) * MSE_obs + gamma * MSE_int
#
#  Three-phase estimator (Section 5 of the thesis):
#    Phase 1  Alternating block coordinate descent over the projection kernels
#             (C1 for X1, C2 for X2), with BFGS and fast cubic-polynomial
#             profiling of the link function (OLS on the pooled data).
#    Phase 2  Indices fixed, expanded in a B-spline basis, and passed to a
#             (cross-validated) Group Lasso that removes unstable variables.
#    Phase 3  Unpenalized natural spline on the X1 index to recover the shape
#             of the link g(x1 -> y).
# =============================================================================

library(splines)
library(grpreg)

#' Fit the non-linear functional worst-risk minimizer.
#'
#' @param dat        output of simulate_fsem().
#' @param tm         time grid.
#' @param gamma      worst-risk trade-off in [0, 1] (weight of the
#'                   interventional environment).
#' @param d          number of basis functions for the kernels.
#' @param spline_df  degrees of freedom of the B-spline / natural spline bases.
#' @param max_iter   number of alternating coordinate descent rounds.
#' @param lambda     Group Lasso penalty: a number, 0 for no penalty (OLS),
#'                   or "cv" to choose it by cross-validation.
#' @param on_improve optional callback(C1, C2, label) called whenever the loss
#'                   improves during Phase 1 (used to record an animation).
#' @return list with the estimated kernels `beta1`, `beta2` and the fitted
#'   link model `g1_model` (an lm on the variable `p1_matching`).
fit_fwrm_nonlinear <- function(dat, tm, gamma = 0.5, d = 5, spline_df = 5,
                               max_iter = 5, lambda = 0.05, on_improve = NULL) {
  dt <- tm[2] - tm[1]
  phi <- make_on_basis(d, tm)
  kernel <- function(C) phi %*% C %*% t(phi)
  normalize <- function(C) C / sqrt(sum(C^2))

  Y_obs <- as.vector(dat$yxo[, 1, ])
  Y_int <- as.vector(dat$yxa[, 1, ])
  Y_pool <- c(Y_obs, Y_int)
  X1_obs <- dat$yxo[, 2, ]
  X1_int <- dat$yxa[, 2, ]
  X2_obs <- dat$yxo[, 3, ]
  X2_int <- dat$yxa[, 3, ]
  n_obs <- length(Y_obs)

  # Pooled index of a kernel over both environments.
  index1 <- function(b) c(project(b, X1_obs, dt), project(b, X1_int, dt))
  index2 <- function(b) c(project(b, X2_obs, dt), project(b, X2_int, dt))

  # Worst-risk loss with the link profiled out by a cubic polynomial in each index.
  profiled_loss <- function(p1, p2) {
    X_design <- cbind(1, p1, p1^2, p1^3, p2, p2^2, p2^3)
    coefs <- .lm.fit(X_design, Y_pool)$coefficients
    coefs[is.na(coefs)] <- 0
    y_hat <- as.vector(X_design %*% coefs)
    mse_obs <- mean((Y_obs - y_hat[1:n_obs])^2)
    mse_int <- mean((Y_int - y_hat[(n_obs + 1):length(y_hat)])^2)
    (1 - gamma) * mse_obs + gamma * mse_int
  }

  # Random unit-norm initialization (C1 first, then C2).
  C1 <- normalize(matrix(rnorm(d * d), d, d))
  C2 <- normalize(matrix(rnorm(d * d), d, d))

  # Square-root environment weights for the weighted Group Lasso in Phase 2.
  sqrt_w <- sqrt(c(rep(1 - gamma, n_obs), rep(gamma, length(Y_int))))
  y_weighted <- Y_pool * sqrt_w

  best_loss <- c(Inf, Inf)

  # ---------------------------------------------------------------------------
  # Phase 1: alternating block coordinate descent with BFGS
  # ---------------------------------------------------------------------------
  for (iter in seq_len(max_iter)) {
    p2_curr <- index2(kernel(C2))

    loss_C1 <- function(params) {
      C1_tmp <- normalize(matrix(params, d, d))
      val <- profiled_loss(index1(kernel(C1_tmp)), p2_curr)
      if (!is.null(on_improve) && val < best_loss[1] - 1e-4) {
        best_loss[1] <<- val
        on_improve(C1_tmp, C2, sprintf("Iter %d, Opt C1", iter))
      }
      val
    }
    opt1 <- optim(par = as.vector(C1), fn = loss_C1, method = "BFGS",
                  control = list(maxit = 300))
    C1 <- normalize(matrix(opt1$par, d, d))
    p1_curr <- index1(kernel(C1))

    loss_C2 <- function(params) {
      C2_tmp <- normalize(matrix(params, d, d))
      val <- profiled_loss(p1_curr, index2(kernel(C2_tmp)))
      if (!is.null(on_improve) && val < best_loss[2] - 1e-4) {
        best_loss[2] <<- val
        on_improve(C1, C2_tmp, sprintf("Iter %d, Opt C2", iter))
      }
      val
    }
    opt2 <- optim(par = as.vector(C2), fn = loss_C2, method = "BFGS",
                  control = list(maxit = 300))
    C2 <- normalize(matrix(opt2$par, d, d))

    cat(sprintf("Iteration %d: C1 loss = %f | C2 loss = %f\n",
                iter, opt1$value, opt2$value))
  }

  # ---------------------------------------------------------------------------
  # Phase 2: variable selection with the Group Lasso
  # ---------------------------------------------------------------------------
  b1 <- kernel(C1)
  b2 <- kernel(C2)
  B1_final <- suppressWarnings(bs(index1(b1), df = spline_df))
  B2_final <- suppressWarnings(bs(index2(b2), df = spline_df))

  X_weighted <- cbind(sqrt_w, B1_final * sqrt_w, B2_final * sqrt_w)
  groups <- c(0, rep(1, spline_df), rep(2, spline_df))

  # Standardize the spline columns so both groups are penalized fairly.
  X_model <- function() cbind(X_weighted[, 1], scale(X_weighted[, -1]))

  final_model <- tryCatch(
    {
      if (identical(lambda, "cv")) {
        cv.grpreg(X_model(), y_weighted, group = groups,
                  penalty = "grLasso", intercept = FALSE)
      } else if (lambda == 0) {
        lm(y_weighted ~ 0 + X_weighted)
      } else {
        grpreg(X_model(), y_weighted, group = groups,
               penalty = "grLasso", intercept = FALSE, lambda = lambda)
      }
    },
    error = function(e) lm(y_weighted ~ 0 + X_weighted)
  )

  final_coefs <- as.vector(coef(final_model))
  if (length(final_coefs) > ncol(X_weighted)) final_coefs <- final_coefs[-1]
  final_coefs[is.na(final_coefs)] <- 0

  # Scale of each estimated pathway: a group shrunk to zero gives a flat kernel.
  idx1 <- 1:(spline_df + 1)
  idx2 <- c(1, (spline_df + 2):length(final_coefs))
  s1 <- sd(cbind(1, scale(B1_final)) %*% final_coefs[idx1])
  s2 <- sd(cbind(1, scale(B2_final)) %*% final_coefs[idx2])

  # Fix the sign ambiguity of the kernels.
  if (sum(diag(C1)) < 0) b1 <- -b1
  if (sum(diag(C2)) < 0) b2 <- -b2

  # ---------------------------------------------------------------------------
  # Phase 3: unpenalized natural spline for the shape of g(x1 -> y)
  # ---------------------------------------------------------------------------
  p1_matching <- index1(b1 * s1)
  g1_model <- lm(Y_pool ~ ns(p1_matching, df = spline_df))

  list(beta1 = b1 * s1, beta2 = b2 * s2, g1_model = g1_model)
}

#' Callback that saves one animation frame per loss improvement
#' (kernel heatmap and current link fit). Frames go to `frame_dir`.
make_frame_recorder <- function(dat, tm, d, frame_dir) {
  dir.create(frame_dir, showWarnings = FALSE, recursive = TRUE)
  unlink(file.path(frame_dir, "*.png"))
  dt <- tm[2] - tm[1]
  phi <- make_on_basis(d, tm)
  Y_pool <- c(as.vector(dat$yxo[, 1, ]), as.vector(dat$yxa[, 1, ]))
  palette <- colorRampPalette(c("dodgerblue", "white", "firebrick"))(100)
  frame <- 0

  function(C1, C2, label) {
    frame <<- frame + 1
    png(file.path(frame_dir, sprintf("frame_%04d.png", frame)),
        width = 1200, height = 500, res = 100)
    on.exit(dev.off())
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

    b1 <- phi %*% C1 %*% t(phi)
    b_max <- max(abs(b1)) + 1e-6
    image(tm, tm, b1, col = palette, zlim = c(-b_max, b_max),
          xlab = "t", ylab = "tau", main = paste("Kernel X1->Y |", label),
          cex.main = 1.5)
    contour(tm, tm, b1, add = TRUE, col = "black")

    p1 <- c(project(b1, dat$yxo[, 2, ], dt), project(b1, dat$yxa[, 2, ], dt))
    plot(p1, Y_pool, col = rgb(0.5, 0.5, 0.5, 0.2), pch = 16,
         xlab = "Projected Index", ylab = "Y", main = "Link Recovery X1->Y",
         cex.main = 1.5)
    fit <- tryCatch(lm(Y_pool ~ ns(p1, df = 5)), error = function(e) NULL)
    if (!is.null(fit)) {
      grid <- seq(min(p1), max(p1), length.out = 100)
      lines(grid, predict(fit, newdata = data.frame(p1 = grid)),
            col = "red", lwd = 4)
    }
  }
}
