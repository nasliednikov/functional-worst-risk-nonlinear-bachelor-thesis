# =============================================================================
#  Plotting functions for the figures of the thesis and the poster.
#  Each draw_* function draws on the current device; save_plot() writes it
#  as both PNG (high resolution) and PDF.
# =============================================================================

save_plot <- function(draw, path_no_ext, width = 10, height = 8, res = 500) {
  png(paste0(path_no_ext, ".png"), width = width * res, height = height * res, res = res)
  draw()
  dev.off()
  pdf(paste0(path_no_ext, ".pdf"), width = width, height = height)
  draw()
  dev.off()
  invisible(path_no_ext)
}

# Pairwise normalization: scale two kernels by their joint maximum.
normalize_pair <- function(b1, b2) {
  m <- max(abs(c(b1, b2)))
  list(b1 = b1 / m, b2 = b2 / m)
}

# Figure 3: one observational and one interventional sample curve.
draw_environments <- function(dat, tms) {
  function() {
    par(mfrow = c(1, 2))
    for (env in list(list(x = dat$yxo, title = "Observational Environment"),
                     list(x = dat$yxa, title = "Interventional Environment"))) {
      plot(tms, env$x[, 1, 1], type = "l", ylim = range(env$x[, , 1]),
           ylab = "Value", xlab = "Time (t)", main = env$title)
      lines(tms, env$x[, 2, 1], lty = 2, col = "blue")
      lines(tms, env$x[, 3, 1], lty = 3, col = "red")
      legend("bottomleft", legend = c("Y", "X1", "X2"), lty = 1:3,
             col = c("black", "blue", "red"), cex = 0.7)
    }
  }
}

# Figure 4: true, pooled (gamma = 0.5), regularized population and
# regularized small-sample kernels, pairwise normalized, X1 row and X2 row.
draw_kernel_heatmaps <- function(tms, true_b1, fits) {
  true_pair <- normalize_pair(true_b1, matrix(0, length(tms), length(tms)))
  pairs <- c(list(true_pair),
             lapply(fits, function(f) normalize_pair(f$beta1, f$beta2)))
  titles_x1 <- list(expression("True " * beta[x[1] * y]),
                    expression("Pooled population " * beta[x[1] * y]),
                    expression("Regularized population " * beta[x[1] * y]),
                    expression("Regularized estimate of " * beta[x[1] * y]))
  titles_x2 <- list(expression("True " * beta[x[2] * y]),
                    expression("Pooled population " * beta[x[2] * y]),
                    expression("Regularized population " * beta[x[2] * y]),
                    expression("Regularized estimate of " * beta[x[2] * y]))
  panel <- function(mat, title) {
    image(tms, tms, mat, zlim = c(-1, 1), axes = FALSE, main = title, cex.main = 1.0)
    box(col = "darkgray")
    axis(1, at = seq(0, 1, 0.2), cex.axis = 0.8)
    axis(2, at = seq(0, 1, 0.2), cex.axis = 0.8, las = 0)
  }
  function() {
    par(mfrow = c(2, 4), mar = c(3, 3, 3, 1), oma = c(4, 0, 0, 0))
    for (k in 1:4) panel(pairs[[k]]$b1, titles_x1[[k]])
    for (k in 1:4) panel(pairs[[k]]$b2, titles_x2[[k]])
    for (k in 1:4) {
      mtext(sprintf("(%s)", letters[k]), side = 1, line = 1, outer = TRUE,
            at = (2 * k - 1) / 8, cex = 1.2)
    }
  }
}

# Figures 5-8: true vs estimated kernels as 3D surfaces (X1 row, X2 row).
draw_kernels_3d <- function(tms, true_b1, est) {
  est_n <- normalize_pair(est$beta1, est$beta2)
  true_n <- list(b1 = true_b1 / max(abs(true_b1)),
                 b2 = matrix(0, length(tms), length(tms)))
  surface <- function(mat, col, title) {
    persp(tms, tms, mat, col = col, border = "gray70", zlim = c(-1, 1),
          xlab = "t", ylab = "tau", zlab = "Beta", main = title,
          theta = 35, phi = 25, shade = 0.15, ticktype = "detailed", cex.main = 1.3)
  }
  function() {
    par(mfrow = c(2, 2), mar = c(1, 2, 3, 1))
    surface(true_n$b1, "dodgerblue", expression("TRUE Kernel (" * X[1] %->% Y * ")"))
    surface(est_n$b1, "limegreen", expression("ESTIMATED Kernel (" * X[1] %->% Y * ")"))
    surface(true_n$b2, "dodgerblue", expression("TRUE Kernel (" * X[2] %->% Y * ")"))
    surface(est_n$b2, "limegreen", expression("ESTIMATED Kernel (" * X[2] %->% Y * ")"))
  }
}

# Figure 9: true link g(x1 -> y) vs the recovered natural spline.
draw_recovered_link <- function(true_index, actual_Y, curve_x, curve_y, g_true) {
  function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(true_index, actual_Y, col = rgb(0.5, 0.5, 0.5, 0.4), pch = 16,
         xlab = "True Projected Index of X1", ylab = "Observed Y",
         main = "Recovered X1 -> Y Non-Linear Relationship",
         cex.main = 1.5, cex.lab = 1.2)
    grid_x <- seq(min(true_index), max(true_index), length.out = 200)
    lines(grid_x, g_true(grid_x), col = "blue", lwd = 5, lty = 2)
    lines(curve_x, curve_y, col = "red", lwd = 5)
    legend("topleft", legend = c("Noisy Data", "Estimated g(x) Spline", "True g(x)"),
           col = c("gray", "red", "blue"), pch = c(16, NA, NA),
           lty = c(NA, 1, 2), lwd = c(NA, 5, 5), cex = 1.2, bg = "white")
  }
}

# Figure 10: true vs predicted Y on the shifted test set.
draw_prediction_tracking <- function(index, y_true, y_pred) {
  function() {
    plot(index, y_true, pch = 16, col = rgb(0, 0, 1, 0.3),
         xlab = "Estimated Projected Index of X1 (Unseen Shifted Test Set)",
         ylab = "Y Value", main = "True Y vs Predicted Y (Interventional Out-of-Sample)",
         xlim = range(index, na.rm = TRUE),
         ylim = range(c(y_true, y_pred), na.rm = TRUE))
    points(index, y_pred, pch = 16, col = rgb(1, 0, 0, 0.3))
    legend("topleft", legend = c("True Y (Shifted Env)", "Predicted Y"),
           col = c("blue", "red"), pch = 16)
  }
}

# Figure 11: 95% prediction interval of the recovered link on the test set.
draw_prediction_interval <- function(index, y_true, intervals) {
  ord <- order(index)
  function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(index, y_true, pch = 16, col = rgb(0, 0, 1, 0.2),
         xlab = "Estimated Projected Index of X1", ylab = "Y Value",
         main = "Out-of-Sample Prediction Interval (95%)",
         cex.main = 1.5, cex.lab = 1.2)
    lines(index[ord], intervals[ord, "lwr"], col = "red", lty = 2, lwd = 3)
    lines(index[ord], intervals[ord, "upr"], col = "red", lty = 2, lwd = 3)
    lines(index[ord], intervals[ord, "fit"], col = "red", lwd = 5)
    legend("topleft", legend = c("True Y (Test Set)", "Predicted Trend", "95% Prediction Bounds"),
           col = c("blue", "red", "red"), lty = c(NA, 1, 2), pch = c(16, NA, NA),
           lwd = c(NA, 5, 3), cex = 1.2, bg = "white")
  }
}

# Figure 12: actual vs predicted on the test set.
draw_actual_vs_predicted <- function(y_pred, y_true) {
  function() {
    par(mar = c(4.5, 4.5, 3, 1))
    plot(y_pred, y_true, pch = 16, col = rgb(0, 0, 1, 0.3),
         xlab = "Predicted Y (Model Expectation)", ylab = "Actual Y (Test Set Data)",
         main = "Actual vs. Predicted")
    abline(a = 0, b = 1, col = "black", lwd = 3, lty = 2)
    legend("topleft", legend = c("Test Samples", "Ideal Fit (Y = X)"),
           col = c("blue", "black"), pch = c(16, NA), lty = c(NA, 2),
           lwd = c(NA, 3), cex = 0.8)
  }
}

# Figure 13: residuals on the test set.
draw_residuals <- function(y_pred, residuals) {
  function() {
    par(mar = c(4.5, 4.5, 3, 1))
    plot(y_pred, residuals, pch = 16, col = rgb(0, 0, 1, 0.3),
         xlab = "Predicted Y (Model Expectation)", ylab = "Residuals (Error)",
         main = "Residuals Analysis")
    abline(h = 0, col = "black", lwd = 3, lty = 2)
    legend("topleft", legend = c("Prediction Errors", "Zero Error Line"),
           col = c("blue", "black"), pch = c(16, NA), lty = c(NA, 2),
           lwd = c(NA, 3), cex = 0.8)
  }
}

# Presentation figure: both environments along the true X1 index, with densities.
draw_environment_scatter <- function(index_obs, y_obs, index_int, y_int, g_true) {
  function() {
    layout(matrix(c(1, 2), nrow = 2), heights = c(1, 3))
    par(mar = c(0, 5, 4, 2) + 0.1)
    d_obs <- density(index_obs)
    d_int <- density(index_int)
    plot(d_obs, xlim = range(c(index_obs, index_int)),
         ylim = c(0, max(c(d_obs$y, d_int$y)) * 1.1),
         main = "Observational vs Interventional Environments", xlab = "",
         ylab = "Density", col = "blue", lwd = 3, axes = FALSE, cex.main = 1.8)
    polygon(d_obs, col = rgb(0, 0.4, 1, 0.4), border = NA)
    lines(d_int, col = "red", lwd = 3)
    polygon(d_int, col = rgb(1, 0.2, 0.2, 0.4), border = NA)
    legend("topleft", legend = c("Obs. Density", "Int. Density"),
           fill = c(rgb(0, 0.4, 1, 0.6), rgb(1, 0.2, 0.2, 0.6)),
           border = NA, bty = "n", cex = 1.4)

    par(mar = c(5, 5, 0, 2) + 0.1)
    plot(index_int, y_int, col = rgb(1, 0.2, 0.2, 0.3), pch = 16,
         xlab = expression("True Projected Index (" * beta^T * X[1] * ")"),
         ylab = "Target Y", ylim = range(c(y_obs, y_int)),
         xlim = range(c(index_obs, index_int)), cex.lab = 1.4, cex.axis = 1.2)
    points(index_obs, y_obs, col = rgb(0, 0.4, 1, 0.4), pch = 16)
    grid_x <- seq(min(c(index_obs, index_int)), max(c(index_obs, index_int)), length.out = 200)
    lines(grid_x, g_true(grid_x), col = "black", lwd = 6, lty = 2)
    legend("topleft",
           legend = c("Observational Env.", "Interventional Env.", "Invariant True g(x)"),
           col = c(rgb(0, 0.4, 1), rgb(1, 0.2, 0.2), "black"), pch = c(16, 16, NA),
           lty = c(NA, NA, 2), lwd = c(NA, NA, 6), cex = 1.4, bg = "white")
  }
}

# Presentation figure: 2x2 heatmaps of true and estimated kernels.
draw_kernel_heatmaps_2x2 <- function(tms, true_b1, est) {
  est_n <- normalize_pair(est$beta1, est$beta2)
  true_n <- list(b1 = true_b1 / max(abs(true_b1)),
                 b2 = matrix(0, length(tms), length(tms)))
  colors <- colorRampPalette(c("dodgerblue", "white", "firebrick"))(100)
  panel <- function(mat, title) {
    image(tms, tms, mat, col = colors, zlim = c(-1, 1),
          xlab = "t", ylab = "tau", main = title, cex.main = 1.5)
    contour(tms, tms, mat, add = TRUE, col = "black")
  }
  function() {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))
    panel(true_n$b1, expression("TRUE Kernel (" * X[1] %->% Y * ")"))
    panel(est_n$b1, expression("ESTIMATED Kernel (" * X[1] %->% Y * ")"))
    panel(true_n$b2, expression("TRUE Kernel (" * X[2] %->% Y * ")"))
    panel(est_n$b2, expression("ESTIMATED Kernel (" * X[2] %->% Y * ")"))
    mtext("Causal Kernels 2D Heatmaps", side = 3, outer = TRUE, cex = 2, font = 2)
  }
}
