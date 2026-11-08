# =============================================================================
#  Illustration of the shift budget sets C_A^gamma(A) (poster figure).
#
#  C_A^gamma(A) = { A' : gamma * K_A - K_A' is positive semi-definite }
#  In 2D this means the covariance ellipse of a new shift A' must lie inside
#  the ellipse of gamma * K_A.
#
#  Usage (from the `code/` folder):  Rscript scripts/shift_budget_regions.R
#  Output: results/shift_budget_regions.png and .pdf
# =============================================================================

source("R/plots.R")  # for save_plot()

#' Coordinates of the `level` confidence ellipse of a 2D covariance matrix.
ellipse_points <- function(cov_matrix, center = c(0, 0), level = 0.95, npoints = 200) {
  radius <- sqrt(qchisq(level, df = 2))
  eig <- eigen(cov_matrix)
  theta <- seq(0, 2 * pi, length.out = npoints)
  circle <- cbind(cos(theta), sin(theta))
  t(center + radius * eig$vectors %*% diag(sqrt(eig$values)) %*% t(circle))
}

K_A <- matrix(c(1.5, 0.8,
                0.8, 1.2), nrow = 2)                 # observed shift covariance
K_valid <- matrix(c(0.6, 0.1, 0.1, 0.4), nrow = 2)   # A' inside C^1(A)
K_invalid <- matrix(c(1.9, -0.5, -0.5, 0.8), nrow = 2) # A' outside C^1(A)

gammas <- c(2.0, 1.0, 0.5)
fills <- c(rgb(0.7, 0.8, 0.9, 0.6), rgb(0.4, 0.5, 0.8, 0.6), rgb(0.1, 0.2, 0.6, 0.7))

draw_shift_budget <- function() {
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(NA, xlim = c(-4.5, 4.5), ylim = c(-4.5, 4.5),
       xlab = expression("Shift Magnitude on " * X[1]),
       ylab = expression("Shift Magnitude on " * X[2]),
       main = expression("Distributional Shift Budget Areas " * C[A]^gamma * "(A)"),
       cex.main = 1.4, cex.lab = 1.2)
  grid(col = "lightgray", lty = "dotted", lwd = 1)

  for (i in seq_along(gammas)) {
    ell <- ellipse_points(gammas[i] * K_A)
    polygon(ell[, 1], ell[, 2], col = fills[i], border = "darkblue", lwd = 2)
  }
  ell <- ellipse_points(K_valid)
  lines(ell[, 1], ell[, 2], col = "green3", lwd = 3, lty = 2)
  ell <- ellipse_points(K_invalid)
  lines(ell[, 1], ell[, 2], col = "red2", lwd = 3, lty = 2)

  legend("topleft",
         legend = c(expression(gamma == 2.0 * " (Expanded Budget)"),
                    expression(gamma == 1.0 * " (Baseline Shift " * K[A] * ")"),
                    expression(gamma == 0.5 * " (Restricted Budget)"),
                    expression("Valid " * A * "'" %in% C[A]^"1.0" * "(A)"),
                    expression("Invalid " * A * "'" %notin% C[A]^"1.0" * "(A)")),
         fill = c(fills, NA, NA), border = c(rep("darkblue", 3), NA, NA),
         col = c(NA, NA, NA, "green3", "red2"), lty = c(NA, NA, NA, 2, 2),
         lwd = c(NA, NA, NA, 3, 3), cex = 1.0, bg = "white")
  text(0, -4.0, expression("Condition: " * gamma * K[A] - K["A'"] * " is positive semi-definite"),
       cex = 1.1, font = 3)
}

dir.create("results", showWarnings = FALSE)
save_plot(draw_shift_budget, "results/shift_budget_regions", width = 10, height = 8.6)
