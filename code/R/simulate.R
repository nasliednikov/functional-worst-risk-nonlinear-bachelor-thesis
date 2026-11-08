# =============================================================================
#  Data generation for the non-linear functional structural equation model
#
#  DAG:  X1 -> Y -> X2
#    X1 is the invariant cause of Y; X2 is a spurious variable (a consequence
#    of Y). Curves live on t in [0, 1] and are built from an orthonormal sine
#    basis with d functions.
#
#  Two environments are simulated:
#    - observational:  exogenous noise only
#    - interventional: an additive shift A ~ N(mu_A, s_A^2) is injected into
#                      X1 and X2, while the mechanism that generates Y is
#                      left unchanged.
# =============================================================================

#' Orthonormal sine basis on the grid `tms`: phi_k(t) = sqrt(2) sin(k pi t).
make_on_basis <- function(d, tms) {
  sapply(seq_len(d), function(k) sqrt(2) * sin(k * pi * tms))
}

#' True non-linear link functions for the two scenarios used in the thesis.
link_functions <- function(scenario = c("sinusoidal", "cubic")) {
  scenario <- match.arg(scenario)
  switch(scenario,
    sinusoidal = list(
      x1_to_y = function(x) x + sin(2 * x),
      y_to_x2 = function(x) -x + sin(x)
    ),
    cubic = list(
      x1_to_y = function(x) x + 0.1 * x^3,
      y_to_x2 = function(x) -x + 0.1 * x^3
    )
  )
}

#' Simulate `n_samples` curves per environment.
#'
#' @param links  list with functions `x1_to_y` and `y_to_x2` (see link_functions()).
#' @return list with arrays `yxo` (observational) and `yxa` (interventional),
#'   each of dimension [time, variable, sample] with variables (Y, X1, X2),
#'   and the true kernels `beta_x1y`, `beta_yx2`.
simulate_fsem <- function(n_samples, tm, links, d = 5, mu_A = 0, s_A = 1) {
  m <- length(tm)
  dt <- tm[2] - tm[1]
  phi <- make_on_basis(d, tm)

  beta_x1y <- phi %*% diag(d) %*% t(phi)
  beta_yx2 <- phi %*% diag(d) %*% t(phi)

  yxo <- array(dim = c(m, 3, n_samples))
  yxa <- array(dim = c(m, 3, n_samples))

  for (i in seq_len(n_samples)) {
    # Exogenous noise curves for (Y, X1, X2) in each environment.
    # The order of the random draws below fixes the reproducible data set.
    e0_curves <- phi %*% matrix(rnorm(d * 3), ncol = 3)
    e1_curves <- phi %*% matrix(rnorm(d * 3), ncol = 3)

    # Shift A acts on X1 and X2 only; Y is left untouched.
    A_mat <- matrix(0, nrow = d, ncol = 3)
    A_mat[, 2] <- rnorm(d, mean = mu_A, sd = s_A)
    A_mat[, 3] <- rnorm(d, mean = mu_A, sd = s_A)
    A_curves <- phi %*% A_mat

    # Observational environment
    X1_o <- e0_curves[, 2]
    Y_o <- links$x1_to_y(t(beta_x1y) %*% X1_o * dt) + e0_curves[, 1]
    X2_o <- links$y_to_x2(t(beta_yx2) %*% Y_o * dt) + e0_curves[, 3]
    yxo[, , i] <- cbind(Y_o, X1_o, X2_o)

    # Interventional environment (same mechanisms, shifted inputs)
    X1_a <- e1_curves[, 2] + A_curves[, 2]
    Y_a <- links$x1_to_y(t(beta_x1y) %*% X1_a * dt) + e1_curves[, 1]
    X2_a <- links$y_to_x2(t(beta_yx2) %*% Y_a * dt) + e1_curves[, 3] + A_curves[, 3]
    yxa[, , i] <- cbind(Y_a, X1_a, X2_a)
  }

  list(yxo = yxo, yxa = yxa, beta_x1y = beta_x1y, beta_yx2 = beta_yx2)
}

#' Project curves X (time x samples) with kernel beta: <beta, X>(t).
project <- function(beta, X, dt) {
  as.vector(t(beta) %*% X * dt)
}
