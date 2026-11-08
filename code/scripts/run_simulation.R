# =============================================================================
#  Reproduces the simulation study of the thesis (Section 6) for one scenario.
#
#  Usage (from the `code/` folder):
#    Rscript scripts/run_simulation.R sinusoidal
#    Rscript scripts/run_simulation.R cubic
#    Rscript scripts/run_simulation.R sinusoidal --gif   # also records the
#                                                        # optimization animation
#
#  Output: results/<scenario>/  (figures as PNG + PDF, and metrics.txt)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else "sinusoidal"
make_gif <- "--gif" %in% args

source("R/simulate.R")
source("R/estimate.R")
source("R/plots.R")

out_dir <- file.path("results", scenario)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(name) file.path(out_dir, name)

links <- link_functions(scenario)
tms <- seq(0, 1, length = 50)
dt <- tms[2] - tms[1]

# Simulation settings used in the thesis
D_BASIS <- 5        # basis functions for the kernels
MU_SHIFT <- 3.0     # mean of the interventional shift A
SD_SHIFT <- 1.0     # sd of the interventional shift A
GAMMA <- 0.8        # worst-risk trade-off of the main fit
LAMBDA <- 0.1       # Group Lasso penalty of the main fit

# -----------------------------------------------------------------------------
# 1. Main data set and fit
#    Note: all random draws below come from one stream started by set.seed(42).
#    Their order (data, fits, test set) fixes the numbers reported in the thesis.
# -----------------------------------------------------------------------------
set.seed(42)
cat(sprintf("Scenario: %s\nGenerating data...\n", scenario))
dat <- simulate_fsem(n_samples = 1000, tm = tms, links = links, d = D_BASIS,
                     mu_A = MU_SHIFT, s_A = SD_SHIFT)

recorder <- if (make_gif) make_frame_recorder(dat, tms, D_BASIS, out("gif_frames")) else NULL

cat("Fitting (coordinate descent with BFGS + Group Lasso)...\n")
fit <- fit_fwrm_nonlinear(dat, tms, gamma = GAMMA, d = D_BASIS, spline_df = 5,
                          max_iter = 5, lambda = LAMBDA, on_improve = recorder)

if (make_gif) {
  if (requireNamespace("gifski", quietly = TRUE)) {
    frames <- list.files(out("gif_frames"), pattern = "\\.png$", full.names = TRUE)
    gifski::gifski(frames, gif_file = out("algorithm_steps.gif"),
                   delay = 0.1, width = 1200, height = 500)
  } else {
    warning("Package 'gifski' is not installed: frames saved, GIF not created.")
  }
}

# -----------------------------------------------------------------------------
# 2. Figures from the main fit
# -----------------------------------------------------------------------------
X1_obs <- dat$yxo[, 2, ]
X1_int <- dat$yxa[, 2, ]
Y_obs <- as.vector(dat$yxo[, 1, ])
Y_int <- as.vector(dat$yxa[, 1, ])

true_index <- project(dat$beta_x1y, X1_obs, dt)
estimated_index <- project(fit$beta1, X1_obs, dt)

# Recovered link: evaluate the spline on the estimated index range and map
# that range onto the true index scale, so both curves share one x-axis.
grid_est <- seq(min(estimated_index), max(estimated_index), length.out = 200)
curve_y <- predict(fit$g1_model, newdata = data.frame(p1_matching = grid_est))
to_true_scale <- lm(true_index ~ estimated_index)
curve_x <- predict(to_true_scale, newdata = data.frame(estimated_index = grid_est))

save_plot(draw_environments(dat, tms), out("fig03_environments"), width = 10, height = 5)
save_plot(draw_kernels_3d(tms, dat$beta_x1y, fit), out("fig05-08_kernels_3d"), width = 10, height = 9)
save_plot(draw_recovered_link(true_index, Y_obs, curve_x, curve_y, links$x1_to_y),
          out("fig09_recovered_link"))
save_plot(draw_environment_scatter(true_index, Y_obs, project(dat$beta_x1y, X1_int, dt),
                                   Y_int, links$x1_to_y),
          out("extra_environment_scatter"))
save_plot(draw_kernel_heatmaps_2x2(tms, dat$beta_x1y, fit), out("extra_kernel_heatmaps_2x2"),
          width = 7, height = 10)

# -----------------------------------------------------------------------------
# 3. Figure 4: pooled vs regularized estimation, large vs small sample
# -----------------------------------------------------------------------------
cat("Estimating the models for Figure 4 (population n = 400, estimate n = 100)...\n")
set.seed(42)
dat_pop <- simulate_fsem(n_samples = 400, tm = tms, links = links, d = D_BASIS,
                         mu_A = MU_SHIFT, s_A = SD_SHIFT)
dat_emp <- simulate_fsem(n_samples = 100, tm = tms, links = links, d = D_BASIS,
                         mu_A = MU_SHIFT, s_A = SD_SHIFT)
invisible(capture.output({
  fit_pooled_pop <- fit_fwrm_nonlinear(dat_pop, tms, gamma = 0.5, d = D_BASIS, lambda = 0)
  fit_reg_pop <- fit_fwrm_nonlinear(dat_pop, tms, gamma = 0.9, d = D_BASIS, lambda = "cv")
  fit_reg_emp <- fit_fwrm_nonlinear(dat_emp, tms, gamma = 0.9, d = D_BASIS, lambda = "cv")
}))
save_plot(draw_kernel_heatmaps(tms, dat_pop$beta_x1y,
                               list(fit_pooled_pop, fit_reg_pop, fit_reg_emp)),
          out("fig04_kernel_heatmaps"), width = 12, height = 6.5)

# -----------------------------------------------------------------------------
# 4. Out-of-sample evaluation on a new shifted test set
# -----------------------------------------------------------------------------
cat("Evaluating on an out-of-sample test set (n = 200)...\n")
dat_test <- simulate_fsem(n_samples = 200, tm = tms, links = links, d = D_BASIS,
                          mu_A = MU_SHIFT, s_A = SD_SHIFT)
Y_obs_test <- as.vector(dat_test$yxo[, 1, ])
Y_int_test <- as.vector(dat_test$yxa[, 1, ])
p1_obs_test <- project(fit$beta1, dat_test$yxo[, 2, ], dt)
p1_int_test <- project(fit$beta1, dat_test$yxa[, 2, ], dt)

pred_obs_test <- predict(fit$g1_model, newdata = data.frame(p1_matching = p1_obs_test))
intervals <- predict(fit$g1_model, newdata = data.frame(p1_matching = p1_int_test),
                     interval = "prediction", level = 0.95)
pred_int_test <- intervals[, "fit"]

mse_obs <- mean((Y_obs_test - pred_obs_test)^2)
mse_int <- mean((Y_int_test - pred_int_test)^2)

metrics <- c(
  sprintf("Scenario: %s", scenario),
  sprintf("L2 distance, true vs estimated beta_1:     %.4f",
          sqrt(sum((as.vector(dat$beta_x1y) - as.vector(fit$beta1))^2))),
  sprintf("Correlation, true vs estimated X1 index:  %.4f", cor(true_index, estimated_index)),
  sprintf("Out-of-sample MSE, observational:         %.4f", mse_obs),
  sprintf("Out-of-sample MSE, interventional:        %.4f", mse_int),
  sprintf("Invariance gap |MSE_obs - MSE_int|:       %.4f", abs(mse_obs - mse_int))
)
writeLines(metrics, out("metrics.txt"))
cat(metrics, sep = "\n")

save_plot(draw_prediction_tracking(p1_int_test, Y_int_test, pred_int_test),
          out("fig10_prediction_tracking"))
save_plot(draw_prediction_interval(p1_int_test, Y_int_test, intervals),
          out("fig11_prediction_interval"))
save_plot(draw_actual_vs_predicted(pred_int_test, Y_int_test),
          out("fig12_actual_vs_predicted"))
save_plot(draw_residuals(pred_int_test, Y_int_test - pred_int_test),
          out("fig13_residuals"))

cat(sprintf("Done. Results in %s/\n", out_dir))
