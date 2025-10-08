# =========================================================
# Simulation under H0 (Type I error estimation)
# =========================================================

library(functintermethsnp)
library(mgcv)

# Load utility functions
source("scripts/functions_sim.R")

# =========================================================
# 1. Choose sample size
# =========================================================
sample_size <- 200   # <-- change to 400 for large-sample scenario

if (sample_size == 400 && file.exists("data/mat_est_meth_n400.rds")) {
  mat_methyl <- readRDS("data/mat_est_meth_n400.rds")
  message("Loaded: mat_est_meth_n400 (", nrow(mat_methyl), " individuals)")
} else if (sample_size == 200 && file.exists("data/mat_est_meth_n200.rds")) {
  mat_methyl <- readRDS("data/mat_est_meth_n200.rds")
  message("Loaded: mat_est_meth_n200 (", nrow(mat_methyl), " individuals)")
} else {
  stop("Please set sample_size = 200 or 400 and ensure the corresponding file exists in 'data/'.")
}

N <- nrow(mat_methyl)

# =========================================================
# 2. Model parameters
# =========================================================
nb_snp <- 20
nb_fct <- 10
rho_true <- 0.1
rho_fit  <- 0.1
kernel_name <- "convex"
region_range <- c(11190000, 11460000)
RGN <- seq(0, 1, length.out = 1000)
MC_H0 <- 5000

# =========================================================
# 3. Design components (same across repetitions)
# =========================================================
set.seed(1234)
W <- matrix(rnorm(N, 0, 0.1), ncol = 1)
colnames(W) <- "W1"

frq_all <- runif(nb_snp, 0.05, 0.2)
G <- sapply(frq_all, function(p) rbinom(N, 2, p))
colnames(G) <- paste0("snp", 1:nb_snp)

# Keep SNP order identical to G columns (no sort)
pos_snp <- sample(seq(region_range[1], region_range[2]), nb_snp)

# Normalize SNP positions into [0, 1]
pos_snp_norm <- (pos_snp - region_range[1]) / (region_range[2] - region_range[1])

# Compute integral of beta(t) * methylation(t)
Z <- intg_beta_methyl(id = 1:N, mat_methyl = mat_methyl)

# Compute matrix of integrals between kernel and methylation curves
Omega_rho_true <- mat_intg_kernel_methyl(
  pos_snp = pos_snp_norm,
  id = 1:N,
  rho = rho_true,
  mat_methyl = mat_methyl,
  kernel_name = "convex"
)

# Compute matrix of interaction terms
K_rho_true <- Omega_rho_true * G


zeta0 <- 0
zeta  <- 0.3
alpha <- c(1.2, 1.0, 0.8, 0.5, 0.3, 0.9, 1.6, 1.3, 0.67, 0.89,
           1.45, 1.4, 0.45, 0.7, 1.35, 0.95, 0.55, 0.88, 1.3, 0.5)

# =========================================================
# 4. Monte Carlo simulations under H0
# =========================================================
pval_H0 <- numeric(MC_H0)
message("\n===== Running ", MC_H0, " simulations under H0 =====")

for (i in 1:MC_H0) {
  set.seed(1000 + i)
  eta_H0 <- rep(0, nb_snp)
  Y_H0 <- round(simulate_y_h1(N, zeta0, zeta, alpha, eta_H0, W, G, Z, K_rho_true), 5)
  
  fit <- try(
    fitMethSNPfunctInter_continu(
      Y = Y_H0, W = W, G = G, mat_methyl = mat_methyl,
      pos_snp = pos_snp, rho = rho_fit, nb_fct = nb_fct,
      kernel_name = kernel_name, region_range = region_range, RGN = RGN
    ),
    silent = TRUE
  )
  
  pval_H0[i] <- if (inherits(fit, "try-error")) NA else fit$pvalue
  if (i %% 500 == 0) message("... ", i, " / ", MC_H0)
}

# =========================================================
# 5. Summary results
# =========================================================
alpha_level <- 0.05
type_I_error <- mean(pval_H0 < alpha_level, na.rm = TRUE)

cat("\n========================================================\n")
cat("Sample size:", N, "\n")
cat("Type I error (H0):", round(type_I_error, 4), "\n")
cat("========================================================\n")

# =========================================================
# 6. Save results
# =========================================================
results_H0 <- list(
  N = N,
  nb_snp = nb_snp,
  rho_true = rho_true,
  rho_fit = rho_fit,
  kernel = kernel_name,
  nb_fct = nb_fct,
  alpha_level = alpha_level,
  type_I_error = type_I_error,
  pval_H0 = pval_H0
)

dir.create("results", showWarnings = FALSE)
saveRDS(results_H0, file = sprintf("results/sim_H0_N%d_%s.rds", N, format(Sys.Date(), "%Y%m%d")))
message("Results saved to: results/sim_H0_N", N, "_", format(Sys.Date(), "%Y%m%d"), ".rds")
