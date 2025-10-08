# =========================================================
# Power analysis varying eta values
# Requires: mat_est_meth200.rds or mat_est_meth400.rds in /data
# =========================================================

library(functintermethsnp)
library(mgcv)

# Load simulation utility functions
source("scripts/functions_sim.R")

# =========================================================
# 1. Load methylation data
# =========================================================
sample_size <- 400   # <-- change to 200 if needed

if (sample_size == 400 && file.exists("data/mat_est_meth_n400.rds")) {
  mat_methyl <- readRDS("data/mat_est_meth_n400.rds")
  message("Loaded: mat_est_meth_n400 (", nrow(mat_methyl), " individuals)")
} else if (sample_size == 200 && file.exists("data/mat_est_meth_n200.rds")) {
  mat_methyl <- readRDS("data/mat_est_meth_n200.rds")
  message("Loaded: mat_est_meth_n200 (", nrow(mat_methyl), " individuals)")
}


N <- nrow(mat_methyl)

# =========================================================
# 2. Simulation parameters
# =========================================================
nb_snp <- 20
nb_fct <- 10
rho_true <- 0.1
rho_fit  <- 0.1
kernel_name <- "convex"
region_range <- c(11190000, 11460000)
RGN <- seq(0, 1, length.out = 1000)
MC <- 1000

# Values of η to test
eta_values <- seq(1, 40, by = 2)
mmc <- length(eta_values)

# =========================================================
# 3. Design components
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
alpha <- c(1.2,1,0.8,0.5,0.3,0.9,1.6,1.3,0.67,0.89,
           1.45,1.4,0.45,0.7,1.35,0.95,0.55,0.88,1.3,0.5)

# =========================================================
# 4. Monte Carlo simulation loop over η values
# =========================================================
pval_eta <- matrix(NA, nrow = MC, ncol = mmc)
power_eta <- numeric(mmc)

message("\n===== Running power simulation across eta values =====")

for (k in 1:mmc) {
  eta_val <- rep(eta_values[k], nb_snp)
  message("\n--- η = ", eta_values[k], " ---")
  
  for (m in 1:MC) {
    set.seed(1000 + m)
    Y <- round(simulate_y_h1(N, zeta0, zeta, alpha, eta_val, W, G, Z, K_rho_true), 5)
    
    fit <- try(
      fitMethSNPfunctInter_continu(
        Y = Y, W = W, G = G, mat_methyl = mat_methyl,
        pos_snp = pos_snp, rho = rho_fit, nb_fct = nb_fct,
        kernel_name = kernel_name, region_range = region_range, RGN = RGN
      ),
      silent = TRUE
    )
    
    pval_eta[m, k] <- if (inherits(fit, "try-error")) NA else fit$pvalue
  }
  
  power_eta[k] <- mean(pval_eta[, k] < 0.05, na.rm = TRUE)
  message("Power(η=", eta_values[k], ") = ",power_eta[k])
}

# =========================================================
# 5. Save results
# =========================================================
dir.create("results", showWarnings = FALSE)
res <- list(
  N = N,
  nb_snp = nb_snp,
  rho_true = rho_true,
  rho_fit = rho_fit,
  eta_values = eta_values,
  power_eta = power_eta,
  pval_eta = pval_eta
)

saveRDS(res, file = sprintf("results/power_eta_grid_N%d_%s.rds", N, format(Sys.Date(), "%Y%m%d")))
message("Results saved to: results/power_eta_grid_N", N, "_", format(Sys.Date(), "%Y%m%d"), ".rds")
