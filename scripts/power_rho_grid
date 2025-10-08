# =========================================================
# Power analysis across (rho, eta) grid
# Requires: data/mat_est_meth_n200.rds or _n400.rds
# =========================================================

library(functintermethsnp)
library(mgcv)

# Utils (true_beta, intg_beta_methyl, kernels, mat_intg_kernel_methyl, simulate_y_h1)
source("scripts/functions_sim.R")

# =========================================================
# 1) Load methylation curves
# =========================================================
sample_size <- 400  # <-- mettre 200 si besoin

if (sample_size == 400 && file.exists("data/mat_est_meth_n400.rds")) {
  mat_methyl <- readRDS("data/mat_est_meth_n400.rds")
  message("Loaded: mat_est_meth_n400 (", nrow(mat_methyl), " individuals)")
} else if (sample_size == 200 && file.exists("data/mat_est_meth_n200.rds")) {
  mat_methyl <- readRDS("data/mat_est_meth_n200.rds")
  message("Loaded: mat_est_meth_n200 (", nrow(mat_methyl), " individuals)")
} else {
  stop("Please set sample_size = 200 or 400 and ensure the RDS file exists in data/.")
}
N <- nrow(mat_methyl)

# =========================================================
# 2) Parameters
# =========================================================
nb_snp      <- 20
nb_fct      <- 10
kernel_name <- "convex"
region_range <- c(11190000, 11460000)
RGN <- seq(0, 1, length.out = 1000)
MC  <- 1000

# Grid over rho (fit) and eta (truth)
val_rho    <- c(0.1, 0.2, 0.4, 0.9, 1, 1.2, 2, 8, 8.5, 9)
eta_values <- seq(1, 40, by = 2)
nb_rho <- length(val_rho)
mmc    <- length(eta_values)

# =========================================================
# 3) Design (constant across reps)
# =========================================================
set.seed(1234)
W <- matrix(rnorm(N, 0, 0.1), ncol = 1); colnames(W) <- "W1"

frq_all <- runif(nb_snp, 0.05, 0.2)
G <- sapply(frq_all, function(p) rbinom(N, 2, p))
colnames(G) <- paste0("snp", 1:nb_snp)

# Positions SNP (ordre = colonnes de G), puis normalisation [0,1]
pos_snp <- sample(seq(region_range[1], region_range[2]), nb_snp)
pos_snp_norm <- (pos_snp - region_range[1]) / (region_range[2] - region_range[1])

# Intégrale ∫ beta(t)*Pi_i(t) dt
Z <- intg_beta_methyl(id = 1:N, mat_methyl = mat_methyl)

# K_rho utilisé pour la génération (rho_true fixe)
rho_true <- 0.1
Omega_rho_true <- mat_intg_kernel_methyl(
  pos_snp = pos_snp_norm, id = 1:N, rho = rho_true,
  mat_methyl = mat_methyl, kernel_name = kernel_name
)
K_rho_true <- Omega_rho_true * G

# Effets fixes
zeta0 <- 0
zeta  <- 0.3
alpha <- c(1.2,1,0.8,0.5,0.3,0.9,1.6,1.3,0.67,0.89,
           1.45,1.4,0.45,0.7,1.35,0.95,0.55,0.88,1.3,0.5)

# =========================================================
# 4) MC loops: for each rho_fit, sweep eta_values
# =========================================================
pval_grid  <- array(NA_real_, dim = c(MC, nb_rho, mmc))  # [rep, rho, eta]
power_grid <- matrix(NA_real_, nrow = nb_rho, ncol = mmc,
                     dimnames = list(paste0("rho=", val_rho),
                                     paste0("eta=", eta_values)))

message("\n===== Power over (rho, eta) grid =====")
for (kk in 1:nb_rho) {
  rho_fit <- val_rho[kk]
  message("\n--- Fitting with rho = ", rho_fit, " ---")
  
  for (k in 1:mmc) {
    eta_vec <- rep(eta_values[k], nb_snp)  # η varie (comme dans power_eta_grid.R)
    for (m in 1:MC) {
      set.seed(3000 + m + 10000*k + 100000*kk)
      Y <- round(simulate_y_h1(N, zeta0, zeta, alpha, eta_vec, W, G, Z, K_rho_true), 5)
      
      fit <- try(
        fitMethSNPfunctInter_continu(
          Y = Y, W = W, G = G, mat_methyl = mat_methyl,
          pos_snp = pos_snp, rho = rho_fit, nb_fct = nb_fct,
          kernel_name = kernel_name, region_range = region_range, RGN = RGN
        ),
        silent = TRUE
      )
      
      pval_grid[m, kk, k] <- if (inherits(fit, "try-error")) NA_real_ else fit$pvalue
    }
    power_grid[kk, k] <- mean(pval_grid[, kk, k] < 0.05, na.rm = TRUE)
    message(sprintf("  -> Power(rho=%.3g, eta=%g) = %.4f",
                    rho_fit, eta_values[k], power_grid[kk, k]))
  }
}

# =========================================================
# 5) Save
# =========================================================
dir.create("results", showWarnings = FALSE)
res <- list(
  N = N,
  nb_snp = nb_snp,
  region_range = region_range,
  RGN = RGN,
  val_rho = val_rho,
  eta_values = eta_values,
  power_grid = power_grid,
  pval_grid = pval_grid
)

outfile <- sprintf("results/power_rho_eta_grid_N%d_%s.rds",
                   N, format(Sys.Date(), "%Y%m%d"))
saveRDS(res, file = outfile)
message("Saved: ", outfile)
