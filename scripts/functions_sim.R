# =========================================================
# Utility functions for simulation scripts
# Author: Yvelin Gansou
# =========================================================

# ---------- True functional coefficient ----------
true_beta <- function(t) cos(3 * pi * t)

# ---------- Kernel definitions ----------
kernel_convex <- function(t, pos_snp, rho) exp(-rho * abs(t - pos_snp))
kernel_gaussian <- function(t, pos_snp, rho) exp(-(rho^2) * (t - pos_snp)^2)
Kernel_linear <- function(t, pos_snp, rho) pmax(1 - rho * abs(t - pos_snp), 0)

# ---------- Simulate continuous response ----------
simulate_y_h1 <- function(N, zeta_0, zeta, alpha, eta, W, G, Z, K_rho) {
  sig <- sqrt((
    var(
      zeta_0 +
        rowSums((rep(zeta, each = N)) * W) +
        rowSums((rep(alpha, each = N)) * G) +
        Z +
        rowSums((rep(eta, each = N)) * K_rho)
    ) / 10
  ))
  Y <- zeta_0 +
    rowSums((rep(zeta, each = N)) * W) +
    rowSums((rep(alpha, each = N)) * G) +
    Z +
    rowSums((rep(eta, each = N)) * K_rho) +
    rnorm(N, mean = 0, sd = sig)
  return(Y)
}

# =========================================================
# Approximations of integrals
# =========================================================

# ---- Integral of beta(t) * methylation(t) ----
intg_beta_methyl <- function(id, mat_methyl, RGN = seq(0, 1, length.out = 1000)) {
  intg_beta_methyl_elm <- function(mat_methyl, id, RGN) {
    (sum(true_beta(t = RGN) * mat_methyl[id, ])) / length(RGN)
  }
  sapply(id, intg_beta_methyl_elm, mat_methyl = mat_methyl, RGN = RGN)
}

# ---- Integral of kernel_convex * methylation ----
intg_kernelvex_methyl <- function(id, pos_snp, rho, mat_methyl, RGN = seq(0, 1, length.out = 1000)) {
  intg_kernelvex_methyl_elm <- function(pos_snp, rho, id, mat_methyl, RGN) {
    (sum(kernel_convex(t = RGN, pos_snp = pos_snp, rho = rho) * mat_methyl[id, ])) / length(RGN)
  }
  sapply(id, intg_kernelvex_methyl_elm, pos_snp = pos_snp, rho = rho, mat_methyl = mat_methyl, RGN = RGN)
}

# ---- Integral of kernel_gaussian * methylation ----
intg_kernelgaus_methyl <- function(id, pos_snp, rho, mat_methyl, RGN = seq(0, 1, length.out = 1000)) {
  intg_kernelgaus_methyl_elm <- function(pos_snp, rho, id, mat_methyl, RGN) {
    (sum(kernel_gaussian(t = RGN, pos_snp = pos_snp, rho = rho) * mat_methyl[id, ])) / length(RGN)
  }
  sapply(id, intg_kernelgaus_methyl_elm, pos_snp = pos_snp, rho = rho, mat_methyl = mat_methyl, RGN = RGN)
}

# ---- Integral of kernel_linear * methylation ----
intg_kernellin_methyl <- function(id, pos_snp, rho, mat_methyl, RGN = seq(0, 1, length.out = 1000)) {
  intg_kernellin_methyl_elm <- function(pos_snp, rho, id, mat_methyl, RGN) {
    (sum(Kernel_linear(t = RGN, pos_snp = pos_snp, rho = rho) * mat_methyl[id, ])) / length(RGN)
  }
  sapply(id, intg_kernellin_methyl_elm, pos_snp = pos_snp, rho = rho, mat_methyl = mat_methyl, RGN = RGN)
}

# ---- Matrix of integrals kernel * methylation ----
mat_intg_kernel_methyl <- function(pos_snp, id, rho, mat_methyl, kernel_name, RGN = seq(0, 1, length.out = 1000)) {
  if (kernel_name == "convex") {
    all_intg <- as.vector(sapply(pos_snp, intg_kernelvex_methyl, id = id, rho = rho, mat_methyl = mat_methyl, RGN = RGN))
  } else if (kernel_name == "gaussian") {
    all_intg <- as.vector(sapply(pos_snp, intg_kernelgaus_methyl, id = id, rho = rho, mat_methyl = mat_methyl, RGN = RGN))
  } else if (kernel_name == "linear") {
    all_intg <- as.vector(sapply(pos_snp, intg_kernellin_methyl, id = id, rho = rho, mat_methyl = mat_methyl, RGN = RGN))
  } else {
    stop("Unsupported kernel")
  }
  nbr_snp <- length(pos_snp)
  mat_intg <- matrix(all_intg, ncol = nbr_snp)
  return(mat_intg)
}
