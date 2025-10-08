 # funInterMethSNP_SimCont

This repository contains simulation scripts used to reproduce the experiments presented in the paper  
**“A functional approach for testing overall effects between DNA methylation and SNPs”**,  
implemented in the R package [`functintermethsnp`](https://github.com/YvelinGansou/functintermethsnp).



---

##  Context

The model aims to test the overall interaction effect between functional DNA methylation signals and discrete SNP genotypes on a continuous phenotype (e.g., weight or gene expression).

The initial methylation data consist of **8 B-cell samples**, which were used as templates to generate larger sample sizes by duplication (\(N = 200\) or \(N = 400\)) for simulation purposes.  
The **original 8 B-cell methylation profiles** can be made available upon request.

To obtain smooth methylation trajectories for each individual, we applied the function `estim_methyl_curve()` from the [`functintermethsnp`](https://github.com/YvelinGansou/functintermethsnp) package.  
The estimated methylation matrices are stored in the `/data` directory as `.rds` files (`mat_est_meth_n200.rds`, `mat_est_meth_n400.rds`) to avoid recomputation and reduce runtime.

---

##  Repository structure

```
├── scripts/
│ ├── functions_sim.R # Utility functions (kernels, integrals, simulators)
│ ├── sim_H0.R # Type I error simulation (null hypothesis)
│ ├── power_eta_grid.R # Power simulation across η values
│ ├── power_rho_eta_grid.R # Power simulation across (ρ, η) grid
│
├── data/
│ ├── mat_est_meth_n200.rds # Estimated methylation curves for N = 200
│ ├── mat_est_meth_n400.rds # Estimated methylation curves for N = 400
│
├── results/ # Generated automatically (power, p-values, etc.)
└── README.md
```
