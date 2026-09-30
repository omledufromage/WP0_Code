# Template for a simulation study.
# Copy this file to simulations/NN_<name>.R and edit the sections below.
# Run from the project root (opening Code.Rproj does this for you):
#   Rscript simulations/00_template.R

source("R/utils.R")

# ---- Settings ---------------------------------------------------------------
study  <- "template"
seed   <- 2026
n_reps <- 100

grid <- expand.grid(
    n  = c(50, 200),
    mu = c(0, 0.5)
)

# ---- One replication --------------------------------------------------------
# Generate data and compute whatever you want to record for one replication.
sim_once <- function(params, rep) {
    x  <- rnorm(params$n, mean = params$mu)
    ci <- t.test(x)$conf.int
    list(
        estimate = mean(x),
        covered  = ci[1] <= params$mu && params$mu <= ci[2]
    )
}

# ---- Run and save -----------------------------------------------------------
results <- run_simulation(sim_once, grid, n_reps, seed)

save_results(results, study,
             settings = list(seed = seed, n_reps = n_reps, grid = grid))

# ---- Quick summary ----------------------------------------------------------
print(aggregate(cbind(estimate, covered) ~ n + mu, data = results, FUN = mean))
