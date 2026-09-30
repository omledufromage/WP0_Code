# Shared helpers for all simulation studies.
# Source from any script with: source("R/utils.R")

# Run `sim_fun` for every row of `grid`, `n_reps` times each.
#
# sim_fun: function(params, rep) returning a named list or one-row data.frame
#          of results for a single replication. `params` is a one-row list
#          with the scenario's parameter values.
# grid:    data.frame of scenarios (e.g. built with expand.grid()).
# seed:    base seed. Scenario i uses seed + i, so each scenario is
#          reproducible on its own, independent of which others are run.
#
# Returns a data.frame with one row per (scenario, replication).
run_simulation <- function(sim_fun, grid, n_reps, seed) {
    out <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
        set.seed(seed + i)
        params <- as.list(grid[i, , drop = FALSE])
        reps <- lapply(seq_len(n_reps), function(r) {
            res <- as.data.frame(sim_fun(params, r))
            cbind(scenario = i, rep = r, grid[i, , drop = FALSE], res,
                  row.names = NULL)
        })
        out[[i]] <- do.call(rbind, reps)
        message(sprintf("Scenario %d/%d done", i, nrow(grid)))
    }
    do.call(rbind, out)
}

# Save results to results/<study>/<study>_<timestamp>.rds together with
# the settings used, so every file records how it was produced.
save_results <- function(results, study, settings = list()) {
    dir <- file.path("results", study)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(dir, sprintf("%s_%s.rds", study,
                                     format(Sys.time(), "%Y%m%d_%H%M%S")))
    saveRDS(list(results = results, settings = settings,
                 session = sessionInfo(), created = Sys.time()), path)
    message("Saved: ", path)
    invisible(path)
}

# Load the most recent results file for a study.
load_latest_results <- function(study) {
    files <- list.files(file.path("results", study), pattern = "\\.rds$",
                        full.names = TRUE)
    if (length(files) == 0) stop("No results found for study: ", study)
    readRDS(sort(files, decreasing = TRUE)[1])
}
