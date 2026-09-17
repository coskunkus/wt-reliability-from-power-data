## =================================================================
##  run_simulation.R
##
##  The sampling design of Step 3 replicated on the panel of the site.
##
##  Steps 3 and 4 of Section 2.2 prescribe random instants, balanced
##  over the twelve months and over day and night and at least 48 hours
##  apart, at which the n WTs are inspected and m0 is counted.  Here
##  that selection is repeated R times on the panel, so that the
##  sampling variability of the estimators in the field can be set
##  against the exact expressions (9) and (13) of the paper.
##
##  Only the maximum likelihood estimator is examined here, since the
##  moment estimator is not admissible on these data.
##
##  Two designs:
##    A  all six WTs are read at each selected instant, which is Step 4
##       as written; the six readings share the wind of that instant
##    B  one WT, chosen at random, is read at each selected instant, so
##       that the N observations are taken at N different times
##
##  Parallel through the base package parallel, with one L'Ecuyer
##  stream per replication, so that the results do not depend on the
##  number of workers.
##
##  Input  : data/kelmarsh_panel.rds, data/run_all_objects.rds
##  Output : results/simulation_report.txt, which contains the numbers
##           quoted in the last paragraph of Section 4
##
##  Base R only.  Run time is about a minute.
## =================================================================

## ---- locate the project root ------------------------------------
local({
  find_root <- function(d) {
    if (file.exists(file.path(d, "code", "core_wtrel.R"))) return(d)
    if (file.exists(file.path(d, "core_wtrel.R")))
      return(if (basename(d) == "code") dirname(d) else d)
    NA_character_
  }
  cand <- character(0)
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) cand <- c(cand, dirname(normalizePath(of, mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    pth <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(pth, "try-error") && length(pth) == 1L && nzchar(pth))
      cand <- c(cand, dirname(normalizePath(pth, mustWork = FALSE)))
  }
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg))
    cand <- c(cand, dirname(normalizePath(sub("^--file=", "", arg[1]),
                                          mustWork = FALSE)))
  cand <- c(cand, getwd())
  for (d in cand) {
    r <- find_root(d)
    if (!is.na(r)) { setwd(r); return(invisible(NULL)) }
  }
  stop("core_wtrel.R was not found; setwd() to the project root first.")
})
source(file.path("code", "core_wtrel.R"))
library(parallel)
if (!dir.exists(res_dir)) dir.create(res_dir)

## ---- the master seed, quoted in the paper -----------------------
MASTER_SEED <- 20260909L
R_DESIGN    <- 1000L       # replications of each design
PER_STRATUM <- 22L         # instants per (month x day/night) stratum
GAP_SLOTS   <- 288L        # 48 hours in ten-minute slots

new_report(file.path(res_dir, paste0(tag("simulation"), "_report.txt")))
say("run_simulation.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)
say("master seed: ", MASTER_SEED, " ,  RNG: L'Ecuyer-CMRG, ",
    "one stream per replication")

need_file(file.path("data", paste0(tag("run_all_objects"), ".rds")), "code/run_all.R")
need_file(PANEL_RDS, "code/prepare_data.R")
O   <- readRDS(file.path("data", paste0(tag("run_all_objects"), ".rds")))
D   <- readRDS(PANEL_RDS)
need_panel(D)
pan <- D$panel
## The manufacturer's cut-in and cut-out speeds, and the two quantities
## that Eq. (12) needs, exactly as the section uses them: d is the
## proportion of the records of the interval, read off the wind speed
## record of the sample, and Eq. (13) is evaluated at the estimate that
## the section reports.
v_ci <- O$v_ci; v_co <- O$v_co
d_use <- O$d_rec; p_ref <- O$p_rec

## ---- L'Ecuyer streams -------------------------------------------
make_streams <- function(n, seed) {
  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old[1]), add = TRUE)
  set.seed(seed)
  s <- .Random.seed
  out <- vector("list", n)
  for (i in seq_len(n)) { out[[i]] <- s; s <- nextRNGStream(s) }
  out
}

## ---- the balanced panel as one matrix, instants by WTs -----------
## This is the one place where the balanced panel is used: the design of
## Step 4 reads every WT at the selected instant, so that only instants
## at which every WT is retained can serve.  The records are those of
## the estimation sample of the core, restricted to such instants, and
## the zero-power indicator is the one of the core, the assumption of
## Eq. (2) included.
sub  <- estimation_sample(pan)
sub  <- sub[sub$full_instant, ]
inst <- sort(unique(sub$slot))
n_inst <- length(inst)
n_slot <- D$n_slot
row_of <- integer(n_slot); row_of[inst + 1L] <- seq_len(n_inst)
sub    <- sub[order(sub$slot, sub$turbine), ]
zero_mat <- matrix(as.numeric(zero_of(sub, v_ci, v_co)),
                   nrow = n_inst, byrow = TRUE)
stratum  <- make_strata(inst * 600 + D$t_origin)
## the panel itself is of no further use and is large, so it goes now
rm(pan, sub, D); invisible(gc(FALSE))

say_rule("Step 3 replicated on the panel")
say("admissible instants, all WTs retained     : ", fmti(n_inst))
say("strata, month x day/night                 : ", length(unique(stratum)))
say("instants per stratum                      : ", PER_STRATUM)
say("minimum separation                        : ", GAP_SLOTS / 6, " hours")
say("replications of each design               : ", fmti(R_DESIGN))
say("d used in Eq. (12), wind speed record     : ", fmt5(d_use))
say("p at which Eq. (13) is evaluated          : ", fmt5(p_ref),
    "  (the estimate of the section)")
say("zero-power rule                           : the one of the core, ",
    "a record outside [", v_ci, ", ", v_co, ") counting as zero")

one_rep <- function(i, mode, streams) {
  assign(".Random.seed", streams[[i]], envir = globalenv())
  sl  <- draw_slots(inst, stratum, PER_STRATUM, GAP_SLOTS, n_slot)
  idx <- row_of[sl + 1L]
  if (mode == "A") {
    m0 <- sum(zero_mat[idx, ]); N <- length(idx) * ncol(zero_mat)
  } else {
    j  <- sample.int(ncol(zero_mat), length(idx), replace = TRUE)
    m0 <- sum(zero_mat[cbind(idx, j)]); N <- length(idx)
  }
  ph <- p_hat_wt(m0, N, d_use)
  c(T = length(sl), N = N, m0 = m0, phat = ph, pM = min(ph, 1))
}

n_core <- max(1L, min(detectCores(logical = FALSE), 8L))
say("workers                                   : ", n_core)
## The cluster is stopped explicitly below.  on.exit() is deliberately
## not used, because at the top level of a file that is sourced it may
## be executed at once, which would stop the cluster before it is used.
cl <- makeCluster(n_core)
## A worker cannot source the core when the project path contains
## non-ASCII characters, so the objects are exported instead.  Only the
## objects that one replication needs are exported: exporting the whole
## workspace would send the panel to every worker, which is large enough
## to break the connection to the workers when several scripts have been
## run in the same session.
clusterExport(cl, c("inst", "stratum", "PER_STRATUM", "GAP_SLOTS", "n_slot",
                    "row_of", "zero_mat", "d_use", "one_rep", "draw_slots",
                    "p_hat_wt"),
              envir = environment())

run_design <- function(mode, tag) {
  streams <- make_streams(R_DESIGN, MASTER_SEED + 7919L * (mode == "B"))
  one <- one_rep(1L, mode, streams)          # one replication on the master
  stopifnot(all(is.finite(one)))
  chunks <- splitIndices(R_DESIGN, n_core)
  res <- do.call(rbind, clusterApply(cl, chunks, function(ii)
    do.call(rbind, lapply(ii, function(i)
      tryCatch(one_rep(i, mode, streams),
               error = function(e) rep(NA_real_, 5L))))))
  colnames(res) <- names(one)
  res <- res[stats::complete.cases(res), , drop = FALSE]
  N   <- mean(res[, "N"])
  data.frame(design = tag,
             T = mean(res[, "T"]), N = N,
             mean_phat = mean(res[, "phat"]), sd_phat = sd(res[, "phat"]),
             rmse13 = se_mle_wt(p_ref, N, d_use),
             mean_pM = mean(res[, "pM"]), sd_pM = sd(res[, "pM"]),
             reps = nrow(res), row.names = NULL)
}
## rmse13 is the square root of Eq. (13) evaluated with the d that the
## estimator actually uses.

tab <- rbind(run_design("A", "A"), run_design("B", "B"))
stopCluster(cl)

say("")
say("A: all six WTs read at each instant (Step 4 as written)")
say("B: one WT read at each instant")
say_block(format(tab, digits = 5))
say("")
say("ratio of the simulated variance of p_hat to Eq. (13)")
say("  design A : ", fmt4((tab$sd_phat[1] / tab$rmse13[1])^2))
say("  design B : ", fmt4((tab$sd_phat[2] / tab$rmse13[2])^2))

saveRDS(list(tab = tab, master_seed = MASTER_SEED, R_design = R_DESIGN),
        file.path("data", paste0(tag("design_study"), ".rds")))

say("")
say("written: data/", tag("design_study"), ".rds")
close_report()
