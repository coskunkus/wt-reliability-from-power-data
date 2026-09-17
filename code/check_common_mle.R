## =================================================================
##  check_common_mle.R
##
##  The verification of the additions of Section 2.2 and of the new
##  Section 2.5: the sample of N periods in which n_j WTs remain, the
##  clustering constant C, the two mean squared errors that carry the
##  factor 1 + rho C, and the maximum likelihood estimator obtained
##  when the WTs of one period see the same wind.
##
##    group 1  the additions reduce to the published formulas at C = 0
##    group 2  the distribution of m_0 and the variance it implies
##    group 3  the new estimator reduces to the modified MLE at n = 1
##    group 4  the Fisher information, in closed form and numerically
##    group 5  the exact MSE of the new estimator against a simulation
##    group 6  the new estimator never leaves the unit interval
##
##  Everything is computed from the functions of code/core_wtrel.R, so
##  that the paper, the tables and this verification cannot part
##  company.  Output : results/check_common_mle_report.txt
##
##  Base R only.  Run time is about two minutes.
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
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg))
    cand <- c(cand, dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = FALSE)))
  cand <- c(cand, getwd())
  for (d in cand) { r <- find_root(d); if (!is.na(r)) { setwd(r); return(invisible(NULL)) } }
  stop("core_wtrel.R was not found; setwd() to the project root first.")
})
source(file.path("code", "core_wtrel.R"))
new_report(file.path(res_dir, "check_common_mle_report.txt"))

W    <- wind_weibull(PAPER$k, PAPER$alpha)
v_ci <- PAPER$v_ci; v_r <- PAPER$v_r; v_co <- PAPER$v_co; P_r <- PAPER$P_r
d    <- d_wt(W, v_ci, v_co)
m1   <- Eg_wt(W, v_ci, v_r, v_co, P_r, 1)
m2   <- Eg_wt(W, v_ci, v_r, v_co, P_r, 2)
ps   <- PAPER$p_tab
M    <- 30L
n_set<- c(1L, 2L, 3L, 6L)
ok   <- TRUE
check <- function(txt, cond) { ok <<- ok && isTRUE(cond)
  say(if (isTRUE(cond)) "  PASS  " else "  FAIL  ", txt) }

say_rule("the setting")
say("Weibull k = ", PAPER$k, " , alpha = ", PAPER$alpha,
    " ; v_ci = ", v_ci, " , v_r = ", v_r, " , v_co = ", v_co, " , P_r = ", P_r)
say("d = ", fmt5(d), " , m1 = ", fmt4(m1), " , m2 = ", fmt4(m2),
    " , M = ", M, " readings")

## ---- group 1 -----------------------------------------------------
say_rule("group 1: at C = 0 the additions return the published formulas")
e <- c(max(abs(mse_mle_C(ps, M, 0, d)    - mse_mle_wt(ps, M, W, v_ci, v_co))),
       max(abs(mse_moment_C(ps, M, 0, m1, m2) -
               mse_moment_wt(ps, M, W, v_ci, v_r, v_co, P_r))),
       max(abs(sapply(ps, function(p) mse_mmle_C(p, rep(1L, M), d)) -
               sapply(ps, function(p) mse_mmle_wt(p, M, W, v_ci, v_co)))))
say("largest difference, MLE / moment / modified MLE : ",
    paste(format(e, digits = 3), collapse = "  "))
check("the three MSEs of Section 2 are the C = 0 case of the new ones",
      max(e) < 1e-12)
say("C of a sample of five periods of six WTs : ", fmt4(C_design(rep(6L, 5L))),
    "  (n - 1 = 5)")
say("C of a sample with n_j = 6,4,1,1         : ", fmt4(C_design(c(6,4,1,1))))
check("C is n - 1 when every period keeps the same number of WTs",
      abs(C_design(rep(6L, 5L)) - 5) < 1e-12)

## ---- group 2 -----------------------------------------------------
say_rule("group 2: the distribution of m_0 and its variance")
say("the closed form M p d (1-pd) (1 + rho_w C) against the moments of")
say("the distribution of m_0, and against the mechanism itself")
set.seed(20260917, kind = "L'Ecuyer-CMRG")
tab <- NULL
for (n in n_set) {
  N <- M / n; nj <- rep(n, N); p <- 0.95
  pk <- pmf_m0_wt(p, nj, d); k <- 0:M
  v_pmf <- sum(pk * k^2) - sum(pk * k)^2
  v_cf  <- M * p * d * (1 - p * d) * (1 + rho_w_wt(p, d) * C_design(nj))
  reps <- 100000L
  vv <- matrix(qweibull(runif(reps * N), W$k, W$c), nrow = reps)
  b  <- (vv >= v_ci & vv < v_co)
  x  <- matrix(rbinom(reps * N * n, 1L, p), nrow = reps)
  m0 <- M - rowSums(b[, rep(seq_len(N), each = n), drop = FALSE] * x)
  tab <- rbind(tab, data.frame(n = n, N = N, C = C_design(nj),
                               var_pmf = v_pmf, var_cf = v_cf,
                               var_sim = var(m0)))
}
say_block(data.frame(n = tab$n, N = tab$N, C = fmt4(tab$C),
                     from_pmf = fmt4(tab$var_pmf),
                     closed_form = fmt4(tab$var_cf),
                     simulated = fmt4(tab$var_sim)))
check("the closed form agrees with the distribution of m_0",
      max(abs(tab$var_pmf - tab$var_cf)) < 1e-8)
check("and with the mechanism, within two per cent",
      max(abs(tab$var_sim / tab$var_cf - 1)) < 0.02)

## ---- group 3 -----------------------------------------------------
say_rule("group 3: at n = 1 the new estimator is the modified MLE")
nj <- rep(1L, M)
e3 <- max(sapply(0:M, function(m0)
  abs(p_ml_wt(c(rep(1L, m0), rep(0L, M - m0)), nj, d) - p_M_wt(m0, M, d))))
say("largest difference over all m_0 = 0 .. ", M, " : ", format(e3, digits = 3))
check("p_ML equals p_M at every possible value of m_0", e3 < 1e-6)
e3b <- max(abs(sapply(ps, function(p) mse_ml_C(p, M, 1L, d)["mse"]) -
               sapply(ps, function(p) mse_mmle_wt(p, M, W, v_ci, v_co))))
say("largest difference of the two MSEs at n = 1 : ", format(e3b, digits = 3))
check("and so do the two mean squared errors", e3b < 1e-8)

## ---- group 4 -----------------------------------------------------
say_rule("group 4: the Fisher information")
inum <- function(p, n, h = 1e-6) {
  P <- pmf_period(p, n, d)[c(seq_len(n), n + 1L)]
  D <- (pmf_period(p + h, n, d) - pmf_period(p - h, n, d))[c(seq_len(n), n + 1L)] / (2 * h)
  sum(D^2 / P)
}
e4 <- max(abs(outer(c(.90, ps), n_set, Vectorize(function(p, n) info_period_wt(p, n, d))) -
              outer(c(.90, ps), n_set, Vectorize(inum))))
say("largest difference, closed form against numerical : ", format(e4, digits = 3))
check("the closed-form information is right", e4 < 1e-6)
e4b <- max(abs(sapply(ps, function(p) se_ml_wt(p, rep(1L, M), d)^2) -
               mse_mle_wt(ps, M, W, v_ci, v_co)))
say("largest difference, 1 / sum I at n = 1 against Eq. (13) : ",
    format(e4b, digits = 3))
check("at n = 1 the information returns the variance of the MLE", e4b < 1e-12)

## ---- group 5 -----------------------------------------------------
say_rule("group 5: the exact MSE of the new estimator against a simulation")
set.seed(20260917, kind = "L'Ecuyer-CMRG")
reps <- 50000L
sim_ml <- function(p, N, n) {
  vv <- matrix(qweibull(runif(reps * N), W$k, W$c), nrow = reps)
  b  <- (vv >= v_ci & vv < v_co)
  m0 <- matrix(rbinom(reps * N, n, 1 - p), nrow = reps)
  m0 <- ifelse(b, m0, n)                       # out of the band: all idle
  a  <- rowSums(m0 == n)
  S  <- rowSums(m0) - a * n
  key <- paste(a, S); u <- !duplicated(key)
  lut <- vapply(which(u), function(i)
           p_ml_stats(n * (N - a[i]) - S[i], S[i], rep(n, a[i]), d), numeric(1))
  mean((lut[match(key, key[u])] - p)^2)
}

g5 <- do.call(rbind, lapply(n_set, function(n) {
  N <- M / n; p <- 0.95
  data.frame(n = n, N = N,
             exact = unname(mse_ml_C(p, N, n, d)["mse"]),
             simulated = sim_ml(p, N, n),
             bound = n / (M * info_period_wt(p, n, d)))
}))
say_block(data.frame(n = g5$n, N = g5$N, exact = fmt5(g5$exact),
                     simulated = fmt5(g5$simulated), bound = fmt5(g5$bound)))
check("the exact MSE agrees with the simulation, within five per cent",
      max(abs(g5$simulated / g5$exact - 1)) < 0.05)

## ---- group 6 -----------------------------------------------------
say_rule("group 6: the estimator stays in the unit interval")
rng <- range(sapply(n_set, function(n) {
  N <- M / n
  cm <- compositions_wt(as.integer(N), n + 1L)
  range(apply(cm, 1L, function(r) {
    m0j <- c(rep(0:(n - 1L), r[seq_len(n)]), rep(n, r[n + 1L]))
    p_ml_wt(m0j, rep(n, N), d)
  }))
}))
say("over every sample the four designs can produce, p_ML runs from ",
    fmt5(rng[1]), " to ", fmt5(rng[2]))
check("no sample sends p_ML outside [0,1], so no modification is needed",
      rng[1] >= 0 && rng[2] <= 1)

say_rule()
say(if (ok) "all checks passed" else "SOME CHECKS FAILED")
close_report()
if (!ok) stop("check_common_mle.R: some checks failed", call. = FALSE)
