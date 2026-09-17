## =================================================================
##  check_formulas.R
##
##  Independent numerical verification of every closed-form expression
##  of the paper.  The script tests the core directly.  A second,
##  independent implementation of the same formulas is in
##  code/check_formulas.py, so that an error in the transcription from
##  the paper to the code cannot survive both.
##
##  Two layers, in this order.
##    Layer 1, exact checks.  Every expression is compared with an
##      independent route: numerical integration in the other variable,
##      finite differences, exact enumeration of all samples for a tiny
##      sample size, a root of the score equation against the closed
##      form, and the Cramer-Rao bound against Eq. (13).
##    Layer 2, empirical checks.  The distributional claims are checked
##      by Monte Carlo with a fixed seed, each reported as a z-score
##      against the Monte Carlo standard error, with |z| < 4 passing.
##      The standard error of a variance is computed from the fourth
##      central moment and not under normality.
##
##  Output : results/check_formulas_report.txt
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
new_report(file.path(res_dir, "check_formulas_report.txt"))

CHK <- new.env()
CHK$n_pass <- 0L; CHK$n_fail <- 0L; CHK$fails <- character(0)

check_rel <- function(name, value, target, tol = 1e-8) {
  rel <- if (abs(target) > 0) abs(value - target) / abs(target)
         else abs(value - target)
  ok <- is.finite(rel) && rel < tol
  if (ok) CHK$n_pass <- CHK$n_pass + 1L else {
    CHK$n_fail <- CHK$n_fail + 1L; CHK$fails <- c(CHK$fails, name)
  }
  say(sprintf("%-58s %-4s  rel.err = %.3e", name, if (ok) "PASS" else "FAIL",
              rel))
  invisible(ok)
}

check_z <- function(name, est, target, se, tol = 4) {
  z  <- (est - target) / se
  ok <- is.finite(z) && abs(z) < tol
  if (ok) CHK$n_pass <- CHK$n_pass + 1L else {
    CHK$n_fail <- CHK$n_fail + 1L; CHK$fails <- c(CHK$fails, name)
  }
  say(sprintf("%-58s %-4s  z = %+7.3f", name, if (ok) "PASS" else "FAIL", z))
  invisible(ok)
}

say("check_formulas.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

## ---- two configurations -----------------------------------------
## C1 is the configuration of Section 3 of the paper, C2 that of the
## real data application (manufacturer's characteristics, Weibull fit).
## the setting of Section 3, taken from the core so that this check and
## the script that produces the table and the figures cannot drift apart
C1 <- list(k = PAPER$k, c = PAPER$alpha, v_ci = PAPER$v_ci,
           v_r = PAPER$v_r, v_co = PAPER$v_co, P_r = PAPER$P_r)
C2 <- list(k = 2.4204, c = 6.8501, v_ci = 3, v_r = 12.5,
           v_co = 24, P_r = 2050)

with_cfg <- function(C) {
  W <- wind_weibull(C$k, C$c)
  list(W = W, v_ci = C$v_ci, v_r = C$v_r, v_co = C$v_co, P_r = C$P_r,
       d = d_wt(W, C$v_ci, C$v_co))
}

## =================================================================
##  Layer 1.  Exact checks
## =================================================================
for (nm in c("C1", "C2")) {
  C <- get(nm); E <- with_cfg(C)
  W <- E$W; v_ci <- E$v_ci; v_r <- E$v_r; v_co <- E$v_co; P_r <- E$P_r
  say_rule(paste0("Layer 1, exact checks, configuration ", nm,
                  ":  k = ", C$k, " , c = ", C$c, " , v_ci = ", C$v_ci,
                  " , v_r = ", C$v_r, " , v_co = ", C$v_co,
                  " , P_r = ", C$P_r))

  ## --- Eq. (2) -------------------------------------------------
  check_rel(paste0(nm, "  Eq. (2): g(v_ci) = 0"),
            g_wt(v_ci, v_ci, v_r, v_co, P_r), 0, 1e-12)
  check_rel(paste0(nm, "  Eq. (2): g is continuous at v_r"),
            g_wt(v_r - 1e-9, v_ci, v_r, v_co, P_r), P_r, 1e-8)
  check_rel(paste0(nm, "  Eq. (2): g(v_co) = 0"),
            g_wt(v_co, v_ci, v_r, v_co, P_r), 0, 1e-12)

  ## --- Eq. (3), H1, against integration in the wind speed scale --
  xs <- P_r * c(0.05, 0.3, 0.6, 0.9)
  for (x in xs) {
    v_x <- ((x / P_r) * (v_r^3 - v_ci^3) + v_ci^3)^(1 / 3)
    ind <- integrate(W$f, 0, v_x, rel.tol = 1e-12)$value +
      integrate(W$f, v_co, Inf, rel.tol = 1e-12)$value
    check_rel(sprintf("%s  Eq. (3): H1(%.4g) by integration", nm, x),
              H1_wt(x, W, v_ci, v_r, v_co, P_r), ind, 1e-8)
  }

  ## --- the two atoms of Eq. (4) --------------------------------
  for (p in c(0.5, 0.9, 0.99)) {
    ## P{P_WT = 0} = p P{V < v_ci or V >= v_co} + (1 - p), by integration
    tgt <- p * (integrate(W$f, 0, v_ci, rel.tol = 1e-12)$value +
                  integrate(W$f, v_co, Inf, rel.tol = 1e-12)$value) + (1 - p)
    check_rel(sprintf("%s  Eq. (4): P{P_WT = 0}, p = %.2f", nm, p),
              P_zero_wt(p, W, v_ci, v_co), tgt, 1e-8)
    tgt2 <- p * integrate(W$f, v_r, v_co, rel.tol = 1e-12)$value
    check_rel(sprintf("%s  Eq. (4): P{P_WT = P_r}, p = %.2f", nm, p),
              P_rated_wt(p, W, v_r, v_co), tgt2, 1e-8)
    ## total mass one, the absolutely continuous part by integration
    h1 <- function(x) {
      v_x <- ((x / P_r) * (v_r^3 - v_ci^3) + v_ci^3)^(1 / 3)
      W$f(v_x) * (v_r^3 - v_ci^3) / (3 * P_r * v_x^2)
    }
    mass <- P_zero_wt(p, W, v_ci, v_co) + P_rated_wt(p, W, v_r, v_co) +
      p * integrate(h1, 0, P_r, rel.tol = 1e-12)$value
    check_rel(sprintf("%s  Eq. (4): total mass one, p = %.2f", nm, p),
              mass, 1, 1e-7)
    ## h1 against the finite difference of H1, Eq. (6)
    x0 <- 0.4 * P_r; e <- 1e-5 * P_r
    fd <- (H1_wt(x0 + e, W, v_ci, v_r, v_co, P_r) -
             H1_wt(x0 - e, W, v_ci, v_r, v_co, P_r)) / (2 * e)
    check_rel(sprintf("%s  Eq. (6): h1 by finite difference", nm),
              h1(x0), fd, 1e-6)
  }

  ## --- Eqs. (5) and (10) in the power scale ---------------------
  h1 <- function(x) {
    v_x <- ((x / P_r) * (v_r^3 - v_ci^3) + v_ci^3)^(1 / 3)
    W$f(v_x) * (v_r^3 - v_ci^3) / (3 * P_r * v_x^2)
  }
  m1_x <- integrate(function(x) x * h1(x), 0, P_r, rel.tol = 1e-12)$value +
    P_r * (W$F(v_co) - W$F(v_r))
  m2_x <- integrate(function(x) x^2 * h1(x), 0, P_r, rel.tol = 1e-12)$value +
    P_r^2 * (W$F(v_co) - W$F(v_r))
  check_rel(paste0(nm, "  Eq. (5): mean, power scale against wind scale"),
            Eg_wt(W, v_ci, v_r, v_co, P_r, 1), m1_x, 1e-7)
  check_rel(paste0(nm, "  Eq. (10): second moment, both scales"),
            Eg_wt(W, v_ci, v_r, v_co, P_r, 2), m2_x, 1e-7)

  ## --- Eqs. (9) and (13) by exact enumeration -------------------
  for (n in c(5L, 8L)) for (p in c(0.6, 0.9, 0.98)) {
    d <- E$d; k <- 0:n
    pk <- dbinom(k, n, 1 - p * d)
    ph <- (1 - k / n) / d
    check_rel(sprintf("%s  Eq. (12): E[p_hat] = p, n = %d, p = %.2f", nm, n, p),
              sum(ph * pk), p, 1e-10)
    check_rel(sprintf("%s  Eq. (13): MSE by enumeration, n = %d, p = %.2f",
                      nm, n, p),
              sum((ph - p)^2 * pk), mse_mle_wt(p, n, W, v_ci, v_co), 1e-10)
    pm <- pmin(ph, 1)
    check_rel(sprintf("%s  MSE of p_M by enumeration, n = %d, p = %.2f",
                      nm, n, p),
              sum((pm - p)^2 * pk), mse_mmle_wt(p, n, W, v_ci, v_co), 1e-10)
    ## Eq. (9) by the variance of the power output, an independent route
    mse9 <- (mu2_wt(p, W, v_ci, v_r, v_co, P_r) -
               mu_wt(p, W, v_ci, v_r, v_co, P_r)^2) /
      (n * Eg_wt(W, v_ci, v_r, v_co, P_r, 1)^2)
    check_rel(sprintf("%s  Eq. (9): MSE from Var(P_WT), n = %d, p = %.2f",
                      nm, n, p),
              mse_moment_wt(p, n, W, v_ci, v_r, v_co, P_r), mse9, 1e-10)
  }

  ## --- the score equation and the information -------------------
  n <- 40L; p0 <- 0.95; d <- E$d
  m0 <- 4L
  loglik <- function(p) m0 * log(1 - p * d) + (n - m0) * log(p * d)
  ## the score of Eq. (11), written analytically
  score  <- function(p) -m0 * d / (1 - p * d) + (n - m0) / p
  ph_cf  <- p_hat_wt(m0, n, d)
  check_rel(paste0(nm, "  Eq. (12): root of the score equation"),
            uniroot(score, c(1e-6, 1 / d - 1e-6), tol = 1e-14)$root,
            ph_cf, 1e-8)
  ## the same score by a central difference of the log-likelihood
  h <- 1e-6
  check_rel(paste0(nm, "  Eq. (11): numerical score at the MLE is zero"),
            (loglik(ph_cf + h) - loglik(ph_cf - h)) / (2 * h * n), 0, 1e-5)
  d2 <- (loglik(ph_cf + 1e-5) - 2 * loglik(ph_cf) + loglik(ph_cf - 1e-5)) /
    1e-10
  check_rel(paste0(nm, "  Eq. (11): observed information is positive"),
            as.numeric(-d2 > 0), 1, 1e-12)
  ## Cramer-Rao bound: the inverse of the Fisher information equals Eq. (13)
  fisher <- function(p) n * d^2 / ((1 - p * d) * (p * d))
  check_rel(paste0(nm, "  Eq. (13) attains the Cramer-Rao bound"),
            1 / fisher(p0), mse_mle_wt(p0, n, W, v_ci, v_co), 1e-12)

}

## --- Table 1 of the paper, reproduced exactly --------------------
say_rule("Layer 1, Table 1 of the paper reproduced from Eqs. (13) and the MSE of p_M")
W1 <- wind_weibull(PAPER$k, PAPER$alpha)
tab1 <- data.frame(
  p = rep(PAPER$p_tab, times = length(PAPER$n_tab)),
  n = rep(PAPER$n_tab, each = length(PAPER$p_tab)))
tab1$mse_pM   <- mapply(function(p, n)
  mse_mmle_wt(p, n, W1, PAPER$v_ci, PAPER$v_co), tab1$p, tab1$n)
tab1$mse_phat <- mapply(function(p, n)
  mse_mle_wt(p, n, W1, PAPER$v_ci, PAPER$v_co), tab1$p, tab1$n)
say_block(data.frame(p = tab1$p, n = tab1$n,
                     MSE_pM = fmt4(tab1$mse_pM),
                     MSE_phat = fmt4(tab1$mse_phat)))
## the values printed in the paper, held in the core
paper <- PAPER_TABLE1
ours <- as.numeric(t(cbind(round(tab1$mse_pM, 4), round(tab1$mse_phat, 4))))
check_rel("Table 1: all 32 entries agree to four decimals",
          max(abs(ours - paper)), 0, 1e-12)

## =================================================================
##  Layer 2.  Empirical checks by Monte Carlo
## =================================================================
say_rule("Layer 2, Monte Carlo checks, |z| < 4 passes")
set.seed(20260909L)
B <- 200000L
E <- with_cfg(C2)
W <- E$W; v_ci <- E$v_ci; v_r <- E$v_r; v_co <- E$v_co; P_r <- E$P_r; d <- E$d
se_var <- function(x) {
  n <- length(x); m4 <- mean((x - mean(x))^4); v <- var(x)
  sqrt(max(m4 - v^2, 0) / n)
}

for (p in c(0.90, 0.97)) for (n in c(100L, 400L)) {
  m0 <- rbinom(B, n, 1 - p * d)
  ph <- p_hat_wt(m0, n, d)
  pm <- pmin(ph, 1)
  check_z(sprintf("MC  E[p_hat] = p, n = %d, p = %.2f", n, p),
          mean(ph), p, sd(ph) / sqrt(B))
  check_z(sprintf("MC  Var(p_hat) = Eq. (13), n = %d, p = %.2f", n, p),
          var(ph), mse_mle_wt(p, n, W, v_ci, v_co), se_var(ph))
  check_z(sprintf("MC  MSE(p_M) formula, n = %d, p = %.2f", n, p),
          mean((pm - p)^2), mse_mmle_wt(p, n, W, v_ci, v_co),
          sd((pm - p)^2) / sqrt(B))
}

## the atoms of Eq. (4) and the moment estimator, by simulating the model
for (p in c(0.90, 0.97)) {
  n <- 400L; Bm <- 40000L
  V <- matrix(rweibull(n * Bm, W$k, W$c), nrow = Bm)
  X <- matrix(rbinom(n * Bm, 1L, p), nrow = Bm)
  PW <- g_wt(as.vector(V), v_ci, v_r, v_co, P_r) * as.vector(X)
  dim(PW) <- c(Bm, n)
  z0 <- mean(PW <= 0)
  check_z(sprintf("MC  P{P_WT = 0} = 1 - pd, p = %.2f", p),
          z0, P_zero_wt(p, W, v_ci, v_co),
          sqrt(z0 * (1 - z0) / (Bm * n)))
  zr <- mean(PW >= P_r - 1e-9)
  check_z(sprintf("MC  P{P_WT = P_r}, p = %.2f", p),
          zr, P_rated_wt(p, W, v_r, v_co), sqrt(zr * (1 - zr) / (Bm * n)))
  pt <- rowMeans(PW) / Eg_wt(W, v_ci, v_r, v_co, P_r, 1)
  check_z(sprintf("MC  E[p_tilde] = p, n = %d, p = %.2f", n, p),
          mean(pt), p, sd(pt) / sqrt(Bm))
  check_z(sprintf("MC  Var(p_tilde) = Eq. (9), n = %d, p = %.2f", n, p),
          var(pt), mse_moment_wt(p, n, W, v_ci, v_r, v_co, P_r), se_var(pt))
  rm(V, X, PW)
}

## the rate claimed in Section 2.2: p_hat and p_M are asymptotically
## equivalent, so sqrt(n) E|p_hat - p_M| tends to zero
say("")
say("rate check: mean of sqrt(n) |p_hat - p_M| must fall as n doubles")
p <- 0.97
prev <- NA_real_
ok_rate <- TRUE
for (n in c(200L, 400L, 800L, 1600L)) {
  m0 <- rbinom(B, n, 1 - p * d)
  ph <- p_hat_wt(m0, n, d)
  val <- mean(sqrt(n) * abs(ph - pmin(ph, 1)))
  say(sprintf("   n = %5d   sqrt(n) E|p_hat - p_M| = %.6e", n, val))
  if (!is.na(prev) && val > prev) ok_rate <- FALSE
  prev <- val
}
if (ok_rate) CHK$n_pass <- CHK$n_pass + 1L else {
  CHK$n_fail <- CHK$n_fail + 1L; CHK$fails <- c(CHK$fails, "rate of p_hat - p_M")
}
say(sprintf("%-58s %s", "rate: the sequence is decreasing",
            if (ok_rate) "PASS" else "FAIL"))

## =================================================================
say_rule("Summary")
say("checks passed : ", CHK$n_pass)
say("checks failed : ", CHK$n_fail)
if (CHK$n_fail > 0L) say("failed: ", paste(CHK$fails, collapse = " ; "))
say(if (CHK$n_fail == 0L) "ALL CHECKS PASS" else "SOME CHECKS FAIL")
close_report()
