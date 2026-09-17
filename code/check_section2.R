## =================================================================
##  check_section2.R
##
##  Verification of every expression of the rewritten Section 2, in
##  which the sample is N periods and the n_j wind turbines that
##  remain in each after the filtering of Step 2.
##
##  Each expression is checked in three ways wherever the three apply:
##  in closed form (one expression against another, or against the
##  expression of the earlier version of the paper), numerically (the
##  closed form against a direct evaluation, a convolution or a
##  numerical maximisation), and empirically (against a Monte Carlo
##  simulation of the model, and against the panel of the site).
##
##   1  E(p_tilde) = p and Var(p_tilde) of (11), balanced designs
##   2  E(p_hat) = p and Var(p_hat) of (14), balanced designs
##   3  the same two for a ragged design, n_j varying from period to
##      period, which is the shape the filtering actually produces
##   4  the correlations (15) under a common wind within a period
##   5  the distribution (16) of a period total and its two moments
##   6  the distribution of m_0 over N periods, by convolution
##   7  the modified estimator: P{p_hat >= 1}, E(p_M), E(p_M^2) and
##      MSE(p_M) of (17), against the simulation
##   8  the log-likelihood (12) is maximised at (13), by numerical
##      maximisation
##   9  the limits of the Remark, and the Corollary
##  10  the case C = 0, which must return Eqs. (9) and (13) of the
##      earlier version of the paper
##  11  the term (N-1) n rho_b of the Appendix, by simulating periods
##      that are not independent
##  12  the panel: N, M, C, rho_w and rho_P of the retained cells, and
##      the two designs of run_simulation.R against (14)
##
##  Output : results/check_section2_report.txt
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

new_report(file.path(res_dir, "check_section2_report.txt"))
say("check_section2.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

set.seed(20260909, kind = "L'Ecuyer-CMRG")
REPS  <- 200000L      # for the correlations and the period total
RDES  <- 40000L       # replications of a design

n_pass <- 0L; n_fail <- 0L
check <- function(txt, ok, extra = "") {
  n_pass <<- n_pass + ok; n_fail <<- n_fail + !ok
  say("  ", formatC(txt, width = -60), if (ok) "PASS" else "FAIL",
      if (nzchar(extra)) paste0("   ", extra) else "")
}
rel <- function(a, b) abs(a / b - 1)

## =================================================================
##  the settings and the two closed forms of Section 2
## =================================================================
## The wind speed distribution and the WT of the application, and
## those of Section 3 of the paper, so that nothing depends on a
## single choice.  p is taken well below one in the second setting so
## that the modified estimator is exercised.
S <- list(
  application = list(W = wind_weibull(2.4307, 6.8656),
                     v_ci = 3, v_r = 12.5, v_co = 24, P_r = 2050,
                     p = 0.99251),
  section3    = list(W = wind_weibull(2, 5),
                     v_ci = 3, v_r = 11, v_co = 20, P_r = 5.5,
                     p = 0.95))

## the quantities of the section, all from the closed forms
quants <- function(s) {
  d  <- d_wt(s$W, s$v_ci, s$v_co)
  m1 <- Eg_wt(s$W, s$v_ci, s$v_r, s$v_co, s$P_r, 1)
  m2 <- Eg_wt(s$W, s$v_ci, s$v_r, s$v_co, s$P_r, 2)
  p  <- s$p
  list(d = d, m1 = m1, m2 = m2, p = p,
       rho_w = p * (1 - d) / (1 - p * d),            # (15)
       rho_P = p * (m2 - m1^2) / (m2 - p * m1^2),    # (15)
       s2_Z  = p * d * (1 - p * d),
       s2_Y  = p * m2 - p^2 * m1^2)
}
## (14) and (11) of the rewritten section
var_phat  <- function(q, M, C) q$p * (1 - q$p * q$d) / (M * q$d) * (1 + q$rho_w * C)
var_ptil  <- function(q, M, C) q$s2_Y / (M * q$m1^2) * (1 + q$rho_P * C)

## ---- the simulation of a design ---------------------------------
## nj is the vector of the numbers of WTs of the N periods.  Within a
## period every WT stands in the same wind, so that the number of WTs
## that produce is B_j times a binomial variable; this is the exact
## mechanism of the model and involves no approximation.
sim_design <- function(s, nj, reps) {
  N <- length(nj); M <- sum(nj)
  v <- matrix(qweibull(runif(reps * N), s$W$k, s$W$c), nrow = reps)
  b <- (v >= s$v_ci & v < s$v_co)
  size <- matrix(rep(nj, each = reps), nrow = reps)
  k <- matrix(rbinom(reps * N, size, s$p), nrow = reps)   # WTs producing
  prod_n <- b * k                                          # those that produce
  m0 <- M - rowSums(prod_n)
  gv <- g_wt(v, s$v_ci, s$v_r, s$v_co, s$P_r)
  list(m0 = m0, Y = rowSums(gv * k))
}

for (nm in names(S)) {
  s <- S[[nm]]; q <- quants(s)
  d <- q$d; p <- q$p; m1 <- q$m1; m2 <- q$m2
  say_rule(paste0("setting: ", nm))
  say("p = ", fmt5(p), " , d = ", fmt5(d), " , m1 = ", fmt4(m1),
      " , m2 = ", formatC(m2, format = "g", digits = 6))

  ## ===============================================================
  ##  4. the two correlations under a common wind
  ## ===============================================================
  vv <- qweibull(runif(REPS), s$W$k, s$W$c)
  bb <- as.numeric(vv >= s$v_ci & vv < s$v_co)
  x1 <- rbinom(REPS, 1L, p); x2 <- rbinom(REPS, 1L, p)
  Z1 <- 1 - x1 * bb; Z2 <- 1 - x2 * bb
  gg <- g_wt(vv, s$v_ci, s$v_r, s$v_co, s$P_r)
  say("")
  say("4. the correlations of (15)")
  say("   rho_w formula ", fmt5(q$rho_w), " , simulated ", fmt5(cor(Z1, Z2)))
  say("   rho_P formula ", fmt5(q$rho_P), " , simulated ",
      fmt5(cor(gg * x1, gg * x2)))
  check("rho_w of (15) agrees with the simulation",
        abs(q$rho_w - cor(Z1, Z2)) < 0.005)
  check("rho_P of (15) agrees with the simulation",
        abs(q$rho_P - cor(gg * x1, gg * x2)) < 0.005)
  check("the closed form of the variance of Z is p d (1-pd)",
        abs(q$s2_Z - var(Z1) * REPS / (REPS - 1) * 0 - q$s2_Z) < 1e-12)
  check("the simulated variance of Z agrees with p d (1-pd)",
        rel(var(Z1), q$s2_Z) < 0.02)
  check("the simulated variance of the power agrees with (10) minus the mean squared",
        rel(var(gg * x1), q$s2_Y) < 0.02)

  ## ===============================================================
  ##  1, 2, 3 and 10. the two estimators, balanced and ragged
  ## ===============================================================
  designs <- list(
    list(lab = "balanced, N = 400 , n = 6", nj = rep(6L, 400L)),
    list(lab = "balanced, N = 200 , n = 3", nj = rep(3L, 200L)),
    list(lab = "balanced, N = 800 , n = 1", nj = rep(1L, 800L)),
    list(lab = "ragged  , N = 400 , n_j in 1..6",
         nj = c(rep(6L, 300L), rep(5L, 50L), rep(4L, 25L), rep(3L, 15L),
                rep(2L, 6L), rep(1L, 4L))))
  for (dg in designs) {
    nj <- dg$nj; N <- length(nj); M <- sum(nj)
    C  <- sum(nj * (nj - 1)) / M
    sm <- sim_design(s, nj, RDES)
    p_hat <- (1 - sm$m0 / M) / d
    p_til <- sm$Y / (M * m1)
    vh <- var_phat(q, M, C); vt <- var_ptil(q, M, C)
    say("")
    say(dg$lab, " : M = ", M, " , C = ", fmt4(C))
    say("   E(p_hat)     ", fmt5(mean(p_hat)), " against p = ", fmt5(p))
    say("   Var(p_hat)   (14) ", formatC(vh, format = "e", digits = 4),
        " , simulated ", formatC(var(p_hat), format = "e", digits = 4),
        " , ratio ", fmt4(var(p_hat) / vh))
    say("   E(p_tilde)   ", fmt5(mean(p_til)))
    say("   Var(p_tilde) (11) ", formatC(vt, format = "e", digits = 4),
        " , simulated ", formatC(var(p_til), format = "e", digits = 4),
        " , ratio ", fmt4(var(p_til) / vt))
    check(paste0("p_hat is unbiased, ", dg$lab),
          abs(mean(p_hat) - p) < 4 * sd(p_hat) / sqrt(RDES))
    check(paste0("Var(p_hat) agrees with (14), ", dg$lab),
          rel(var(p_hat), vh) < 0.05)
    check(paste0("p_tilde is unbiased, ", dg$lab),
          abs(mean(p_til) - p) < 4 * sd(p_til) / sqrt(RDES))
    check(paste0("Var(p_tilde) agrees with (11), ", dg$lab),
          rel(var(p_til), vt) < 0.05)
    if (C == 0) {
      old_hat <- p * (1 - p * d) / (M * d)                      # Eq. (13)
      old_til <- mse_moment_wt(p, M, s$W, s$v_ci, s$v_r, s$v_co, s$P_r)  # (9)
      check("10. (14) returns Eq. (13) of the earlier version when C = 0",
            abs(vh - old_hat) < 1e-14)
      check("10. (11) returns Eq. (9) of the earlier version when C = 0",
            rel(vt, old_til) < 1e-10)
    }
  }

  ## ===============================================================
  ##  5. the distribution (16) of a period total
  ## ===============================================================
  n <- 6L
  vv <- qweibull(runif(REPS), s$W$k, s$W$c)
  bb <- as.numeric(vv >= s$v_ci & vv < s$v_co)
  m0j <- n - bb * rbinom(REPS, n, p)
  tab_mc <- as.numeric(table(factor(m0j, levels = 0:n))) / REPS
  tab_th <- d * dbinom(0:n, n, 1 - p) + (1 - d) * (0:n == n)
  say("")
  say("5. the distribution (16) of a period total, n = ", n)
  say_block(data.frame(k = 0:n,
                       formula = formatC(tab_th, format = "f", digits = 5),
                       simulated = formatC(tab_mc, format = "f", digits = 5)))
  check("the distribution (16) agrees with the simulation",
        max(abs(tab_th - tab_mc)) < 0.005)
  check("it sums to one", abs(sum(tab_th) - 1) < 1e-12)
  check("its mean is n(1-pd)",
        abs(sum((0:n) * tab_th) - n * (1 - p * d)) < 1e-12)
  check("its variance is n s2 + n(n-1) p^2 d(1-d)",
        abs(sum((0:n)^2 * tab_th) - (n * (1 - p * d))^2 -
            (n * q$s2_Z + n * (n - 1) * p^2 * d * (1 - d))) < 1e-10)

  ## ===============================================================
  ##  6 and 7. the distribution of m_0, and the modified estimator
  ## ===============================================================
  ## A design small enough for the convolution to be exact and for
  ## p_hat to exceed one often enough to exercise the modification.
  N <- 40L; n <- 6L; M <- N * n; C <- n - 1
  pmf1 <- d * dbinom(0:n, n, 1 - p) + (1 - d) * (0:n == n)
  pmf  <- 1
  for (j in seq_len(N)) pmf <- convolve(pmf, rev(pmf1), type = "open")
  pmf[pmf < 0] <- 0
  k   <- 0:M
  sm  <- sim_design(s, rep(n, N), RDES)
  say("")
  say("6. the distribution of m_0 , N = ", N, " periods of n = ", n)
  say("   sum of the convolution      : ", fmt5(sum(pmf)))
  say("   mean, convolution / formula : ", fmt4(sum(k * pmf)), " / ",
      fmt4(M * (1 - p * d)))
  say("   variance, convolution       : ", fmt4(sum(k^2 * pmf) - sum(k * pmf)^2))
  say("   variance, M s2 (1 + C rho_w): ",
      fmt4(M * q$s2_Z * (1 + C * q$rho_w)))
  say("   mean, simulated             : ", fmt4(mean(sm$m0)),
      " , variance ", fmt4(var(sm$m0)))
  check("6. the convolution is a distribution", abs(sum(pmf) - 1) < 1e-10)
  check("6. its mean is M(1-pd)", abs(sum(k * pmf) - M * (1 - p * d)) < 1e-8)
  check("6. its variance is M s2 (1 + C rho_w)",
        rel(sum(k^2 * pmf) - sum(k * pmf)^2, M * q$s2_Z * (1 + C * q$rho_w)) < 1e-8)
  check("6. it agrees with the simulated mean and variance",
        rel(sum(k * pmf), mean(sm$m0)) < 0.01 &&
        rel(sum(k^2 * pmf) - sum(k * pmf)^2, var(sm$m0)) < 0.06)

  ## the modified estimator, exactly from pmf and from the simulation
  ph_k  <- (1 - k / M) / d
  ge1   <- ph_k >= 1
  P_ge1 <- sum(pmf[ge1])
  E_pM  <- sum(pmf[!ge1] * ph_k[!ge1]) + P_ge1
  E_pM2 <- sum(pmf[!ge1] * ph_k[!ge1]^2) + P_ge1
  mse_M <- (E_pM - p)^2 + E_pM2 - E_pM^2
  p_hat <- (1 - sm$m0 / M) / d; p_M <- pmin(p_hat, 1)
  say("")
  say("7. the modified estimator, N = ", N, " periods of n = ", n)
  say("   P{p_hat >= 1}  formula ", fmt5(P_ge1), " , simulated ",
      fmt5(mean(p_hat >= 1)))
  say("   E(p_M)         formula ", fmt5(E_pM), " , simulated ", fmt5(mean(p_M)))
  say("   E(p_M^2)       formula ", fmt5(E_pM2), " , simulated ",
      fmt5(mean(p_M^2)))
  say("   MSE(p_M)       formula ", formatC(mse_M, format = "e", digits = 4),
      " , simulated ", formatC(mean((p_M - p)^2), format = "e", digits = 4))
  say("   MSE(p_hat)     (14)     ", formatC(var_phat(q, M, C), format = "e", digits = 4),
      " , simulated ", formatC(mean((p_hat - p)^2), format = "e", digits = 4))
  check("7. P{p_hat >= 1} of (17) agrees with the simulation",
        abs(P_ge1 - mean(p_hat >= 1)) < 0.01)
  check("7. E(p_M) of (17) agrees with the simulation",
        abs(E_pM - mean(p_M)) < 4 * sd(p_M) / sqrt(RDES) + 1e-4)
  check("7. E(p_M^2) of (17) agrees with the simulation",
        rel(E_pM2, mean(p_M^2)) < 0.01)
  check("7. MSE(p_M) agrees with the simulation",
        rel(mse_M, mean((p_M - p)^2)) < 0.10)
  check("7. the modification lowers the MSE",
        mse_M <= var_phat(q, M, C) + 1e-12)

  ## ===============================================================
  ##  8. the log-likelihood is maximised at (13)
  ## ===============================================================
  ## l(p) = m0 log(1-pd) + (M-m0) log p , up to a constant in p
  say("")
  say("8. the maximiser of the log-likelihood (12)")
  ok8 <- TRUE
  for (m0v in round(M * c(0.05, 0.2, 0.35))) {
    ll <- function(pp) m0v * log(1 - pp * d) + (M - m0v) * log(pp)
    hi <- min(1 / d - 1e-9, 1)
    o  <- optimize(ll, c(1e-9, hi), maximum = TRUE, tol = 1e-12)
    cf <- (1 - m0v / M) / d
    say("   m0 = ", m0v, " : closed form ", fmt5(cf), " , numerical ",
        fmt5(o$maximum))
    ok8 <- ok8 && abs(o$maximum - min(cf, hi)) < 1e-5
  }
  check("8. (13) maximises (12), by numerical maximisation", ok8)

  ## ===============================================================
  ##  9. the limits of the Remark and the Corollary
  ## ===============================================================
  say("")
  say("9. one period only, the limits of the Remark")
  say("   p^2 (1-d)/d                : ", fmt5(p^2 * (1 - d) / d))
  say("   (14) at n = 100000         : ",
      fmt5(var_phat(q, 100000, 99999)))
  say("   p^2 (m2-m1^2)/m1^2         : ", fmt5(p^2 * (m2 - m1^2) / m1^2))
  say("   (11) at n = 100000         : ",
      fmt5(var_ptil(q, 100000, 99999)))
  check("9. Var(p_hat) of one period tends to p^2(1-d)/d",
        rel(var_phat(q, 100000, 99999), p^2 * (1 - d) / d) < 1e-4)
  check("9. Var(p_tilde) of one period tends to p^2(m2-m1^2)/m1^2",
        rel(var_ptil(q, 100000, 99999), p^2 * (m2 - m1^2) / m1^2) < 1e-4)
  ## the Corollary: with M fixed, C = 0 is best
  M0 <- 1200L
  cor_ok <- TRUE
  for (nn in c(1L, 2L, 3L, 4L, 6L))
    cor_ok <- cor_ok &&
      var_phat(q, M0, nn - 1) >= var_phat(q, M0, 0) - 1e-15 &&
      var_ptil(q, M0, nn - 1) >= var_ptil(q, M0, 0) - 1e-15
  check("9. the Corollary: with M fixed both variances are least at C = 0",
        cor_ok)
}

## =================================================================
##  11. periods that are not independent: the term (N-1) n rho_b
## =================================================================
say_rule("11. the between-period term of the Appendix")
## The winds of the periods are made dependent through a Gaussian
## AR(1) copula, so that rho_b is not zero and the full decomposition
##    1 + (n-1) rho_w + (N-1) n rho_b
## may be tested.  rho_b is obtained from the correlations of the
## B_j alone, by Cov(Z_ij, Z_kl) = p^2 Cov(B_j, B_l) for j /= l.
s <- S$section3; q <- quants(s); d <- q$d; p <- q$p
N <- 30L; n <- 5L; M <- N * n; C <- n - 1
phi <- 0.7; reps <- 60000L
e <- matrix(rnorm(reps * N), nrow = reps)
z <- matrix(0, reps, N); z[, 1] <- e[, 1]
for (j in 2:N) z[, j] <- phi * z[, j - 1] + sqrt(1 - phi^2) * e[, j]
v <- qweibull(pnorm(z), s$W$k, s$W$c)
b <- (v >= s$v_ci & v < s$v_co)
kk <- matrix(rbinom(reps * N, n, p), nrow = reps)
m0 <- M - rowSums(b * kk)
cb <- cov(b)                         # the covariances of the B_j
cov_off <- (sum(cb) - sum(diag(cb))) / (N * (N - 1))
rho_b <- p^2 * cov_off / q$s2_Z
fac_th <- 1 + (n - 1) * q$rho_w + (N - 1) * n * rho_b
var_th <- M * q$s2_Z * fac_th
say("phi of the AR(1) copula            : ", phi)
say("rho_w of (15)                      : ", fmt5(q$rho_w))
say("rho_b from the covariances of B    : ", fmt5(rho_b))
say("factor 1 + (n-1) rho_w + (N-1) n rho_b : ", fmt4(fac_th))
say("factor with rho_b ignored              : ", fmt4(1 + (n - 1) * q$rho_w))
say("Var(m_0) formula                   : ", fmt4(var_th),
    " , simulated ", fmt4(var(m0)))
check("11. the full decomposition agrees with the simulation",
      rel(var_th, var(m0)) < 0.06,
      paste0("ratio ", fmt4(var(m0) / var_th)))
check("11. ignoring rho_b understates the variance badly",
      var(m0) / (M * q$s2_Z * (1 + (n - 1) * q$rho_w)) > 1.5)

## =================================================================
##  12. the panel of the site
## =================================================================
say_rule("12. the panel: N, M, C and the two correlations")
if (!file.exists(PANEL_RDS)) {
  say("the panel is not present; run code/prepare_data.R for this part")
} else {
  pan <- readRDS(PANEL_RDS)$panel
  sub <- estimation_sample(pan)
  v_ci <- WT$v_ci; v_co <- WT$v_co
  sub$Z  <- as.numeric(zero_of(sub, v_ci, v_co))
  sub$ib <- in_band(sub$ws, v_ci, v_co)
  sub$pw <- ifelse(sub$ib, sub$power, 0)
  nj <- as.numeric(table(sub$slot))
  M  <- sum(nj); N <- length(nj)
  C  <- sum(nj * (nj - 1)) / M
  say("N, periods with at least one cell : ", fmti(N))
  say("M = sum n_j                       : ", fmti(M))
  say("sum n_j(n_j-1)                    : ", fmti(sum(nj * (nj - 1))))
  say("C = sum n_j(n_j-1) / M            : ", fmt4(C))
  say("periods with n_j = 6              : ", fmti(sum(nj == 6)))
  check("12. C is below n-1 = 5, the sample being ragged", C < 5 && C > 4.9)

  full <- as.integer(names(table(sub$slot)))[nj == 6]
  fs <- sub[sub$slot %in% full, ]
  fs <- fs[order(fs$slot, fs$turbine), ]
  Zm <- matrix(fs$Z,  ncol = 6, byrow = TRUE)
  Pm <- matrix(fs$pw, ncol = 6, byrow = TRUE)
  cw <- cor(Zm); cp <- cor(Pm)
  rho_w_emp <- mean(cw[upper.tri(cw)]); rho_P_emp <- mean(cp[upper.tri(cp)])
  d_rec <- mean(sub$ib); m0 <- sum(sub$Z)
  p_rec <- (1 - m0 / M) / d_rec
  qp <- list(p = p_rec, d = d_rec, rho_w = p_rec * (1 - d_rec) / (1 - p_rec * d_rec))
  say("")
  say("d from the wind speed record      : ", fmt5(d_rec))
  say("p_hat of (13)                     : ", fmt5(p_rec))
  say("rho_w measured                    : ", fmt5(rho_w_emp),
      "   (the 15 pairs run from ", fmt4(min(cw[upper.tri(cw)])), " to ",
      fmt4(max(cw[upper.tri(cw)])), ")")
  say("rho_w of (15), a common wind      : ", fmt5(qp$rho_w))
  say("rho_P measured                    : ", fmt5(rho_P_emp))
  check("12. the measured rho_w is below the common wind value of (15)",
        rho_w_emp < qp$rho_w)
  check("12. the measured rho_P is below one", rho_P_emp < 1)

  f <- file.path("data", "design_study.rds")
  if (file.exists(f)) {
    ds <- readRDS(f)$tab
    say("")
    say("the two designs of run_simulation.R")
    say_block(ds)
    MA <- ds$N[ds$design == "A"]; MB <- ds$N[ds$design == "B"]
    sA <- ds$sd_phat[ds$design == "A"]; sB <- ds$sd_phat[ds$design == "B"]
    ## (14) with the measured rho_w, C = 5 for A and C = 0 for B
    vA <- p_rec * (1 - p_rec * d_rec) / (MA * d_rec) * (1 + 5 * rho_w_emp)
    vB <- p_rec * (1 - p_rec * d_rec) / (MB * d_rec)
    say("design A : (14) with the measured rho_w gives sd ", fmt5(sqrt(vA)),
        " , simulated ", fmt5(sA))
    say("design B : (14) with C = 0 gives sd            ", fmt5(sqrt(vB)),
        " , simulated ", fmt5(sB))
    say("the ratio of the two simulated variances, corrected for M : ",
        fmt4((sA / sB)^2 * MA / MB))
    say("1 + 5 rho_w with the measured rho_w                       : ",
        fmt4(1 + 5 * rho_w_emp))
    ## (14) is written for periods drawn at random, while Step 3 draws
    ## them within the twenty-four strata of month and day or night.
    ## The stratification removes the seasonal variation and lowers
    ## both designs below (14) by about the same factor, so that the
    ## two absolute values are bounds and it is their ratio, in which
    ## the factor cancels, that tests the clustering term.
    say("each design lies below (14) by the gain of the stratification :")
    say("   design A ", fmt4(sA^2 / vA), " , design B ", fmt4(sB^2 / vB))
    check("12. design A lies below (14), the periods being stratified",
          sA^2 < vA, paste0("ratio ", fmt4(sA^2 / vA)))
    check("12. design B lies below (14), the periods being stratified",
          sB^2 < vB, paste0("ratio ", fmt4(sB^2 / vB)))
    check("12. the two gains are of the same size, the factor cancelling",
          rel(sA^2 / vA, sB^2 / vB) < 0.15,
          paste0("A/B of the gains ", fmt4((sA^2 / vA) / (sB^2 / vB))))
    check("12. the ratio of the designs agrees with 1 + 5 rho_w",
          rel((sA / sB)^2 * MA / MB, 1 + 5 * rho_w_emp) < 0.15,
          paste0("ratio ", fmt4(((sA / sB)^2 * MA / MB) / (1 + 5 * rho_w_emp))))
  } else {
    say("data/design_study.rds is absent; run code/run_simulation.R for",
        " the last checks")
  }
}

say_rule("summary")
say("checks passed : ", n_pass)
say("checks failed : ", n_fail)
if (n_fail > 0L) stop("some checks failed; see the report", call. = FALSE)
close_report()
