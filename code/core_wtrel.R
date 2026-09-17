## =================================================================
##  core_wtrel.R
##
##  Shared core for
##    S. Eryilmaz and C. Kus,
##    "On the estimation of wind turbine reliability from turbine
##     power output data".
##
##  Every other script in code/ sources this file, so that the real
##  data application, the Monte Carlo study and the two formula
##  checks use exactly one implementation of the model and of the
##  estimators.
##
##  Base R only. No contributed package is required.
##
##  Notation is that of the paper:
##    p       long-run probability that the WT is in an operating state
##    n       number of observations of the state of a WT
##    m0      number of observations with zero power output
##    v_ci    cut-in wind speed          v_r   rated wind speed
##    v_co    cut-out wind speed         P_r   nominal (rated) power
##    F, f    cdf and pdf of the wind speed V
##    g       wind speed to power output function, Eq. (2)
##    d       F(v_co) - F(v_ci)
##    mu(p)   mean power output, Eq. (5)
##    p_hat   MLE, Eq. (12)              p_M   modified MLE
##    p_tilde moment estimator, Eq. (8)
##  Two symbols belong to a figure that the code draws and the paper
##  does not carry; they do not appear in the paper and are index free:
##    q(v)    probability that an operating WT produces power at the
##            wind speed v, which the model takes to be one between
##            v_ci and v_co
##    u       wind speed at which q equals one half, that is the speed
##            at which an operating WT actually ceases to produce
## =================================================================

## ---- the folders written by the scripts -------------------------
## Reports go to results/, the figure to figures/ as .eps and as .pdf,
## and the LaTeX of every table and of the new section to
## "Claude outputs/".  Every table is
## in addition printed into the report of the script that produced it,
## so that the manuscript can be brought up to date from results/
## alone.
## =================================================================
##  The folder of one run
## =================================================================
##  Every report, table and figure of a run is written under
##      runs/<date>_<time>/
##  so that two runs never mix.  The folder is made by the first script
##  of the run and its name is kept in runs/LATEST.txt, which the
##  scripts that follow read, so that they all write to the same place.
##  code/new_run.R starts a new one.  The caches of the data stay in
##  data/ and are shared by every run, the reading of the archives being
##  the slow step.
run_root <- local({
  env <- Sys.getenv("WT_RUN", "")
  if (nzchar(env)) return(env)
  latest <- file.path("runs", "LATEST.txt")
  if (file.exists(latest)) {
    p <- trimws(readLines(latest, warn = FALSE))[1]
    if (!is.na(p) && nzchar(p)) return(p)
  }
  p <- file.path("runs", format(Sys.time(), "%Y-%m-%d_%H%M%S"))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  writeLines(p, latest)
  p
})
for (d in c(run_root, file.path(run_root, c("results", "figures", "tables"))))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

res_dir <- file.path(run_root, "results")
fig_dir <- file.path(run_root, "figures")

## =================================================================
##  The site
## =================================================================
## The real data application is carried out on the Kelmarsh wind
## farm.  The same code prepares and analyses the Penmanshiel wind farm
## of the same data set, which serves as a second site; the site is
## named in the environment variable WT_SITE and the characteristics
## below are in every case those published by the manufacturer.
WT_SITES <- list(
  Kelmarsh = list(site = "Kelmarsh", country = "United Kingdom",
                  model = "Senvion MM92", P_r = 2050, v_ci = 3,
                  v_r = 12.5, v_co = 24, years = 2017:2019,
                  doi = "10.5281/zenodo.5841833"),
  ## one complete year, which is what Step 1 asks for; the earlier
  ## archives of this site do not carry every signal that is needed
  Penmanshiel = list(site = "Penmanshiel", country = "United Kingdom",
                  model = "Senvion MM82", P_r = 2050, v_ci = 3.5,
                  v_r = 14.5, v_co = 25, years = 2019L,
                  doi = "10.5281/zenodo.5946808"))
SITE <- Sys.getenv("WT_SITE", "Kelmarsh")
if (!SITE %in% names(WT_SITES))
  stop("WT_SITE must be one of ", paste(names(WT_SITES), collapse = " or "),
       ", not ", SITE, call. = FALSE)
WT <- WT_SITES[[SITE]]

## Every file that a script writes carries the name of the site, except
## at Kelmarsh, whose file names are those of the paper.
tag <- function(x) if (SITE == "Kelmarsh") x else
  paste0(x, "_", tolower(SITE))
## The folders in which the archives of the site are looked for.  They are
## large and are often kept once for several working folders, so the
## folder data/ of the project is searched first, then the folder data/
## beside it, and then whatever WT_DATA names.  Everything a run writes
## goes to data/ of the project whatever the archives are read from.
data_dirs <- function() {
  d <- c(Sys.getenv("WT_DATA", ""), "data", file.path("..", "data"))
  d <- d[nzchar(d)]
  d[dir.exists(d)]
}
find_data_file <- function(pattern) {
  for (d in data_dirs()) {
    hit <- list.files(d, pattern, full.names = TRUE)
    if (length(hit)) return(hit)
  }
  character(0)
}

## Everything a run caches is written here, whatever folder the archives
## were read from.
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

PANEL_RDS  <- file.path("data", paste0(tolower(SITE), "_panel.rds"))
STATUS_RDS <- file.path("data", paste0(tolower(SITE), "_status.rds"))
tex_dir <- file.path(run_root, "tables")

## =================================================================
##  The setting of the numerical illustrations of the paper
## =================================================================
## Section 3 of the paper assumes a Weibull wind speed with k = 2 and
## alpha = 5 and a WT with v_ci = 3, v_r = 11, v_co = 20 m/s and
## P_r = 5.5, and its Table 1 and Figures 1 and 2 are computed from
## these.  The values are held here, so that the script that produces
## them and the two formula checks read one and the same setting.
PAPER <- list(k = 2, alpha = 5, v_ci = 3, v_r = 11, v_co = 20, P_r = 5.5,
              p_tab = c(.95, .96, .97, .98),   # the rows of Table 1
              n_tab = c(10L, 15L, 20L, 30L),   # its blocks, and Figure 1
              n_fig = 10:30)                   # the abscissa of Figure 2
## The 32 values printed in Table 1 of the paper, read off the paper and
## used to verify what the code produces.
PAPER_TABLE1 <- c(0.0257, 0.0459, 0.0251, 0.0454, 0.0245, 0.0449, 0.0240,
                  0.0444, 0.0170, 0.0306, 0.0166, 0.0303, 0.0161, 0.0300,
                  0.0158, 0.0296, 0.0131, 0.0230, 0.0126, 0.0227, 0.0122,
                  0.0225, 0.0119, 0.0222, 0.0090, 0.0153, 0.0085, 0.0151,
                  0.0082, 0.0150, 0.0079, 0.0148)

## ---- scaling of the variances in the tables ---------------------
var_scale     <- 1e4
var_scale_tex <- "10^{4}"

## ---- formatting -------------------------------------------------
fmt4 <- function(x) formatC(as.numeric(x), format = "f", digits = 4)
fmt5 <- function(x) formatC(as.numeric(x), format = "f", digits = 5)
fmt7 <- function(x) formatC(as.numeric(x), format = "f", digits = 7)
fmti <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")

## =================================================================
##  1. The model
## =================================================================

## Eq. (2): the power output of a WT as a function of the wind speed
g_wt <- function(v, v_ci, v_r, v_co, P_r) {
  out <- numeric(length(v))
  mid <- which(v >= v_ci & v < v_r)
  out[mid] <- P_r * (v[mid]^3 - v_ci^3) / (v_r^3 - v_ci^3)
  out[which(v >= v_r & v < v_co)] <- P_r
  out
}

## A wind speed distribution is carried as a list with the cdf F and
## the pdf f, so that any continuous distribution may be supplied.
wind_weibull <- function(k, c_scale)
  list(F = function(v) pweibull(v, k, c_scale),
       f = function(v) dweibull(v, k, c_scale),
       k = k, c = c_scale, name = "Weibull")

## H1 of Section 2 and the two atoms of the distribution (3) - (4)
H1_wt <- function(x, W, v_ci, v_r, v_co, P_r)
  1 - W$F(v_co) + W$F(((x / P_r) * (v_r^3 - v_ci^3) + v_ci^3)^(1 / 3))

d_wt      <- function(W, v_ci, v_co) W$F(v_co) - W$F(v_ci)
P_zero_wt <- function(p, W, v_ci, v_co) 1 - p * d_wt(W, v_ci, v_co)
P_rated_wt<- function(p, W, v_r, v_co) p * (W$F(v_co) - W$F(v_r))

## E[g(V)^j], j = 1, 2.  The integral over (v_ci, v_r) is taken in the
## wind speed scale, which is equivalent to the integral of x^j h1(x)
## over (0, P_r) appearing in Eqs. (5) and (10).
Eg_wt <- function(W, v_ci, v_r, v_co, P_r, j = 1) {
  int <- integrate(function(v) g_wt(v, v_ci, v_r, v_co, P_r)^j * W$f(v),
                   lower = v_ci, upper = v_r, rel.tol = 1e-10)$value
  int + P_r^j * (W$F(v_co) - W$F(v_r))
}

mu_wt  <- function(p, W, v_ci, v_r, v_co, P_r) p * Eg_wt(W, v_ci, v_r, v_co, P_r, 1)
mu2_wt <- function(p, W, v_ci, v_r, v_co, P_r) p * Eg_wt(W, v_ci, v_r, v_co, P_r, 2)

## =================================================================
##  2. The estimators and their mean squared errors
## =================================================================

## Eq. (12) and its modification
p_hat_wt   <- function(m0, n, d) (1 - m0 / n) / d
p_M_wt     <- function(m0, n, d) pmin(p_hat_wt(m0, n, d), 1)
## Eq. (8)
p_tilde_wt <- function(power, W, v_ci, v_r, v_co, P_r)
  mean(power) / Eg_wt(W, v_ci, v_r, v_co, P_r, 1)

## Eq. (9)
mse_moment_wt <- function(p, n, W, v_ci, v_r, v_co, P_r) {
  m1 <- Eg_wt(W, v_ci, v_r, v_co, P_r, 1)
  m2 <- Eg_wt(W, v_ci, v_r, v_co, P_r, 2)
  (p * m2 - p^2 * m1^2) / (n * m1^2)
}

## Eq. (13)
mse_mle_wt <- function(p, n, W, v_ci, v_co) {
  d <- d_wt(W, v_ci, v_co)
  p * (1 - p * d) / (n * d)
}

## MSE of the modified MLE, from the two moments given in Section 2.2
mse_mmle_wt <- function(p, n, W, v_ci, v_co) {
  d    <- d_wt(W, v_ci, v_co)
  k    <- 0:n
  pk   <- dbinom(k, n, 1 - p * d)
  thr  <- floor(n * (1 - d))
  low  <- k > thr
  Pge1 <- sum(pk[!low])
  e1 <- sum((1 - k[low] / n) / d * pk[low]) + Pge1
  e2 <- sum(((1 - k[low] / n) / d)^2 * pk[low]) + Pge1
  (e1 - p)^2 + (e2 - e1^2)
}

## Standard error implied by Eq. (13)
se_mle_wt <- function(p, n, d) sqrt(pmax(p, 0) * (1 - pmax(p, 0) * d) / (n * d))

## =================================================================
##  2b. The period sample of Section 2.2: the clustering constant,
##      the two correlations, and the maximum likelihood estimator
##      under a wind speed common to the WTs of one period
## =================================================================
##  The sample is N periods; n_j of the n WTs remain in period j after
##  the filtering of Step 2, so that
##      M = sum_i N_i = sum_j n_j        and       C = sum_j n_j(n_j-1)/M .
##  Two WTs of one period see the same wind and are therefore dependent;
##  two WTs of two periods are independent.  Everything below is exact.

## C of Section 2.2, from the n_j
C_design <- function(nj) sum(nj * (nj - 1)) / sum(nj)

## The two correlations of Section 2.2
rho_w_wt <- function(p, d)          p * (1 - d) / (1 - p * d)
rho_P_wt <- function(p, m1, m2)     p * (m2 - m1^2) / (m2 - p * m1^2)

## The MSE of the moment estimator and of the MLE over such a sample
mse_moment_C <- function(p, M, C, m1, m2)
  (p * m2 - p^2 * m1^2) / (M * m1^2) * (1 + rho_P_wt(p, m1, m2) * C)
mse_mle_C <- function(p, M, C, d)
  p * (1 - p * d) / (M * d) * (1 + rho_w_wt(p, d) * C)

## The distribution of the number of zeros of one period, and of m_0
pmf_period <- function(p, n, d)
  d * dbinom(0:n, n, 1 - p) + (1 - d) * (0:n == n)

pmf_m0_wt <- function(p, nj, d) {
  out <- 1
  for (n in nj) out <- convolve(out, rev(pmf_period(p, n, d)), type = "open")
  out[out < 0] <- 0
  out / sum(out)
}

## The MSE of the modified MLE over the period sample
mse_mmle_C <- function(p, nj, d) {
  M  <- sum(nj)
  pk <- pmf_m0_wt(p, nj, d)
  k  <- 0:M
  ph <- (1 - k / M) / d
  ge <- ph >= 1
  e1 <- sum(pk[!ge] * ph[!ge]) + sum(pk[ge])
  e2 <- sum(pk[!ge] * ph[!ge]^2) + sum(pk[ge])
  (e1 - p)^2 + e2 - e1^2
}

## ---- the likelihood of the period sample ------------------------
##  A period in which every retained WT is idle is uninformative about
##  the reason; a period in which at least one WT produces power tells
##  that the wind was between the cut-in and the cut-out speeds, so that
##  the zeros of that period are failures.  With A the set of the first
##  kind, S the zeros and R the positive cells of the periods outside A,
##      l(p) = R log p + S log(1-p)
##             + sum_{j in A} log{ (1-d) + d (1-p)^{n_j} } .
##  The statistics the likelihood needs: R positive cells and S zeros
##  among the periods that are not idle throughout, and the sizes nA of
##  the periods that are.
stats_common <- function(m0j, nj) {
  A <- m0j >= nj
  list(R = sum(nj[!A] - m0j[!A]), S = sum(m0j[!A]), nA = nj[A])
}

loglik_common <- function(p, m0j, nj, d) {
  st <- stats_common(m0j, nj)
  loglik_stats(p, st$R, st$S, st$nA, d)
}

loglik_stats <- function(p, R, S, nA, d) {
  v <- numeric(length(p))
  if (R > 0) v <- v + R * log(p)
  if (S > 0) v <- v + S * log(1 - p)
  if (length(nA)) {
    tb <- table(nA)                       # the distinct period sizes only
    nn <- as.numeric(names(tb)); cc <- as.numeric(tb)
    u  <- 1 - p
    v  <- v + as.vector(cc %*% log((1 - d) + d * outer(nn, u,
                                     function(a, b) b^a)))
  }
  v
}

##  The estimator itself.  The likelihood is a likelihood, so that its
##  maximiser lies in [0,1] by construction and no modification of the
##  kind of Section 2.4 is called for.  A grid locates the maximum and
##  optimize refines it.
p_ml_stats <- function(R, S, nA, d, grid = 401L) {
  if (R == 0) return(0)                      # every period idle throughout
  g  <- seq(1e-10, 1 - 1e-10, length.out = grid)
  lg <- loglik_stats(g, R, S, nA, d)
  i  <- which.max(lg)
  if (i == grid) return(1)
  lo <- g[max(i - 1L, 1L)]; hi <- g[min(i + 1L, grid)]
  optimize(function(q) loglik_stats(q, R, S, nA, d), c(lo, hi),
           maximum = TRUE, tol = 1e-12)$maximum
}

p_ml_wt <- function(m0j, nj, d, grid = 401L) {
  st <- stats_common(m0j, nj)
  p_ml_stats(st$R, st$S, st$nA, d, grid)
}

##  The Fisher information of one period of n WTs, in closed form, and
##  the standard error it implies over the whole sample.
info_period_wt <- function(p, n, d)
  n * d / (p * (1 - p)) - n^2 * d * (1 - p)^(n - 2) +
  n^2 * d^2 * (1 - p)^(2 * n - 2) / (1 - d + d * (1 - p)^n)

se_ml_wt <- function(p, nj, d) 1 / sqrt(sum(info_period_wt(p, nj, d)))

##  The exact MSE of the estimator when every period keeps n WTs: the
##  N periods fall into the n+1 cells of pmf_period, and every one of
##  the compositions of N into those cells is enumerated.
compositions_wt <- function(N, m) {
  if (m == 1L) return(matrix(N, nrow = 1L))
  out <- list()
  for (k in 0:N) {
    rest <- compositions_wt(N - k, m - 1L)
    out[[length(out) + 1L]] <- cbind(k, rest, deparse.level = 0)
  }
  do.call(rbind, out)
}

mse_ml_C <- function(p, N, n, d) {
  N <- as.integer(N); n <- as.integer(n)
  cm  <- compositions_wt(N, n + 1L)
  pr  <- pmf_period(p, n, d)
  lw  <- apply(cm, 1L, function(r) lfactorial(N) - sum(lfactorial(r)) +
                                   sum(ifelse(r > 0, r * log(pr), 0)))
  w   <- exp(lw)
  a   <- cm[, n + 1L]
  S   <- if (n == 1L) rep(0, nrow(cm)) else
           as.vector(cm[, seq_len(n), drop = FALSE] %*% (0:(n - 1L)))
  key <- paste(a, S)
  u   <- !duplicated(key)
  lut <- vapply(which(u), function(i)
           p_ml_stats(n * (N - a[i]) - S[i], S[i], rep(n, a[i]), d), numeric(1))
  est <- lut[match(key, key[u])]
  c(mse = sum(w * (est - p)^2), bias = sum(w * est) - p,
    P0 = sum(w[est <= 0]), P1 = sum(w[est >= 1]))
}

## =================================================================
##  3. The cut-in speed in the field and the fit of the Weibull F
## =================================================================

## q(v), the probability that an operating WT produces power at the wind
## speed v, and u, the speed at which q equals one half.  The records
## are grouped into narrow wind speed classes of width `width`, classes
## with fewer than `min_count` records being discarded, and u is read
## off by linear interpolation.  Only records of an operating WT are to
## be passed, since a WT that is down produces nothing at any speed.
q_positive <- function(v, power, lo = 1.5, hi = 4.5, width = 0.05,
                       min_count = 200L) {
  b    <- round(v / width) * width
  keep <- b >= lo & b <= hi
  qq   <- tapply(power[keep] > 0, b[keep], mean)
  nn   <- tapply(power[keep] > 0, b[keep], length)
  ok   <- nn >= min_count
  cv   <- data.frame(v = as.numeric(names(qq))[ok], q = as.numeric(qq)[ok],
                     n = as.numeric(nn)[ok])
  cv   <- cv[order(cv$v), ]
  list(u = as.numeric(approx(cv$q, cv$v, xout = 0.5, ties = mean)$y),
       curve = cv)
}

## Weibull fit by maximum likelihood, base R, with the moment estimates
## as starting values
fit_weibull_ml <- function(v) {
  v  <- v[is.finite(v) & v > 0]
  st <- fit_weibull_mom(v)
  ## optim may probe extreme parameters during the line search; the
  ## warning that dweibull then emits is harmless and is suppressed
  nll <- function(th) {
    val <- suppressWarnings(-sum(dweibull(v, exp(th[1]), exp(th[2]), log = TRUE)))
    if (!is.finite(val)) 1e100 else val
  }
  o  <- optim(log(c(st$k, st$c)), nll, method = "BFGS",
              control = list(reltol = 1e-14, maxit = 500))
  c(k = exp(o$par[1]), c = exp(o$par[2]), convergence = o$convergence)
}

## Starting values for the maximum likelihood fit, from the usual
## moment approximation; not reported anywhere
fit_weibull_mom <- function(v) {
  v <- v[is.finite(v) & v > 0]
  m <- mean(v); s <- sd(v)
  k <- (s / m)^(-1.086)
  list(k = k, c = m / gamma(1 + 1 / k))
}

## -----------------------------------------------------------------
##  The choice of the wind speed distribution
## -----------------------------------------------------------------
## The paper takes F to be Weibull.  That choice is not assumed here
## but verified: the families below, all of them in common use for
## wind speed, are fitted to the same wind speed record by maximum
## likelihood and compared by the Akaike information criterion.  Every
## family is supported on the positive half line and is parametrised
## on the log scale, so that the maximisation is unconstrained.
WIND_FAMILIES <- list(
  list(name = "Weibull", npar = 2L,
       start = function(x) { s <- fit_weibull_mom(x); log(c(s$k, s$c)) },
       dens  = function(x, th) dweibull(x, exp(th[1]), exp(th[2]), log = TRUE),
       cdf   = function(q, th) pweibull(q, exp(th[1]), exp(th[2])),
       pars  = function(th) c(k = exp(th[1]), alpha = exp(th[2]))),
  list(name = "Gamma", npar = 2L,
       start = function(x) log(c(mean(x)^2 / var(x), mean(x) / var(x))),
       dens  = function(x, th) dgamma(x, exp(th[1]), exp(th[2]), log = TRUE),
       cdf   = function(q, th) pgamma(q, exp(th[1]), exp(th[2])),
       pars  = function(th) c(shape = exp(th[1]), rate = exp(th[2]))),
  list(name = "Lognormal", npar = 2L,
       start = function(x) c(mean(log(x)), log(sd(log(x)))),
       dens  = function(x, th) dlnorm(x, th[1], exp(th[2]), log = TRUE),
       cdf   = function(q, th) plnorm(q, th[1], exp(th[2])),
       pars  = function(th) c(meanlog = th[1], sdlog = exp(th[2]))),
  list(name = "Log-logistic", npar = 2L,
       start = function(x) log(c(median(x), 2)),
       dens  = function(x, th) {
         a <- exp(th[1]); b <- exp(th[2])
         log(b / a) + (b - 1) * log(x / a) - 2 * log1p((x / a)^b)
       },
       cdf   = function(q, th) {
         a <- exp(th[1]); b <- exp(th[2])
         1 / (1 + (q / a)^(-b))
       },
       pars  = function(th) c(alpha = exp(th[1]), beta = exp(th[2]))),
  list(name = "Rayleigh", npar = 1L,
       start = function(x) log(sqrt(mean(x^2) / 2)),
       dens  = function(x, th) dweibull(x, 2, sqrt(2) * exp(th[1]), log = TRUE),
       cdf   = function(q, th) pweibull(q, 2, sqrt(2) * exp(th[1])),
       pars  = function(th) c(sigma = exp(th[1]))))

fit_family_ml <- function(f, x) {
  nll <- function(th) {
    val <- suppressWarnings(-sum(f$dens(x, th)))
    if (!is.finite(val)) 1e100 else val
  }
  o <- optim(f$start(x), nll,
             method = if (f$npar == 1L) "Brent" else "BFGS",
             lower  = if (f$npar == 1L) -10 else -Inf,
             upper  = if (f$npar == 1L)  10 else  Inf,
             control = list(reltol = 1e-14, maxit = 1000))
  list(par = o$par, nll = o$value, conv = o$convergence)
}

## The comparison itself.  A wind speed recorded as exactly zero, if
## there is any, is carried as an atom at zero and is common to every
## family, F(v) = pi0 + (1 - pi0) G(v); it lies below the cut-in speed
## in any case.  The table reports, beside the criterion, the only two
## features of F that the estimators use, F(v_ci) and d.
fit_wind_families <- function(v_all, v_ci, v_co) {
  v   <- v_all[is.finite(v_all) & v_all > 0]
  pi0 <- mean(v_all == 0)
  fits <- list()
  tab  <- data.frame()
  for (f in WIND_FAMILIES) {
    fit <- fit_family_ml(f, v)
    fits[[f$name]] <- list(f = f, fit = fit)
    Fci <- pi0 + (1 - pi0) * f$cdf(v_ci, fit$par)
    Fco <- pi0 + (1 - pi0) * f$cdf(v_co, fit$par)
    tab <- rbind(tab, data.frame(
      family = f$name, npar = f$npar, logLik = -fit$nll,
      AIC = 2 * f$npar + 2 * fit$nll, F_vci = Fci, d = Fco - Fci,
      conv = fit$conv, stringsAsFactors = FALSE))
  }
  tab      <- tab[order(tab$AIC), ]
  tab$dAIC <- tab$AIC - tab$AIC[1]
  row.names(tab) <- NULL
  list(tab = tab, fits = fits, pi0 = pi0, n = length(v),
       best = tab$family[1])
}

## =================================================================
##  3b. The archives of the site
## =================================================================
## Returns the folder that holds the csv files of one year, unpacking
## the archive when necessary.  The second element of the returned list
## says whether that folder is scratch and may be removed afterwards.
year_folder <- function(y, keep_raw = FALSE) {
  kept <- file.path("data", "raw", y)
  if (dir.exists(kept) && length(list.files(kept, "^Turbine_Data_")))
    return(list(dir = kept, scratch = FALSE))
  z <- find_data_file(sprintf("^%s_SCADA_%d_.*\\.zip$", SITE, y))
  if (!length(z))
    stop("\n\n  The data of ", y, " was not found: there is neither a folder\n",
         "  data/raw/", y, "/ nor an archive ", SITE, "_SCADA_", y,
         "_*.zip in data/.\n\n",
         "  Open  https://doi.org/", WT$doi, "  in a browser and download\n",
         "  the archive of ", y, ", named ", SITE, "_SCADA_", y,
         "_<number>.zip , into\n",
         "  the folder data/ of this project, beside code/, with its name\n",
         "  unchanged.  A folder data/ beside this one is searched as well,\n",
         "  and so is the folder named in the variable WT_DATA, so that the\n",
         "  archives need not be copied twice.  Then run this script again.  The paper uses the three\n",
         "  archives of 2017, 2018 and 2019 and nothing else.  The Data\n",
         "  section of README.md gives their exact names and sizes.\n",
         call. = FALSE)
  target <- if (keep_raw) kept else
    file.path(tempdir(), paste0(tolower(SITE), "_raw_", y))
  dir.create(target, showWarnings = FALSE, recursive = TRUE)
  ## the target is a folder under the session temporary directory; its
  ## path is not printed, so that the report does not carry the name of
  ## the machine on which it was produced
  ## A year may come in more than one archive, as it does at Penmanshiel.
  ## unzip() of base R reports a failure by a warning only, and a silent
  ## failure would leave the panel incomplete, so the table of contents
  ## of every archive is compared with what reached the folder.
  for (zz in z) {
    say("unpacking ", basename(zz))
    toc  <- unzip(zz, list = TRUE)$Name
    toc  <- toc[grepl("\\.csv$", toc, ignore.case = TRUE)]
    w <- NULL
    withCallingHandlers(unzip(zz, exdir = target),
                        warning = function(cond) {
                          w <<- c(w, conditionMessage(cond))
                          invokeRestart("muffleWarning")
                        })
    miss <- toc[!file.exists(file.path(target, toc))]
    if (length(miss))
      stop("\n\n  ", length(miss), " of the ", length(toc),
           " csv files of\n  ", basename(zz),
           " were not extracted, the first being\n    ", miss[1],
           "\n  The archive may be incomplete or damaged; download it again",
           " and\n  run this script anew.\n", call. = FALSE)
    if (length(w))
      say("  note: unzip reported ", length(w), " warning(s), but every one",
          " of the ", length(toc), " csv files of the archive is present")
  }
  list(dir = target, scratch = !keep_raw)
}

## =================================================================
##  4. Step 3 of the protocol: the sampling design
## =================================================================
## The 24 strata defined by the twelve months and by day and night
## receive the same number of instants, and two selected instants are
## at least gap slots apart, one slot being ten minutes.

make_strata <- function(ts) {
  mo <- as.integer(format(ts, "%m"))
  hr <- as.integer(format(ts, "%H"))
  paste0(mo, ifelse(hr >= 6 & hr < 18, "D", "N"))
}

draw_slots <- function(slot, stratum, per, gap, n_slot) {
  blocked <- logical(n_slot + 2L * gap + 2L)
  off <- gap + 1L
  acc <- integer(0)
  for (s in sample(unique(stratum))) {
    cand <- slot[stratum == s]
    cand <- cand[sample.int(length(cand))]
    got <- 0L
    for (cc in cand) {
      if (got >= per) break
      if (!blocked[cc + off]) {
        acc <- c(acc, cc)
        blocked[(cc + off - gap):(cc + off + gap)] <- TRUE
        got <- got + 1L
      }
    }
  }
  sort(acc)
}

## =================================================================
##  5. Figures: one device helper, so that every figure has the same
##     margins and fonts, and is written both as .eps and as .pdf
## =================================================================
## The rendering is done in a temporary folder and the result copied,
## because a synchronised drive or a non-ASCII path makes postscript()
## and pdf() fail on some systems whereas file.copy() handles both.
fig_device <- function(name, expr, width = 5.2, height = 3.6,
                       mfrow = c(1, 1)) {
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  draw <- function() {
    par(mfrow = mfrow, mar = c(4.2, 4.4, 1.0, 0.8), mgp = c(2.6, 0.8, 0),
        cex = 0.85, cex.main = 1.0, font.main = 1, las = 0)
    eval(expr)
  }
  tmp_eps <- file.path(tempdir(), paste0(name, ".eps"))
  postscript(tmp_eps, width = width, height = height, horizontal = FALSE,
             onefile = FALSE, paper = "special", family = "Helvetica")
  draw(); dev.off()
  tmp_pdf <- file.path(tempdir(), paste0(name, ".pdf"))
  pdf(tmp_pdf, width = width, height = height, family = "Helvetica")
  draw(); dev.off()
  file.copy(tmp_eps, file.path(fig_dir, paste0(name, ".eps")), overwrite = TRUE)
  file.copy(tmp_pdf, file.path(fig_dir, paste0(name, ".pdf")), overwrite = TRUE)
  invisible(file.path(fig_dir, paste0(name, c(".eps", ".pdf"))))
}

## =================================================================
##  6. The single LaTeX table emitter
## =================================================================
## body   : character matrix or data frame, already formatted
## header : character vector, one entry per header row, cells separated
##          by " & " already, or a plain vector of column names
## rules  : row numbers after which a \hline is placed
tex_table <- function(body, colspec, header, caption, label,
                      rules = integer(0), file = NULL, note = NULL,
                      label_row = NULL) {
  body <- as.matrix(body)
  ## an optional full-width label row inserted before body row `at`
  if (!is.null(label_row)) {
    lab <- c(label_row$text, rep(NA_character_, ncol(body) - 1L))
    body <- rbind(body[seq_len(label_row$at - 1L), , drop = FALSE], lab,
                  body[seq(label_row$at, nrow(body)), , drop = FALSE])
  }
  ## A cell equal to NA is dropped, so that a row consisting of a single
  ## \multicolumn spanning the whole width may be placed inside the body.
  rows <- apply(body, 1, function(r) paste(r[!is.na(r)], collapse = " & "))
  out <- c("\\begin{table}[htbp]",
           "\\centering",
           paste0("\\caption{", caption, "}\\label{", label, "}"),
           "\\footnotesize",
           "\\setlength{\\tabcolsep}{2pt}",
           paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}",
                  colspec, "}"),
           "\\hline")
  for (h in header) out <- c(out, paste0(h, " \\\\"))
  out <- c(out, "\\hline")
  for (i in seq_along(rows)) {
    out <- c(out, paste0(rows[i], " \\\\"))
    if (i %in% rules) out <- c(out, "\\hline")
  }
  out <- c(out, "\\hline", "\\end{tabular*}")
  ## The note stands below the closing rule of the table and outside the
  ## tabular.  Inside it, a \\multicolumn of a fixed width would fix the
  ## width of the widest row, and the columns of the body, being narrower,
  ## would then be spread by \\extracolsep over that width rather than
  ## over the width of the text, which throws the last column of the body
  ## against the right margin.
  if (!is.null(note))
    out <- c(out, "", paste0("\\noindent\\parbox{\\textwidth}{\\footnotesize\\raggedright ",
                             note, "\\par}"))
  out <- c(out, "\\end{table}")
  if (!is.null(file)) {
    if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE)
    writeLines(out, file)
  }
  invisible(out)
}

## Writes the table to "Claude outputs/" and prints the identical
## LaTeX into the report, so that the manuscript row and the report row
## can be compared character by character.
emit_table <- function(name, ...) {
  file <- file.path(tex_dir, paste0(name, ".tex"))
  lines <- tex_table(..., file = file)
  say("")
  say("LaTeX of this table, also written to ", file)
  say(strrep("=", 70))
  for (l in lines) say(l)
  say(strrep("=", 70))
  invisible(lines)
}

## =================================================================
##  The sample of the application, defined once
## =================================================================
## Every script that estimates p works on the same records and counts
## m_0 in the same way.  The two rules are here and nowhere else.
##
##   estimation_sample()  the records that Steps 1 and 2 leave, all of
##                        them; no instant is required to carry every
##                        WT, since such a requirement would select the
##                        instants at which no WT has a gap in its
##                        signals.  The replication of Step 3 is the one
##                        exception and says so, the design of Step 4
##                        asking for instants with every WT present.
##
##   zero_of()            the zero-power indicator under the assumption
##                        of Eq. (2): a record whose wind speed lies
##                        outside [v_ci, v_co) counts as zero-power
##                        whatever its measured power, and a measured
##                        power that is not positive counts as zero.
estimation_sample <- function(pan) pan[pan$retained, ]

in_band <- function(v, v_ci, v_co) v >= v_ci & v < v_co

zero_of <- function(z, v_ci = WT$v_ci, v_co = WT$v_co)
  z$zero | !in_band(z$ws, v_ci, v_co)

## -----------------------------------------------------------------
##  Splicing a table into the section
## -----------------------------------------------------------------
## The section is one file, tables included, as a journal asks.  The
## table is nevertheless never typed by hand: emit_table() writes it,
## and splice_table() puts that very text into the section between two
## markers, so that the file carries the table and the code remains its
## only author.  The markers stay in the file and mark the block that
## the next run replaces.
splice_table <- function(name, file = file.path(tex_dir,
                                                "section_real_data.tex")) {
  src <- file.path(tex_dir, paste0(name, ".tex"))
  if (!file.exists(file) || !file.exists(src)) return(invisible(FALSE))
  txt <- readLines(file, warn = FALSE)
  a <- grep(paste0("^% <<< ", name, "$"), txt)
  b <- grep(paste0("^% >>> ", name, "$"), txt)
  if (length(a) != 1L || length(b) != 1L || b <= a) {
    say("note: the markers of ", name, " were not found in ",
        basename(file), "; the table was not spliced")
    return(invisible(FALSE))
  }
  out <- c(txt[seq_len(a)], readLines(src, warn = FALSE), txt[b:length(txt)])
  writeLines(out, file)
  say("the table ", name, " was spliced into ", basename(file))
  invisible(TRUE)
}

## =================================================================
##  7. Guard: a missing input stops the script with an instruction
## =================================================================
need_file <- function(path, produced_by) {
  if (file.exists(path)) return(invisible(TRUE))
  stop("\n\n  The file  ", path, "  is missing.\n",
       "  It is produced by  ", produced_by, "\n",
       "  Run that script first, then run this one again.\n",
       call. = FALSE)
}

## The columns that data/kelmarsh_panel.rds must carry.  A panel built
## by an earlier version of prepare_data.R lacks some of them, and every
## script then stops with an instruction instead of failing obscurely.
PANEL_VERSION <- 4L
PANEL_COLUMNS <- c("ts", "ws", "power", "power_max", "e_exp", "lp_curt",
                   "av_sys", "av_b32", "turbine", "slot",
                   "excl_event", "ev_fo", "ev_ts", "excl_curt", "excluded",
                   "observed", "retained", "zero", "year", "full_instant")
need_panel <- function(D) {
  ver <- if (is.null(D$panel_version)) 0L else D$panel_version
  miss <- setdiff(PANEL_COLUMNS, names(D$panel))
  if (ver < PANEL_VERSION || length(miss))
    stop("\n\n  data/kelmarsh_panel.rds was built by an earlier version of",
         " code/prepare_data.R\n  (panel version ", ver, ", this code needs ",
         PANEL_VERSION, if (length(miss)) paste0("; missing columns: ",
         paste(miss, collapse = ", ")) else "", ").\n",
         "  Run  source(\"code/prepare_data.R\")  once more, then run this",
         " script again.\n", call. = FALSE)
  invisible(TRUE)
}

## =================================================================
##  8. Reporting
## =================================================================
new_report <- function(path) {
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  con <- file(path, "w")
  assign(".wt_report_con", con, envir = .GlobalEnv)
  invisible(con)
}

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  if (exists(".wt_report_con", envir = .GlobalEnv))
    cat(txt, "\n", sep = "", file = get(".wt_report_con", envir = .GlobalEnv))
}

say_block <- function(x) {
  txt <- paste(capture.output(print(x)), collapse = "\n")
  say(txt)
}

close_report <- function() {
  if (exists(".wt_report_con", envir = .GlobalEnv)) {
    close(get(".wt_report_con", envir = .GlobalEnv))
    rm(".wt_report_con", envir = .GlobalEnv)
  }
}

say_rule <- function(title = NULL) {
  say(strrep("-", 70))
  if (!is.null(title)) { say(title); say(strrep("-", 70)) }
}
