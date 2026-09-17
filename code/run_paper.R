## =================================================================
##  run_paper.R
##
##  The numerical illustrations of the paper: Table 1 and Figures 1
##  and 2 of Section 3, computed from the exact expressions of the
##  paper and from nothing else.
##
##    Table 1   MSE(p_M) and MSE(p_hat) for p = .95 to .98 and
##              n = 10, 15, 20, 30
##    Figure 1  the MSE of p_hat, Eq. (13), and of p_tilde, Eq. (9),
##              as functions of p, one panel per sample size
##    Figure 2  the same two, as functions of n, one panel per p
##
##  The setting is the paper's own, held in PAPER of the core: a
##  Weibull wind speed with k = 2 and alpha = 5, and a WT with
##  v_ci = 3, v_r = 11, v_co = 20 m/s and P_r = 5.5.  Nothing is
##  simulated here: every value is exact.  The 32 entries of Table 1
##  are in addition compared with the values printed in the paper, so
##  that the script fails if the code and the paper ever part company.
##
##  Output : results/run_paper_report.txt
##           "Claude outputs"/table1_mse.tex
##           figures/figure1_mse_p.{eps,pdf}
##           figures/figure2_mse_n.{eps,pdf}
##
##  Base R only.  Run time is a few seconds.
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
## data/ is created as well: this script needs no data, and in a fresh
## clone of the repository the folder does not exist yet
for (dd in c(res_dir, fig_dir, tex_dir, "data"))
  if (!dir.exists(dd)) dir.create(dd, recursive = TRUE)

new_report(file.path(res_dir, "run_paper_report.txt"))
say("run_paper.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

W    <- wind_weibull(PAPER$k, PAPER$alpha)
v_ci <- PAPER$v_ci; v_r <- PAPER$v_r; v_co <- PAPER$v_co; P_r <- PAPER$P_r

say_rule("the setting of Section 3 of the paper")
say("wind speed distribution : Weibull, k = ", PAPER$k,
    " , alpha = ", PAPER$alpha)
say("WT                      : v_ci = ", v_ci, " , v_r = ", v_r,
    " , v_co = ", v_co, " m/s , P_r = ", P_r)
say("d = F(v_co) - F(v_ci)   : ", fmt5(d_wt(W, v_ci, v_co)))
say("mu(1) of Eq. (5)        : ", fmt5(Eg_wt(W, v_ci, v_r, v_co, P_r, 1)))

## =================================================================
##  Table 1: the MSE of the MLE and of the modified MLE
## =================================================================
## Eq. (13) for p_hat, and the two moments of Section 2.2 for p_M.
say_rule("Table 1 of the paper")
tab <- expand.grid(p = PAPER$p_tab, n = PAPER$n_tab)
tab$mse_pM   <- mapply(function(p, n) mse_mmle_wt(p, n, W, v_ci, v_co),
                       tab$p, tab$n)
tab$mse_phat <- mapply(function(p, n) mse_mle_wt(p, n, W, v_ci, v_co),
                       tab$p, tab$n)
say_block(data.frame(p = tab$p, n = tab$n,
                     MSE_pM = fmt4(tab$mse_pM),
                     MSE_phat = fmt4(tab$mse_phat)))

## the same 32 values as the paper prints, in the order of its table,
## set beside them so that the agreement can be read off entry by entry
ours  <- as.numeric(t(cbind(round(tab$mse_pM, 4), round(tab$mse_phat, 4))))
prnt  <- matrix(PAPER_TABLE1, ncol = 2, byrow = TRUE)
say("")
say("the printed table of the paper against this computation")
say_block(data.frame(
  p = formatC(tab$p, format = "f", digits = 2),
  n = formatC(tab$n, format = "d"),
  MSE_pM_paper = formatC(prnt[, 1], format = "f", digits = 4),
  MSE_pM_code  = formatC(tab$mse_pM, format = "f", digits = 6),
  MSE_phat_paper = formatC(prnt[, 2], format = "f", digits = 4),
  MSE_phat_code  = formatC(tab$mse_phat, format = "f", digits = 6),
  agree = ifelse(round(tab$mse_pM, 4) == prnt[, 1] &
                 round(tab$mse_phat, 4) == prnt[, 2], "yes", "NO")))
say("")
say("largest departure from the 32 values printed in the paper : ",
    formatC(max(abs(ours - PAPER_TABLE1)), format = "e", digits = 2))
stopifnot(max(abs(ours - PAPER_TABLE1)) < 1e-12)
say("the table of the paper is reproduced exactly")

## the LaTeX of the table, in the two-block layout of the paper
half <- nrow(tab) / 2L
emit_table(
  "table1_mse",
  body = cbind(formatC(tab$p[1:half], format = "f", digits = 2),
               formatC(tab$n[1:half], format = "d"),
               fmt4(tab$mse_pM[1:half]), fmt4(tab$mse_phat[1:half]),
               formatC(tab$n[half + 1:half], format = "d"),
               fmt4(tab$mse_pM[half + 1:half]),
               fmt4(tab$mse_phat[half + 1:half])),
  colspec = "ccrrcrr",
  header = paste0("$p$ & $n$ & $MSE(\\hat p_M)$ & $MSE(\\hat p)$ & $n$ & ",
                  "$MSE(\\hat p_M)$ & $MSE(\\hat p)$"),
  caption = "The MSEs of the maximum likelihood estimators.",
  label = "tab:mse", rules = 4L,
  note = paste0("Note: the wind speed is Weibull with $k = ", PAPER$k,
                "$ and $\\alpha = ", PAPER$alpha, "$, and the wind turbine ",
                "has $v_{ci} = ", v_ci, "$, $v_r = ", v_r, "$, $v_{co} = ",
                v_co, "$ m/s and $P_r = ", P_r, "$. The values are exact, ",
                "(15) with $C = 0$ giving $MSE(\\hat p)$ and the two ",
                "moments (18) giving $MSE(\\hat p_M)$."))

## =================================================================
##  Figure 1: the MSEs as functions of p
## =================================================================
say_rule("Figure 1 of the paper, the MSEs as functions of p")
pp <- seq(0.01, 1, by = 0.01)
f1 <- lapply(PAPER$n_tab, function(n)
  list(n = n,
       mle = mse_mle_wt(pp, n, W, v_ci, v_co),
       mom = sapply(pp, function(p)
         mse_moment_wt(p, n, W, v_ci, v_r, v_co, P_r))))
for (z in f1)
  say("n = ", formatC(z$n, width = 2), " : at p = 1 the MSE of p_hat is ",
      fmt5(z$mle[length(pp)]), " and that of p_tilde ",
      fmt5(z$mom[length(pp)]), " , their ratio ",
      fmt4(z$mom[length(pp)] / z$mle[length(pp)]))

fig_device("figure1_mse_p", quote({
  for (z in f1) {
    plot(pp, z$mom, type = "l", lty = 2, lwd = 2, col = "grey35",
         ylim = c(0, max(z$mom)), xlab = "p", ylab = "MSE",
         main = paste0("n = ", z$n), cex.main = 1)
    lines(pp, z$mle, lty = 1, lwd = 2)
    legend("topleft", bty = "n", lwd = 2, lty = c(1, 2), cex = 0.9,
           col = c("black", "grey35"),
           legend = c(expression(hat(p)), expression(tilde(p))))
  }
}), width = 6.0, height = 5.0, mfrow = c(2, 2))

## =================================================================
##  Figure 2: the MSEs as functions of n
## =================================================================
say_rule("Figure 2 of the paper, the MSEs as functions of n")
nn <- PAPER$n_fig
f2 <- lapply(PAPER$p_tab, function(p)
  list(p = p,
       mle = sapply(nn, function(n) mse_mle_wt(p, n, W, v_ci, v_co)),
       mom = sapply(nn, function(n)
         mse_moment_wt(p, n, W, v_ci, v_r, v_co, P_r))))
for (z in f2)
  say("p = ", formatC(z$p, format = "f", digits = 2),
      " : at n = ", min(nn), " the MSE of p_hat is ", fmt5(z$mle[1]),
      " and that of p_tilde ", fmt5(z$mom[1]),
      " ; at n = ", max(nn), " they are ", fmt5(z$mle[length(nn)]),
      " and ", fmt5(z$mom[length(nn)]))

fig_device("figure2_mse_n", quote({
  for (z in f2) {
    plot(nn, z$mom, type = "l", lty = 2, lwd = 2, col = "grey35",
         ylim = c(0, max(z$mom)), xlab = "n", ylab = "MSE",
         main = paste0("p = ", formatC(z$p, format = "f", digits = 2)),
         cex.main = 1)
    lines(nn, z$mle, lty = 1, lwd = 2)
    legend("topright", bty = "n", lwd = 2, lty = c(1, 2), cex = 0.9,
           col = c("black", "grey35"),
           legend = c(expression(hat(p)), expression(tilde(p))))
  }
}), width = 6.0, height = 5.0, mfrow = c(2, 2))

## =================================================================
##  Table 2: the effect of the clustering, exactly
## =================================================================
## The sample of Section 2.2 is N periods of n WTs, M = Nn cells and
## C = n - 1.  Nothing is simulated here either: MSE(p_hat) is (14),
## and MSE(p_M) is computed from the exact distribution of m_0, which
## is the N-fold convolution of (16).  The column n = 1 is the case
## C = 0 and must return the entries of Table 1.
say_rule("Table 2, the effect of reading n WTs in one period")
## Everything here comes from code/core_wtrel.R, so that the table of
## the paper and the formulas of Section 2 cannot part company.
d_p  <- d_wt(W, v_ci, v_co)
m1_p <- Eg_wt(W, v_ci, v_r, v_co, P_r, 1)
m2_p <- Eg_wt(W, v_ci, v_r, v_co, P_r, 2)

M_fix  <- 30L                       # the same number of readings throughout
n_set  <- c(1L, 2L, 3L, 6L)         # WTs read in one period
tab2 <- expand.grid(p = PAPER$p_tab, n = n_set)
tab2$N       <- M_fix / tab2$n
tab2$C       <- tab2$n - 1L
tab2$mse_pt  <- mse_moment_C(tab2$p, M_fix, tab2$C, m1_p, m2_p)
tab2$mse_ph  <- mse_mle_C(tab2$p, M_fix, tab2$C, d_p)
tab2$mse_pM  <- mapply(function(p, N, n) mse_mmle_C(p, rep(n, N), d_p),
                       tab2$p, tab2$N, tab2$n)
tab2$mse_ML  <- mapply(function(p, N, n) unname(mse_ml_C(p, N, n, d_p)["mse"]),
                       tab2$p, tab2$N, tab2$n)

say("M = ", M_fix, " readings in every design, N = M/n periods of n WTs")
say("rho_w of (16) at p = ", fmt4(PAPER$p_tab[1]), " is ",
    fmt5(rho_w_wt(PAPER$p_tab[1], d_p)), " , at p = ",
    fmt4(PAPER$p_tab[4]), " it is ", fmt5(rho_w_wt(PAPER$p_tab[4], d_p)))
say("rho_P of (16) at the same two values : ",
    fmt5(rho_P_wt(PAPER$p_tab[1], m1_p, m2_p)), " and ",
    fmt5(rho_P_wt(PAPER$p_tab[4], m1_p, m2_p)))
say_block(data.frame(p = formatC(tab2$p, format = "f", digits = 2),
                     n = formatC(tab2$n, format = "d"),
                     N = formatC(tab2$N, format = "d"),
                     C = formatC(tab2$C, format = "d"),
                     MSE_ptilde = fmt4(tab2$mse_pt),
                     MSE_phat   = fmt4(tab2$mse_ph),
                     MSE_pM     = fmt4(tab2$mse_pM),
                     MSE_pML    = fmt4(tab2$mse_ML)))

## the column n = 1 must reproduce the n = 30 block of Table 1
one <- tab2[tab2$n == 1L, ]
ref <- tab[tab$n == 30L, ]
say("")
say("the column n = 1 against the n = 30 block of Table 1")
say_block(data.frame(p = formatC(one$p, format = "f", digits = 2),
                     MSE_pM_new = fmt4(one$mse_pM),
                     MSE_pM_tab1 = fmt4(ref$mse_pM),
                     MSE_phat_new = fmt4(one$mse_ph),
                     MSE_phat_tab1 = fmt4(ref$mse_phat),
                     MSE_pML_new = fmt4(one$mse_ML)))
stopifnot(max(abs(one$mse_ph - ref$mse_phat)) < 1e-12)
stopifnot(max(abs(one$mse_pM - ref$mse_pM)) < 1e-10)
stopifnot(max(abs(one$mse_ML - ref$mse_pM)) < 1e-8)
say("the case C = 0 returns Table 1 exactly, and p_ML is then p_M itself")

say("")
say("the MSE of p_hat and of p_ML, relative to one WT per period")
for (p0 in PAPER$p_tab) {
  z <- tab2[tab2$p == p0, ]
  say("p = ", formatC(p0, format = "f", digits = 2),
      " : p_hat ", paste0("n=", z$n, " -> ", fmt4(z$mse_ph / z$mse_ph[1]),
                          collapse = " , "))
  say("            p_ML   ",
      paste0("n=", z$n, " -> ", fmt4(z$mse_ML / z$mse_ML[1]), collapse = " , "))
}

emit_table(
  "table2_design",
  body = cbind(formatC(tab2$p, format = "f", digits = 2),
               formatC(tab2$n, format = "d"),
               formatC(tab2$N, format = "d"),
               formatC(tab2$C, format = "d"),
               fmt4(tab2$mse_pt), fmt4(tab2$mse_ph),
               fmt4(tab2$mse_pM), fmt4(tab2$mse_ML)),
  colspec = "cccc rrrr",
  header = paste0("$p$ & $n$ & $N$ & $C$ & $MSE(\\tilde p)$ & ",
                  "$MSE(\\hat p)$ & $MSE(\\hat p_M)$ & $MSE(\\hat p_{ML})$"),
  caption = paste0("The MSEs of the four estimators when the same ",
                   "$M = ", M_fix, "$ readings are collected as $N$ periods ",
                   "of $n$ wind turbines."),
  label = "tab:design", rules = 4L,
  note = paste0("Note: the wind speed is Weibull with $k = ", PAPER$k,
                "$ and $\\alpha = ", PAPER$alpha, "$, and the wind turbine ",
                "has $v_{ci} = ", v_ci, "$, $v_r = ", v_r, "$, $v_{co} = ",
                v_co, "$ m/s and $P_r = ", P_r, "$, so that $d = ",
                fmt5(d_p), "$. Every value is exact and nothing is ",
                "simulated. The row $n = 1$ is the case $C = 0$: it ",
                "returns the $n = 30$ block of Table 1, and there ",
                "$\\hat p_{ML}$ is $\\hat p_M$ itself."))

## a figure of the same thing: the two MSEs against n at a fixed M
say("")
say("Figure 3, MSE against n at M = ", M_fix, " , for p_hat and p_ML")
nn3 <- c(1L, 2L, 3L, 5L, 6L, 10L, 15L)      # the divisors of M
fig3 <- lapply(PAPER$p_tab, function(p0) {
  data.frame(n = nn3,
             ph = mse_mle_C(p0, M_fix, nn3 - 1L, d_p),
             ml = sapply(nn3, function(n)
                    unname(mse_ml_C(p0, M_fix / n, n, d_p)["mse"])))
})
names(fig3) <- formatC(PAPER$p_tab, format = "f", digits = 2)
for (i in seq_along(fig3))
  say("p = ", names(fig3)[i], " : p_hat ",
      paste(fmt4(fig3[[i]]$ph), collapse = " "), "  |  p_ML ",
      paste(fmt4(fig3[[i]]$ml), collapse = " "))
fig_device("figure3_mse_cluster", quote({
  for (i in seq_along(fig3)) {
    z <- fig3[[i]]
    plot(z$n, z$ph, type = "b", lwd = 2, pch = 16, cex = 0.6,
         ylim = c(0, max(z$ph)),
         xlab = "n, wind turbines read in one period", ylab = "MSE",
         main = paste0("p = ", names(fig3)[i]), cex.main = 1)
    lines(z$n, z$ml, lwd = 2, lty = 2)
    points(z$n, z$ml, pch = 1, cex = 0.6)
    if (i == 1L)
      legend("topleft", c(expression(hat(p)), expression(hat(p)[ML])),
             lty = c(1, 2), lwd = 2, bty = "n", cex = 0.9)
  }
}), width = 6.0, height = 5.0, mfrow = c(2, 2))

say("")
say("figures written to ", fig_dir, "/ , tables written to \"", tex_dir, "\"/")
saveRDS(list(tab = tab, tab2 = tab2, fig3 = fig3, f1 = f1, f2 = f2, setting = PAPER,
             d = d_p, M_fix = M_fix),
        file.path("data", "paper_numbers.rds"))
close_report()
