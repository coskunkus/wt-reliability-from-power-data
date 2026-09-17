## =================================================================
##  run_example.R
##
##  The application of the rewritten Section 2 to the Kelmarsh data:
##  the sample of N periods and the n_j wind turbines that remain in
##  each, the two estimators, their variances with the clustering
##  constant C, and the availability indicators of the same records.
##
##  It produces the three tables of the example and splices them into
##  the section between the markers, so that nothing has to be \input:
##
##    table_example_sample     the sample after Steps 1 and 2: N, the
##                             distribution of n_j, M, C, and the N_i
##    table_example_estimates  the quantities the estimators require,
##                             the estimates, and four availability
##                             indicators
##    table_example_variance   the variance of (15): C and rho_w of the
##                             record against (16), and the two designs
##                             of the replication of Step 3
##
##  Output : results/run_example_report.txt
##           "Claude outputs"/table_example_{sample,estimates,variance}.tex
##           the tables spliced into "Claude outputs"/section_example.tex
##
##  Run after code/prepare_data.R and code/run_simulation.R.
##  Base R only.  Run time is under a minute.
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
for (dd in c(res_dir, fig_dir, tex_dir))
  if (!dir.exists(dd)) dir.create(dd, recursive = TRUE)

if (!file.exists(PANEL_RDS))
  stop("\n\n  The file  ", PANEL_RDS, "  is missing.\n",
       "  It is produced by  code/prepare_data.R\n",
       "  Run that script first, then run this one again.\n", call. = FALSE)

new_report(file.path(res_dir, "run_example_report.txt"))
say("run_example.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

D    <- readRDS(PANEL_RDS)
pan  <- D$panel
v_ci <- WT$v_ci; v_r <- WT$v_r; v_co <- WT$v_co; P_r <- WT$P_r

## =================================================================
##  1. The sample of Section 2.2: N, n_j, M and C
## =================================================================
say_rule("1. the sample after Steps 1 and 2")
sub       <- estimation_sample(pan)
sub$ib    <- in_band(sub$ws, v_ci, v_co)
sub$Z     <- as.numeric(zero_of(sub, v_ci, v_co))
sub$pw    <- ifelse(sub$ib, sub$power, 0)   # the power output of the model

slots_all <- length(unique(pan$slot))       # every ten-minute period
nj_tab    <- table(sub$slot)
nj        <- as.numeric(nj_tab)
N         <- length(nj)                     # periods with at least one cell
M         <- sum(nj)
Csum      <- sum(nj * (nj - 1))
C         <- Csum / M
Ni        <- as.numeric(table(sub$turbine))
n_wt      <- length(Ni)

say("periods in the record, N in the widest sense : ", fmti(slots_all))
say("periods with at least one WT retained, N     : ", fmti(N))
say("periods with none                            : ", fmti(slots_all - N))
say("cells retained, M = sum n_j                  : ", fmti(M))
say("                  = sum N_i                  : ", fmti(sum(Ni)))
stopifnot(sum(Ni) == M)
say("")
say("the distribution of n_j")
cnt <- table(factor(nj, levels = 1:n_wt))
say_block(data.frame(n_j = names(cnt), periods = fmti(as.numeric(cnt)),
                     cells = fmti(as.numeric(names(cnt)) * as.numeric(cnt))))
say("")
say("sum n_j(n_j-1)                               : ", fmti(Csum))
say("C = sum n_j(n_j-1) / M                       : ", fmt4(C))
say("C of a balanced sample of six WTs would be   : 5.0000")
say("")
say("N_i, the periods each WT contributes")
say_block(data.frame(WT = seq_len(n_wt), N_i = fmti(Ni)))
say("sum N_i(N_i-1) / M, the constant of the other grouping : ",
    fmti(sum(Ni * (Ni - 1)) / M))

## =================================================================
##  2. The quantities the estimators require
## =================================================================
say_rule("2. the quantities of Section 2.1")
fw    <- fit_weibull_ml(sub$ws)
k_ml  <- unname(fw["k"]); a_ml <- unname(fw["c"])
W_ml  <- wind_weibull(k_ml, a_ml)
d_ml  <- d_wt(W_ml, v_ci, v_co)
d_rec <- mean(sub$ib)
se_d  <- sqrt(d_rec * (1 - d_rec) / M)
m1    <- Eg_wt(W_ml, v_ci, v_r, v_co, P_r, 1)
m2    <- Eg_wt(W_ml, v_ci, v_r, v_co, P_r, 2)
say("Weibull by maximum likelihood : k = ", fmt4(k_ml),
    " , alpha = ", fmt4(a_ml), " m/s")
say("d of the fitted F             : ", fmt5(d_ml))
say("d of the wind speed record    : ", fmt5(d_rec),
    "   standard error ", fmt5(se_d))
say("m1 = E g(V)                   : ", fmt4(m1), " kW")
say("m2 = E g(V)^2                 : ", formatC(m2, format = "g", digits = 6))

## =================================================================
##  3. The estimates of Section 2.3 and 2.4
## =================================================================
say_rule("3. the estimates")
m0     <- sum(sub$Z)
p_rec  <- p_hat_wt(m0, M, d_rec)
p_wml  <- p_hat_wt(m0, M, d_ml)
pM_wml <- p_M_wt(m0, M, d_ml)
p_til  <- mean(sub$pw) / m1
say("m0                            : ", fmti(m0),
    "   m0/M = ", fmt5(m0 / M))
say("p_hat of (14), d of the record : ", fmt5(p_rec))
say("p_hat of (14), d of the fitted F : ", fmt5(p_wml),
    "  ->  p_M = ", fmt5(pM_wml))
say("p_tilde of (10)                : ", fmt5(p_til),
    "   (inadmissible, the cubic law understating the power curve)")
say("observed mean power output    : ", fmt4(mean(sub$power)), " kW")
say("mean of the model's power     : ", fmt4(mean(sub$pw)), " kW")

## ---- why the d of the record is the one to use -------------------
## With d estimated by the fraction B/M of the cells whose wind speed
## lies in the band, (14) reduces to the ratio of the cells that produce
## power to the cells that could have produced it.  A cell that produces
## power has its wind speed in the band, so that the numerator never
## exceeds the denominator and p_hat <= 1 holds for every sample.  The
## d of the fitted F carries no such guarantee.
say("")
Bib <- sum(sub$ib)                      # cells with the wind in the band
Rpos <- M - m0                          # cells that produce power
say("cells whose wind lies in the band, B          : ", fmti(Bib))
say("cells that produce power, M - m0              : ", fmti(Rpos))
say("d of the record is B/M                        : ", fmt5(Bib / M))
say("(14) is then (M - m0)/B                       : ", fmt7(Rpos / Bib))
stopifnot(Rpos <= Bib)                  # a producing cell is in the band
stopifnot(abs(p_rec - Rpos / Bib) < 1e-12)
say("a producing cell is in the band, so that the ratio cannot exceed")
say("one: with this d the estimate is admissible for every sample, and")
say("the modification of Section 2.4 is never called for")
say("with the d of the fitted F the estimate is         ", fmt5(p_wml),
    " , above one by ", fmt5(p_wml - 1))
say("the two values of d differ by                      ", fmt5(d_rec - d_ml),
    " , that is ", fmt4(100 * (d_rec - d_ml) / d_rec), " per cent")

## its standard error: a proportion over the B cells, inflated by the
## clustering of the wind farm, both measured on the same records
ibs   <- sub[sub$ib, ]
njb   <- as.numeric(table(ibs$slot))
Cb    <- sum(njb * (njb - 1)) / sum(njb)
sb    <- as.numeric(tapply(ibs$Z, ibs$slot, sum))
se_bin <- sqrt(p_rec * (1 - p_rec) / Bib)
deff   <- (sum((sb - njb * (1 - p_rec))^2) / sum(njb)) / (p_rec * (1 - p_rec))
se_p   <- se_bin * sqrt(deff)
say("")
say("binomial standard error of (M - m0)/B         : ",
    formatC(se_bin, format = "e", digits = 3))
say("C of the cells in the band                    : ", fmt4(Cb))
say("design effect measured on the same cells      : ", fmt4(deff),
    "   (rho = ", fmt5((deff - 1) / Cb), ")")
say("standard error of p_hat                       : ",
    formatC(se_p, format = "e", digits = 3))
say("p_hat and twice its standard error            : ", fmt5(p_rec),
    " +- ", fmt5(2 * se_p))
## what the residual dependence is: the failures of the cells in the band
## do not fall singly, a site-wide outage taking several WTs at once
fullb <- njb == n_wt
obs   <- as.numeric(table(factor(sb[fullb], levels = 0:n_wt)))
exp0  <- sum(fullb) * dbinom(0:n_wt, n_wt, 1 - p_rec)
say("")
say("periods in which all six WTs have their wind in the band : ",
    fmti(sum(fullb)))
say("the number of them that fail together, against what independence")
say("would give")
say_block(data.frame(failures = 0:n_wt, observed = fmti(obs),
                     independent = formatC(exp0, format = "f", digits = 1)))
say("single failures fall short and joint failures are far in excess:")
say("nine periods lose all six at once, which independence makes")
say("impossible, so that the residual correlation is the signature of")
say("the outages that strike the whole site at one time")
say("")
say("the correlation of the failure indicators of the cells in the band")
say("is ", fmt5((deff - 1) / Cb), " ; section 4 below measures that of the ",
    "zero indicators themselves, which is far larger, conditioning on a ",
    "wind in the band removing most of the dependence")

## =================================================================
##  4. The correlations of (16) and the variance of (15)
## =================================================================
say_rule("4. the clustering, and the variance of (15)")
full <- as.integer(names(nj_tab))[nj == n_wt]
fs   <- sub[sub$slot %in% full, ]
fs   <- fs[order(fs$slot, fs$turbine), ]
Zm   <- matrix(fs$Z,  ncol = n_wt, byrow = TRUE)
Pm   <- matrix(fs$pw, ncol = n_wt, byrow = TRUE)
cw   <- cor(Zm); cp <- cor(Pm)
rho_w <- mean(cw[upper.tri(cw)])
rho_P <- mean(cp[upper.tri(cp)])
rho_w_th <- p_rec * (1 - d_rec) / (1 - p_rec * d_rec)
rho_P_th <- p_rec * (m2 - m1^2) / (m2 - p_rec * m1^2)
say("rho_w of the record           : ", fmt5(rho_w),
    "   (the ", n_wt * (n_wt - 1) / 2, " pairs run from ",
    fmt4(min(cw[upper.tri(cw)])), " to ", fmt4(max(cw[upper.tri(cw)])), ")")
say("rho_w of (16), a common wind  : ", fmt5(rho_w_th))
say("rho_P of the record           : ", fmt5(rho_P))
say("rho_P of (16), a common wind  : ", fmt5(rho_P_th))
say("")
## Eq. (15) with the measured rho_w.  The periods of the whole panel are
## ten minutes apart and are not independent, so that this is a bound
## and not the standard error of the estimate of the whole record.
var14 <- function(p, d, M, C, rw) p * (1 - p * d) / (M * d) * (1 + rw * C)
sd_all <- sqrt(var14(p_rec, d_rec, M, C, rho_w))
say("(15) with C = ", fmt4(C), " and the measured rho_w, over the whole ",
    "record : sd = ", formatC(sd_all, format = "e", digits = 3))
say("the periods of the whole record are ten minutes apart and are not")
say("independent, so that this value is a bound and not a standard error")

## the two designs of the replication of Step 3
des <- NULL
f <- file.path("data", "design_study.rds")
if (file.exists(f)) {
  des <- readRDS(f)$tab
  MA <- des$N[des$design == "A"]; MB <- des$N[des$design == "B"]
  sA <- des$sd_phat[des$design == "A"]; sB <- des$sd_phat[des$design == "B"]
  vA <- var14(p_rec, d_rec, MA, n_wt - 1, rho_w)
  vB <- var14(p_rec, d_rec, MB, 0, rho_w)
  say("")
  say("the replication of Step 3, 1000 times")
  say("  design A, six WTs at each period, C = 5")
  say("    M = ", fmt4(MA), " , (15) gives sd ", fmt5(sqrt(vA)),
      " , simulated ", fmt5(sA))
  say("  design B, one WT at each period, C = 0")
  say("    M = ", fmt4(MB), " , (15) gives sd ", fmt5(sqrt(vB)),
      " , simulated ", fmt5(sB))
  say("  the ratio of the two simulated variances, corrected for M : ",
      fmt4((sA / sB)^2 * MA / MB))
  say("  1 + (n-1) rho_w of the record                             : ",
      fmt4(1 + (n_wt - 1) * rho_w))
  say("  each design falls below (15) by the gain of the stratification:")
  say("    design A ", fmt4(sA^2 / vA), " , design B ", fmt4(sB^2 / vB))
} else {
  say("data/design_study.rds is absent; run code/run_simulation.R for the",
      " two designs")
}

## =================================================================
##  5. The likelihood of Remark 1, the wind of a period being common
## =================================================================
say_rule("5. the estimator of Remark 1, a wind common to the period")
## the zeros of each period, in the order of nj
m0j <- as.numeric(tapply(sub$Z, sub$slot, sum))
stopifnot(length(m0j) == N, sum(m0j) == m0)
idle <- m0j >= nj                       # periods idle throughout
S_st <- sum(m0j[!idle]); R_st <- sum(nj[!idle] - m0j[!idle])
say("periods in which every retained WT is idle   : ", fmti(sum(idle)),
    "   (", fmt4(100 * mean(idle)), " per cent)")
say("cells of the other periods, zeros R and S    : ",
    fmti(R_st), " positive , ", fmti(S_st), " idle")
say("in those periods the wind is known to lie in the band, so that")
say("their zeros are failures and their share of idle cells is ",
    fmt5(S_st / (R_st + S_st)))

p_ML_rec <- p_ml_wt(m0j, nj, d_rec)
p_ML_wml <- p_ml_wt(m0j, nj, d_ml)
se_ML    <- se_ml_wt(p_ML_wml, nj, d_ml)
say("")
say("p_ML, d of the wind speed record             : ", fmt5(p_ML_rec))
say("p_ML, d of the fitted F                      : ", fmt5(p_ML_wml))
say("its standard error, 1/sqrt(sum I_{n_j})      : ",
    formatC(se_ML, format = "e", digits = 3))
say("p_hat of (14) for comparison                 : ", fmt5(p_wml),
    "   p_M = ", fmt5(pM_wml))
say("the three estimates differ by at most ",
    fmt5(max(abs(outer(c(p_ML_wml, p_wml, pM_wml), c(p_ML_wml, p_wml, pM_wml),
                       "-")))))
say("")
say("the standard error of (15) with C = ", fmt4(C),
    " for the same sample : ", formatC(sd_all, format = "e", digits = 3))
say("the ratio of the two standard errors         : ",
    fmt4(sd_all / se_ML))
say("both are bounds rather than standard errors, the periods of the")
say("whole record being ten minutes apart and not independent")

## ---- why p_ML falls below the other estimates --------------------
## The estimator reads one and the same wind speed over a period.  Each
## WT carries its own nacelle anemometer, so that a WT may be below the
## cut-in speed while its neighbour is above it; the estimator charges
## such a zero to a failure.  The record says how often that happens.
say("")
say_rule("5b. the wind of a period is common only to a degree")
inp  <- sub[sub$slot %in% as.integer(names(nj_tab))[!idle], ]
z_out <- sum(inp$Z == 1 & !inp$ib)
z_in  <- sum(inp$Z == 1 &  inp$ib)
say("cells of the periods that are not idle throughout : ", fmti(nrow(inp)))
say("  their zeros                                     : ", fmti(z_out + z_in))
say("    the WT's own wind speed outside the band      : ", fmti(z_out),
    "   (", fmt4(100 * z_out / (z_out + z_in)), " per cent)")
say("    the WT's own wind speed inside the band       : ", fmti(z_in))
say("a common wind would make every one of these zeros a failure, which")
say("is what sends p_ML down to ", fmt5(p_ML_wml), " ; counting as failures")
say("only the zeros of the cells whose own wind lies in the band gives")
say("  1 - (zeros in band) / (cells in band)           : ",
    fmt5(1 - sum(sub$Z == 1 & sub$ib) / sum(sub$ib)))
say("which is (14) read cell by cell, and agrees with p_hat of ",
    fmt5(p_rec), " to ", fmt5(abs(p_rec - (1 - sum(sub$Z == 1 & sub$ib) /
                                           sum(sub$ib)))))
Wm  <- matrix(fs$ws, ncol = n_wt, byrow = TRUE)
cws <- cor(Wm); cor_ws_mean <- mean(cws[upper.tri(cws)])
say("the six nacelle anemometers correlate ",
    fmt4(cor_ws_mean), " with one another, so that the wind of a period is")
say("common only to a degree, and (14) is the estimator to use in the field")

## =================================================================
##  6. The availability indicators of the same cells
## =================================================================
say_rule("6. the availability indicators, same cells")
A_FO   <- 1 - mean(sub$ev_fo)
A_FOTS <- 1 - mean(sub$ev_fo | sub$ev_ts)
A_sys  <- mean(sub$av_sys, na.rm = TRUE)
A_b32  <- mean(sub$av_b32, na.rm = TRUE)
say("1 - FOR, forced outage                       : ", fmt5(A_FO))
say("forced outage or technical standby           : ", fmt5(A_FOTS))
say("operator's time-based system availability    : ", fmt5(A_sys))
say("operator's time-based IEC 61400-26 B.3.2     : ", fmt5(A_b32))
dif <- p_rec - c(A_FO, A_FOTS, A_sys, A_b32)
say("p_hat minus each, in percentage points       : ",
    paste(formatC(100 * dif, format = "f", digits = 2), collapse = " , "))
say("largest departure, in percentage points      : ",
    formatC(100 * max(abs(dif)), format = "f", digits = 2))

## =================================================================
##  7. The tables of the example
## =================================================================
say_rule("7. the tables of the example")

## ---- Table: the sample -------------------------------------------
lab <- names(cnt); per <- as.numeric(cnt)
body1 <- rbind(
  c("Ten-minute periods in the record", fmti(slots_all)),
  c("Periods with at least one WT retained, $N$", fmti(N)),
  c("Cells retained, $M=\\sum_j n_j=\\sum_i N_i$", fmti(M)),
  c("$\\sum_j n_j(n_j-1)$", fmti(Csum)),
  c("$C=\\sum_j n_j(n_j-1)/M$", fmt4(C)),
  cbind(paste0("Periods with $n_j = ", rev(lab), "$"), fmti(rev(per))),
  cbind(paste0("Periods contributed by WT ", seq_len(n_wt), ", $N_i$"),
        fmti(Ni)))
emit_table(
  tag("table_example_sample"),
  body = body1, colspec = "lr",
  header = "Quantity & Value",
  caption = paste0("The sample of the wind farm after Steps 1 and 2 ",
                   "(Kelmarsh, 2017 to 2019)."),
  label = "tab:sample", rules = 5L,
  note = paste0("Note: a cell is one wind turbine in one ten-minute period. ",
                "Step 1 removes the cells without a valid power or wind ",
                "speed signal and Step 2 those falling inside an event of ",
                "the categories scheduled maintenance, out of electrical ",
                "specification and requested shutdown, one cell at a time, ",
                "so that the number $n_j$ of wind turbines that remain ",
                "varies from period to period. The constant $C$ would be ",
                "five for a sample in which every period kept all six."))

## ---- Table: the estimates ----------------------------------------
body2 <- rbind(
  c("Nominal power $P_r$", "2050 kW, manufacturer"),
  c("Cut-in wind speed $v_{ci}$", paste0(fmt4(v_ci), " m/s, manufacturer")),
  c("Rated wind speed $v_r$", paste0(fmt4(v_r), " m/s, manufacturer")),
  c("Cut-out wind speed $v_{co}$", paste0(fmt4(v_co), " m/s, manufacturer")),
  c("Wind speed distribution $F$",
    paste0("Weibull, $\\hat k = ", fmt4(k_ml), "$, $\\hat\\alpha = ",
           fmt4(a_ml), "$ m/s")),
  c("$d=F(v_{co})-F(v_{ci})$ of the fitted $F$", fmt5(d_ml)),
  c("$d$ of the wind speed record", fmt5(d_rec)),
  c("$m_1=E\\{g(V)\\}$", paste0(fmt4(m1), " kW")),
  c("$\\hat p$ of (14), $d$ of the wind speed record", fmt5(p_rec)),
  c("$\\hat p_M$ of (14), $d$ of the fitted $F$", fmt5(pM_wml)),
  c("$1-FOR$, forced outage", fmt5(A_FO)),
  c("Forced outage or technical standby", fmt5(A_FOTS)),
  c("Operator's time-based system availability", fmt5(A_sys)),
  c("Operator's IEC 61400-26 B.3.2 availability", fmt5(A_b32)))
emit_table(
  tag("table_example_estimates"),
  body = body2, colspec = "lr",
  header = "Quantity & Value",
  caption = paste0("The quantities that the estimators of Section 2 require, ",
                   "the two estimates of $p$, and four availability ",
                   "indicators of the same wind turbines over the same ",
                   "records (Kelmarsh, 2017 to 2019)."),
  label = "tab:estimates", rules = 8L,
  note = paste0("Note: $M = ", fmti(M), "$, $m_0 = ", fmti(m0),
                "$ and $m_0/M = ", fmt5(m0 / M),
                "$. The first estimate is $(M-m_0)/B = ", fmti(M - m0), "/",
                fmti(Bib), "$, the ratio of the cells that produce power to ",
                "the cells whose wind speed lies in the band, and its ",
                "standard error is ", formatC(se_p, format = "f", digits = 6),
                ", the clustering inflating the binomial value ",
                formatC(se_bin, format = "f", digits = 6), " by ", fmt4(deff),
                ". A cell that produces power has its wind speed in the ",
                "band, so that this estimate cannot exceed one, whereas the ",
                "estimate that the fitted $F$ supplies carries no such ",
                "guarantee and is ", fmt5(p_wml), ". ",
                "The parameters of the wind speed distribution are ",
                "estimated by maximum likelihood from the wind speeds of ",
                "the same cells, and the standard error of the proportion ",
                "that gives the second value of $d$ is ", fmt5(se_d),
                ". A cell whose wind speed lies outside $[v_{ci}, v_{co})$ ",
                "counts as a zero-power observation, as $g$ of (2) ",
                "prescribes. The moment estimator (10) takes the ",
                "inadmissible value ", fmt4(p_til),
                ". The first two indicators are one minus the fraction of ",
                "the cells that fall inside an event of the category named, ",
                "a technical standby being charged only when it stops the ",
                "wind turbine; the last two are the availability signals ",
                "that the operator's own system reports, averaged over the ",
                "same cells."))

## ---- Table: the variance -----------------------------------------
if (!is.null(des)) {
  body3 <- rbind(
    c("$\\rho_w$ of the wind speed record", fmt5(rho_w), ""),
    c("$\\rho_w$ of (16), one and the same wind", fmt5(rho_w_th), ""),
    c("$\\rho_P$ of the wind speed record", fmt5(rho_P), ""),
    c("$\\rho_P$ of (16), one and the same wind", fmt5(rho_P_th), ""),
    c("$1+(n-1)\\rho_w$ with the measured $\\rho_w$",
      fmt4(1 + (n_wt - 1) * rho_w), ""),
    c(paste0("Design A, six WTs in each period, $M = ", fmt4(MA),
             "$"), fmt5(sqrt(vA)), fmt5(sA)),
    c(paste0("Design B, one WT in each period, $M = ", fmt4(MB),
             "$"), fmt5(sqrt(vB)), fmt5(sB)))
  emit_table(
    tag("table_example_variance"),
    body = body3, colspec = "lrr",
    header = "Quantity & Eq. (15) or (16) & Simulated",
    caption = paste0("The clustering of the wind farm and the standard ",
                     "deviation of $\\hat p$ under the two designs of the ",
                     "replication of Step 3 (Kelmarsh, 2017 to 2019)."),
    label = "tab:variance", rules = 5L,
    note = paste0("Note: $\\rho_w$ and $\\rho_P$ are averaged over the ",
                  n_wt * (n_wt - 1) / 2, " pairs of wind turbines of the ",
                  fmti(length(full)), " periods in which all six remain, ",
                  "and run from ", fmt4(min(cw[upper.tri(cw)])), " to ",
                  fmt4(max(cw[upper.tri(cw)])), " for $\\rho_w$. Step 3 is ",
                  "replicated 1000 times, with $22$ periods in each of the ",
                  "twenty-four strata of month and of day and night and at ",
                  "least forty-eight hours between two periods. The two ",
                  "designs differ in $n$ alone, so that the ratio of their ",
                  "variances is $1+(n-1)\\rho_w$; it is ",
                  fmt4((sA / sB)^2 * MA / MB), " against the ",
                  fmt4(1 + (n_wt - 1) * rho_w), " of (15). Each design ",
                  "falls below (15) taken alone, by ", fmt4(sA^2 / vA),
                  " and ", fmt4(sB^2 / vB),
                  ", the periods being drawn within the strata rather than ",
                  "at random."))
}

## ---- Table: what a further WT of the same period is worth --------
## At a fixed number of periods, reading one more WT adds a cell but
## not an independent one.  The variance of (15) relative to the
## design that reads one WT in each of the same periods is
##    {1 + (n-1) rho_w} / n ,
## and n divided by that factor is the number of independent readings
## that a period carries.  Nothing is simulated here.
say_rule("8. the value of reading n WTs in one period")
nn   <- 1:n_wt
fac  <- function(rw, n) (1 + (n - 1) * rw) / n
r_emp <- fac(rho_w, nn);     r_th <- fac(rho_w_th, nn)
e_emp <- nn / (1 + (nn - 1) * rho_w)
e_th  <- nn / (1 + (nn - 1) * rho_w_th)
say("rho_w measured ", fmt5(rho_w), " , rho_w of (16) ", fmt5(rho_w_th))
say_block(data.frame(n = nn,
                     ratio_measured = fmt4(r_emp),
                     readings_measured = fmt4(e_emp),
                     ratio_common_wind = fmt4(r_th),
                     readings_common_wind = fmt4(e_th)))
say("as n grows the ratio tends to rho_w itself, ", fmt4(rho_w),
    " measured and ", fmt4(rho_w_th), " under a common wind,")
say("and the number of independent readings to 1/rho_w, ",
    fmt4(1 / rho_w), " and ", fmt4(1 / rho_w_th))

body4 <- rbind(
  cbind(formatC(nn, format = "d"), fmt4(r_emp), fmt4(e_emp)),
  c("$\\infty$", fmt4(rho_w), fmt4(1 / rho_w)))
emit_table(
  tag("table_example_design"),
  body = body4, colspec = "crr",
  header = "$n$ & Variance & Independent readings",
  caption = paste0("The value of reading $n$ wind turbines in each period ",
                   "instead of one, at the same number of periods ",
                   "(Kelmarsh, 2017 to 2019)."),
  label = "tab:value", rules = 6L,
  note = paste0("Note: the variance is that of (15) relative to the design ",
                "that reads one wind turbine in each of the same periods, ",
                "namely $\\{1+(n-1)\\rho_w\\}/n$ with the correlation ",
                "$\\rho_w = ", fmt5(rho_w), "$ measured on these data, and ",
                "the second column is $n$ divided by it, that is, the number ",
                "of independent readings that one period carries. Six wind ",
                "turbines read together therefore carry the information of ",
                fmt4(e_emp[n_wt]), ", and no number of them carries more ",
                "than ", fmt4(1 / rho_w), ". Under the common wind of (16), ",
                "for which $\\rho_w = ", fmt5(rho_w_th), "$, the same two ",
                "numbers are ", fmt4(e_th[n_wt]), " and ", fmt4(1 / rho_w_th),
                "."))

## ---- splice them into the section --------------------------------
say("")
for (nm in c("table_example_sample", "table_example_estimates",
             "table_example_variance", "table_example_design"))
  splice_table(tag(nm), file.path(tex_dir, "section_example.tex"))

saveRDS(list(N = N, M = M, C = C, Csum = Csum, nj = nj, Ni = Ni,
             m0 = m0, d_rec = d_rec, d_ml = d_ml, se_d = se_d,
             p_rec = p_rec, p_wml = p_wml, pM_wml = pM_wml, p_til = p_til,
             rho_w = rho_w, rho_w_th = rho_w_th,
             rho_P = rho_P, rho_P_th = rho_P_th,
             A = c(A_FO = A_FO, A_FOTS = A_FOTS, A_sys = A_sys, A_b32 = A_b32),
             k = k_ml, alpha = a_ml, m1 = m1, m2 = m2),
        file.path("data", tag2 <- paste0(tolower(SITE), "_example.rds")))
say("objects written to data/", tag2)
close_report()
