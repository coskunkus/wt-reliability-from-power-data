## =================================================================
##  run_all.R
##
##  Everything in the real data application of the paper except the
##  Monte Carlo study of the sampling design, which is in
##  code/run_simulation.R.
##
##  Input   : data/kelmarsh_panel.rds, data/kelmarsh_status.rds
##            (written by code/prepare_data.R)
##  Output  : results/run_all_report.txt, which contains every number
##            of the section and the LaTeX of its single table
##            "Claude outputs"/table2_quantities.tex, the table of the
##            section, and table_wind_fit.tex, which is not in it
##            figures/figure_cut_in.{eps,pdf} and
##            figures/figure_wind_fit.{eps,pdf}, neither of them in the
##            paper: they are the evidence for two of its statements
##
##  Nothing outside the paper is used: the WT characteristics are the
##  manufacturer's (Section 2), the wind speed distribution is Weibull
##  with k and alpha estimated by maximum likelihood (Sections 3 and 5),
##  the Weibull being the family that the Akaike information criterion
##  selects out of five families in common use for wind speed,
##  and p is estimated by Eqs. (8) and (12) with p_M = min(p_hat, 1).
##
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
for (dd in c(res_dir, fig_dir, tex_dir))
  if (!dir.exists(dd)) dir.create(dd, recursive = TRUE)

new_report(file.path(res_dir, paste0(tag("run_all"), "_report.txt")))
say("run_all.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

need_file(PANEL_RDS, "code/prepare_data.R")
need_file(STATUS_RDS, "code/prepare_data.R")
D   <- readRDS(PANEL_RDS)
need_panel(D)
st  <- readRDS(STATUS_RDS)
pan <- D$panel
n_turb <- D$n_turb

## The characteristics are the manufacturer's, as listed by Bauer and
## Matysik, wind-turbine-models.com: the Senvion MM92 of Kelmarsh has a
## cut-in speed of 3.0 m/s, a rated speed of 12.5 m/s, a cut-out speed
## of 24.0 m/s and a nominal power of 2050 kW, and the Senvion MM82 of
## Penmanshiel has 3.5, 14.5, 25.0 m/s and 2050 kW.  They are held in
## WT_SITES of the core.
P_r  <- WT$P_r          # nominal power, in kW
v_ci <- WT$v_ci         # manufacturer's cut-in wind speed, in m/s
v_r  <- WT$v_r          # manufacturer's rated wind speed, in m/s
v_co <- WT$v_co         # manufacturer's cut-out wind speed, in m/s

## The sample and the zero-power indicator are those of the core, which
## defines them once for every script that estimates p.
ret  <- estimation_sample(pan)
sub  <- ret

## =================================================================
##  1. The data and the filtering
## =================================================================
say_rule("1. Data and filtering")
say("site                                 : ", D$site, ", ", n_turb,
    " x ", WT$model, ", P_r = ", P_r, " kW")
say("period                               : ", format(min(pan$ts)), " to ",
    format(max(pan$ts)))
say("ten-minute WT records                : ", fmti(nrow(pan)))
say("  with power and wind speed signal   : ", fmti(sum(pan$observed)),
    "  (", fmt4(100 * mean(pan$observed)), " per cent)")
say("  without one of the two, Step 1     : ", fmti(sum(!pan$observed)),
    "  (", fmt4(100 * mean(!pan$observed)), " per cent)")
say("  removed by the event log, Step 2   : ", fmti(sum(pan$excl_event)),
    "  (", fmt4(100 * mean(pan$excl_event)), " per cent)")
say("  retained                           : ", fmti(sum(pan$retained)),
    "  (", fmt4(100 * mean(pan$retained)), " per cent)")
say("sample size n                        : ", fmti(nrow(sub)))
say("(instants at which every WT is retained, used in Step 3 only : ",
    fmti(sum(table(ret$slot) == n_turb)), ")")

## =================================================================
##  2. The WT characteristics
## =================================================================
say_rule("2. WT characteristics, manufacturer's values (Section 2 of the paper)")
say("nominal power        P_r  : ", fmt4(P_r), " kW")
say("cut-in wind speed    v_ci : ", fmt4(v_ci), " m/s")
say("rated wind speed     v_r  : ", fmt4(v_r), " m/s")
say("cut-out wind speed   v_co : ", fmt4(v_co), " m/s")
## the records of an operating WT, for the figure of the cut-in
avl <- sub[sub$av_sys >= 0.999, ]

## =================================================================
##  3. The wind speed distribution
## =================================================================
say_rule("3. Wind speed distribution")
## The paper takes F to be Weibull.  That choice is not assumed here:
## five families in common use for wind speed are fitted to the same
## record by maximum likelihood and compared by the Akaike information
## criterion.  Beside the criterion the table gives the only two
## features of F that the estimators use, F(v_ci) and d, against the
## proportions of the wind speed record itself.  The comparison is not
## in the section, which takes F to be Weibull as Section 3 of the
## paper does; it is kept here, with its table, against the question
## why that family and no other.
FAM   <- fit_wind_families(sub$ws, v_ci, v_co)
d_rec <- mean(sub$ws >= v_ci & sub$ws < v_co)
F_rec <- mean(sub$ws < v_ci)
fam   <- FAM$tab

say("five families fitted by maximum likelihood on ", fmti(FAM$n),
    " wind speed records")
one <- function(x) formatC(x, format = "f", digits = 1)
say_block(data.frame(family = fam$family, npar = fam$npar,
                     logLik = one(fam$logLik), AIC = one(fam$AIC),
                     dAIC = one(fam$dAIC),
                     F_vci = fmt5(fam$F_vci), d = fmt5(fam$d),
                     err_d = fmt5(fam$d - d_rec)))
say("")
for (nm in fam$family) {
  z <- FAM$fits[[nm]]
  say("  ", nm, " : ",
      paste(names(z$f$pars(z$fit$par)), fmt4(z$f$pars(z$fit$par)),
            sep = " = ", collapse = " ,  "))
}
say("")
say("smallest AIC                   : ", FAM$best, ", smaller by ",
    fmt4(fam$dAIC[2]), " than the next")
say("the family of the paper is the one with the smallest AIC : ",
    identical(FAM$best, "Weibull"))
stopifnot(identical(FAM$best, "Weibull"))
say("from the wind speed record, without a family")
say("  F(v_ci) = P{V < ", v_ci, "}           : ", fmt5(F_rec))
say("  d       = P{", v_ci, " <= V < ", v_co, "}  : ", fmt5(d_rec))
say("the family whose d is closest to the record : ",
    fam$family[which.min(abs(fam$d - d_rec))])
say("")

## Weibull with k and alpha estimated by maximum likelihood, as the paper
## prescribes; the historical wind record is the retained sample itself
ws <- sub$ws[sub$ws > 0]
ml  <- fit_weibull_ml(ws)
W_ml <- wind_weibull(ml[["k"]], ml[["c"]])
d_ml <- d_wt(W_ml, v_ci, v_co)
Eg1  <- Eg_wt(W_ml, v_ci, v_r, v_co, P_r, 1)
## the two routes to the Weibull agree to the last reported digit
stopifnot(abs(d_ml - fam$d[fam$family == "Weibull"]) < 1e-8)

say("Weibull, maximum likelihood    : k = ", fmt4(ml[["k"]]),
    " , alpha = ", fmt4(ml[["c"]]), " m/s  (convergence code ",
    ml[["convergence"]], ")")
say("d = F(v_co) - F(v_ci)          : ", fmt5(d_ml))
say("F(v_ci) = P{V < ", v_ci, "}           : ", fmt4(W_ml$F(v_ci)))
say("mu(1), the denominator of Eq. (8), Weibull ML : ", fmt4(Eg1), " kW")
say("observed mean power output       : ", fmt4(mean(sub$power)), " kW")

## =================================================================
##  4. The benchmark from the failure and maintenance records
## =================================================================
say_rule("4. Availability from the failure and maintenance records")
## Two families of indicators, both over the retained sample, so that
## they are comparable with the estimators.
##  (i)  from the event log directly, by the union of the event
##       intervals of a category on the ten-minute grid:
##       A_FO    charges the forced outages only, which is the state
##               that 1 - p describes in Section 2;
##       A_FOTS  charges in addition the technical standby stops, the
##               other state in which a WT with sufficient wind does not
##               produce for a reason internal to the WT.
##  (ii) the operator's own indicators, as computed by the SCADA system
##       from the same event log: the time-based system availability and
##       the IEC 61400-26 B.3.2 availability.
A_FO   <- 1 - mean(sub$ev_fo)
A_FOTS <- 1 - mean(sub$ev_fo | sub$ev_ts)
A_sys  <- mean(sub$av_sys, na.rm = TRUE)
A_b32  <- mean(sub$av_b32, na.rm = TRUE)

say("from the event log, union of the event intervals, retained sample")
say("  1 - P{forced outage} = 1 - FOR                         : ", fmt5(A_FO))
say("  1 - P{forced outage or technical standby stop} = FOTS   : ", fmt5(A_FOTS))
say("operator's indicators, retained sample")
say("  operator's time-based system availability              : ", fmt5(A_sys))
say("  operator's time-based IEC 61400-26 B.3.2 availability  : ", fmt5(A_b32))

## =================================================================
##  5. The estimators on the full retained sample
## =================================================================
## Section 2 assumes, through g of (2), that a WT produces no power
## when the wind speed lies outside [v_ci, v_co).  The sample is made
## to conform to that assumption: a record whose wind speed lies
## outside the interval is counted as a zero-power observation,
## whatever its measured power.  The number of records that this
## affects is reported, since they are records in which the WT did
## produce.
say_rule("5. Estimates of p, Eqs. (8) and (12) of the paper")
N      <- nrow(sub)
inband <- in_band(sub$ws, v_ci, v_co)
m0_obs <- sum(sub$zero)                       # zero as measured
n_adj  <- sum(!inband & !sub$zero)            # produced outside the interval
m0     <- sum(zero_of(sub, v_ci, v_co))       # zero under the assumption
pw     <- ifelse(inband, sub$power, 0)        # the power output of the model

## d_rec, the value read off the wind speed record of the site, was
## obtained in part 3 above
stopifnot(abs(d_rec - mean(inband)) < 1e-12)
p_rec  <- p_hat_wt(m0, N, d_rec)
p_wml  <- p_hat_wt(m0, N, d_ml) # d from the fitted distribution
pM_wml <- p_M_wt(m0, N, d_ml)
p_til  <- p_tilde_wt(pw, W_ml, v_ci, v_r, v_co, P_r)

say("n                                : ", fmti(N))
say("zero-power records as measured   : ", fmti(m0_obs), "  (", fmt5(m0_obs / N), ")")
say("records outside [v_ci, v_co)     : ", fmti(sum(!inband)), "  (", fmt5(mean(!inband)), ")")
say("  below v_ci                     : ", fmti(sum(sub$ws < v_ci)),
    "  , of which the WT produced : ", fmti(sum(sub$ws < v_ci & !sub$zero)))
say("  at or above v_co               : ", fmti(sum(sub$ws >= v_co)),
    "  , of which the WT produced : ", fmti(sum(sub$ws >= v_co & !sub$zero)))
say("  produced outside the interval  : ", fmti(n_adj), "  (", fmt5(n_adj / N), ")")
say("m0 under the assumption of (2)   : ", fmti(m0), "  ,  m0/n = ", fmt5(m0 / N))
say("")
say("d from the wind speed record     : ", fmt5(d_rec))
## Section 2 takes F, and with it d, to be known.  It is not known
## here, and the proportion above is put in its place; the standard
## error of that proportion says what the substitution costs.
say("  standard error of that proportion : ",
    formatC(sqrt(d_rec * (1 - d_rec) / N), format = "f", digits = 5))
say("(12) with that d                 : ", fmt5(p_rec))
say("d from the fitted Weibull        : ", fmt5(d_ml))
say("(12) with that d                 : ", fmt5(p_wml), "  ->  p_M = min(p_hat, 1) = ", fmt5(pM_wml))
say("(8)  p_tilde                     : ", fmt5(p_til), "   (mu(1) = ", fmt4(Eg1), " kW)")
say("")
say("with d of the wind speed record, (12) is the proportion of the records")
say("of the interval in which the WT produced, and cannot exceed one:")
say("  proportion                     : ", fmt5(mean(pw[inband] > 0)))
say("")
say("estimate minus the indicators, d from the wind speed record:")
say("  1 - FOR : ", formatC(p_rec - A_FO,   format = "f", digits = 6),
    "   FO or TS : ", formatC(p_rec - A_FOTS, format = "f", digits = 6),
    "   operator system : ", formatC(p_rec - A_sys, format = "f", digits = 6),
    "   operator B3.2 : ", formatC(p_rec - A_b32, format = "f", digits = 6))
say("largest |estimate - indicator|   : ",
    fmt4(max(abs(p_rec - c(A_FO, A_FOTS, A_sys, A_b32)))))
say("")
say("the assumption matters because the WTs do produce below v_ci; the")
say("probability that an operating WT produces, by wind speed class, is")
say("drawn in the figure of this section")

## What the estimate is made of: the records of the interval in which
## the WT did not produce, which are the only ones that 1 - p_hat
## counts, sorted by what the event log says of them.
ib  <- sub[inband, ]
zib <- ib[ib$zero, ]
say("")
say("the records of the interval in which the WT did not produce")
say("  records of the interval        : ", fmti(nrow(ib)))
say("  of them with zero power        : ", fmti(nrow(zib)),
    "  (", fmt5(nrow(zib) / nrow(ib)), ")")
say("    inside a forced outage       : ", fmti(sum(zib$ev_fo)))
say("    otherwise                    : ", fmti(sum(!zib$ev_fo)))
say("    otherwise, and below v_ci + 1 m/s : ",
    fmti(sum(!zib$ev_fo & zib$ws < v_ci + 1)),
    "  (", fmt4(mean(zib$ws[!zib$ev_fo] < v_ci + 1)), " of the rest)")
say("  the last group is worth, in the estimate : ",
    fmt5(sum(!zib$ev_fo & zib$ws < v_ci + 1) / nrow(ib)))
say("  they are the records of an operating WT that has not yet started,")
say("  the transition of the figure being smooth rather than a step")

## =================================================================
##  5b. One WT and one year at a time
## =================================================================
## The estimator is applied to every WT and year separately, the wind
## speed record of that same WT and year supplying d, so that no WT
## year is singled out.  The whole table is printed here; the section
## quotes the WT year whose estimate is the smallest, which is the one
## in which the WT lost most of its time to failures.
say_rule("5b. The estimator applied to one WT and one year at a time")
ret$year <- as.integer(format(ret$ts, "%Y"))
cell <- do.call(rbind, lapply(split(ret, list(ret$turbine, ret$year)),
  function(z) {
    if (!nrow(z)) return(NULL)
    inb <- in_band(z$ws, v_ci, v_co)
    data.frame(turbine = z$turbine[1], year = z$year[1], n = nrow(z),
               d = mean(inb), p_hat = p_hat_wt(sum(zero_of(z, v_ci, v_co)),
                                               nrow(z), mean(inb)),
               A_FO = 1 - mean(z$ev_fo), A_FOTS = 1 - mean(z$ev_fo | z$ev_ts),
               A_sys = mean(z$av_sys, na.rm = TRUE),
               A_b32 = mean(z$av_b32, na.rm = TRUE))
  }))
cell <- cell[order(cell$p_hat), ]
row.names(cell) <- NULL
say_block(data.frame(turbine = cell$turbine, year = cell$year,
                     n = fmti(cell$n), d = fmt5(cell$d),
                     p_hat = fmt5(cell$p_hat), A_FO = fmt5(cell$A_FO),
                     A_FOTS = fmt5(cell$A_FOTS), A_sys = fmt5(cell$A_sys),
                     A_b32 = fmt5(cell$A_b32)))
say("")
say("WT years                             : ", nrow(cell))
for (nm in c("A_FO", "A_FOTS", "A_sys", "A_b32"))
  say("largest |estimate - ", nm, "| over the WT years : ",
      fmt5(max(abs(cell$p_hat - cell[[nm]]))))
say("the smallest estimate and the smallest value of every indicator",
    " belong to the same WT year : ",
    all(vapply(c("p_hat", "A_FO", "A_FOTS", "A_sys", "A_b32"),
               function(nm) which.min(cell[[nm]]) == which.min(cell$p_hat),
               logical(1))))

## the WT year of the section, and what the event log says of it
cs   <- cell[which.min(cell$p_hat), ]
sfo  <- st[st$turbine == cs$turbine & !is.na(st$t0) & !is.na(st$t1) &
             as.integer(format(st$t0, "%Y")) == cs$year &
             grepl("forced outage", st$iec, ignore.case = TRUE), ]
sfo$hrs <- as.numeric(difftime(sfo$t1, sfo$t0, units = "hours"))
sfo <- sfo[order(-sfo$hrs), ]
say("")
say("the WT year with the smallest estimate : WT ", cs$turbine, " in ", cs$year)
say("  n                                  : ", fmti(cs$n))
say("  d of that WT year                  : ", fmt5(cs$d))
say("  (12)                               : ", fmt5(cs$p_hat))
say("  1 - FOR                            : ", fmt5(cs$A_FO))
say("  forced outage or technical standby : ", fmt5(cs$A_FOTS))
say("  operator's system availability     : ", fmt5(cs$A_sys))
say("  operator's B.3.2 availability      : ", fmt5(cs$A_b32))
say("  forced outage events of that year  : ", nrow(sfo))
say("  hours lost to forced outages       : ",
    formatC(sum(sfo$hrs), format = "f", digits = 2))
say("  longest single forced outage       : ",
    formatC(sfo$hrs[1], format = "f", digits = 2), " hours, ",
    ## the date is written in figures, not in words: the name of the
    ## month would be that of the locale of the machine, and the report
    ## must not depend on it
    format(sfo$t0[1], "%Y-%m-%d"), " , ", sfo$message[1])

## =================================================================
##  5c. Without Step 2: every stop counted as a stop
## =================================================================
## Step 2 removes the periods of scheduled maintenance and of the other
## stops that are not failures, and it needs the event log to do so.
## The event log may be dispensed with altogether, if a WT that is
## stopped is counted as not being in the operating state whatever the
## cause of the stop.  Nothing is then removed, and the estimation uses
## nothing but the power output and the wind speed.  The benchmark is
## then the fraction of the time in which no stop of any kind is
## registered in the event log.
say_rule("5c. Without Step 2, every stop counted as a stop")
al  <- pan[pan$observed, ]
tb  <- table(al$slot)
alb <- al[as.character(al$slot) %in% names(tb)[tb == n_turb], ]
inb2 <- in_band(alb$ws, v_ci, v_co)
N2   <- nrow(alb)
m02  <- sum(zero_of(alb, v_ci, v_co))
d2   <- mean(inb2)
p2   <- p_hat_wt(m02, N2, d2)
A_any <- 1 - mean(alb$ev_fo | alb$ev_ts | alb$excl_event)
say("records with a power and a wind speed signal : ", fmti(nrow(al)))
say("instants with every WT present               : ",
    fmti(length(unique(alb$slot))))
say("n                                            : ", fmti(N2))
say("m0/n                                         : ", fmt5(m02 / N2))
say("d from the wind speed record                 : ", fmt5(d2))
say("(12)                                         : ", fmt5(p2))
say("share of the time under scheduled maintenance, out of electrical")
say("  specification or requested shutdown         : ",
    fmt5(mean(alb$excl_event)))
say("1 - P{a stop of any kind in the event log}   : ", fmt5(A_any))
say("operator's time-based system availability    : ",
    fmt5(mean(alb$av_sys, na.rm = TRUE)))
say("operator's time-based IEC 61400-26 B.3.2     : ",
    fmt5(mean(alb$av_b32, na.rm = TRUE)))
say("estimate minus that benchmark, in percentage points : ",
    formatC(100 * (p2 - A_any), format = "f", digits = 2))
say("the two differ, in percentage points         : ",
    formatC(abs(100 * (p2 - A_any)), format = "f", digits = 2))
say("")
say("the two readings differ by definition: with Step 2 the estimate is")
say("of 1 - FOR, without it of the availability of the WT")

## =================================================================
##  6. The figure of the paper
## =================================================================
## The wind speed record and the five fitted densities.  The figure is
## not in the paper, the information of it being in the table of the
## families, but it is drawn here because it is what the table asserts.
fig_device(tag("figure_wind_fit"), quote({
  hist(sub$ws[sub$ws > 0], breaks = 80, freq = FALSE, col = "grey92",
       border = "grey75", xlab = "Wind speed v (m/s)", ylab = "Density",
       main = "", xlim = c(0, ceiling(v_co / 2) * 2))
  xs <- seq(0.01, ceiling(v_co / 2) * 2, by = 0.02)
  lt <- c(1, 2, 3, 4, 5)
  cl <- c("black", "grey30", "grey45", "grey55", "grey65")
  for (i in seq_len(nrow(fam))) {
    z <- FAM$fits[[fam$family[i]]]
    lines(xs, exp(z$f$dens(xs, z$fit$par)) * (1 - FAM$pi0),
          lwd = 2, lty = lt[i], col = cl[i])
  }
  legend("topright", bty = "n", lwd = 2, lty = lt[seq_len(nrow(fam))],
         col = cl[seq_len(nrow(fam))], cex = 0.75,
         legend = paste0(fam$family, ", AIC ",
                         formatC(fam$AIC, format = "d", big.mark = ",")))
}), width = 5.6, height = 3.8)

## The transition that q describes, as it is observed, against the step
## that Eq. (2) prescribes at the nominal cut-in speed.  Neither this
## figure nor the one above is in the paper, the section reporting in
## numbers what they show; they are drawn because they are the evidence
## for two of its statements.
Qc <- q_positive(avl$ws, avl$power)$curve
fig_device(tag("figure_cut_in"), quote({
  plot(Qc$v, Qc$q, type = "n", xlim = c(1.5, 4.5), ylim = c(0, 1),
       xlab = "Wind speed v (m/s)",
       ylab = "Probability of a positive power output")
  lines(c(1.5, v_ci, v_ci, 4.5), c(0, 0, 1, 1), lwd = 2, lty = 2,
        col = "grey35")
  lines(Qc$v, Qc$q, lwd = 2)
  text(v_ci, 0.20, expression(v[ci]), pos = 4, cex = 0.95, col = "grey35")
  legend("topleft", inset = c(0.02, 0.02), bty = "n", lwd = 2,
         lty = c(1, 2), cex = 0.8, col = c("black", "grey35"),
         legend = c("Observed", "Step of Eq. (2)"))
}))

## =================================================================
##  7. The LaTeX table of the paper
## =================================================================
emit_table(
  tag("table_wind_fit"),
  body = cbind(fam$family,
               formatC(fam$npar, format = "d"),
               formatC(fam$logLik, format = "f", digits = 1, big.mark = ","),
               formatC(fam$AIC, format = "f", digits = 1, big.mark = ","),
               formatC(fam$dAIC, format = "f", digits = 1, big.mark = ","),
               fmt5(fam$F_vci), fmt5(fam$d)),
  colspec = "lcrrrrr",
  header = paste0("Family & Parameters & Log-likelihood & AIC & ",
                  "$\\Delta$AIC & $F(v_{ci})$ & $d$"),
  caption = paste0("Five distributions in common use for wind speed, ",
                   "fitted to the ", fmti(FAM$n), " wind speed records of ",
                   "the sample by maximum likelihood and ordered by the ",
                   "Akaike information criterion (", SITE, ", ",
                   min(WT$years), " to ", max(WT$years), "). The last two ",
                   "columns are the only two features of $F$ that the ",
                   "estimators of Section 2 use; the wind speed record ",
                   "itself gives $", fmt5(F_rec), "$ and $", fmt5(d_rec),
                   "$ for them."),
  label = "tab:windfit", rules = integer(0))

emit_table(
  tag("table2_quantities"),
  body = cbind(
    c("Nominal power $P_r$", "Cut-in wind speed $v_{ci}$",
      "Rated wind speed $v_r$", "Cut-out wind speed $v_{co}$",
      "Wind speed distribution $F$",
      "$d = F(v_{co}) - F(v_{ci})$ of the fitted $F$",
      "$\\mu(1)$ of Eq. (5)",
      paste0("$\\hat p$ of Eq. (12), $d = ", fmt5(d_rec),
             "$ from the wind speed record"),
      paste0("$\\hat p_M$ of Eq. (12), $d = ", fmt5(d_ml),
             "$ from the fitted $F$"),
      "$1 - FOR$, forced outage", "Forced outage or technical standby",
      "Operator's time-based system availability",
      "Operator's IEC 61400-26 B.3.2 availability"),
    c(paste0(P_r, " kW, manufacturer"),
      paste0(fmt4(v_ci), " m/s, manufacturer"),
      paste0(fmt4(v_r), " m/s, manufacturer"),
      paste0(fmt4(v_co), " m/s, manufacturer"),
      paste0("Weibull, $\\hat k = ", fmt4(ml[["k"]]), "$, $\\hat\\alpha = ",
             fmt4(ml[["c"]]), "$ m/s"),
      fmt5(d_ml), paste0(fmt4(Eg1), " kW"),
      fmt5(p_rec), fmt5(pM_wml),
      fmt5(A_FO), fmt5(A_FOTS), fmt5(A_sys), fmt5(A_b32))),
  label_row = list(at = 8L,
                   text = "\\multicolumn{2}{l}{Estimates, and the failure and maintenance records}"),
  colspec = "lr",
  header = c("\\multicolumn{2}{l}{Quantities required by the estimators} ",
             "Quantity & Value"),
  caption = paste0("The quantities that the estimators of Section 2 require, ",
                   "the two estimates of $p$, and four availability ",
                   "indicators of the same wind turbines over the same period",
                   " (", SITE, ", ", min(WT$years), " to ", max(WT$years), ")."),
  label = "tab:quantities", rules = c(7L, 8L),
  note = paste0("Note: $n = ", fmti(N), "$, $m_0 = ", fmti(m0), "$, $m_0/n = ",
                fmt5(m0 / N), "$. The parameters of the wind speed ",
                "distribution are estimated by maximum likelihood. A record ",
                "whose wind speed lies outside $[v_{ci}, v_{co})$ counts as a ",
                "zero-power observation, as $g$ of Eq. (2) prescribes; this ",
                "concerns ", fmti(n_adj), " records in which the wind turbine ",
                "did produce. Before the modification the second estimate is ",
                fmt4(p_wml), ", and the moment estimator (8) takes the ",
                "inadmissible value ", fmt4(p_til), ". The first two ",
                "indicators are one minus the fraction of the records that ",
                "fall inside an event of the category named, a technical ",
                "standby being charged only when it stops the wind turbine; ",
                "the last two are the availability signals that the ",
                "operator's own system reports, averaged over the same ",
                "records."))

say("")
say("")
splice_table(tag("table2_quantities"))
say("figure written to ", fig_dir, "/ , table written to \"", tex_dir, "\"/")

saveRDS(list(P_r = P_r, v_ci = v_ci, v_r = v_r, v_co = v_co,
             fam = fam, F_rec = F_rec, ml = ml, d_ml = d_ml, Eg1 = Eg1, N = N, m0 = m0,
             m0_obs = m0_obs, n_adj = n_adj, d_rec = d_rec, p_rec = p_rec,
             p_wml = p_wml, pM_wml = pM_wml, p_til = p_til,
             A_FO = A_FO, A_FOTS = A_FOTS, A_sys = A_sys, A_b32 = A_b32),
        file.path("data", paste0(tag("run_all_objects"), ".rds")))
close_report()
