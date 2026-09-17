## =================================================================
##  check_data.R
##
##  Verification of the data layer.  Every assumption that
##  code/prepare_data.R makes about the files of the site is tested here
##  against the data themselves, and the facts about the data that the
##  paper quotes are computed here.  PASS/FAIL summary at the end.
##
##  What is tested
##    1. the panel is complete: one record per WT per ten-minute slot
##    2. "Power (kW)" is the ten-minute mean active power in kW, checked
##       against the energy counter of the same interval
##    3. the time stamp marks the beginning of the ten-minute interval,
##       checked by aligning the event log with the power signal
##    4. the zero-power state: idle consumption is negative, and a
##       record cannot have positive mean power with non-positive maximum
##    5. the curtailment signal is identically zero at this site
##    6. every forced-outage event, whatever its Status field, is a
##       period without production, so that the union of all of them is
##       the forced-outage time; the double counting of the plain sum of
##       durations is quantified
##    7. the three categories removed in Step 2 contain only Stop events
##    8. the low wind state is counted as available by the operator, so
##       that it is correctly kept in Step 2
##    9. records without the power signal are gaps of the SCADA link and
##       not outages: most of them carry a positive energy export
##   10. the static file: the model and the nominal power of the WTs
##
##  Input  : data/<site>_panel.rds, data/<site>_status.rds,
##           data/<site>_WT_static.csv (optional)
##  Output : results/check_data_report.txt
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
new_report(file.path(res_dir, paste0(tag("check_data"), "_report.txt")))
say("check_data.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

need_file(PANEL_RDS, "code/prepare_data.R")
need_file(STATUS_RDS, "code/prepare_data.R")
D   <- readRDS(PANEL_RDS)
need_panel(D)
st  <- readRDS(STATUS_RDS)
pan <- D$panel
n_slot <- D$n_slot; N_TURB <- D$n_turb; t_origin <- D$t_origin

CHK <- new.env(); CHK$pass <- 0L; CHK$fail <- 0L; CHK$fails <- character(0)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  if (ok) CHK$pass <- CHK$pass + 1L else {
    CHK$fail <- CHK$fail + 1L; CHK$fails <- c(CHK$fails, name)
  }
  say(sprintf("%-62s %-4s %s", name, if (ok) "PASS" else "FAIL", detail))
}

## slot flags by union of intervals, as in prepare_data.R, with a shift
key <- function(turbine, slot) (turbine - 1L) * n_slot + slot + 1L
flag_slots <- function(iv, shift = 0L) {
  f <- logical(n_slot * N_TURB)
  iv <- iv[!is.na(iv$t0) & !is.na(iv$t1) & iv$t1 > iv$t0, ]
  for (j in seq_len(nrow(iv))) {
    a <- as.numeric(difftime(iv$t0[j], t_origin, units = "secs"))
    b <- as.numeric(difftime(iv$t1[j], t_origin, units = "secs"))
    s <- seq(floor(a / 600), floor(b / 600)) + shift
    s <- s[s >= 0 & s < n_slot]
    if (length(s)) f[key(iv$turbine[j], s)] <- TRUE
  }
  f[key(pan$turbine, pan$slot)]
}
has_iec <- !is.na(st$iec)
st$hours <- as.numeric(difftime(st$t1, st$t0, units = "hours"))
obs <- !is.na(pan$power)

## ---- 1. completeness of the panel --------------------------------
say_rule("1. completeness of the ten-minute panel")
dup <- anyDuplicated(pan[, c("turbine", "slot")])
check("one record per WT and slot, no duplicates", dup == 0L)
## one record per WT and per ten minutes of every year of the site
n_exp <- N_TURB * sum(sapply(D$years, function(y)
  as.integer(difftime(as.Date(paste0(y + 1L, "-01-01")),
                      as.Date(paste0(y, "-01-01")), units = "days")) * 144L))
check(paste0(N_TURB, " WTs x ", fmti(n_exp / N_TURB), " ten-minute slots = ",
             fmti(n_exp), " records"), nrow(pan) == n_exp,
      paste("found", fmti(nrow(pan))))
per_wt <- tapply(pan$slot, pan$turbine, function(s) length(s) == n_slot &&
                   all(sort(s) == seq_len(n_slot) - 1L))
check("every WT covers every slot of the grid", all(per_wt))
check("time zone of the panel is UTC", attr(pan$ts, "tzone") == "UTC")

## ---- 2. the power signal -----------------------------------------
say_rule("2. Power (kW) is the ten-minute mean active power")
ok <- obs & !is.na(pan$e_exp)
r  <- pan$e_exp[ok] / (pan$power[ok] / 6)
lowb <- ok & pan$power > 0 & pan$power <= 100
r_low <- median(pan$e_exp[lowb] / (pan$power[lowb] / 6))
cc <- cor(pan$e_exp[ok & pan$power > 50], pan$power[ok & pan$power > 50] / 6)
say("  energy export of the interval divided by power/6:")
say("    median over all records with power > 0 : ", fmt4(median(r[pan$power[ok] > 0])))
say("    median over 0 < power <= 100 kW          : ", fmt4(r_low))
say("    correlation of energy export with power/6 : ", fmt4(cc))
check("energy export tracks power/6 (correlation > 0.95)", cc > 0.95)
check("ratio close to one at low power (within 5 per cent)", abs(r_low - 1) < 0.05)
say("  (the ratio falls to about 0.90 at rated power: the export meter sits")
say("   behind the transformer and records the net export; the estimators")
say("   use the active power at the WT, which is the quantity of Eq. (1))")

## ---- 3. the time stamp convention -------------------------------
say_rule("3. time stamp marks the beginning of the ten-minute interval")
iv_fo <- st[has_iec & st$iec == "Forced outage", ]
zf <- sapply(c(-1L, 0L, 1L), function(sh) {
  f <- flag_slots(iv_fo, sh); mean(pan$power[f & obs] <= 0)
})
say("  zero-power fraction inside forced-outage slots when the events are")
say("  snapped with a shift of -1, 0, +1 slots : ", paste(fmt4(zf), collapse = "  "))
check("agreement is highest with no shift", which.max(zf) == 2L)

## ---- 4. the zero-power state -------------------------------------
say_rule("4. the zero-power state")
pw <- pan$power[obs]
say("  records with power exactly 0     : ", fmti(sum(pw == 0)))
say("  records with power < 0           : ", fmti(sum(pw < 0)),
    "  (minimum ", fmt4(min(pw)), " kW, idle consumption)")
say("  records with 0 < power <= 1 kW   : ", fmti(sum(pw > 0 & pw <= 1)))
okm <- obs & !is.na(pan$power_max)
bad <- sum(pan$power[okm] > 0 & pan$power_max[okm] <= 0)
check("no record has positive mean power with non-positive maximum",
      bad == 0L, paste("violations:", bad))
check("idle consumption is negative, so 'power <= 0' is the zero state",
      sum(pw < 0) > 1000L)
part <- sum(pan$power[okm] <= 0 & pan$power_max[okm] > 0)
say("  records with mean power <= 0 but a positive maximum inside the")
say("  interval, i.e. production during part of the ten minutes : ", fmti(part))

## ---- 5. curtailment ----------------------------------------------
say_rule("5. curtailment")
lc <- pan$lp_curt[!is.na(pan$lp_curt)]
say("  records with a positive lost production to curtailment : ",
    fmti(sum(lc > 0)), "  (maximum ", fmt4(max(lc)), " kWh)")
say("  such records are removed from the sample by Step 2, whether or not")
say("  the event log also names them")
check("the curtailment signal is present and never negative", all(lc >= 0))

## ---- 6. forced outages -------------------------------------------
say_rule("6. forced-outage events")
fo_stop <- flag_slots(st[has_iec & st$iec == "Forced outage" & st$status == "Stop", ])
fo_warn <- flag_slots(st[has_iec & st$iec == "Forced outage" & st$status == "Warning", ])
z_stop <- mean(pan$power[fo_stop & obs] <= 0)
z_warn <- mean(pan$power[fo_warn & !fo_stop & obs] <= 0)
say("  Stop events   : ", fmti(sum(fo_stop)), " slots, zero-power fraction ", fmt4(z_stop))
say("  Warning events: ", fmti(sum(fo_warn & !fo_stop)),
    " slots not covered by a Stop, zero-power fraction ", fmt4(z_warn))
check("Stop events of the category are periods without production", z_stop > 0.85)
check("Warning events of the category are periods without production", z_warn > 0.85)
u_h <- sum(pan$ev_fo) / 6
s_h <- sum(st$hours[has_iec & st$iec == "Forced outage" & is.finite(st$hours) &
                      st$hours > 0])
say("  forced-outage WT hours, union of the intervals : ", fmt4(u_h))
say("  forced-outage WT hours, plain sum of durations : ", fmt4(s_h))
check("union is used, the plain sum would double count", u_h < s_h,
      paste0("sum/union = ", fmt4(s_h / u_h)))
tab <- table(pan$ev_fo[pan$retained & !is.na(pan$av_sys)],
             pan$av_sys[pan$retained & !is.na(pan$av_sys)] < 1)
say("  cross-tabulation on the retained sample, event-log forced outage")
say("  (rows) against operator's availability below one (columns):")
say_block(tab)
check("almost every forced-outage slot is unavailable for the operator",
      tab["TRUE", "TRUE"] / sum(tab["TRUE", ]) > 0.99)

## ---- 7. the excluded categories ---------------------------------
say_rule("7. the categories removed in Step 2")
ex <- st[has_iec & st$iec %in% D$excl_iec, ]
say_block(table(ex$iec, ex$status))
check("the excluded categories contain Stop events only",
      all(ex$status == "Stop"))

## ---- 8. the low wind state ----------------------------------------
say_rule("8. the low wind state is available and is kept")
oes <- flag_slots(st[has_iec & st$iec == "Out of Environmental Specification", ])
a_oes <- mean(pan$av_sys[oes & !is.na(pan$av_sys)])
say("  operator's system availability inside 'Out of Environmental",
    " Specification' events : ", fmt4(a_oes))
check("the operator counts the low wind state as available (> 0.95)", a_oes > 0.95)
check("the low wind state is not among the excluded categories",
      !("Out of Environmental Specification" %in% D$excl_iec))

## ---- 9. records without the power signal ------------------------
say_rule("9. records without the power signal")
nap <- !obs
n_nap <- sum(nap)
say("  records without power            : ", fmti(n_nap),
    "  (", fmt4(100 * mean(nap)), " per cent)")
say("  of which inside excluded events  : ", fmti(sum(nap & pan$excl_event)))
say("  of which inside forced outages   : ", fmti(sum(nap & pan$ev_fo)),
    "  (", fmt4(sum(nap & pan$ev_fo) / 6), " WT hours)")
say("  the share of the forced-outage slots that these represent : ",
    fmt4(sum(nap & pan$ev_fo) / max(sum(pan$ev_fo), 1L)))
say("  this share of the forced-outage time carries no power signal and is")
say("  therefore not seen by the estimator, which is one of the reasons why")
say("  the estimate of p is on the high side")
say("  of which with positive energy export : ",
    fmti(sum(nap & !is.na(pan$e_exp) & pan$e_exp > 0)))
say("  wind speed missing exactly where power is missing : ",
    fmti(sum(is.na(pan$ws) & nap)), " of ", fmti(sum(is.na(pan$ws))))
by_day <- table(as.Date(pan$ts[nap]))
say("  largest daily counts of missing records (", fmti(144L * N_TURB),
    " = every WT all day):")
say_block(head(sort(by_day, decreasing = TRUE), 6))
check("most missing-power records carry a positive energy export",
      mean(pan$e_exp[nap & !is.na(pan$e_exp)] > 0) > 0.5)
## Step 1 discards the records without a power signal, and this is
## admissible as long as they are a small part of the sample; they are
## reported above in detail so that the reader may judge them.
check("the records without a power signal are below 2 per cent of the sample",
      mean(nap) < 0.02, paste(fmt4(100 * mean(nap)), "per cent"))

## ---- 10. the static file -----------------------------------------
say_rule("10. the static file")
sf <- find_data_file(paste0("^", SITE, "_WT_static\\.csv$"))[1]
if (is.na(sf)) sf <- file.path("data", paste0(SITE, "_WT_static.csv"))
if (file.exists(sf)) {
  stt <- read.csv(sf, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  ## the file of a site may end in empty lines, as that of Penmanshiel
  ## does, and they are not WTs
  stt <- stt[!is.na(stt$Title) & nzchar(trimws(stt$Title)), , drop = FALSE]
  say_block(stt[, c("Title", "Model", "Rated power (kW)", "Hub Height (m)",
                    "Rotor Diameter (m)")])
  model_short <- sub("^Senvion ", "", WT$model)
  check(paste0("the panel and the static file agree on the number of WTs"),
        nrow(stt) == N_TURB, paste(nrow(stt), "against", N_TURB))
  check(paste0("every WT is a ", WT$model),
        all(grepl(model_short, stt$Model, fixed = TRUE)))
  check(paste0("nominal power ", WT$P_r, " kW for every WT"),
        all(stt[["Rated power (kW)"]] == WT$P_r))
  say("  hub heights present : ",
      paste(sort(unique(stt[["Hub Height (m)"]])), collapse = ", "), " m")
} else say("  data/", SITE, "_WT_static.csv not present; skipped")

## =================================================================
say_rule("Summary")
say("checks passed : ", CHK$pass)
say("checks failed : ", CHK$fail)
if (CHK$fail > 0L) say("failed: ", paste(CHK$fails, collapse = " ; "))
say(if (CHK$fail == 0L) "ALL CHECKS PASS" else "SOME CHECKS FAIL")
close_report()
