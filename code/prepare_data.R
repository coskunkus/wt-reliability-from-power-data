## =================================================================
##  prepare_data.R
##
##  Steps 1 and 2 of the estimation protocol of Section 2.2.
##
##  Reads the open wind farm data of Plumley (2022; Zenodo
##  record 8252025, CC-BY-4.0) and writes the ten-minute panel that
##  every other script uses.  This is the only script that touches the
##  raw archives, and it is run once.
##
##  Expected input, either of the two:
##    data/Kelmarsh_SCADA_2017_3083.zip   (and 2018, 2019), or
##    data/raw/2017/Turbine_Data_Kelmarsh_1_....csv  and the rest
##
##  The archives hold about 2.3 GB of csv files, which are of no further
##  use once the panel has been built.  They are therefore unpacked one
##  year at a time into a scratch folder under tempdir(), read, and
##  deleted again, so that nothing large is written into the project
##  folder.  A project folder that lives on a synchronised drive would
##  otherwise upload 2.3 GB of intermediate csv files.  Set
##  KEEP_RAW <- TRUE below to unpack into data/raw/ and keep the files,
##  and an existing data/raw/<year>/ is always used as it stands.
##
##  Output:
##    data/<site>_panel.rds
##    data/<site>_status.rds
##    results/prepare_data_report.txt
##
##  Base R only.  Run time is about four minutes, peak scratch use
##  about 900 MB.
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

new_report(file.path(res_dir, paste0(tag("prepare_data"), "_report.txt")))
say("prepare_data.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

## ---- the site -----------------------------------------------------
## SITE, WT and the two file names come from the core; the number of WTs
## is not fixed here but read off the archives, since the numbering of a
## site may have gaps.
DOI      <- WT$doi
YEARS    <- as.integer(WT$years)
SLOT_SEC <- 600L

## The columns retained from the 299 columns of a Turbine_Data file.
## The four availability signals are Greenbyte's IEC 61400-26 indicators
## and are used only to build the benchmark of the application.
SCADA_WANTED <- c("# Date and time",
                  "Wind speed (m/s)",
                  "Power (kW)",
                  "Power, Maximum (kW)",
                  "Energy Export (kWh)",
                  "Lost Production to Curtailment (Total) (kWh)",
                  "Time-based System Avail.",
                  "Time-based IEC B.3.2 (Manufacturers View)")
SCADA_SHORT  <- c("ts", "ws", "power", "power_max", "e_exp", "lp_curt",
                  "av_sys", "av_b32")

## Step 2: the IEC 61400-26 categories whose time stamps are removed,
## because a zero power output in these periods is not caused by a
## failure of the WT.  The low wind state, which the model itself
## describes through d, is not among them.
EXCL_IEC <- c("Scheduled Maintenance",
              "Out of Electrical Specification",
              "Requested Shutdown")

KEEP_RAW <- FALSE     # TRUE unpacks into data/raw/ and keeps the csv files
## year_folder() itself is in the core, so that every script that reads
## the archives of the site unpacks them in the same way.

## =================================================================
##  Reading, one year at a time
## =================================================================
turbine_of <- function(f)
  as.integer(sub(paste0(".*_", SITE, "_([0-9]+)_.*"), "\\1", basename(f)))

read_scada <- function(f) {
  hdr <- scan(f, what = "", sep = ",", skip = 9, nlines = 1, quiet = TRUE)
  cc  <- rep("NULL", length(hdr))
  idx <- match(SCADA_WANTED, hdr)
  if (anyNA(idx))
    stop("\n\n  The following signal is not carried by ", basename(f), ":\n    ",
         paste(SCADA_WANTED[is.na(idx)], collapse = "\n    "),
         "\n  Not every signal is available over the whole period of this data",
         "\n  set. Prepare a later year, by setting years in WT_SITES of",
         "\n  code/core_wtrel.R, and run this script again.\n", call. = FALSE)
  cc[idx] <- c("character", rep("numeric", length(SCADA_WANTED) - 1L))
  d <- read.csv(f, skip = 9, header = TRUE, check.names = FALSE,
                colClasses = cc, na.strings = c("NaN", "NA", ""))
  d <- d[, SCADA_WANTED, drop = FALSE]
  names(d) <- SCADA_SHORT
  d$ts <- as.POSIXct(d$ts, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  d$turbine <- turbine_of(f)
  d
}

read_status <- function(f) {
  d <- read.csv(f, skip = 9, header = TRUE, check.names = FALSE,
                colClasses = "character", na.strings = c("NaN", "NA", ""))
  keep <- c("Timestamp start", "Timestamp end", "Status", "Code",
            "Message", "Service contract category", "IEC category")
  d <- d[, keep, drop = FALSE]
  names(d) <- c("t0", "t1", "status", "code", "message", "svc", "iec")
  to_time <- function(x) {
    x <- trimws(x); x[x == "" | x == "NaN" | x == "NA"] <- NA_character_
    as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  d$t0 <- to_time(d$t0); d$t1 <- to_time(d$t1)
  d$turbine <- turbine_of(f)
  d
}

pan_list <- vector("list", length(YEARS))
st_list  <- vector("list", length(YEARS))
for (i in seq_along(YEARS)) {
  y  <- YEARS[i]
  yf <- year_folder(y, KEEP_RAW)
  fs <- sort(list.files(yf$dir, paste0("^Turbine_Data_", SITE, "_.*\\.csv$"),
                        full.names = TRUE, recursive = TRUE))
  ss <- sort(list.files(yf$dir, paste0("^Status_", SITE, "_.*\\.csv$"),
                        full.names = TRUE, recursive = TRUE))
  say(y, ": ", length(fs), " SCADA files and ", length(ss), " status files")
  if (!length(fs))
    stop("\n\n  No SCADA file was found for ", y,
         ".\n  The archive of ", y, " seems to be incomplete.\n",
         call. = FALSE)
  pan_list[[i]] <- do.call(rbind, lapply(fs, read_scada))
  st_list[[i]]  <- do.call(rbind, lapply(ss, read_status))
  ## the csv files are of no further use, so the scratch copy goes now
  if (yf$scratch) unlink(yf$dir, recursive = TRUE)
}

pan <- do.call(rbind, pan_list)
st  <- do.call(rbind, st_list)
rm(pan_list, st_list)

## The WTs are numbered as the files name them, and the numbering may
## have a gap, as it has at Penmanshiel where there is no WT 3.  They
## are therefore reindexed from one, the file numbers being kept for the
## report; at a site whose numbering is complete nothing changes.
TURB_ID <- sort(unique(pan$turbine))
N_TURB  <- length(TURB_ID)
pan$turbine <- match(pan$turbine, TURB_ID)
st$turbine  <- match(st$turbine,  TURB_ID)
st <- st[!is.na(st$turbine), ]

pan <- pan[order(pan$turbine, pan$ts), ]
rownames(pan) <- NULL

## ---- the ten-minute grid ----------------------------------------
t_origin  <- min(pan$ts)
pan$slot  <- as.integer(as.numeric(difftime(pan$ts, t_origin, units = "secs")) /
                          SLOT_SEC)
n_slot    <- max(pan$slot) + 1L

## =================================================================
##  Step 2: the filtering, and the event flags of the benchmark
## =================================================================
## Every event interval is snapped onto the ten-minute grid and the
## slots are flagged by the UNION of the intervals of a category, so
## that overlapping events are never counted twice.  The time stamp of
## a SCADA record marks the beginning of its ten-minute interval, which
## is what the slot arithmetic below assumes; code/check_data.R tests
## this convention against the event log.
key <- function(turbine, slot) (turbine - 1L) * n_slot + slot + 1L
flag_slots <- function(iv) {
  f <- logical(n_slot * N_TURB)
  iv <- iv[!is.na(iv$t0) & !is.na(iv$t1) & iv$t1 > iv$t0, ]
  for (j in seq_len(nrow(iv))) {
    a <- as.numeric(difftime(iv$t0[j], t_origin, units = "secs"))
    b <- as.numeric(difftime(iv$t1[j], t_origin, units = "secs"))
    s <- seq(floor(a / SLOT_SEC), floor(b / SLOT_SEC))
    s <- s[s >= 0 & s < n_slot]
    if (length(s)) f[key(iv$turbine[j], s)] <- TRUE
  }
  f[key(pan$turbine, pan$slot)]
}
has_iec <- !is.na(st$iec)
pan$excl_event <- flag_slots(st[has_iec & st$iec %in% EXCL_IEC, ])
## forced outages: every event of that category, whatever its Status.
## Both the Stop and the Warning events of this category correspond to
## a WT that is not producing (checked in code/check_data.R).
pan$ev_fo <- flag_slots(st[has_iec & st$iec == "Forced outage", ])
## technical standby that actually stops the WT
pan$ev_ts <- flag_slots(st[has_iec & st$iec == "Technical Standby" &
                             st$status == "Stop", ])

## Curtailment is reported explicitly by the data set.  At this site the
## signal is identically zero, so that the Requested Shutdown category
## of the event log is the only source of curtailment and is already
## contained in EXCL_IEC.
pan$excl_curt <- !is.na(pan$lp_curt) & pan$lp_curt > 0
pan$excluded  <- pan$excl_event | pan$excl_curt
pan$observed  <- !is.na(pan$power) & !is.na(pan$ws)
pan$retained  <- pan$observed & !pan$excluded
pan$zero      <- pan$power <= 0
pan$year      <- as.integer(format(pan$ts, "%Y"))

## instants at which every WT is retained
tab <- tapply(pan$retained, pan$slot, sum)
ok  <- as.integer(names(tab))[tab == N_TURB]
pan$full_instant <- pan$slot %in% ok

## =================================================================
##  Report
## =================================================================
say_rule("Step 1 and Step 2")
say("WT numbers in the files              : ",
    paste(TURB_ID, collapse = ", "))
say("site                                 : ", SITE, ", ", N_TURB,
    " x ", WT$model)
say("period                               : ", format(min(pan$ts)), " to ",
    format(max(pan$ts)))
say("ten-minute WT records                : ", fmti(nrow(pan)))
say("  with power and wind speed signal   : ", fmti(sum(pan$observed)),
    "  (", fmt4(100 * mean(pan$observed)), " per cent)")
say("  removed by the event log, Step 2   : ", fmti(sum(pan$excl_event)),
    "  (", fmt4(100 * mean(pan$excl_event)), " per cent)")
say("  removed as curtailed               : ", fmti(sum(pan$excl_curt)))
say("  retained                           : ", fmti(sum(pan$retained)),
    "  (", fmt4(100 * mean(pan$retained)), " per cent)")
say("instants with every WT retained      : ", fmti(length(ok)))
say("observations at those instants       : ",
    fmti(sum(pan$retained & pan$full_instant)))

st$hours <- as.numeric(difftime(st$t1, st$t0, units = "hours"))
say_rule("records without the power signal")
nap <- is.na(pan$power)
say("records with missing power         : ", fmti(sum(nap)),
    "  (", fmt4(100 * mean(nap)), " per cent)")
say("  of which inside an excluded event : ", fmti(sum(nap & pan$excl_event)))
say("  of which inside a forced outage   : ", fmti(sum(nap & pan$ev_fo)),
    "  (", fmt4(sum(nap & pan$ev_fo) / 6), " WT hours)")
say("  of which with positive energy export, i.e. the WT was producing")
say("  while the SCADA power signal was missing : ",
    fmti(sum(nap & !is.na(pan$e_exp) & pan$e_exp > 0)))
say("  => these are gaps of the SCADA link and not outages; they are")
say("     removed from the estimator and from the benchmark alike")

say_rule("forced outage time, by union of the event intervals")
say("WT hours in the period              : ", fmti(nrow(pan) / 6))
say("WT hours in forced outage           : ", fmt4(sum(pan$ev_fo) / 6))
say("1 - FOR over the whole period       : ", fmt5(1 - mean(pan$ev_fo)))
say("(the plain sum of the event durations would give ",
    fmt4(sum(st$hours[has_iec & st$iec == "Forced outage" &
                        is.finite(st$hours) & st$hours > 0])),
    " hours, because overlapping events would be counted twice)")

say_rule("event log, hours by IEC 61400-26 category, plain sum of durations")
hh <- st[is.finite(st$hours) & st$hours > 0, ]
agg <- data.frame(category = names(tapply(hh$hours, hh$iec, sum)),
                  events = as.integer(tapply(hh$hours, hh$iec, length)),
                  hours = round(as.numeric(tapply(hh$hours, hh$iec, sum)), 1))
say_block(agg[order(-agg$hours), ])

saveRDS(list(panel = pan, t_origin = t_origin, n_slot = n_slot,
             n_turb = N_TURB, turb_id = TURB_ID, years = YEARS, site = SITE,
             excl_iec = EXCL_IEC, panel_version = PANEL_VERSION),
        PANEL_RDS)
saveRDS(st, STATUS_RDS)
say("")
say("written: ", PANEL_RDS, " and ", STATUS_RDS)
close_report()
