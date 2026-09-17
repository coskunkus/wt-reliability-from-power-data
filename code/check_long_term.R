## =================================================================
##  check_long_term.R
##
##  The estimator (12) uses the wind speed distribution through the
##  single number d = F(v_co) - F(v_ci).  In the section that number is
##  read off the wind speed record of the sample itself.  It is a
##  feature of the site rather than of the period, so that it may be
##  taken from a record that is independent of the sample, and this
##  script does so: every archive of the site that is present in data/
##  and does not belong to the period of the sample is read, the wind
##  speed alone, and d is obtained from those years.  The estimate of p
##  is then recomputed with that value, the count m_0 and the size n of
##  the sample being unchanged.
##
##  It reports in addition how far d moves from one year to another,
##  which is the quantity that decides whether a value taken from
##  another period may be used at all.
##
##  Input  : data/<Site>_SCADA_<year>_*.zip  for years outside those of
##           WT_SITES, together with data/<site>_panel.rds
##  Output : results/<site>check_long_term_report.txt
##
##  Download the further years from the deposit named in the report.
##  Base R only.  Run time is about one minute per year.
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
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

new_report(file.path(res_dir, paste0(tag("check_long_term"), "_report.txt")))
say("check_long_term.R  -  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)
say("site: ", SITE, " , deposit https://doi.org/", WT$doi)

v_ci <- WT$v_ci; v_co <- WT$v_co
IN_SAMPLE <- as.integer(WT$years)

## ---- which years are available ----------------------------------
## the archives are looked for wherever they are kept, as the core says
zips  <- basename(find_data_file(sprintf("^%s_SCADA_[0-9]{4}_.*\\.zip$", SITE)))
years <- sort(unique(as.integer(sub(sprintf("^%s_SCADA_([0-9]{4})_.*$", SITE),
                                    "\\1", zips))))
extra <- setdiff(years, IN_SAMPLE)
say_rule("the archives that were found")
say("folders searched                     : ",
    paste(data_dirs(), collapse = " , "))
say("years of the sample                  : ", paste(IN_SAMPLE, collapse = ", "))
say("years present                        : ", paste(years, collapse = ", "))
say("years outside the sample             : ",
    if (length(extra)) paste(extra, collapse = ", ") else "none")
if (!length(extra))
  stop("\n\n  No archive of a year outside ", paste(IN_SAMPLE, collapse = ", "),
       " was found.\n  The folders searched are ", paste(data_dirs(),
       collapse = " , "), ".\n  Download one or more further years of ", SITE,
       " from\n  https://doi.org/", WT$doi, "  into data/ and run this script",
       " again.\n", call. = FALSE)

## ---- the wind speed of those years ------------------------------
## Only the time stamp and the wind speed are read, so that a year is
## read in a few seconds.  No filtering of the event log is applied:
## the wind blows whether the WT is available or not, and it is the
## wind speed alone that d describes.
read_ws <- function(f) {
  hdr <- scan(f, what = "", sep = ",", skip = 9, nlines = 1, quiet = TRUE)
  want <- c("# Date and time", "Wind speed (m/s)")
  idx  <- match(want, hdr)
  if (anyNA(idx)) return(NULL)
  cc <- rep("NULL", length(hdr))
  cc[idx] <- c("character", "numeric")
  d <- read.csv(f, skip = 9, header = TRUE, check.names = FALSE,
                colClasses = cc, na.strings = c("NaN", "NA", ""))
  names(d) <- c("ts", "ws")
  d$ts <- as.POSIXct(d$ts, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  d
}

say_rule("the wind speed of the years outside the sample")
per_year <- data.frame()
ws_all   <- numeric(0)
for (y in extra) {
  yf <- year_folder(y)
  fs <- sort(list.files(yf$dir, paste0("^Turbine_Data_", SITE, "_.*\\.csv$"),
                        full.names = TRUE, recursive = TRUE))
  if (!length(fs)) stop("no SCADA file of ", y, " was found after unpacking",
                        call. = FALSE)
  w <- unlist(lapply(fs, function(f) { z <- read_ws(f); if (is.null(z))
    NULL else z$ws }), use.names = FALSE)
  if (yf$scratch) unlink(yf$dir, recursive = TRUE)
  w <- w[is.finite(w)]
  per_year <- rbind(per_year, data.frame(
    year = y, files = length(fs), records = length(w),
    mean = mean(w), F_vci = mean(w < v_ci),
    d = mean(w >= v_ci & w < v_co)))
  ws_all <- c(ws_all, w)
}
say_block(data.frame(year = per_year$year, files = per_year$files,
                     records = fmti(per_year$records),
                     mean = fmt4(per_year$mean),
                     F_vci = fmt5(per_year$F_vci), d = fmt5(per_year$d)))

d_long <- mean(ws_all >= v_ci & ws_all < v_co)
F_long <- mean(ws_all < v_ci)
say("")
say("the years outside the sample, taken together")
say("  records                            : ", fmti(length(ws_all)))
say("  F(v_ci) = P{V < ", v_ci, "}                 : ", fmt5(F_long))
say("  d       = P{", v_ci, " <= V < ", v_co, "}        : ", fmt5(d_long))
say("  largest difference between two years of d : ",
    fmt5(diff(range(per_year$d))))

## the same five families, fitted to this independent record
FAML <- fit_wind_families(ws_all, v_ci, v_co)
say("")
say("the five families fitted to this record, by maximum likelihood")
say_block(data.frame(family = FAML$tab$family,
                     AIC = formatC(FAML$tab$AIC, format = "f", digits = 1),
                     dAIC = formatC(FAML$tab$dAIC, format = "f", digits = 1),
                     d = fmt5(FAML$tab$d),
                     err_d = fmt5(FAML$tab$d - d_long)))
say("smallest AIC : ", FAML$best)

## ---- what it does to the estimate -------------------------------
need_file(PANEL_RDS, "code/prepare_data.R")
D   <- readRDS(PANEL_RDS)
need_panel(D)
pan <- D$panel
## the sample and the zero-power rule of the core, as every script uses
sub    <- estimation_sample(pan)
N      <- nrow(sub)
inband <- in_band(sub$ws, v_ci, v_co)
m0     <- sum(zero_of(sub, v_ci, v_co))
d_rec  <- mean(inband)

say_rule("the estimate of p, with d of the other years")
say("n                                    : ", fmti(N))
say("m0 under the assumption of (2)       : ", fmti(m0))
say("d of the sample                      : ", fmt5(d_rec),
    "  ->  (12) = ", fmt5(p_hat_wt(m0, N, d_rec)))
say("d of the years outside the sample    : ", fmt5(d_long),
    "  ->  (12) = ", fmt5(p_hat_wt(m0, N, d_long)))
say("difference between the two estimates : ",
    formatC(p_hat_wt(m0, N, d_long) - p_hat_wt(m0, N, d_rec),
            format = "f", digits = 5))
say("the estimate from the independent d is admissible : ",
    p_hat_wt(m0, N, d_long) <= 1)

saveRDS(list(per_year = per_year, d_long = d_long, F_long = F_long,
             fam = FAML$tab, d_rec = d_rec, N = N, m0 = m0),
        file.path("data", paste0(tag("long_term"), ".rds")))
close_report()
