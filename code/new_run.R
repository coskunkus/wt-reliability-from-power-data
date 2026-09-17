## =================================================================
##  new_run.R
##
##  Starts a new run.  It makes the folder
##      runs/<date>_<time>/
##  with results/ , tables/ and figures/ inside it, and writes its name
##  into runs/LATEST.txt.  Every script that is run afterwards reads
##  that file and writes there, so that the whole of a run lands in one
##  folder and two runs never mix.
##
##  Run this first, then the scripts in the order given in README.md.
##  Running it again starts another folder and leaves the earlier ones
##  untouched.
##
##  The caches of the data are not touched: they stay in data/ and are
##  shared by every run, the reading of the three archives being the
##  slow step.
##
##  Base R only.  It takes no time at all.
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
  for (d in cand) { r <- find_root(d); if (!is.na(r)) { setwd(r); return(invisible(NULL)) } }
  stop("core_wtrel.R was not found; setwd() to the project root first.")
})

stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
root  <- file.path("runs", stamp)
for (d in c(root, file.path(root, c("results", "tables", "figures"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
writeLines(root, file.path("runs", "LATEST.txt"))

writeLines(c(
  paste0("run          : ", stamp),
  paste0("started      : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("R            : ", R.version.string),
  paste0("platform     : ", R.version$platform),
  paste0("folder       : ", normalizePath(root, winslash = "/"))),
  file.path(root, "run_info.txt"))

cat("\nthe run folder is  ", root, "\n",
    "  results/   the reports of every script\n",
    "  tables/    the LaTeX of every table\n",
    "  figures/   every figure, as eps and as pdf\n\n",
    "Every script run from now on writes there.  Run them in the order\n",
    "given in README.md, and start the next run with this script again.\n\n",
    sep = "")
