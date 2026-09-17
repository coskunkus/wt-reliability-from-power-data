## =================================================================
##  push_to_github.R
##
##  Replaces the whole of the GitHub repository with this folder.  The
##  repository is meant to hold the code of one paper and nothing else,
##  so a fresh history is built here and pushed over whatever is on
##  GitHub: every file that is on the remote and not in this folder
##  disappears.
##
##  What is pushed: code/ , README.md , LICENSE , .gitignore ,
##  .gitattributes .  The data, the reports, the tables and the figures
##  are not pushed; they are what a run produces.
##
##  Run it from the folder that holds this file:
##
##      setwd("<the folder that holds this file>")
##      source("push_to_github.R")
##
##  Authentication is left to git.  Git on Windows and on macOS keeps
##  the GitHub login in its credential helper and a push to an https
##  remote then authenticates by itself, asking in a browser the first
##  time.  A personal access token is needed only where there is no
##  helper, and is then given in GITHUB_PAT or GITHUB_TOKEN.
##
##  Base R only.  Nothing outside this folder is read or written.
## =================================================================

GH_USER <- "coskunkus"
GH_REPO <- "wt-reliability-from-power-data"
BRANCH  <- "main"

## ---- the folder that holds this file ----------------------------
local({
  here <- NA_character_
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) here <- dirname(normalizePath(of, mustWork = FALSE))
  }
  if (is.na(here)) {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg))
      here <- dirname(normalizePath(sub("^--file=", "", arg[1]),
                                    mustWork = FALSE))
  }
  if (is.na(here)) here <- getwd()
  if (!file.exists(file.path(here, "code", "core_wtrel.R")))
    stop("code/core_wtrel.R was not found beside this script; ",
         "setwd() to the folder that holds it first.", call. = FALSE)
  setwd(here)
})
PROJECT <- getwd()
cat("\nthe folder to publish : ", PROJECT, "\n", sep = "")

## ---- git and the credentials ------------------------------------
if (!nzchar(Sys.which("git"))) stop("git was not found on the PATH.")
pat <- Sys.getenv("GITHUB_PAT", Sys.getenv("GITHUB_TOKEN", ""))
helper <- suppressWarnings(system2("git", c("config", "--get",
                                            "credential.helper"),
                                   stdout = TRUE, stderr = FALSE))
cat("authentication        : ",
    if (nzchar(pat)) "the token of the environment" else
      if (length(helper)) paste0("the credential helper of git (",
                                 helper[1], ")") else
        "none found, git will ask", "\n", sep = "")

## ---- a staging folder whose path is ASCII only ------------------
## git receives its arguments in the system encoding, so that a folder
## name with a Turkish letter reaches it corrupted and every command
## fails.  No path is therefore passed to git: the files are copied
## under the temporary folder and the working directory is changed into
## it before each call.
is_ascii <- function(x) !grepl("[^\x01-\x7f]", x)
stage <- file.path(tempdir(), paste0("ghpush_", format(Sys.time(), "%H%M%S")))
if (!is_ascii(stage)) stage <- file.path("C:/tmp", basename(stage))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
if (!is_ascii(stage)) stop("No ASCII staging path is available.", call. = FALSE)

PUBLISH <- c("README.md", "LICENSE", ".gitignore", ".gitattributes",
             "push_to_github.R")
file.copy(file.path(PROJECT, "code"), stage, recursive = TRUE)
for (f in PUBLISH)
  if (file.exists(file.path(PROJECT, f)))
    file.copy(file.path(PROJECT, f), stage, overwrite = TRUE)

## nothing that a run produces is ever published
for (d in c("data", "runs", "results", "figures", "tables", "Claude outputs"))
  unlink(file.path(stage, d), recursive = TRUE)

files <- list.files(stage, recursive = TRUE, all.files = TRUE, no.. = TRUE)
cat("files to publish      : ", length(files), "\n", sep = "")
cat("\n", paste0("  ", sort(files), collapse = "\n"), "\n\n", sep = "")

## ---- the fresh history ------------------------------------------
old <- getwd(); on.exit(setwd(old), add = TRUE)
setwd(stage)
## An argument that carries a space has to reach git as one argument, and
## system2() hands them to the shell, so every such argument is quoted in
## the way the shell of the platform expects.
qq <- function(x) shQuote(x, type = if (.Platform$OS.type == "windows")
                                      "cmd" else "sh")
run <- function(...) {
  st <- system2("git", c(...), stdout = TRUE, stderr = TRUE)
  code <- attr(st, "status")
  if (length(st)) cat(paste0("  ", st, collapse = "\n"), "\n", sep = "")
  if (!is.null(code) && code != 0L)
    stop("git ", paste(c(...), collapse = " "), " failed.", call. = FALSE)
  invisible(st)
}
remote <- if (nzchar(pat))
  sprintf("https://%s@github.com/%s/%s.git", pat, GH_USER, GH_REPO) else
  sprintf("https://github.com/%s/%s.git", GH_USER, GH_REPO)

cat("building a fresh history\n")
run("init", "-q", "-b", BRANCH)
run("config", "user.name",  qq("Coskun Kus"))
run("config", "user.email", "coskun@selcuk.edu.tr")
run("add", "-A")
run("commit", "-q", "-m",
    qq(paste0("The code of the paper, rebuilt on ", format(Sys.Date()))))
run("remote", "add", "origin", qq(remote))

cat("\npushing over ", GH_USER, "/", GH_REPO, ", branch ", BRANCH,
    "\nevery file on the remote that is not in this folder will be removed\n",
    sep = "")
run("push", "--force", "origin", BRANCH)

cat("\ndone.  https://github.com/", GH_USER, "/", GH_REPO, "\n", sep = "")
