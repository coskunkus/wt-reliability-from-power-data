# Wind turbine reliability from power output data

The code of

> S. Eryilmaz and C. Kus, *On the estimation of wind turbine reliability from
> turbine power output data*.

Every number, table and figure of the paper is produced by these scripts and by
nothing else. A clone reproduces the whole of the numerical work: the exact
expressions of Section 2, the tables and figures of Section 3, and the real data
application of Section 4 on the open Kelmarsh wind farm records.

The scripts are base R throughout, with one Python file that repeats the formula
check independently. No package has to be installed.

## Starting from a clone

    git clone https://github.com/coskunkus/wt-reliability-from-power-data.git
    cd wt-reliability-from-power-data
    mkdir data

The folder `data/` is not in the repository, the archives of the site being far
too large to publish and being published already by their own deposit, so it
has to be made. Open

    https://doi.org/10.5281/zenodo.5841833

in a browser and download these three files into the `data/` folder that was
just made, with their names unchanged:

    Kelmarsh_SCADA_2017_3083.zip
    Kelmarsh_SCADA_2018_3084.zip
    Kelmarsh_SCADA_2019_3085.zip

The folder then looks like this, and nothing else has to be arranged:

    wt-reliability-from-power-data/
        code/                       the scripts, from the repository
        data/                       the three archives, from the deposit
        README.md
        runs/                       made by the first run

The archives are never unpacked by hand, never renamed and never modified. Then
run the scripts in the order given below.

## Where the output of a run goes

Every report, table and figure is written under

    runs/<date>_<time>/results/     the report of every script, as plain text
    runs/<date>_<time>/tables/      the LaTeX of every table
    runs/<date>_<time>/figures/     every figure, as eps and as pdf

`code/new_run.R` makes such a folder and writes its name into
`runs/LATEST.txt`; every script run afterwards reads that file and writes
there, so that the whole of a run lands together and two runs never mix.
Running `new_run.R` again starts another folder and leaves the earlier ones
untouched. A script run without any run folder makes one for itself.

The caches of the data are not part of a run and stay in `data/`, which every
run shares, the reading of the three archives being the slow step.

## What runs without any data

Four scripts need no data at all and may be run on a fresh clone:

    Rscript code/new_run.R            # the run folder, first of all
    Rscript code/run_paper.R          # Tables 1 and 2, Figures 1, 2 and 3
    Rscript code/check_formulas.R     # 120 checks of the formulas
    Rscript code/check_section2.R     # 89 checks of the theory of Section 2
    Rscript code/check_common_mle.R   # 11 checks of the clustering and of
                                      # the likelihood under a common wind
    python3 code/check_formulas.py    # the formulas again, coded independently

The formula check is written twice, once in R and once in Python, from the
paper rather than from the other implementation, so that an error in the
transcription cannot survive both.

## The data, and where to get them

The data are the open Kelmarsh wind farm records of C. Plumley (2022): six
Senvion MM92 wind turbines at one site in Northamptonshire, United Kingdom,
released by Cubico Sustainable Investments Ltd under a CC-BY-4.0 licence. They
are held in the Zenodo deposit

    https://doi.org/10.5281/zenodo.5841833

The record is a ten-minute one and carries, for each wind turbine and each
ten-minute interval, the power output, the wind speed of the nacelle
anemometer and 297 further signals, together with an event log in which every
change of the state of a wind turbine is registered with its category in the
sense of IEC 61400-26. The paper uses the three complete years 2017, 2018 and
2019 and nothing else.

Open the deposit in a browser and download these three archives into `data/`,
beside `code/`, with their names unchanged:

    data/Kelmarsh_SCADA_2017_3083.zip      167 MB
    data/Kelmarsh_SCADA_2018_3084.zip      251 MB
    data/Kelmarsh_SCADA_2019_3085.zip      297 MB

The number at the end of a name is the file number of the deposit and it is
part of the name. A script that needs an archive and does not find it stops
with a message that names the archive, the deposit and the folder to put it in,
so that nothing has to be guessed.

The archives are large and are often kept once for several working folders, so
they are looked for in three places, in this order: the folder `data/` of the
project, a folder `data/` beside the project, and the folder named in the
environment variable `WT_DATA`. Keeping the project inside a folder whose `data/`
already holds the archives is therefore enough, and nothing has to be copied:

    <somewhere>/data/Kelmarsh_SCADA_2017_3083.zip      the archives, once
    <somewhere>/<this project>/code/                   the scripts
    <somewhere>/<this project>/data/                   what the run caches

Whatever the archives are read from, everything a run caches is written to
`data/` of the project, which is made when it is missing.

The archives are not unpacked by hand and are not modified: `prepare_data.R`
reads them as they are, keeps the eight signals it needs and writes the
ten-minute panel into `data/`. The archives themselves are never published with
the code, and `.gitignore` keeps `data/` out of the repository.

The characteristics of the wind turbine are the manufacturer's and are not in
the data set: the cut-in and the cut-out wind speeds are those of the technical
brochure of the model, and the nominal power and the rated wind speed those of
the uprated version installed at the site. They are held in `code/core_wtrel.R`,
in the list `WT_SITES`, with the source of each value named beside it, and no
script takes them from anywhere else.

## The order of the run

    setwd("<the folder that holds this README>")

    ## the folder of this run, always first
    source("code/new_run.R")          # makes runs/<date>_<time>/

    ## these need no data
    source("code/run_paper.R")        # Tables 1 and 2, Figures 1, 2 and 3
    source("code/check_formulas.R")   # the formulas, about two minutes
    source("code/check_section2.R")   # the theory of Section 2, about two minutes
    source("code/check_common_mle.R") # the additions, about two minutes

    ## these need the three archives
    source("code/prepare_data.R")     # Steps 1 and 2, about three minutes
    source("code/check_data.R")       # the data layer, under a minute
    source("code/run_all.R")          # the estimates of the record
    source("code/run_simulation.R")   # Step 3 replicated, about a minute
    source("code/run_example.R")      # Section 4 and its four tables

    ## optional, and needs a further year of the site in data/
    source("code/check_long_term.R")

`prepare_data.R` writes the ten-minute panel that the four scripts after it
read, and `run_simulation.R` reads what `run_all.R` writes, so the order within
the second block matters. The panel carries a version number, and a script that
meets an out-of-date panel says so and stops. `check_section2.R` reads the panel
as well, for the part of its checks that is carried out on the record itself, so
it is better run after `prepare_data.R`; without the panel it says so and skips
that part alone.

Nothing is written by hand and nothing is copied between scripts by hand: every
number of the paper is in one of the reports of a run, and every table of the
paper is one of the files under `tables/`.

## Which script produces what

    code/new_run.R           makes the folder of a run and names it in
                             runs/LATEST.txt

    code/core_wtrel.R        the shared core: the power curve g of (2), the
                             distribution of the power output, the two
                             estimators and their mean squared errors, the
                             clustering constant C and the two correlations,
                             the distribution of m_0 over a sample of periods,
                             the likelihood under a wind common to a period
                             and the Fisher information it carries, the fit of
                             the Weibull F, the reading of the SCADA files,
                             the sampling design of Step 3, the figure device
                             and the LaTeX table emitter. Every other script
                             sources this one and defines no formula of its
                             own

    code/run_paper.R         Section 3. Table 1, the mean squared errors of
                             the two maximum likelihood estimators for a
                             sample of independent readings; Table 2, the four
                             estimators when the same M = 30 readings are
                             collected as N = M/n periods of n wind turbines;
                             Figures 1, 2 and 3. Everything is exact and
                             nothing is simulated. The 32 entries of Table 1
                             are compared with the values printed in the
                             paper, and the script stops if they ever differ.
                             The row n = 1 of Table 2 must return the n = 30
                             block of Table 1, and that too is checked

    code/prepare_data.R      Steps 1 and 2 of the protocol: unpacks the three
                             archives, reads the eight signals needed out of
                             the 299 columns of each SCADA file, applies the
                             filtering of the event log one cell at a time,
                             and writes the ten-minute panel

    code/run_all.R           the estimates over the whole record and over each
                             of the eighteen wind turbine years, the four
                             availability indicators of the same cells, and
                             the reading that uses no event log

    code/run_simulation.R    Step 3 replicated a thousand times on the panel,
                             once with all six wind turbines read at each
                             instant and once with one, the two designs being
                             set against the variance of the paper

    code/run_example.R       Section 4. The sample of N periods and the n_j
                             wind turbines that remain in each, M and C, the
                             fit of F and the two values of d, the estimates
                             and their standard errors, the correlations of
                             the record against the closed form, the estimator
                             obtained under a common wind, and the four tables
                             of the section

    code/check_long_term.R   optional. Reads the wind speed of every archive
                             of the site in data/ whose year is not one of the
                             three of the sample, obtains d from those years
                             and recomputes the estimate with it, so that d
                             comes from a record independent of the sample

## The verification

Four scripts verify rather than produce. They are what makes the numbers of the
paper checkable without trusting any single implementation.

    code/check_formulas.R    every expression of the paper against an
                             independent route, exactly where an exact route
                             exists and by Monte Carlo as a z-score otherwise.
                             120 checks

    code/check_formulas.py   the same checks, coded in Python from the paper
                             and not from the R code

    code/check_section2.R    the theory of Section 2 in twelve groups: the
                             closed forms against the exact distribution of
                             m_0, the numerical maximisation of the
                             log-likelihood, a Monte Carlo of the sampling
                             mechanism itself, and a panel of periods with an
                             autocorrelated wind. 89 checks

    code/check_common_mle.R  the additions: the two mean squared errors at
                             C = 0 return the published formulas, the closed
                             form of the variance of m_0 agrees both with the
                             distribution of m_0 and with the mechanism, the
                             estimator obtained under a common wind is the
                             modified maximum likelihood estimator when one
                             wind turbine is read per period, the Fisher
                             information in closed form agrees with a
                             numerical one, the exact mean squared error
                             agrees with a simulation, and no sample sends
                             that estimator outside the unit interval.
                             11 checks

    code/check_data.R        every assumption made about the Kelmarsh files,
                             against the files: the unit and the meaning of
                             the power signal, the time stamp convention, the
                             zero-power state, the categories of the event
                             log, the missing records

A verification script prints one line per check and stops with a non-zero
status if any of them fails, so that it may be run unattended.

## Reproducibility

The random number generator is the L'Ecuyer-CMRG stream and the master seed is
20260909, so that every simulation returns the same numbers on every machine.
The reports carry the version of R that produced them.

## Publishing this folder

    setwd("<the folder that holds this README>")
    source("push_to_github.R")

The script builds a fresh history from this folder and pushes it over the
repository, so that every file on the remote that is not here is removed. Only
`code/`, `README.md`, `LICENSE`, `.gitignore`, `.gitattributes` and the script
itself are published; the data, the runs and everything a run produces are not.
Authentication is left to git, which keeps the login of GitHub in its credential
helper; a personal access token is needed only where there is no helper, and is
then given in `GITHUB_PAT`.

## Data availability

The data are the open Kelmarsh records named above. The code is released under
the MIT licence, which the file LICENSE carries.
