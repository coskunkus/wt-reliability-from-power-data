#!/usr/bin/env python3
# =================================================================
#  check_formulas.py
#
#  Second, independent implementation of the closed-form expressions
#  of
#     S. Eryilmaz and C. Kus,
#     "On the estimation of wind turbine reliability from turbine
#      power output data".
#
#  The formulas are coded here directly from the equations of the
#  paper and not translated from code/core_wtrel.R, so that an error
#  in the transcription from the paper to the code cannot survive
#  both scripts.  The two layers and the pass criteria are those of
#  code/check_formulas.R: exact checks against an independent route
#  first, then Monte Carlo checks reported as z-scores with |z| < 4
#  passing.
#
#  Output : <run folder>/results/check_formulas_py_report.txt
#
#  Requires numpy and scipy.  Run time is about one minute.
#
#      python3 code/check_formulas.py
# =================================================================
import os
import sys
from datetime import datetime

import numpy as np
from scipy import integrate, optimize, stats

# ---- locate the project root ------------------------------------
def find_root():
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (here, os.path.dirname(here), os.getcwd(),
              os.path.dirname(os.getcwd())):
        if os.path.exists(os.path.join(d, "code", "core_wtrel.R")):
            return d
        if os.path.exists(os.path.join(d, "core_wtrel.R")):
            return os.path.dirname(d) if os.path.basename(d) == "code" else d
    raise SystemExit("core_wtrel.R was not found; run from the project root.")


ROOT = find_root()
os.chdir(ROOT)
# the run folder that the R scripts made, so that both reports of a run
# land together; a folder is made here when this script is run first
_latest = os.path.join("runs", "LATEST.txt")
if os.environ.get("WT_RUN"):
    RUN_ROOT = os.environ["WT_RUN"]
elif os.path.exists(_latest):
    RUN_ROOT = open(_latest).read().strip().splitlines()[0]
else:
    import datetime
    RUN_ROOT = os.path.join("runs",
                            datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S"))
    os.makedirs(RUN_ROOT, exist_ok=True)
    os.makedirs("runs", exist_ok=True)
    open(_latest, "w").write(RUN_ROOT + "\n")
os.makedirs(os.path.join(RUN_ROOT, "results"), exist_ok=True)
REPORT = open(os.path.join(RUN_ROOT, "results",
                           "check_formulas_py_report.txt"), "w")


def say(*parts):
    txt = "".join(str(x) for x in parts)
    print(txt)
    REPORT.write(txt + "\n")


def rule(title=None):
    say("-" * 70)
    if title is not None:
        say(title)
        say("-" * 70)


N_PASS = 0
N_FAIL = 0
FAILS = []


def check_rel(name, value, target, tol=1e-8):
    global N_PASS, N_FAIL
    rel = abs(value - target) / abs(target) if target != 0 else abs(value - target)
    ok = np.isfinite(rel) and rel < tol
    if ok:
        N_PASS += 1
    else:
        N_FAIL += 1
        FAILS.append(name)
    say("%-58s %-4s  rel.err = %.3e" % (name, "PASS" if ok else "FAIL", rel))
    return ok


def check_z(name, est, target, se, tol=4.0):
    global N_PASS, N_FAIL
    z = (est - target) / se
    ok = np.isfinite(z) and abs(z) < tol
    if ok:
        N_PASS += 1
    else:
        N_FAIL += 1
        FAILS.append(name)
    say("%-58s %-4s  z = %+7.3f" % (name, "PASS" if ok else "FAIL", z))
    return ok


# =================================================================
#  The model, coded from the equations of the paper
# =================================================================
def g_wt(v, v_ci, v_r, v_co, P_r):
    """Equation (2)."""
    v = np.asarray(v, dtype=float)
    out = np.zeros_like(v)
    mid = (v >= v_ci) & (v < v_r)
    out[mid] = P_r * (v[mid] ** 3 - v_ci ** 3) / (v_r ** 3 - v_ci ** 3)
    out[(v >= v_r) & (v < v_co)] = P_r
    return out


def F_weib(v, k, c):
    return 1.0 - np.exp(-((np.asarray(v, dtype=float) / c) ** k))


def f_weib(v, k, c):
    v = np.asarray(v, dtype=float)
    return (k / c) * (v / c) ** (k - 1.0) * np.exp(-((v / c) ** k))


def H1(x, k, c, v_ci, v_r, v_co, P_r):
    """The function H1 of Equation (3)."""
    v_x = ((x / P_r) * (v_r ** 3 - v_ci ** 3) + v_ci ** 3) ** (1.0 / 3.0)
    return 1.0 - F_weib(v_co, k, c) + F_weib(v_x, k, c)


def h1(x, k, c, v_ci, v_r, v_co, P_r):
    """Equation (6), the derivative of H1."""
    v_x = ((x / P_r) * (v_r ** 3 - v_ci ** 3) + v_ci ** 3) ** (1.0 / 3.0)
    return f_weib(v_x, k, c) * (v_r ** 3 - v_ci ** 3) / (3.0 * P_r * v_x ** 2)


def d_par(k, c, v_ci, v_co):
    return F_weib(v_co, k, c) - F_weib(v_ci, k, c)


def P_zero(p, k, c, v_ci, v_co):
    """From Equation (4)."""
    return p * (F_weib(v_ci, k, c) + 1.0 - F_weib(v_co, k, c)) + 1.0 - p


def P_rated(p, k, c, v_r, v_co):
    """From Equation (4)."""
    return p * (F_weib(v_co, k, c) - F_weib(v_r, k, c))


def moment_g(k, c, v_ci, v_r, v_co, P_r, j=1):
    """The bracket of Equations (5) and (10), in the power scale."""
    val = integrate.quad(lambda x: x ** j * h1(x, k, c, v_ci, v_r, v_co, P_r),
                         0.0, P_r, epsabs=1e-13, epsrel=1e-13)[0]
    return val + P_r ** j * (F_weib(v_co, k, c) - F_weib(v_r, k, c))


def mse_moment(p, n, k, c, v_ci, v_r, v_co, P_r):
    """Equation (9)."""
    m1 = moment_g(k, c, v_ci, v_r, v_co, P_r, 1)
    m2 = moment_g(k, c, v_ci, v_r, v_co, P_r, 2)
    return (p * m2 - p ** 2 * m1 ** 2) / (n * m1 ** 2)


def mse_mle(p, n, k, c, v_ci, v_co):
    """Equation (13)."""
    d = d_par(k, c, v_ci, v_co)
    return p * (1.0 - p * d) / (n * d)


def mse_mmle(p, n, k, c, v_ci, v_co):
    """Mean squared error of the modified MLE of Section 2.2."""
    d = d_par(k, c, v_ci, v_co)
    kk = np.arange(n + 1)
    pk = stats.binom.pmf(kk, n, 1.0 - p * d)
    thr = np.floor(n * (1.0 - d))
    low = kk > thr
    Pge1 = pk[~low].sum()
    ph = (1.0 - kk[low] / n) / d
    e1 = float((ph * pk[low]).sum() + Pge1)
    e2 = float((ph ** 2 * pk[low]).sum() + Pge1)
    return (e1 - p) ** 2 + (e2 - e1 ** 2)


def p_hat(m0, n, d):
    """Equation (12)."""
    return (1.0 - m0 / n) / d


# =================================================================
say("check_formulas.py  -  ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
say("python ", sys.version.split()[0], " numpy ", np.__version__,
    " scipy ", __import__("scipy").__version__)

C1 = dict(k=2.0, c=5.0, v_ci=3.0, v_r=11.0, v_co=20.0, P_r=5.5)
C2 = dict(k=2.4204, c=6.8501, v_ci=3.0, v_r=12.5, v_co=24.0, P_r=2050.0)

# =================================================================
#  Layer 1.  Exact checks
# =================================================================
for nm, C in (("C1", C1), ("C2", C2)):
    k, c = C["k"], C["c"]
    v_ci, v_r, v_co, P_r = C["v_ci"], C["v_r"], C["v_co"], C["P_r"]
    d = d_par(k, c, v_ci, v_co)
    rule("Layer 1, exact checks, configuration %s:  k = %g , c = %g , "
         "v_ci = %g , v_r = %g , v_co = %g , P_r = %g"
         % (nm, k, c, v_ci, v_r, v_co, P_r))

    # --- Equation (2) --------------------------------------------
    check_rel("%s  Eq. (2): g(v_ci) = 0" % nm,
              float(g_wt([v_ci], v_ci, v_r, v_co, P_r)[0]), 0.0, 1e-12)
    check_rel("%s  Eq. (2): g is continuous at v_r" % nm,
              float(g_wt([v_r - 1e-9], v_ci, v_r, v_co, P_r)[0]), P_r, 1e-8)
    check_rel("%s  Eq. (2): g(v_co) = 0" % nm,
              float(g_wt([v_co], v_ci, v_r, v_co, P_r)[0]), 0.0, 1e-12)

    # --- Equation (3), H1 against integration of f ---------------
    for x in (0.05 * P_r, 0.3 * P_r, 0.6 * P_r, 0.9 * P_r):
        v_x = ((x / P_r) * (v_r ** 3 - v_ci ** 3) + v_ci ** 3) ** (1.0 / 3.0)
        ind = (integrate.quad(lambda t: f_weib(t, k, c), 0.0, v_x,
                              epsabs=1e-14, epsrel=1e-14)[0]
               + integrate.quad(lambda t: f_weib(t, k, c), v_co, np.inf,
                                epsabs=1e-14, epsrel=1e-14)[0])
        check_rel("%s  Eq. (3): H1(%.4g) by integration" % (nm, x),
                  float(H1(x, k, c, v_ci, v_r, v_co, P_r)), ind, 1e-8)

    # --- the atoms and the total mass of Equation (4) ------------
    for p in (0.5, 0.9, 0.99):
        tgt = p * (integrate.quad(lambda t: f_weib(t, k, c), 0.0, v_ci,
                                  epsabs=1e-14, epsrel=1e-14)[0]
                   + integrate.quad(lambda t: f_weib(t, k, c), v_co, np.inf,
                                    epsabs=1e-14, epsrel=1e-14)[0]) + (1.0 - p)
        check_rel("%s  Eq. (4): P{P_WT = 0}, p = %.2f" % (nm, p),
                  float(P_zero(p, k, c, v_ci, v_co)), tgt, 1e-8)
        tgt2 = p * integrate.quad(lambda t: f_weib(t, k, c), v_r, v_co,
                                  epsabs=1e-14, epsrel=1e-14)[0]
        check_rel("%s  Eq. (4): P{P_WT = P_r}, p = %.2f" % (nm, p),
                  float(P_rated(p, k, c, v_r, v_co)), tgt2, 1e-8)
        mass = (P_zero(p, k, c, v_ci, v_co) + P_rated(p, k, c, v_r, v_co)
                + p * integrate.quad(lambda x: h1(x, k, c, v_ci, v_r, v_co, P_r),
                                     0.0, P_r, epsabs=1e-13, epsrel=1e-13)[0])
        check_rel("%s  Eq. (4): total mass one, p = %.2f" % (nm, p),
                  float(mass), 1.0, 1e-7)

    # --- Equation (6) by a central difference of H1 --------------
    x0 = 0.4 * P_r
    e = 1e-5 * P_r
    fd = (H1(x0 + e, k, c, v_ci, v_r, v_co, P_r)
          - H1(x0 - e, k, c, v_ci, v_r, v_co, P_r)) / (2.0 * e)
    check_rel("%s  Eq. (6): h1 by finite difference" % nm,
              float(h1(x0, k, c, v_ci, v_r, v_co, P_r)), float(fd), 1e-6)

    # --- Equations (5) and (10): power scale against wind scale --
    m1_v = (integrate.quad(lambda v: g_wt([v], v_ci, v_r, v_co, P_r)[0]
                           * f_weib(v, k, c), v_ci, v_r,
                           epsabs=1e-13, epsrel=1e-13)[0]
            + P_r * (F_weib(v_co, k, c) - F_weib(v_r, k, c)))
    m2_v = (integrate.quad(lambda v: g_wt([v], v_ci, v_r, v_co, P_r)[0] ** 2
                           * f_weib(v, k, c), v_ci, v_r,
                           epsabs=1e-13, epsrel=1e-13)[0]
            + P_r ** 2 * (F_weib(v_co, k, c) - F_weib(v_r, k, c)))
    check_rel("%s  Eq. (5): mean, power scale against wind scale" % nm,
              moment_g(k, c, v_ci, v_r, v_co, P_r, 1), m1_v, 1e-7)
    check_rel("%s  Eq. (10): second moment, both scales" % nm,
              moment_g(k, c, v_ci, v_r, v_co, P_r, 2), m2_v, 1e-7)

    # --- Equations (9), (12), (13) by exact enumeration ----------
    for n in (5, 8):
        for p in (0.6, 0.9, 0.98):
            kk = np.arange(n + 1)
            pk = stats.binom.pmf(kk, n, 1.0 - p * d)
            ph = (1.0 - kk / n) / d
            check_rel("%s  Eq. (12): E[p_hat] = p, n = %d, p = %.2f"
                      % (nm, n, p), float((ph * pk).sum()), p, 1e-10)
            check_rel("%s  Eq. (13): MSE by enumeration, n = %d, p = %.2f"
                      % (nm, n, p), float(((ph - p) ** 2 * pk).sum()),
                      mse_mle(p, n, k, c, v_ci, v_co), 1e-10)
            pm = np.minimum(ph, 1.0)
            check_rel("%s  MSE of p_M by enumeration, n = %d, p = %.2f"
                      % (nm, n, p), float(((pm - p) ** 2 * pk).sum()),
                      mse_mmle(p, n, k, c, v_ci, v_co), 1e-10)
            m1 = moment_g(k, c, v_ci, v_r, v_co, P_r, 1)
            m2 = moment_g(k, c, v_ci, v_r, v_co, P_r, 2)
            var_p = p * m2 - (p * m1) ** 2
            check_rel("%s  Eq. (9): MSE from Var(P_WT), n = %d, p = %.2f"
                      % (nm, n, p),
                      mse_moment(p, n, k, c, v_ci, v_r, v_co, P_r),
                      var_p / (n * m1 ** 2), 1e-10)

    # --- the score equation, the information and the bound -------
    n, m0 = 40, 4
    def loglik(p):
        return m0 * np.log(1.0 - p * d) + (n - m0) * np.log(p * d)

    def score(p):
        return -m0 * d / (1.0 - p * d) + (n - m0) / p

    ph_cf = p_hat(m0, n, d)
    check_rel("%s  Eq. (12): root of the score equation" % nm,
              optimize.brentq(score, 1e-6, 1.0 / d - 1e-6, xtol=1e-15),
              ph_cf, 1e-8)
    h = 1e-6
    check_rel("%s  Eq. (11): numerical score at the MLE is zero" % nm,
              (loglik(ph_cf + h) - loglik(ph_cf - h)) / (2.0 * h * n),
              0.0, 1e-5)
    d2 = (loglik(ph_cf + 1e-5) - 2.0 * loglik(ph_cf)
          + loglik(ph_cf - 1e-5)) / 1e-10
    check_rel("%s  Eq. (11): observed information is positive" % nm,
              1.0 if -d2 > 0 else 0.0, 1.0, 1e-12)
    p0 = 0.95
    fisher = n * d ** 2 / ((1.0 - p0 * d) * (p0 * d))
    check_rel("%s  Eq. (13) attains the Cramer-Rao bound" % nm,
              1.0 / fisher, mse_mle(p0, n, k, c, v_ci, v_co), 1e-12)

# --- Table 1 of the paper ----------------------------------------
rule("Layer 1, Table 1 of the paper")
paper = [0.0257, 0.0459, 0.0251, 0.0454, 0.0245, 0.0449, 0.0240, 0.0444,
         0.0170, 0.0306, 0.0166, 0.0303, 0.0161, 0.0300, 0.0158, 0.0296,
         0.0131, 0.0230, 0.0126, 0.0227, 0.0122, 0.0225, 0.0119, 0.0222,
         0.0090, 0.0153, 0.0085, 0.0151, 0.0082, 0.0150, 0.0079, 0.0148]
ours = []
for n in (10, 15, 20, 30):
    for p in (0.95, 0.96, 0.97, 0.98):
        ours.append(round(mse_mmle(p, n, 2.0, 5.0, 3.0, 20.0), 4))
        ours.append(round(mse_mle(p, n, 2.0, 5.0, 3.0, 20.0), 4))
say("  p      n    MSE(p_M)   MSE(p_hat)")
i = 0
for n in (10, 15, 20, 30):
    for p in (0.95, 0.96, 0.97, 0.98):
        say("  %.2f  %3d    %.4f      %.4f" % (p, n, ours[i], ours[i + 1]))
        i += 2
check_rel("Table 1: all 32 entries agree to four decimals",
          float(np.max(np.abs(np.array(ours) - np.array(paper)))), 0.0, 1e-12)

# =================================================================
#  Layer 2.  Monte Carlo checks
# =================================================================
rule("Layer 2, Monte Carlo checks, |z| < 4 passes")
rng = np.random.default_rng(20260909)
B = 200000
k, c = C2["k"], C2["c"]
v_ci, v_r, v_co, P_r = C2["v_ci"], C2["v_r"], C2["v_co"], C2["P_r"]
d = d_par(k, c, v_ci, v_co)


def se_var(x):
    n = len(x)
    m4 = np.mean((x - x.mean()) ** 4)
    v = x.var(ddof=1)
    return np.sqrt(max(m4 - v ** 2, 0.0) / n)


for p in (0.90, 0.97):
    for n in (100, 400):
        m0 = rng.binomial(n, 1.0 - p * d, size=B)
        ph = (1.0 - m0 / n) / d
        pm = np.minimum(ph, 1.0)
        check_z("MC  E[p_hat] = p, n = %d, p = %.2f" % (n, p),
                ph.mean(), p, ph.std(ddof=1) / np.sqrt(B))
        check_z("MC  Var(p_hat) = Eq. (13), n = %d, p = %.2f" % (n, p),
                ph.var(ddof=1), mse_mle(p, n, k, c, v_ci, v_co), se_var(ph))
        sq = (pm - p) ** 2
        check_z("MC  MSE(p_M) formula, n = %d, p = %.2f" % (n, p),
                sq.mean(), mse_mmle(p, n, k, c, v_ci, v_co),
                sq.std(ddof=1) / np.sqrt(B))

for p in (0.90, 0.97):
    n, Bm = 400, 40000
    V = rng.weibull(k, size=(Bm, n)) * c
    X = rng.binomial(1, p, size=(Bm, n))
    PW = g_wt(V.ravel(), v_ci, v_r, v_co, P_r).reshape(Bm, n) * X
    z0 = float((PW <= 0).mean())
    check_z("MC  P{P_WT = 0} = 1 - pd, p = %.2f" % p,
            z0, P_zero(p, k, c, v_ci, v_co),
            np.sqrt(z0 * (1.0 - z0) / (Bm * n)))
    zr = float((PW >= P_r - 1e-9).mean())
    check_z("MC  P{P_WT = P_r}, p = %.2f" % p,
            zr, P_rated(p, k, c, v_r, v_co),
            np.sqrt(zr * (1.0 - zr) / (Bm * n)))
    pt = PW.mean(axis=1) / moment_g(k, c, v_ci, v_r, v_co, P_r, 1)
    check_z("MC  E[p_tilde] = p, n = %d, p = %.2f" % (n, p),
            pt.mean(), p, pt.std(ddof=1) / np.sqrt(Bm))
    check_z("MC  Var(p_tilde) = Eq. (9), n = %d, p = %.2f" % (n, p),
            pt.var(ddof=1), mse_moment(p, n, k, c, v_ci, v_r, v_co, P_r),
            se_var(pt))
    del V, X, PW

say("")
say("rate check: mean of sqrt(n) |p_hat - p_M| must fall as n doubles")
prev, ok_rate = None, True
p = 0.97
for n in (200, 400, 800, 1600):
    m0 = rng.binomial(n, 1.0 - p * d, size=B)
    ph = (1.0 - m0 / n) / d
    val = float(np.mean(np.sqrt(n) * np.abs(ph - np.minimum(ph, 1.0))))
    say("   n = %5d   sqrt(n) E|p_hat - p_M| = %.6e" % (n, val))
    if prev is not None and val > prev:
        ok_rate = False
    prev = val
if ok_rate:
    N_PASS += 1
else:
    N_FAIL += 1
    FAILS.append("rate of p_hat - p_M")
say("%-58s %s" % ("rate: the sequence is decreasing",
                  "PASS" if ok_rate else "FAIL"))

# =================================================================
rule("Summary")
say("checks passed : ", N_PASS)
say("checks failed : ", N_FAIL)
if N_FAIL:
    say("failed: ", " ; ".join(FAILS))
say("ALL CHECKS PASS" if N_FAIL == 0 else "SOME CHECKS FAIL")
REPORT.close()
sys.exit(0 if N_FAIL == 0 else 1)
