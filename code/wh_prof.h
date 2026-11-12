/* wh_prof.h (11-11-2026) -- the 2024 warp shell's radial construction (Fuchs et al., arXiv:2405.02709, section 3, and the Warp Factory
 * code metricGet_WarpShellComoving.m / TOVconstDensity.m / alphaNumericSolver.m / compactSigmoid.m) as finite integer sums in Q.40, and
 * the DECLARED metric functions between the grid points. Shared by wh_warp.c (and the shell-only program's successors).
 *
 * UNITS: R2 = 1 (20 m), G = c = 1; m_total = 1/6 (their M = R2 c^2/6G); R1 = 1/2; rho = 1/(7 pi). Energies in c^4/(G R2^2) = 3.0256e41 J/m^3.
 * THE TABLES (N = 100000 points, r_i = i dr, dr = rmax/(N-1), rmax = 1.2 sqrt(29.9^2 + 29.9^2 + 0.4^2) m / 20 m, their grid):
 *   rho_s  the density 1/(7 pi) on (R1, R2), MATLAB-smoothed (four passes, span 7159)
 *   f_     4 pi r^2 rho_s;  m_ its cumulative trapezoid (the mass function)
 *   P_s    the constant-density TOV pressure (closed form, total mass) for r < R2 -- variant AS-PUBLISHED (1) zeroes it inside R1 as the
 *          paper says AND as the public code does (TOVconstDensity.m multiplies by the step density, zero for r <= R1); variant
 *          HYPOTHETICAL (0) keeps it there (neither the paper nor its code) -- MATLAB-smoothed (four passes, span 3999)
 *   da_    a' = (m + 4 pi r^3 P_s) / (r^2 - 2 m r);  a_ its cumulative trapezoid pinned to (1/2) ln(1 - 2m/r) at the last point
 *   S0 -> smoothed twice (span 3999) -> bp_ = its central difference (the unit shift's derivative table), sh_ = the trapezoid of bp_
 *          from sh_0 = 1 (so the unit shift S is C^1 and within O(dr^2) of the smoothed sigmoid; the max deviation is printed)
 * THE DECLARED OBJECT: between grid points f, da and bp are linear, so m, a and S are piecewise quadratic and C^1; m', a', S' are
 * piecewise linear; m'', a'', S'' are piecewise constant. prof_mid(i) returns everything at the midpoint of interval i, where it all
 * exists. The shift used by a metric is beta(r) = beta0 S(r) (their vWarp = 0.02 at the start).
 */
#ifndef WH_PROF_H
#define WH_PROF_H
#include "whq.h"

#define N 100000
static q rho_s[N], f_[N], m_[N], P_s[N], da_[N], a_[N], sh_[N], bp_[N];
static int64_t S_[N + 1];
static q RHO, RMAX, DR, S23, SHDEV;                   /* 1/(7 pi), rmax, dr, sqrt(2/3), max |S - smoothed sigmoid| */
static const q R1 = ONE / 2, R2 = ONE;

static void smooth_passes(q *y, int span, int passes) /* MATLAB smooth(): centred mean, span rounded down to odd, shrinking windows at the ends */
{
    const int w = (span % 2 ? span : span - 1) / 2;
    for (int pass = 0; pass < passes; pass++) {
        S_[0] = 0; for (int i = 0; i < N; i++) S_[i + 1] = S_[i] + y[i];
        for (int i = 0; i < N; i++) {
            int wi = w; if (i < wi) wi = i; if (N - 1 - i < wi) wi = N - 1 - i;
            y[i] = rdiv(S_[i + wi + 1] - S_[i - wi], 2 * wi + 1);
        }
    }
}
static q unit_sigmoid(q r)                            /* compactSigmoid.m with sigma = 0, Rbuff = 0: 1 inside R1, 0 beyond R2, smooth between */
{
    if (r <= R1) return ONE;
    if (r >= R2) return 0;
    const q X = mulq(R2 - R1, divq(ONE, r - R2) + divq(ONE, r - R1));   /* f = 1/(1 + e^{-X}) */
    if (X >= 0) { const q e = expq(-X); return divq(ONE, ONE + e); }
    const q e = expq(X); return divq(e, ONE + e);
}
static void prof_build(int variant)                   /* variant 0 = HYPOTHETICAL (pressure kept inside R1), 1 = AS-PUBLISHED (zeroed there: paper and code) */
{
    RHO = divq(ONE, 7 * PI);
    RMAX = divq(6 * sqrtq(divq(178818 * ONE, 100 * ONE)), 100 * ONE);
    DR = divq(RMAX, (q)(N - 1) * ONE);
    S23 = sqrtq(divq(2 * ONE, 3 * ONE));
    for (int i = 0; i < N; i++) { const q r = (q)i * DR; rho_s[i] = (r > R1 && r < R2) ? RHO : 0; }
    smooth_passes(rho_s, 7159, 4);
    for (int i = 0; i < N; i++) { const q r = (q)i * DR; f_[i] = mulq(mulq(4 * PI, mulq(r, r)), rho_s[i]); }
    m_[0] = 0; for (int i = 1; i < N; i++) m_[i] = m_[i - 1] + rdiv(mulq(f_[i - 1] + f_[i], DR), 2);
    for (int i = 0; i < N; i++) {
        const q r = (q)i * DR; q P = 0;
        if (r < R2 && !(variant && r < R1)) {
            const q b = sqrtq(ONE - divq(mulq(r, r), 3 * ONE));
            P = divq(mulq(RHO, b - S23), 3 * S23 - b);
        }
        P_s[i] = P;
    }
    smooth_passes(P_s, 3999, 4);
    da_[0] = 0;
    for (int i = 1; i < N; i++) {
        const q r = (q)i * DR, r2 = mulq(r, r);
        da_[i] = divq(m_[i] + mulq(mulq(4 * PI, mulq(r2, r)), P_s[i]), r2 - mulq(2 * m_[i], r));
    }
    a_[0] = 0; for (int i = 1; i < N; i++) a_[i] = a_[i - 1] + rdiv(mulq(da_[i - 1] + da_[i], DR), 2);
    {
        const q rend = (q)(N - 1) * DR;
        const q pin = rdiv(lnq(ONE - divq(2 * m_[N - 1], rend)), 2) - a_[N - 1];
        for (int i = 0; i < N; i++) a_[i] += pin;
    }
    /* the unit shift: the compact sigmoid smoothed twice, then its derivative table and the C^1 integral of that table */
    for (int i = 0; i < N; i++) sh_[i] = unit_sigmoid((q)i * DR);
    smooth_passes(sh_, 3999, 2);
    bp_[0] = 0; bp_[N - 1] = 0;
    for (int i = 1; i + 1 < N; i++) bp_[i] = divq(sh_[i + 1] - sh_[i - 1], 2 * DR);
    {
        q S = ONE; SHDEV = 0;
        for (int i = 0; i < N; i++) {
            const q d = S - sh_[i]; if (d > SHDEV) SHDEV = d; if (-d > SHDEV) SHDEV = -d;
            sh_[i] = S;                                                        /* sh_ now holds the C^1 unit shift at the nodes */
            if (i + 1 < N) S += rdiv(mulq(bp_[i] + bp_[i + 1], DR), 2);
        }
    }
}
typedef struct { q r, m, mp, mpp, a, ap, app, S, Sp, Spp; } wh_mid;
static void prof_mid(int i, wh_mid *o)                /* everything at the midpoint of interval i (0 <= i <= N-2) */
{
    const q h2 = DR / 2, h8 = DR / 8;
    o->r = (q)i * DR + h2;
    o->m = m_[i] + mulq(f_[i], h2) + mulq(f_[i + 1] - f_[i], h8);
    o->mp = rdiv(f_[i] + f_[i + 1], 2);
    o->mpp = divq(f_[i + 1] - f_[i], DR);
    o->a = a_[i] + mulq(da_[i], h2) + mulq(da_[i + 1] - da_[i], h8);
    o->ap = rdiv(da_[i] + da_[i + 1], 2);
    o->app = divq(da_[i + 1] - da_[i], DR);
    o->S = sh_[i] + mulq(bp_[i], h2) + mulq(bp_[i + 1] - bp_[i], h8);
    o->Sp = rdiv(bp_[i] + bp_[i + 1], 2);
    o->Spp = divq(bp_[i + 1] - bp_[i], DR);
}
#endif
