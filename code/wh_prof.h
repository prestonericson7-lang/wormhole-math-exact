/* wh_prof.h (11-11-2026; v2 9:00 PM after the math buddy's review; v3 11-12 1:00 AM) -- the 2024 warp shell's radial construction (Fuchs et al.,
 * VERSION 3 (11-12, the reviewer's whole file, verified by three builds): two changes against v2 and nothing else. (1) The second derivative's
 * dr^2 is kept exact in Q.80 as 2 DR^2 (v2 truncated dr^2 to 707 units of 2^-40 where its value is 707.819049, which divided every S''
 * by a number 0.116 percent too small and lowered every crossing of the shift limit by 3.25e-5: the correction published 12:03 AM).
 * (2) The mass table's trapezoid is kept exact in Q.80 (M2_) and m is rounded once per node, and da is formed from the exact numerator
 * (M2_ + 2^41 X) / (2 Dn) rounded once (v2 rounded twice, enough for a false DEC violation at beta0 = 0 on the 61-cell theta grid).
 * arXiv:2405.02709, section 3, and the Warp Factory code metricGet_WarpShellComoving.m / TOVconstDensity.m / alphaNumericSolver.m /
 * compactSigmoid.m) as finite integer sums in Q.40, and the DECLARED metric functions between the grid points. Shared by wh_warp.c.
 *
 * UNITS: R2 = 1 (20 m), G = c = 1; m_total = 1/6 (their M = R2 c^2/6G); R1 = 1/2; rho = 1/(7 pi). Energies in c^4/(G R2^2) = 3.0256e41 J/m^3.
 * THE TABLES (N = 100000 points, r_i = i dr, dr = rmax/(N-1), rmax = 1.2 sqrt(29.9^2 + 29.9^2 + 0.4^2) m / 20 m, their grid):
 *   rho_s  the density 1/(7 pi) on (R1, R2), MATLAB-smoothed (four passes, span 7159)
 *   f_     4 pi r^2 rho_s;  m_ its cumulative trapezoid (the mass function)
 *   MT     the trapezoid of the UNSMOOTHED step density (their TOVconstDensity uses M(end) of that cumtrapz: 0.16666023, not 1/6; v2 does the same)
 *   P_s    the constant-density TOV pressure with MT (their closed form) for r < R2 -- variant AS-PUBLISHED (1) zeroes it inside R1 as the
 *          paper says AND as the public code does (TOVconstDensity.m multiplies by the step density, zero for r <= R1); variant
 *          HYPOTHETICAL (0) keeps it there (neither the paper nor its code) -- MATLAB-smoothed (four passes, span 3999)
 *   da_    a' = (m + 4 pi r^3 P_s) / (r^2 - 2 m r);  a_ its cumulative trapezoid pinned to (1/2) ln(1 - 2m/r) at the last point
 *   THE SHIFT (v2): y = compactSigmoid at the nodes (sigma = 0, Rbuff = 0; the paper's eq. 27 has Rb > 0 and the example's value is not
 *          extracted: FIND), smoothed twice (span 3999) EXACTLY: the two box passes are one triangle kernel of integer weights, so
 *          (2w+1)^2 S_i is an exact 128-bit integer sum; S_i, S'_i (the central difference of those sums) and S''_i (their second
 *          difference) are each divided ONCE (S'' in two stages, the first exact to 1 part in 1e28). v1 re-integrated a rounded central
 *          difference and carried 4e-4 of rounding noise in S''; v2's S'' noise is one unit of 2^-40.
 * THE DECLARED OBJECT: between grid points f and da are linear, so m and a are piecewise quadratic and C^1; m', a' piecewise linear; m'', a''
 * piecewise constant. The shift: S'' on interval i is the exact kernel second difference; S' = bp_i + S''_i (r - r_i); S = S_i + bp_i (r - r_i)
 * + S''_i (r - r_i)^2 / 2; S' jumps by at most one unit of 2^-40 at a node (the roundings of bp). prof_mid(i) returns everything at the
 * midpoint of interval i. The shift used by a metric is beta(r) = beta0 S(r) (their vWarp = 0.02 at the start).
 */
#ifndef WH_PROF_H
#define WH_PROF_H
#include "whq.h"

#define N 100000
static q rho_s[N], f_[N], m_[N], P_s[N], da_[N], a_[N], sh_[N], bp_[N], spp_[N];
static int64_t S_[N + 1];
static q RHO, RMAX, DR, MT, SQA, SHDEV;              /* 1/(7 pi), rmax, dr, the step mass, sqrt(1 - 2 MT), max |S - the twice-smoothed sigmoid| */
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
    const q X = mulq(R2 - R1, divq(ONE, r - R2) + divq(ONE, r - R1));   /* f = 1/(1 + e^{-X}); r - R2 < 0: whq.h v2's divq takes it */
    if (X >= 0) { const q e = expq(-X); return divq(ONE, ONE + e); }
    const q e = expq(X); return divq(e, ONE + e);
}
/* ---- 128-bit helpers for the exact shift ---- */
typedef struct { uint64_t hi, lo; } u128;
static inline u128 u128_add(u128 a, u128 b) { u128 r; r.lo = a.lo + b.lo; r.hi = a.hi + b.hi + (r.lo < a.lo); return r; }
static inline u128 u128_sub(u128 a, u128 b) { u128 r; r.lo = a.lo - b.lo; r.hi = a.hi - b.hi - (a.lo < b.lo); return r; }
static inline int u128_lt(u128 a, u128 b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
static inline u128 u128_shl(u128 a, int s) { u128 r; if (s == 0) return a; if (s >= 64) { r.hi = a.lo << (s - 64); r.lo = 0; return r; } r.hi = (a.hi << s) | (a.lo >> (64 - s)); r.lo = a.lo << s; return r; }
/* (num) / den for a 128-bit num and a 64-bit den, the quotient as 128 bits (long division, 128 steps), remainder out */
static u128 udiv128_64(u128 num, uint64_t den, uint64_t *rem)
{
    u128 qt = {0, 0}; uint64_t r = 0;
    for (int i = 127; i >= 0; i--) {
        const int bit = i >= 64 ? (num.hi >> (i - 64)) & 1 : (num.lo >> i) & 1;
        const int top = (r >> 63) & 1; r = (r << 1) | (uint64_t)bit;
        if (top || r >= den) { r -= den; if (i >= 64) qt.hi |= (uint64_t)1 << (i - 64); else qt.lo |= (uint64_t)1 << i; }
    }
    *rem = r; return qt;
}
/* a signed Q.40 value from an exact 128-bit magnitude numerator over a 64-bit divisor, rounded to nearest (needs the quotient < 2^63) */
static q qround_128_64(u128 mag, int neg, uint64_t den)
{
    uint64_t rem; u128 qt = udiv128_64(mag, den, &rem);
    uint64_t r = qt.lo; if (2 * rem >= den) r++;
    return neg ? -(q)r : (q)r;
}
static void prof_build(int variant)                   /* variant 0 = HYPOTHETICAL (pressure kept inside R1), 1 = AS-PUBLISHED (zeroed there: paper and code) */
{
    RHO = divq(ONE, 7 * PI);
    RMAX = divq(6 * sqrtq(divq(178818 * ONE, 100 * ONE)), 100 * ONE);
    DR = divq(RMAX, (q)(N - 1) * ONE);
    for (int i = 0; i < N; i++) { const q r = (q)i * DR; rho_s[i] = (r > R1 && r < R2) ? RHO : 0; }
    {   /* MT: the trapezoid of 4 pi r^2 times the unsmoothed step, as their cumtrapz; one rounding at the end */
        u128 acc = {0, 0}; q prev = 0;
        for (int i = 0; i < N; i++) {
            const q r = (q)i * DR; const q fi = mulq(mulq(4 * PI, mulq(r, r)), rho_s[i]);
            if (i) { const q s = prev + fi; uint64_t hi, lo; mul128((uint64_t)s, (uint64_t)DR, &hi, &lo); u128 t = {hi, lo}; acc = u128_add(acc, t); }
            prev = fi;
        }
        MT = qround_128_64(acc, 0, (uint64_t)2 << FB);                     /* / (2 * ONE): the trapezoid's half and the Q.40 product's scale */
    }
    SQA = sqrtq(ONE - 2 * MT);
    smooth_passes(rho_s, 7159, 4);
    for (int i = 0; i < N; i++) { const q r = (q)i * DR; f_[i] = mulq(mulq(4 * PI, mulq(r, r)), rho_s[i]); }
    static u128 M2_[N];                                                       /* 2 x the mass trapezoid in Q.80, exact (as MT above); m_ rounded once per node */
    M2_[0].hi = 0; M2_[0].lo = 0; m_[0] = 0;
    for (int i = 1; i < N; i++) { uint64_t hi, lo; mul128((uint64_t)(f_[i - 1] + f_[i]), (uint64_t)DR, &hi, &lo); u128 t = {hi, lo};
        M2_[i] = u128_add(M2_[i - 1], t); m_[i] = qround_128_64(M2_[i], 0, (uint64_t)2 << FB); }
    for (int i = 0; i < N; i++) {
        const q r = (q)i * DR; q P = 0;
        if (r < R2 && !(variant && r < R1)) {
            const q b = sqrtq(ONE - mulq(2 * MT, mulq(r, r)));              /* their closed form with M = MT, R = 1: P = rho (sqrt(1-2M) - sqrt(1-2M r^2)) / (sqrt(1-2M r^2) - 3 sqrt(1-2M)) */
            P = divq(mulq(RHO, SQA - b), b - 3 * SQA);
        }
        P_s[i] = P;
    }
    smooth_passes(P_s, 3999, 4);
    da_[0] = 0;
    for (int i = 1; i < N; i++) {
        const q r = (q)i * DR, r2 = mulq(r, r);
        { const q X = mulq(mulq(4 * PI, mulq(r2, r)), P_s[i]), Dn = r2 - mulq(2 * m_[i], r);   /* the numerator from the exact M2_: (M2_ + 2^41 X) / (2 Dn), rounded once */
          u128 xs = {0, (uint64_t)X}; const u128 N2 = u128_add(M2_[i], u128_shl(xs, FB + 1)); da_[i] = qround_128_64(N2, 0, 2 * (uint64_t)Dn); }
    }
    a_[0] = 0; for (int i = 1; i < N; i++) a_[i] = a_[i - 1] + rdiv(mulq(da_[i - 1] + da_[i], DR), 2);
    {
        const q rend = (q)(N - 1) * DR;
        const q pin = rdiv(lnq(ONE - divq(2 * m_[N - 1], rend)), 2) - a_[N - 1];
        for (int i = 0; i < N; i++) a_[i] += pin;
    }
    /* THE SHIFT, exactly: y at the nodes, two box passes as one triangle kernel in 128-bit, each derived table divided once */
    {
        static q y[N]; static int64_t P1[N + 1]; static u128 P2[N + 1];
        const int w = 1999;                                                  /* span 3999 */
        for (int i = 0; i < N; i++) y[i] = unit_sigmoid((q)i * DR);
        P1[0] = 0; for (int i = 0; i < N; i++) P1[i + 1] = P1[i] + y[i];
        /* B1_i = the first box sum (full window; y is constant where the window would shrink, so the full window is exact there too) */
        static int64_t B1[N];
        for (int i = 0; i < N; i++) { const int lo = i - w < 0 ? 0 : i - w, hi = i + w > N - 1 ? N - 1 : i + w; B1[i] = P1[hi + 1] - P1[lo]; if (hi - lo + 1 < 2 * w + 1) B1[i] = (int64_t)y[i] * (2 * w + 1); }
        P2[0].hi = 0; P2[0].lo = 0;
        for (int i = 0; i < N; i++) { u128 t = {0, (uint64_t)B1[i]}; P2[i + 1] = u128_add(P2[i], t); }
        static u128 B2[N];                                                   /* (2w+1)^2 S_i, exact */
        for (int i = 0; i < N; i++) { const int lo = i - w < 0 ? 0 : i - w, hi = i + w > N - 1 ? N - 1 : i + w; if (hi - lo + 1 < 2 * w + 1) { uint64_t h, l; mul128((uint64_t)B1[i], (uint64_t)(2 * w + 1), &h, &l); B2[i].hi = h; B2[i].lo = l; } else B2[i] = u128_sub(P2[hi + 1], P2[lo]); }
        const uint64_t K2 = (uint64_t)(2 * w + 1) * (uint64_t)(2 * w + 1);
        SHDEV = 0;
        for (int i = 0; i < N; i++) {
            sh_[i] = qround_128_64(B2[i], 0, K2);                            /* S_i, one rounding */
            q d = sh_[i] - y[i]; if (d < 0) d = -d; /* against the raw sigmoid: the smoothing's own change, printed as a sanity number */
            if (d > SHDEV) SHDEV = d;
        }
        /* bp_i = S'_i = (B2_{i+1} - B2_{i-1}) * ONE / (2 dr K2), one division; the ends by the same rule with a one-sided difference */
        const uint64_t den1 = 2 * (uint64_t)DR; /* times K2 below via a two-step: den fits 64 bits as 2 dr_q K2 (9e14) */
        uint64_t den_bp_hi, den_bp_lo; mul128(den1, K2, &den_bp_hi, &den_bp_lo); const uint64_t den_bp = den_bp_lo;   /* 2 dr K2 < 2^63 */
        (void)den_bp_hi;
        for (int i = 0; i < N; i++) {
            const int ip = i + 1 < N ? i + 1 : i, im = i > 0 ? i - 1 : i;
            u128 dB; int neg; if (u128_lt(B2[ip], B2[im])) { dB = u128_sub(B2[im], B2[ip]); neg = 1; } else { dB = u128_sub(B2[ip], B2[im]); neg = 0; }
            const int steps = ip - im;                                        /* 2 in the interior, 1 at the ends */
            u128 num = u128_shl(dB, FB);
            bp_[i] = steps == 2 ? qround_128_64(num, neg, den_bp) : qround_128_64(num, neg, den_bp / 2);
        }
        /* spp_i = S'' on interval i = (bp_{i+1} - bp_i) / dr taken from the EXACT sums: ((B2_{i+2} - B2_{i+1}) - (B2_i - B2_{i-1})) ONE^2 / (2 dr^2 K2)
           in two stages: first / K2 with a 128-bit quotient (exact to the remainder), then / (2 dr^2) rounded once */
        uint64_t d2hi, d2lo; mul128((uint64_t)DR, (uint64_t)DR, &d2hi, &d2lo); (void)d2hi; /* VERIFY FIX: DR^2 kept in Q.80 (7.78e14, d2hi = 0) */
        const uint64_t den2 = 2 * d2lo;
        for (int i = 0; i + 1 < N; i++) {
            const int a0 = i - 1 >= 0 ? i - 1 : 0, a1 = i, a2 = i + 1, a3 = i + 2 < N ? i + 2 : N - 1;
            /* dd = (B2[a3] - B2[a2]) - (B2[a1] - B2[a0]) as a signed 128-bit */
            u128 p, m; int pn, mn;
            if (u128_lt(B2[a3], B2[a2])) { p = u128_sub(B2[a2], B2[a3]); pn = 1; } else { p = u128_sub(B2[a3], B2[a2]); pn = 0; }
            if (u128_lt(B2[a1], B2[a0])) { m = u128_sub(B2[a0], B2[a1]); mn = 1; } else { m = u128_sub(B2[a1], B2[a0]); mn = 0; }
            u128 dd; int neg;
            if (pn == mn) { if (u128_lt(p, m)) { dd = u128_sub(m, p); neg = !pn; } else { dd = u128_sub(p, m); neg = pn; } }
            else { dd = u128_add(p, m); neg = pn; }
            /* the end intervals have one-sided steps: scale by the step counts (2 and 2 in the interior) */
            const int s1 = a3 - a2, s0 = a1 - a0; (void)s1; (void)s0;
            /* S''_q = (dd / K2) * ONE / (2 dr2_q), dr2_q = dr^2 in Q.40 units: ONE factor of 2^40 (dd / K2 is the second difference of S in Q.40 units) */
            u128 num = u128_shl(dd, FB);                                      /* dd * ONE: < 2^128 since dd < 2^88 */
            uint64_t rem1; u128 q1 = udiv128_64(num, K2, &rem1);             /* stage 1: / K2, exact to the remainder (relative 1e-16 at most) */
            spp_[i] = qround_128_64(u128_shl(q1, FB), neg, den2);  /* VERIFY FIX: q1 2^40 / (2 DR^2) */                           /* stage 2: / (2 dr^2), rounded once */
            if (a0 == a1 || a3 == a2) spp_[i] = 0;                            /* the grid's first and last interval: S is constant there */
        }
        spp_[N - 1] = 0;
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
    o->Spp = spp_[i];
    o->Sp = bp_[i] + mulq(spp_[i], h2);
    o->S = sh_[i] + mulq(bp_[i], h2) + mulq(spp_[i], mulq(h2, h2) / 2);
}
#endif
