/* wh_relaxt.c (11-12-2026, 3:20 AM) -- WORMHOLE MATH 1b, THE P0 KERNEL: the relaxation of several coupled fields on the (x, t) grid of the
 * patched geometry's first problem, in the form the reviewer's reconciled sizing (3:14 AM) prescribes, built on an exact test ahead of the
 * equations so that P0's three sources drop into one callback when its specification lands.
 *
 * THE OPERATOR. Rows uniform in the log depth t = ln(-y) from the cut t0 (Dirichlet data there, the conformal boundary's) down to t1 = ln 14
 * (Neumann with a given normal derivative, the deep end); columns uniform in x with spacing 1 (free for an x-independent problem), Neumann
 * with given values on both x sides. For each field F_k:
 *     F_tt - F_t + e^{2t} F_xx = S_k(F, t, x),
 * discretized with centred differences (exact for quadratics in t and in x) and multiplied through by Dt^2, so that every coefficient is
 * below the Q.60 ceiling: (1 + Dt/2) F_S + (1 - Dt/2) F_N + c_j (F_W + F_E) - (2 + 2 c_j) F = Dt^2 S_k, with c_j = e^{2 t_j} Dt^2 / Dx^2 (at
 * t1 = ln 14: 1.14 at N = 128, 4.4 at N = 64), and then divided by its diagonal: the Gauss-Seidel correction is the weighted sum
 * d = w_S F_S + w_N F_N + w_X (F_W + F_E) - w_D S - F with the per-row weights w_S = (1 + Dt/2)/(2 + 2 c_j), w_N = (1 - Dt/2)/(2 + 2 c_j),
 * w_X = c_j/(2 + 2 c_j), w_D = Dt^2/(2 + 2 c_j), all below one at any N, rounded once at setup; no division per point. The residual of the
 * operator in units of 2^-60 is d (2 + 2 c_j)/Dt^2, formed in 128 bits and saturating. Ghost points carry the Neumann values:
 * F_{-1} = F_1 - 2 Dt g on a t side, F_{-1} = F_1 - 2 Dx g on an x side.
 * THE SWEEP. Red-black over the points; at each point all fields are updated from the same state (never one field alone: the reviewer's
 * spectra show the warp factor's equation alone has four negative eigenvalues and diverges in a few sweeps, while the coupled column is held
 * by its first-order cross-couplings). The plain diagonal with a per-row reciprocal r_j = 1/(2 + 2 c_j), rounded once at setup: the
 * Gauss-Seidel correction d_k = r_j R_k with R_k the row-scaled residual, the update omega d_k; the sources lagged (the sizing: the full
 * Newton diagonal costs 22-180 percent near the ceiling and buys nothing; the one-sided one is equally good at 1.7 times the cost).
 * omega = 2 - 90/N by default (the linear iteration grows the error 100-750 times before it decays, so the safe value sits well under the
 * linear edge), with THE CHANGE GUARD: a checkpoint of the fields every 64 sweeps when the sweep's largest change is at or below the last
 * checkpoint's; if a sweep's largest change exceeds 32 times the smallest seen, the checkpoint is restored and omega drops by 0.1.
 * THE STOPS. The largest Gauss-Seidel correction at or below tol (2^-40 by default, in field units); the largest residual, exact in 128
 * bits and divided back by Dt^2, at or below tol, checked every 16 sweeps; and, when an exact solution is known, its largest deviation.
 * THE TEST (this file): F_k = P_k(t) Q(x), P_k quadratic in t, Q quadratic in x, with S_k the operator applied to it symbolically and then
 * evaluated on the grid in the same arithmetic (e^{2t} rounded once), so the discrete solution equals F_k at the nodes to the arithmetic,
 * plus an optional lagged nonlinear coupling lam (F_k^2 + SUM_{m != k} F_m) to exercise the plain update as P0's sources will. MEASURED
 * (3:37 AM): the operator, weights, ghosts, guard and stops land on the exact solution from zero at N = 64-256 (errors 4e-14 to 5e-13 at the
 * 2^-40 stop; started at the exact solution the correction is 3 units); lam = 0.1 converges likewise; lam = 1 from zero leaves the basin
 * (the test's solution is 0.04 in size and the transient overshoots to 0.3, where a quadratic source has no nearby solution), which no
 * guard cures: a nonlinear problem starts near its solution, as P0 will from its reference.
 * MODEL=ads3 (3:45 AM): THE P0 SOURCES, the reviewer's three equations in the t-form (stage 1's verified system, y = -e^t, times e^{2t}):
 *   U_tt - U_t + e^{2t} U_xx = -[U_t (U_t + 2 V_t) + e^{2t} U_x (U_x + 2 V_x)] + e^{2(B+t)} [4/L^2 + (Q2/3) e^{-4V}]
 *   V_tt - V_t + e^{2t} V_xx = -[V_t (U_t + 2 V_t) + e^{2t} V_x (U_x + 2 V_x)] + e^{2(B+t)} [4/L^2 + kappa e^{-2V} - (2/3) Q2 e^{-4V}]
 *   B_tt - B_t + e^{2t} B_xx =  [2 U_t V_t + V_t^2 + e^{2t} (2 U_x V_x + V_x^2)] + e^{2(B+t)} [-2/L^2 - kappa e^{-2V} + (5/6) Q2 e^{-4V}]
 * in DEVIATION VARIABLES u, v, b around a reference (Ur(t), Vr, Br(t)), because the raw fields pass Q.60 (U near 7 at the cut, B near -13
 * at the deep end): the reference's node values are held in Q.56, its t-derivative Ur' and the factor e^{2(Br + t)} are formed per row
 * through exponent arithmetic in 128 bits (e^t is split into an integer and a fraction, since it reaches 14), and the operator's action on
 * the reference is the per-row constant g_j = w_S (Ur_{j-1} - Ur_j) + w_N (Ur_{j+1} - Ur_j) (likewise for Br). The weighted source splits as
 * w_D S_a + w_X S_b + w_E e^{2b} S_c with S_a the t-derivative terms, S_b the x-derivative terms (their e^{2t} is w_X's), S_c the bracket,
 * w_E = Dt^2 e^{2(Br+t)}/(2 + 2 c_j); the exponentials of the deviations (e^{2b}, e^{-2v}, e^{-4v}) are the only ones evaluated per point.
 * THE TEST is global AdS3 x S2 with flux (the reviewer's, 3:22 AM): L = 1, kappa = 1, w = 9/8 = 1/S0^2, Q2 = (3w + 12)/(2w^2), L3 =
 * (3 + w/4)^(-1/2), U = ln(L3 coth e^t), B = ln(L3 / sinh e^t), V = ln S0, which satisfies the three equations exactly (the reviewer: to
 * 1e-124). The deviations are then zero, and to make zero the EXACT discrete solution the reference's discretization defect D_k = g_k -
 * (its weighted source) is subtracted once (the manufactured-solution form; in the real P0 that defect is the source that drives the
 * deviation from a reference that is not a solution). Dirichlet zero at the cut, Neumann zero at the deep end and on the x sides. A start
 * pert= away from zero must relax back to zero: that checks the sources, the stencil, the ghosts and the guard in P0's own variables.
 * Usage: whrelaxt [model=quad|ads3|p0] [N=128] [Gx=4] [lam=0] [pert=0.01] [defect=1] [Q2=6] [a=0.02] [fu4=0] [fb4=0] [fv4=0] [ramp=10] [scale=1] [omega=auto|value] [tol=1048576] [t0=-7.0] [guard=1]    N 4..1024 (ads3: 64 or more; defect=0 leaves the reference's discretization defect in); Gx up to 64. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "whq60.h"
#include "par2.h"

#define MAXN 1024
#define MAXGX 64
#define NF 3
static int N, GX; static int64_t DT, DX, T0, OMEGA, LAM;
#define ONE56 ((int64_t)1 << 56)
static int64_t T056, DT56;   /* the row coordinate t_j = t0 + Dt j in Q.56: t itself passes the Q.60 ceiling at the real cut (t0 = -10.85) and the span is 13.5 */
static int64_t UR56[MAXN + 2], BR56[MAXN + 2], LNS0, LNL3, CZ[MAXN + 2];   /* CZ: coth(z)/z per row where z >= 0.4 (else -1), for the first-integral rows */   /* the reference's node values in Q.56 (U about 11 and B about -13 at the ends), one ghost row past the deep end; ln S0 and ln L3 */
static int64_t F[NF][MAXGX + 1][MAXN + 1], SRC[NF][MAXGX + 1][MAXN + 1], EX[NF][MAXGX + 1][MAXN + 1], CKP[NF][MAXGX + 1][MAXN + 1];
static int64_t CJ[MAXN + 1], WS[MAXN + 1], WN[MAXN + 1], WX[MAXN + 1], WD[MAXN + 1], DIAG1[MAXN + 1];   /* per row: c_j, the weights, and 1 + c_j */
static inline void split60(s128 x, int64_t *n, int64_t *f)         /* a 128-bit Q.60 exponent into integer + fraction for exp60n */
{ *n = (int64_t)(((uint64_t)x.hi << 4) | (x.lo >> 60)); *f = (int64_t)(x.lo & (((uint64_t)1 << 60) - 1)); }
#define TROW56(j) (T056 + DT56 * (int64_t)(j))   /* t_j = t0 + Dt j in Q.56 (|t| below 128: no wrap) */
static inline s128 s128_neg(s128 a) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); return a; }
/* ln(1 - e^{-c z}) in 128 bits, c = 2 or 4, z = zi + zf (Q.60), em = e^{-c z}: for z >= 1 the direct ln60; below 1 as ln c + t + ln q with q = (1 - e^{-x})/x,
   x = c z (the series for x < 1/8, else the quotient), since near the cut the value is about t + ln c, past the Q.60 ceiling */
static s128 ln1mexp(int c, int64_t zi, int64_t zf, int64_t em, int64_t t56)
{
    s128 r = {0, 0};
    if (zi >= 1) { s128_addv(&r, ln60(ONE60 - em)); return r; }
    const int64_t x = c * zf; int64_t q;
    if (x < ONE60 / 8) { int64_t term = ONE60; q = ONE60; for (int k = 1; k <= 14; k++) { term = -rdiv64(mul60(term, x), k + 1); if (!term) break; q += term; } }
    else q = div60(ONE60 - em, x);
    s128_addv(&r, (c == 2 ? 1 : 2) * ln2_60()); s128_addp(&r, 16, t56); s128_addv(&r, ln60(q));
    return r;
}
static int64_t parseq(const char *s, int fb)   /* a decimal into Q.fb (fb = 56: |v| below 128; fb = 48: below 32768): the integer part, then the fraction through parse60 */
{
    int neg = 0; if (*s == '-') { neg = 1; s++; } else if (*s == '+') s++;
    int64_t ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (*s++ - '0');
    int64_t fr = 0; if (*s == '.') { char buf[40] = "0."; strncat(buf, s + 1, 30); fr = (parse60(buf) + ((int64_t)1 << (59 - fb))) >> (60 - fb); }
    const int64_t v = (ip << fb) + fr; return neg ? -v : v;
}
static int64_t parse56(const char *s) { return parseq(s, 56); }
static void fmt56(int64_t v, char *out)   /* twelve decimals of a Q.56 value */
{
    const uint64_t u = v < 0 ? -(uint64_t)v : (uint64_t)v; const uint64_t ip = u >> 56, fr = u & ((uint64_t)ONE56 - 1);
    uint64_t h, l; mul128(fr, 1000000000000ull, &h, &l); uint64_t q = (h << 8) | (l >> 56); if (l & ((uint64_t)1 << 55)) q++;
    uint64_t ipr = ip; if (q >= 1000000000000ull) { q -= 1000000000000ull; ipr++; }
    sprintf(out, "%s%llu.%012llu", v < 0 ? "-" : "", (unsigned long long)ipr, (unsigned long long)q);
}
static int64_t wofq2(int64_t Q2)   /* w = 1/S0^2 from 2 Q2 w^2 - 3 w - 12 = 0: w = 3/(4 Q2) + 4 s/Q2 with s = sqrt(9/256 + 3 Q2/8) (every piece below 8) */
{ const int64_t s = sqrt60(9 * (ONE60 / 256) + mul60(Q2, 3 * (ONE60 / 8))); return div60(3 * (ONE60 / 4), Q2) + div60(4 * s, Q2); }
static int64_t GT1[NF][MAXGX + 1], GX0[NF][MAXN + 1], GX1[NF][MAXN + 1];   /* the Neumann values: d_t F at t1 per column; d_x F at x0 and x1 per row */
static int MODEL = 0;                                                     /* 0 quad, 1 ads3 */
static int DEFECT = 1;   /* ads3: 1 subtracts the manufactured defect (zero deviation is then the exact discrete solution); 0 leaves it in, so the end deviation is the reference's discretization error on this grid */
static int64_t URP[MAXN + 1], WE[MAXN + 1], GU[MAXN + 1], GB[MAXN + 1], DEF[NF][MAXN + 1];   /* ads3 per row: Ur', w_E, the operator on Ur and Br, the manufactured defects */
static int64_t INV2DT56, K_C1, K_C2, K_C3, K_C4, K_C5, K_C6;             /* 1/(2 Dt) in Q.56; the bracket constants (Q2/3) w^2, kappa w, (2/3) Q2 w^2, (5/6) Q2 w^2 */
static inline int64_t exp_of128(s128 x) { int64_t n, f; split60(x, &n, &f); return n < -62 ? 0 : exp60n(n, f); }
static int64_t AS, AN;                                             /* 1 + Dt/2 and 1 - Dt/2 */
/* the exact test: P_k(t) = (p0_k + p1_k s + p2_k s^2)/16 with s = t/8 (so |s| < 1 and nothing passes the ceiling: t^2 itself would be 49 at the
   cut), Q(x) = 1 + x/8 - (x/8)^2 (Q' = 1/8 - x/32, Q'' = -1/32); the sources below */
static const int64_t PC[NF][3] = { { 1, 1, -1 }, { 2, -1, 1 }, { -1, 2, 1 } };     /* times ONE60 / 16 */
static inline int64_t Pk(int k, int64_t sv) { const int64_t u = ONE60 / 16; return PC[k][0] * u + mul60(PC[k][1] * u, sv) + mul60(PC[k][2] * u, mul60(sv, sv)); }   /* sv = t/8 in Q.60 = 2 t in Q.56 */
static inline int64_t Pkp(int k, int64_t sv) { const int64_t u = ONE60 / 16; return (PC[k][1] * u + 2 * mul60(PC[k][2] * u, sv)) / 8; }   /* dP/dt = (p1 + 2 p2 s)/(16 x 8) */
static inline int64_t Pkpp(int k) { return 2 * PC[k][2] * (ONE60 / 16) / 64; }                                                                  /* d2P/dt2 = 2 p2/(16 x 64) */
static inline int64_t Qx(int64_t x) { const int64_t x8 = x / 8; return ONE60 + x8 - mul60(x8, x8); }
static inline int64_t Qxp(int64_t x) { return ONE60 / 8 - x / 32; }
static inline int64_t Qxpp(void) { return -ONE60 / 32; }
/* the field value with the ghosts: j = -1 does not occur (Dirichlet at the cut); j = N + 1 mirrors with the t1 value; i = -1 and GX + 1 mirror with the x values */
static inline int64_t Fv(int k, int i, int j)
{
    if (j > N) return F[k][i][N - 1] + 2 * mul60(DT, GT1[k][i]);
    if (i < 0) return F[k][1][j] - 2 * mul60(DX, GX0[k][j]);
    if (i > GX) return F[k][GX - 1][j] + 2 * mul60(DX, GX1[k][j]);
    return F[k][i][j];
}
static inline int64_t source(int k, int i, int j)                  /* S_k = the test's operator value + lam (F_k^2 + SUM_{m != k} F_m), lagged */
{
    int64_t s = SRC[k][i][j];
    if (LAM) { s += mul60(LAM, mul60(F[k][i][j], F[k][i][j])); for (int m = 0; m < NF; m++) if (m != k) s += mul60(LAM, F[m][i][j]); }
    return s;
}
static inline int64_t dT(int k, int i, int j) { s128 r = {0, 0}; s128_addp(&r, Fv(k, i, j + 1) - Fv(k, i, j - 1), INV2DT56); return s128_round_sh(r, 56); }   /* centred d_t of a deviation */
static inline int64_t dX(int k, int i, int j) { return (Fv(k, i + 1, j) - Fv(k, i - 1, j)) / 2; }                                                         /* centred d_x (Dx = 1) */
/* The first integral of the x-independent problem (the reviewer's stage 2): Q = e^{U + 2V} (d_y U - d_y B) is constant in y; in t = ln(-y),
   Q = -e^{U + 2V - t} (U_t - B_t), whose reference part is exactly 1, so Q/(-L3 S0^2) = e^{u + 2v} [1 + (coth z / z)(u_t - b_t)]: the same on
   every row for the exact solution and e^{u_c + 2 v_c} = c2_eff/c2 at the deep end. Near the cut coth z / z = 1/z^2 is 1e9 and the row is ill-conditioned. */
/* ads3: the three weighted sources at a point: w_D S_a + w_X S_b + w_E e^{2b} S_c (the manufactured defect is subtracted in correction()) */
static void ads3_wsource(int i, int j, int64_t out[NF])
{
    const int64_t ut = URP[j] + dT(0, i, j), vt = dT(1, i, j), ux = dX(0, i, j), vx = dX(1, i, j), b = F[2][i][j], v = F[1][i][j];
    const int64_t e2b = exp60(2 * b); s128 mv = {0, 0}; s128_addp(&mv, -2, v); const int64_t em2v = exp_of128(mv), em4v = mul60(em2v, em2v);   /* v reaches 7 at the real cut: -2v and -4v pass the Q.60 ceiling, so the exponent is formed in 128 bits and e^{-4v} is the square */
    const int64_t sa[NF] = { -mul60(ut, ut + 2 * vt), -mul60(vt, ut + 2 * vt), 2 * mul60(ut, vt) + mul60(vt, vt) };
    const int64_t sb[NF] = { -mul60(ux, ux + 2 * vx), -mul60(vx, ux + 2 * vx), 2 * mul60(ux, vx) + mul60(vx, vx) };
    const int64_t sc[NF] = { 4 * ONE60 + mul60(K_C1, em4v), 4 * ONE60 + mul60(K_C2, em2v) - mul60(K_C3, em4v), -2 * ONE60 - mul60(K_C2, em2v) + mul60(K_C4, em4v) };
    for (int k = 0; k < NF; k++) {
        s128 r = {0, 0}; s128_addp(&r, WD[j], sa[k]); s128_addp(&r, WX[j], sb[k]); s128_addp(&r, mul60(WE[j], e2b), sc[k]);
        out[k] = s128_round_sh(r, F60);
    }
}
static inline int64_t correction(int k, int i, int j)            /* d_k = w_S F_S + w_N F_N + w_X (F_W + F_E) - w_D S_k - F_k, rounded once */
{
    s128 r = {0, 0};
    s128_addp(&r, WS[j], Fv(k, i, j - 1)); s128_addp(&r, WN[j], Fv(k, i, j + 1));
    s128_addp(&r, WX[j], Fv(k, i - 1, j) + Fv(k, i + 1, j));
    if (MODEL == 0) { s128_addp(&r, -WD[j], source(k, i, j)); return s128_round_sh(r, F60) - F[k][i][j]; }
    int64_t ws[NF]; ads3_wsource(i, j, ws);                          /* the operator's action on the reference, minus its weighted source */
    return s128_round_sh(r, F60) - ws[k] + (k == 0 ? GU[j] : k == 2 ? GB[j] : 0) - (DEFECT ? DEF[k][j] : 0) - F[k][i][j];
}
static int64_t DT2HALF;                                          /* Dt^2 / 2 in Q.60 */
static inline int64_t residual_units(int64_t d, int j)           /* |F_tt - F_t + e^{2t} F_xx - S| = |d| (2 + 2 c_j) / Dt^2, in units of 2^-60, saturating */
{
    const uint64_t ad = (uint64_t)(d < 0 ? -d : d);                   /* |d| / w_D = |d| (2 + 2 c_j) / Dt^2, in units: (|d| 2^60) / w_D */
    const uint64_t hi = ad >> 4, lo = ad << 60;
    if (hi >= (uint64_t)WD[j]) return INT64_MAX;
    uint64_t rem; uint64_t q = div128(hi, lo, (uint64_t)WD[j], &rem); if (2 * rem >= (uint64_t)WD[j]) q++;
    return q > (uint64_t)INT64_MAX ? INT64_MAX : (int64_t)q;
}
typedef struct { int color, i_begin, i_end; int64_t maxchg, maxcorr, maxres; } rjob;
static void sweep_rows(void *p)
{
    rjob *jb = (rjob *)p; const int color = jb->color; int64_t maxchg = 0, maxcorr = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j <= N; j++) if (((i + j) & 1) == color) {
        int64_t d[NF];
        for (int k = 0; k < NF; k++) d[k] = correction(k, i, j);
        for (int k = 0; k < NF; k++) if (d[k] > 3 * ONE60 || d[k] < -3 * ONE60) { printf("WHRELAXT: a correction passed 3 at column %d row %d: the start is outside the basin -- stop\n", i, j); exit(2); }
        for (int k = 0; k < NF; k++) {
            const int64_t c = d[k] < 0 ? -d[k] : d[k]; if (c > maxcorr) maxcorr = c;
            const int64_t step = mul60(OMEGA, d[k]); const int64_t ch = step < 0 ? -step : step; if (ch > maxchg) maxchg = ch;
            F[k][i][j] += step;
        }
    }
    jb->maxchg = maxchg; jb->maxcorr = maxcorr;
}
static void resid_rows(void *p)
{
    rjob *jb = (rjob *)p; int64_t maxres = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j <= N; j++) for (int k = 0; k < NF; k++) {
        const int64_t rv = residual_units(correction(k, i, j), j); if (rv > maxres) maxres = rv;
    }
    jb->maxres = maxres;
}
static int64_t P0_W, P0_Q2V, P0_FU4, P0_FB4, P0_FV4, P0_SCALE = ONE60;   /* the three free v^4 coefficients (U, B, V) in Q.48: thousands at a = 1/8 */
static s128 ln_s128(int64_t x)   /* ln x for any positive Q.60 x, in 128 bits: x = 2^-k m with m in [1/2, 1), ln x = ln m - k ln 2 (ln60 alone wraps below e^-8) */
{ int k = 0; while (x < ONE60 / 2) { x <<= 1; k++; } s128 r = {0, 0}; s128_addv(&r, ln60(x)); s128_addp(&r, -k, ln2_60()); return r; }
/* the brane's data at the cut for the S2 radius a, as deviations from the reference (the large parts cancel in 128 bits first); the small FG terms
   (y^2/(18 a^2), 7 y^2/(36 a^2), the v^4 ln v term, the free v^4 terms) only when full (at a ramp's smaller radii the series is not valid and
   v^2/a^4 would pass the ceiling); returns 0 if V's deviation passes the Q.60 ceiling */
static int p0_cutdata(int64_t a, int full, int64_t *u0, int64_t *v0, int64_t *b0, int64_t *logt_out)
{
    int64_t fg_ub = 0, fg_v = 0, logt = 0, f4u = 0, f4b = 0, f4v = 0;
    if (full) {
        const int64_t a2 = mul60(a, a), a4 = mul60(a2, a2);
        s128 ev = {0, 0}; s128_addp(&ev, 16, T056); const int64_t v = exp_of128(ev), v2 = mul60(v, v), v4 = mul60(v2, v2);   /* v = -y = e^{t0} */
        fg_ub = div60(v2, 18 * a2); fg_v = div60(7 * v2, 36 * a2);                 /* v^2/(18 a^2) and 7 v^2/(36 a^2) */
        const int64_t v2lnv = 16 * mul60(v2, T056);                                 /* v^2 ln v, with ln v = t0 in Q.56 (hence the 16) */
        const int64_t cl_v2 = mul60(rdiv64(P0_Q2V, 10) - ONE60 / 60, div60(v2, a4));   /* (6 Q2 - 1)/60 = Q2/10 - 1/60, times v^2/a^4 (6 Q2 itself, and 10 as a Q.60 divisor, would pass the ceiling) */
        logt = mul60(cl_v2, v2lnv);                                                 /* (6 Q2 - 1)/(60 a^4) v^4 ln v: U and B carry it, V carries -3/2 of it */
        { s128 r = {0, 0}; s128_addp(&r, P0_FU4, v4); f4u = s128_round_sh(r, 48); }     /* a Q.48 coefficient times a Q.60 v^4: rounded once to Q.60 */
        { s128 r = {0, 0}; s128_addp(&r, P0_FB4, v4); f4b = s128_round_sh(r, 48); }
        { s128 r = {0, 0}; s128_addp(&r, P0_FV4, v4); f4v = s128_round_sh(r, 48); }
    }
    s128 vs = ln_s128(a); s128_addp(&vs, -16, T056); s128_addv(&vs, -LNS0);          /* V0 - ln S0 without the small terms: ln a - t0 - ln S0 (7.00 at the reconciled cut and a = 0.02) */
    if (vs.hi > 0 || vs.hi < -1 || (vs.hi == 0 && vs.lo >= (uint64_t)(8 * ONE60)) || (vs.hi == -1 && vs.lo < (uint64_t)(-8 * ONE60))) return 0;
    const int64_t vlarge = (int64_t)vs.lo;                                           /* within (-8, 8): the low word is the value */
    int64_t u = ((-T056 - UR56[0]) << 4) + fg_ub + logt + f4u, b = ((-T056 - BR56[0]) << 4) + fg_ub + logt + f4b, vv = vlarge - fg_v - rdiv64(3 * logt, 2) + f4v;
    if (P0_SCALE != ONE60) { u = mul60(u, P0_SCALE); vv = mul60(vv, P0_SCALE); b = mul60(b, P0_SCALE); }
    *u0 = u; *v0 = vv; *b0 = b; *logt_out = logt; return 1;
}
/* the start profile of the V deviation for the radius a at row j: (1/2) ln(1 + a^2 w e^{-2t}), the shape V = (1/2) ln(S0^2 + a^2 e^{-2t}) has
   (ln a - t - ln S0 near the cut, slope -1 down to t = (1/2) ln(a^2 w), then zero); a straight line left the basin in 63 sweeps */
static int64_t p0_profile(int64_t a, int j)
{
    s128 xx = ln_s128(a); xx.hi = xx.hi * 2 + (int64_t)(xx.lo >> 63); xx.lo <<= 1;   /* 2 ln a */
    s128_addv(&xx, ln60(P0_W)); s128_addp(&xx, -32, TROW56(j));                         /* x = ln(a^2 w) - 2t */
    int64_t r;
    if (xx.hi > 0 || (xx.hi == 0 && xx.lo > 0)) r = s128_round_sh(xx, 1) + ln60(ONE60 + exp_of128(s128_neg(xx))) / 2;   /* x/2 + (1/2) ln(1 + e^-x) */
    else r = ln60(ONE60 + exp_of128(xx)) / 2;                                                                      /* (1/2) ln(1 + e^x) */
    return P0_SCALE != ONE60 ? mul60(r, P0_SCALE) : r;
}
int main(int argc, char **argv)
{
    N = 128; GX = 4; int guard = 1, start_exact = 0, dbg = 0; int64_t lam = 0, omega_in = 0, tol = (int64_t)1 << 20, t0 = INT64_MIN, pert = 0, Q2in = 6 * ONE60, ain = parse60("0.02"), fu4 = 0, fb4 = 0, fv4 = 0, scale = ONE60; int ramp = 0;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "Gx=", 3)) GX = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "lam=", 4)) lam = parse60(argv[i] + 4);
        else if (!strncmp(argv[i], "tol=", 4)) tol = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "t0=", 3)) t0 = parse56(argv[i] + 3);
        else if (!strncmp(argv[i], "Q2=", 3)) Q2in = parse60(argv[i] + 3);
        else if (!strncmp(argv[i], "a=", 2)) ain = parse60(argv[i] + 2);
        else if (!strncmp(argv[i], "fu4=", 4)) fu4 = parseq(argv[i] + 4, 48);   /* the free v^4 coefficients of U, B and V (the energy and the pressure of the brane's state), Q.48 */
        else if (!strncmp(argv[i], "fb4=", 4)) fb4 = parseq(argv[i] + 4, 48);
        else if (!strncmp(argv[i], "fv4=", 4)) fv4 = parseq(argv[i] + 4, 48);
        else if (!strncmp(argv[i], "scale=", 6)) scale = parse60(argv[i] + 6);
        else if (!strncmp(argv[i], "ramp=", 5)) ramp = atoi(argv[i] + 5);
        else if (!strncmp(argv[i], "guard=", 6)) guard = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "start=", 6)) start_exact = !strcmp(argv[i] + 6, "exact");
        else if (!strncmp(argv[i], "dbg=", 4)) dbg = atoi(argv[i] + 4);
        else if (!strncmp(argv[i], "model=", 6)) MODEL = !strcmp(argv[i] + 6, "ads3") ? 1 : !strcmp(argv[i] + 6, "p0") ? 2 : 0;
        else if (!strncmp(argv[i], "pert=", 5)) pert = parse60(argv[i] + 5);
        else if (!strncmp(argv[i], "defect=", 7)) DEFECT = atoi(argv[i] + 7);
        else if (!strncmp(argv[i], "omega=", 6)) { if (strcmp(argv[i] + 6, "auto")) omega_in = parse60(argv[i] + 6); }
    }
    if (N < 4 || N > MAXN || (N & (N - 1))) { printf("WHRELAXT: N must be a power of two in 4..%d\n", MAXN); return 1; }
    if (GX < 2 || GX > MAXGX) { printf("WHRELAXT: Gx must be 2..%d\n", MAXGX); return 1; }
    char a[48], b[48], c[48], d[48], e[48];
    if (MODEL == 2 && ramp == 0) ramp = 10;                            /* the real P0 climbs to its brane radius in ten doublings unless told otherwise */
    if (t0 == INT64_MIN) t0 = MODEL == 2 ? parse56("-10.85") : -7 * ONE56;   /* the cut: the exact tests at t0 = -7; the real P0 at y = -1.94e-5 (the reviewer's reconciled depth) */
    const int64_t t1 = ln60(14 * (ONE60 / 16)) + 4 * ln2_60();          /* ln 14 = ln(14/16) + 4 ln 2 (14 itself is past the ceiling) */
    T056 = t0; DT56 = rdiv64(((t1 + 8) >> 4) - t0, N); DT = DT56 << 4; T0 = 0; DX = ONE60; LAM = lam;   /* the rows in Q.56: the span t1 - t0 (9.6, or 13.5 at the real cut) passes the Q.60 ceiling; Dt itself is below 8 and goes to Q.60 exactly */
    AS = ONE60 + DT / 2; AN = ONE60 - DT / 2; DT2HALF = mul60(DT, DT) / 2;
    for (int j = 0; j <= N; j++) {                                    /* per row: t_j, e^{2 t_j}, c_j, the reciprocal */
        const int64_t tj56 = TROW56(j);
        /* c_j = e^{2t} Dt^2 / Dx^2 = e^{2t + 2 ln Dt} (Dx = 1): one exponential, below 2 at N >= 128 (e^{2t} itself reaches 196, past the ceiling) */
        s128 ex2 = {0, 0}; s128_addp(&ex2, 32, tj56); s128_addp(&ex2, 2, ln60(DT)); int64_t n2, f2; split60(ex2, &n2, &f2);
        if (ex2.hi < 0) {                                                 /* c_j = e^g below 1: the direct form; the diagonal 2 (1 + c_j) */
            CJ[j] = exp60n(n2, f2); DIAG1[j] = ONE60 + CJ[j];
            WS[j] = div60(AS / 2, DIAG1[j]); WN[j] = div60(AN / 2, DIAG1[j]); WX[j] = div60(CJ[j] / 2, DIAG1[j]); WD[j] = div60(mul60(DT, DT) / 2, DIAG1[j]);
        } else {                                                          /* c_j at or above 1 (coarse grids): through 1/c = e^-g, which is below 1 */
            s128 neg = {0, 0}; s128_addp(&neg, -32, tj56); s128_addp(&neg, -2, ln60(DT)); int64_t n3, f3; split60(neg, &n3, &f3);
            const int64_t rc = exp60n(n3, f3), den = ONE60 + rc;           /* 1/c and 1 + 1/c */
            CJ[j] = rc > ONE60 / 8 ? div60(ONE60, rc) : -1; DIAG1[j] = 0;   /* c itself only where it fits (the banner shows -1 past the ceiling); the weights never need it */
            WS[j] = div60(mul60(AS / 2, rc), den); WN[j] = div60(mul60(AN / 2, rc), den); WX[j] = div60(ONE60 / 2, den); WD[j] = div60(mul60(mul60(DT, DT) / 2, rc), den);
        }
    }
    { uint64_t rem; INV2DT56 = (int64_t)div128((uint64_t)1 << (56 + 60 - 64), 0, (uint64_t)(2 * DT), &rem); if (2 * rem >= (uint64_t)(2 * DT)) INV2DT56++; }   /* 2^116 / (2 Dt 2^60) = 2^56/(2 Dt) */
    if (MODEL >= 1) {                                                /* the AdS3 x S2 reference, with the manufactured defects (ads3) or the real P0's data (p0) */
        if (N < 64) { printf("WHRELAXT: model=ads3|p0 needs N >= 64\n"); return 1; }
        const int64_t w = MODEL == 2 ? wofq2(Q2in) : parse60("1.125"), w2 = mul60(w, w);   /* w = 1/S0^2: 9/8 for the exact test; for the real P0 the root of 2 Q2 w^2 - 3 w - 12 = 0 (1.1328 at Q2 = 6), so the reference carries the real sources */
        const int64_t Q2v = div60(3 * w / 2 + 6 * ONE60, w2);              /* Q2 = (3w + 12)/(2 w^2) = 6.0741, formed as (3w/2 + 6)/w^2: 3w + 12 would pass the ceiling */
        const int64_t L3 = div60(ONE60, sqrt60(3 * ONE60 + w / 4));       /* L3 = (3 + w/4)^(-1/2) = 0.5521 */
        const int64_t lnL3 = ln60(L3), lnw = ln60(w), lnS0 = -lnw / 2; LNS0 = lnS0; LNL3 = lnL3;   /* V = ln S0 = -ln(w)/2 */
        K_C1 = mul60(Q2v / 3, w2); K_C2 = w; K_C3 = mul60(mul60(Q2v, 2 * (ONE60 / 3)), w2); K_C4 = mul60(mul60(Q2v, 5 * (ONE60 / 6)), w2);   /* (Q2/3) w^2, kappa w, (2/3) Q2 w^2, (5/6) Q2 w^2 (2 Q2 and 5 Q2 would pass the ceiling) */
        for (int j = 0; j <= N + 1; j++) {
            const int64_t tj56 = TROW56(j);                       /* e^t as integer + fraction: z = e^t (up to 14, and past 16 on a coarse grid's ghost row): m halvings until t - m ln 2 < 2, then z = 2^m e^{t - m ln 2} */
            s128 e = {0, 0}; s128_addp(&e, 16, tj56); int64_t zi, zf; int m = 0;
            while (e.hi > 0 || (e.hi == 0 && e.lo >= (uint64_t)(2 * ONE60))) { s128_addp(&e, -1, ln2_60()); m++; }
            { const int64_t z = exp_of128(e); zi = z >> (60 - m); zf = (int64_t)(((uint64_t)z << m) & ((((uint64_t)1) << 60) - 1)); }   /* the low 60 bits survive the shift whatever the size */
            /* e^{-2z}, e^{-4z} with z = zi + zf: exponents in 128 bits */
            s128 m2 = {0, 0}; s128_addv(&m2, -2 * zf); s128_add(&m2, 0, 0); { s128 zz = { (-2 * zi) >> 4, (uint64_t)(-2 * zi) << 60 }; s128_add(&m2, zz.hi, zz.lo); }
            const int64_t em2z = exp_of128(m2);
            s128 m4 = {0, 0}; s128_addv(&m4, -4 * zf); { s128 zz = { (-4 * zi) >> 4, (uint64_t)(-4 * zi) << 60 }; s128_add(&m4, zz.hi, zz.lo); }
            const int64_t em4z = exp_of128(m4);
            CZ[j] = -1;                                                       /* coth(z)/z = (1 + e^-2z)/((1 - e^-2z)) e^-t, kept where z >= 0.4 (below 4.8; near the cut it is 1/z^2) */
            if (zi >= 1 || zf >= 2 * (ONE60 / 5)) { s128 mt = {0, 0}; s128_addp(&mt, -16, tj56); CZ[j] = mul60(div60(ONE60 + em2z, ONE60 - em2z), exp_of128(mt)); }
            const int64_t ln1p = ln60(ONE60 + em2z);                                  /* ln(1 + e^{-2z}), in [0, ln 2] */
            const s128 ln1m = ln1mexp(2, zi, zf, em2z, tj56), ln1m4 = ln1mexp(4, zi, zf, em4z, tj56);   /* ln(1 - e^{-2z}), ln(1 - e^{-4z}) in 128 bits (about t + ln 2 near the cut, past the Q.60 ceiling) */
            /* ln coth z = ln1p - ln1m;  ln sinh z = z + ln1m - ln 2;  Ur = ln L3 + ln coth z;  Br = ln L3 - ln sinh z (all in 128 bits, stored in Q.56) */
            s128 lncoth = {0, 0}; s128_addv(&lncoth, ln1p); { const s128 ng = s128_neg(ln1m); s128_add(&lncoth, ng.hi, ng.lo); }
            s128 lnsinh = { zi >> 4, (uint64_t)zi << 60 }; s128_addv(&lnsinh, zf - ln2_60()); s128_add(&lnsinh, ln1m.hi, ln1m.lo);   /* z + ln1m - ln2 in Q.60, 128 bits */
            { s128 ur = lncoth; s128_addv(&ur, lnL3); UR56[j] = s128_round_sh(ur, 4); }
            { s128 br = { (lnL3) < 0 ? -1 : 0, (uint64_t)lnL3 }; s128 neg = lnsinh; neg.lo = ~neg.lo + 1; neg.hi = ~neg.hi + (neg.lo == 0); s128_add(&br, neg.hi, neg.lo); BR56[j] = s128_round_sh(br, 4); }
            if (j <= N) {
                /* Ur' = -2 e^t / sinh(2 e^t) = -exp(t + ln 2 - ln sinh 2z),  ln sinh 2z = 2z + ln(1 - e^{-4z}) - ln 2 */
                s128 ex = {0, 0}; s128_addp(&ex, 16, tj56); s128_addv(&ex, 2 * ln2_60() - 2 * zf); { const s128 ng = s128_neg(ln1m4); s128_add(&ex, ng.hi, ng.lo); } { s128 zz = { (-2 * zi) >> 4, (uint64_t)(-2 * zi) << 60 }; s128_add(&ex, zz.hi, zz.lo); }
                URP[j] = -exp_of128(ex);
                /* e^{2(Br + t)} = exp(2 ln L3 + 2t - 2 ln sinh z) */
                s128 eb = {0, 0}; s128_addp(&eb, 2, lnL3); s128_addp(&eb, 32, tj56); {   /* 2 ln L3 + 2t in 128 bits: it is -15.19 on the cut row, past the 64-bit ceiling */ s128 two = lnsinh; two.hi = two.hi * 2 + (int64_t)(two.lo >> 63); two.lo <<= 1; two.lo = ~two.lo + 1; two.hi = ~two.hi + (two.lo == 0); s128_add(&eb, two.hi, two.lo); }
                WE[j] = mul60(WD[j], exp_of128(eb));                 /* w_E = w_D e^{2(Br + t)} */
            }
        }
        for (int j = 1; j <= N; j++) {                                  /* g_j for U and B from the Q.56 node differences; the deep end's ghost row mirrors (Neumann 0 on the reference's own derivative is NOT assumed: use the true node past the end) */
            s128 gu = {0, 0}; s128_addp(&gu, WS[j], UR56[j - 1] - UR56[j]); s128_addp(&gu, WN[j], UR56[j + 1] - UR56[j]); GU[j] = s128_round_sh(gu, 56);
            s128 gb = {0, 0}; s128_addp(&gb, WS[j], BR56[j - 1] - BR56[j]); s128_addp(&gb, WN[j], BR56[j + 1] - BR56[j]); GB[j] = s128_round_sh(gb, 56);
        }
        /* the deviation fields start at zero (Dirichlet zero at the cut); the x-side and deep-end Neumann values are zero */
        memset(F, 0, sizeof F); memset(EX, 0, sizeof EX); memset(GT1, 0, sizeof GT1); memset(GX0, 0, sizeof GX0); memset(GX1, 0, sizeof GX1); memset(DEF, 0, sizeof DEF);
        for (int j = 1; j <= N; j++) { int64_t ws[NF]; ads3_wsource(0, j, ws); DEF[0][j] = GU[j] - ws[0]; DEF[1][j] = -ws[1]; DEF[2][j] = GB[j] - ws[2]; }   /* the manufactured defects: at zero deviation the correction is then zero */
        if (MODEL == 2) {                                            /* the real P0: the brane's data at the cut as deviations from the reference, the same on every column, reached by a ramp in a */
            P0_W = w; P0_Q2V = Q2v; P0_FU4 = fu4; P0_FB4 = fb4; P0_FV4 = fv4; P0_SCALE = scale;
            const int nst = ramp > 1 ? ramp : 1; const int64_t a1 = ain >> (nst - 1);     /* stage k has the radius a 2^{k - nst}: the V deviation at the cut grows by ln 2 per stage */
            int64_t u0, v0, b0, logt; if (!p0_cutdata(a1, nst == 1, &u0, &v0, &b0, &logt)) { printf("WHRELAXT p0: the V deviation at the cut passes the Q.60 ceiling -- stop\n"); return 1; }
            for (int i = 0; i <= GX; i++) { F[0][i][0] = u0; F[1][i][0] = v0; F[2][i][0] = b0; }
            for (int i = 0; i <= GX; i++) for (int j = 1; j <= N; j++) { F[0][i][j] = mul60(u0, ONE60 - (ONE60 / N) * j); F[2][i][j] = mul60(b0, ONE60 - (ONE60 / N) * j); F[1][i][j] = p0_profile(a1, j); }   /* the start */
            char p1[48], p2[48], p3[48], p4[48], p5[48]; fmt60(ain, p1); fmt60(a1, p5); fmt60(u0, p2); fmt60(v0, p3); fmt60(b0, p4);
            printf("WHRELAXT p0: the brane's data at the cut (a = %s; U0 = B0 = -t0 + v^2/(18 a^2) + (6 Q2 - 1)/(60 a^4) v^4 ln v, V0 = ln a - t0 - 7 v^2/(36 a^2) - 3/2 of the log term, plus the free v^4 terms fu4, fb4, fv4: %s), reached in %d stages of the radius doubling from %s (the small terms only at the full radius); stage 1: deviations u0 = %s, v0 = %s, b0 = %s; the start: u and b linear from the cut values to zero at the deep end, v = (1/2) ln(1 + a^2 w e^-2t); later stages start from the previous solution plus the profile's increment\n",
                   p1, (fu4 || fb4 || fv4) ? "given" : "zero", nst, p5, p2, p3, p4);
        }
        if (pert && MODEL == 1) for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 1; j <= N; j++) {
            const int64_t s1 = ONE60 / N * j, bump = mul60(s1, ONE60 - s1);     /* a smooth bump (j/N)(1 - j/N), zero at the cut */
            F[k][i][j] = mul60(pert, bump) * (1 + (i & 1)) / 2 * (k + 1) / 2;
        }
        char q1[48], q2[48], q3[48], q4[48], qw[48], q4b[80]; fmt60(Q2v, q1); fmt60(L3, q2); fmt60(lnS0, q3); fmt60(pert, q4); fmt60(w, qw); sprintf(q4b, "a bump of amplitude %s", q4);
        printf("WHRELAXT %s: the AdS3 x S2 reference with w = %s: Q2 = %s, L3 = %s, V = ln S0 = %s; deviations (u, v, b) from (ln(L3 coth e^t), ln S0, ln(L3/sinh e^t)); the discretization defect %s; start = %s\n",
               MODEL == 2 ? "p0" : "ads3", qw, q1, q2, q3, DEFECT ? "subtracted so zero is the exact discrete solution" : "left in", MODEL == 2 ? "the brane's data at the cut (see the p0 line)" : q4b);
    } else
    /* the exact test fields, their Neumann values and their sources */
    for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 0; j <= N; j++) {
        const int64_t tj56 = TROW56(j), sv = 2 * tj56, xi = DX * i, P = Pk(k, sv), Q = Qx(xi);   /* sv = t/8 in Q.60 */
        EX[k][i][j] = mul60(P, Q);
        F[k][i][j] = (j == 0) ? EX[k][i][j] : 0;                        /* Dirichlet at the cut; the interior and the other sides start at 0 */
        /* S = P'' Q - P' Q + e^{2t} P Q'' : the operator on the exact field, evaluated with the same roundings */
        const int64_t ppq = mul60(Pkpp(k) - Pkp(k, sv), Q);
        s128 ex = {0, 0}; s128_addp(&ex, 32, tj56); s128_addp(&ex, -8, ln2_60()); int64_t n, f; split60(ex, &n, &f);
        const int64_t e2over256 = exp60n(n, f);                          /* e^{2t}/256 (below 0.77): e^{2t} itself would pass the ceiling */
        s128 term = {0, 0}; s128_addp(&term, e2over256, mul60(P, Qxpp()));   /* e^{2t} P Q'' / 256 in Q.120 */
        SRC[k][i][j] = ppq + s128_round_sh(term, F60 - 8);              /* times 256, rounded once */
        if (LAM) { SRC[k][i][j] -= mul60(LAM, mul60(EX[k][i][j], EX[k][i][j])); }
        if (j == N) GT1[k][i] = mul60(Pkp(k, sv), Q);                   /* d_t F at the deep end */
        if (i == 0) GX0[k][j] = mul60(P, Qxp(0)); if (i == GX) GX1[k][j] = mul60(P, Qxp(xi));   /* d_x F on the x sides */
    }
    if (LAM) for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 0; j <= N; j++) for (int m = 0; m < NF; m++) if (m != k) SRC[k][i][j] -= mul60(LAM, EX[m][i][j]);
    int64_t omega = N > 90 ? 2 * ONE60 - (ONE60 / N) * 90 : ONE60;  /* 2 - 90/N (N a power of two: exact); 1 below N = 91 */
    if (omega_in) omega = omega_in;
    OMEGA = omega;
    char lamb[48]; fmt60(lam, lamb); fmt56(T056, a); fmt60(t1, b); fmt60(DT, c); fmt60(CJ[N], d); fmt60(omega, e);
    printf("WHRELAXT: Q.60 fixed point; %d fields on rows uniform in t = ln(-y) from t0 = %s (Dirichlet, the cut) to t1 = ln 14 = %s (Neumann with values), N = %d (Dt = %s), Gx = %d columns of spacing 1 (Neumann with values); the operator F_tt - F_t + e^{2t} F_xx normalized by its row diagonal (c_j at the deep end = %s); red-black, all fields from one state, the plain diagonal as per-row weights, omega = %s%s, the change guard %s; lam = %s\n",
           NF, a, b, N, c, GX, d, e, omega_in ? " (given)" : " (2 - 90/N)", guard ? "on" : "off", lamb);
    const int64_t worst = whq60_selftest("WHRELAXT"); printf("WHRELAXT selftest: largest deviation %lld units of 2^-60\n", (long long)worst);
    if (start_exact) memcpy(F, EX, sizeof F);
    if (dbg) {                                                       /* the correction at the exact solution, per row of column 1, and a few values */
        for (int j = 1; j <= N; j += (N > 16 ? N / 16 : 1)) { char q1[48], q2[48]; const int64_t c0 = correction(0, 1, j); fmt60(EX[0][1][j], q1); fmt60(F[0][1][j], q2);
            printf("WHRELAXT dbg row %d: exact F_1 = %s, F_1 = %s, correction at this state = %lld units, w_S %lld w_N %lld w_X %lld w_D %lld c_j %lld (units)\n", j, q1, q2, (long long)c0, (long long)WS[j], (long long)WN[j], (long long)WX[j], (long long)WD[j], (long long)CJ[j]); }
    }
    uint64_t chain = 0xcbf29ce484222325ull;
    int sweep; int64_t maxres = -1, maxchg = 0, maxcorr = 0, smallest = -1, ckp_chg = -1; int restores = 0, last_ckp = 0;
    const int mid = (GX + 1) / 2;                                      /* the columns 0..GX in two halves (the x sides are unknowns too) */
    const int nstages = (MODEL == 2 && ramp > 1) ? ramp : 1; int total_sweeps = 0, stopped = 0;
    for (int stage = 1; stage <= nstages; stage++) {
    if (MODEL == 2 && stage > 1) {                                   /* the next radius: new data on the cut row, the previous solution plus the profile's increment as the start */
        const int64_t ak = ain >> (nstages - stage), akp = ain >> (nstages - stage + 1);
        int64_t u0, v0, b0, lt; if (!p0_cutdata(ak, stage == nstages, &u0, &v0, &b0, &lt)) { printf("WHRELAXT p0: the V deviation at the cut passes the Q.60 ceiling -- stop\n"); return 1; }
        for (int i = 0; i <= GX; i++) { F[0][i][0] = u0; F[1][i][0] = v0; F[2][i][0] = b0; for (int j = 1; j <= N; j++) F[1][i][j] += p0_profile(ak, j) - p0_profile(akp, j); }
        OMEGA = omega; smallest = -1; ckp_chg = -1; last_ckp = 0; maxres = -1; stopped = 0;
        char s1[48], s2[48], s3[48]; fmt60(ak, s1); fmt60(v0, s2); fmt60(lt, s3);
        printf("WHRELAXT p0 stage %d of %d: a = %s, the V deviation at the cut %s%s%s\n", stage, nstages, s1, s2, stage == nstages ? "; the small FG terms in, the log term " : "", stage == nstages ? s3 : "");
    }
    memcpy(CKP, F, sizeof F);
    for (sweep = 1; sweep <= 200000; sweep++) {
        maxchg = 0; maxcorr = 0;
        for (int color = 0; color < 2; color++) {
            rjob x = { color, 0, mid, 0, 0, 0 }, y = { color, mid, GX + 1, 0, 0, 0 };
            par2(sweep_rows, &x, &y);
            if (x.maxchg > maxchg) maxchg = x.maxchg; if (y.maxchg > maxchg) maxchg = y.maxchg;
            if (x.maxcorr > maxcorr) maxcorr = x.maxcorr; if (y.maxcorr > maxcorr) maxcorr = y.maxcorr;
        }
        chain = fnv(chain, &maxcorr, 8);
        if (smallest < 0 || maxchg < smallest) smallest = maxchg;
        if (guard && smallest > 0 && maxchg / 32 > smallest && sweep - last_ckp > 1) {   /* the change guard (compared by division: 32 x could pass the int64) */
            memcpy(F, CKP, sizeof F); OMEGA -= ONE60 / 10; if (OMEGA < ONE60) OMEGA = ONE60; restores++;
            fmt60(OMEGA, a); printf("WHRELAXT guard: the change %lld units exceeded 32 x the smallest %lld at sweep %d: checkpoint restored (sweep %d), omega -> %s\n", (long long)maxchg, (long long)smallest, sweep, last_ckp, a);
            smallest = -1; continue;
        }
        if (sweep % 64 == 0 && (ckp_chg < 0 || maxchg <= ckp_chg)) { memcpy(CKP, F, sizeof F); ckp_chg = maxchg; last_ckp = sweep; }
        if (sweep % 16 == 0 || maxcorr <= tol) {
            rjob x = { 0, 0, mid, 0, 0, 0 }, y = { 0, mid, GX + 1, 0, 0, 0 };
            par2(resid_rows, &x, &y); maxres = x.maxres > y.maxres ? x.maxres : y.maxres;
        }
        if ((sweep & (sweep - 1)) == 0 || sweep % 1000 == 0) printf("WHRELAXT sweep %d: max correction = %lld units of 2^-60, max change = %lld, max |residual| = %lld\n", sweep, (long long)maxcorr, (long long)maxchg, (long long)maxres);
        if (maxcorr <= tol && maxres >= 0 && maxres <= tol) { printf("WHRELAXT: stop at sweep %d: the largest correction %lld and residual %lld are at or below tol %lld units\n", sweep, (long long)maxcorr, (long long)maxres, (long long)tol); stopped = 1; break; }
        if (maxchg == 0) { stopped = 1; break; }
    }
    if (sweep > 200000) sweep = 200000; total_sweeps += sweep;
    if (nstages > 1) { fmt60(OMEGA, a); printf("WHRELAXT p0 stage %d of %d: %s after %d sweeps (omega at the end %s, %d guard restores so far); %d sweeps in all\n", stage, nstages, stopped ? "converged" : "NOT converged", sweep, a, restores, total_sweeps); }
    }
    if (nstages > 1) sweep = total_sweeps;
    int64_t maxerr[NF] = {0, 0, 0};
    for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 0; j <= N; j++) { int64_t er = F[k][i][j] - EX[k][i][j]; if (er < 0) er = -er; if (er > maxerr[k]) maxerr[k] = er; }
    fmt60(F[0][GX / 2][N / 2], a); fmt60(EX[0][GX / 2][N / 2], b); fmt60(OMEGA, c);
    printf("WHRELAXT: after %d sweeps (omega at the end %s, %d guard restores): max |F - exact| = %lld, %lld, %lld units of 2^-60 for the three fields; F_1 at the middle = %s (exact %s)\n",
           sweep, c, restores, (long long)maxerr[0], (long long)maxerr[1], (long long)maxerr[2], a, b);
    if (MODEL >= 1) {                                                /* the continuum reference's defect on this grid, as residuals: the interior rows and the deep (Neumann) row */
        int64_t din[NF] = {0, 0, 0}, dnn[NF];
        for (int k = 0; k < NF; k++) { for (int j = 1; j < N; j++) { const int64_t rr = residual_units(DEF[k][j], j); if (rr > din[k]) din[k] = rr; } dnn[k] = residual_units(DEF[k][N], N); }
        printf("WHRELAXT ads3: the reference's defect on this grid as residuals: interior rows max = %lld, %lld, %lld units of 2^-60, the deep Neumann row = %lld, %lld, %lld; the defect was %s\n",
               (long long)din[0], (long long)din[1], (long long)din[2], (long long)dnn[0], (long long)dnn[1], (long long)dnn[2], DEFECT ? "subtracted (zero deviation is the exact discrete solution)" : "left in (the end deviation is the discretization error)");
    }
    if (MODEL == 2) {                                                /* the real P0's answer: the full fields on five rows (U, B in Q.56; V in Q.60) and the column spread (zero for x-independent data) */
        int64_t spread = 0; for (int k = 0; k < NF; k++) for (int i = 1; i <= GX; i++) for (int j = 0; j <= N; j++) { int64_t dd = F[k][i][j] - F[k][0][j]; if (dd < 0) dd = -dd; if (dd > spread) spread = dd; }
        for (int j = 0; j <= N; j += N / 4) {
            char r1[48], r2[48], r3[48], r4[48], r5[48], r6[48], r7[48];
            fmt56(TROW56(j), r1); fmt56(UR56[j] + ((F[0][0][j] + 8) >> 4), r2); fmt56(BR56[j] + ((F[2][0][j] + 8) >> 4), r3); fmt60(LNS0 + F[1][0][j], r4); fmt60(F[0][0][j], r5); fmt60(F[1][0][j], r6); fmt60(F[2][0][j], r7);
            char r8[48] = "(near the cut: ill-conditioned, not printed)";
            if (CZ[j] >= 0) { const int64_t qr = mul60(exp60(F[0][0][j] + 2 * F[1][0][j]), ONE60 + mul60(CZ[j], dT(0, 0, j) - dT(2, 0, j))) - ONE60; fmt60(qr, r8); }
            printf("WHRELAXT p0 row %d (t = %s): U = %s, B = %s, V = %s; deviations u = %s, v = %s, b = %s; the first integral Q/(-L3 S0^2) - 1 = %s\n", j, r1, r2, r3, r4, r5, r6, r7, r8);
        }
        { char c1[48], c2[48], c3[48]; fmt56(UR56[N] + ((F[0][0][N] + 8) >> 4), c1); fmt60(LNS0 + F[1][0][N], c2); fmt60(LNL3 + F[2][0][N], c3);   /* B + ln sinh(-y) = (Br + b) + ln sinh z = ln L3 + b */
          const int64_t c2r = exp60(F[0][0][N] + 2 * F[1][0][N]); char c4[48]; fmt60(c2r - ONE60, c4);   /* c2_eff/c2 - 1 = e^{u_c + 2 v_c} - 1 */
          printf("WHRELAXT p0 centre (the deep end, y = -14): U_c = %s, V_c = %s, B + ln sinh(-y) = %s; the central charge against the AdS3 x S2 fixed point's: c2_eff/c2 - 1 = e^{u_c + 2 v_c} - 1 = %s (c2_eff = 3 L_c e^{U_c + 2 V_c}/G5 with L_c = 2 pi, c2 = 6 pi L3 S*^2/G5; the reviewer's law 1 + 5.19 (l_B/L_x)^2 with l_B = lhat a, L_x = 2 pi)\n", c1, c2, c3, c4); }
        printf("WHRELAXT p0: the column spread (x-independent data must give identical columns) = %lld units of 2^-60\n", (long long)spread);
    }
    printf("WHRELAXT: chain %016llx\n", (unsigned long long)chain);
    printf("WHRELAXT: done\n");
    return 0;
}
