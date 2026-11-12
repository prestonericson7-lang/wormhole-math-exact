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
 * Usage: whrelaxt [model=quad|ads3] [N=128] [Gx=4] [lam=0] [pert=0.01] [defect=1] [omega=auto|value] [tol=1048576] [t0=-7.0] [guard=1]    N 4..1024 (ads3: 64 or more; defect=0 leaves the reference's discretization defect in); Gx up to 64. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "whq60.h"
#include "par2.h"

#define MAXN 1024
#define MAXGX 64
#define NF 3
static int N, GX; static int64_t DT, DX, T0, OMEGA, LAM;
static int64_t F[NF][MAXGX + 1][MAXN + 1], SRC[NF][MAXGX + 1][MAXN + 1], EX[NF][MAXGX + 1][MAXN + 1], CKP[NF][MAXGX + 1][MAXN + 1];
static int64_t CJ[MAXN + 1], WS[MAXN + 1], WN[MAXN + 1], WX[MAXN + 1], WD[MAXN + 1], DIAG1[MAXN + 1];   /* per row: c_j, the weights, and 1 + c_j */
static inline void split60(s128 x, int64_t *n, int64_t *f)         /* a 128-bit Q.60 exponent into integer + fraction for exp60n */
{ *n = (int64_t)(((uint64_t)x.hi << 4) | (x.lo >> 60)); *f = (int64_t)(x.lo & (((uint64_t)1 << 60) - 1)); }
#define TROW(j) ((int64_t)((uint64_t)T0 + (uint64_t)DT * (uint64_t)(j)))   /* t_j = t0 + Dt j: Dt j alone reaches 9.6, past the Q.60 ceiling, so the sum is formed modulo 2^64 (defined) and lands back below 8 */
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
static inline int64_t Pk(int k, int64_t t) { const int64_t u = ONE60 / 16, sv = t / 8; return PC[k][0] * u + mul60(PC[k][1] * u, sv) + mul60(PC[k][2] * u, mul60(sv, sv)); }
static inline int64_t Pkp(int k, int64_t t) { const int64_t u = ONE60 / 16, sv = t / 8; return (PC[k][1] * u + 2 * mul60(PC[k][2] * u, sv)) / 8; }   /* dP/dt = (p1 + 2 p2 s)/(16 x 8) */
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
/* ads3: the three weighted sources at a point: w_D S_a + w_X S_b + w_E e^{2b} S_c (the manufactured defect is subtracted in correction()) */
static void ads3_wsource(int i, int j, int64_t out[NF])
{
    const int64_t ut = URP[j] + dT(0, i, j), vt = dT(1, i, j), ux = dX(0, i, j), vx = dX(1, i, j), b = F[2][i][j], v = F[1][i][j];
    const int64_t e2b = exp60(2 * b), em2v = exp60(-2 * v), em4v = exp60(-4 * v);
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
int main(int argc, char **argv)
{
    N = 128; GX = 4; int guard = 1, start_exact = 0, dbg = 0; int64_t lam = 0, omega_in = 0, tol = (int64_t)1 << 20, t0 = -7 * ONE60, pert = 0;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "Gx=", 3)) GX = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "lam=", 4)) lam = parse60(argv[i] + 4);
        else if (!strncmp(argv[i], "tol=", 4)) tol = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "t0=", 3)) t0 = parse60(argv[i] + 3);
        else if (!strncmp(argv[i], "guard=", 6)) guard = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "start=", 6)) start_exact = !strcmp(argv[i] + 6, "exact");
        else if (!strncmp(argv[i], "dbg=", 4)) dbg = atoi(argv[i] + 4);
        else if (!strncmp(argv[i], "model=", 6)) MODEL = !strcmp(argv[i] + 6, "ads3");
        else if (!strncmp(argv[i], "pert=", 5)) pert = parse60(argv[i] + 5);
        else if (!strncmp(argv[i], "defect=", 7)) DEFECT = atoi(argv[i] + 7);
        else if (!strncmp(argv[i], "omega=", 6)) { if (strcmp(argv[i] + 6, "auto")) omega_in = parse60(argv[i] + 6); }
    }
    if (N < 4 || N > MAXN || (N & (N - 1))) { printf("WHRELAXT: N must be a power of two in 4..%d\n", MAXN); return 1; }
    if (GX < 2 || GX > MAXGX) { printf("WHRELAXT: Gx must be 2..%d\n", MAXGX); return 1; }
    char a[48], b[48], c[48], d[48], e[48];
    const int64_t t1 = ln60(14 * (ONE60 / 16)) + 4 * ln2_60();          /* ln 14 = ln(14/16) + 4 ln 2 (14 itself is past the ceiling) */
    T0 = t0; DT = rdiv64(t1, N) - rdiv64(t0, N); DX = ONE60; LAM = lam;   /* t1 - t0 (9.6) is past the Q.60 ceiling: divide the ends separately */
    AS = ONE60 + DT / 2; AN = ONE60 - DT / 2; DT2HALF = mul60(DT, DT) / 2;
    for (int j = 0; j <= N; j++) {                                    /* per row: t_j, e^{2 t_j}, c_j, the reciprocal */
        const int64_t tj = TROW(j);
        /* c_j = e^{2t} Dt^2 / Dx^2 = e^{2t + 2 ln Dt} (Dx = 1): one exponential, below 2 at N >= 128 (e^{2t} itself reaches 196, past the ceiling) */
        s128 ex2 = {0, 0}; s128_addp(&ex2, 2, tj); s128_addp(&ex2, 2, ln60(DT)); int64_t n2, f2; split60(ex2, &n2, &f2);
        if (ex2.hi < 0) {                                                 /* c_j = e^g below 1: the direct form; the diagonal 2 (1 + c_j) */
            CJ[j] = exp60n(n2, f2); DIAG1[j] = ONE60 + CJ[j];
            WS[j] = div60(AS / 2, DIAG1[j]); WN[j] = div60(AN / 2, DIAG1[j]); WX[j] = div60(CJ[j] / 2, DIAG1[j]); WD[j] = div60(mul60(DT, DT) / 2, DIAG1[j]);
        } else {                                                          /* c_j at or above 1 (coarse grids): through 1/c = e^-g, which is below 1 */
            s128 neg = {0, 0}; s128_addp(&neg, -2, tj); s128_addp(&neg, -2, ln60(DT)); int64_t n3, f3; split60(neg, &n3, &f3);
            const int64_t rc = exp60n(n3, f3), den = ONE60 + rc;           /* 1/c and 1 + 1/c */
            CJ[j] = rc > ONE60 / 8 ? div60(ONE60, rc) : -1; DIAG1[j] = 0;   /* c itself only where it fits (the banner shows -1 past the ceiling); the weights never need it */
            WS[j] = div60(mul60(AS / 2, rc), den); WN[j] = div60(mul60(AN / 2, rc), den); WX[j] = div60(ONE60 / 2, den); WD[j] = div60(mul60(mul60(DT, DT) / 2, rc), den);
        }
    }
    { uint64_t rem; INV2DT56 = (int64_t)div128((uint64_t)1 << (56 + 60 - 64), 0, (uint64_t)(2 * DT), &rem); if (2 * rem >= (uint64_t)(2 * DT)) INV2DT56++; }   /* 2^116 / (2 Dt 2^60) = 2^56/(2 Dt) */
    if (MODEL == 1) {                                                /* the AdS3 x S2 reference and the manufactured defects */
        if (N < 64) { printf("WHRELAXT: model=ads3 needs N >= 64\n"); return 1; }
        const int64_t w = parse60("1.125"), w2 = mul60(w, w);         /* w = 1/S0^2 = 9/8 */
        const int64_t Q2v = div60(3 * w / 2 + 6 * ONE60, w2);              /* Q2 = (3w + 12)/(2 w^2) = 6.0741, formed as (3w/2 + 6)/w^2: 3w + 12 would pass the ceiling */
        const int64_t L3 = div60(ONE60, sqrt60(3 * ONE60 + w / 4));       /* L3 = (3 + w/4)^(-1/2) = 0.5521 */
        const int64_t lnL3 = ln60(L3), lnw = ln60(w), lnS0 = -lnw / 2;   /* V = ln S0 = -ln(w)/2 */
        K_C1 = mul60(Q2v / 3, w2); K_C2 = w; K_C3 = mul60(mul60(Q2v, 2 * (ONE60 / 3)), w2); K_C4 = mul60(mul60(Q2v, 5 * (ONE60 / 6)), w2);   /* (Q2/3) w^2, kappa w, (2/3) Q2 w^2, (5/6) Q2 w^2 (2 Q2 and 5 Q2 would pass the ceiling) */
        static int64_t UR56[MAXN + 2], BR56[MAXN + 2];                 /* the reference's node values in Q.56 (U near 7, B near -13), one ghost row past the deep end */
        for (int j = 0; j <= N + 1; j++) {
            const int64_t tj = TROW(j);                           /* e^t as integer + fraction: z = e^t (up to 14, and past 16 on a coarse grid's ghost row): m halvings until t - m ln 2 < 2, then z = 2^m e^{t - m ln 2} */
            s128 e = {0, 0}; s128_addv(&e, tj); int64_t zi, zf; int m = 0;
            while (tj - m * ln2_60() >= 2 * ONE60) m++;
            s128_addp(&e, -m, ln2_60());
            { const int64_t z = exp_of128(e); zi = z >> (60 - m); zf = (int64_t)(((uint64_t)z << m) & ((((uint64_t)1) << 60) - 1)); }   /* the low 60 bits survive the shift whatever the size */
            /* e^{-2z}, e^{-4z} with z = zi + zf: exponents in 128 bits */
            s128 m2 = {0, 0}; s128_addv(&m2, -2 * zf); s128_add(&m2, 0, 0); { s128 zz = { (-2 * zi) >> 4, (uint64_t)(-2 * zi) << 60 }; s128_add(&m2, zz.hi, zz.lo); }
            const int64_t em2z = exp_of128(m2);
            s128 m4 = {0, 0}; s128_addv(&m4, -4 * zf); { s128 zz = { (-4 * zi) >> 4, (uint64_t)(-4 * zi) << 60 }; s128_add(&m4, zz.hi, zz.lo); }
            const int64_t em4z = exp_of128(m4);
            const int64_t ln1m = ln60(ONE60 - em2z), ln1p = ln60(ONE60 + em2z);   /* ln(1 - e^{-2z}), ln(1 + e^{-2z}) */
            /* ln coth z = ln1p - ln1m;  ln sinh z = z + ln1m - ln 2;  Ur = ln L3 + ln coth z;  Br = ln L3 - ln sinh z (Q.56) */
            const int64_t lncoth = ln1p - ln1m;
            s128 lnsinh = { zi >> 4, (uint64_t)zi << 60 }; s128_addv(&lnsinh, zf + ln1m - ln2_60());                       /* z + ln1m - ln2 in Q.60, 128 bits */
            UR56[j] = s128_round_sh((s128){ (lnL3 + lncoth) < 0 ? -1 : 0, (uint64_t)(lnL3 + lncoth) }, 4);
            { s128 br = { (lnL3) < 0 ? -1 : 0, (uint64_t)lnL3 }; s128 neg = lnsinh; neg.lo = ~neg.lo + 1; neg.hi = ~neg.hi + (neg.lo == 0); s128_add(&br, neg.hi, neg.lo); BR56[j] = s128_round_sh(br, 4); }
            if (j <= N) {
                /* Ur' = -2 e^t / sinh(2 e^t) = -exp(t + ln 2 - ln sinh 2z),  ln sinh 2z = 2z + ln(1 - e^{-4z}) - ln 2 */
                s128 ex = {0, 0}; s128_addv(&ex, tj + 2 * ln2_60() - 2 * zf - ln60(ONE60 - em4z)); { s128 zz = { (-2 * zi) >> 4, (uint64_t)(-2 * zi) << 60 }; s128_add(&ex, zz.hi, zz.lo); }
                URP[j] = -exp_of128(ex);
                /* e^{2(Br + t)} = exp(2 ln L3 + 2t - 2 ln sinh z) */
                s128 eb = {0, 0}; s128_addp(&eb, 2, lnL3); s128_addp(&eb, 2, tj); {   /* 2 ln L3 + 2t in 128 bits: it is -15.19 on the cut row, past the 64-bit ceiling */ s128 two = lnsinh; two.hi = two.hi * 2 + (int64_t)(two.lo >> 63); two.lo <<= 1; two.lo = ~two.lo + 1; two.hi = ~two.hi + (two.lo == 0); s128_add(&eb, two.hi, two.lo); }
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
        if (pert) for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 1; j <= N; j++) {
            const int64_t s1 = ONE60 / N * j, bump = mul60(s1, ONE60 - s1);     /* a smooth bump (j/N)(1 - j/N), zero at the cut */
            F[k][i][j] = mul60(pert, bump) * (1 + (i & 1)) / 2 * (k + 1) / 2;
        }
        char q1[48], q2[48], q3[48], q4[48]; fmt60(Q2v, q1); fmt60(L3, q2); fmt60(lnS0, q3); fmt60(pert, q4);
        printf("WHRELAXT ads3: the AdS3 x S2 reference with w = 9/8: Q2 = %s, L3 = %s, V = ln S0 = %s; deviations (u, v, b) from (ln(L3 coth e^t), ln S0, ln(L3/sinh e^t)); the discretization defect subtracted so zero is the exact discrete solution; start = a bump of amplitude %s\n", q1, q2, q3, q4);
    } else
    /* the exact test fields, their Neumann values and their sources */
    for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 0; j <= N; j++) {
        const int64_t tj = TROW(j), xi = DX * i, P = Pk(k, tj), Q = Qx(xi);
        EX[k][i][j] = mul60(P, Q);
        F[k][i][j] = (j == 0) ? EX[k][i][j] : 0;                        /* Dirichlet at the cut; the interior and the other sides start at 0 */
        /* S = P'' Q - P' Q + e^{2t} P Q'' : the operator on the exact field, evaluated with the same roundings */
        const int64_t ppq = mul60(Pkpp(k) - Pkp(k, tj), Q);
        s128 ex = {0, 0}; s128_addp(&ex, 2, tj); s128_addp(&ex, -8, ln2_60()); int64_t n, f; split60(ex, &n, &f);
        const int64_t e2over256 = exp60n(n, f);                          /* e^{2t}/256 (below 0.77): e^{2t} itself would pass the ceiling */
        s128 term = {0, 0}; s128_addp(&term, e2over256, mul60(P, Qxpp()));   /* e^{2t} P Q'' / 256 in Q.120 */
        SRC[k][i][j] = ppq + s128_round_sh(term, F60 - 8);              /* times 256, rounded once */
        if (LAM) { SRC[k][i][j] -= mul60(LAM, mul60(EX[k][i][j], EX[k][i][j])); }
        if (j == N) GT1[k][i] = mul60(Pkp(k, tj), Q);                   /* d_t F at the deep end */
        if (i == 0) GX0[k][j] = mul60(P, Qxp(0)); if (i == GX) GX1[k][j] = mul60(P, Qxp(xi));   /* d_x F on the x sides */
    }
    if (LAM) for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 0; j <= N; j++) for (int m = 0; m < NF; m++) if (m != k) SRC[k][i][j] -= mul60(LAM, EX[m][i][j]);
    int64_t omega = N > 90 ? 2 * ONE60 - (ONE60 / N) * 90 : ONE60;  /* 2 - 90/N (N a power of two: exact); 1 below N = 91 */
    if (omega_in) omega = omega_in;
    OMEGA = omega;
    char lamb[48]; fmt60(lam, lamb); fmt60(T0, a); fmt60(t1, b); fmt60(DT, c); fmt60(CJ[N], d); fmt60(omega, e);
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
        if (maxcorr <= tol && maxres >= 0 && maxres <= tol) { printf("WHRELAXT: stop at sweep %d: the largest correction %lld and residual %lld are at or below tol %lld units\n", sweep, (long long)maxcorr, (long long)maxres, (long long)tol); break; }
        if (maxchg == 0) break;
    }
    int64_t maxerr[NF] = {0, 0, 0};
    for (int k = 0; k < NF; k++) for (int i = 0; i <= GX; i++) for (int j = 0; j <= N; j++) { int64_t er = F[k][i][j] - EX[k][i][j]; if (er < 0) er = -er; if (er > maxerr[k]) maxerr[k] = er; }
    fmt60(F[0][GX / 2][N / 2], a); fmt60(EX[0][GX / 2][N / 2], b); fmt60(OMEGA, c);
    printf("WHRELAXT: after %d sweeps (omega at the end %s, %d guard restores): max |F - exact| = %lld, %lld, %lld units of 2^-60 for the three fields; F_1 at the middle = %s (exact %s)\n",
           sweep, c, restores, (long long)maxerr[0], (long long)maxerr[1], (long long)maxerr[2], a, b);
    if (MODEL == 1) {                                                /* the continuum reference's defect on this grid, as residuals: the interior rows and the deep (Neumann) row */
        int64_t din[NF] = {0, 0, 0}, dnn[NF];
        for (int k = 0; k < NF; k++) { for (int j = 1; j < N; j++) { const int64_t rr = residual_units(DEF[k][j], j); if (rr > din[k]) din[k] = rr; } dnn[k] = residual_units(DEF[k][N], N); }
        printf("WHRELAXT ads3: the reference's defect on this grid as residuals: interior rows max = %lld, %lld, %lld units of 2^-60, the deep Neumann row = %lld, %lld, %lld; the defect was %s\n",
               (long long)din[0], (long long)din[1], (long long)din[2], (long long)dnn[0], (long long)dnn[1], (long long)dnn[2], DEFECT ? "subtracted (zero deviation is the exact discrete solution)" : "left in (the end deviation is the discretization error)");
    }
    printf("WHRELAXT: chain %016llx\n", (unsigned long long)chain);
    printf("WHRELAXT: done\n");
    return 0;
}
