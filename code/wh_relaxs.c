/* wh_relaxs.c (11-12-2026, 2:00 AM) -- WORMHOLE MATH 1b, THE SYSTEM KERNEL: the relaxation of several coupled fields at once, each with a
   nonlinear source, with a Newton step per point over the fields' local Jacobian block, in exact fixed point; built and verified on a
   problem with a known answer before the 1b equations arrive. The single-field kernel (wh_relax.c, versions 1-4) proved the sweep, the
   residual accounting, the two cores, a Robin boundary and stretched grids. The 1b system has three unknown functions in the conformal
   gauge and five in the Einstein-DeTurck gauge, with sources that depend on all of them (the exponential of the brane's warp factor among
   them), so the kernel must relax a VECTOR of fields and take the Newton step over the whole local block, not field by field: where a
   source's own derivative is comparable to the Laplacian's diagonal, the reviewer's float tests showed, a point update that lags it
   overshoots and diverges on coarse grids.

   THE TEST PROBLEM: NF = 3 fields on the unit square, zero on the boundary,
       Laplacian u_k = f_k + lam u_k^2 + SUM_{m != k} u_m,       k = 1, 2, 3,
   with the exact solution u_k = c_k x(1-x) y(1-y), c = (1, 2, 3)/scale (scale = 4 by default), and f_k defined from it on the grid with the
   same products, so the discrete solution equals u_k at the nodes to the arithmetic (the five-point Laplacian is exact for this u; the
   products are rounded once). lam, an integer, sets the self-nonlinearity: lam = 1 is mild; THE STRESS CASE is scale = 40 with
   lam = 150000, where the source lam u^2 stays below the Q.60 ceiling (3.3) while its own derivative 2 lam h^2 u reaches a third of the
   Laplacian's diagonal 4 at G = 32 (the regime where the reviewer's float tests overshot with a lagged source); lam = 0 leaves the
   linear coupling only.

   THE KERNEL: red-black Newton-SOR over the vector of fields. At each point the residuals R_k = 4 u_k - S_k + h^2 (f_k + lam u_k^2 +
   SUM_{m != k} u_m) (S_k the four neighbours of field k) and the Jacobian J_km = dR_k/du_m = (4 + 2 lam h^2 u_k) delta_km + h^2 (1 - delta_km).
   jac=block (the default) solves J d = -R by Gaussian elimination in Q.60 (a few roundings; the solve only steers the iteration, the
   converged solution is fixed by R = 0, which is measured exactly) and updates all fields at once with the over-relaxation omega;
   jac=diag takes the Newton step field by field with the other fields lagged (the diagonal of J only); jac=none takes the plain
   Gauss-Seidel step 4 u_k = S_k - h^2 (...) with every source term lagged, the update a float code would write first. The residual
   (S_k - 4 u_k) G^2 - f_k - lam u_k^2 - SUM u_m is measured in 128 bits with the products rounded once, and its maximum over the fields and
   the points is the convergence measure and the chain; omega = 2 / (1 + sin(pi/G)) or omega= given. Two cores as wh_relax.c version 3:
   the rows in halves through par2, the same bits either way.

   Divergence is declared when a sweep's largest change exceeds 2 (the fields are below 1); the residual measure saturates at the int64
   through the transient without stopping anything.
   Usage: whrelaxs [G=64] [lam=1] [scale=4] [jac=block|diag|none] [omega=auto|value] [thr=1]      G a power of two, 4..512. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "whq.h"
#include "../par2.h"

#define MAXG 512
#define NF 3
#define F60 60
static const int64_t ONE60 = (int64_t)1 << F60;
typedef struct { int64_t hi; uint64_t lo; } s128;
static int64_t U[NF][MAXG + 1][MAXG + 1], Fv[NF][MAXG + 1][MAXG + 1], Ex[NF][MAXG + 1][MAXG + 1];

static inline void s128_add(s128 *a, int64_t hi, uint64_t lo) { uint64_t l = a->lo + lo; a->hi += hi + (l < lo); a->lo = l; }
static inline void mulprod(int64_t a, int64_t b, int64_t *hi, uint64_t *lo)
{
    int neg = (a < 0) ^ (b < 0); uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b, h, l;
    mul128(ua, ub, &h, &l); if (neg) { l = ~l + 1; h = ~h + (l == 0); } *hi = (int64_t)h; *lo = l;
}
static inline void s128_addp(s128 *a, int64_t x, int64_t y) { int64_t h; uint64_t l; mulprod(x, y, &h, &l); s128_add(a, h, l); }
static inline void s128_addv(s128 *a, int64_t v) { s128_add(a, v < 0 ? -1 : 0, (uint64_t)v); }
static inline int64_t s128_round_sh(s128 a, int sh)
{
    int neg = a.hi < 0; if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t r = sh ? ((a.lo >> sh) | ((uint64_t)a.hi << (64 - sh))) : a.lo;
    if (sh && (a.lo & ((((uint64_t)1) << sh) - 1)) >= (((uint64_t)1) << (sh - 1))) r++;
    if (r > (uint64_t)INT64_MAX) r = (uint64_t)INT64_MAX;                        /* the residual measure saturates */
    return neg ? -(int64_t)r : (int64_t)r;
}
static inline int64_t mul60(int64_t a, int64_t b) { s128 s = {0, 0}; s128_addp(&s, a, b); return s128_round_sh(s, F60); }
static inline int64_t rdiv64(int64_t x, int64_t n) { return x >= 0 ? (x + n / 2) / n : -((-x + n / 2) / n); }
static int64_t div60(int64_t a, int64_t b)             /* round(a 2^60 / b), b > 0, |a| < 8 b */
{
    int neg = a < 0; uint64_t ua = neg ? -(uint64_t)a : (uint64_t)a; const uint64_t hi = ua >> 4, lo = ua << 60;
    if (hi >= (uint64_t)b) { printf("WHRELAXS: div60 overflow -- stop\n"); exit(1); }
    uint64_t rem; uint64_t qq = div128(hi, lo, (uint64_t)b, &rem); if (2 * rem >= (uint64_t)b) qq++;
    return neg ? -(int64_t)qq : (int64_t)qq;
}
static int64_t parse60(const char *s)
{
    int neg = 0; if (*s == '-') { neg = 1; s++; }
    uint64_t ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (uint64_t)(*s++ - '0');
    uint64_t fp = 0, scale = 1; if (*s == '.') { s++; while (*s >= '0' && *s <= '9' && scale < 1000000000000ull) { fp = fp * 10 + (uint64_t)(*s++ - '0'); scale *= 10; } }
    if (ip >= 7) { printf("WHRELAXS: a value at or above 7 cannot be held in Q.60\n"); exit(1); }
    int64_t v = (int64_t)(ip << F60);
    if (fp) { uint64_t h, l, rem; mul128(fp, (uint64_t)1 << F60, &h, &l); uint64_t qq = div128(h, l, scale, &rem); if (2 * rem >= scale) qq++; v += (int64_t)qq; }
    return neg ? -v : v;
}
static const int64_t PI60 = 3622009729038561421ll;
static int64_t sin60_small(int64_t x)                  /* sin x for |x| < 1 by its Taylor series in Q.60 */
{
    const int64_t x2 = mul60(x, x); int64_t term = x, sum = x;
    for (int k = 1; k < 40; k++) { term = -rdiv64(mul60(term, x2), (int64_t)(2 * k) * (2 * k + 1)); if (!term) break; sum += term; }
    return sum;
}
static void fmt60(int64_t v, char *out)
{
    int neg = v < 0; uint64_t u = neg ? -(uint64_t)v : (uint64_t)v;
    uint64_t ip = u >> F60, frac = u & ((((uint64_t)1) << F60) - 1), h, l;
    mul128(frac, 1000000000000ull, &h, &l);
    uint64_t dg = (l >> F60) | (h << (64 - F60)); if ((l >> (F60 - 1)) & 1) dg++;
    if (dg >= 1000000000000ull) { dg -= 1000000000000ull; ip++; }
    sprintf(out, "%s%llu.%012llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)dg);
}

static int Gg, JAC; static int64_t G2g, OMEGA, LAM, H2;
/* the source of field k at a point from the current fields: f_k + lam u_k^2 + SUM_{m != k} u_m (Q.60, the square rounded once) */
static inline int64_t source(int k, int i, int j)
{
    int64_t s = Fv[k][i][j] + LAM * mul60(U[k][i][j], U[k][i][j]);
    for (int m = 0; m < NF; m++) if (m != k) s += U[m][i][j];
    return s;
}
/* J d = -R for the 3 x 3 local block by Gaussian elimination in Q.60 (the diagonal 4 dominates: no pivoting) */
static void solve_block(int64_t J[NF][NF], int64_t R[NF], int64_t d[NF])
{
    for (int k = 0; k < NF; k++) for (int i = k + 1; i < NF; i++) {
        const int64_t f = div60(J[i][k], J[k][k]);
        for (int j = k; j < NF; j++) J[i][j] -= mul60(f, J[k][j]);
        R[i] -= mul60(f, R[k]);
    }
    for (int i = NF - 1; i >= 0; i--) {
        int64_t s = -R[i]; for (int j = i + 1; j < NF; j++) s -= mul60(J[i][j], d[j]);
        d[i] = div60(s, J[i][i]);
    }
}
typedef struct { int color, i_begin, i_end; int64_t maxchg, maxres; } rjob;
static void sweep_rows(void *p)
{
    rjob *jb = (rjob *)p; const int color = jb->color; int64_t maxchg = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1 + ((i + color) & 1); j < Gg; j += 2) {
        int64_t R[NF], J[NF][NF], d[NF];
        for (int k = 0; k < NF; k++) {
            const int64_t S = U[k][i - 1][j] + U[k][i + 1][j] + U[k][i][j - 1] + U[k][i][j + 1];
            R[k] = 4 * U[k][i][j] - S + mul60(H2, source(k, i, j));
            for (int m = 0; m < NF; m++) J[k][m] = (m == k) ? 4 * ONE60 + 2 * LAM * mul60(H2, U[k][i][j]) : H2;
        }
        if (JAC == 2) solve_block(J, R, d);
        else if (JAC == 1) { for (int k = 0; k < NF; k++) d[k] = -div60(R[k], J[k][k]); }
        else { for (int k = 0; k < NF; k++) d[k] = -rdiv64(R[k], 4); }
        for (int k = 0; k < NF; k++) {
            const int64_t step = mul60(OMEGA, d[k]); int64_t ch = step < 0 ? -step : step; if (ch > maxchg) maxchg = ch;
            U[k][i][j] += step;
        }
    }
    jb->maxchg = maxchg;
}
static void resid_rows(void *p)
{
    rjob *jb = (rjob *)p; int64_t maxres = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j < Gg; j++) for (int k = 0; k < NF; k++) {
        const int64_t lap = U[k][i - 1][j] + U[k][i + 1][j] + U[k][i][j - 1] + U[k][i][j + 1] - 4 * U[k][i][j];
        s128 r = {0, 0}; s128_addp(&r, lap, G2g); s128_addv(&r, -source(k, i, j));
        int64_t rv = s128_round_sh(r, 0); if (rv < 0) rv = -rv; if (rv > maxres) maxres = rv;
    }
    jb->maxres = maxres;
}

int main(int argc, char **argv)
{
    int G = 64, jac = 2, scale = 4; int64_t lam = 1, omega_in = 0, thr = 1;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "G=", 2)) G = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "lam=", 4)) lam = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "scale=", 6)) scale = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "thr=", 4)) thr = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "jac=", 4)) jac = !strcmp(argv[i] + 4, "block") ? 2 : !strcmp(argv[i] + 4, "diag") ? 1 : 0;
        else if (!strncmp(argv[i], "omega=", 6)) { if (strcmp(argv[i] + 6, "auto")) omega_in = parse60(argv[i] + 6); }
    }
    if (G < 4 || G > MAXG || (G & (G - 1))) { printf("WHRELAXS: G must be a power of two in 4..%d\n", MAXG); return 1; }
    char A[48], B[48], C[48], E[48];
    const int64_t hq = ONE60 / G; H2 = ONE60 / ((int64_t)G * G); LAM = lam;
    for (int k = 0; k < NF; k++) {                                   /* the exact solution and its sources on the grid */
        const int64_t ck = (ONE60 / scale) * (k + 1);
        for (int i = 0; i <= G; i++) for (int j = 0; j <= G; j++) {
            const int64_t x = hq * i, y = hq * j, xa = mul60(x, ONE60 - x), ya = mul60(y, ONE60 - y);
            Ex[k][i][j] = mul60(ck, mul60(xa, ya)); U[k][i][j] = 0;
            Fv[k][i][j] = -2 * mul60(ck, xa + ya);                     /* the Laplacian of u_k, exactly */
        }
    }
    for (int k = 0; k < NF; k++) for (int i = 0; i <= G; i++) for (int j = 0; j <= G; j++) {   /* f_k = Lap u_k - lam u_k^2 - SUM u_m */
        Fv[k][i][j] -= LAM * mul60(Ex[k][i][j], Ex[k][i][j]);
        for (int m = 0; m < NF; m++) if (m != k) Fv[k][i][j] -= Ex[m][i][j];
    }
    const int64_t s = sin60_small(PI60 / G);
    int64_t omega = div60(2 * ONE60, ONE60 + s); if (omega_in) omega = omega_in;
    fmt60(omega, A); sprintf(B, "%lld", (long long)lam); fmt60(H2, C);
    printf("WHRELAXS: Q.%d fixed point; %d coupled fields on the unit square, G = %d (h^2 = %s), Laplacian u_k = f_k + lam u_k^2 + SUM_(m != k) u_m with lam = %s; red-black Newton-SOR, the Newton step over the %s, omega = %s%s; exact solution u_k = (k/%d) x(1-x) y(1-y)\n",
           F60, NF, G, C, B, jac == 2 ? "local 3 x 3 Jacobian block" : jac == 1 ? "Jacobian's diagonal only (fields lagged)" : "plain Gauss-Seidel diagonal 4 (sources lagged)", A, omega_in ? " (given)" : "", scale);
    const int64_t G2 = (int64_t)G * G;
    uint64_t chain = 0xcbf29ce484222325ull;
    int sweep; int64_t maxres = 0, maxchg = 0, best = -1; int lastimp = 0, floor_sweep = 0, blew = 0;
    Gg = G; JAC = jac; G2g = G2; OMEGA = omega;
    const int mid = 1 + (G - 1) / 2;
    for (sweep = 1; sweep <= 100000; sweep++) {
        maxchg = 0;
        for (int color = 0; color < 2; color++) {
            rjob a = { color, 1, mid, 0, 0 }, b = { color, mid, G, 0, 0 };
            par2(sweep_rows, &a, &b);
            if (a.maxchg > maxchg) maxchg = a.maxchg; if (b.maxchg > maxchg) maxchg = b.maxchg;
        }
        { rjob a = { 0, 1, mid, 0, 0 }, b = { 0, mid, G, 0, 0 };
          par2(resid_rows, &a, &b);
          maxres = a.maxres > b.maxres ? a.maxres : b.maxres; }
        chain = fnv(chain, &maxres, 8);
        if (best < 0 || maxres < best) { best = maxres; lastimp = sweep; floor_sweep = sweep; }
        if ((sweep & (sweep - 1)) == 0 || sweep % 1000 == 0) printf("WHRELAXS sweep %d: max |residual| = %lld units of 2^-60 (%.3e), max change = %lld units\n", sweep, (long long)maxres, (double)maxres / 1152921504606846976.0, (long long)maxchg);
        if (maxchg > 2 * ONE60) { blew = 1; break; }                         /* a change above 2 in fields below 1: diverged */
        if (maxres <= thr || maxchg == 0 || (sweep > 8 * G && sweep - lastimp >= G)) break;
    }
    if (blew) printf("WHRELAXS: DIVERGED (the largest change passed 2 at sweep %d)\n", sweep);
    else printf("WHRELAXS: the residual reached its floor (%lld units, about %lld x G^2) at sweep %d\n", (long long)best, (long long)(best / G2), floor_sweep);
    int64_t maxerr[NF] = {0, 0, 0};
    for (int k = 0; k < NF; k++) for (int i = 0; i <= G; i++) for (int j = 0; j <= G; j++) { int64_t e = U[k][i][j] - Ex[k][i][j]; if (e < 0) e = -e; if (e > maxerr[k]) maxerr[k] = e; }
    fmt60(U[0][G / 2][G / 2], A); fmt60(Ex[0][G / 2][G / 2], B); fmt60(maxerr[0], C); fmt60(maxerr[1], E);
    printf("WHRELAXS: %s after %d sweeps: max |residual| = %lld units, last max change %lld units; u_1(1/2, 1/2) = %s (exact %s); max |u_h - u| = %s, %s, %lld units for the three fields\n",
           blew ? "stopped" : "converged", sweep, (long long)maxres, (long long)maxchg, A, B, C, E, (long long)maxerr[2]);
    printf("WHRELAXS: chain %016llx\n", (unsigned long long)chain);
    printf("WHRELAXS: done\n");
    return 0;
}
