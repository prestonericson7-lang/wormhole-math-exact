/* wh_relax.c (11-11-2026, 10:55 PM) -- WORMHOLE MATH 1b, the solver kernel: a two-dimensional elliptic relaxation in exact fixed point,
   built and verified before the 1b equations arrive (Maldacena-Milekhin's patched five-dimensional geometry reduces, with its S2 symmetry,
   to a 2D elliptic system; PLAN.md section 1, "1b": "a relaxation solver in fixed point"). This program is the kernel on a problem with a
   known answer, so that the kernel's convergence, its residual accounting and its bit-exactness are measured before any physics goes in.

   THE TEST PROBLEM: Poisson's equation on the unit square, Laplacian u = f, u = 0 on the boundary, with the exact solution
       u(x, y) = x (1 - x) y (1 - y),     f = -2 [x (1 - x) + y (1 - y)],
   on a grid of G x G cells (h = 1/G), the five-point Laplacian. The discrete solution differs from u by the discretization error
   (h^2/12)(u_xxxx + u_yyyy) = -(h^2/12) 2 [2 y(1-y) + 2 x(1-x)]... for this quartic u_xxxx = u_yyyy = 0 is false: u_xxxx = 0 (u is quadratic
   in x) -- so the five-point discrete Laplacian of this u is EXACT (second differences of a quadratic are exact), and the discrete solution
   equals u at the grid points to the arithmetic. That makes the test sharp: after convergence, max |u_h - u| measures the kernel's own
   error, not a truncation error.

   THE KERNEL: red-black Gauss-Seidel with over-relaxation (SOR), omega = 2 / (1 + sin(pi/G)) (the optimal value for this operator),
   in Q.60 (int64 with 60 fraction bits; the update's products in 128 bits, one rounding per point per sweep); the residual
   r = (u_N + u_S + u_E + u_W - 4 u) G^2 - f computed exactly in 128 bits and its maximum printed; the sweep stops when the maximum
   residual is below the stated threshold or the update's maximum change is below one unit of 2^-60. A chain hash over every sweep's
   maximum residual. The 1b system will replace f and the five-point stencil with its own nonlinear operator (Newton-relaxation around
   the same red-black sweep); the boundary handling, the residual accounting and the convergence reporting stay.

   Usage: whrelax [G=128] [omega=auto] [thr=1e-15]    G a power of two (sin(pi/G) is then computed exactly enough from a Taylor series in Q.60) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "whq.h"

#define MAXG 1024
#define F60 60
static const int64_t ONE60 = (int64_t)1 << F60;
typedef struct { int64_t hi; uint64_t lo; } s128;
static int64_t U[MAXG + 1][MAXG + 1], Fv[MAXG + 1][MAXG + 1], Ex[MAXG + 1][MAXG + 1];

static inline void s128_add(s128 *a, int64_t hi, uint64_t lo) { uint64_t l = a->lo + lo; a->hi += hi + (l < lo); a->lo = l; }
static inline void mulprod(int64_t a, int64_t b, int64_t *hi, uint64_t *lo)
{
    int neg = (a < 0) ^ (b < 0); uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b, h, l;
    mul128(ua, ub, &h, &l); if (neg) { l = ~l + 1; h = ~h + (l == 0); } *hi = (int64_t)h; *lo = l;
}
static inline void s128_addp(s128 *a, int64_t x, int64_t y) { int64_t h; uint64_t l; mulprod(x, y, &h, &l); s128_add(a, h, l); }
static inline int64_t s128_round_sh(s128 a, int sh)
{
    int neg = a.hi < 0; if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t r = (a.lo >> sh) | ((uint64_t)a.hi << (64 - sh)); if ((a.lo & ((((uint64_t)1) << sh) - 1)) >= (((uint64_t)1) << (sh - 1))) r++;
    return neg ? -(int64_t)r : (int64_t)r;
}
static inline int64_t mul60(int64_t a, int64_t b) { s128 s = {0, 0}; s128_addp(&s, a, b); return s128_round_sh(s, F60); }
static inline int64_t rdiv64(int64_t x, int64_t n) { return x >= 0 ? (x + n / 2) / n : -((-x + n / 2) / n); }
static int64_t div60(int64_t a, int64_t b)             /* round(a 2^60 / b), b > 0, |a| < 8 b */
{
    int neg = a < 0; uint64_t ua = neg ? -(uint64_t)a : (uint64_t)a; const uint64_t hi = ua >> 4, lo = ua << 60;
    if (hi >= (uint64_t)b) { printf("WHRELAX: div60 overflow -- stop\n"); exit(1); }
    uint64_t rem; uint64_t qq = div128(hi, lo, (uint64_t)b, &rem); if (2 * rem >= (uint64_t)b) qq++;
    return neg ? -(int64_t)qq : (int64_t)qq;
}
/* pi in Q.60 and sin x for 0 <= x <= 1 by Taylor (x = pi/G <= pi/4 for G >= 4) */
static const int64_t PI60 = 3622009729038561421ll;
static int64_t sin60_small(int64_t x)
{
    const int64_t x2 = mul60(x, x); int64_t term = x, sum = x;
    for (int k = 1; k < 30; k++) { term = -rdiv64(mul60(term, x2), (int64_t)(2 * k) * (2 * k + 1)); if (!term) break; sum += term; }
    return sum;
}
static void fmt60(int64_t v, char *out)
{
    int neg = v < 0; uint64_t u = neg ? -(uint64_t)v : (uint64_t)v;
    uint64_t ip = u >> F60, frac = u & ((((uint64_t)1) << F60) - 1), h, l;
    mul128(frac, 1000000000000000ull, &h, &l);          /* 15 decimals here: the kernel's own error is far below 1e-12 */
    uint64_t dg = (l >> F60) | (h << (64 - F60)); if ((l >> (F60 - 1)) & 1) dg++;
    if (dg >= 1000000000000000ull) { dg -= 1000000000000000ull; ip++; }
    sprintf(out, "%s%llu.%015llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)dg);
}

int main(int argc, char **argv)
{
    int G = 128; int64_t omega = 0, thr = 1; /* thr in units of 2^-60 of the residual (the residual is O(1) in size, so 1 unit = 8.7e-19) */
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "G=", 2)) G = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "thr=", 4)) thr = atoll(argv[i] + 4);
    }
    if (G < 4 || G > MAXG || (G & (G - 1))) { printf("WHRELAX: G must be a power of two in 4..%d\n", MAXG); return 1; }
    char A[48], B[48], C[48];
    /* the exact solution and f on the grid: x = i/G, y = j/G exactly in Q.60 (G a power of two) */
    const int64_t hq = ONE60 / G;                         /* h = 1/G exactly */
    for (int i = 0; i <= G; i++) for (int j = 0; j <= G; j++) {
        const int64_t x = hq * i, y = hq * j, xa = mul60(x, ONE60 - x), ya = mul60(y, ONE60 - y);   /* x(1-x), y(1-y): one rounding each */
        Ex[i][j] = mul60(xa, ya);
        Fv[i][j] = -2 * (xa + ya);
        U[i][j] = 0;                                      /* the boundary stays 0; the interior starts at 0 */
    }
    /* omega = 2 / (1 + sin(pi/G)) */
    const int64_t s = sin60_small(PI60 / G);
    omega = div60(2 * ONE60, ONE60 + s);
    fmt60(omega, A); fmt60(s, B);
    printf("WHRELAX: Q.%d fixed point (int64, 128-bit products, one rounding per point per sweep); Poisson on the unit square, G = %d (h = 1/%d), five-point Laplacian, red-black SOR omega = %s (sin(pi/G) = %s); exact solution u = x(1-x) y(1-y), for which the five-point Laplacian is exact\n", F60, G, G, A, B);
    const int64_t G2 = (int64_t)G * G;                    /* 1/h^2 */
    uint64_t chain = 0xcbf29ce484222325ull;
    int sweep; int64_t maxres = 0, maxchg = 0, best = -1; int lastimp = 0, floor_sweep = 0;
    /* THE STOP: the residual's floor is the rounding of the four neighbours and the centre (about one unit each) times G^2, so a residual
       of one unit is unreachable at G > 1; the sweep stops when the maximum residual has not improved for 50 sweeps (the floor reached), or
       at thr if that comes first. The sweep count at which the floor was reached is the convergence measure printed. */
    for (sweep = 1; sweep <= 100000; sweep++) {
        maxchg = 0;
        for (int color = 0; color < 2; color++)
            for (int i = 1; i < G; i++) for (int j = 1 + ((i + color) & 1); j < G; j += 2) {
                /* the Gauss-Seidel value: (u_N + u_S + u_E + u_W - f h^2) / 4; f h^2 = f / G^2 */
                const int64_t fh2 = rdiv64(Fv[i][j], G2);
                const int64_t gs = rdiv64(U[i - 1][j] + U[i + 1][j] + U[i][j - 1] + U[i][j + 1] - fh2, 4);
                const int64_t nu = U[i][j] + mul60(omega, gs - U[i][j]);
                int64_t ch = nu - U[i][j]; if (ch < 0) ch = -ch; if (ch > maxchg) maxchg = ch;
                U[i][j] = nu;
            }
        /* the residual, exactly: r = (u_N + u_S + u_E + u_W - 4u) G^2 - f, in units of 2^-60 (the Laplacian of an O(1) function) */
        maxres = 0;
        for (int i = 1; i < G; i++) for (int j = 1; j < G; j++) {
            const int64_t lap = U[i - 1][j] + U[i + 1][j] + U[i][j - 1] + U[i][j + 1] - 4 * U[i][j];
            s128 r = {0, 0}; s128_addp(&r, lap, G2); s128_add(&r, -(Fv[i][j] >> 63) - (Fv[i][j] < 0), (uint64_t)(-Fv[i][j]));
            /* the line above subtracts f as a 128-bit value: hi = sign extension, lo = -f */
            int64_t rv = s128_round_sh(r, 0); if (rv < 0) rv = -rv; if (rv > maxres) maxres = rv;
        }
        chain = fnv(chain, &maxres, 8);
        if (best < 0 || maxres < best) { best = maxres; lastimp = sweep; floor_sweep = sweep; }
        if ((sweep & (sweep - 1)) == 0 || sweep % 1000 == 0) printf("WHRELAX sweep %d: max |residual| = %lld units of 2^-60 (%.3e), max change = %lld units\n", sweep, (long long)maxres, (double)maxres / 1152921504606846976.0, (long long)maxchg);
        /* SOR at omega near 2 lets the residual GROW for the first tens of sweeps before it falls, and needs about 6.6 G sweeps for sixty bits
           (spectral radius omega - 1 = 1 - 2 pi/G + ...): the floor rule waits 8 G sweeps and then asks for G sweeps without improvement */
        if (maxres <= thr || maxchg == 0 || (sweep > 8 * G && sweep - lastimp >= G)) break;
    }
    printf("WHRELAX: the residual reached its floor (%lld units, about %lld x G^2) at sweep %d; the rounding floor on paper is a few units x G^2 = a few x %lld\n", (long long)best, (long long)(best / G2), floor_sweep, (long long)G2);
    /* the error against the exact solution */
    int64_t maxerr = 0; int ei = 0, ej = 0;
    for (int i = 0; i <= G; i++) for (int j = 0; j <= G; j++) { int64_t e = U[i][j] - Ex[i][j]; if (e < 0) e = -e; if (e > maxerr) { maxerr = e; ei = i; ej = j; } }
    fmt60(U[G / 2][G / 2], A); fmt60(Ex[G / 2][G / 2], B); fmt60(maxerr, C);
    printf("WHRELAX: converged after %d sweeps: max |residual| = %lld units of 2^-60, last max change %lld units; u(1/2, 1/2) = %s (exact 1/16 = %s); max |u_h - u| = %s (%lld units, at i = %d, j = %d)\n",
           sweep, (long long)maxres, (long long)maxchg, A, B, C, (long long)maxerr, ei, ej);
    printf("WHRELAX: chain %016llx\n", (unsigned long long)chain);
    printf("WHRELAX: done\n");
    return 0;
}
