/* wh_relax.c (11-11-2026, 10:55 PM; v2 nl=1; v3 two cores 11:47 PM; v4 11-12 1:35 AM) -- WORMHOLE MATH 1b, the solver kernel:
   VERSION 4 (what the reviewer's float tests of the 1b boundary problems asked for; the default path is v3 bit for bit, its chains
   unchanged): (1) omega= a decimal or 'auto' (the formula, the default). An adaptive omega keyed to the residual was tried and
   dropped: SOR's residual grows through its first tens of sweeps at any omega, which such a rule reads as divergence, and the
   rule then backed off to 1.0 and crawled (1744 sweeps to a poor floor at G = 128); the production solver takes omega from measured
   sweeps instead;
   (2) robin=1: a Robin condition on the y = 0 side (a brane's Israel condition has this form), u + c u_y = c g with g = x(1 - x) the
   exact u_y there and c HALF the first cell, discretized with the second-order one-sided derivative (exact for the quadratic-in-y test
   solution, so the discrete solution still equals u at the nodes) and ELIMINATED into the first interior row (the ghost-point
   method; c = the whole first cell would cancel that row's diagonal, a singular system, and an explicit boundary update diverges):
   on the uniform grid u_0 = u_2 - 4 u_1 + 2 h g and the row-1 equation (u_E + u_W + 2 u_2 - 8 u_1 + 2 h g) G^2 = f; on the stretched grid
   u_0 = (12 g/G^2 - 16 u_1 + u_2)/9 and 18 G^2 (u_E + u_W) + 4 G^4 u_2 - (36 G^2 + 28 G^4) u_1 + 12 G^2 g = 18 f; u_0 is set from the
   converged row 1 at the end. The residual is measured in units of 2^-60 as before but SATURATES at the int64 instead of stopping:
   on the stretched grid the weights near y = 0 reach G^4, so an early error of 2^-20 in u is a residual of 2^64 units; (3) stretch=1: the rows at y_j = (j/G)^2, dense toward y = 0 (the
   brane side), with the three-point second difference on the non-uniform spacing, exact for quadratics; the equation at row j is
   multiplied by D_j = 2j (2j - 1)(2j + 1) so that every stencil weight is an INTEGER (G^4 (2j + 1), G^4 (2j - 1), G^2 D_j) and the
   Gauss-Seidel value is one 128-by-64-bit division rounded once: no weight is ever rounded. With stretch=1 the Robin row reads
   u_0 = (16 u_1 - u_2)/3 - 4 g/G^2. The nonlinear option (nl=1) stays on the uniform grid.
   The original header follows. WORMHOLE MATH 1b, the solver kernel: a two-dimensional elliptic relaxation in exact fixed point,
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
#include "../par2.h"                           /* v3 (11:47 PM): the second A9 core takes half the rows of every sweep and of every residual pass;
                                                  without it (the PC, or a machine with one core) both halves run here: the same bits */

#define MAXG 1024
#define F60 60
static const int64_t ONE60 = (int64_t)1 << F60;
typedef struct { int64_t hi; uint64_t lo; } s128;
static int64_t U[MAXG + 1][MAXG + 1], Fv[MAXG + 1][MAXG + 1], Ex[MAXG + 1][MAXG + 1];
static int64_t Yq[MAXG + 1], Gx[MAXG + 1];                              /* v4: y of each row (Q.60); g_i = x_i (1 - x_i) for the Robin row */
static int64_t WS[MAXG + 1], WN[MAXG + 1], WX[MAXG + 1], WD[MAXG + 1], DJ[MAXG + 1];   /* v4 stretch: the integer stencil weights of row j */
static int STRETCH = 0, ROBIN = 0; static int64_t HQ;

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
static int64_t parse60(const char *s)                 /* v4: a decimal into Q.60, round(v 2^60) */
{
    int neg = 0; if (*s == '-') { neg = 1; s++; }
    uint64_t ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (uint64_t)(*s++ - '0');
    uint64_t fp = 0, scale = 1; if (*s == '.') { s++; while (*s >= '0' && *s <= '9' && scale < 1000000000000ull) { fp = fp * 10 + (uint64_t)(*s++ - '0'); scale *= 10; } }
    if (ip >= 7) { printf("WHRELAX: a value at or above 7 cannot be held in Q.60 here\n"); exit(1); }
    int64_t v = (int64_t)(ip << F60);
    if (fp) { uint64_t h, l, rem; mul128(fp, (uint64_t)1 << F60, &h, &l); uint64_t qq = div128(h, l, scale, &rem); if (2 * rem >= scale) qq++; v += (int64_t)qq; }
    return neg ? -v : v;
}
static int64_t div128s(s128 num, int64_t den)         /* v4: round(num / den), num signed 128-bit, 0 < den < 2^63, the quotient an int64 */
{
    int neg = num.hi < 0; if (neg) { num.lo = ~num.lo + 1; num.hi = ~num.hi + (num.lo == 0); }
    if ((uint64_t)num.hi >= (uint64_t)den) { printf("WHRELAX: a quotient does not fit -- stop\n"); exit(1); }
    uint64_t rem; uint64_t q = div128((uint64_t)num.hi, num.lo, (uint64_t)den, &rem); if (2 * rem >= (uint64_t)den) q++;
    return neg ? -(int64_t)q : (int64_t)q;
}
static int64_t div128sat(s128 num, int64_t den)       /* v4: as div128s but a quotient that does not fit saturates (the residual measure) */
{
    int neg = num.hi < 0; if (neg) { num.lo = ~num.lo + 1; num.hi = ~num.hi + (num.lo == 0); }
    if ((uint64_t)num.hi >= (uint64_t)den) return INT64_MAX;
    uint64_t rem; uint64_t q = div128((uint64_t)num.hi, num.lo, (uint64_t)den, &rem); if (2 * rem >= (uint64_t)den) q++;
    if (q > (uint64_t)INT64_MAX) return INT64_MAX;
    return neg ? -(int64_t)q : (int64_t)q;
}
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

/* the sweep and the residual pass over a range of rows: pure computation (reads U's other colour and the tables, writes its own points
   or its own job record), so a half can run on core 1 */
static int Gg, NL; static int64_t G2g, OMEGA;
typedef struct { int color, i_begin, i_end; int64_t maxchg, maxres; } rjob;
static void sweep_rows(void *p)
{
    rjob *jb = (rjob *)p; const int color = jb->color; int64_t maxchg = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1 + ((i + color) & 1); j < Gg; j += 2) {
        const int64_t fh2 = rdiv64(Fv[i][j], G2g);
        const int64_t S = U[i - 1][j] + U[i + 1][j] + U[i][j - 1] + U[i][j + 1];
        int64_t gs;
        if (ROBIN && j == 1 && STRETCH) {                /* v4: row 1 with the Robin row eliminated, stretched: 18 G^2 (E + W) + 4 G^4 u_2 + 12 G^2 g - 18 f over 36 G^2 + 28 G^4 */
            const int64_t G4 = G2g * G2g; s128 num = {0, 0};
            s128_addp(&num, 18 * G2g, U[i - 1][1] + U[i + 1][1]); s128_addp(&num, 4 * G4, U[i][2]); s128_addp(&num, 12 * G2g, Gx[i]); s128_addp(&num, -18, Fv[i][1]);
            gs = div128s(num, 36 * G2g + 28 * G4);
        }
        else if (ROBIN && j == 1) {                      /* v4: row 1 with the Robin row eliminated, uniform: (E + W + 2 u_2 + 2 h g - f h^2) / 8 */
            gs = rdiv64(U[i - 1][1] + U[i + 1][1] + 2 * U[i][2] + mul60(2 * HQ, Gx[i]) - fh2, 8);
        }
        else if (STRETCH) {                              /* v4: D_j (u_xx + u_yy) = D_j f with integer weights; one rounding */
            s128 num = {0, 0};
            s128_addp(&num, WX[j], U[i - 1][j] + U[i + 1][j]); s128_addp(&num, WS[j], U[i][j - 1]); s128_addp(&num, WN[j], U[i][j + 1]);
            s128_addp(&num, -DJ[j], Fv[i][j]);
            gs = div128s(num, WD[j]);
        }
        else if (!NL) gs = rdiv64(S - fh2, 4);
        else { const int64_t u = U[i][j], u2h2 = rdiv64(mul60(u, u), G2g);
               const int64_t res = 4 * u + u2h2 - S + fh2, dres = 4 * ONE60 + 2 * rdiv64(u, G2g); gs = u - div60(res, dres); }
        const int64_t nu = U[i][j] + mul60(OMEGA, gs - U[i][j]);
        int64_t ch = nu - U[i][j]; if (ch < 0) ch = -ch; if (ch > maxchg) maxchg = ch;
        U[i][j] = nu;
    }
    jb->maxchg = maxchg;
}
static void resid_rows(void *p)
{
    rjob *jb = (rjob *)p; int64_t maxres = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j < Gg; j++) {
        if (ROBIN && j == 1 && STRETCH) {                /* v4: the eliminated row-1 equation's residual, divided back by 18 */
            const int64_t G4 = G2g * G2g; s128 r = {0, 0};
            s128_addp(&r, 18 * G2g, U[i - 1][1] + U[i + 1][1]); s128_addp(&r, 4 * G4, U[i][2]); s128_addp(&r, 12 * G2g, Gx[i]); s128_addp(&r, -18, Fv[i][1]);
            s128_addp(&r, -(36 * G2g + 28 * G4), U[i][1]);
            int64_t rv = div128sat(r, 18); if (rv < 0) rv = -rv; if (rv > maxres) maxres = rv;
            continue;
        }
        if (ROBIN && j == 1) {                           /* v4: (E + W + 2 u_2 - 8 u_1 + 2 h g) G^2 - f */
            const int64_t lap1 = U[i - 1][1] + U[i + 1][1] + 2 * U[i][2] - 8 * U[i][1] + mul60(2 * HQ, Gx[i]);
            s128 r = {0, 0}; s128_addp(&r, lap1, G2g); s128_add(&r, -(Fv[i][1] >> 63) - (Fv[i][1] < 0), (uint64_t)(-Fv[i][1]));
            int64_t rv = div128sat(r, 1); if (rv < 0) rv = -rv; if (rv > maxres) maxres = rv;
            continue;
        }
        if (STRETCH) {                                   /* v4: the weighted residual divided back by D_j, in units of 2^-60 (saturating) */
            s128 r = {0, 0};
            s128_addp(&r, WX[j], U[i - 1][j] + U[i + 1][j]); s128_addp(&r, WS[j], U[i][j - 1]); s128_addp(&r, WN[j], U[i][j + 1]);
            s128_addp(&r, -WD[j], U[i][j]); s128_addp(&r, -DJ[j], Fv[i][j]);
            int64_t rv = div128sat(r, DJ[j]); if (rv < 0) rv = -rv; if (rv > maxres) maxres = rv;
            continue;
        }
        const int64_t lap = U[i - 1][j] + U[i + 1][j] + U[i][j - 1] + U[i][j + 1] - 4 * U[i][j];
        s128 r = {0, 0}; s128_addp(&r, lap, G2g); s128_add(&r, -(Fv[i][j] >> 63) - (Fv[i][j] < 0), (uint64_t)(-Fv[i][j]));
        if (NL) { const int64_t u2 = mul60(U[i][j], U[i][j]); s128_add(&r, -(u2 >> 63) - (u2 < 0), (uint64_t)(-u2)); }
        int64_t rv = s128_round_sh(r, 0); if (rv < 0) rv = -rv; if (rv > maxres) maxres = rv;
    }
    jb->maxres = maxres;
}

int main(int argc, char **argv)
{
    int G = 128, nl = 0, robin = 0, stretch = 0; int64_t omega = 0, thr = 1, omega_in = 0; /* thr in units of 2^-60 of the residual (the residual is O(1) in size, so 1 unit = 8.7e-19) */
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "G=", 2)) G = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "thr=", 4)) thr = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "nl=", 3)) nl = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "omega=", 6)) { if (strcmp(argv[i] + 6, "auto")) omega_in = parse60(argv[i] + 6); }
        else if (!strncmp(argv[i], "robin=", 6)) robin = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "stretch=", 8)) stretch = atoi(argv[i] + 8);
    }
    if (stretch && nl) { printf("WHRELAX: stretch=1 is the linear test only\n"); return 1; }
    STRETCH = stretch; ROBIN = robin;
    if (robin && nl) { printf("WHRELAX: robin=1 is the linear test only\n"); return 1; }
    /* nl=1: THE NONLINEAR TEST, the shape of 1b's equations (a Laplacian equal to a nonlinear function of the unknown): Laplacian u - u^2 = f - u_e^2
       with the same exact solution u_e; each visit takes one Newton step of the local scalar equation 4u + h^2 u^2 - S + h^2 g = 0 (S the
       four neighbours, g = f - u_e^2): u <- u - (4u + h^2 u^2 - S + h^2 g) / (4 + 2 h^2 u), then the same over-relaxation. The linear path
       (nl=0) is unchanged bit for bit. */
    if (G < 4 || G > MAXG || (G & (G - 1))) { printf("WHRELAX: G must be a power of two in 4..%d\n", MAXG); return 1; }
    char A[48], B[48], C[48];
    /* the exact solution and f on the grid: x = i/G, y = j/G exactly in Q.60 (G a power of two) */
    const int64_t hq = ONE60 / G;                         /* h = 1/G exactly */
    HQ = hq;
    for (int j = 0; j <= G; j++) Yq[j] = stretch ? (ONE60 / ((int64_t)G * G)) * j * j : hq * j;   /* v4: y_j = (j/G)^2 exactly when stretched */
    if (stretch) {                                        /* v4: the integer weights of row j (see the header); G^4 (2j+1) < 2^52 at G = 1024 */
        const int64_t G4 = (int64_t)G * G * (int64_t)G * G;
        for (int j = 1; j < G; j++) { DJ[j] = 2 * (int64_t)j * (2 * (int64_t)j - 1) * (2 * (int64_t)j + 1); WS[j] = G4 * (2 * (int64_t)j + 1); WN[j] = G4 * (2 * (int64_t)j - 1);
                                      WX[j] = (int64_t)G * G * DJ[j]; WD[j] = 2 * WX[j] + WS[j] + WN[j]; }
    }
    for (int i = 0; i <= G; i++) for (int j = 0; j <= G; j++) {
        const int64_t x = hq * i, y = Yq[j], xa = mul60(x, ONE60 - x), ya = mul60(y, ONE60 - y);   /* x(1-x), y(1-y): one rounding each */
        if (j == 0) Gx[i] = xa;
        Ex[i][j] = mul60(xa, ya);
        Fv[i][j] = -2 * (xa + ya);
        if (nl) Fv[i][j] -= mul60(Ex[i][j], Ex[i][j]);   /* g = f - u_e^2 for the nonlinear test */
        U[i][j] = 0;                                      /* the boundary stays 0; the interior starts at 0 */
    }
    /* omega = 2 / (1 + sin(pi/G)) */
    const int64_t s = sin60_small(PI60 / G);
    omega = div60(2 * ONE60, ONE60 + s);
    if (omega_in) omega = omega_in;                       /* v4 */
    fmt60(omega, A); fmt60(s, B);
    printf("WHRELAX: Q.%d fixed point (int64, 128-bit products, one rounding per point per sweep); Poisson on the unit square, G = %d (h = 1/%d), five-point Laplacian, red-black SOR omega = %s (sin(pi/G) = %s); exact solution u = x(1-x) y(1-y), for which the five-point Laplacian is exact\n", F60, G, G, A, B);
    if (omega_in || robin || stretch) { fmt60(omega, A); printf("WHRELAX v4 options: omega %s%s; robin=%d (y = 0 side: u + c u_y = c g, second-order one-sided); stretch=%d (rows at y = (j/G)^2, integer weights)\n", A, omega_in ? " (given)" : "", robin, stretch); }
    const int64_t G2 = (int64_t)G * G;                    /* 1/h^2 */
    uint64_t chain = 0xcbf29ce484222325ull;
    int sweep; int64_t maxres = 0, maxchg = 0, best = -1; int lastimp = 0, floor_sweep = 0;
    /* THE STOP: the residual's floor is the rounding of the four neighbours and the centre (about one unit each) times G^2, so a residual
       of one unit is unreachable at G > 1; the sweep stops when the maximum residual has not improved for 50 sweeps (the floor reached), or
       at thr if that comes first. The sweep count at which the floor was reached is the convergence measure printed. */
    Gg = G; NL = nl; G2g = G2; OMEGA = omega;
    const int mid = 1 + (G - 1) / 2;                      /* the rows split in two halves, one per core */
    for (sweep = 1; sweep <= 100000; sweep++) {
        maxchg = 0;

        for (int color = 0; color < 2; color++) {           /* the Gauss-Seidel value (u_N + u_S + u_E + u_W - f h^2) / 4, or one Newton step for nl=1; over-relaxed */
            rjob a = { color, 1, mid, 0, 0 }, b = { color, mid, G, 0, 0 };
            par2(sweep_rows, &a, &b);
            if (a.maxchg > maxchg) maxchg = a.maxchg; if (b.maxchg > maxchg) maxchg = b.maxchg;
        }
        /* the residual, exactly: r = (u_N + u_S + u_E + u_W - 4u) G^2 - f (- u^2 for nl=1), in units of 2^-60 */
        { rjob a = { 0, 1, mid, 0, 0 }, b = { 0, mid, G, 0, 0 };
          par2(resid_rows, &a, &b);
          maxres = a.maxres > b.maxres ? a.maxres : b.maxres; }
        chain = fnv(chain, &maxres, 8);
        if (best < 0 || maxres < best) { best = maxres; lastimp = sweep; floor_sweep = sweep; }
        if ((sweep & (sweep - 1)) == 0 || sweep % 1000 == 0) printf("WHRELAX sweep %d: max |residual| = %lld units of 2^-60 (%.3e), max change = %lld units\n", sweep, (long long)maxres, (double)maxres / 1152921504606846976.0, (long long)maxchg);
        /* SOR at omega near 2 lets the residual GROW for the first tens of sweeps before it falls, and needs about 6.6 G sweeps for sixty bits
           (spectral radius omega - 1 = 1 - 2 pi/G + ...): the floor rule waits 8 G sweeps and then asks for G sweeps without improvement */
        if (maxres <= thr || maxchg == 0 || (sweep > 8 * G && sweep - lastimp >= G)) break;
    }
    printf("WHRELAX: the residual reached its floor (%lld units, about %lld x G^2) at sweep %d; the rounding floor on paper is a few units x G^2 = a few x %lld\n", (long long)best, (long long)(best / G2), floor_sweep, (long long)G2);
    if (robin) for (int i = 1; i < G; i++)                   /* v4: the Robin row from the converged row 1 (u_e = 0 there) */
        U[i][0] = stretch ? rdiv64(rdiv64(12 * Gx[i], G2) - 16 * U[i][1] + U[i][2], 9) : U[i][2] - 4 * U[i][1] + mul60(2 * hq, Gx[i]);
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
