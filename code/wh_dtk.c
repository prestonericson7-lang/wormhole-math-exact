/* wh_dtk.c (11-12-2026, 6:50 AM) -- WORMHOLE MATH 1b, the kernel for the two-dimensional production problems (P1, the mouths; P4) in the
 * Einstein-DeTurck form (G2) of the reviewer's stage 1 (b.5): five unknowns on a grid in (eta, zeta), each equation's principal part the
 * SAME variable-coefficient operator h^ab d_a d_b (the inverse metric of the orbit space), the rest first order and lagged. Built ahead of the
 * equations, as the P0 kernel was: this file holds the machinery the equations will not change and proves it on a manufactured problem.
 *   - Five fields F_k(i, j), i along eta (periodic: the circle; or Neumann on symmetry lines), j along zeta = ln z (Dirichlet at both ends
 *     in the test; the production problem has Dirichlet at the cutoff and the regular centre at the other end).
 *   - The principal part at a point uses the REFERENCE metric's inverse (fixed per point, so the weights w_E, w_W, w_N, w_S, w_X and the
 *     row scale 1/D are formed once, with the only divisions of the run): h^ee F_ee/De^2 + 2 h^ez F_ez/(4 De Dz) + h^zz F_zz/Dz^2, with the
 *     diagonal D = 2 h^ee/De^2 + 2 h^zz/Dz^2; the cross term couples the four diagonal neighbours (+NE +SW -NW -SE). The difference between
 *     the current and the reference inverse metric times the second differences goes to the lagged source in the production problem.
 *   - Red-black Gauss-Seidel over (i + j) parity, all five fields updated from one state at a point (the P0 rule: never a field alone);
 *     the diagonal neighbours have the point's own colour, so the cross term is Jacobi-like, which is fine while it is small.
 *   - A fixed omega with the change guard, the stops of the P0 kernel (correction and 128-bit residual at or below tol).
 *   - Q.60 end to end, every product and quotient through 128 bits and rounded once; the row coordinates below 8.
 * THE MANUFACTURED TEST: the reference inverse metric h^ee = 1 + p sin(eta) cos(zeta), h^zz = 1 + p cos(eta) sin(zeta), h^ez = q sin(eta) sin(zeta)
 * (positive definite for p + q < 1), the exact fields F_k = A_k cos(m_k eta + c_k) (1 - zeta^2) / 8 with different m_k, c_k, and the source
 * S_k = h^ab d_a d_b F_k computed from the continuum derivatives in the same arithmetic (so the discrete error is the stencil's h^2, measured
 * by grid doubling) plus an optional lagged coupling lam (F_k^2 + SUM_{m != k} F_m) as P0's test had. Domain eta in [0, 2 pi) periodic, zeta in
 * [-1, 1]. Usage: wh_dtk [Ge=32] [N=32] [p=0.3] [q=0.2] [lam=0] [omega=auto|value] [tol=1048576] [pert=0.01] [guard=1] */
#include "whq60.h"
#include <string.h>
#include <stdlib.h>
#include "par2.h"

#define MAXGE 256
#define MAXN 256
#define NF 5
static int GE, N; static int64_t DE, DZ, OMEGA, LAM, PP, QQ;
static int64_t F[NF][MAXGE + 1][MAXN + 1], EX[NF][MAXGE + 1][MAXN + 1], SRC[NF][MAXGE + 1][MAXN + 1], CKP[NF][MAXGE + 1][MAXN + 1];
static int64_t WE[MAXGE + 1][MAXN + 1], WN[MAXGE + 1][MAXN + 1], WX[MAXGE + 1][MAXN + 1], WD[MAXGE + 1][MAXN + 1];   /* per point: the weights (E = W, N = S), the cross weight, 1/D */
static int64_t HEE[MAXGE + 1][MAXN + 1], HZZ[MAXGE + 1][MAXN + 1], HEZ[MAXGE + 1][MAXN + 1];                        /* the reference inverse metric */
static inline int ie(int i) { return i < 0 ? i + GE : (i >= GE ? i - GE : i); }                  /* periodic eta: i = GE is i = 0 */
static inline int64_t Fv(int k, int i, int j) { return F[k][ie(i)][j]; }                            /* j stays in 0..N (Dirichlet ends) */
static inline int64_t etaof(int i) { return (int64_t)i * DE; }                                      /* eta_i = i De, De = 2 pi / GE */
static inline int64_t zof(int j) { return -ONE60 + (int64_t)j * DZ; }                               /* zeta_j = -1 + j Dz, Dz = 2 / N */
/* the manufactured exact fields and their continuum derivatives: F_k = A_k cos(m_k eta + c_k) (1 - zeta^2)/8 */
static const int MK[NF] = { 1, 2, 1, 3, 2 }; static const int64_t CK[NF] = { 0, ONE60 / 3, 2 * (ONE60 / 3), ONE60 / 5, ONE60 / 7 }; static const int AK[NF] = { 3, 2, -2, 1, -1 };
static void exact(int k, int i, int64_t z, int64_t *f, int64_t *fee, int64_t *fez, int64_t *fzz)   /* the angle m eta reduced on the grid: (i m mod Ge) De stays below 2 pi (m eta itself would pass the Q.60 ceiling) */
{
    int64_t c, s; cossin60((int64_t)((i * MK[k]) % GE) * DE + CK[k], &c, &s);
    const int64_t amp = AK[k] * (ONE60 / 8), g = ONE60 - mul60(z, z), gz = -2 * z, gzz = -2 * ONE60;
    *f = mul60(mul60(amp, c), g);
    *fee = -MK[k] * MK[k] * mul60(mul60(amp, c), g);
    *fez = -MK[k] * mul60(mul60(amp, s), gz);
    *fzz = mul60(mul60(amp, c), gzz);
}
static inline int64_t source(int k, int i, int j)
{
    int64_t s = SRC[k][i][j];
    if (LAM) { s += mul60(LAM, mul60(F[k][i][j], F[k][i][j])); for (int m = 0; m < NF; m++) if (m != k) s += mul60(LAM, F[m][i][j]); }
    return s;
}
static inline int64_t correction(int k, int i, int j)       /* d = w_E (F_E + F_W) + w_N (F_N + F_S) + w_X (F_NE + F_SW - F_NW - F_SE) - S/D - F, rounded once */
{
    s128 r = {0, 0};
    s128_addp(&r, WE[i][j], Fv(k, i + 1, j) + Fv(k, i - 1, j));
    s128_addp(&r, WN[i][j], Fv(k, i, j + 1) + Fv(k, i, j - 1));
    s128_addp(&r, WX[i][j], Fv(k, i + 1, j + 1) + Fv(k, i - 1, j - 1) - Fv(k, i - 1, j + 1) - Fv(k, i + 1, j - 1));
    s128_addp(&r, -WD[i][j], source(k, i, j));
    return s128_round_sh(r, F60) - F[k][i][j];
}
typedef struct { int color, i_begin, i_end; int64_t maxchg, maxcorr, maxres; } rjob;
static void sweep_rows(void *p)
{
    rjob *jb = (rjob *)p; const int color = jb->color; int64_t maxchg = 0, maxcorr = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j < N; j++) if (((i + j) & 1) == color) {
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
static void resid_rows(void *p)                                 /* the residual |h^ab d_a d_b F - S| in units of 2^-60 = |d| D, saturating */
{
    rjob *jb = (rjob *)p; int64_t maxres = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j < N; j++) for (int k = 0; k < NF; k++) {
        const int64_t d = correction(k, i, j); const uint64_t ad = (uint64_t)(d < 0 ? -d : d);
        const uint64_t hi = ad >> 4, lo = ad << 60; if (hi >= (uint64_t)WD[i][j]) { maxres = INT64_MAX; continue; }
        uint64_t rem; const uint64_t q = div128(hi, lo, (uint64_t)WD[i][j], &rem); const int64_t r = q > (uint64_t)INT64_MAX ? INT64_MAX : (int64_t)q;
        if (r > maxres) maxres = r;
    }
    jb->maxres = maxres;
}
int main(int argc, char **argv)
{
    GE = 32; N = 32; int guard = 1; int64_t lam = 0, omega_in = 0, tol = (int64_t)1 << 20, pert = ONE60 / 100; PP = parse60("0.3"); QQ = parse60("0.2");
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "Ge=", 3)) GE = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "p=", 2)) PP = parse60(argv[i] + 2);
        else if (!strncmp(argv[i], "q=", 2)) QQ = parse60(argv[i] + 2);
        else if (!strncmp(argv[i], "lam=", 4)) lam = parse60(argv[i] + 4);
        else if (!strncmp(argv[i], "tol=", 4)) tol = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "pert=", 5)) pert = parse60(argv[i] + 5);
        else if (!strncmp(argv[i], "guard=", 6)) guard = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "omega=", 6)) { if (strcmp(argv[i] + 6, "auto")) omega_in = parse60(argv[i] + 6); }
    }
    if (GE < 4 || GE > MAXGE || (GE & 1) || N < 4 || N > MAXN || (N & 1)) { printf("WHDTK: Ge and N must be even, 4..%d and 4..%d\n", MAXGE, MAXN); return 1; }
    DE = rdiv64(2 * PI60, GE); DZ = rdiv64(2 * ONE60, N); LAM = lam;
    char a[48], b[48], c[48], d[48];
    /* the reference inverse metric and the weights, per point */
    int64_t hmin = ONE60 * 8;
    for (int i = 0; i < GE; i++) for (int j = 0; j <= N; j++) {
        int64_t ce, se, cz, sz; cossin60(etaof(i), &ce, &se); cossin60(zof(j), &cz, &sz);
        HEE[i][j] = ONE60 + mul60(PP, mul60(se, cz)); HZZ[i][j] = ONE60 + mul60(PP, mul60(ce, sz)); HEZ[i][j] = mul60(QQ, mul60(se, sz));
        const int64_t det = mul60(HEE[i][j], HZZ[i][j]) - mul60(HEZ[i][j], HEZ[i][j]); if (det < hmin) hmin = det;
        /* D = 2 h^ee/De^2 + 2 h^zz/Dz^2: formed as a 128-bit value scaled by De^2 Dz^2 to stay below the ceiling, then the weights as ratios */
        const int64_t de2 = mul60(DE, DE), dz2 = mul60(DZ, DZ);                               /* De^2 = 0.0385 at Ge = 32, Dz^2 = 0.0039 at N = 32 */
        const int64_t tee = 2 * mul60(HEE[i][j], dz2), tzz = 2 * mul60(HZZ[i][j], de2);         /* D De^2 Dz^2 = 2 h^ee Dz^2 + 2 h^zz De^2 (below 8) */
        const int64_t dd = tee + tzz;
        WE[i][j] = div60(mul60(HEE[i][j], dz2), dd);                                              /* w_E = (h^ee/De^2)/D = h^ee Dz^2 / (D De^2 Dz^2) */
        WN[i][j] = div60(mul60(HZZ[i][j], de2), dd);
        WX[i][j] = div60(mul60(HEZ[i][j], mul60(DE, DZ)) / 2, dd);                               /* w_X = (2 h^ez/(4 De Dz))/D = h^ez De Dz / (2 D De^2 Dz^2) */
        WD[i][j] = div60(mul60(de2, dz2), dd);                                                    /* 1/D = De^2 Dz^2 / (D De^2 Dz^2) */
    }
    /* the exact fields, the Dirichlet ends, the sources */
    for (int k = 0; k < NF; k++) for (int i = 0; i < GE; i++) for (int j = 0; j <= N; j++) {
        int64_t f, fee, fez, fzz; exact(k, i, zof(j), &f, &fee, &fez, &fzz);
        EX[k][i][j] = f; F[k][i][j] = (j == 0 || j == N) ? f : 0;
        s128 s = {0, 0}; s128_addp(&s, HEE[i][j], fee); s128_addp(&s, 2 * HEZ[i][j], fez); s128_addp(&s, HZZ[i][j], fzz);
        SRC[k][i][j] = s128_round_sh(s, F60);
        if (LAM) { SRC[k][i][j] -= mul60(LAM, mul60(f, f)); }
    }
    if (LAM) for (int k = 0; k < NF; k++) for (int i = 0; i < GE; i++) for (int j = 0; j <= N; j++) for (int m = 0; m < NF; m++) if (m != k) SRC[k][i][j] -= mul60(LAM, EX[m][i][j]);
    for (int k = 0; k < NF; k++) for (int i = 0; i < GE; i++) for (int j = 1; j < N; j++) F[k][i][j] = mul60(pert, EX[k][i][j]);   /* the start: a fraction of the exact fields */
    for (int k = 0; k < NF; k++) for (int j = 0; j <= N; j++) { F[k][GE][j] = F[k][0][j]; EX[k][GE][j] = EX[k][0][j]; }           /* the periodic copy at i = Ge for the printers */
    const int nn = GE > N ? GE : N; int64_t omega = nn > 90 ? 2 * ONE60 - (ONE60 / nn) * 90 : ONE60; if (omega_in) omega = omega_in; OMEGA = omega;
    fmt60(PP, a); fmt60(QQ, b); fmt60(hmin, c); fmt60(omega, d); char lamb[48]; fmt60(lam, lamb);
    printf("WHDTK: Q.60 fixed point; %d fields on a periodic eta (Ge = %d, De = 2 pi/Ge) by zeta in [-1, 1] (N = %d, Dirichlet ends); the principal part h^ab d_a d_b with the reference inverse metric h^ee = 1 + p sin eta cos zeta, h^zz = 1 + p cos eta sin zeta, h^ez = q sin eta sin zeta (p = %s, q = %s; the smallest determinant %s); per-point weights formed once; red-black with all fields from one state; omega %s%s; the guard %s; lam = %s\n",
           NF, GE, N, a, b, c, d, omega_in ? " (given)" : " (2 - 90/max(Ge, N) or 1)", guard ? "on" : "off", lamb);
    const int64_t worst = whq60_selftest("WHDTK"); printf("WHDTK selftest: largest deviation %lld units of 2^-60\n", (long long)worst);
    uint64_t chain = 0xcbf29ce484222325ull;
    int sweep; int64_t maxres = -1, maxchg = 0, maxcorr = 0, smallest = -1, ckp_chg = -1; int restores = 0, last_ckp = 0;
    const int mid = GE / 2;
    memcpy(CKP, F, sizeof F);
    for (sweep = 1; sweep <= 200000; sweep++) {
        maxchg = 0; maxcorr = 0;
        for (int color = 0; color < 2; color++) {
            rjob x = { color, 0, mid, 0, 0, 0 }, y = { color, mid, GE, 0, 0, 0 };
            par2(sweep_rows, &x, &y);
            if (x.maxchg > maxchg) maxchg = x.maxchg; if (y.maxchg > maxchg) maxchg = y.maxchg;
            if (x.maxcorr > maxcorr) maxcorr = x.maxcorr; if (y.maxcorr > maxcorr) maxcorr = y.maxcorr;
        }
        chain = fnv(chain, &maxcorr, 8);
        if (smallest < 0 || maxchg < smallest) smallest = maxchg;
        if (guard && smallest > 0 && maxchg / 32 > smallest && sweep - last_ckp > 1) {
            memcpy(F, CKP, sizeof F); OMEGA -= ONE60 / 10; if (OMEGA < ONE60) OMEGA = ONE60; restores++;
            fmt60(OMEGA, a); printf("WHDTK guard: the change %lld units exceeded 32 x the smallest %lld at sweep %d: checkpoint restored (sweep %d), omega -> %s\n", (long long)maxchg, (long long)smallest, sweep, last_ckp, a);
            smallest = -1; continue;
        }
        if (sweep % 64 == 0 && (ckp_chg < 0 || maxchg <= ckp_chg)) { memcpy(CKP, F, sizeof F); ckp_chg = maxchg; last_ckp = sweep; }
        if (sweep % 16 == 0 || maxcorr <= tol) { rjob x = { 0, 0, mid, 0, 0, 0 }, y = { 0, mid, GE, 0, 0, 0 }; par2(resid_rows, &x, &y); maxres = x.maxres > y.maxres ? x.maxres : y.maxres; }
        if ((sweep & (sweep - 1)) == 0 || sweep % 1000 == 0) printf("WHDTK sweep %d: max correction = %lld units of 2^-60, max change = %lld, max |residual| = %lld\n", sweep, (long long)maxcorr, (long long)maxchg, (long long)maxres);
        if (maxcorr <= tol && maxres >= 0 && maxres <= tol) { printf("WHDTK: stop at sweep %d: the largest correction %lld and residual %lld are at or below tol %lld units\n", sweep, (long long)maxcorr, (long long)maxres, (long long)tol); break; }
        if (maxchg == 0) break;
    }
    int64_t maxerr[NF] = {0, 0, 0, 0, 0};
    for (int k = 0; k < NF; k++) for (int i = 0; i < GE; i++) for (int j = 0; j <= N; j++) { int64_t er = F[k][i][j] - EX[k][i][j]; if (er < 0) er = -er; if (er > maxerr[k]) maxerr[k] = er; }
    fmt60(F[0][GE / 4][N / 2], a); fmt60(EX[0][GE / 4][N / 2], b); fmt60(OMEGA, c);
    printf("WHDTK: after %d sweeps (omega at the end %s, %d guard restores): max |F - exact| = %lld, %lld, %lld, %lld, %lld units of 2^-60 for the five fields; F_1 at (eta = pi/2, zeta = 0) = %s (exact %s)\n",
           sweep, c, restores, (long long)maxerr[0], (long long)maxerr[1], (long long)maxerr[2], (long long)maxerr[3], (long long)maxerr[4], a, b);
    printf("WHDTK: chain %016llx\n", (unsigned long long)chain);
    printf("WHDTK: done\n");
    return 0;
}
