/* wh_dtk.c (11-12-2026, 6:50 AM; v2 7:10 AM) -- WORMHOLE MATH 1b, the kernel for the two-dimensional production problems (P1, the mouths;
 * P4) in the Einstein-DeTurck form (G2) of the reviewer's stage 1 (b.5): five unknowns on a grid in (eta, zeta), each equation's principal
 * part the SAME variable-coefficient operator h^ab d_a d_b (the inverse metric of the orbit space), the rest first order and lagged. Built
 * ahead of the equations, as the P0 kernel was: this file holds the machinery the equations will not change and proves it on a manufactured
 * problem.
 *   - Five fields on i along eta (periodic, the circle; or two reflection lines at eta = 0 and pi with a parity per field: Neumann for the
 *     even fields, zero for the odd ones) and j along zeta = ln z (Dirichlet at the cutoff j = 0; Dirichlet or Neumann at the deep end j = N).
 *   - The principal part at a point uses the REFERENCE metric's inverse (fixed per point, so the weights w_E, w_N, w_X and the row scale 1/D
 *     are formed once, with the only divisions of the run): h^ee F_ee/De^2 + 2 h^ez F_ez/(4 De Dz) + h^zz F_zz/Dz^2, D = 2 h^ee/De^2 +
 *     2 h^zz/Dz^2; the cross term couples the four diagonal neighbours (+NE +SW -NW -SE). The difference between the current and the reference
 *     inverse metric times the second differences goes to the lagged source in the production problem.
 *   - Deviation variables (ref=1): the unknowns are D = F - R for a reference table R (any smooth table, not an exact solution); the
 *     operator on R is formed once per point (LR = SUM w R_nb + w_X cross(R) - R), the sources are evaluated on the full field R + D, and the
 *     Dirichlet rows hold the data minus the reference. With defect=1 the reference's own discrete defect is subtracted (the P0 form: zero
 *     deviation is then the exact discrete solution of the reference's problem; only right when the reference solves the continuum problem).
 *   - Red-black Gauss-Seidel over (i + j) parity, all five fields updated from one state at a point; the diagonal neighbours have the point's
 *     own colour, so the cross term is Jacobi-like, which is fine while it is small. A fixed omega with the change guard; the P0 stops.
 *   - Q.60 end to end, every product and quotient through 128 bits and rounded once; the coordinates below 8 (eta up to 2 pi, zeta in [-1, 1]).
 * THE MANUFACTURED TEST: the reference inverse metric h^ee = 1 + p sin(eta) cos(zeta), h^zz = 1 + p cos(eta) sin(zeta), h^ez = q sin(eta) sin(zeta)
 * (positive definite for p + q < 1); the exact fields F_k = A_k cos(m_k eta + c_k) g(zeta)/8 (periodic) or A_k cos(m_k eta) g/8 for the even
 * and A_k sin(m_k eta) g/8 for the odd fields (reflect), g = 1 - zeta^2 (Dirichlet ends) or (1 - zeta^2)^2 / 2 (Neumann deep end; the full square's g'' = 8 at the end is past the ceiling); the source
 * S_k = h^ab d_a d_b F_k from the continuum derivatives in the same arithmetic, so the discrete error is the stencil's h^2, measured by grid
 * doubling; an optional lagged coupling lam (F_k^2 + SUM_{m != k} F_m). In the test with ref=1 the reference is R = (1 + rp) F + rp g/8 (rp given),
 * so the deviation to find is -rp F - rp g/8.
 * Usage: wh_dtk [Ge=32] [N=32] [bc=periodic|reflect] [top=dirichlet|neumann] [ref=0|1] [rp=0.2] [defect=0|1] [p=0.3] [q=0.2] [lam=0]
 *                [omega=auto|value] [tol=1048576] [pert=0.01] [guard=1] */
#include "whq60.h"
#include <string.h>
#include <stdlib.h>
#include "par2.h"

#define MAXGE 256
#define MAXN 256
#define NF 5
static int GE, N, REFLECT = 0, TOPNEU = 0, REF = 0, DEFECT = 0; static int64_t DE, DZ, OMEGA, LAM, PP, QQ, RP;
static const int PAR[NF] = { 1, 1, -1, 1, -1 };                     /* the parity of each field under the reflections (reflect mode): + Neumann, - zero */
static int64_t F[NF][MAXGE + 1][MAXN + 1], EX[NF][MAXGE + 1][MAXN + 1], SRC[NF][MAXGE + 1][MAXN + 1], CKP[NF][MAXGE + 1][MAXN + 1];
static int64_t R[NF][MAXGE + 1][MAXN + 1], LR[NF][MAXGE + 1][MAXN + 1], DEF[NF][MAXGE + 1][MAXN + 1];   /* the reference, the operator on it (correction units), its defect */
static int64_t WE[MAXGE + 1][MAXN + 1], WN[MAXGE + 1][MAXN + 1], WX[MAXGE + 1][MAXN + 1], WD[MAXGE + 1][MAXN + 1];   /* per point: the weights (E = W, N = S), the cross weight, 1/D */
static int64_t HEE[MAXGE + 1][MAXN + 1], HZZ[MAXGE + 1][MAXN + 1], HEZ[MAXGE + 1][MAXN + 1];                        /* the reference inverse metric */
static inline int ie(int i) { return i < 0 ? i + GE : (i >= GE ? i - GE : i); }                  /* periodic eta: i = GE is i = 0 */
static inline int64_t Av(int64_t (*A)[MAXGE + 1][MAXN + 1], int k, int i, int j)                   /* a field's value with the ghosts of both modes */
{
    if (j < 0) j = 0; else if (j > N) j = N - 1;                                                   /* j = N + 1 appears only with the Neumann deep end: the mirror of N - 1 */
    if (!REFLECT) return A[k][ie(i)][j];
    if (i < 0) return PAR[k] * A[k][-i][j];
    if (i > GE) return PAR[k] * A[k][2 * GE - i][j];
    return A[k][i][j];
}
static inline int64_t Fv(int k, int i, int j) { return Av(F, k, i, j); }
static inline int64_t zof(int j) { return -ONE60 + (int64_t)j * DZ; }                               /* zeta_j = -1 + j Dz, Dz = 2/N */
/* the manufactured exact fields and their continuum derivatives */
static const int MK[NF] = { 1, 2, 1, 3, 2 }; static const int64_t CK[NF] = { 0, ONE60 / 3, 2 * (ONE60 / 3), ONE60 / 5, ONE60 / 7 }; static const int AK[NF] = { 3, 2, -2, 1, -1 };
static void gfun(int64_t z, int64_t *g, int64_t *gz, int64_t *gzz)
{
    const int64_t h = ONE60 - mul60(z, z);
    if (!TOPNEU) { *g = h; *gz = -2 * z; *gzz = -2 * ONE60; }
    else { *g = mul60(h, h) / 2; *gz = -2 * mul60(z, h); *gzz = -2 * ONE60 + 6 * mul60(z, z); }   /* (1 - zeta^2)^2 / 2: with the full square g'' is exactly 8 at the deep end, the one value Q.60 cannot hold (found 7:25 AM: the source wrapped to -2.12) */
}
static void exact(int k, int i, int64_t z, int64_t *f, int64_t *fee, int64_t *fez, int64_t *fzz)   /* the angle m eta reduced on the grid: (i m mod period) De stays below 2 pi */
{
    const int per = REFLECT ? 2 * GE : GE;                                                         /* the full period in grid steps */
    int64_t c, s; cossin60((int64_t)((i * MK[k]) % per) * DE + (REFLECT ? 0 : CK[k]), &c, &s);
    const int odd = REFLECT && PAR[k] < 0;
    const int64_t amp = AK[k] * (ONE60 / 8), ang = odd ? s : c, dang = odd ? c : -s;                 /* the angular factor and its derivative over m */
    int64_t g, gz, gzz; gfun(z, &g, &gz, &gzz);
    *f = mul60(mul60(amp, ang), g);
    *fee = -MK[k] * MK[k] * mul60(mul60(amp, ang), g);
    *fez = MK[k] * mul60(mul60(amp, dang), gz);
    *fzz = mul60(mul60(amp, ang), gzz);
}
static inline int64_t source(int k, int i, int j)       /* the lagged source on the full field */
{
    int64_t s = SRC[k][i][j];
    if (LAM) { const int64_t fk = F[k][i][j] + (REF ? R[k][i][j] : 0); s += mul60(LAM, mul60(fk, fk)); for (int m = 0; m < NF; m++) if (m != k) s += mul60(LAM, F[m][i][j] + (REF ? R[m][i][j] : 0)); }
    return s;
}
static inline int64_t opw(int64_t (*A)[MAXGE + 1][MAXN + 1], int k, int i, int j)   /* SUM w A_nb + w_X cross(A), rounded once: the operator's neighbour part in correction units */
{
    s128 r = {0, 0};
    s128_addp(&r, WE[i][j], Av(A, k, i + 1, j) + Av(A, k, i - 1, j));
    s128_addp(&r, WN[i][j], Av(A, k, i, j + 1) + Av(A, k, i, j - 1));
    s128_addp(&r, WX[i][j], Av(A, k, i + 1, j + 1) + Av(A, k, i - 1, j - 1) - Av(A, k, i - 1, j + 1) - Av(A, k, i + 1, j - 1));
    return s128_round_sh(r, F60);
}
static inline int64_t correction(int k, int i, int j)       /* d = [SUM w F_nb + w_X cross] - S/D + LR - DEF - F, each piece rounded once */
{
    s128 r = {0, 0};
    s128_addp(&r, WE[i][j], Fv(k, i + 1, j) + Fv(k, i - 1, j));
    s128_addp(&r, WN[i][j], Fv(k, i, j + 1) + Fv(k, i, j - 1));
    s128_addp(&r, WX[i][j], Fv(k, i + 1, j + 1) + Fv(k, i - 1, j - 1) - Fv(k, i - 1, j + 1) - Fv(k, i + 1, j - 1));
    s128_addp(&r, -WD[i][j], source(k, i, j));
    return s128_round_sh(r, F60) + (REF ? LR[k][i][j] : 0) - (DEFECT ? DEF[k][i][j] : 0) - F[k][i][j];
}
static inline int fixed_point(int k, int i, int j)         /* a point that is data, not an unknown */
{
    if (j == 0) return 1; if (j == N && !TOPNEU) return 1;
    if (REFLECT && (i == 0 || i == GE) && PAR[k] < 0) return 1;
    return 0;
}
typedef struct { int color, i_begin, i_end; int64_t maxchg, maxcorr, maxres; } rjob;
static void sweep_rows(void *p)
{
    rjob *jb = (rjob *)p; const int color = jb->color; int64_t maxchg = 0, maxcorr = 0;
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j <= N; j++) if (((i + j) & 1) == color) {
        int64_t d[NF];
        for (int k = 0; k < NF; k++) d[k] = fixed_point(k, i, j) ? 0 : correction(k, i, j);
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
    for (int i = jb->i_begin; i < jb->i_end; i++) for (int j = 1; j <= N; j++) for (int k = 0; k < NF; k++) {
        if (fixed_point(k, i, j)) continue;
        const int64_t d = correction(k, i, j); const uint64_t ad = (uint64_t)(d < 0 ? -d : d);
        const uint64_t hi = ad >> 4, lo = ad << 60; if (hi >= (uint64_t)WD[i][j]) { maxres = INT64_MAX; continue; }
        uint64_t rem; const uint64_t q = div128(hi, lo, (uint64_t)WD[i][j], &rem); const int64_t r = q > (uint64_t)INT64_MAX ? INT64_MAX : (int64_t)q;
        if (r > maxres) maxres = r;
    }
    jb->maxres = maxres;
}
int main(int argc, char **argv)
{
    GE = 32; N = 32; int guard = 1; int64_t lam = 0, omega_in = 0, tol = (int64_t)1 << 20, pert = ONE60 / 100; PP = parse60("0.3"); QQ = parse60("0.2"); RP = parse60("0.2");
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "Ge=", 3)) GE = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "bc=", 3)) REFLECT = !strcmp(argv[i] + 3, "reflect");
        else if (!strncmp(argv[i], "top=", 4)) TOPNEU = !strcmp(argv[i] + 4, "neumann");
        else if (!strncmp(argv[i], "ref=", 4)) REF = atoi(argv[i] + 4);
        else if (!strncmp(argv[i], "rp=", 3)) RP = parse60(argv[i] + 3);
        else if (!strncmp(argv[i], "defect=", 7)) DEFECT = atoi(argv[i] + 7);
        else if (!strncmp(argv[i], "p=", 2)) PP = parse60(argv[i] + 2);
        else if (!strncmp(argv[i], "q=", 2)) QQ = parse60(argv[i] + 2);
        else if (!strncmp(argv[i], "lam=", 4)) lam = parse60(argv[i] + 4);
        else if (!strncmp(argv[i], "tol=", 4)) tol = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "pert=", 5)) pert = parse60(argv[i] + 5);
        else if (!strncmp(argv[i], "guard=", 6)) guard = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "omega=", 6)) { if (strcmp(argv[i] + 6, "auto")) omega_in = parse60(argv[i] + 6); }
    }
    if (GE < 4 || GE > MAXGE || (GE & 1) || N < 4 || N > MAXN || (N & 1)) { printf("WHDTK: Ge and N must be even, 4..%d and 4..%d\n", MAXGE, MAXN); return 1; }
    if (DEFECT && !REF) { printf("WHDTK: defect=1 needs ref=1\n"); return 1; }
    const int IMAX = REFLECT ? GE : GE - 1;                                                           /* the last unknown column */
    DE = rdiv64(REFLECT ? PI60 : 2 * PI60, GE); DZ = rdiv64(2 * ONE60, N); LAM = lam;
    char a[48], b[48], c[48], d[48];
    /* the reference inverse metric and the weights, per point (the lines i = 0 and i = Ge included) */
    int64_t hmin = 7 * ONE60;
    for (int i = 0; i <= GE; i++) for (int j = 0; j <= N; j++) {
        int64_t ce, se, cz, sz; cossin60((int64_t)i * DE, &ce, &se); cossin60(zof(j), &cz, &sz);
        HEE[i][j] = ONE60 + mul60(PP, mul60(se, cz)); HZZ[i][j] = ONE60 + mul60(PP, mul60(ce, sz)); HEZ[i][j] = REFLECT ? 0 : mul60(QQ, mul60(se, sz));   /* reflect: h^ez is odd and zero on the lines; the test keeps it zero */
        const int64_t det = mul60(HEE[i][j], HZZ[i][j]) - mul60(HEZ[i][j], HEZ[i][j]); if (det < hmin) hmin = det;
        const int64_t de2 = mul60(DE, DE), dz2 = mul60(DZ, DZ);
        const int64_t tee = 2 * mul60(HEE[i][j], dz2), tzz = 2 * mul60(HZZ[i][j], de2), dd = tee + tzz;   /* D De^2 Dz^2 (below 8) */
        WE[i][j] = div60(mul60(HEE[i][j], dz2), dd); WN[i][j] = div60(mul60(HZZ[i][j], de2), dd);
        WX[i][j] = div60(mul60(HEZ[i][j], mul60(DE, DZ)) / 2, dd); WD[i][j] = div60(mul60(de2, dz2), dd);
    }
    /* the exact fields, the sources, the reference, the Dirichlet rows, the start */
    for (int k = 0; k < NF; k++) for (int i = 0; i <= GE; i++) for (int j = 0; j <= N; j++) {
        int64_t f, fee, fez, fzz; exact(k, i, zof(j), &f, &fee, &fez, &fzz);
        EX[k][i][j] = f;
        s128 s = {0, 0}; s128_addp(&s, HEE[i][j], fee); s128_addp(&s, 2 * HEZ[i][j], fez); s128_addp(&s, HZZ[i][j], fzz);
        SRC[k][i][j] = s128_round_sh(s, F60);
        if (LAM) SRC[k][i][j] -= mul60(LAM, mul60(f, f));
        int64_t g, gz, gzz; gfun(zof(j), &g, &gz, &gzz);
        R[k][i][j] = REF ? f + mul60(RP, f) + mul60(RP, g / 8) : 0;                                   /* the reference: (1 + rp) F + rp g/8 */
    }
    if (LAM) for (int k = 0; k < NF; k++) for (int i = 0; i <= GE; i++) for (int j = 0; j <= N; j++) for (int m = 0; m < NF; m++) if (m != k) SRC[k][i][j] -= mul60(LAM, EX[m][i][j]);
    for (int k = 0; k < NF; k++) for (int i = 0; i <= GE; i++) for (int j = 0; j <= N; j++) {
        const int64_t target = EX[k][i][j] - R[k][i][j];                                              /* the exact deviation (the field itself when ref = 0) */
        F[k][i][j] = fixed_point(k, i, j) ? target : mul60(pert, target);
    }
    if (REF) for (int k = 0; k < NF; k++) for (int i = 0; i <= IMAX; i++) for (int j = 1; j <= N; j++) {
        LR[k][i][j] = opw(R, k, i, j) - R[k][i][j];                                                    /* the operator on the reference, correction units */
        if (DEFECT) { s128 r = {0, 0}; s128_addp(&r, -WD[i][j], SRC[k][i][j]); DEF[k][i][j] = LR[k][i][j] + s128_round_sh(r, F60); }   /* at zero deviation: the correction the reference leaves */
    }
    const int nn = GE > N ? GE : N; int64_t omega = nn > 90 ? 2 * ONE60 - (ONE60 / nn) * 90 : ONE60; if (omega_in) omega = omega_in; OMEGA = omega;
    fmt60(PP, a); fmt60(QQ, b); fmt60(hmin, c); fmt60(omega, d); char lamb[48], rpb[48]; fmt60(lam, lamb); fmt60(RP, rpb);
    printf("WHDTK: Q.60 fixed point; %d fields on eta %s (Ge = %d) by zeta in [-1, 1] (N = %d; the cutoff Dirichlet, the deep end %s); the principal part h^ab d_a d_b with the reference inverse metric h^ee = 1 + p sin eta cos zeta, h^zz = 1 + p cos eta sin zeta, h^ez = q sin eta sin zeta (p = %s, q = %s; the smallest determinant %s); per-point weights formed once; red-black with all fields from one state; omega %s%s; the guard %s; lam = %s; %s\n",
           NF, REFLECT ? "in [0, pi] with reflection lines (parity + + - + -)" : "periodic (De = 2 pi/Ge)", GE, N, TOPNEU ? "Neumann" : "Dirichlet", a, b, c, d, omega_in ? " (given)" : " (2 - 90/max(Ge, N) or 1)", guard ? "on" : "off", lamb,
           REF ? (DEFECT ? "deviations from the reference (1 + rp) F + rp g/8 with its defect subtracted" : "deviations from the reference (1 + rp) F + rp g/8") : "the fields themselves");
    if (REF) printf("WHDTK: rp = %s\n", rpb);
    if (getenv("WHDTK_DBG")) for (int j = 1; j <= N; j++) { const int64_t c0 = correction(0, 1, j); char q1[48]; fmt60(c0, q1); fprintf(stderr, "dbg start row %d: correction of field 0 at column 1 = %s (%lld units); w_N %lld w_E %lld w_X %lld 1/D %lld; S %lld; F %lld F_N %lld F_S %lld\n", j, q1, (long long)c0, (long long)WN[1][j], (long long)WE[1][j], (long long)WX[1][j], (long long)WD[1][j], (long long)SRC[0][1][j], (long long)F[0][1][j], (long long)Fv(0, 1, j + 1), (long long)Fv(0, 1, j - 1)); }
    const int64_t worst = whq60_selftest("WHDTK"); printf("WHDTK selftest: largest deviation %lld units of 2^-60\n", (long long)worst);
    uint64_t chain = 0xcbf29ce484222325ull;
    int sweep; int64_t maxres = -1, maxchg = 0, maxcorr = 0, smallest = -1, ckp_chg = -1; int restores = 0, last_ckp = 0;
    const int mid = (IMAX + 1) / 2;
    memcpy(CKP, F, sizeof F);
    for (sweep = 1; sweep <= 200000; sweep++) {
        maxchg = 0; maxcorr = 0;
        for (int color = 0; color < 2; color++) {
            rjob x = { color, 0, mid, 0, 0, 0 }, y = { color, mid, IMAX + 1, 0, 0, 0 };
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
        if (sweep % 16 == 0 || maxcorr <= tol) { rjob x = { 0, 0, mid, 0, 0, 0 }, y = { 0, mid, IMAX + 1, 0, 0, 0 }; par2(resid_rows, &x, &y); maxres = x.maxres > y.maxres ? x.maxres : y.maxres; }
        if ((sweep & (sweep - 1)) == 0 || sweep % 1000 == 0) printf("WHDTK sweep %d: max correction = %lld units of 2^-60, max change = %lld, max |residual| = %lld\n", sweep, (long long)maxcorr, (long long)maxchg, (long long)maxres);
        if (maxcorr <= tol && maxres >= 0 && maxres <= tol) { printf("WHDTK: stop at sweep %d: the largest correction %lld and residual %lld are at or below tol %lld units\n", sweep, (long long)maxcorr, (long long)maxres, (long long)tol); break; }
        if (maxchg == 0) break;
    }
    int64_t maxerr[NF] = {0, 0, 0, 0, 0};
    for (int k = 0; k < NF; k++) for (int i = 0; i <= IMAX; i++) for (int j = 0; j <= N; j++) { int64_t er = F[k][i][j] + (REF ? R[k][i][j] : 0) - EX[k][i][j]; if (er < 0) er = -er; if (er > maxerr[k]) maxerr[k] = er; }
    fmt60(F[0][GE / 4][N / 2] + (REF ? R[0][GE / 4][N / 2] : 0), a); fmt60(EX[0][GE / 4][N / 2], b); fmt60(OMEGA, c);
    printf("WHDTK: after %d sweeps (omega at the end %s, %d guard restores): max |F - exact| = %lld, %lld, %lld, %lld, %lld units of 2^-60 for the five fields; F_1 at (eta = Ge/4, zeta = 0) = %s (exact %s)\n",
           sweep, c, restores, (long long)maxerr[0], (long long)maxerr[1], (long long)maxerr[2], (long long)maxerr[3], (long long)maxerr[4], a, b);
    if (getenv("WHDTK_DBG")) for (int j = 0; j <= N; j++) { int64_t m = 0, fx = 0, ff = 0; for (int i = 0; i <= IMAX; i++) { int64_t er = F[0][i][j] + (REF ? R[0][i][j] : 0) - EX[0][i][j]; if (er < 0) er = -er; if (er > m) { m = er; fx = EX[0][i][j]; ff = F[0][i][j] + (REF ? R[0][i][j] : 0); } }
        char q1[48], q2[48]; fmt60(fx, q1); fmt60(ff, q2); fprintf(stderr, "dbg row %d: max |F_0 - exact| = %lld units (exact %s, got %s)\n", j, (long long)m, q1, q2); }
    printf("WHDTK: chain %016llx\n", (unsigned long long)chain);
    printf("WHDTK: done\n");
    return 0;
}
