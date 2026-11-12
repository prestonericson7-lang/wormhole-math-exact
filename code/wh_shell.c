/* wh_shell.c (11-11-2026) -- WORMHOLE MATH 3b, part 1: the 2024 "warp shell" WITHOUT its shift: Fuchs et al.'s matter shell, its
 * construction reproduced as finite integer sums, the metric DECLARED as the piecewise interpolant of those tables, the Einstein tensor
 * of that object in closed form at every interval's midpoint, the four energy conditions, and two rules for the cavity.
 *
 * Fuchs, Helmerich, Bobrick, Sellers, Melcher, Martire, "Constant velocity physical warp drive solution", Class. Quantum Grav. 41, 095013
 * (2024), arXiv:2405.02709, section 3; their public code (github.com/NerdsWithAttitudes/WarpFactory): metricGet_WarpShellComoving.m,
 * TOVconstDensity.m, alphaNumericSolver.m; their example W1_Warp_Shell.mlx fixes the numbers the paper leaves out.
 *
 * UNITS: R2 = 1 (20 m), G = c = 1. Their mass M = R2 c^2/(6G) makes 2GM/c^2 = R2/3, so m_total = 1/6; the shell is R1 = 1/2 < r < 1;
 * the density rho = m_total / (4/3 pi (1 - 1/8)) = 1/(7 pi) = 0.045473; energies and pressures are in units of c^4/(G R2^2) =
 * 3.0256e41 J/m^3 (their ~1e39 J/m^3 is ~3.3e-3 here).
 *
 * THE CONSTRUCTION (theirs, in their order):
 *   1. r_i = i dr, i = 0..N-1, N = 100000, dr = rmax/(N-1), rmax = 1.2 sqrt(29.9^2 + 29.9^2 + 0.4^2) m / 20 m = 2.537212643828 (their grid).
 *   2. rho_i = 1/(7 pi) for R1 < r_i < R2 (strict, as their code), else 0.
 *   3. P_i = TOVconstDensity's closed form with the total mass, for r_i < R2:  P/rho = [sqrt(2/3) - sqrt(1 - r^2/3)] / [sqrt(1 - r^2/3) - 3 sqrt(2/3)]
 *      variant CODE : for every r_i < 1 (what metricGet_WarpShellComoving.m does: nothing zeroes it inside R1);
 *      variant PAPER: additionally P_i = 0 for r_i < R1 (the paper's step (ii): "this pressure is set to zero for r < R1 to enforce a vacuum interior").
 *   4. MATLAB smooth(y, span) = the centred moving average with span rounded down to odd and symmetric shrinking windows at the ends
 *      (y~_i = mean of y_j over |j - i| <= w_i, w_i = min(w, i, N-1-i)); FOUR passes; density span 7159 (1.79 x 4000 as the code), pressure
 *      span 3999 (4000); each mean rounded to nearest in Q.40.
 *   5. m_i = cumulative trapezoid of f = 4 pi r^2 rho~ (m_0 = 0).   6. da_i = (m_i + 4 pi r_i^3 P~_i) / (r_i^2 - 2 m_i r_i), da_0 = 0;
 *      a_i = cumulative trapezoid of da, shifted so that a_{N-1} = (1/2) ln(1 - 2 m_{N-1}/r_{N-1}) (alphaNumericSolver.m).
 *
 * THE DECLARED OBJECT: between grid points the integrands f and da are linear, so m and a are piecewise quadratic and C^1, m' = f and
 * a' = da are piecewise linear, m'' = f' and a'' = da' are piecewise constant. The metric ds^2 = -e^{2a} dt^2 + dr^2/(1 - 2m/r) + r^2 dOmega^2
 * is evaluated at every interval's midpoint, where all of these are defined; nothing is differentiated numerically.
 * Its Einstein tensor in the orthonormal frame (Carroll, Spacetime and Geometry, eq. 5.154-5.156, with e^{-2b} = 1 - 2m/r and so
 * b' = (m' r - m)/(r^2 - 2 m r)):
 *   8 pi rho_T = 2 m'/r^2 ;   8 pi P_r = -2m/r^3 + 2 a' (1 - 2m/r)/r ;   8 pi P_t = (1 - 2m/r) [a'' + a'^2 - a' b' + (a' - b')/r].
 * The four energy conditions for T = diag(rho, P_r, P_t, P_t): NEC rho + P_r >= 0, rho + P_t >= 0; WEC adds rho >= 0; SEC adds
 * rho + P_r + 2 P_t >= 0; DEC: rho - |P_r| >= 0, rho - |P_t| >= 0. Six margins; their minima over the cavity (r < 1/2), the shell
 * (1/2 <= r <= 1) and the exterior (r > 1), with where, for both variants; a table every 0.05 in r; the vacuum residuals.
 *
 * EXACT: whq.h (Q.40, 128-bit intermediates, round to nearest). Build: cc -O2 -std=gnu11 -o whshell wh_shell.c  (whq.h beside it)
 */
#include "whq.h"
#include <stdlib.h>

#define N 100000
static q rho_s[N], P_s[N], f_[N], m_[N], da_[N];
static int64_t S_[N + 1];

static void smooth4(q *y, int span)                   /* MATLAB smooth(), four passes, in place */
{
    const int w = (span % 2 ? span : span - 1) / 2;   /* span rounded down to odd; w = half window */
    for (int pass = 0; pass < 4; pass++) {
        S_[0] = 0; for (int i = 0; i < N; i++) S_[i + 1] = S_[i] + y[i];
        for (int i = 0; i < N; i++) {
            int wi = w; if (i < wi) wi = i; if (N - 1 - i < wi) wi = N - 1 - i;
            y[i] = rdiv(S_[i + wi + 1] - S_[i - wi], 2 * wi + 1);
        }
    }
}
typedef struct { q v, r; } minrec;
static void upd(minrec *m, q v, q r) { if (v < m->v) { m->v = v; m->r = r; } }

int main(void)
{
    char A[48], B[48], C[48], D[48], E[48], F[48], G[48], H[48], I[48];
    const q RHO = divq(ONE, 7 * PI);                                                         /* 1/(7 pi) */
    const q rmax = divq(6 * sqrtq(divq(178818 * ONE, 100 * ONE)), 100 * ONE);               /* 1.2 sqrt(1788.18)/20 */
    const q dr = divq(rmax, (q)(N - 1) * ONE), R1 = ONE / 2, R2 = ONE;
    const q s23 = sqrtq(divq(2 * ONE, 3 * ONE));                                             /* sqrt(2/3) */
    fmt(RHO, A); fmt(rmax, B); fmt(dr, C); fmt(s23, D);
    printf("WHSHELL: Q.%d (whq.h); R2 = 1 = 20 m, G = c = 1, m_total = 1/6; rho = 1/(7 pi) = %s; grid N = %d, rmax = %s, dr = %s; sqrt(2/3) = %s; energies in c^4/(G R2^2) = 3.0256e41 J/m^3\n", FB, A, N, B, C, D);

    /* the density, smoothed once for both variants */
    for (int i = 0; i < N; i++) { const q r = (q)i * dr; rho_s[i] = (r > R1 && r < R2) ? RHO : 0; }
    smooth4(rho_s, 7159);
    {
        int lo = -1, hi = -1; for (int i = 0; i < N; i++) if (rho_s[i] != 0) { if (lo < 0) lo = i; hi = i; }
        fmt((q)lo * dr, A); fmt((q)hi * dr, B);
        /* the smoothed mass: trapezoid of 4 pi r^2 rho~ */
        for (int i = 0; i < N; i++) { const q r = (q)i * dr; f_[i] = mulq(mulq(4 * PI, mulq(r, r)), rho_s[i]); }
        m_[0] = 0; for (int i = 1; i < N; i++) m_[i] = m_[i - 1] + rdiv(mulq(f_[i - 1] + f_[i], dr), 2);
        fmt(m_[N - 1], C); fmt(divq(ONE, 6 * ONE), D);
        printf("WHSHELL density: four passes of span 7159 (3.63 m) spread the shell's matter over r = %s .. %s (R1 = 0.5, R2 = 1); the smoothed total mass m = %s (unsmoothed 1/6 = %s)\n", A, B, C, D);
    }

    uint64_t hash = 0xcbf29ce484222325ull;
    for (int var = 0; var < 2; var++) {
        const char *name = var ? "PAPER (P = 0 inside R1)" : "CODE (P kept inside R1)";
        /* the pressure */
        for (int i = 0; i < N; i++) {
            const q r = (q)i * dr; q P = 0;
            if (r < R2 && !(var && r < R1)) {
                const q b = sqrtq(ONE - divq(mulq(r, r), 3 * ONE));                          /* sqrt(1 - r^2/3) */
                P = divq(mulq(RHO, b - s23), 3 * s23 - b);                                   /* rho (b - a)/(3a - b) = rho (a - b)/(b - 3a) */
            }
            P_s[i] = P;
        }
        smooth4(P_s, 3999);
        int plo = -1, phi = -1; for (int i = 0; i < N; i++) if (P_s[i] != 0) { if (plo < 0) plo = i; phi = i; }
        /* a' = da and a by the trapezoid, pinned at the end */
        da_[0] = 0;
        for (int i = 1; i < N; i++) {
            const q r = (q)i * dr, r2 = mulq(r, r);
            const q num = m_[i] + mulq(mulq(4 * PI, mulq(r2, r)), P_s[i]);
            const q den = r2 - mulq(2 * m_[i], r);
            da_[i] = divq(num, den);
        }
        q a = 0; for (int i = 1; i < N; i++) a += rdiv(mulq(da_[i - 1] + da_[i], dr), 2);   /* a_{N-1} before the pin */
        const q rend = (q)(N - 1) * dr;
        const q pin = rdiv(lnq(ONE - divq(2 * m_[N - 1], rend)), 2);
        const q a0 = pin - a;                                                                  /* a_0 after the pin */
        fmt((q)plo * dr, A); fmt((q)phi * dr, B); fmt(P_s[0], C); fmt(a0, D); fmt(pin, E);
        printf("WHSHELL %s: smoothed pressure nonzero over r = %s .. %s; P~(0) = %s; a(0) = %s, a(rmax) = %s (= (1/2) ln(1 - 2m/r) there)\n", name, A, B, C, D, E);

        /* the Einstein tensor at every midpoint; the margins' minima per region; the table every 0.05 */
        minrec mn[3][6]; for (int g = 0; g < 3; g++) for (int k = 0; k < 6; k++) { mn[g][k].v = (q)1 << 62; mn[g][k].r = 0; }
        q maxres_cav = 0, maxres_ext = 0, acur = a0; int lastbin = -1;
        printf("WHSHELL %s: r | rho~ | P~ | m | a | B | rho_T | P_r | P_t   (midpoints; every 0.05 in r)\n", name);
        for (int i = 0; i + 1 < N; i++) {
            const q r = (q)i * dr + dr / 2, r2 = mulq(r, r), r3 = mulq(r2, r);
            const q m = rdiv(m_[i] + m_[i + 1], 2) - rdiv(mulq(f_[i + 1] - f_[i], dr), 8);   /* the quadratic's midpoint value */
            const q mp = rdiv(f_[i] + f_[i + 1], 2);                                         /* m' */
            const q ap = rdiv(da_[i] + da_[i + 1], 2);                                       /* a' */
            const q app = divq(da_[i + 1] - da_[i], dr);                                     /* a'' */
            const q amid = acur + rdiv(mulq(da_[i], dr), 2) + rdiv(mulq(da_[i + 1] - da_[i], dr), 8);
            const q u = ONE - divq(2 * m, r);                                                /* 1 - 2m/r = 1/B */
            const q bp = divq(mulq(mp, r) - m, r2 - mulq(2 * m, r));                         /* b' */
            const q rhoT = divq(mp, mulq(4 * PI, r2));
            const q Pr = divq(-divq(2 * m, r3) + divq(mulq(2 * ap, u), r), 8 * PI);
            const q Pt = divq(mulq(u, app + mulq(ap, ap) - mulq(ap, bp) + divq(ap - bp, r)), 8 * PI);
            const q absPr = Pr < 0 ? -Pr : Pr, absPt = Pt < 0 ? -Pt : Pt;
            const q marg[6] = { rhoT + Pr, rhoT + Pt, rhoT, rhoT + Pr + 2 * Pt, rhoT - absPr, rhoT - absPt };
            const int g = r < R1 ? 0 : (r <= R2 ? 1 : 2);
            for (int k = 0; k < 6; k++) upd(&mn[g][k], marg[k], r);
            if (var && i < 5000 && absPr > maxres_cav) maxres_cav = absPr;                     /* deep cavity r < 0.127: vacuum in the PAPER variant */
            if (i > 60000 && (absPr > maxres_ext || absPt > maxres_ext)) maxres_ext = absPr > absPt ? absPr : absPt;   /* r > 1.52: Schwarzschild vacuum */
            hash = fnv(hash, &rhoT, 8); hash = fnv(hash, &Pr, 8); hash = fnv(hash, &Pt, 8);
            const int bin = (int)((r * 20) >> FB);                                           /* every 0.05 */
            if (bin != lastbin && bin <= 50) {
                lastbin = bin;
                fmt(r, A); fmt(rdiv(rho_s[i] + rho_s[i + 1], 2), B); fmt(rdiv(P_s[i] + P_s[i + 1], 2), C); fmt(m, D); fmt(amid, E);
                fmt(divq(ONE, u), F); fmt(rhoT, G); fmt(Pr, H); fmt(Pt, I);
                printf("WHSHELL %s: %s | %s | %s | %s | %s | %s | %s | %s | %s\n", name, A, B, C, D, E, F, G, H, I);
            }
            acur += rdiv(mulq(da_[i] + da_[i + 1], dr), 2);
        }
        const char *gn[3] = { "cavity r < 1/2", "shell 1/2 <= r <= 1", "exterior r > 1" };
        const char *kn[6] = { "rho + P_r (NEC)", "rho + P_t (NEC)", "rho (WEC)", "rho + P_r + 2 P_t (SEC)", "rho - |P_r| (DEC)", "rho - |P_t| (DEC)" };
        for (int g = 0; g < 3; g++) for (int k = 0; k < 6; k++) {
            fmt(mn[g][k].v, A); fmt(mn[g][k].r, B);
            printf("WHSHELL %s: min over the %s of %s = %s at r = %s %s\n", name, gn[g], kn[k], A, B, mn[g][k].v < 0 ? "VIOLATED" : "holds");
        }
        fmt(maxres_ext, A); fmt(maxres_cav, B);
        printf("WHSHELL %s: vacuum residuals: max |P_r|, |P_t| for r > 1.52 (Schwarzschild, must be 0) = %s; max |P_r| for r < 0.127 (flat cavity in PAPER, must be 0) = %s\n", name, A, B);
    }
    printf("WHSHELL: chain %016llx\n", (unsigned long long)hash);
    printf("WHSHELL: done\n");
    return 0;
}
