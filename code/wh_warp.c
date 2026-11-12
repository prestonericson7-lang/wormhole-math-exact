/* wh_warp.c (11-11-2026; v2 9:05 PM, after the math buddy's review W1-W10) -- WORMHOLE MATH 3b, part 2: the 2024 warp shell WITH its
 * shift, judged exactly; the shift limit its authors left open.
 *
 * Fuchs et al., "Constant velocity physical warp drive solution", Class. Quantum Grav. 41, 095013 (2024), arXiv:2405.02709, section 4: the
 * matter shell of section 3 (wh_prof.h reproduces its construction as integer sums) with a shift vector along x inside it, their eq. (26):
 * g_tx = -beta0 S(r), beta0 = 0.02 (their vWarp), S the compact sigmoid of eq. (27)-(28) smoothed twice; everything else unchanged. Their
 * claim: all four energy conditions hold at every point for beta0 = 0.02, "likely not an upper limit"; the limit is "a future direction
 * of work". This program evaluates the four conditions exactly on a grid in (r, theta) for a list of beta0 and brackets the limit.
 *
 * THE METRIC, with the shift direction as the polar axis (x = r cos theta), in coordinates (t, r, theta, phi):
 *   g_tt = -e^{2a(r)}, g_tr = -beta cos theta, g_t theta = beta r sin theta, g_rr = B(r) = 1/(1 - 2m/r), g_theta theta = r^2, g_phi phi = r^2 sin^2 theta,
 * beta(r) = beta0 S(r). Axisymmetric and stationary. The review confirmed this is Warp Factory's metric (to 1.2e-41) and that v1's
 * curvature code was right off the axis; v2 keeps that algebra and changes HOW it is evaluated:
 * v2 (W2): THE WARPED PRODUCT. g = h + f^2 dphi^2 with h the (t, r, theta) block and f = r sin theta. Then, exactly (O'Neill, fiber dimension 1),
 *   R_ab = R[h]_ab - (nabla_a nabla_b f) / f,   R_phi phi = -f Laplacian_h f,   R_a phi = 0,   R = h^{ab} R[h]_ab - 2 Laplacian_h f / f,
 * so g^{phi phi} and its derivatives (1/(r^2 sin^3 theta), which overflowed Q.40 next to the axis in v1) are never formed; the only
 * angular factor is 1/f. The 3D Ricci of h is the generic Gamma / d Gamma / R_mn code of v1 on three indices. The v1 four-index code is
 * kept as a REFERENCE and the two are compared at off-axis points (printed: must agree to ~1e-11).
 * v2 (W4): the shift's S, S', S'' from the exact kernel sums (wh_prof.h v2): S'' noise one unit of 2^-40, not 4e-4.
 * v2 (W3): A WRITTEN FLOOR under every verdict. The curvature code's rounding envelope, from the review's measurement against 50-digit
 *   evaluations with the axis terms removed, is 2^-40 [0.1/r^4 + 10] plus 2^-40 [2/(r sin theta)] for the 1/f terms; the tables' rounding
 *   enters through a'' = (da_{i+1} - da_i)/dr: one unit of 2^-40 over dr, times u/(8 pi): 3 x 2^-41/(8 pi dr). F = twice the first plus
 *   the second. A margin is VIOLATED only below -F; between -F and 0 it is "within the floor"; vacuum cells (matter exactly zero at both
 *   nodes of the interval) are reported apart, as residuals of the declared object, never as violations (the review: the declared
 *   exterior has an exact NEC of about -4e-10 from the a'' table, independent of beta0).
 * v2 (W10): THE BISECTION resolves the limit: the beta0 list is {0, 0.02, 0.025, 0.03, 0.035, 0.05, 0.1}; from the last pass and the
 *   first failure, 16 halvings over the CANDIDATE cells (every cell whose worst margin at the failing beta0 is below 20 F, plus a 5 x 5
 *   block around the minimum; the count is printed), the verdict at each halving = the minimum over k of (margin_k + F) over them.
 *   The DEC's |Y| uses a 128-bit norm when |Y| is small (W5 addition).
 * THE FRAME: the Eulerian observer n^mu = -alpha g^{mu t}, alpha = 1/sqrt(-g^{tt}); the triad e1 = d_r / sqrt(B), e2 = d_theta / r,
 * e3 = d_phi / f (g_r theta = 0, so they are orthonormal and normal to n). rho = T(n,n); p_a = -T(n,e_a); P_ab = T(e_a,e_b); p_3 = P_13 = P_23 = 0.
 * THE CONDITIONS over observers V = n + v m, |m| = 1 in the triad, 0 <= v <= 1 (v = 1: null): NEC: v = 1; WEC: all v; SEC: T(V,V) +
 * (1/2) T (1 - v^2); DEC: Y = -T^a_b V^b future-pointing and causal: margin Y^0 - |Y|. The minimum over m is taken on 100 Fibonacci
 * directions x 10 speeds IN THE (r, theta, phi) TRIAD with V^0 = 1 (the formulas of Warp Factory's getEvenPointsOnSphere /
 * generateUniformField, but not its Cartesian Cholesky triad nor its (1+v^2)^-1/2 normalization: the review, W8), then refined by a local
 * descent in the angles. Values in c^4/(G R2^2) = 3.0256e41 J/m^3 (their 1e39 is 3.3e-3 here).
 * THE SELF-TESTS: (1) beta0 = 0 against the shell's closed form at five radii; (2) flat space with a CONSTANT shift (m = 0, a = 0, S = 1,
 * S' = S'' = 0): every frame component must vanish to the floor; (3) the warped-product code against the v1 four-index reference at
 * off-axis points with beta0 = 0.2.
 * THE LIGHT RAYS (their section 5.1): along the axis a null ray has E (dt/dr)^2 -+ 2 beta (dt/dr) - B = 0, so the ray with the shift
 * crosses the diameter in 2 INT (S - beta)/E dr and the ray against it in 2 INT (S + beta)/E dr, S = sqrt(beta^2 + E B): the difference
 * is 4 INT beta e^{-2a} dr, printed in coordinate time and in the proper time of an emitter at rmax (times e^{a(rmax)}). The paper's
 * 7.6 ns is reproduced by neither (the review: its Table 1 is captioned 0.04 c; the same definition gives their Alcubierre 8.0 ns); stated as open.
 *
 * EXACT: whq.h v2 + wh_prof.h v2 (Q.40, 128-bit intermediates, round to nearest). ENV: WH_RSTRIDE (midpoints per r step, default 40), WH_NTH
 * (theta cells on [0, pi], default 61), WH_VARIANT (0 = HYPOTHETICAL, 1 = AS-PUBLISHED, default 1), WH_BETAS (how many of the beta0 values are used,
 * default 7), WH_HALVINGS (default 16). Build: cc -O2 -std=gnu11 -o whwarp wh_warp.c
 */
#include "wh_prof.h"
#include <stdlib.h>

typedef struct { q g[4][4], gi[4][4], d[2][4][4], dd[3][4][4]; } met_t;   /* the 4D metric (v1 layout: kept for the reference code) */

static q D1(const met_t *M, int k, int a, int b) { return k == 1 ? M->d[0][a][b] : k == 2 ? M->d[1][a][b] : 0; }
static q D2(const met_t *M, int k, int l, int a, int b)
{
    if (k == 1 && l == 1) return M->dd[0][a][b];
    if (k == 2 && l == 2) return M->dd[2][a][b];
    if ((k == 1 && l == 2) || (k == 2 && l == 1)) return M->dd[1][a][b];
    return 0;
}
/* the metric at a point: the 4D arrays (phi row/column too, used only by the reference), the (t, r, theta) inverse by the adjugate */
static void metric_at(const wh_mid *p, q beta0, q c, q s, met_t *M)
{
    const q r = p->r, r2 = mulq(r, r), r3 = mulq(r2, r);
    const q E = expq(2 * p->a), ap = p->ap, app = p->app;
    const q u = ONE - divq(2 * p->m, r);
    const q up = -divq(2 * p->mp, r) + divq(2 * p->m, r2);
    const q upp = -divq(2 * p->mpp, r) + divq(4 * p->mp, r2) - divq(4 * p->m, r3);
    const q u2 = mulq(u, u);
    const q B = divq(ONE, u), B1 = -divq(up, u2), B2 = -divq(upp, u2) + divq(2 * mulq(up, up), mulq(u, u2));
    const q b = mulq(beta0, p->S), bp = mulq(beta0, p->Sp), bpp = mulq(beta0, p->Spp);
    const q s2 = mulq(s, s), c2 = mulq(c, c);
    memset(M, 0, sizeof *M);
    M->g[0][0] = -E; M->g[0][1] = -mulq(b, c); M->g[0][2] = mulq(mulq(b, r), s); M->g[1][1] = B; M->g[2][2] = r2; M->g[3][3] = mulq(r2, s2);
    M->d[0][0][0] = -mulq(2 * ap, E); M->d[0][0][1] = -mulq(bp, c); M->d[0][0][2] = mulq(mulq(bp, r) + b, s); M->d[0][1][1] = B1; M->d[0][2][2] = 2 * r; M->d[0][3][3] = mulq(2 * r, s2);
    M->d[1][0][1] = mulq(b, s); M->d[1][0][2] = mulq(mulq(b, r), c); M->d[1][3][3] = mulq(mulq(2 * r2, s), c);
    M->dd[0][0][0] = -mulq(2 * app + mulq(4 * ap, ap), E); M->dd[0][0][1] = -mulq(bpp, c); M->dd[0][0][2] = mulq(mulq(bpp, r) + 2 * bp, s); M->dd[0][1][1] = B2; M->dd[0][2][2] = 2 * ONE; M->dd[0][3][3] = 2 * s2;
    M->dd[1][0][1] = mulq(bp, s); M->dd[1][0][2] = mulq(mulq(bp, r) + b, c); M->dd[1][3][3] = mulq(mulq(4 * r, s), c);
    M->dd[2][0][1] = mulq(b, c); M->dd[2][0][2] = -mulq(mulq(b, r), s); M->dd[2][3][3] = mulq(2 * r2, c2 - s2);
    for (int a = 0; a < 4; a++) for (int bb = a + 1; bb < 4; bb++) {
        M->g[bb][a] = M->g[a][bb];
        for (int k = 0; k < 2; k++) M->d[k][bb][a] = M->d[k][a][bb];
        for (int k = 0; k < 3; k++) M->dd[k][bb][a] = M->dd[k][a][bb];
    }
    const q g00 = M->g[0][0], g01 = M->g[0][1], g02 = M->g[0][2], g11 = M->g[1][1], g22 = M->g[2][2];
    const q det = mulq(mulq(g00, g11), g22) - mulq(mulq(g01, g01), g22) - mulq(mulq(g02, g02), g11);   /* < 0: whq.h v2's divq takes it */
    M->gi[0][0] = divq(mulq(g11, g22), det);
    M->gi[0][1] = M->gi[1][0] = -divq(mulq(g01, g22), det);
    M->gi[0][2] = M->gi[2][0] = -divq(mulq(g02, g11), det);
    M->gi[1][1] = divq(mulq(g00, g22) - mulq(g02, g02), det);
    M->gi[1][2] = M->gi[2][1] = divq(mulq(g01, g02), det);
    M->gi[2][2] = divq(mulq(g00, g11) - mulq(g01, g01), det);
    M->gi[3][3] = M->g[3][3] ? divq(ONE, M->g[3][3]) : 0;                   /* only the reference uses it */
}
/* ---- v2: the warped-product Einstein tensor. T[a][b] for a, b in (t, r, theta) covariant; T33 = T_phi phi / g_phi phi (the e3 e3 frame value) ---- */
static void einstein_wp(const met_t *M, q r, q c, q s, q T[3][3], q *T33)
{
    q Gm[3][3][3], dGm[2][3][3][3], dgi[2][3][3];
    for (int k = 0; k < 2; k++) for (int l = 0; l < 3; l++) for (int ss = 0; ss < 3; ss++) {
        q acc = 0;
        for (int a = 0; a < 3; a++) { if (!M->gi[l][a]) continue; for (int b = 0; b < 3; b++) { const q dab = D1(M, k + 1, a, b); if (dab && M->gi[b][ss]) acc -= mulq(mulq(M->gi[l][a], dab), M->gi[b][ss]); } }
        dgi[k][l][ss] = acc;
    }
    for (int l = 0; l < 3; l++) for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++) {
        q acc = 0;
        for (int ss = 0; ss < 3; ss++) { if (!M->gi[l][ss]) continue; const q w = D1(M, m, ss, n) + D1(M, n, ss, m) - D1(M, ss, m, n); if (w) acc += mulq(M->gi[l][ss], w); }
        Gm[l][m][n] = rdiv(acc, 2);
        for (int k = 0; k < 2; k++) {
            q acc2 = 0;
            for (int ss = 0; ss < 3; ss++) {
                const q w = D1(M, m, ss, n) + D1(M, n, ss, m) - D1(M, ss, m, n);
                if (w && dgi[k][l][ss]) acc2 += mulq(dgi[k][l][ss], w);
                if (M->gi[l][ss]) { const q w2 = D2(M, k + 1, m, ss, n) + D2(M, k + 1, n, ss, m) - D2(M, k + 1, ss, m, n); if (w2) acc2 += mulq(M->gi[l][ss], w2); }
            }
            dGm[k][l][m][n] = rdiv(acc2, 2);
        }
    }
    q Rh[3][3];                                                           /* the 3D Ricci of h */
    for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++) {
        q acc = 0;
        for (int l = 1; l <= 2; l++) acc += dGm[l - 1][l][m][n];
        if (n == 1 || n == 2) for (int l = 0; l < 3; l++) acc -= dGm[n - 1][l][m][l];
        for (int l = 0; l < 3; l++) for (int ss = 0; ss < 3; ss++) {
            if (Gm[l][l][ss] && Gm[ss][m][n]) acc += mulq(Gm[l][l][ss], Gm[ss][m][n]);
            if (Gm[l][n][ss] && Gm[ss][m][l]) acc -= mulq(Gm[l][n][ss], Gm[ss][m][l]);
        }
        Rh[m][n] = acc;
    }
    /* the Hessian of f = r sin theta on h: H_ab = d_a d_b f - Gamma^c_ab d_c f, with d_r f = s, d_theta f = r c, d_rtheta f = c, d_thetatheta f = -r s */
    const q df[3] = {0, s, mulq(r, c)};
    q H[3][3];
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
        q v = 0;
        if ((a == 1 && b == 2) || (a == 2 && b == 1)) v = c; else if (a == 2 && b == 2) v = -mulq(r, s);
        for (int cc = 1; cc < 3; cc++) if (Gm[cc][a][b] && df[cc]) v -= mulq(Gm[cc][a][b], df[cc]);
        H[a][b] = v;
    }
    const q f = mulq(r, s);
    q lap = 0; for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) if (M->gi[a][b] && H[a][b]) lap += mulq(M->gi[a][b], H[a][b]);
    const q lapf = divq(lap, f);                                           /* Laplacian f / f: the one 1/(r sin theta) in the whole computation */
    q R4[3][3], Rs = 0;
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) { R4[a][b] = Rh[a][b] - divq(H[a][b], f); }
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) if (M->gi[a][b] && Rh[a][b]) Rs += mulq(M->gi[a][b], Rh[a][b]);
    Rs -= 2 * lapf;                                                        /* R = h^{ab} R[h]_ab - 2 lap f / f */
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) T[a][b] = divq(R4[a][b] - rdiv(mulq(M->g[a][b], Rs), 2), 8 * PI);
    *T33 = divq(-lapf - rdiv(Rs, 2), 8 * PI);                            /* T_phi phi / g_phi phi = (R_phi phi / f^2 - R/2) / 8 pi */
}
/* ---- the v1 four-index reference (unchanged algebra; overflows next to the axis, used off-axis only for the comparison) ---- */
static void einstein_ref(const met_t *M, q T[4][4])
{
    static q Gm[4][4][4], dGm[2][4][4][4], dgi[2][4][4];
    for (int k = 0; k < 2; k++) for (int l = 0; l < 4; l++) for (int ss = 0; ss < 4; ss++) {
        q acc = 0;
        for (int a = 0; a < 4; a++) { if (!M->gi[l][a]) continue; for (int b = 0; b < 4; b++) { const q dab = D1(M, k + 1, a, b); if (dab && M->gi[b][ss]) acc -= mulq(mulq(M->gi[l][a], dab), M->gi[b][ss]); } }
        dgi[k][l][ss] = acc;
    }
    for (int l = 0; l < 4; l++) for (int m = 0; m < 4; m++) for (int n = 0; n < 4; n++) {
        q acc = 0;
        for (int ss = 0; ss < 4; ss++) { if (!M->gi[l][ss]) continue; const q w = D1(M, m, ss, n) + D1(M, n, ss, m) - D1(M, ss, m, n); if (w) acc += mulq(M->gi[l][ss], w); }
        Gm[l][m][n] = rdiv(acc, 2);
        for (int k = 0; k < 2; k++) {
            q acc2 = 0;
            for (int ss = 0; ss < 4; ss++) {
                const q w = D1(M, m, ss, n) + D1(M, n, ss, m) - D1(M, ss, m, n);
                if (w && dgi[k][l][ss]) acc2 += mulq(dgi[k][l][ss], w);
                if (M->gi[l][ss]) { const q w2 = D2(M, k + 1, m, ss, n) + D2(M, k + 1, n, ss, m) - D2(M, k + 1, ss, m, n); if (w2) acc2 += mulq(M->gi[l][ss], w2); }
            }
            dGm[k][l][m][n] = rdiv(acc2, 2);
        }
    }
    q R[4][4];
    for (int m = 0; m < 4; m++) for (int n = 0; n < 4; n++) {
        q acc = 0;
        for (int l = 1; l <= 2; l++) acc += dGm[l - 1][l][m][n];
        if (n == 1 || n == 2) for (int l = 0; l < 4; l++) acc -= dGm[n - 1][l][m][l];
        for (int l = 0; l < 4; l++) for (int ss = 0; ss < 4; ss++) {
            if (Gm[l][l][ss] && Gm[ss][m][n]) acc += mulq(Gm[l][l][ss], Gm[ss][m][n]);
            if (Gm[l][n][ss] && Gm[ss][m][l]) acc -= mulq(Gm[l][n][ss], Gm[ss][m][l]);
        }
        R[m][n] = acc;
    }
    q Rs = 0;
    for (int m = 0; m < 4; m++) for (int n = 0; n < 4; n++) if (M->gi[m][n] && R[m][n]) Rs += mulq(M->gi[m][n], R[m][n]);
    for (int m = 0; m < 4; m++) for (int n = 0; n < 4; n++) T[m][n] = divq(R[m][n] - rdiv(mulq(M->g[m][n], Rs), 2), 8 * PI);
}
typedef struct { q rho, p1, p2, P11, P12, P22, P33; } frame_t;
static void frame_wp(const met_t *M, const q T[3][3], q T33, q r, frame_t *F)
{
    const q alpha = divq(ONE, sqrtq(-M->gi[0][0]));
    q n[3]; for (int m = 0; m < 3; m++) n[m] = -mulq(alpha, M->gi[m][0]);
    const q sB = sqrtq(M->g[1][1]);
    q rho = 0, t1 = 0, t2 = 0;
    for (int m = 0; m < 3; m++) { if (!n[m]) continue; for (int k = 0; k < 3; k++) if (n[k] && T[m][k]) rho += mulq(mulq(n[m], n[k]), T[m][k]); t1 += mulq(n[m], T[m][1]); t2 += mulq(n[m], T[m][2]); }
    F->rho = rho; F->p1 = -divq(t1, sB); F->p2 = -divq(t2, r);
    F->P11 = divq(T[1][1], M->g[1][1]); F->P12 = divq(T[1][2], mulq(r, sB)); F->P22 = divq(T[2][2], mulq(r, r)); F->P33 = T33;
}
static void frame_ref(const met_t *M, const q T[4][4], q r, q s, frame_t *F)
{
    const q alpha = divq(ONE, sqrtq(-M->gi[0][0]));
    q n[4]; for (int m = 0; m < 4; m++) n[m] = -mulq(alpha, M->gi[m][0]);
    const q sB = sqrtq(M->g[1][1]);
    q rho = 0, t1 = 0, t2 = 0;
    for (int m = 0; m < 4; m++) { if (!n[m]) continue; for (int k = 0; k < 4; k++) if (n[k] && T[m][k]) rho += mulq(mulq(n[m], n[k]), T[m][k]); t1 += mulq(n[m], T[m][1]); t2 += mulq(n[m], T[m][2]); }
    F->rho = rho; F->p1 = -divq(t1, sB); F->p2 = -divq(t2, r);
    F->P11 = divq(T[1][1], M->g[1][1]); F->P12 = divq(T[1][2], mulq(r, sB)); F->P22 = divq(T[2][2], mulq(r, r)); F->P33 = divq(T[3][3], M->g[3][3]);
    (void)s;
}
/* |Y| with a 128-bit norm when Y is small (W5 addition): Y in Q.40; returns Q.40 */
static q norm3(q Y1, q Y2, q Y3)
{
    const q sq = mulq(Y1, Y1) + mulq(Y2, Y2) + mulq(Y3, Y3);
    if (sq >= (ONE >> 20)) return sqrtq(sq);                               /* |Y| >= 1e-3: the Q.40 route is exact enough (error below 2^-41 / 2|Y|) */
    uint64_t h, l; u128 acc = {0, 0};
    const q ys[3] = {Y1, Y2, Y3};
    for (int i = 0; i < 3; i++) { const uint64_t a = ys[i] < 0 ? -(uint64_t)ys[i] : (uint64_t)ys[i]; mul128(a, a, &h, &l); u128 t = {h, l}; acc = u128_add(acc, t); }
    /* acc = |Y|^2 in Q.80; its square root in Q.40 by 64 bit-setting steps with a 128-bit square */
    uint64_t rr = 0;
    for (int bit = 63; bit >= 0; bit--) {
        const uint64_t cand = rr | ((uint64_t)1 << bit);
        if (bit >= 44) continue;                                           /* |Y| < 2^3 here: cand^2 < 2^88 */
        mul128(cand, cand, &h, &l); u128 c2 = {h, l};
        if (!u128_lt(acc, c2)) rr = cand;
    }
    return (q)rr;
}
/* the four margins for a direction m (unit, triad) and a speed v */
static void margins(const frame_t *F, const q m[3], q v, q out[4])
{
    const q pm = mulq(F->p1, m[0]) + mulq(F->p2, m[1]);
    const q Pm0 = mulq(F->P11, m[0]) + mulq(F->P12, m[1]), Pm1 = mulq(F->P12, m[0]) + mulq(F->P22, m[1]), Pm2 = mulq(F->P33, m[2]);
    const q mPm = mulq(m[0], Pm0) + mulq(m[1], Pm1) + mulq(m[2], Pm2);
    const q v2 = mulq(v, v);
    const q TVV = F->rho - 2 * mulq(v, pm) + mulq(v2, mPm);
    const q tr = -F->rho + F->P11 + F->P22 + F->P33;
    out[0] = F->rho - 2 * pm + mPm;                                   /* NEC (v = 1) */
    out[1] = TVV;                                                     /* WEC */
    out[2] = TVV + rdiv(mulq(tr, ONE - v2), 2);                       /* SEC */
    const q Y0 = F->rho - mulq(v, pm), Y1 = F->p1 - mulq(v, Pm0), Y2 = F->p2 - mulq(v, Pm1), Y3 = -mulq(v, Pm2);
    out[3] = Y0 - norm3(Y1, Y2, Y3);                                  /* DEC */
}
static q GOLD;
static void fib_dir(int i, int n, q m[3])
{
    const q th = divq(mulq(2 * PI, (q)i * ONE), GOLD);
    const q cph = ONE - divq((2 * (q)i + 1) * ONE, (q)n * ONE), sph = sqrtq(ONE - mulq(cph, cph));
    m[0] = mulq(cosq(th), sph); m[1] = mulq(sinq(th), sph); m[2] = cph;
}
static void angles_dir(q th, q ph, q m[3]) { const q sph = sinq(ph); m[0] = mulq(cosq(th), sph); m[1] = mulq(sinq(th), sph); m[2] = cosq(ph); }
static void point_minima(const frame_t *F, q coarse[4], q fine[4])
{
    const int ND = 100, NV = 10;
    q bestth[4] = {0, 0, 0, 0}, bestph[4] = {0, 0, 0, 0}, bestv[4] = {ONE, ONE, ONE, ONE};
    for (int k = 0; k < 4; k++) coarse[k] = (q)1 << 62;
    for (int i = 0; i < ND; i++) {
        q m[3]; fib_dir(i, ND, m);
        const q th = divq(mulq(2 * PI, (q)i * ONE), GOLD), cph = ONE - divq((2 * (q)i + 1) * ONE, (q)ND * ONE);
        q ph = PI / 2; for (int it = 0; it < 30; it++) { const q cs = cosq(ph), sn = sinq(ph); if (sn == 0) break; ph += divq(cs - cph, sn); }
        for (int j = 0; j < NV; j++) {
            const q v = ONE - divq((q)j * ONE, 9 * ONE);
            q o[4]; margins(F, m, v, o);
            for (int k = 0; k < 4; k++) if (k == 0 ? j == 0 : 1) if (o[k] < coarse[k]) { coarse[k] = o[k]; bestth[k] = th; bestph[k] = ph; bestv[k] = v; }
        }
    }
    for (int k = 0; k < 4; k++) {
        q th = bestth[k], ph = bestph[k], v = bestv[k], best = coarse[k];
        q step = divq(PI, 10 * ONE), vstep = divq(ONE, 9 * ONE);
        for (int round = 0; round < 16; round++) {
            int moved = 0;
            for (int dth = -1; dth <= 1; dth++) for (int dph = -1; dph <= 1; dph++) for (int dv = -1; dv <= 1; dv++) {
                if (k == 0 && dv) continue;
                if (!dth && !dph && !dv) continue;
                const q t2 = th + dth * step, p2 = ph + dph * step; q v2 = v + dv * vstep;
                if (v2 < 0) v2 = 0; if (v2 > ONE) v2 = ONE;
                q m[3]; angles_dir(t2, p2, m); q o[4]; margins(F, m, v2, o);
                if (o[k] < best) { best = o[k]; th = t2; ph = p2; v = v2; moved = 1; }
            }
            if (!moved) { step = rdiv(step, 2); vstep = rdiv(vstep, 2); }
        }
        fine[k] = best;
    }
}
static q beta_list(int i) { static const int num[7] = {0, 20, 25, 30, 35, 50, 100}; return divq((q)num[i] * ONE, 1000 * ONE); }
/* THE FLOOR at a cell */
static q floor_at(q r, q s)
{
    const q r2 = mulq(r, r), r4 = mulq(r2, r2);
    const q curv = divq(ONE, 10 * r4) + 10 * ONE + divq(2 * ONE, mulq(r, s));     /* 2^-40 [0.1/r^4 + 10 + 2/(r sin theta)] in units of 2^-40 (a Q.40 count) */
    const q table = divq(3 * ONE, mulq(16 * PI, DR));                            /* 3 x 2^-41 / (8 pi dr) = 3/(16 pi dr) units */
    return 3 * (curv >> FB) + (table >> FB) + 2;                                  /* in Q.40 units (2^-40 each): three times the review's envelope (one vacuum cell at r = 0.1 sat at 1.35 x it with the factor 2), plus the table term, +2 for the margin's own roundings */
}
/* one cell: the frame and the minima, with the floor */
typedef struct { frame_t F; q fine[4], coarse[4], floor; int vac; } cell_t;
static void eval_cell(int i, q beta0, q c, q s, cell_t *o)
{
    wh_mid p; prof_mid(i, &p); met_t M; metric_at(&p, beta0, c, s, &M); q T[3][3], T33; einstein_wp(&M, p.r, c, s, T, &T33); frame_wp(&M, T, T33, p.r, &o->F);
    point_minima(&o->F, o->coarse, o->fine);
    o->floor = floor_at(p.r, s);
    o->vac = (f_[i] == 0 && f_[i + 1] == 0 && P_s[i] == 0 && P_s[i + 1] == 0);   /* the declared matter is exactly zero on this interval */
}

int main(void)
{
    char A[48], B[48], C[48], D[48], E[48], F6[48];
    const char *ev;
    const int rstride = (ev = getenv("WH_RSTRIDE")) ? atoi(ev) : 40, nth = (ev = getenv("WH_NTH")) ? atoi(ev) : 61;
    const int variant = (ev = getenv("WH_VARIANT")) ? atoi(ev) : 1, nb = (ev = getenv("WH_BETAS")) ? atoi(ev) : 7, nhalf = (ev = getenv("WH_HALVINGS")) ? atoi(ev) : 16;
    GOLD = rdiv(ONE + sqrtq(5 * ONE), 2);
    prof_build(variant);
    fmt(RMAX, A); fmt(DR, B); fmt(SHDEV, C); fmt(MT, D);
    printf("WHWARP v2: Q.%d (whq.h v2 + wh_prof.h v2); the 2024 shell with its shift, variant %s; grid r from 0.1 to 1.5 every %d midpoints (dr = %s), %d theta cells on [0, pi]; rmax %s; the step mass MT = %s (their M(end)); the twice-smoothed shift differs from the raw sigmoid by at most %s\n",
           FB, variant ? "AS-PUBLISHED (P = 0 inside R1)" : "HYPOTHETICAL (P kept inside R1)", rstride, B, nth, A, D, C);
    const int i0 = (int)((ONE / 10) / DR), i1 = (int)((3 * ONE / 2) / DR);
    uint64_t hash = 0xcbf29ce484222325ull;
    /* self-test 1: beta0 = 0 against the shell's closed form at five radii, theta = pi/3 */
    {
        const q th = divq(PI, 3 * ONE), c = cosq(th), s = sinq(th);
        printf("WHWARP selftest 1 (beta0 = 0, theta = pi/3): r | rho generic | rho closed | P_r generic | P_r closed | P_t generic | P_t closed | p1 p2 P12 (must be 0) | floor\n");
        for (int t = 0; t < 5; t++) {
            const int i = i0 + (i1 - i0) * t / 4; wh_mid p; prof_mid(i, &p);
            met_t M; metric_at(&p, 0, c, s, &M); q T[3][3], T33; einstein_wp(&M, p.r, c, s, T, &T33); frame_t Fr; frame_wp(&M, T, T33, p.r, &Fr);
            const q r = p.r, r2 = mulq(r, r), r3 = mulq(r2, r), u = ONE - divq(2 * p.m, r);
            const q bprime = divq(mulq(p.mp, r) - p.m, r2 - mulq(2 * p.m, r));
            const q rhoc = divq(p.mp, mulq(4 * PI, r2));
            const q Prc = divq(-divq(2 * p.m, r3) + divq(mulq(2 * p.ap, u), r), 8 * PI);
            const q Ptc = divq(mulq(u, p.app + mulq(p.ap, p.ap) - mulq(p.ap, bprime) + divq(p.ap - bprime, r)), 8 * PI);
            fmt(r, A); fmt(Fr.rho, B); fmt(rhoc, C); fmt(Fr.P11, D); fmt(Prc, E); fmt(Fr.P22, F6);
            char G1[48], G2[48], G3[48], G4[48], G5[48]; fmt(Ptc, G1); fmt(Fr.p1, G2); fmt(Fr.p2, G3); fmt(Fr.P12, G4); fmt(floor_at(r, s), G5);
            printf("WHWARP selftest 1: %s | %s | %s | %s | %s | %s | %s | %s %s %s | %s\n", A, B, C, D, E, F6, G1, G2, G3, G4, G5);
            hash = fnv(hash, &Fr.rho, 8); hash = fnv(hash, &Fr.P22, 8);
        }
    }
    /* self-test 2: flat space with a constant shift: m = a = 0, S = 1, S' = S'' = 0: every frame component must be 0 to the floor */
    {
        q worst = 0;
        for (int t = 0; t < 3; t++) for (int bi = 0; bi < 2; bi++) for (int tj = 0; tj < 3; tj++) {
            wh_mid p; memset(&p, 0, sizeof p); p.r = (t + 1) * (3 * ONE / 10); p.S = ONE;        /* r = 0.3, 0.6, 0.9 */
            const q beta0 = bi ? 3 * ONE / 10 : ONE / 10;
            const q th = divq(mulq(PI, (2 * (q)tj + 1) * ONE), 6 * ONE), c = cosq(th), s = sinq(th);   /* pi/6, pi/2, 5pi/6 */
            met_t M; metric_at(&p, beta0, c, s, &M); q T[3][3], T33; einstein_wp(&M, p.r, c, s, T, &T33); frame_t Fr; frame_wp(&M, T, T33, p.r, &Fr);
            const q vals[7] = {Fr.rho, Fr.p1, Fr.p2, Fr.P11, Fr.P12, Fr.P22, Fr.P33};
            for (int k = 0; k < 7; k++) { q v = vals[k] < 0 ? -vals[k] : vals[k]; if (v > worst) worst = v; }
        }
        fmt(worst, A); printf("WHWARP selftest 2 (flat space with a constant shift, beta0 = 0.1 and 0.3, three radii, three angles): max |frame component| = %s (must be 0 to the floor, about 1e-11)\n", A);
        hash = fnv(hash, &worst, 8);
    }
    /* self-test 3: the warped-product code against the v1 four-index reference at off-axis points, beta0 = 0.2 */
    {
        q worst = 0; const q beta0 = ONE / 5;
        for (int t = 0; t < 5; t++) for (int tj = 0; tj < 2; tj++) {
            const int i = i0 + (i1 - i0) * t / 4; wh_mid p; prof_mid(i, &p);
            const q th = tj ? divq(PI, 3 * ONE) : divq(2 * PI, 3 * ONE), c = cosq(th), s = sinq(th);
            met_t M; metric_at(&p, beta0, c, s, &M);
            q T[3][3], T33; einstein_wp(&M, p.r, c, s, T, &T33); frame_t Fa; frame_wp(&M, T, T33, p.r, &Fa);
            q T4[4][4]; einstein_ref(&M, T4); frame_t Fb; frame_ref(&M, T4, p.r, s, &Fb);
            const q da[7] = {Fa.rho - Fb.rho, Fa.p1 - Fb.p1, Fa.p2 - Fb.p2, Fa.P11 - Fb.P11, Fa.P12 - Fb.P12, Fa.P22 - Fb.P22, Fa.P33 - Fb.P33};
            for (int k = 0; k < 7; k++) { q v = da[k] < 0 ? -da[k] : da[k]; if (v > worst) worst = v; }
        }
        fmt(worst, A); printf("WHWARP selftest 3 (warped product vs the four-index reference, beta0 = 0.2, five radii x two angles): max |difference| = %s (must be ~1e-11)\n", A);
        hash = fnv(hash, &worst, 8);
    }
    /* the light rays along the axis */
    for (int bi = 1; bi < nb; bi++) {
        const q beta0 = beta_list(bi);
        q diff = 0, twith = 0, tagainst = 0, flat = 0;
        for (int i = 0; i + 1 < N; i++) {
            wh_mid p; prof_mid(i, &p);
            const q E2 = expq(2 * p.a), Bv = divq(ONE, ONE - divq(2 * p.m, p.r)), b = mulq(beta0, p.S);
            const q S = sqrtq(mulq(b, b) + mulq(E2, Bv));
            diff += mulq(4 * divq(b, E2), DR); twith += mulq(2 * divq(S - b, E2), DR); tagainst += mulq(2 * divq(S + b, E2), DR); flat += 2 * DR;
        }
        const q ns = 73351520255790;                                       /* 20 m / c in ns, Q.40 */
        const q ea = expq(a_[N - 1]);                                      /* the emitter's lapse at rmax: proper time = e^a x coordinate time */
        fmt(beta0, A); fmt(diff, B); fmt(mulq(diff, ns), C); fmt(mulq(mulq(diff, ns), ea), D); fmt(twith - flat, E); fmt(tagainst - flat, F6);
        printf("WHWARP light rays, beta0 = %s: round-trip difference 4 INT beta e^-2a dr = %s R2 = %s ns coordinate time, %s ns in the proper time of an emitter at rmax (R2 = 20 m); the ray with the shift crosses the diameter %s R2 later than flat space, the ray against it %s R2 later (delays). The paper's 7.6 ns matches neither definition (open).\n", A, B, C, D, E, F6);
        hash = fnv(hash, &diff, 8);
    }
    /* THE SCAN over beta0 (from 0): the four minima over the grid with the floor; vacuum cells apart */
    const char *kn[4] = {"NEC", "WEC", "SEC", "DEC"};
    const int nr = (i1 - i0) / rstride + 1;
    cell_t *cells = malloc(sizeof(cell_t) * (size_t)nr * nth);
    q firstfail = -1, lastpass = -1; int ifail = 0, jfail = 0, kfail = -1;
    for (int bi = 0; bi < nb; bi++) {
        const q beta0 = beta_list(bi);
        q gmin[4], gr[4], gth[4], gfl[4]; int gi_[4] = {0, 0, 0, 0}, gj_[4] = {0, 0, 0, 0}; for (int k = 0; k < 4; k++) { gmin[k] = (q)1 << 62; gr[k] = gth[k] = gfl[k] = 0; }
        q vacres = 0, vacfl = 0; int nvac = 0, nviol[4] = {0, 0, 0, 0};
        for (int j = 0; j < nth; j++) {
            const q th = divq(mulq(PI, (2 * (q)j + 1) * ONE), 2 * (q)nth * ONE), c = cosq(th), s = sinq(th);
            for (int ii = 0; ii < nr; ii++) {
                const int i = i0 + ii * rstride; cell_t *o = &cells[(size_t)ii * nth + j]; eval_cell(i, beta0, c, s, o);
                wh_mid p; prof_mid(i, &p);
                for (int k = 0; k < 4; k++) hash = fnv(hash, &o->fine[k], 8);
                if (o->vac) { nvac++; for (int k = 0; k < 4; k++) { q v = o->fine[k] < 0 ? -o->fine[k] : o->fine[k]; if (v > vacres) { vacres = v; vacfl = o->floor; } } continue; }
                for (int k = 0; k < 4; k++) { if (o->fine[k] < gmin[k]) { gmin[k] = o->fine[k]; gr[k] = p.r; gth[k] = th; gfl[k] = o->floor; gi_[k] = ii; gj_[k] = j; } if (o->fine[k] < -o->floor) nviol[k]++; }
            }
        }
        fmt(beta0, A); fmt(vacres, B); fmt(vacfl, C);
        printf("WHWARP beta0 = %s: %d vacuum cells (the declared matter exactly zero): max |residual| = %s against their floor %s (the declared object's own, not a violation)\n", A, nvac, B, C);
        int fail = 0;
        for (int k = 0; k < 4; k++) {
            fmt(gmin[k], B); fmt(gfl[k], C); fmt(gr[k], D); fmt(gth[k], E);
            const char *verdict = gmin[k] < -gfl[k] ? "VIOLATED" : gmin[k] < 0 ? "within the floor" : "holds";
            printf("WHWARP beta0 = %s: min %s over the matter cells = %s (floor %s) at r = %s, theta = %s: %s (%d cells below the floor)\n", A, kn[k], B, C, D, E, verdict, nviol[k]);
            if (gmin[k] < -gfl[k] && !fail) { fail = 1; if (firstfail < 0) { firstfail = beta0; ifail = gi_[k]; jfail = gj_[k]; kfail = k; } }
        }
        if (!fail) lastpass = beta0;
        if (firstfail >= 0) break;
    }
    if (firstfail < 0) { fmt(lastpass, A); printf("WHWARP shift limit: every condition holds on the grid up to beta0 = %s (the last value tried); extend WH_BETAS\n", A); }
    else {
        /* the candidates: cells whose worst (margin + floor) at the failing beta0 is below 20 floors, plus the 5 x 5 block around the first failure */
        static int cand_i[200000], cand_j[200000]; int nc = 0;
        for (int j = 0; j < nth; j++) for (int ii = 0; ii < nr; ii++) {
            const cell_t *o = &cells[(size_t)ii * nth + j]; if (o->vac) continue;
            q worst = (q)1 << 62; for (int k = 0; k < 4; k++) if (o->fine[k] < worst) worst = o->fine[k];
            const int near = abs(ii - ifail) <= 2 && abs(j - jfail) <= 2;
            if (worst < 20 * o->floor || near) { if (nc < 200000) { cand_i[nc] = ii; cand_j[nc] = j; nc++; } }
        }
        fmt(lastpass, A); fmt(firstfail, B);
        printf("WHWARP bisection: between beta0 = %s (all hold) and %s (%s first fails); %d candidate cells of %d; %d halvings; the verdict at each = min over the candidates and the four conditions of (margin + floor)\n", A, B, kn[kfail], nc, nr * nth, nhalf);
        q lo = lastpass, hi = firstfail; int klim = kfail; q rlim = 0, thlim = 0;
        for (int step = 0; step < nhalf; step++) {
            const q mid = rdiv(lo + hi, 2); q worst = (q)1 << 62; int kw = 0; q rw = 0, thw = 0;
            for (int cidx = 0; cidx < nc; cidx++) {
                const int ii = cand_i[cidx], j = cand_j[cidx], i = i0 + ii * rstride;
                const q th = divq(mulq(PI, (2 * (q)j + 1) * ONE), 2 * (q)nth * ONE), c = cosq(th), s = sinq(th);
                cell_t o; eval_cell(i, mid, c, s, &o); wh_mid p; prof_mid(i, &p);
                for (int k = 0; k < 4; k++) { const q v = o.fine[k] + o.floor; if (v < worst) { worst = v; kw = k; rw = p.r; thw = th; } }
            }
            if (worst < 0) { hi = mid; klim = kw; rlim = rw; thlim = thw; } else lo = mid;
        }
        fmt(lo, A); fmt(hi, B); fmt(rlim, C); fmt(thlim, D);
        printf("WHWARP shift limit: %s first fails (below its floor) between beta0 = %s and %s, at r = %s, theta = %s; this bracket, the r grid and the floor bound it (the review's resolution: about 1e-6 in beta0 each)\n", kn[klim], A, B, C, D);
        hash = fnv(hash, &lo, 8); hash = fnv(hash, &hi, 8);
    }
    free(cells);
    printf("WHWARP: chain %016llx\n", (unsigned long long)hash);
    printf("WHWARP: done\n");
    return 0;
}
