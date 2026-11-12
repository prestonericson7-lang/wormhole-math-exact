/* wh_syk3.c (11-11-2026, 9:40 PM) -- WORMHOLE MATH 2, v3: the SYK traversable-wormhole teleportation protocol (Gao-Jafferis, arXiv:1911.07416)
   EXACT AT FINITE N FOR ONE REALIZATION OF THE COUPLINGS: the full reduced density matrix rho_TR of the reference qubit R (half of an EPR
   pair, untouched after the swap) and the right Dirac mode T (the readout; the paper's letters) from its six independent matrix elements (the paper's (2.15)-(2.20)). v1 and v2 printed the paper's (2.63)-(2.66), which
   are the ENSEMBLE-AVERAGED forms: they set rho_22 = rho_33 = 1 - rho_11 and rho_23 = 0 by the SO(N) symmetry of the average, and hold
   for one seed only approximately (an independent review found I_RT from them off by -12.9% to +5.2% at N = 8). v3 prints both: the exact
   matrix, and with formulas=1 (the default) the averaged forms on the same seed, so the two stand side by side on one line.

   THE MATRIX. The message qubit is the left Dirac mode chi_l = (psi_l^1 + i psi_l^2)/sqrt2 at time -t', entangled with the reference R
   (R is half of an EPR pair whose other half was swapped into the mode); the readout, the paper's T, is the right mode
   chi_r = (psi_r^1 + i psi_r^2)/sqrt2 at time t, after the coupling U = e^{i mu V}. With M_ab = |a><b| on the left mode (M_00 = chi chi^dag, M_11 = chi^dag chi, M_10 = chi^dag,
   M_01 = chi) and N_cc' = |c><c'| on the right mode (the same with chi_r), the elements in the basis |a c> (a = R the reference, c = T the readout; index k = 2a + c)
   are
       rho_(ac),(a'c') = (1/2) SUM_b < M_ba'(-t') U^dag N_c'c(t) U M_ab(-t') >_TFD,
   nonzero only for a + a' + c + c' even (fermion parity): the four diagonal elements, rho_03 and rho_12 (six, the paper's (2.15)-(2.20)).
   It is Hermitian with trace 1 by construction, and rho_00 + rho_11 = rho_22 + rho_33 = 1/2 exactly (the reference's marginal: R is half of an
   EPR pair and nothing touches it after the swap; the readout's marginal rho_00 + rho_22, rho_11 + rho_33 is NOT 1/2 for one seed when
   mu != 0, as the review checked), which is printed as a check. I_RT = S(R) + S(T) - S(RT): S(RT) from the two 2 x 2 blocks {0,3} and
   {1,2}, S(R) from (rho_00 + rho_11, rho_22 + rho_33) = (1/2, 1/2), S(T) from (rho_00 + rho_22, rho_11 + rho_33). The averaged forms are the special case rho_00 = rho_33,
   rho_11 = rho_22, rho_12 = 0, in which rho_11(paper) = 2 rho_00 and rho_14(paper) = 2 rho_03.

   THE COMPUTATION, in the factorized layout of v2 (left Majoranas on qubits 0..N/2-1, right on N/2..N-1; Psi a D x D matrix, D = 2^{N/2};
   e^{-iHt} Psi = U_L Psi U_R^T). With the bare Majoranas M_j = sqrt2 psi_l^j and R_j = sqrt2 psi_r^j (Pauli strings): 2 chi = M_1 + i M_2,
   2 chi^dag = M_1 - i M_2, 2 chi chi^dag = 1 - i M_1 M_2, 2 chi^dag chi = 1 + i M_1 M_2, and the same on the right, so every 2M and 2N is
   the sum of two Pauli strings with phases, applied exactly (no rounding). The four states Z_ab = e^{-iHt} U e^{-iHt'} (2 M_ab) e^{iHt'} |TFD>
   cost 5t' + 4t evolution steps together (v2's formulas cost 16t); each element is (1/16) SUM_b <Z_a'b | 2 N_c'c | Z_ab>, its twelve
   inner-product sums accumulated in signed 128-bit and rounded ONCE per element after the exact division by 16. K = <{A1, B1}>, the
   paper's large-N variable, is printed as before.

   Everything else is v2's, bit for bit: whq.h (Q.40 with 128-bit intermediates), the exact integer Hamiltonians, the Taylor matrix
   exponentials, the couplings J ~ N(0, 6/N^3) from xorshift64* with the declared seed (so v1, v2 and v3 share them), |I> from the proven
   phase rule e(L) = SUM_k (k even ? 3 : 1) bit_k + 2 C(|L|, 2) mod 4 (derived by the review 9:26 PM; verified exactly by the self-test).
   Command line as v2, plus formulas=0|1:
       whsyk3 N=14 seeds=1 beta=4 mu=-12,12,0 t=1,2,3,4,5,6 [tp=...] [pair=1,2] [vnorm=N] [formulas=1]
   v3 takes INTEGER t and t' (U(1)^t). Self-tests printed: |I>'s annihilation, exactly; H_L Psi_I = Psi_I H_R^T exactly; the two EPR
   identities rho_00 + rho_11 = rho_22 + rho_33 = 1/2 and the imaginary parts of the diagonal (all to the rounding); mu = 0 -> I_RT = 0. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "whq.h"

#define MAXN 16
#define MAXD (1 << (MAXN / 2))
#define MAXTERMS 1820                       /* C(16, 4) */

typedef struct { q re, im; } cq;
typedef struct { int64_t hi; uint64_t lo; } s128;
typedef struct { uint32_t x, z; int ph; } pstr;

static int N, H, D, NM, PA = 1, PB = 2;
static pstr ML[MAXN + 1], MR[MAXN + 1];     /* the bare Majoranas (sqrt2 psi) as Pauli strings over all N qubits */
static int64_t HLre[MAXD][MAXD], HLim[MAXD][MAXD], HRre[MAXD][MAXD], HRim[MAXD][MAXD];   /* exact integer entries (Q.40 units) */
static cq UL[2][MAXD][MAXD], UR[2][MAXD][MAXD];    /* [0] = e^{-iH}, [1] = e^{+iH} for dt = 1, each side */
static cq EL[MAXD][MAXD], ER[MAXD][MAXD];          /* e^{-beta H/4} each side */
static cq PsiI[MAXD][MAXD], T[MAXD][MAXD];
static cq WA[MAXD][MAXD], WB[MAXD][MAXD], WC[MAXD][MAXD], WD[MAXD][MAXD], WE[MAXD][MAXD], WF[MAXD][MAXD], WG[MAXD][MAXD], WH[MAXD][MAXD], WI[MAXD][MAXD], WJ[MAXD][MAXD];
static cq ZS[4][MAXD][MAXD], XS[2][MAXD][MAXD];    /* v3: the four states Z_ab and two work matrices */
static s128 AR[MAXD][MAXD], AI[MAXD][MAXD];
static uint64_t chain = 0xcbf29ce484222325ull;
static int taylor_max = 0;

/* ---- 128-bit ---- */
static inline void s128_add(s128 *a, int64_t hi, uint64_t lo) { uint64_t l = a->lo + lo; a->hi += hi + (l < lo); a->lo = l; }
static inline void s128_adds(s128 *a, s128 b) { s128_add(a, b.hi, b.lo); }
static inline void s128_subs(s128 *a, s128 b) { b.lo = ~b.lo + 1; b.hi = ~b.hi + (b.lo == 0); s128_add(a, b.hi, b.lo); }
static inline void mulprod(q a, q b, int64_t *hi, uint64_t *lo)
{
    int neg = (a < 0) ^ (b < 0); uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b, h, l;
    mul128(ua, ub, &h, &l); if (neg) { l = ~l + 1; h = ~h + (l == 0); } *hi = (int64_t)h; *lo = l;
}
static inline void s128_addq(s128 *a, q v) { s128_add(a, v >> (64 - FB), (uint64_t)v << FB); }
/* round(a / 2^sh) to an int64; sh = FB rounds a Q.80 product sum to Q.40, sh = FB + 4 also divides by 16 exactly: ONE rounding */
static inline q s128_round_sh(s128 a, int sh)
{
    int neg = a.hi < 0; if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t r = (a.lo >> sh) | ((uint64_t)a.hi << (64 - sh)); if ((a.lo & ((((uint64_t)1) << sh) - 1)) >= (((uint64_t)1) << (sh - 1))) r++;
    return neg ? -(q)r : (q)r;
}
static inline q s128_round(s128 a) { return s128_round_sh(a, FB); }
static inline void cmac(s128 *re, s128 *im, cq a, cq b)      /* (re, im) += a b */
{
    int64_t h; uint64_t l;
    mulprod(a.re, b.re, &h, &l); s128_add(re, h, l); mulprod(-a.im, b.im, &h, &l); s128_add(re, h, l);
    mulprod(a.re, b.im, &h, &l); s128_add(im, h, l); mulprod(a.im, b.re, &h, &l); s128_add(im, h, l);
}
/* C = A B (D x D), one rounding per entry */
static void mmul(cq A[MAXD][MAXD], cq B[MAXD][MAXD], cq C[MAXD][MAXD])
{
    for (int i = 0; i < D; i++) {
        s128 re[MAXD], im[MAXD]; memset(re, 0, sizeof(s128) * D); memset(im, 0, sizeof(s128) * D);
        for (int k = 0; k < D; k++) { const cq a = A[i][k]; if (!a.re && !a.im) continue; for (int j = 0; j < D; j++) cmac(&re[j], &im[j], a, B[k][j]); }
        for (int j = 0; j < D; j++) { C[i][j].re = s128_round(re[j]); C[i][j].im = s128_round(im[j]); }
    }
}
/* C = A B^T */
static void mmulT(cq A[MAXD][MAXD], cq B[MAXD][MAXD], cq C[MAXD][MAXD])
{
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
        s128 re = {0, 0}, im = {0, 0};
        for (int k = 0; k < D; k++) cmac(&re, &im, A[i][k], B[j][k]);
        C[i][j].re = s128_round(re); C[i][j].im = s128_round(im);
    }
}
/* the exact integer Hamiltonian as a cq matrix (entries are exact Q.40 values) */
static void hmat(int64_t Hre[MAXD][MAXD], int64_t Him[MAXD][MAXD], cq M[MAXD][MAXD]) { for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) { M[i][j].re = (q)Hre[i][j]; M[i][j].im = (q)Him[i][j]; } }
/* E = exp(s M) by Taylor: term_{k+1} = s M term_k / (k+1), each product rounded once; the sum in 128-bit, rounded once */
static void mexp(cq M[MAXD][MAXD], q sre, q sim, cq E[MAXD][MAXD], cq term[MAXD][MAXD], cq u[MAXD][MAXD])
{
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) { AR[i][j].hi = AR[i][j].lo = 0; AI[i][j].hi = AI[i][j].lo = 0; term[i][j].re = term[i][j].im = 0; }
    for (int i = 0; i < D; i++) { term[i][i].re = ONE; s128_addq(&AR[i][i], ONE); }
    int small = 0, k;
    for (k = 1; k <= 80; k++) {
        mmul(M, term, u);
        q mx = 0;
        for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
            q re = mulq(u[i][j].re, sre) - mulq(u[i][j].im, sim), im = mulq(u[i][j].re, sim) + mulq(u[i][j].im, sre);
            term[i][j].re = rdiv(re, k); term[i][j].im = rdiv(im, k);
            s128_addq(&AR[i][j], term[i][j].re); s128_addq(&AI[i][j], term[i][j].im);
            q a = term[i][j].re < 0 ? -term[i][j].re : term[i][j].re, b = term[i][j].im < 0 ? -term[i][j].im : term[i][j].im; if (a > mx) mx = a; if (b > mx) mx = b;
        }
        if (mx <= 1) { if (++small >= 2) break; } else small = 0;
    }
    if (k > taylor_max) taylor_max = k;
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) { E[i][j].re = s128_round(AR[i][j]); E[i][j].im = s128_round(AI[i][j]); }
}

/* ---- Pauli strings and the Majoranas ---- */
static inline int parity32(uint32_t v) { v ^= v >> 16; v ^= v >> 8; v ^= v >> 4; v ^= v >> 2; v ^= v >> 1; return v & 1; }
static int pop32(uint32_t v) { int c = 0; while (v) { v &= v - 1; c++; } return c; }
static pstr pmul(pstr a, pstr b) { pstr r; r.x = a.x ^ b.x; r.z = a.z ^ b.z; r.ph = (a.ph + b.ph + 2 * pop32(a.z & b.x)) & 3; return r; }
static pstr PX(int k) { pstr p = { 1u << k, 0, 0 }; return p; }
static pstr PY(int k) { pstr p = { 1u << k, 1u << k, 1 }; return p; }
static pstr PZ(int k) { pstr p = { 0, 1u << k, 0 }; return p; }
static pstr PID(void) { pstr p = { 0, 0, 0 }; return p; }
static void build_majoranas(void)
{
    for (int j = 1; j <= N; j++) {
        const int site = (j - 1) / 2; pstr s = PID();
        for (int k = 0; k < site; k++) s = pmul(s, PZ(k));
        ML[j] = pmul(s, (j & 1) ? PX(site) : PY(site));
        const int rs = H + site; pstr t = PID();
        for (int k = 0; k < rs; k++) t = pmul(t, PZ(k));
        MR[j] = pmul(t, (j & 1) ? PX(rs) : PY(rs));
    }
}
/* apply a Pauli string (over all N qubits) to Psi: out[L'][R'] = i^ph (-1)^{parity(z & (L | R<<h))} Psi[L][R], L' = L ^ xL, R' = R ^ xR */
static void papply(const pstr *p, cq in[MAXD][MAXD], cq out[MAXD][MAXD])
{
    const uint32_t xL = p->x & (uint32_t)(D - 1), xR = p->x >> H, zL = p->z & (uint32_t)(D - 1), zR = p->z >> H;
    for (int L = 0; L < D; L++) {
        const int pL = parity32(zL & (uint32_t)L);
        for (int R = 0; R < D; R++) {
            const int rot = (p->ph + 2 * (pL ^ parity32(zR & (uint32_t)R))) & 3; const cq v = in[L][R]; cq w;
            switch (rot) { case 0: w = v; break; case 1: w.re = -v.im; w.im = v.re; break; case 2: w.re = -v.re; w.im = -v.im; break; default: w.re = v.im; w.im = -v.re; }
            out[L ^ xL][R ^ xR] = w;
        }
    }
}
/* out = A + i^k B, exact */
static void addrot(cq A[MAXD][MAXD], cq B[MAXD][MAXD], int k, cq out[MAXD][MAXD])
{
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) {
        const cq a = A[L][R], b = B[L][R]; cq r;
        switch (k & 3) { case 0: r = b; break; case 1: r.re = -b.im; r.im = b.re; break; case 2: r.re = -b.re; r.im = -b.im; break; default: r.re = b.im; r.im = -b.re; }
        out[L][R].re = a.re + r.re; out[L][R].im = a.im + r.im;
    }
}
/* out = (2 M_ab) in, M_ab = |a><b| on the left Dirac mode chi = (M_PA + i M_PB)/2: 2M_00 = 1 - i M1 M2, 2M_11 = 1 + i M1 M2, 2M_10 = 2 chi^dag = M1 - i M2, 2M_01 = 2 chi = M1 + i M2 */
static void apply2M(int a, int b, cq in[MAXD][MAXD], cq out[MAXD][MAXD])
{
    if (a == b) { papply(&ML[PB], in, XS[0]); papply(&ML[PA], XS[0], XS[1]); addrot(in, XS[1], a ? 1 : 3, out); }
    else { papply(&ML[PA], in, XS[0]); papply(&ML[PB], in, XS[1]); addrot(XS[0], XS[1], a ? 3 : 1, out); }
}
/* (re, im) += i^k <A|B>, the products exact in 128-bit (no rounding here) */
static void inner_acc(cq A[MAXD][MAXD], cq B[MAXD][MAXD], int k, s128 *re, s128 *im)
{
    s128 pr = {0, 0}, pi = {0, 0}; int64_t h; uint64_t l;
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) {
        mulprod(A[L][R].re, B[L][R].re, &h, &l); s128_add(&pr, h, l); mulprod(A[L][R].im, B[L][R].im, &h, &l); s128_add(&pr, h, l);
        mulprod(A[L][R].re, B[L][R].im, &h, &l); s128_add(&pi, h, l); mulprod(-A[L][R].im, B[L][R].re, &h, &l); s128_add(&pi, h, l);
    }
    switch (k & 3) {
        case 0: s128_adds(re, pr); s128_adds(im, pi); break;
        case 1: s128_subs(re, pi); s128_adds(im, pr); break;      /*  i (x + iy) = -y + ix */
        case 2: s128_subs(re, pr); s128_subs(im, pi); break;
        default: s128_adds(re, pi); s128_subs(im, pr); break;     /* -i (x + iy) =  y - ix */
    }
}
static void build_I(cq Psi[MAXD][MAXD])
{
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { Psi[L][R].re = Psi[L][R].im = 0; }
    for (int L = 0; L < D; L++) {
        int e = 0, pc = 0;
        for (int k = 0; k < H; k++) if ((L >> k) & 1) { e += (k & 1) ? 1 : 3; pc++; }
        e = (e + 2 * (pc * (pc - 1) / 2)) & 3;
        const int R = (D - 1) ^ L;
        switch (e) { case 0: Psi[L][R].re = ONE; break; case 1: Psi[L][R].im = ONE; break; case 2: Psi[L][R].re = -ONE; break; default: Psi[L][R].im = -ONE; }
    }
}
static uint64_t rng;
static uint64_t xs64(void) { rng ^= rng >> 12; rng ^= rng << 25; rng ^= rng >> 27; return rng * 0x2545F4914F6CDD1Dull; }
static q gauss40(void) { q s = 0; for (int i = 0; i < 12; i++) s += (q)(xs64() >> 24); return s - 6 * ONE; }
/* add coef x (the bare 4-product Pauli string) into the side's exact integer matrix: the string must act within that side */
static void addterm(pstr P, q coef, int side)
{
    const uint32_t x = side ? P.x >> H : P.x, z = side ? P.z >> H : P.z;
    if (side && (P.x & (uint32_t)(D - 1))) { printf("WHSYK3: a right term touches the left half -- stop\n"); exit(1); }
    if (side && (P.z & (uint32_t)(D - 1))) { printf("WHSYK3: a right term carries a left parity string -- stop\n"); exit(1); }
    if (!side && ((P.x >> H) || (P.z >> H))) { printf("WHSYK3: a left term touches the right half -- stop\n"); exit(1); }
    int64_t (*Hre)[MAXD] = side ? HRre : HLre, (*Him)[MAXD] = side ? HRim : HLim;
    for (int b = 0; b < D; b++) {
        const int rot = (P.ph + 2 * parity32(z & (uint32_t)b)) & 3; const int bp = b ^ (int)x;
        switch (rot) { case 0: Hre[bp][b] += coef; break; case 1: Him[bp][b] += coef; break; case 2: Hre[bp][b] -= coef; break; default: Him[bp][b] -= coef; }
    }
}
static void build_syk(uint64_t seed)
{
    rng = seed ? seed : 0x9E3779B97F4A7C15ull;
    const q var = divq(6 * ONE, (q)N * N * N * ONE), sd = sqrtq(var);
    memset(HLre, 0, sizeof HLre); memset(HLim, 0, sizeof HLim); memset(HRre, 0, sizeof HRre); memset(HRim, 0, sizeof HRim);
    int nt = 0; q first = 0;
    for (int i = 1; i <= N; i++) for (int j = i + 1; j <= N; j++) for (int k = j + 1; k <= N; k++) for (int l = k + 1; l <= N; l++) {
        const q J = mulq(gauss40(), sd), c = -rdiv(J, 4); if (!nt) first = c;
        addterm(pmul(pmul(pmul(ML[i], ML[j]), ML[k]), ML[l]), c, 0);
        addterm(pmul(pmul(pmul(MR[i], MR[j]), MR[k]), MR[l]), c, 1);
        nt++;
    }
    NM = N; { char fb[32]; fmt(first, fb); printf("WHSYK3 seed %llu: %d terms, first coefficient %s (v1 and v2 print the same two numbers: the same couplings)\n", (unsigned long long)seed, nt, fb); }
}
static void build_model1(void)
{
    static const int idx[5][4] = { {1, 2, 4, 5}, {1, 3, 4, 7}, {1, 3, 5, 6}, {2, 3, 4, 6}, {2, 3, 5, 7} };
    static const double cf[5] = { -0.36, 0.19, -0.71, 0.22, 0.49 };
    memset(HLre, 0, sizeof HLre); memset(HLim, 0, sizeof HLim); memset(HRre, 0, sizeof HRre); memset(HRim, 0, sizeof HRim);
    for (int t = 0; t < 5; t++) {
        const double c = cf[t] * (double)ONE / 4.0; const q cq40 = (q)(c < 0 ? c - 0.5 : c + 0.5);
        addterm(pmul(pmul(pmul(ML[idx[t][0]], ML[idx[t][1]]), ML[idx[t][2]]), ML[idx[t][3]]), cq40, 0);
        addterm(pmul(pmul(pmul(MR[idx[t][0]], MR[idx[t][1]]), MR[idx[t][2]]), MR[idx[t][3]]), cq40, 1);
    }
    NM = 7;
}
/* the coupling: Psi <- PROD_j (c + i ss W_j) Psi, W_j = i bare_l^j bare_r^j (a Pauli string over both halves) */
static void apply_U(cq Psi[MAXD][MAXD], q c, q s, int sign, cq w[MAXD][MAXD])
{
    for (int j = 1; j <= NM; j++) {
        pstr W = pmul(ML[j], MR[j]); W.ph = (W.ph + 1) & 3;
        papply(&W, Psi, w);
        const q ss = sign > 0 ? s : -s;
        for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) {
            const q wr = w[L][R].re, wi = w[L][R].im;
            const q re = mulq(Psi[L][R].re, c) - mulq(wi, ss), im = mulq(Psi[L][R].im, c) + mulq(wr, ss);
            Psi[L][R].re = re; Psi[L][R].im = im;
        }
    }
}
/* Psi <- U_L(t) Psi U_R(t)^T with U(t) = U(1)^t, sign 0 = e^{-iHt}, 1 = e^{+iHt}; t a non-negative integer */
static void evolve(cq Psi[MAXD][MAXD], int t, int sign, cq w[MAXD][MAXD])
{
    for (int step = 0; step < t; step++) { mmul(UL[sign], Psi, w); mmulT(w, UR[sign], Psi); }
}
static q norm2(cq Psi[MAXD][MAXD])
{
    s128 a = {0, 0}; int64_t h; uint64_t l;
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { mulprod(Psi[L][R].re, Psi[L][R].re, &h, &l); s128_add(&a, h, l); mulprod(Psi[L][R].im, Psi[L][R].im, &h, &l); s128_add(&a, h, l); }
    return s128_round(a);
}
static cq inner(cq A[MAXD][MAXD], cq B[MAXD][MAXD])
{
    s128 re = {0, 0}, im = {0, 0}; inner_acc(A, B, 0, &re, &im);
    cq r = { s128_round(re), s128_round(im) }; return r;
}
static void normalize(cq Psi[MAXD][MAXD]) { const q n = sqrtq(norm2(Psi)); for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { Psi[L][R].re = divq(Psi[L][R].re, n); Psi[L][R].im = divq(Psi[L][R].im, n); } }
static void copy(cq A[MAXD][MAXD], cq B[MAXD][MAXD]) { memcpy(B, A, sizeof(cq) * MAXD * MAXD); }
static q maxabs(cq A[MAXD][MAXD]) { q m = 0; for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { q a = A[L][R].re < 0 ? -A[L][R].re : A[L][R].re, b = A[L][R].im < 0 ? -A[L][R].im : A[L][R].im; if (a > m) m = a; if (b > m) m = b; } return m; }
static void hashq(q v) { chain = fnv(chain, &v, 8); }
static q ent(q lam, int *clamp) { if (lam <= 0) { *clamp = 1; return 0; } return -mulq(lam, lnq(lam)); }
static q qabs(q v) { return v < 0 ? -v : v; }
static q parseq(const char *s)
{
    int neg = 0; if (*s == '-') { neg = 1; s++; }
    q ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (*s++ - '0');
    q fp = 0, scale = 1; if (*s == '.') { s++; while (*s >= '0' && *s <= '9' && scale < 1000000000000ll) { fp = fp * 10 + (*s++ - '0'); scale *= 10; } }
    q v = ip * ONE + (fp ? divq(fp * ONE, scale * ONE) : 0); return neg ? -v : v;
}
static int parselist(const char *s, q *out, int max) { int n = 0; char buf[256]; strncpy(buf, s, 255); buf[255] = 0; for (char *tok = strtok(buf, ","); tok && n < max; tok = strtok(NULL, ",")) out[n++] = parseq(tok); return n; }

int main(int argc, char **argv)
{
    N = 8; int model = 0, vnormN = 0, formulas = 1; q beta = 4 * ONE; q mus[16], ts[32], tps[32]; int nmu = 0, nt = 0, ntp = 0; uint64_t seeds[16]; int nseed = 0;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "model=", 6)) model = !strcmp(argv[i] + 6, "m1");
        else if (!strncmp(argv[i], "beta=", 5)) beta = parseq(argv[i] + 5);
        else if (!strncmp(argv[i], "mu=", 3)) nmu = parselist(argv[i] + 3, mus, 16);
        else if (!strncmp(argv[i], "t=", 2)) nt = parselist(argv[i] + 2, ts, 32);
        else if (!strncmp(argv[i], "tp=", 3)) ntp = parselist(argv[i] + 3, tps, 32);
        else if (!strncmp(argv[i], "pair=", 5)) sscanf(argv[i] + 5, "%d,%d", &PA, &PB);
        else if (!strncmp(argv[i], "vnorm=", 6)) vnormN = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "formulas=", 9)) formulas = atoi(argv[i] + 9);
        else if (!strncmp(argv[i], "seeds=", 6)) { char b[256]; strncpy(b, argv[i] + 6, 255); b[255] = 0; for (char *t = strtok(b, ","); t && nseed < 16; t = strtok(NULL, ",")) seeds[nseed++] = strtoull(t, 0, 10); }
    }
    if (N < 4 || N > MAXN || (N & 1)) { printf("WHSYK3: N must be even, 4..%d\n", MAXN); return 1; }
    if (!nmu) { mus[0] = -12 * ONE; mus[1] = 12 * ONE; nmu = 2; }
    if (!nt) { ts[0] = ONE; ts[1] = 2 * ONE; ts[2] = 3 * ONE; nt = 3; }
    if (!nseed) { seeds[0] = 1; nseed = 1; }
    if (model) { N = 8; nseed = 1; if (!vnormN) vnormN = 10; }
    if (!vnormN) vnormN = N;
    if (PA < 1 || PA > N || PB < 1 || PB > N || PA == PB) { printf("WHSYK3: pair must be two different Majorana indices in 1..N\n"); return 1; }
    H = N / 2; D = 1 << H;
    char A[32], B[32], C[32], Dd[32], E[32], F[32], G[32], S1[32], S2[32], S3[32], S4[32], S5[32], S6[32], S7[32], S8[32];
    fmt(beta, A);
    printf("WHSYK3: Q.%d (whq.h); %s; N = %d Majoranas a side, FACTORIZED (left on qubits 0..%d, right on %d..%d; D = %d a side); beta = %s; V normalized with N = %d; the Dirac modes chi = (psi^%d + i psi^%d)/sqrt2; the evolution U(1)^t by Taylor matrix exponentials; the EXACT rho_TR from its six elements (GJ (2.15)-(2.20)) %s\n",
           FB, model ? "MODEL 1 of the 2022 experiment (Kobrin-Schuster-Yao eq. 1: 7 Majoranas, 5 terms)" : "the FULL q = 4 SYK, <J^2> = 6/N^3, J = 1", N, H - 1, H, N - 1, D, A, vnormN, PA, PB,
           formulas ? "beside the averaged forms (2.63)-(2.66) on the same seed" : "alone (formulas=0)");
    build_majoranas();
    build_I(PsiI);
    /* self-test 1: (bare_l^j + i bare_r^j) |I> = 0 exactly for every j (unnormalized |I>) */
    q worst = 0;
    for (int j = 1; j <= N; j++) {
        papply(&ML[j], PsiI, WA); papply(&MR[j], PsiI, WB);
        for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { q re = WA[L][R].re - WB[L][R].im, im = WA[L][R].im + WB[L][R].re; re = qabs(re); im = qabs(im); if (re > worst) worst = re; if (im > worst) worst = im; }
    }
    fmt(worst, A); printf("WHSYK3 selftest 1: max |(psi_l^j + i psi_r^j)|I>| = %s (must be 0.000000000000)\n", A);
    for (int si = 0; si < nseed; si++) {
        if (model) build_model1(); else build_syk(seeds[si]);
        /* self-test 2: H_L Psi_I = Psi_I H_R^T exactly (the matrices are exact integers, |I> is exact: no rounding anywhere) */
        {
            hmat(HLre, HLim, WC); hmat(HRre, HRim, WD);
            mmul(WC, PsiI, WA); mmulT(PsiI, WD, WB);
            q w2 = 0; for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { q re = qabs(WA[L][R].re - WB[L][R].re), im = qabs(WA[L][R].im - WB[L][R].im); if (re > w2) w2 = re; if (im > w2) w2 = im; }
            fmt(w2, A); fmt(maxabs(WA), B);
            printf("WHSYK3 seed %llu selftest 2: max |H_L Psi_I - Psi_I H_R^T| = %s (must be 0); |H_L Psi_I| max entry %s\n", (unsigned long long)seeds[si], A, B);
            hashq(w2);
        }
        /* the matrix exponentials: U(+-1) each side, e^{-beta H/4} each side */
        hmat(HLre, HLim, WC); hmat(HRre, HRim, WD);
        mexp(WC, 0, -ONE, UL[0], WE, WF); mexp(WC, 0, ONE, UL[1], WE, WF);
        mexp(WD, 0, -ONE, UR[0], WE, WF); mexp(WD, 0, ONE, UR[1], WE, WF);
        const q bq = rdiv(beta, 4);
        mexp(WC, -bq, 0, EL, WE, WF); mexp(WD, -bq, 0, ER, WE, WF);
        /* the thermofield double: E_L Psi_I E_R^T, normalized */
        mmul(EL, PsiI, WA); mmulT(WA, ER, T);
        const q nT = norm2(T), nI = norm2(PsiI); normalize(T);
        fmt(divq(nT, nI), A); fmt(norm2(T), B);
        printf("WHSYK3 seed %llu: |e^{-beta H/4} I|^2 / |I|^2 = %s (Z_beta / Z_inf); normalized to %s; Taylor terms at most %d\n", (unsigned long long)seeds[si], A, B, taylor_max);
        hashq(divq(nT, nI));
        for (int mi = 0; mi < nmu; mi++) {
            const q arg = divq(mus[mi], (q)(2 * 4 * vnormN) * ONE), c = cosq(arg), s = sinq(arg);
            for (int ti = 0; ti < nt; ti++) {
                const q tq = ts[ti], tpq = ntp ? tps[ti < ntp ? ti : ntp - 1] : tq;
                const int t = (int)(tq / ONE), tp = (int)(tpq / ONE);
                if ((q)t * ONE != tq || (q)tp * ONE != tpq || t < 0 || tp < 0) { printf("WHSYK3: v3 takes non-negative integer t and t' (U(1)^t); got otherwise -- skipped\n"); continue; }
                /* ---- THE EXACT MATRIX ---- */
                copy(T, WG); evolve(WG, tp, 1, WE);                                   /* WG = e^{iHt'} T */
                for (int ab = 0; ab < 4; ab++) {                                      /* Z_ab = e^{-iHt} U e^{-iHt'} (2 M_ab) e^{iHt'} T */
                    apply2M(ab >> 1, ab & 1, WG, ZS[ab]); evolve(ZS[ab], tp, 0, WE); apply_U(ZS[ab], c, s, 1, WE); evolve(ZS[ab], t, 0, WE);
                }
                s128 ar[6], ai[6]; memset(ar, 0, sizeof ar); memset(ai, 0, sizeof ai);   /* 0: rho00, 1: rho11, 2: rho22, 3: rho33, 4: rho03, 5: rho12 */
                for (int b = 0; b < 2; b++) {
                    cq (*Z0)[MAXD] = ZS[b], (*Z1)[MAXD] = ZS[2 + b];
                    papply(&MR[PB], Z0, XS[0]); papply(&MR[PA], XS[0], XS[1]);                   /* XS[1] = R1 R2 Z0 */
                    inner_acc(Z0, Z0, 0, &ar[0], &ai[0]); inner_acc(Z0, XS[1], 3, &ar[0], &ai[0]);   /* rho00 += <Z0|(1 - i R1 R2)|Z0> */
                    inner_acc(Z0, Z0, 0, &ar[1], &ai[1]); inner_acc(Z0, XS[1], 1, &ar[1], &ai[1]);   /* rho11 += <Z0|(1 + i R1 R2)|Z0> */
                    papply(&MR[PA], Z0, XS[0]); papply(&MR[PB], Z0, XS[1]);                      /* R1 Z0, R2 Z0 */
                    inner_acc(Z1, XS[0], 0, &ar[4], &ai[4]); inner_acc(Z1, XS[1], 3, &ar[4], &ai[4]);  /* rho03 += <Z1|(R1 - i R2)|Z0> = <Z1| 2 chi_r^dag |Z0> */
                    inner_acc(Z1, XS[0], 0, &ar[5], &ai[5]); inner_acc(Z1, XS[1], 1, &ar[5], &ai[5]);  /* rho12 += <Z1|(R1 + i R2)|Z0> = <Z1| 2 chi_r |Z0> */
                    papply(&MR[PB], Z1, XS[0]); papply(&MR[PA], XS[0], XS[1]);                   /* XS[1] = R1 R2 Z1 */
                    inner_acc(Z1, Z1, 0, &ar[2], &ai[2]); inner_acc(Z1, XS[1], 3, &ar[2], &ai[2]);   /* rho22 */
                    inner_acc(Z1, Z1, 0, &ar[3], &ai[3]); inner_acc(Z1, XS[1], 1, &ar[3], &ai[3]);   /* rho33 */
                }
                q rr[6], ri[6]; for (int e = 0; e < 6; e++) { rr[e] = s128_round_sh(ar[e], FB + 4); ri[e] = s128_round_sh(ai[e], FB + 4); }
                const q r00 = rr[0], r11 = rr[1], r22 = rr[2], r33 = rr[3];
                const q epr1 = r00 + r11 - (ONE >> 1), epr2 = r22 + r33 - (ONE >> 1);
                q imdiag = 0; for (int e = 0; e < 4; e++) if (qabs(ri[e]) > imdiag) imdiag = qabs(ri[e]);
                /* the entropies: S(RT) from the blocks {0,3} and {1,2}; S(R) and S(T) from the diagonal marginals */
                int clamp = 0; q SRT = 0, SR = 0, ST = 0;
                { const q m = rdiv(r00 + r33, 2), d = rdiv(r00 - r33, 2), off2 = mulq(rr[4], rr[4]) + mulq(ri[4], ri[4]), disc = sqrtq(mulq(d, d) + off2);
                  SRT += ent(m + disc, &clamp) + ent(m - disc, &clamp); }
                { const q m = rdiv(r11 + r22, 2), d = rdiv(r11 - r22, 2), off2 = mulq(rr[5], rr[5]) + mulq(ri[5], ri[5]), disc = sqrtq(mulq(d, d) + off2);
                  SRT += ent(m + disc, &clamp) + ent(m - disc, &clamp); }
                ST = ent(r00 + r22, &clamp) + ent(r11 + r33, &clamp);   /* the readout T's marginal (not 1/2, 1/2 for one seed) */
                SR = ent(r00 + r11, &clamp) + ent(r22 + r33, &clamp);   /* the reference R's marginal: (1/2, 1/2), so S(R) = ln 2 */
                const q I = SR + ST - SRT, frac = divq(I, 2 * LN2);
                /* K = <{A1, B1}> = Re <M1(-t') T | U^dag R1(t) U T> (the bare Majoranas carry the two 1/sqrt2) */
                q K, rho11f = 0, r14re = 0, r14im = 0, fracf = 0; int clampf = 0;
                if (formulas) {
                    /* v2's block, unchanged: the averaged forms on this seed */
                    papply(&ML[PA], WG, WA); evolve(WA, tp, 0, WE);                  /* a1 */
                    papply(&ML[PB], WG, WB); evolve(WB, tp, 0, WE);                  /* a2 */
                    copy(WA, WH); evolve(WH, tp, 1, WE); papply(&ML[PB], WH, WC); evolve(WC, tp, 0, WE);   /* a21 */
                    copy(WB, WH); evolve(WH, tp, 1, WE); papply(&ML[PA], WH, WD); evolve(WD, tp, 0, WE);   /* a12 */
                    copy(T, WH); apply_U(WH, c, s, 1, WE); evolve(WH, t, 0, WE);     /* up */
                    papply(&MR[PA], WH, WF); evolve(WF, t, 1, WE); apply_U(WF, c, s, -1, WE);              /* b1 = WF */
                    const cq Y1 = inner(WA, WF);
                    papply(&MR[PB], WH, WI); papply(&MR[PA], WI, WJ); evolve(WJ, t, 1, WE); apply_U(WJ, c, s, -1, WE);   /* bb = B12 T = WJ */
                    const cq X1 = inner(WC, WJ);
                    copy(WB, WH); apply_U(WH, c, s, 1, WE); evolve(WH, t, 0, WE);
                    papply(&MR[PB], WH, WI); papply(&MR[PA], WI, WJ); evolve(WJ, t, 1, WE); apply_U(WJ, c, s, -1, WE);   /* B12 a2 */
                    const cq X2 = inner(WA, WJ);
                    copy(WA, WH); apply_U(WH, c, s, 1, WE); evolve(WH, t, 0, WE);
                    papply(&MR[PB], WH, WI); evolve(WI, t, 1, WE); apply_U(WI, c, s, -1, WE);                         /* B2 a1 */
                    const cq Y2 = inner(WD, WI);
                    copy(WD, WH); apply_U(WH, c, s, 1, WE); evolve(WH, t, 0, WE);
                    papply(&MR[PB], WH, WI); evolve(WI, t, 1, WE); apply_U(WI, c, s, -1, WE);                         /* B2 a12 */
                    const cq Y3 = inner(WA, WI);
                    rho11f = (ONE >> 1) - rdiv(X1.re - X2.re, 4);
                    K = Y1.re;
                    r14re = rdiv(Y1.re, 2) - rdiv(Y2.re + Y3.re, 4); r14im = -rdiv(Y2.im + Y3.im, 4);
                    const q m14 = sqrtq(mulq(r14re, r14re) + mulq(r14im, r14im));
                    const q l1 = rdiv(rho11f + m14, 2), l2 = rdiv(rho11f - m14, 2), l3 = rdiv(ONE - rho11f, 2);
                    q Sf = 0; const q lam[4] = { l1, l2, l3, l3 };
                    for (int e = 0; e < 4; e++) Sf += ent(lam[e], &clampf);
                    fracf = divq(2 * LN2 - Sf, 2 * LN2);
                } else {
                    papply(&ML[PA], WG, WA); evolve(WA, tp, 0, WE);                  /* a1 = M1(-t') T */
                    copy(T, WH); apply_U(WH, c, s, 1, WE); evolve(WH, t, 0, WE);
                    papply(&MR[PA], WH, WF); evolve(WF, t, 1, WE); apply_U(WF, c, s, -1, WE);              /* b1 = U^dag R1(t) U T */
                    K = inner(WA, WF).re;
                }
                char HB[32], HT[32]; fmt(beta, HB); fmt(tpq, HT);
                fmt(tq, A); fmt(mus[mi], B); fmt(r00, C); fmt(r11, Dd); fmt(r22, E); fmt(r33, F); fmt(rr[4], G); fmt(ri[4], S1); fmt(rr[5], S2); fmt(ri[5], S3);
                fmt(epr1, S4); fmt(epr2, S5); fmt(imdiag, S6); fmt(K, S7); fmt(frac, S8);
                printf("WHSYK3 N=%d seed=%llu beta=%s mu=%s t=%s tp=%s: EXACT rho00 %s rho11 %s rho22 %s rho33 %s | rho03 %s %s i | rho12 %s %s i | checks (the reference R's marginal, exact) rho00+rho11-1/2 %s rho22+rho33-1/2 %s max|Im diag| %s | K %s | I_RT/(2 ln 2) EXACT %s%s",
                       N, (unsigned long long)seeds[si], HB, B, A, HT, C, Dd, E, F, G, S1, S2, S3, S4, S5, S6, S7, S8, clamp ? " (an eigenvalue clamped at 0)" : "");
                hashq(r00); hashq(r11); hashq(r22); hashq(r33); hashq(rr[4]); hashq(ri[4]); hashq(rr[5]); hashq(ri[5]); hashq(K); hashq(frac);
                if (formulas) {
                    char F1[32], F2[32], F3[32], F4[32], F5[32]; fmt(rho11f, F1); fmt(r14re, F2); fmt(r14im, F3); fmt(fracf, F4); fmt(frac - fracf, F5);
                    printf(" | AVERAGED FORMS (2.63)-(2.66) on this seed: rho11 %s rho14 %s %s i I_RT/(2 ln 2) %s%s | exact - averaged %s", F1, F2, F3, F4, clampf ? " (clamped)" : "", F5);
                    hashq(rho11f); hashq(r14re); hashq(r14im); hashq(fracf);
                }
                printf(" | Taylor terms <= %d\n", taylor_max);
            }
        }
    }
    printf("WHSYK3: chain %016llx\n", (unsigned long long)chain);
    printf("WHSYK3: done\n");
    return 0;
}
