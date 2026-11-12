/* wh_syk4.c (11-11-2026, 10:30 PM) -- WORMHOLE MATH 2, v4: the SYK traversable-wormhole teleportation protocol (Gao-Jafferis, arXiv:1911.07416)
   exact at finite N for ONE realization of the couplings, the full six-element density matrix rho_TR of v3, NOW IN Q.60 END TO END, and
   with the time steps built once per seed. Why v4 (the review's bit-exact emulation, 11-11 10:20 PM): v3's Q.40 chain U(1)^t left up to 21
   units of 2^-40 in the evolved states at t = 6, and the Q.40 cosine and sine of the coupling (c^2 + s^2 - 1 = -2.4 units) moved K by up to
   27 units, so v2/v3's twelve printed places were good to ten. v4 carries everything in Q.60 (int64 with 60 fraction bits, products in
   128 bits, one rounding per entry per product) and rounds nothing to Q.40 at all: the printed twelve places are rounded from Q.60 values
   whose own error is about 1e-16.

   THE PROTOCOL AND THE MATRIX: exactly v3's (the derivation is in wh_syk3.c's header): the message is the left Dirac mode
   chi_l = (psi_l^1 + i psi_l^2)/sqrt2 at -t', entangled with the reference R (half of an EPR pair); the readout is the right mode chi_r at t
   after the coupling U = e^{i mu V}; in the basis |a c> (a = R, c = T, index k = 2a + c) the six elements are
       rho_(ac),(a'c') = (1/2) SUM_b < M_ba'(-t') U^dag N_c'c(t) U M_ab(-t') >,
   built from the four states Z_ab = e^{-iHt} U e^{-iHt'} (2 M_ab) e^{iHt'} |TFD> as (1/16) SUM_b <Z_a'b | 2N_c'c | Z_ab>; the reference's
   marginal rho_00 + rho_11 = rho_22 + rho_33 = 1/2 exactly (printed as a check); I_RT = S(R) + S(T) - S(RT); K = <{A1, B1}>.

   THE ARITHMETIC, v4:
   - H_L and H_R: exact integers in Q.60 (the Q.40 coefficients -J/4 of v1-v3 shifted by 20 bits: the same couplings, bit for bit).
   - ONE set of powers per seed and side: P_k = H P_{k-1} / k in Q.60 (one rounding per entry per k), summed twice: S_W = SUM (-i)^k P_k =
     e^{-iH} (the rotation by (-i)^k is exact) and S_E = SUM (-beta/4)^k P_k = e^{-beta H/4} (exact for beta = 4), both accumulated in
     128 bits and rounded once; the series stops when a term is below one unit of 2^-60 twice (about 22 terms at N = 14).
   - W(t) = W(1) W(t-1) for t = 2..tmax, once per seed, in Q.60; e^{+iHt} is W(t)^dag, the conjugate transpose, no new rounding; an
     evolution by t is then ONE pair of products W_L(t) Psi W_R(t)^T whatever t is (v3 took t pairs and recomputed them for every point).
   - Products skip exact zeros, so the block structure of H, W and E (fermion parity) is used without bookkeeping: half the multiplies.
   - The coupling's cosine and sine of mu/(8 N) in Q.60 (Taylor, 128-bit products; c^2 + s^2 - 1 printed in units of 2^-60), the coupling
     applied in Q.60 with one rounding per entry per Majorana pair.
   - The 2M and 2N operators are sums of two Pauli strings with phases: exact. The six elements: twelve inner products in 128 bits, divided
     by 16 and rounded once to Q.60. The eigenvalues by a 128-bit integer square root; the entropies by a Q.60 logarithm (range reduction
     by powers of two, then 2 artanh((m-1)/(m+1)) to 2^-60; ln 2 itself computed the same way and printed as a self-test).
   - |I> is built normalized exactly (amplitude 2^{-H/2} for H = N/2 even, (1 + i) 2^{-(H+1)/2} for H odd: a global phase, the review's
     trick), so |I|^2 = 1 and no sum of squares approaches 2^127 (the first build's unnormalized |I> had |I|^2 = D 2^120 = 2^127 at D = 128
     and overflowed the 128-bit norm there: found at N = 14 before any run); the thermofield double is then normalized in Q.60 (a 128-bit
     norm, an integer square root, a 128-bit division per entry).
   Self-tests printed: ln 2; c^2 + s^2 - 1; (psi_l^j + i psi_r^j)|I> = 0 exactly; H_L Psi_I = Psi_I H_R^T exactly; e^{-beta H_L/4} Psi_I =
   Psi_I (e^{-beta H_R/4})^T to the series' rounding; W(1)^dag W(1) - 1 (the unitarity defect, in units of 2^-60); the reference's marginal.
   Command line: whsyk4 N=14 seeds=1,2,3 beta=4 mu=-12,12,0 t=1,2,3,4,5,6 [tp=...] [pair=1,2] [vnorm=N]; integer t and t' (t <= 8). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "whq.h"

#define MAXN 16
#define MAXD (1 << (MAXN / 2))
#define TMAX 8
#define F60 60
static const int64_t ONE60 = (int64_t)1 << F60;

typedef struct { int64_t re, im; } c60;
typedef struct { int64_t hi; uint64_t lo; } s128;
typedef struct { uint32_t x, z; int ph; } pstr;

static int N, H, D, NM, PA = 1, PB = 2;
static pstr ML[MAXN + 1], MR[MAXN + 1];
static c60 HL[MAXD][MAXD], HR[MAXD][MAXD];              /* exact Q.60 integers */
static c60 WL[TMAX + 1][MAXD][MAXD], WR[TMAX + 1][MAXD][MAXD];   /* W(t) = e^{-iHt}, t = 0..tmax, each side */
static c60 EL[MAXD][MAXD], ER[MAXD][MAXD];              /* e^{-beta H/4} */
static c60 PsiI[MAXD][MAXD], T[MAXD][MAXD];
static c60 WA[MAXD][MAXD], WB[MAXD][MAXD], WC[MAXD][MAXD], WG[MAXD][MAXD], WX[MAXD][MAXD];
static c60 ZS[4][MAXD][MAXD], XS[2][MAXD][MAXD];
static s128 SWre[MAXD][MAXD], SWim[MAXD][MAXD], SEre[MAXD][MAXD], SEim[MAXD][MAXD];   /* the two series' accumulators */
static uint64_t chain = 0xcbf29ce484222325ull;
static int taylor_max = 0; static int64_t unit_defect_max = 0;
static int64_t LN2_60;
static int EOFL[MAXD];                                   /* e(L), the phase exponent of |I>'s amplitude on |L>|~L>, kept for self-test 5 */

/* ---- 128-bit ---- */
static inline void s128_add(s128 *a, int64_t hi, uint64_t lo) { uint64_t l = a->lo + lo; a->hi += hi + (l < lo); a->lo = l; }
static inline void s128_adds(s128 *a, s128 b) { s128_add(a, b.hi, b.lo); }
static inline void s128_subs(s128 *a, s128 b) { b.lo = ~b.lo + 1; b.hi = ~b.hi + (b.lo == 0); s128_add(a, b.hi, b.lo); }
static inline void mulprod(int64_t a, int64_t b, int64_t *hi, uint64_t *lo)
{
    int neg = (a < 0) ^ (b < 0); uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b, h, l;
    mul128(ua, ub, &h, &l); if (neg) { l = ~l + 1; h = ~h + (l == 0); } *hi = (int64_t)h; *lo = l;
}
static inline void s128_addp(s128 *a, int64_t x, int64_t y) { int64_t h; uint64_t l; mulprod(x, y, &h, &l); s128_add(a, h, l); }
/* round(a / 2^sh), sh in 56..64, result must fit int64 */
static inline int64_t s128_round_sh(s128 a, int sh)
{
    int neg = a.hi < 0; if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t r = (a.lo >> sh) | ((uint64_t)a.hi << (64 - sh)); if ((a.lo & ((((uint64_t)1) << sh) - 1)) >= (((uint64_t)1) << (sh - 1))) r++;
    return neg ? -(int64_t)r : (int64_t)r;
}
static inline int64_t mul60(int64_t a, int64_t b) { s128 s = {0, 0}; s128_addp(&s, a, b); return s128_round_sh(s, F60); }   /* round(a b / 2^60) */
static inline int64_t rdiv64(int64_t x, int64_t n) { return x >= 0 ? (x + n / 2) / n : -((-x + n / 2) / n); }
/* round(a 2^60 / b) for b > 0, |a| < 8 b */
static int64_t div60(int64_t a, int64_t b)
{
    int neg = a < 0; uint64_t ua = neg ? -(uint64_t)a : (uint64_t)a;
    const uint64_t hi = ua >> 4, lo = ua << 60;
    if (hi >= (uint64_t)b) { printf("WHSYK4: div60 overflow -- stop\n"); exit(1); }
    uint64_t rem; uint64_t qq = div128(hi, lo, (uint64_t)b, &rem); if (2 * rem >= (uint64_t)b) qq++;
    return neg ? -(int64_t)qq : (int64_t)qq;
}
/* the integer square root of a 128-bit radicand (hi:lo), digit by digit; the result < 2^64 */
static uint64_t isqrt128(uint64_t hi, uint64_t lo)
{
    uint64_t rhi = 0, rlo = 0, root = 0;
    for (int i = 0; i < 64; i++) {
        const int top = i < 32 ? (int)((hi >> (62 - 2 * i)) & 3) : (int)((lo >> (126 - 2 * i)) & 3);
        rhi = (rhi << 2) | (rlo >> 62); rlo = (rlo << 2) | (uint64_t)top;
        const uint64_t thi = root >> 62, tlo = (root << 2) | 1;
        root <<= 1;
        if (rhi > thi || (rhi == thi && rlo >= tlo)) { const uint64_t nlo = rlo - tlo; rhi = rhi - thi - (rlo < tlo); rlo = nlo; root |= 1; }
    }
    return root;
}
/* sqrt of a non-negative Q.120 value held in an s128 -> Q.60 (floor) */
static int64_t sqrt120(s128 a) { if (a.hi < 0) return 0; return (int64_t)isqrt128((uint64_t)a.hi, a.lo); }
/* ln x for x > 0 in Q.60: x = 2^k m, m in [1/2, 1); ln m = 2 artanh((m - 1)/(m + 1)), |z| <= 1/3, 24 terms */
static int64_t artanh2_60(int64_t z)                   /* 2 (z + z^3/3 + z^5/5 + ...) */
{
    const int64_t z2 = mul60(z, z); int64_t term = z, sum = 0;
    for (int k = 0; k < 40; k++) { const int64_t add = rdiv64(term, 2 * k + 1); if (!add && k > 2) break; sum += add; term = mul60(term, z2); if (!term) break; }
    return 2 * sum;
}
static int64_t ln60(int64_t x)
{
    if (x <= 0) { printf("WHSYK4: ln of a non-positive value -- stop\n"); exit(1); }
    int k = 0; int64_t m = x;
    while (m < (ONE60 >> 1)) { m <<= 1; k--; }
    while (m >= ONE60) { m >>= 1; k++; }                /* m in [1/2, 1); x = m 2^k; division by 2 is exact only when m is even: carry the lost bit */
    /* the shifts right above lose at most one bit each: redo exactly for x >= 1 by dividing with rounding (x < 8 here, so at most 3 shifts) */
    m = x; k = 0;
    while (m < (ONE60 >> 1)) { m <<= 1; k--; }
    while (m >= ONE60) { m = rdiv64(m, 2); k++; }
    const int64_t z = div60(m - ONE60, m + ONE60);
    return artanh2_60(z) + (int64_t)k * LN2_60;
}
static int64_t ent60(int64_t lam, int *clamp) { if (lam <= 0) { *clamp = 1; return 0; } return -mul60(lam, ln60(lam)); }
/* cos and sin of a small x (|x| < 1) in Q.60 by Taylor, 128-bit products */
static void cossin60(int64_t x, int64_t *c, int64_t *s)
{
    const int64_t x2 = mul60(x, x); int64_t tc = ONE60, ts = x, sc = ONE60, ss = x;
    for (int k = 1; k < 30; k++) {
        tc = -rdiv64(mul60(tc, x2), (int64_t)(2 * k - 1) * (2 * k)); ts = -rdiv64(mul60(ts, x2), (int64_t)(2 * k) * (2 * k + 1));
        sc += tc; ss += ts; if (!tc && !ts) break;
    }
    *c = sc; *s = ss;
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
static void hash64(int64_t v) { chain = fnv(chain, &v, 8); }

/* ---- matrices ---- */
/* C = op(A) B, op = identity (mode 0) or conjugate transpose (mode 1); exact zeros of op(A) skipped; one rounding per entry */
static void mmul(c60 A[MAXD][MAXD], int mode, c60 B[MAXD][MAXD], c60 C[MAXD][MAXD])
{
    for (int i = 0; i < D; i++) {
        s128 re[MAXD], im[MAXD]; memset(re, 0, sizeof(s128) * D); memset(im, 0, sizeof(s128) * D);
        for (int k = 0; k < D; k++) {
            c60 a = mode ? A[k][i] : A[i][k]; if (mode) a.im = -a.im;
            if (!a.re && !a.im) continue;
            for (int j = 0; j < D; j++) {
                const c60 b = B[k][j];
                s128_addp(&re[j], a.re, b.re); s128_addp(&re[j], -a.im, b.im);
                s128_addp(&im[j], a.re, b.im); s128_addp(&im[j], a.im, b.re);
            }
        }
        for (int j = 0; j < D; j++) { C[i][j].re = s128_round_sh(re[j], F60); C[i][j].im = s128_round_sh(im[j], F60); }
    }
}
static void transpose(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD]) { for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) B[j][i] = A[i][j]; }
static void copy(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD]) { memcpy(B, A, sizeof(c60) * MAXD * MAXD); }
/* Psi <- op(A) Psi op(B)^T : X = op(A) Psi; Psi = (op(B) X^T)^T */
static void sandwich(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD], int mode, c60 Psi[MAXD][MAXD])
{
    mmul(A, mode, Psi, WX); transpose(WX, Psi); mmul(B, mode, Psi, WX); transpose(WX, Psi);
}
/* the evolution by t steps: sign 0 = e^{-iHt} (W(t)), sign 1 = e^{+iHt} (W(t)^dag) */
static void evolve(c60 Psi[MAXD][MAXD], int t, int sign) { if (t > 0) sandwich(WL[t], WR[t], sign, Psi); }

/* ---- Pauli strings and the Majoranas (as v2/v3) ---- */
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
static void papply(const pstr *p, c60 in[MAXD][MAXD], c60 out[MAXD][MAXD])
{
    const uint32_t xL = p->x & (uint32_t)(D - 1), xR = p->x >> H, zL = p->z & (uint32_t)(D - 1), zR = p->z >> H;
    for (int L = 0; L < D; L++) {
        const int pL = parity32(zL & (uint32_t)L);
        for (int R = 0; R < D; R++) {
            const int rot = (p->ph + 2 * (pL ^ parity32(zR & (uint32_t)R))) & 3; const c60 v = in[L][R]; c60 w;
            switch (rot) { case 0: w = v; break; case 1: w.re = -v.im; w.im = v.re; break; case 2: w.re = -v.re; w.im = -v.im; break; default: w.re = v.im; w.im = -v.re; }
            out[L ^ xL][R ^ xR] = w;
        }
    }
}
static void addrot(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD], int k, c60 out[MAXD][MAXD])
{
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) {
        const c60 a = A[L][R], b = B[L][R]; c60 r;
        switch (k & 3) { case 0: r = b; break; case 1: r.re = -b.im; r.im = b.re; break; case 2: r.re = -b.re; r.im = -b.im; break; default: r.re = b.im; r.im = -b.re; }
        out[L][R].re = a.re + r.re; out[L][R].im = a.im + r.im;
    }
}
/* out = (2 M_ab) in on the left mode: 2M_00 = 1 - i M1 M2, 2M_11 = 1 + i M1 M2, 2M_10 = M1 - i M2 (2 chi^dag), 2M_01 = M1 + i M2 (2 chi) */
static void apply2M(int a, int b, c60 in[MAXD][MAXD], c60 out[MAXD][MAXD])
{
    if (a == b) { papply(&ML[PB], in, XS[0]); papply(&ML[PA], XS[0], XS[1]); addrot(in, XS[1], a ? 1 : 3, out); }
    else { papply(&ML[PA], in, XS[0]); papply(&ML[PB], in, XS[1]); addrot(XS[0], XS[1], a ? 3 : 1, out); }
}
/* (re, im) += i^k <A|B>, exact in 128 bits */
static void inner_acc(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD], int k, s128 *re, s128 *im)
{
    s128 pr = {0, 0}, pi = {0, 0};
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) {
        s128_addp(&pr, A[L][R].re, B[L][R].re); s128_addp(&pr, A[L][R].im, B[L][R].im);
        s128_addp(&pi, A[L][R].re, B[L][R].im); s128_addp(&pi, -A[L][R].im, B[L][R].re);
    }
    switch (k & 3) {
        case 0: s128_adds(re, pr); s128_adds(im, pi); break;
        case 1: s128_subs(re, pi); s128_adds(im, pr); break;
        case 2: s128_subs(re, pr); s128_subs(im, pi); break;
        default: s128_adds(re, pi); s128_subs(im, pr); break;
    }
}
static void build_I(c60 Psi[MAXD][MAXD])
{
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { Psi[L][R].re = Psi[L][R].im = 0; }
    for (int L = 0; L < D; L++) {
        int e = 0, pc = 0;
        for (int k = 0; k < H; k++) if ((L >> k) & 1) { e += (k & 1) ? 1 : 3; pc++; }
        e = (e + 2 * (pc * (pc - 1) / 2)) & 3; EOFL[L] = e;
        const int R = (D - 1) ^ L;
        /* |I> NORMALIZED EXACTLY, so that |I|^2 = 1 and no sum of squares can reach 2^127 (the unnormalized D x 2^120 did at D = 128):
           H even: amplitude i^e 2^{-H/2}; H odd: i^e (1 + i) 2^{-(H+1)/2} (a global phase (1+i)/sqrt2, the review's trick, changes no expectation value) */
        if ((H & 1) == 0) {
            const int64_t a = ONE60 >> (H / 2);
            switch (e) { case 0: Psi[L][R].re = a; break; case 1: Psi[L][R].im = a; break; case 2: Psi[L][R].re = -a; break; default: Psi[L][R].im = -a; }
        } else {
            const int64_t a = ONE60 >> ((H + 1) / 2);
            switch (e) { case 0: Psi[L][R].re = a; Psi[L][R].im = a; break; case 1: Psi[L][R].re = -a; Psi[L][R].im = a; break;
                         case 2: Psi[L][R].re = -a; Psi[L][R].im = -a; break; default: Psi[L][R].re = a; Psi[L][R].im = -a; }
        }
    }
}
static uint64_t rng;
static uint64_t xs64(void) { rng ^= rng >> 12; rng ^= rng << 25; rng ^= rng >> 27; return rng * 0x2545F4914F6CDD1Dull; }
static q gauss40(void) { q s = 0; for (int i = 0; i < 12; i++) s += (q)(xs64() >> 24); return s - 6 * ONE; }
/* add coef60 x (the bare 4-product Pauli string) into the side's exact integer matrix */
static void addterm(pstr P, int64_t coef, int side)
{
    const uint32_t x = side ? P.x >> H : P.x, z = side ? P.z >> H : P.z;
    if (side && ((P.x & (uint32_t)(D - 1)) || (P.z & (uint32_t)(D - 1)))) { printf("WHSYK4: a right term touches the left half -- stop\n"); exit(1); }
    if (!side && ((P.x >> H) || (P.z >> H))) { printf("WHSYK4: a left term touches the right half -- stop\n"); exit(1); }
    c60 (*Hm)[MAXD] = side ? HR : HL;
    for (int b = 0; b < D; b++) {
        const int rot = (P.ph + 2 * parity32(z & (uint32_t)b)) & 3; const int bp = b ^ (int)x;
        switch (rot) { case 0: Hm[bp][b].re += coef; break; case 1: Hm[bp][b].im += coef; break; case 2: Hm[bp][b].re -= coef; break; default: Hm[bp][b].im -= coef; }
    }
}
static void build_syk(uint64_t seed)
{
    rng = seed ? seed : 0x9E3779B97F4A7C15ull;
    const q var = divq(6 * ONE, (q)N * N * N * ONE), sd = sqrtq(var);
    memset(HL, 0, sizeof HL); memset(HR, 0, sizeof HR);
    int nt = 0; q first = 0;
    for (int i = 1; i <= N; i++) for (int j = i + 1; j <= N; j++) for (int k = j + 1; k <= N; k++) for (int l = k + 1; l <= N; l++) {
        const q J = mulq(gauss40(), sd), c = -rdiv(J, 4); if (!nt) first = c;
        const int64_t c60v = (int64_t)c << 20;                 /* the Q.40 coefficient as an exact Q.60 integer */
        addterm(pmul(pmul(pmul(ML[i], ML[j]), ML[k]), ML[l]), c60v, 0);
        addterm(pmul(pmul(pmul(MR[i], MR[j]), MR[k]), MR[l]), c60v, 1);
        nt++;
    }
    NM = N; { char fb[32]; fmt(first, fb); printf("WHSYK4 seed %llu: %d terms, first coefficient %s (v1, v2, v3 print the same two numbers: the same couplings)\n", (unsigned long long)seed, nt, fb); }
}
/* MODEL 1 of the 2022 experiment (Kobrin-Schuster-Yao eq. 1: 7 Majoranas, 5 terms), the Q.40 coefficients v1-v3 declared, shifted to Q.60 */
static void build_model1(void)
{
    static const int idx[5][4] = { {1, 2, 4, 5}, {1, 3, 4, 7}, {1, 3, 5, 6}, {2, 3, 4, 6}, {2, 3, 5, 7} };
    static const double cf[5] = { -0.36, 0.19, -0.71, 0.22, 0.49 };
    memset(HL, 0, sizeof HL); memset(HR, 0, sizeof HR);
    for (int t = 0; t < 5; t++) {
        const double c = cf[t] * (double)ONE / 4.0; const q cq40 = (q)(c < 0 ? c - 0.5 : c + 0.5);
        addterm(pmul(pmul(pmul(ML[idx[t][0]], ML[idx[t][1]]), ML[idx[t][2]]), ML[idx[t][3]]), (int64_t)cq40 << 20, 0);
        addterm(pmul(pmul(pmul(MR[idx[t][0]], MR[idx[t][1]]), MR[idx[t][2]]), MR[idx[t][3]]), (int64_t)cq40 << 20, 1);
    }
    NM = 7;
    printf("WHSYK4 model 1: 5 terms on 7 Majoranas a side (coefficients -0.36, 0.19, -0.71, 0.22, 0.49 declared in Q.40 as v1 did), V normalized with N = 10\n");
}
/* the two series from one set of powers: W = e^{-i dt H} and E = e^{-(beta/4) H}; bq60 = beta/4, dt60 = dt in Q.60 (dt = 1 is exact) */
static void series(c60 Hm[MAXD][MAXD], int64_t bq60, int64_t dt60, c60 W[MAXD][MAXD], c60 E[MAXD][MAXD], c60 P[MAXD][MAXD], c60 u[MAXD][MAXD])
{
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
        SWre[i][j].hi = SWre[i][j].lo = 0; SWim[i][j].hi = SWim[i][j].lo = 0; SEre[i][j].hi = SEre[i][j].lo = 0; SEim[i][j].hi = SEim[i][j].lo = 0;
        P[i][j].re = P[i][j].im = 0;
    }
    for (int i = 0; i < D; i++) { P[i][i].re = ONE60; s128_add(&SWre[i][i], 0, 0); SWre[i][i].hi = ONE60 >> 4; SWre[i][i].lo = (uint64_t)ONE60 << 60; SEre[i][i] = SWre[i][i]; }
    /* the powers are of s H with s = max(dt, beta/4), so that W's coefficient (dt/s)^k and E's (-beta/(4s))^k are both ratios of modulus <= 1
       and nothing overflows the Q.60 ceiling of 8 (dt^k alone passed it at k = 3 for dt = 2.8); with dt = beta/4 = 1 every scale is exactly 1 */
    const int64_t s60 = dt60 > bq60 ? dt60 : bq60, wr60 = div60(dt60, s60), br60 = div60(bq60, s60);
    int64_t bk = ONE60, wk = ONE60;
    int small = 0, k;
    for (k = 1; k <= 80; k++) {
        mmul(Hm, 0, P, u);                               /* u = H P_{k-1}, one rounding per entry */
        bk = -mul60(bk, br60); wk = mul60(wk, wr60);
        int64_t mx = 0;
        for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
            const c60 p = { rdiv64(mul60(u[i][j].re, s60), k), rdiv64(mul60(u[i][j].im, s60), k) }; P[i][j] = p;   /* P_k = (s H)^k / k! */
            /* S_W += (-i)^k dt^k p: k mod 4 = 0: +p; 1: -i p = (p.im, -p.re); 2: -p; 3: +i p = (-p.im, p.re); times wk (exact when dt = 1) */
            const int64_t pr = mul60(p.re, wk), pi = mul60(p.im, wk);
            int64_t wr, wi; switch (k & 3) { case 0: wr = pr; wi = pi; break; case 1: wr = pi; wi = -pr; break; case 2: wr = -pr; wi = -pi; break; default: wr = -pi; wi = pr; }
            s128_add(&SWre[i][j], wr >> 4, (uint64_t)wr << 60); s128_add(&SWim[i][j], wi >> 4, (uint64_t)wi << 60);
            /* S_E += bk p (a Q.120 product; exact when bk = +-1) */
            s128_addp(&SEre[i][j], bk, p.re); s128_addp(&SEim[i][j], bk, p.im);
            const int64_t a = p.re < 0 ? -p.re : p.re, b = p.im < 0 ? -p.im : p.im; if (a > mx) mx = a; if (b > mx) mx = b;
        }
        if (mx <= 1) { if (++small >= 2) break; } else small = 0;
    }
    if (k > taylor_max) taylor_max = k;
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
        W[i][j].re = s128_round_sh(SWre[i][j], F60); W[i][j].im = s128_round_sh(SWim[i][j], F60);
        E[i][j].re = s128_round_sh(SEre[i][j], F60); E[i][j].im = s128_round_sh(SEim[i][j], F60);
    }
}
/* the coupling: Psi <- PROD_j (c + i ss W_j) Psi, W_j = i bare_l^j bare_r^j; one rounding per entry per j */
static void apply_U(c60 Psi[MAXD][MAXD], int64_t c, int64_t s, int sign, c60 w[MAXD][MAXD])
{
    for (int j = 1; j <= NM; j++) {
        pstr W = pmul(ML[j], MR[j]); W.ph = (W.ph + 1) & 3;
        papply(&W, Psi, w);
        const int64_t ss = sign > 0 ? s : -s;
        for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) {
            s128 re = {0, 0}, im = {0, 0};
            s128_addp(&re, Psi[L][R].re, c); s128_addp(&re, -w[L][R].im, ss);
            s128_addp(&im, Psi[L][R].im, c); s128_addp(&im, w[L][R].re, ss);
            Psi[L][R].re = s128_round_sh(re, F60); Psi[L][R].im = s128_round_sh(im, F60);
        }
    }
}
static s128 norm2_128(c60 Psi[MAXD][MAXD]) { s128 a = {0, 0}; for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { s128_addp(&a, Psi[L][R].re, Psi[L][R].re); s128_addp(&a, Psi[L][R].im, Psi[L][R].im); } return a; }
static void normalize(c60 Psi[MAXD][MAXD])
{
    const int64_t n = sqrt120(norm2_128(Psi));
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { Psi[L][R].re = div60(Psi[L][R].re, n); Psi[L][R].im = div60(Psi[L][R].im, n); }
}
static int64_t maxabs_diff(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD]) { int64_t m = 0; for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { int64_t a = A[L][R].re - B[L][R].re, b = A[L][R].im - B[L][R].im; if (a < 0) a = -a; if (b < 0) b = -b; if (a > m) m = a; if (b > m) m = b; } return m; }
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
    N = 8; int vnormN = 0, model = 0; q beta = 4 * ONE, dt = ONE; q mus[16], ts[32], tps[32]; int nmu = 0, nt = 0, ntp = 0; uint64_t seeds[16]; int nseed = 0;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "model=", 6)) model = !strcmp(argv[i] + 6, "m1");
        else if (!strncmp(argv[i], "dt=", 3)) dt = parseq(argv[i] + 3);          /* the time unit: W(dt) by its own series; t and t' count units */
        else if (!strncmp(argv[i], "beta=", 5)) beta = parseq(argv[i] + 5);
        else if (!strncmp(argv[i], "mu=", 3)) nmu = parselist(argv[i] + 3, mus, 16);
        else if (!strncmp(argv[i], "t=", 2)) nt = parselist(argv[i] + 2, ts, 32);
        else if (!strncmp(argv[i], "tp=", 3)) ntp = parselist(argv[i] + 3, tps, 32);
        else if (!strncmp(argv[i], "pair=", 5)) sscanf(argv[i] + 5, "%d,%d", &PA, &PB);
        else if (!strncmp(argv[i], "vnorm=", 6)) vnormN = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "seeds=", 6)) { char b[256]; strncpy(b, argv[i] + 6, 255); b[255] = 0; for (char *t = strtok(b, ","); t && nseed < 16; t = strtok(NULL, ",")) seeds[nseed++] = strtoull(t, 0, 10); }
    }
    if (N < 4 || N > MAXN || (N & 1)) { printf("WHSYK4: N must be even, 4..%d\n", MAXN); return 1; }
    if (!nmu) { mus[0] = -12 * ONE; mus[1] = 12 * ONE; nmu = 2; }
    if (!nt) { ts[0] = ONE; ts[1] = 2 * ONE; ts[2] = 3 * ONE; nt = 3; }
    if (!nseed) { seeds[0] = 1; nseed = 1; }
    if (model) { N = 8; nseed = 1; if (!vnormN) vnormN = 10; }
    if (!vnormN) vnormN = N;
    if (dt <= 0 || dt > 4 * ONE) { printf("WHSYK4: dt must be in (0, 4]\n"); return 1; }
    if (PA < 1 || PA > N || PB < 1 || PB > N || PA == PB) { printf("WHSYK4: pair must be two different Majorana indices in 1..N\n"); return 1; }
    int tmax = 0; for (int i = 0; i < nt; i++) { int t = (int)(ts[i] / ONE); if ((q)t * ONE != ts[i] || t < 0 || t > TMAX) { printf("WHSYK4: t must be an integer in 0..%d\n", TMAX); return 1; } if (t > tmax) tmax = t; }
    for (int i = 0; i < ntp; i++) { int t = (int)(tps[i] / ONE); if ((q)t * ONE != tps[i] || t < 0 || t > TMAX) { printf("WHSYK4: tp must be an integer in 0..%d\n", TMAX); return 1; } if (t > tmax) tmax = t; }
    H = N / 2; D = 1 << H;
    const int64_t beta60 = (int64_t)beta << 20, bq60 = rdiv64(beta60, 4), dt60 = (int64_t)dt << 20;
    LN2_60 = artanh2_60(div60(ONE60, 3 * ONE60));        /* ln 2 = 2 artanh(1/3) */
    char A[48], B[48], C[48], Dd[48], E[48], F[48], G[48], S1[48], S2[48], S3[48], S4[48], S5[48], S6[48], S7[48], S8[48];
    fmt60(beta60, A); fmt60(LN2_60, B);
    fmt60(dt60, C);
    printf("WHSYK4: Q.%d end to end (int64 with 128-bit products; nothing rounded to Q.40); %s; N = %d Majoranas a side, FACTORIZED (D = %d a side); beta = %s; V normalized with N = %d; the Dirac modes chi = (psi^%d + i psi^%d)/sqrt2; the time unit dt = %s (W(dt) by its own series, W(t) = W(dt)^t) and e^{-beta H/4} from one set of powers per seed; the EXACT rho_TR from its six elements (GJ (2.15)-(2.20)); K = <{A1, B1}>\n",
           F60, model ? "MODEL 1 of the 2022 experiment (Kobrin-Schuster-Yao eq. 1: 7 Majoranas, 5 terms)" : "the FULL q = 4 SYK, <J^2> = 6/N^3, J = 1", N, D, A, vnormN, PA, PB, C);
    printf("WHSYK4 selftest: ln 2 = %s (0.693147180560)\n", B);
    build_majoranas();
    build_I(PsiI);
    int64_t worst = 0;
    for (int j = 1; j <= N; j++) {
        papply(&ML[j], PsiI, WA); papply(&MR[j], PsiI, WB);
        for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { int64_t re = WA[L][R].re - WB[L][R].im, im = WA[L][R].im + WB[L][R].re; if (re < 0) re = -re; if (im < 0) im = -im; if (re > worst) worst = re; if (im > worst) worst = im; }
    }
    printf("WHSYK4 selftest 1: max |(psi_l^j + i psi_r^j)|I>| = %lld units of 2^-60 (must be 0)\n", (long long)worst);
    for (int si = 0; si < nseed; si++) {
        if (model) build_model1(); else build_syk(seeds[si]);
        /* self-test 2: H_L Psi_I = Psi_I H_R^T exactly */
        mmul(HL, 0, PsiI, WA); transpose(PsiI, WX); mmul(HR, 0, WX, WB); transpose(WB, WC);
        printf("WHSYK4 seed %llu selftest 2: max |H_L Psi_I - Psi_I H_R^T| = %lld units of 2^-60 (must be 0)\n", (unsigned long long)seeds[si], (long long)maxabs_diff(WA, WC));
        hash64(maxabs_diff(WA, WC));
        /* the series: W(1) and E on each side; then W(t) = W(1) W(t-1) */
        series(HL, bq60, dt60, WL[1], EL, WA, WB); series(HR, bq60, dt60, WR[1], ER, WA, WB);
        for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) { WL[0][i][j].re = WR[0][i][j].re = (i == j) ? ONE60 : 0; WL[0][i][j].im = WR[0][i][j].im = 0; }
        for (int t = 2; t <= tmax; t++) { mmul(WL[1], 0, WL[t - 1], WL[t]); mmul(WR[1], 0, WR[t - 1], WR[t]); }
        /* self-test 5: THE RIGHT SIDE IS THE LEFT RE-INDEXED. H_L Psi_I = Psi_I H_R^T exactly and Psi_I is a phased permutation (one entry per row,
           on R = ~L, phase i^e(L) up to a global factor), so H_R^T = Psi^-1 H_L Psi and the same for every function of H: W_R(t)[a][b] =
           i^{e(~a) - e(~b)} W_L(t)[~b][~a] and E_R likewise, up to the roundings: the left series multiplies H P_{k-1} and the re-indexed
           right series amounts to P_{k-1} H, equal in exact arithmetic (the powers commute) but rounded differently by a unit or two per
           entry (MEASURED 11:08 PM: 1-3 units at N = 8-14, exactly 0 for the 2022 model). A left recursion written as P_{k-1} H would make
           the identity bit-exact and the right-side series redundant: half the series work and half the stored step operators (the memory
           that stands between this engine and N = 18); that is version 4.1's change. Not hashed: the result chain is unchanged by this test. */
        { int64_t w5 = 0, e5 = 0;
          for (int a = 0; a < D; a++) for (int b = 0; b < D; b++) {
              const int na = (D - 1) ^ a, nb = (D - 1) ^ b, rot = (EOFL[na] - EOFL[nb]) & 3;
              c60 v = WL[1][nb][na], w; c60 ve = EL[nb][na], we;
              switch (rot) { case 0: w = v; we = ve; break; case 1: w.re = -v.im; w.im = v.re; we.re = -ve.im; we.im = ve.re; break;
                             case 2: w.re = -v.re; w.im = -v.im; we.re = -ve.re; we.im = -ve.im; break; default: w.re = v.im; w.im = -v.re; we.re = ve.im; we.im = -ve.re; }
              int64_t d1 = WR[1][a][b].re - w.re, d2 = WR[1][a][b].im - w.im, d3 = ER[a][b].re - we.re, d4 = ER[a][b].im - we.im;
              if (d1 < 0) d1 = -d1; if (d2 < 0) d2 = -d2; if (d3 < 0) d3 = -d3; if (d4 < 0) d4 = -d4;
              if (d1 > w5) w5 = d1; if (d2 > w5) w5 = d2; if (d3 > e5) e5 = d3; if (d4 > e5) e5 = d4;
          }
          printf("WHSYK4 seed %llu selftest 5: |W_R(1) - the re-indexed W_L(1)| = %lld units, |E_R - the re-indexed E_L| = %lld units of 2^-60 (a few units: the two series multiply in opposite orders; the identity is exact in exact arithmetic)\n",
                 (unsigned long long)seeds[si], (long long)w5, (long long)e5); }
        /* self-test 3: the unitarity defect of W_L(1): max |W^dag W - 1| in units of 2^-60 */
        mmul(WL[1], 1, WL[1], WA); for (int i = 0; i < D; i++) WA[i][i].re -= ONE60;
        memset(WB, 0, sizeof(c60) * MAXD * MAXD); const int64_t udef = maxabs_diff(WA, WB); if (udef > unit_defect_max) unit_defect_max = udef;
        /* self-test 4: E_L Psi_I = Psi_I E_R^T to the series' rounding */
        mmul(EL, 0, PsiI, WA); transpose(PsiI, WX); mmul(ER, 0, WX, WB); transpose(WB, WC);
        const int64_t edef = maxabs_diff(WA, WC);
        printf("WHSYK4 seed %llu selftest 3-4: |W_L(1)^dag W_L(1) - 1| <= %lld units of 2^-60; |e^{-beta H_L/4} Psi_I - Psi_I (e^{-beta H_R/4})^T| <= %lld units of 2^-60 (both must be a few units); Taylor terms at most %d\n",
               (unsigned long long)seeds[si], (long long)udef, (long long)edef, taylor_max);
        hash64(udef); hash64(edef);
        /* the thermofield double: E_L Psi_I E_R^T, normalized in Q.60 */
        copy(PsiI, T); sandwich(EL, ER, 0, T);
        { const s128 n2 = norm2_128(T); fmt60(s128_round_sh(n2, F60), B); }     /* |E I|^2 with |I|^2 = 1 exactly: Z_beta / Z_inf */
        normalize(T);
        fmt60(sqrt120(norm2_128(T)), A);
        printf("WHSYK4 seed %llu: |e^{-beta H/4} I|^2 / |I|^2 = %s (Z_beta / Z_inf, the number v2 and v3 print to their digits); the thermofield double normalized to %s\n", (unsigned long long)seeds[si], B, A);
        for (int mi = 0; mi < nmu; mi++) {
            /* arg = mu / (8 vnorm) in Q.60 from the Q.40 mu: round(mu40 2^20 / (8 vnorm)) */
            int64_t arg60; { const q m = mus[mi]; int neg = m < 0; uint64_t um = neg ? -(uint64_t)m : (uint64_t)m, h, l, rem; mul128(um, (uint64_t)1 << 20, &h, &l); uint64_t qq = div128(h, l, (uint64_t)(8 * vnormN), &rem); if (2 * rem >= (uint64_t)(8 * vnormN)) qq++; arg60 = neg ? -(int64_t)qq : (int64_t)qq; }
            int64_t c, s; cossin60(arg60, &c, &s);
            const int64_t csdef = mul60(c, c) + mul60(s, s) - ONE60;
            for (int ti = 0; ti < nt; ti++) {
                const int t = (int)(ts[ti] / ONE), tp = ntp ? (int)(tps[ti < ntp ? ti : ntp - 1] / ONE) : t;
                /* THE EXACT MATRIX: WG = e^{iHt'} T; Z_ab = e^{-iHt} U e^{-iHt'} (2 M_ab) WG */
                copy(T, WG); evolve(WG, tp, 1);
                for (int ab = 0; ab < 4; ab++) { apply2M(ab >> 1, ab & 1, WG, ZS[ab]); evolve(ZS[ab], tp, 0); apply_U(ZS[ab], c, s, 1, WX); evolve(ZS[ab], t, 0); }
                s128 ar[6], ai[6]; memset(ar, 0, sizeof ar); memset(ai, 0, sizeof ai);
                for (int b = 0; b < 2; b++) {
                    c60 (*Z0)[MAXD] = ZS[b], (*Z1)[MAXD] = ZS[2 + b];
                    papply(&MR[PB], Z0, XS[0]); papply(&MR[PA], XS[0], XS[1]);
                    inner_acc(Z0, Z0, 0, &ar[0], &ai[0]); inner_acc(Z0, XS[1], 3, &ar[0], &ai[0]);
                    inner_acc(Z0, Z0, 0, &ar[1], &ai[1]); inner_acc(Z0, XS[1], 1, &ar[1], &ai[1]);
                    papply(&MR[PA], Z0, XS[0]); papply(&MR[PB], Z0, XS[1]);
                    inner_acc(Z1, XS[0], 0, &ar[4], &ai[4]); inner_acc(Z1, XS[1], 3, &ar[4], &ai[4]);
                    inner_acc(Z1, XS[0], 0, &ar[5], &ai[5]); inner_acc(Z1, XS[1], 1, &ar[5], &ai[5]);
                    papply(&MR[PB], Z1, XS[0]); papply(&MR[PA], XS[0], XS[1]);
                    inner_acc(Z1, Z1, 0, &ar[2], &ai[2]); inner_acc(Z1, XS[1], 3, &ar[2], &ai[2]);
                    inner_acc(Z1, Z1, 0, &ar[3], &ai[3]); inner_acc(Z1, XS[1], 1, &ar[3], &ai[3]);
                }
                int64_t rr[6], ri[6]; for (int e = 0; e < 6; e++) { rr[e] = s128_round_sh(ar[e], F60 + 4); ri[e] = s128_round_sh(ai[e], F60 + 4); }
                const int64_t r00 = rr[0], r11 = rr[1], r22 = rr[2], r33 = rr[3];
                const int64_t epr1 = r00 + r11 - (ONE60 >> 1), epr2 = r22 + r33 - (ONE60 >> 1);
                int64_t imdiag = 0; for (int e = 0; e < 4; e++) { int64_t v = ri[e] < 0 ? -ri[e] : ri[e]; if (v > imdiag) imdiag = v; }
                /* the entropies in Q.60 */
                int clamp = 0; int64_t SRT = 0, SR = 0, ST = 0;
                { const int64_t m = rdiv64(r00 + r33, 2), d = rdiv64(r00 - r33, 2); s128 rad = {0, 0}; s128_addp(&rad, d, d); s128_addp(&rad, rr[4], rr[4]); s128_addp(&rad, ri[4], ri[4]);
                  const int64_t disc = sqrt120(rad); SRT += ent60(m + disc, &clamp) + ent60(m - disc, &clamp); }
                { const int64_t m = rdiv64(r11 + r22, 2), d = rdiv64(r11 - r22, 2); s128 rad = {0, 0}; s128_addp(&rad, d, d); s128_addp(&rad, rr[5], rr[5]); s128_addp(&rad, ri[5], ri[5]);
                  const int64_t disc = sqrt120(rad); SRT += ent60(m + disc, &clamp) + ent60(m - disc, &clamp); }
                ST = ent60(r00 + r22, &clamp) + ent60(r11 + r33, &clamp);          /* the readout T's marginal */
                SR = ent60(r00 + r11, &clamp) + ent60(r22 + r33, &clamp);          /* the reference R's: (1/2, 1/2) */
                const int64_t I = SR + ST - SRT, frac = div60(I, 2 * LN2_60);
                /* K = Re <M1(-t') T | U^dag R1(t) U T> */
                papply(&ML[PA], WG, WA); evolve(WA, tp, 0);                      /* a1 */
                copy(T, WB); apply_U(WB, c, s, 1, WX); evolve(WB, t, 0); papply(&MR[PA], WB, WC); evolve(WC, t, 1); apply_U(WC, c, s, -1, WX);   /* b1 */
                s128 kre = {0, 0}, kim = {0, 0}; inner_acc(WA, WC, 0, &kre, &kim); const int64_t K = s128_round_sh(kre, F60);
                fmt60((int64_t)t * dt60, A); fmt(mus[mi], B); fmt60(r00, C); fmt60(r11, Dd); fmt60(r22, E); fmt60(r33, F); fmt60(rr[4], G); fmt60(ri[4], S1); fmt60(rr[5], S2); fmt60(ri[5], S3);
                fmt60(K, S7); fmt60(frac, S8); fmt60((int64_t)tp * dt60, S4); fmt60(beta60, S5); fmt60(epr1, S6);
                printf("WHSYK4 N=%d seed=%llu beta=%s mu=%s t=%s tp=%s: EXACT rho00 %s rho11 %s rho22 %s rho33 %s | rho03 %s %s i | rho12 %s %s i | checks (R's marginal) rho00+rho11-1/2 %s rho22+rho33-1/2 %lld units, max|Im diag| %lld units, c^2+s^2-1 %lld units of 2^-60 | K %s | I_RT/(2 ln 2) EXACT %s%s\n",
                       N, (unsigned long long)seeds[si], S5, B, A, S4, C, Dd, E, F, G, S1, S2, S3, S6, (long long)epr2, (long long)imdiag, (long long)csdef, S7, S8, clamp ? " (an eigenvalue clamped at 0)" : "");
                hash64(r00); hash64(r11); hash64(r22); hash64(r33); hash64(rr[4]); hash64(ri[4]); hash64(rr[5]); hash64(ri[5]); hash64(K); hash64(frac);
            }
        }
    }
    printf("WHSYK4: unitarity defect of W(1) over the seeds at most %lld units of 2^-60; Taylor terms at most %d\n", (long long)unit_defect_max, taylor_max);
    printf("WHSYK4: chain %016llx\n", (unsigned long long)chain);
    printf("WHSYK4: done\n");
    return 0;
}
