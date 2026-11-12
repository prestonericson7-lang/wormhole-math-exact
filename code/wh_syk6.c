/* wh_syk6.c (11-12-2026, 12:45 AM) -- WORMHOLE MATH 2, v6: v5.2 (the exact six-element rho_TR in Q.60 end to end, the 17-matrix pool, the
   two-core products) PLUS THE SIZE-WINDING OBSERVABLE of Kobrin, Schuster and Yao, derived for this layout by the review (item 7, 12:17 AM)
   and computed here from the states v5 already holds, in the same arithmetic. The protocol lines (WHSYK6 N=...: EXACT ...) and their chain
   are v5's bit for bit at dt = 1; the new lines carry their own chain ("winding chain"). One more change: dt= is parsed straight into Q.60
   (round(dt 2^60)); v4 and v5 parsed it through Q.40, so their dt = 2.8 was 2.8 + 1.8e-13 (the review's finding; it moved K by ~3e-14).

   THE OBSERVABLE. For the operator state A = e^{-iHt} g_i e^{iHt} T (sqrt2 times psi_l^i(-t)|tfd> in this layout; ||A||_F = 1):
   - O' = A Pi^dag, i.e. O'[c][L] = A[c][L ^ (D-1)] i^{-e(L)}: the entangled state's phased permutation undone (exact quarter turns).
   - For each x: f_x(c) = O'[c][c ^ x]; F_x(z) = SUM_c (-1)^{popcount(z & c)} f_x(c), a Walsh-Hadamard transform (additions only);
     T_{x,z} = i^{popcount(x & z)} F_x(z) = tr(sigma_{x,z} O'), with sigma_{x,z} the Hermitian Pauli string (X on x\z, Y on x&z, Z on z\x),
     which under the Jordan-Wigner map is +-1 times one Majorana string. Its Majorana size from (x, z): run k from n-1 down to 0 with a
     parity bit p = 0; X or Y on qubit k counts 1, Z counts 2 if p = 0, I counts 2 if p = 1; p flips after each odd count.
   - P(l) = SUM_{|P| = l} |T|^2 / D and q(l) = omega^-2 SUM_{|P| = l} T^2 / D, omega the global phase of the normalized |I> ((1+i)/sqrt2 for
     n odd, so omega^-2 = -i; 1 for n even): KSY's p(l) and q(l) for rho^{1/2} psi_i(t), checked by the review against KSY's own operator.
   ARITHMETIC: each entry of O' is rounded once from Q.60 to Q.56; the transform is int64 additions (|F| <= sqrt(D) ||f||_2 by Cauchy-
   Schwarz, below 2^61 at D = 512); T^2 and |T|^2 are exact 128-bit products summed in 128 bits; the per-size sums are rounded once to Q.60
   at the output (shift 52 + n). FREE CHECKS, printed: the even sizes are exactly 0 (a real check only with windcheck=1, which sums the even-popcount
   diagonals into them; without it an identity, since a string's size has the parity of popcount(x)); SUM_l P(l) - 1 in units; W = SUM_l q(l) (real, t-
   independent); with windcheck=1 the diagonals of even popcount(x) are computed too, must be exactly 0, and feed the even sizes (otherwise
   skipped). Reviewed line by line by the reviewer (1:05 AM): the four functions right; these two text points were its findings.
   THE LINK TO THE PROTOCOL (an exact identity): C(mu) = <tfd| psi_r^i(t) e^{i mu V} psi_l^i(-t) |tfd> = -(i/2) e^{-i mu NM/(8 vnorm)}
   SUM_l q(l) e^{i mu l/(4 vnorm)}, because e^{i mu V} multiplies a size-l string by exactly e^{i mu (l - NM/2)/(4 vnorm)}. So 2|C(mu)| is
   KSY's line-fit function at m = -mu/(4 vnorm): the sign of mu that unwinds arg q(l) is the teleporting sign. v6 prints C(mu) from q(l)
   for every mu (the phases e^{i 2 l arg} as powers of the small rotation e^{i arg}, arg = mu/(8 vnorm), one rounding a step: the Q.60
   series itself only takes angles below 2.8); windcheck=1 also computes C(mu) from the states (<B| e^{i mu V} |A>/2 with B = e^{iHt} g_i^R
   e^{-iHt} T) and prints the difference in units of 2^-60. K itself is not a function of q(l) (U stands on both sides of it); it stays.
   COST: the transform is D^2 log D additions and 3 D^2 products per operator and time, under 0.1 percent of a protocol point; the state A
   of the first operator is the one K already holds; each further operator costs two products (and the state check four more).
   Command line (v5's, plus): wind=0|1|2 (0 off; 1 the pair's two operators, the default; 2 all operators) windcheck=0|1.
       whsyk6 N=14 seeds=1 beta=4 mu=-12,12 dt=2.8 t=1 wind=1 windcheck=1     (KSY's t = 2.8; their t = 4.0 by dt=1 t=4: the pool
       holds W(1..6), so t = 2.8 cannot be reached as seven steps of 0.4)
   THE LAYOUT (v5): one static pool of 17 matrices (68 MB at D = 512), carved by use: W_L(1..6), E_L, PsiI, T, WG, WX, Z_00..Z_11, XS[0..1];
   the winding uses Z_01 (A of the first operator, K's state), Z_10 (A of the others), XS[0] (O'), and with windcheck Z_00/XS[1] (B), Z_11
   (e^{i mu V} A). No right-hand operator is stored (mmul modes 2 and 3 read the left one through the exact re-indexing). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "whq.h"
#include "par2.h"                           /* v5.2 (11:44 PM): the A9's second core takes half the rows of every product (the machines' two-core call); on a
                                                  machine without it, or on the PC, both halves run here one after the other: the same bits either way */

#define MAXN 18
#define MAXD (1 << (MAXN / 2))
#define TMAX 6
#define F60 60
static const int64_t ONE60 = (int64_t)1 << F60;

typedef struct { int64_t re, im; } c60;
typedef struct { int64_t hi; uint64_t lo; } s128;
typedef struct { uint32_t x, z; int ph; } pstr;
typedef c60 Mat[MAXD][MAXD];
typedef s128 Acc[MAXD][MAXD];

static int N, H, D, NM, PA = 1, PB = 2;
static pstr ML[MAXN + 1], MR[MAXN + 1];
/* THE POOL: 17 matrices */
#define NPOOL 17
static c60 POOL[NPOOL][MAXD][MAXD];
#define WL(t) (POOL[(t) - 1])                 /* slots 0..5: W_L(1..6) */
#define EL (POOL[6])
#define PSII (POOL[7])
#define TT (POOL[8])
#define WG (POOL[9])
#define WX (POOL[10])
#define ZS(k) (POOL[11 + (k)])               /* slots 11..14: Z_00, Z_01, Z_10, Z_11 */
#define XS(k) (POOL[15 + (k)])               /* slots 15..16 */
#define HL (POOL[5])                          /* H_L in W_L(6)'s slot: dead after the series */
#define HR (POOL[4])                          /* H_R (check=1 only) in W_L(5)'s slot: dead after the series */
#define SWRE ((Acc *)POOL[11])                /* the series' accumulators in the Z slots */
#define SWIM ((Acc *)POOL[12])
#define SERE ((Acc *)POOL[13])
#define SEIM ((Acc *)POOL[14])
static uint64_t chain = 0xcbf29ce484222325ull;
static int taylor_max = 0; static int64_t unit_defect_max = 0;
static int64_t LN2_60;
static int EOFL[MAXD];

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
static inline int64_t s128_round_sh(s128 a, int sh)
{
    int neg = a.hi < 0; if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t r = (a.lo >> sh) | ((uint64_t)a.hi << (64 - sh)); if ((a.lo & ((((uint64_t)1) << sh) - 1)) >= (((uint64_t)1) << (sh - 1))) r++;
    return neg ? -(int64_t)r : (int64_t)r;
}
static inline int64_t mul60(int64_t a, int64_t b) { s128 s = {0, 0}; s128_addp(&s, a, b); return s128_round_sh(s, F60); }
static inline int64_t rdiv64(int64_t x, int64_t n) { return x >= 0 ? (x + n / 2) / n : -((-x + n / 2) / n); }
static int64_t div60(int64_t a, int64_t b)
{
    int neg = a < 0; uint64_t ua = neg ? -(uint64_t)a : (uint64_t)a;
    const uint64_t hi = ua >> 4, lo = ua << 60;
    if (hi >= (uint64_t)b) { printf("WHSYK6: div60 overflow -- stop\n"); exit(1); }
    uint64_t rem; uint64_t qq = div128(hi, lo, (uint64_t)b, &rem); if (2 * rem >= (uint64_t)b) qq++;
    return neg ? -(int64_t)qq : (int64_t)qq;
}
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
static int64_t sqrt120(s128 a) { if (a.hi < 0) return 0; return (int64_t)isqrt128((uint64_t)a.hi, a.lo); }
static int64_t artanh2_60(int64_t z)
{
    const int64_t z2 = mul60(z, z); int64_t term = z, sum = 0;
    for (int k = 0; k < 40; k++) { const int64_t add = rdiv64(term, 2 * k + 1); if (!add && k > 2) break; sum += add; term = mul60(term, z2); if (!term) break; }
    return 2 * sum;
}
static int64_t ln60(int64_t x)
{
    if (x <= 0) { printf("WHSYK6: ln of a non-positive value -- stop\n"); exit(1); }
    int k = 0; int64_t m = x;
    while (m < (ONE60 >> 1)) { m <<= 1; k--; }
    while (m >= ONE60) { m = rdiv64(m, 2); k++; }
    const int64_t z = div60(m - ONE60, m + ONE60);
    return artanh2_60(z) + (int64_t)k * LN2_60;
}
static int64_t ent60(int64_t lam, int *clamp) { if (lam <= 0) { *clamp = 1; return 0; } return -mul60(lam, ln60(lam)); }
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
/* the entry op(A)[i][k]: mode 0 A[i][k]; 1 conj(A[k][i]); 2 the re-indexed right-hand operator R[i][k] = i^{e(~i) - e(~k)} A[~k][~i];
   3 conj(R[k][i]) = conj(i^{e(~k) - e(~i)} A[~i][~k]) */
static inline c60 entry(c60 A[MAXD][MAXD], int mode, int i, int k)
{
    c60 v;
    if (mode == 0) return A[i][k];
    if (mode == 1) { v = A[k][i]; v.im = -v.im; return v; }
    const int ni = (D - 1) ^ i, nk = (D - 1) ^ k; int rot; c60 w;
    if (mode == 2) { v = A[nk][ni]; rot = (EOFL[ni] - EOFL[nk]) & 3; }
    else { v = A[ni][nk]; rot = (EOFL[nk] - EOFL[ni]) & 3; }
    switch (rot) { case 0: w = v; break; case 1: w.re = -v.im; w.im = v.re; break; case 2: w.re = -v.re; w.im = -v.im; break; default: w.re = v.im; w.im = -v.re; }
    if (mode == 3) w.im = -w.im;
    return w;
}
/* C = op(A) B; exact zeros of op(A) skipped; one rounding per entry. TILED (v5.1, 11:37 PM): the output is formed in 32 x 32 tiles whose
   128-bit accumulators stay in the L1 cache while the k index runs in tiles of 32, so that at D = 512 (a 4 MB operand against the A9's
   512 KB L2) the strided operand is read once per tile instead of once per entry. The sums are exact in 128 bits, so the order of the k
   tiles changes no bit: the result equals the untiled product exactly (checked by the chains at N = 8-12). */
#define TILE 32
typedef struct { c60 (*A)[MAXD]; int mode; c60 (*B)[MAXD]; c60 (*C)[MAXD]; int i_begin, i_end; } mjob;
/* the rows i_begin..i_end-1 of C = op(A) B, tiled; pure computation on its own stack (the accumulators are locals), so it can run on
   core 1: it reads A, B, D and EOFL and writes only its own rows of C */
static void mmul_rows(void *p)
{
    const mjob *jb = (const mjob *)p; c60 (*A)[MAXD] = jb->A, (*B)[MAXD] = jb->B, (*C)[MAXD] = jb->C; const int mode = jb->mode;
    s128 acr[TILE][TILE], aci[TILE][TILE];
    const int T = D < TILE ? D : TILE;
    for (int i0 = jb->i_begin; i0 < jb->i_end; i0 += T) for (int j0 = 0; j0 < D; j0 += T) {
        memset(acr, 0, sizeof acr); memset(aci, 0, sizeof aci);
        for (int k0 = 0; k0 < D; k0 += T)
            for (int i = i0; i < i0 + T; i++) for (int k = k0; k < k0 + T; k++) {
                const c60 a = entry(A, mode, i, k);
                if (!a.re && !a.im) continue;
                s128 *re = acr[i - i0], *im = aci[i - i0]; const c60 *b = &B[k][j0];
                for (int j = 0; j < T; j++) {
                    s128_addp(&re[j], a.re, b[j].re); s128_addp(&re[j], -a.im, b[j].im);
                    s128_addp(&im[j], a.re, b[j].im); s128_addp(&im[j], a.im, b[j].re);
                }
            }
        for (int i = 0; i < T; i++) for (int j = 0; j < T; j++) { C[i0 + i][j0 + j].re = s128_round_sh(acr[i][j], F60); C[i0 + i][j0 + j].im = s128_round_sh(aci[i][j], F60); }
    }
}
/* C = op(A) B: the rows split in two halves, one per core when the machine has two (D >= 64; at D = 16, 32 one core does it all) */
static void mmul(c60 A[MAXD][MAXD], int mode, c60 B[MAXD][MAXD], c60 C[MAXD][MAXD])
{
    mjob j0 = { A, mode, B, C, 0, D }, j1 = { A, mode, B, C, D / 2, D };
    if (D >= 2 * TILE) { j0.i_end = D / 2; par2(mmul_rows, &j0, &j1); }
    else mmul_rows(&j0);
}
static void transpose(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD]) { for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) B[j][i] = A[i][j]; }
static void copy(c60 A[MAXD][MAXD], c60 B[MAXD][MAXD]) { for (int i = 0; i < D; i++) memcpy(B[i], A[i], sizeof(c60) * D); }
/* Psi <- op(L) Psi op(R)^T with R the re-indexed L: sign 0: W_L(t) Psi W_R(t)^T; sign 1: W_L^dag Psi conj(W_R) */
static void sandwich(c60 L[MAXD][MAXD], int sign, c60 Psi[MAXD][MAXD])
{
    mmul(L, sign ? 1 : 0, Psi, WX); transpose(WX, Psi); mmul(L, sign ? 3 : 2, Psi, WX); transpose(WX, Psi);
}
static void evolve(c60 Psi[MAXD][MAXD], int t, int sign) { if (t > 0) sandwich(WL(t), sign, Psi); }

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
static void apply2M(int a, int b, c60 in[MAXD][MAXD], c60 out[MAXD][MAXD])
{
    if (a == b) { papply(&ML[PB], in, XS(0)); papply(&ML[PA], XS(0), XS(1)); addrot(in, XS(1), a ? 1 : 3, out); }
    else { papply(&ML[PA], in, XS(0)); papply(&ML[PB], in, XS(1)); addrot(XS(0), XS(1), a ? 3 : 1, out); }
}
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
static void addterm(pstr P, int64_t coef, int side)
{
    const uint32_t x = side ? P.x >> H : P.x, z = side ? P.z >> H : P.z;
    if (side && ((P.x & (uint32_t)(D - 1)) || (P.z & (uint32_t)(D - 1)))) { printf("WHSYK6: a right term touches the left half -- stop\n"); exit(1); }
    if (!side && ((P.x >> H) || (P.z >> H))) { printf("WHSYK6: a left term touches the right half -- stop\n"); exit(1); }
    c60 (*Hm)[MAXD] = side ? HR : HL;
    for (int b = 0; b < D; b++) {
        const int rot = (P.ph + 2 * parity32(z & (uint32_t)b)) & 3; const int bp = b ^ (int)x;
        switch (rot) { case 0: Hm[bp][b].re += coef; break; case 1: Hm[bp][b].im += coef; break; case 2: Hm[bp][b].re -= coef; break; default: Hm[bp][b].im -= coef; }
    }
}
static void zero(c60 A[MAXD][MAXD]) { for (int i = 0; i < D; i++) memset(A[i], 0, sizeof(c60) * D); }
static void build_syk(uint64_t seed, int check)
{
    rng = seed ? seed : 0x9E3779B97F4A7C15ull;
    const q var = divq(6 * ONE, (q)N * N * N * ONE), sd = sqrtq(var);
    zero(HL); zero(HR);
    int nt = 0; q first = 0;
    for (int i = 1; i <= N; i++) for (int j = i + 1; j <= N; j++) for (int k = j + 1; k <= N; k++) for (int l = k + 1; l <= N; l++) {
        const q J = mulq(gauss40(), sd), c = -rdiv(J, 4); if (!nt) first = c;
        const int64_t c60v = (int64_t)c << 20;
        addterm(pmul(pmul(pmul(ML[i], ML[j]), ML[k]), ML[l]), c60v, 0);
        if (check) addterm(pmul(pmul(pmul(MR[i], MR[j]), MR[k]), MR[l]), c60v, 1);
        nt++;
    }
    NM = N; { char fb[32]; fmt(first, fb); printf("WHSYK6 seed %llu: %d terms, first coefficient %s (v1-v4 print the same two numbers: the same couplings)\n", (unsigned long long)seed, nt, fb); }
}
static void build_model1(int check)
{
    static const int idx[5][4] = { {1, 2, 4, 5}, {1, 3, 4, 7}, {1, 3, 5, 6}, {2, 3, 4, 6}, {2, 3, 5, 7} };
    static const double cf[5] = { -0.36, 0.19, -0.71, 0.22, 0.49 };
    zero(HL); zero(HR);
    for (int t = 0; t < 5; t++) {
        const double c = cf[t] * (double)ONE / 4.0; const q cq40 = (q)(c < 0 ? c - 0.5 : c + 0.5);
        addterm(pmul(pmul(pmul(ML[idx[t][0]], ML[idx[t][1]]), ML[idx[t][2]]), ML[idx[t][3]]), (int64_t)cq40 << 20, 0);
        if (check) addterm(pmul(pmul(pmul(MR[idx[t][0]], MR[idx[t][1]]), MR[idx[t][2]]), MR[idx[t][3]]), (int64_t)cq40 << 20, 1);
    }
    NM = 7;
    printf("WHSYK6 model 1: 5 terms on 7 Majoranas a side (coefficients -0.36, 0.19, -0.71, 0.22, 0.49 declared in Q.40 as v1 did), V normalized with N = 10\n");
}
/* the two series from one set of powers of s H, s = max(dt, beta/4): W = e^{-i dt H}, E = e^{-(beta/4) H}; order 1: P_{k-1} H (the left side), 0: H P_{k-1} */
static void series(c60 Hm[MAXD][MAXD], int64_t bq60, int64_t dt60, int order, c60 W[MAXD][MAXD], c60 E[MAXD][MAXD], c60 P[MAXD][MAXD], c60 u[MAXD][MAXD],
                   Acc *SWre, Acc *SWim, Acc *SEre, Acc *SEim)
{
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
        (*SWre)[i][j].hi = (*SWre)[i][j].lo = 0; (*SWim)[i][j].hi = (*SWim)[i][j].lo = 0; (*SEre)[i][j].hi = (*SEre)[i][j].lo = 0; (*SEim)[i][j].hi = (*SEim)[i][j].lo = 0;
        P[i][j].re = P[i][j].im = 0;
    }
    for (int i = 0; i < D; i++) { P[i][i].re = ONE60; (*SWre)[i][i].hi = ONE60 >> 4; (*SWre)[i][i].lo = (uint64_t)ONE60 << 60; (*SEre)[i][i] = (*SWre)[i][i]; }
    const int64_t s60 = dt60 > bq60 ? dt60 : bq60, wr60 = div60(dt60, s60), br60 = div60(bq60, s60);
    int64_t bk = ONE60, wk = ONE60;
    int small = 0, k;
    for (k = 1; k <= 80; k++) {
        if (order) mmul(P, 0, Hm, u); else mmul(Hm, 0, P, u);
        bk = -mul60(bk, br60); wk = mul60(wk, wr60);
        int64_t mx = 0;
        for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
            const c60 p = { rdiv64(mul60(u[i][j].re, s60), k), rdiv64(mul60(u[i][j].im, s60), k) }; P[i][j] = p;
            const int64_t pr = mul60(p.re, wk), pi = mul60(p.im, wk);
            int64_t wr, wi; switch (k & 3) { case 0: wr = pr; wi = pi; break; case 1: wr = pi; wi = -pr; break; case 2: wr = -pr; wi = -pi; break; default: wr = -pi; wi = pr; }
            s128_add(&(*SWre)[i][j], wr >> 4, (uint64_t)wr << 60); s128_add(&(*SWim)[i][j], wi >> 4, (uint64_t)wi << 60);
            s128_addp(&(*SEre)[i][j], bk, p.re); s128_addp(&(*SEim)[i][j], bk, p.im);
            const int64_t a = p.re < 0 ? -p.re : p.re, b = p.im < 0 ? -p.im : p.im; if (a > mx) mx = a; if (b > mx) mx = b;
        }
        if (mx <= 1) { if (++small >= 2) break; } else small = 0;
    }
    if (k > taylor_max) taylor_max = k;
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
        W[i][j].re = s128_round_sh((*SWre)[i][j], F60); W[i][j].im = s128_round_sh((*SWim)[i][j], F60);
        E[i][j].re = s128_round_sh((*SEre)[i][j], F60); E[i][j].im = s128_round_sh((*SEim)[i][j], F60);
    }
}
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


/* ---- v6: dt straight into Q.60 ---- */
static int64_t parse60(const char *s)
{
    int neg = 0; if (*s == '-') { neg = 1; s++; }
    uint64_t ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (uint64_t)(*s++ - '0');
    uint64_t fp = 0, scale = 1; if (*s == '.') { s++; while (*s >= '0' && *s <= '9' && scale < 1000000000000ull) { fp = fp * 10 + (uint64_t)(*s++ - '0'); scale *= 10; } }
    if (ip >= 7) { printf("WHSYK6: dt must be below 7\n"); exit(1); }
    int64_t v = (int64_t)(ip << F60);
    if (fp) { uint64_t h, l, rem; mul128(fp, (uint64_t)1 << F60, &h, &l); uint64_t qq = div128(h, l, scale, &rem); if (2 * rem >= scale) qq++; v += (int64_t)qq; }
    return neg ? -v : v;
}
/* mu/(8 vnorm) in Q.60 from the Q.40 mu: round(mu40 2^20 / (8 vnorm)) -- the same arithmetic as the protocol loop's */
static int64_t mu_arg60(q m, int vnormN)
{
    int neg = m < 0; uint64_t um = neg ? -(uint64_t)m : (uint64_t)m, h, l, rem; mul128(um, (uint64_t)1 << 20, &h, &l);
    uint64_t qq = div128(h, l, (uint64_t)(8 * vnormN), &rem); if (2 * rem >= (uint64_t)(8 * vnormN)) qq++;
    return neg ? -(int64_t)qq : (int64_t)qq;
}
/* (a + ib)(c + id) in Q.60, one rounding per component */
static void cmul60(int64_t a, int64_t b, int64_t c, int64_t d, int64_t *re, int64_t *im)
{
    s128 r = {0, 0}, i = {0, 0};
    s128_addp(&r, a, c); s128_addp(&r, -b, d); s128_addp(&i, a, d); s128_addp(&i, b, c);
    *re = s128_round_sh(r, F60); *im = s128_round_sh(i, F60);
}

/* ---- v6: THE SIZE-WINDING OBSERVABLE (see the header) ---- */
#define WSHIFT 4                                            /* Q.60 -> Q.56 for the transform's int64 sums */
static int64_t fre[MAXD], fim[MAXD];                        /* one diagonal of O' */
static uint64_t wchain = 0xcbf29ce484222325ull;             /* the winding lines' own chain */
static void whash64(int64_t v) { wchain = fnv(wchain, &v, 8); }
static inline int64_t rsh_round(int64_t v, int sh) { s128 a = { v < 0 ? -1 : 0, (uint64_t)v }; return s128_round_sh(a, sh); }
static int wind_size(uint32_t x, uint32_t z)                /* the Majorana size of sigma_{x,z} */
{
    int size = 0, p = 0;
    for (int k = H - 1; k >= 0; k--) {
        const int xb = (x >> k) & 1, zb = (z >> k) & 1; int c;
        if (xb) c = 1; else if (zb) c = p ? 0 : 2; else c = p ? 2 : 0;
        size += c; if (c & 1) p ^= 1;
    }
    return size;
}
static void build_Oprime(c60 A[MAXD][MAXD], c60 out[MAXD][MAXD])   /* O' = A Pi^dag */
{
    for (int c = 0; c < D; c++) for (int L = 0; L < D; L++) {
        const c60 v = A[c][L ^ (D - 1)]; c60 w;
        switch ((4 - EOFL[L]) & 3) { case 0: w = v; break; case 1: w.re = -v.im; w.im = v.re; break; case 2: w.re = -v.re; w.im = -v.im; break; default: w.re = v.im; w.im = -v.re; }
        out[c][L] = w;
    }
}
/* the transform of O': the per-size sums of T^2 (complex) and |T|^2 in 128 bits; even_max = the largest |T| on an even-popcount diagonal (check only) */
static void wind_transform(c60 Op[MAXD][MAXD], int check, s128 Qre[], s128 Qim[], s128 Pr[], int64_t *even_max)
{
    for (int l = 0; l <= N; l++) { Qre[l].hi = 0; Qre[l].lo = 0; Qim[l].hi = 0; Qim[l].lo = 0; Pr[l].hi = 0; Pr[l].lo = 0; }
    *even_max = 0;
    for (uint32_t x = 0; x < (uint32_t)D; x++) {
        const int evenx = (pop32(x) & 1) == 0;
        if (evenx && !check) continue;                          /* exact zeros: A is fermion-odd */
        for (int c = 0; c < D; c++) { const c60 v = Op[c][c ^ (int)x]; fre[c] = rsh_round(v.re, WSHIFT); fim[c] = rsh_round(v.im, WSHIFT); }
        for (int h = 1; h < D; h <<= 1)                          /* the Walsh-Hadamard butterflies */
            for (int i = 0; i < D; i += 2 * h) for (int j = i; j < i + h; j++) {
                const int64_t ar = fre[j], ai = fim[j], br = fre[j + h], bi = fim[j + h];
                fre[j] = ar + br; fim[j] = ai + bi; fre[j + h] = ar - br; fim[j + h] = ai - bi;
            }
        for (uint32_t z = 0; z < (uint32_t)D; z++) {
            int64_t tr = fre[z], ti = fim[z];
            switch (pop32(x & z) & 3) { case 1: { const int64_t t = tr; tr = -ti; ti = t; } break; case 2: tr = -tr; ti = -ti; break; case 3: { const int64_t t = tr; tr = ti; ti = -t; } break; default: break; }
            if (evenx) { const int64_t a = tr < 0 ? -tr : tr, b = ti < 0 ? -ti : ti; if (a > *even_max) *even_max = a; if (b > *even_max) *even_max = b; }   /* check mode: these strings are summed into their (even) sizes too */
            const int l = wind_size(x, z);
            s128_addp(&Qre[l], tr, tr); s128_addp(&Qre[l], -ti, ti);
            s128_addp(&Qim[l], tr, ti); s128_addp(&Qim[l], tr, ti);
            s128_addp(&Pr[l], tr, tr); s128_addp(&Pr[l], ti, ti);
        }
    }
}
/* the winding lines for one operator state A at the operator time tp: the per-size sums, W, C(mu) from q(l) for every mu, and with
   check the state route for C(mu) */
static void wind_report(int op, uint64_t seed, int64_t beta60, int tp, int64_t dt60, const q *mus, int nmu, int vnormN, int check,
                        c60 A[MAXD][MAXD])
{
    static s128 Qre[MAXN + 1], Qim[MAXN + 1], Pr[MAXN + 1];
    static int64_t Pl[MAXN + 1], qre[MAXN + 1], qim[MAXN + 1];
    char S0[48], S1[48], S2[48], S3[48], S4[48], S5[48], S6[48], SB[48], line[2048]; int64_t even_max;
    build_Oprime(A, XS(0));
    wind_transform(XS(0), check, Qre, Qim, Pr, &even_max);
    const int sh = 2 * (F60 - WSHIFT) - F60 + H;              /* Q.112 sums over D strings -> Q.60: 52 + n */
    s128 psum = {0, 0}; int64_t even_sizes = 0, wre = 0, wim = 0;
    for (int l = 0; l <= N; l++) {
        const int64_t a = s128_round_sh(Qre[l], sh), b = s128_round_sh(Qim[l], sh);
        Pl[l] = s128_round_sh(Pr[l], sh); s128_adds(&psum, Pr[l]);
        if (H & 1) { qre[l] = b; qim[l] = -a; } else { qre[l] = a; qim[l] = b; }   /* omega^-2 = -i for n odd */
        if ((l & 1) == 0) { if (Pl[l] > even_sizes) even_sizes = Pl[l]; }
        else { wre += qre[l]; wim += qim[l]; }
    }
    const int64_t psum60 = s128_round_sh(psum, sh) - ONE60;
    fmt60(beta60, SB); fmt60((int64_t)tp * dt60, S1);
    int n = sprintf(line, "WHSYK6 N=%d seed=%llu beta=%s t=%s op=%d: WINDING", N, (unsigned long long)seed, SB, S1, op);
    for (int l = 1; l <= N; l += 2) {
        fmt60(Pl[l], S2); fmt60(qre[l], S3); fmt60(qim[l], S4);
        n += sprintf(line + n, " | l=%d P %s q %s %s i", l, S2, S3, S4);
        whash64(Pl[l]); whash64(qre[l]); whash64(qim[l]);
    }
    fmt60(wre, S5); fmt60(wim, S6);
    printf("%s || sum P(l) - 1 = %lld units, even sizes (a check only with windcheck=1) at most %lld units, W = %s %s i%s\n", line, (long long)psum60, (long long)even_sizes, S5, S6,
           check ? (even_max ? " -- EVEN-x DIAGONALS NOT ZERO" : ", even-x diagonals exactly 0") : "");
    whash64(psum60); whash64(wre); whash64(wim);
    if (check) {                                               /* B = e^{iHt} g_i^R e^{-iHt} T into XS(1) */
        copy(TT, ZS(0)); evolve(ZS(0), tp, 0); papply(&MR[op], ZS(0), XS(1)); evolve(XS(1), tp, 1);
    }
    for (int mi = 0; mi < nmu; mi++) {
        const int64_t arg60 = mu_arg60(mus[mi], vnormN);
        int64_t c1, s1; cossin60(arg60, &c1, &s1);             /* e^{i arg}, a small angle */
        int64_t e2r, e2i, e4r, e4i; cmul60(c1, s1, c1, s1, &e2r, &e2i); cmul60(e2r, e2i, e2r, e2i, &e4r, &e4i);
        int64_t pr = e2r, pi = e2i;                              /* e^{i 2 l arg}, l = 1, then stepping by e^{i 4 arg} */
        s128 sre = {0, 0}, sim = {0, 0};
        for (int l = 1; l <= N; l += 2) {                     /* S = SUM_l q(l) e^{i mu l/(4 vnorm)} */
            s128_addp(&sre, qre[l], pr); s128_addp(&sre, -qim[l], pi);
            s128_addp(&sim, qre[l], pi); s128_addp(&sim, qim[l], pr);
            cmul60(pr, pi, e4r, e4i, &pr, &pi);
        }
        const int64_t Sre = s128_round_sh(sre, F60), Sim = s128_round_sh(sim, F60);
        int64_t fr = ONE60, fi = 0;                              /* e^{i NM arg}; its conjugate is the prefactor's phase */
        for (int k = 0; k < NM; k++) cmul60(fr, fi, c1, s1, &fr, &fi);
        s128 pre = {0, 0}, pim = {0, 0};                         /* (Sre + i Sim)(fr - i fi) */
        s128_addp(&pre, Sre, fr); s128_addp(&pre, Sim, fi);
        s128_addp(&pim, Sim, fr); s128_addp(&pim, -Sre, fi);
        const int64_t Ere = s128_round_sh(pre, F60), Eim = s128_round_sh(pim, F60);
        const int64_t Cre = rdiv64(Eim, 2), Cim = -rdiv64(Ere, 2);   /* times -i/2 */
        s128 m2 = {0, 0}; s128_addp(&m2, Cre, Cre); s128_addp(&m2, Cim, Cim); const int64_t twoC = 2 * sqrt120(m2);
        fmt(mus[mi], S0); fmt60(Cre, S2); fmt60(Cim, S3); fmt60(twoC, S4);
        whash64(Cre); whash64(Cim);
        if (check) {
            copy(A, ZS(3)); apply_U(ZS(3), c1, s1, 1, WX);
            s128 kre = {0, 0}, kim = {0, 0}; inner_acc(XS(1), ZS(3), 0, &kre, &kim);
            const int64_t Dre = rdiv64(s128_round_sh(kre, F60), 2), Dim = rdiv64(s128_round_sh(kim, F60), 2);
            int64_t d1 = Dre - Cre, d2 = Dim - Cim; if (d1 < 0) d1 = -d1; if (d2 < 0) d2 = -d2; if (d2 > d1) d1 = d2;
            fmt60(Dre, S5); fmt60(Dim, S6);
            printf("WHSYK6 N=%d seed=%llu beta=%s t=%s op=%d mu=%s: C(mu) from q(l) %s %s i, 2|C| %s | from the states %s %s i | differ by %lld units of 2^-60\n",
                   N, (unsigned long long)seed, SB, S1, op, S0, S2, S3, S4, S5, S6, (long long)d1);
        } else
            printf("WHSYK6 N=%d seed=%llu beta=%s t=%s op=%d mu=%s: C(mu) from q(l) %s %s i, 2|C| %s\n", N, (unsigned long long)seed, SB, S1, op, S0, S2, S3, S4);
    }
}
int main(int argc, char **argv)
{
    N = 8; int vnormN = 0, model = 0, check = 0, wind = 1, windcheck = 0; q beta = 4 * ONE, dt = ONE; int64_t dt60_in = ONE60; q mus[16], ts[32], tps[32]; int nmu = 0, nt = 0, ntp = 0; uint64_t seeds[16]; int nseed = 0;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "model=", 6)) model = !strcmp(argv[i] + 6, "m1");
        else if (!strncmp(argv[i], "dt=", 3)) { dt = parseq(argv[i] + 3); dt60_in = parse60(argv[i] + 3); }
        else if (!strncmp(argv[i], "wind=", 5)) wind = atoi(argv[i] + 5);
        else if (!strncmp(argv[i], "windcheck=", 10)) windcheck = atoi(argv[i] + 10);
        else if (!strncmp(argv[i], "check=", 6)) check = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "beta=", 5)) beta = parseq(argv[i] + 5);
        else if (!strncmp(argv[i], "mu=", 3)) nmu = parselist(argv[i] + 3, mus, 16);
        else if (!strncmp(argv[i], "t=", 2)) nt = parselist(argv[i] + 2, ts, 32);
        else if (!strncmp(argv[i], "tp=", 3)) ntp = parselist(argv[i] + 3, tps, 32);
        else if (!strncmp(argv[i], "pair=", 5)) sscanf(argv[i] + 5, "%d,%d", &PA, &PB);
        else if (!strncmp(argv[i], "vnorm=", 6)) vnormN = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "seeds=", 6)) { char b[256]; strncpy(b, argv[i] + 6, 255); b[255] = 0; for (char *t = strtok(b, ","); t && nseed < 16; t = strtok(NULL, ",")) seeds[nseed++] = strtoull(t, 0, 10); }
    }
    if (N < 4 || N > MAXN || (N & 1)) { printf("WHSYK6: N must be even, 4..%d\n", MAXN); return 1; }
    if (!nmu) { mus[0] = -12 * ONE; mus[1] = 12 * ONE; nmu = 2; }
    if (!nt) { ts[0] = ONE; ts[1] = 2 * ONE; ts[2] = 3 * ONE; nt = 3; }
    if (!nseed) { seeds[0] = 1; nseed = 1; }
    if (model) { N = 8; nseed = 1; if (!vnormN) vnormN = 10; }
    if (!vnormN) vnormN = N;
    if (dt <= 0 || dt > 4 * ONE) { printf("WHSYK6: dt must be in (0, 4]\n"); return 1; }
    if (PA < 1 || PA > N || PB < 1 || PB > N || PA == PB) { printf("WHSYK6: pair must be two different Majorana indices in 1..N\n"); return 1; }
    int tmax = 0; for (int i = 0; i < nt; i++) { int t = (int)(ts[i] / ONE); if ((q)t * ONE != ts[i] || t < 0 || t > TMAX) { printf("WHSYK6: t must be an integer in 0..%d\n", TMAX); return 1; } if (t > tmax) tmax = t; }
    for (int i = 0; i < ntp; i++) { int t = (int)(tps[i] / ONE); if ((q)t * ONE != tps[i] || t < 0 || t > TMAX) { printf("WHSYK6: tp must be an integer in 0..%d\n", TMAX); return 1; } if (t > tmax) tmax = t; }
    H = N / 2; D = 1 << H;
    const int64_t beta60 = (int64_t)beta << 20, bq60 = rdiv64(beta60, 4), dt60 = dt60_in;   /* v6: dt in Q.60 from its text */
    LN2_60 = artanh2_60(div60(ONE60, 3 * ONE60));
    char A[48], B[48], C[48], Dd[48], E[48], F[48], G[48], S1[48], S2[48], S3[48], S4[48], S5[48], S6[48], S7[48], S8[48];
    fmt60(beta60, A); fmt60(LN2_60, B); fmt60(dt60, C);
    printf("WHSYK6: v6 = v5 + the size-winding observable (wind=%d, windcheck=%d; dt parsed in Q.60); Q.%d end to end (int64 with 128-bit products; nothing rounded to Q.40); %s; N = %d Majoranas a side, FACTORIZED (D = %d a side; a pool of %d matrices of %d KB); beta = %s; V normalized with N = %d; the Dirac modes chi = (psi^%d + i psi^%d)/sqrt2; the time unit dt = %s (W(dt) by its own series, W(t) = W(t-1) W(1)) and e^{-beta H/4} from one set of powers per seed; the right-hand operators read from the left by the exact re-indexing; the EXACT rho_TR from its six elements (GJ (2.15)-(2.20)); K = <{A1, B1}>\n",
           wind, windcheck, F60, model ? "MODEL 1 of the 2022 experiment (Kobrin-Schuster-Yao eq. 1: 7 Majoranas, 5 terms)" : "the FULL q = 4 SYK, <J^2> = 6/N^3, J = 1", N, D, NPOOL, (int)(sizeof(c60) * D * D / 1024), A, vnormN, PA, PB, C);
    printf("WHSYK6 selftest: ln 2 = %s (0.693147180560)\n", B);
    build_majoranas();
    build_I(PSII);
    int64_t worst = 0;
    for (int j = 1; j <= N; j++) {
        papply(&ML[j], PSII, ZS(0)); papply(&MR[j], PSII, ZS(1));
        for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { int64_t re = ZS(0)[L][R].re - ZS(1)[L][R].im, im = ZS(0)[L][R].im + ZS(1)[L][R].re; if (re < 0) re = -re; if (im < 0) im = -im; if (re > worst) worst = re; if (im > worst) worst = im; }
    }
    printf("WHSYK6 selftest 1: max |(psi_l^j + i psi_r^j)|I>| = %lld units of 2^-60 (must be 0)\n", (long long)worst);
    for (int si = 0; si < nseed; si++) {
        if (model) build_model1(check); else build_syk(seeds[si], check);
        /* self-test 2: H_L Psi_I = Psi_I H_R^T exactly, with H_R read from H_L by the re-indexing (mode 2) */
        mmul(HL, 0, PSII, ZS(0)); transpose(PSII, WX); mmul(HL, 2, WX, ZS(1)); transpose(ZS(1), ZS(2));
        printf("WHSYK6 seed %llu selftest 2: max |H_L Psi_I - Psi_I H_R^T| = %lld units of 2^-60 (must be 0; H_R read from H_L by the re-indexing)\n", (unsigned long long)seeds[si], (long long)maxabs_diff(ZS(0), ZS(2)));
        hash64(maxabs_diff(ZS(0), ZS(2)));
        if (check) { /* self-test 2b: the re-indexed H_L against the H_R built from the right-hand Pauli strings */
            for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) ZS(0)[i][j] = entry(HL, 2, i, j);
            printf("WHSYK6 seed %llu selftest 2b: max |H_R (built) - the re-indexed H_L| = %lld units (must be 0)\n", (unsigned long long)seeds[si], (long long)maxabs_diff(ZS(0), HR)); }
        /* the series on the left (the accumulators in the Z slots; P and u in WG and WX); then W(t) = W(t-1) W(1) */
        series(HL, bq60, dt60, 1, WL(1), EL, WG, WX, SWRE, SWIM, SERE, SEIM);
        if (check) { /* self-test 5: the right-hand series (opposite order) against the re-indexed left: must be 0 units */
            series(HR, bq60, dt60, 0, XS(0), XS(1), WG, WX, SWRE, SWIM, SERE, SEIM);
            int64_t w5 = 0, e5 = 0;
            for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) {
                const c60 w = entry(WL(1), 2, i, j), e = entry(EL, 2, i, j);
                int64_t d1 = XS(0)[i][j].re - w.re, d2 = XS(0)[i][j].im - w.im, d3 = XS(1)[i][j].re - e.re, d4 = XS(1)[i][j].im - e.im;
                if (d1 < 0) d1 = -d1; if (d2 < 0) d2 = -d2; if (d3 < 0) d3 = -d3; if (d4 < 0) d4 = -d4;
                if (d1 > w5) w5 = d1; if (d2 > w5) w5 = d2; if (d3 > e5) e5 = d3; if (d4 > e5) e5 = d4;
            }
            printf("WHSYK6 seed %llu selftest 5: |the right-side series - the re-indexed left| = %lld units (W) and %lld units (E) of 2^-60 (must be 0)\n", (unsigned long long)seeds[si], (long long)w5, (long long)e5);
        }
        for (int t = 2; t <= tmax; t++) mmul(WL(t - 1), 0, WL(1), WL(t));
        /* self-test 3: the unitarity defect of W_L(1) */
        mmul(WL(1), 1, WL(1), ZS(0)); for (int i = 0; i < D; i++) ZS(0)[i][i].re -= ONE60;
        zero(ZS(1)); const int64_t udef = maxabs_diff(ZS(0), ZS(1)); if (udef > unit_defect_max) unit_defect_max = udef;
        /* self-test 4: E_L Psi_I = Psi_I E_R^T */
        mmul(EL, 0, PSII, ZS(0)); transpose(PSII, WX); mmul(EL, 2, WX, ZS(1)); transpose(ZS(1), ZS(2));
        const int64_t edef = maxabs_diff(ZS(0), ZS(2));
        printf("WHSYK6 seed %llu selftest 3-4: |W_L(1)^dag W_L(1) - 1| <= %lld units of 2^-60; |e^{-beta H_L/4} Psi_I - Psi_I (e^{-beta H_R/4})^T| <= %lld units of 2^-60 (both must be a few units); Taylor terms at most %d\n",
               (unsigned long long)seeds[si], (long long)udef, (long long)edef, taylor_max);
        hash64(udef); hash64(edef);
        copy(PSII, TT); sandwich(EL, 0, TT);
        { const s128 n2 = norm2_128(TT); fmt60(s128_round_sh(n2, F60), B); }
        normalize(TT);
        fmt60(sqrt120(norm2_128(TT)), A);
        printf("WHSYK6 seed %llu: |e^{-beta H/4} I|^2 / |I|^2 = %s (Z_beta / Z_inf); the thermofield double normalized to %s\n", (unsigned long long)seeds[si], B, A);
        for (int mi = 0; mi < nmu; mi++) {
            int64_t arg60; { const q m = mus[mi]; int neg = m < 0; uint64_t um = neg ? -(uint64_t)m : (uint64_t)m, h, l, rem; mul128(um, (uint64_t)1 << 20, &h, &l); uint64_t qq = div128(h, l, (uint64_t)(8 * vnormN), &rem); if (2 * rem >= (uint64_t)(8 * vnormN)) qq++; arg60 = neg ? -(int64_t)qq : (int64_t)qq; }
            int64_t c, s; cossin60(arg60, &c, &s);
            const int64_t csdef = mul60(c, c) + mul60(s, s) - ONE60;
            for (int ti = 0; ti < nt; ti++) {
                const int t = (int)(ts[ti] / ONE), tp = ntp ? (int)(tps[ti < ntp ? ti : ntp - 1] / ONE) : t;
                copy(TT, WG); evolve(WG, tp, 1);
                for (int ab = 0; ab < 4; ab++) { apply2M(ab >> 1, ab & 1, WG, ZS(ab)); evolve(ZS(ab), tp, 0); apply_U(ZS(ab), c, s, 1, WX); evolve(ZS(ab), t, 0); }
                s128 ar[6], ai[6]; memset(ar, 0, sizeof ar); memset(ai, 0, sizeof ai);
                for (int b = 0; b < 2; b++) {
                    c60 (*Z0)[MAXD] = ZS(b), (*Z1)[MAXD] = ZS(2 + b);
                    papply(&MR[PB], Z0, XS(0)); papply(&MR[PA], XS(0), XS(1));
                    inner_acc(Z0, Z0, 0, &ar[0], &ai[0]); inner_acc(Z0, XS(1), 3, &ar[0], &ai[0]);
                    inner_acc(Z0, Z0, 0, &ar[1], &ai[1]); inner_acc(Z0, XS(1), 1, &ar[1], &ai[1]);
                    papply(&MR[PA], Z0, XS(0)); papply(&MR[PB], Z0, XS(1));
                    inner_acc(Z1, XS(0), 0, &ar[4], &ai[4]); inner_acc(Z1, XS(1), 3, &ar[4], &ai[4]);
                    inner_acc(Z1, XS(0), 0, &ar[5], &ai[5]); inner_acc(Z1, XS(1), 1, &ar[5], &ai[5]);
                    papply(&MR[PB], Z1, XS(0)); papply(&MR[PA], XS(0), XS(1));
                    inner_acc(Z1, Z1, 0, &ar[2], &ai[2]); inner_acc(Z1, XS(1), 3, &ar[2], &ai[2]);
                    inner_acc(Z1, Z1, 0, &ar[3], &ai[3]); inner_acc(Z1, XS(1), 1, &ar[3], &ai[3]);
                }
                int64_t rr[6], ri[6]; for (int e = 0; e < 6; e++) { rr[e] = s128_round_sh(ar[e], F60 + 4); ri[e] = s128_round_sh(ai[e], F60 + 4); }
                const int64_t r00 = rr[0], r11 = rr[1], r22 = rr[2], r33 = rr[3];
                const int64_t epr1 = r00 + r11 - (ONE60 >> 1), epr2 = r22 + r33 - (ONE60 >> 1);
                int64_t imdiag = 0; for (int e = 0; e < 4; e++) { int64_t v = ri[e] < 0 ? -ri[e] : ri[e]; if (v > imdiag) imdiag = v; }
                int clamp = 0; int64_t SRT = 0, SR = 0, ST = 0;
                { const int64_t m = rdiv64(r00 + r33, 2), d = rdiv64(r00 - r33, 2); s128 rad = {0, 0}; s128_addp(&rad, d, d); s128_addp(&rad, rr[4], rr[4]); s128_addp(&rad, ri[4], ri[4]);
                  const int64_t disc = sqrt120(rad); SRT += ent60(m + disc, &clamp) + ent60(m - disc, &clamp); }
                { const int64_t m = rdiv64(r11 + r22, 2), d = rdiv64(r11 - r22, 2); s128 rad = {0, 0}; s128_addp(&rad, d, d); s128_addp(&rad, rr[5], rr[5]); s128_addp(&rad, ri[5], ri[5]);
                  const int64_t disc = sqrt120(rad); SRT += ent60(m + disc, &clamp) + ent60(m - disc, &clamp); }
                ST = ent60(r00 + r22, &clamp) + ent60(r11 + r33, &clamp);
                SR = ent60(r00 + r11, &clamp) + ent60(r22 + r33, &clamp);
                const int64_t I = SR + ST - SRT, frac = div60(I, 2 * LN2_60);
                /* K = Re <M1(-t') T | U^dag R1(t) U T>, the states in Z_01, Z_10, Z_11 (free now) */
                papply(&ML[PA], WG, ZS(1)); evolve(ZS(1), tp, 0);
                copy(TT, ZS(2)); apply_U(ZS(2), c, s, 1, WX); evolve(ZS(2), t, 0); papply(&MR[PA], ZS(2), ZS(3)); evolve(ZS(3), t, 1); apply_U(ZS(3), c, s, -1, WX);
                s128 kre = {0, 0}, kim = {0, 0}; inner_acc(ZS(1), ZS(3), 0, &kre, &kim); const int64_t K = s128_round_sh(kre, F60);
                fmt60((int64_t)t * dt60, A); fmt(mus[mi], B); fmt60(r00, C); fmt60(r11, Dd); fmt60(r22, E); fmt60(r33, F); fmt60(rr[4], G); fmt60(ri[4], S1); fmt60(rr[5], S2); fmt60(ri[5], S3);
                fmt60(K, S7); fmt60(frac, S8); fmt60((int64_t)tp * dt60, S4); fmt60(beta60, S5); fmt60(epr1, S6);
                printf("WHSYK6 N=%d seed=%llu beta=%s mu=%s t=%s tp=%s: EXACT rho00 %s rho11 %s rho22 %s rho33 %s | rho03 %s %s i | rho12 %s %s i | checks (R's marginal) rho00+rho11-1/2 %s rho22+rho33-1/2 %lld units, max|Im diag| %lld units, c^2+s^2-1 %lld units of 2^-60 | K %s | I_RT/(2 ln 2) EXACT %s%s\n",
                       N, (unsigned long long)seeds[si], S5, B, A, S4, C, Dd, E, F, G, S1, S2, S3, S6, (long long)epr2, (long long)imdiag, (long long)csdef, S7, S8, clamp ? " (an eigenvalue clamped at 0)" : "");
                hash64(r00); hash64(r11); hash64(r22); hash64(r33); hash64(rr[4]); hash64(ri[4]); hash64(rr[5]); hash64(ri[5]); hash64(K); hash64(frac);
                if (wind && mi == 0) {                               /* v6: the winding at the operator time tp (K's Z_01 = g_PA(-tp) T is A for PA) */
                    const int nops = wind >= 2 ? NM : 2;
                    for (int oi = 0; oi < nops; oi++) {
                        const int op = wind >= 2 ? oi + 1 : (oi == 0 ? PA : PB);
                        c60 (*A)[MAXD] = ZS(1);
                        if (op != PA) { papply(&ML[op], WG, ZS(2)); evolve(ZS(2), tp, 0); A = ZS(2); }
                        wind_report(op, seeds[si], beta60, tp, dt60, mus, nmu, vnormN, windcheck, A);
                    }
                }
            }
        }
    }
    printf("WHSYK6: unitarity defect of W(1) over the seeds at most %lld units of 2^-60; Taylor terms at most %d\n", (long long)unit_defect_max, taylor_max);
    printf("WHSYK6: chain %016llx\n", (unsigned long long)chain);
    if (wind) printf("WHSYK6: winding chain %016llx\n", (unsigned long long)wchain);
    printf("WHSYK6: done\n");
    return 0;
}
