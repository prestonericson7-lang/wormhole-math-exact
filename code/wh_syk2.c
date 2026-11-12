/* wh_syk2.c (11-11-2026, 9:10 PM) -- WORMHOLE MATH 2, v2: the SYK traversable-wormhole teleportation protocol (Gao-Jafferis, arXiv:1911.07416,
   sections 2.1-2.3) exact at finite N, with the FACTORIZED evolution, so N = 14-16 Majoranas a side take minutes on the pair instead of weeks.

   THE IDEA. The two SYK systems l and r act on disjoint qubits once the left Majoranas live on qubits 0..h-1 and the right on h..N-1
   (h = N/2, standard Jordan-Wigner on each half; a right Majorana carries the left half's parity string, which cancels in every 4-Majorana
   term of H_r). Then H = H_L (x) 1 + 1 (x) H_R, the state is a D x D matrix Psi[L][R] with D = 2^h, and
       e^{-iHt} Psi = U_L(t) Psi U_R(t)^T,   U = e^{-iHt} on one side, a D x D matrix exponential,
   two matrix products per evolution (2 D^3 multiply-adds: 4e6 at N = 14) instead of a term-by-term Taylor on 2^N amplitudes. The physics is
   representation-independent; v1 (wh_syk.c, the paper's interleaved layout) and v2 must give the same numbers to rounding, and do (printed).

   THE ENTANGLED STATE |I> in this layout (found numerically at N = 4..10 and proven here exactly by the self-test (psi_l^j + i psi_r^j)|I> = 0):
   it lives on the basis states whose right half is the bitwise complement of the left half, with amplitude (unnormalized) i^{e(L)},
       e(L) = SUM over set bits k of L of (k even ? 3 : 1) + 2 x C(popcount L, 2)   (mod 4),
   i.e. a factor -i for each set even site, +i for each set odd site, and -1 for every pair of set bits; the state with L = 0 has phase 1.
   |I> is kept unnormalized (amplitudes exactly 0, +-1, +-i in Q.40); the thermofield double is normalized after e^{-beta H/4}, as before.

   EXACTNESS. H_L and H_R are built as D x D matrices whose entries are exact int64 sums of the terms' coefficients and phases (no rounding
   at all). A matrix exponential e^{sH} (s = -i dt, +i dt, or -beta/4) is the Taylor series with each term a D x D product accumulated in
   signed 128-bit and rounded once per entry, the series summed in 128-bit and rounded once; U(t) for integer t is U(1)^t by products. Every
   matrix product rounds once per entry. Majoranas and the coupling U = PROD_j (cos + i sin W_j) are permutations with phases on the L or R
   index (exact). Inner products in 128-bit, rounded once. The couplings J_{ijkl} ~ N(0, 6/N^3) from xorshift64* with the declared seed in the
   same order as v1, so v1 and v2 share them bit for bit.

   THE PROTOCOL AND THE OUTPUT: exactly v1's (rho11, rho14, K, I_RT/(2 ln 2) from (2.63)-(2.66) and (2.13)); the same command line:
       whsyk2 N=14 model=syk seeds=1,2,3 beta=4 mu=-12,12,0 t=1,2,3,4,5,6 [tp=...] [pair=1,2] [vnorm=N]
       whsyk2 N=8 model=m1 beta=4 mu=-12,12 t=2.8
   Self-tests printed: |I>'s annihilation, exactly; H_L Psi_I = Psi_I H_R^T exactly (the built matrices against the built state); the norm of
   every real-time evolution; mu = 0 -> I_RT = 0. */
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

static int N, H, D, NM;
static pstr ML[MAXN + 1], MR[MAXN + 1];     /* the bare Majoranas (sqrt2 psi) as Pauli strings over all N qubits */
static int64_t HLre[MAXD][MAXD], HLim[MAXD][MAXD], HRre[MAXD][MAXD], HRim[MAXD][MAXD];   /* exact integer entries (Q.40 units) */
static cq UL[2][MAXD][MAXD], UR[2][MAXD][MAXD];    /* [0] = e^{-iH}, [1] = e^{+iH} for dt = 1, each side */
static cq EL[MAXD][MAXD], ER[MAXD][MAXD];          /* e^{-beta H/4} each side */
static cq PsiI[MAXD][MAXD], T[MAXD][MAXD];
static cq WA[MAXD][MAXD], WB[MAXD][MAXD], WC[MAXD][MAXD], WD[MAXD][MAXD], WE[MAXD][MAXD], WF[MAXD][MAXD], WG[MAXD][MAXD], WH[MAXD][MAXD], WI[MAXD][MAXD], WJ[MAXD][MAXD];
static s128 AR[MAXD][MAXD], AI[MAXD][MAXD];
static uint64_t chain = 0xcbf29ce484222325ull;
static int taylor_max = 0; static q drift_max = 0;

/* ---- 128-bit ---- */
static inline void s128_add(s128 *a, int64_t hi, uint64_t lo) { uint64_t l = a->lo + lo; a->hi += hi + (l < lo); a->lo = l; }
static inline void mulprod(q a, q b, int64_t *hi, uint64_t *lo)
{
    int neg = (a < 0) ^ (b < 0); uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b, h, l;
    mul128(ua, ub, &h, &l); if (neg) { l = ~l + 1; h = ~h + (l == 0); } *hi = (int64_t)h; *lo = l;
}
static inline void s128_addq(s128 *a, q v) { s128_add(a, v >> (64 - FB), (uint64_t)v << FB); }
static inline q s128_round(s128 a)
{
    int neg = a.hi < 0; if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t r = (a.lo >> FB) | ((uint64_t)a.hi << (64 - FB)); if ((a.lo & MASK) >= ((uint64_t)1 << (FB - 1))) r++;
    return neg ? -(q)r : (q)r;
}
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
    if (side && (P.x & (uint32_t)(D - 1))) { printf("WHSYK2: a right term touches the left half -- stop\n"); exit(1); }
    if (side && (P.z & (uint32_t)(D - 1))) { printf("WHSYK2: a right term carries a left parity string -- stop\n"); exit(1); }
    if (!side && ((P.x >> H) || (P.z >> H))) { printf("WHSYK2: a left term touches the right half -- stop\n"); exit(1); }
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
    NM = N; { char fb[32]; fmt(first, fb); printf("WHSYK2 seed %llu: %d terms, first coefficient %s (v1 prints the same two numbers when the couplings are the same)\n", (unsigned long long)seed, nt, fb); }
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
    s128 re = {0, 0}, im = {0, 0}; int64_t h; uint64_t l;
    for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) {
        mulprod(A[L][R].re, B[L][R].re, &h, &l); s128_add(&re, h, l); mulprod(A[L][R].im, B[L][R].im, &h, &l); s128_add(&re, h, l);
        mulprod(A[L][R].re, B[L][R].im, &h, &l); s128_add(&im, h, l); mulprod(-A[L][R].im, B[L][R].re, &h, &l); s128_add(&im, h, l);
    }
    cq r = { s128_round(re), s128_round(im) }; return r;
}
static void normalize(cq Psi[MAXD][MAXD]) { const q n = sqrtq(norm2(Psi)); for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { Psi[L][R].re = divq(Psi[L][R].re, n); Psi[L][R].im = divq(Psi[L][R].im, n); } }
static void copy(cq A[MAXD][MAXD], cq B[MAXD][MAXD]) { memcpy(B, A, sizeof(cq) * MAXD * MAXD); }
static q maxabs(cq A[MAXD][MAXD]) { q m = 0; for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { q a = A[L][R].re < 0 ? -A[L][R].re : A[L][R].re, b = A[L][R].im < 0 ? -A[L][R].im : A[L][R].im; if (a > m) m = a; if (b > m) m = b; } return m; }
static void hashq(q v) { chain = fnv(chain, &v, 8); }
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
    N = 8; int model = 0, pa = 1, pb = 2, vnormN = 0; q beta = 4 * ONE; q mus[16], ts[32], tps[32]; int nmu = 0, nt = 0, ntp = 0; uint64_t seeds[16]; int nseed = 0;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "N=", 2)) N = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "model=", 6)) model = !strcmp(argv[i] + 6, "m1");
        else if (!strncmp(argv[i], "beta=", 5)) beta = parseq(argv[i] + 5);
        else if (!strncmp(argv[i], "mu=", 3)) nmu = parselist(argv[i] + 3, mus, 16);
        else if (!strncmp(argv[i], "t=", 2)) nt = parselist(argv[i] + 2, ts, 32);
        else if (!strncmp(argv[i], "tp=", 3)) ntp = parselist(argv[i] + 3, tps, 32);
        else if (!strncmp(argv[i], "pair=", 5)) sscanf(argv[i] + 5, "%d,%d", &pa, &pb);
        else if (!strncmp(argv[i], "vnorm=", 6)) vnormN = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "seeds=", 6)) { char b[256]; strncpy(b, argv[i] + 6, 255); b[255] = 0; for (char *t = strtok(b, ","); t && nseed < 16; t = strtok(NULL, ",")) seeds[nseed++] = strtoull(t, 0, 10); }
    }
    if (N < 4 || N > MAXN || (N & 1)) { printf("WHSYK2: N must be even, 4..%d\n", MAXN); return 1; }
    if (!nmu) { mus[0] = -12 * ONE; mus[1] = 12 * ONE; nmu = 2; }
    if (!nt) { ts[0] = ONE; ts[1] = 2 * ONE; ts[2] = 3 * ONE; nt = 3; }
    if (!nseed) { seeds[0] = 1; nseed = 1; }
    if (model) { N = 8; nseed = 1; if (!vnormN) vnormN = 10; }
    if (!vnormN) vnormN = N;
    H = N / 2; D = 1 << H;
    char A[32], B[32], C[32], Dd[32], E[32], F[32], G[32];
    fmt(beta, A);
    printf("WHSYK2: Q.%d (whq.h); %s; N = %d Majoranas a side, FACTORIZED (left on qubits 0..%d, right on %d..%d; D = %d a side); beta = %s; V normalized with N = %d; the Dirac fermion read = (psi^%d + i psi^%d)/sqrt2; the evolution U(1)^t by Taylor matrix exponentials\n",
           FB, model ? "MODEL 1 of the 2022 experiment (Kobrin-Schuster-Yao eq. 1: 7 Majoranas, 5 terms)" : "the FULL q = 4 SYK, <J^2> = 6/N^3, J = 1", N, H - 1, H, N - 1, D, A, vnormN, pa, pb);
    build_majoranas();
    build_I(PsiI);
    /* self-test 1: (bare_l^j + i bare_r^j) |I> = 0 exactly for every j (unnormalized |I>) */
    q worst = 0;
    for (int j = 1; j <= N; j++) {
        papply(&ML[j], PsiI, WA); papply(&MR[j], PsiI, WB);
        for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { q re = WA[L][R].re - WB[L][R].im, im = WA[L][R].im + WB[L][R].re; if (re < 0) re = -re; if (im < 0) im = -im; if (re > worst) worst = re; if (im > worst) worst = im; }
    }
    fmt(worst, A); printf("WHSYK2 selftest 1: max |(psi_l^j + i psi_r^j)|I>| = %s (must be 0.000000000000)\n", A);
    for (int si = 0; si < nseed; si++) {
        if (model) build_model1(); else build_syk(seeds[si]);
        /* self-test 2: H_L Psi_I = Psi_I H_R^T exactly (the matrices are exact integers, |I> is exact: no rounding anywhere) */
        {
            hmat(HLre, HLim, WC); hmat(HRre, HRim, WD);
            mmul(WC, PsiI, WA); mmulT(PsiI, WD, WB);
            q w2 = 0; for (int L = 0; L < D; L++) for (int R = 0; R < D; R++) { q re = WA[L][R].re - WB[L][R].re, im = WA[L][R].im - WB[L][R].im; if (re < 0) re = -re; if (im < 0) im = -im; if (re > w2) w2 = re; if (im > w2) w2 = im; }
            fmt(w2, A); fmt(maxabs(WA), B);
            printf("WHSYK2 seed %llu selftest 2: max |H_L Psi_I - Psi_I H_R^T| = %s (must be 0); |H_L Psi_I| max entry %s\n", (unsigned long long)seeds[si], A, B);
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
        printf("WHSYK2 seed %llu: |e^{-beta H/4} I|^2 / |I|^2 = %s (Z_beta / Z_inf); normalized to %s; Taylor terms at most %d\n", (unsigned long long)seeds[si], A, B, taylor_max);
        hashq(divq(nT, nI));
        for (int mi = 0; mi < nmu; mi++) {
            const q arg = divq(mus[mi], (q)(2 * 4 * vnormN) * ONE), c = cosq(arg), s = sinq(arg);
            for (int ti = 0; ti < nt; ti++) {
                const q tq = ts[ti], tpq = ntp ? tps[ti < ntp ? ti : ntp - 1] : tq;
                const int t = (int)(tq / ONE), tp = (int)(tpq / ONE);
                if ((q)t * ONE != tq || (q)tp * ONE != tpq) { printf("WHSYK2: v2 takes integer t and t' (U(1)^t); got non-integer -- skipped\n"); continue; }
                drift_max = 0;
                /* a1 = e^{-iHt'} m_l^a e^{iHt'} T; a2 likewise with m_l^b; a21 = A2 a1; a12 = A1 a2 */
                copy(T, WG); evolve(WG, tp, 1, WE);                              /* w = e^{iHt'} T */
                papply(&ML[pa], WG, WA); evolve(WA, tp, 0, WE);                  /* a1 */
                papply(&ML[pb], WG, WB); evolve(WB, tp, 0, WE);                  /* a2 */
                copy(WA, WH); evolve(WH, tp, 1, WE); papply(&ML[pb], WH, WC); evolve(WC, tp, 0, WE);   /* a21 */
                copy(WB, WH); evolve(WH, tp, 1, WE); papply(&ML[pa], WH, WD); evolve(WD, tp, 0, WE);   /* a12 */
                /* the B side on T: up = e^{-iHt} U T; b1 = U^dag e^{iHt} m_r^a up; bb = U^dag e^{iHt} m_r^a m_r^b up */
                copy(T, WG); apply_U(WG, c, s, 1, WE); evolve(WG, t, 0, WE);     /* up */
                papply(&MR[pa], WG, WF); evolve(WF, t, 1, WE); apply_U(WF, c, s, -1, WE);              /* b1 = WF */
                const cq Y1 = inner(WA, WF);
                papply(&MR[pb], WG, WH); papply(&MR[pa], WH, WI); evolve(WI, t, 1, WE); apply_U(WI, c, s, -1, WE);   /* bb = B12 T = WI */
                const cq X1 = inner(WC, WI);
                /* B12 a2 and B2 a1, B2 a12 */
                copy(WB, WG); apply_U(WG, c, s, 1, WE); evolve(WG, t, 0, WE);
                papply(&MR[pb], WG, WH); papply(&MR[pa], WH, WJ); evolve(WJ, t, 1, WE); apply_U(WJ, c, s, -1, WE);   /* B12 a2 */
                const cq X2 = inner(WA, WJ);
                copy(WA, WG); apply_U(WG, c, s, 1, WE); evolve(WG, t, 0, WE);
                papply(&MR[pb], WG, WH); evolve(WH, t, 1, WE); apply_U(WH, c, s, -1, WE);                         /* B2 a1 */
                const cq Y2 = inner(WD, WH);
                copy(WD, WG); apply_U(WG, c, s, 1, WE); evolve(WG, t, 0, WE);
                papply(&MR[pb], WG, WH); evolve(WH, t, 1, WE); apply_U(WH, c, s, -1, WE);                         /* B2 a12 */
                const cq Y3 = inner(WA, WH);
                const q rho11 = (ONE >> 1) - rdiv(X1.re - X2.re, 4);
                const q K = Y1.re;
                const q r14re = rdiv(Y1.re, 2) - rdiv(Y2.re + Y3.re, 4), r14im = -rdiv(Y2.im + Y3.im, 4);
                const q m14 = sqrtq(mulq(r14re, r14re) + mulq(r14im, r14im));
                const q l1 = rdiv(rho11 + m14, 2), l2 = rdiv(rho11 - m14, 2), l3 = rdiv(ONE - rho11, 2);
                q S = 0; int clamp = 0; const q lam[4] = { l1, l2, l3, l3 };
                for (int e = 0; e < 4; e++) { if (lam[e] <= 0) { clamp = 1; continue; } S -= mulq(lam[e], lnq(lam[e])); }
                const q I = 2 * LN2 - S, frac = divq(I, 2 * LN2);
                /* the norm drift: the evolved a1 against T's norm */
                { q n1 = norm2(WA), d = n1 - ONE; if (d < 0) d = -d; if (d > drift_max) drift_max = d; }
                char HB[32], HT[32], HD[32]; fmt(beta, HB); fmt(tpq, HT); fmt(drift_max, HD);
                fmt(tq, A); fmt(mus[mi], B); fmt(rho11, C); fmt(r14re, Dd); fmt(r14im, E); fmt(K, F); fmt(frac, G);
                printf("WHSYK2 N=%d seed=%llu beta=%s mu=%s t=%s tp=%s: rho11 %s | rho14 %s %s i | K %s | I_RT/(2 ln 2) %s%s | norm drift %s | Taylor terms <= %d\n",
                       N, (unsigned long long)seeds[si], HB, B, A, HT, C, Dd, E, F, G, clamp ? " (an eigenvalue clamped at 0)" : "", HD, taylor_max);
                hashq(rho11); hashq(r14re); hashq(r14im); hashq(K); hashq(frac);
            }
        }
    }
    printf("WHSYK2: chain %016llx\n", (unsigned long long)chain);
    printf("WHSYK2: done\n");
    return 0;
}
