/* wh_syk.c (11-11-2026, 7:40 PM) -- WORMHOLE MATH 2: the SYK traversable-wormhole teleportation protocol of Gao and Jafferis
   (arXiv:1911.07416, sections 2.1-2.3) computed EXACTLY at finite N on the two boards, the same program on both, every line to agree
   bit for bit. Nothing here builds a door: it is the number the 2022 "wormhole on a quantum processor" could not compute, the FULL
   q = 4 SYK protocol, and the same protocol on the 2022 experiment's learned five-term model (Kobrin, Schuster and Yao, arXiv:2302.07897,
   eq. 1), so the live dispute gets exact numbers anyone can check.

   WHAT IS COMPUTED. Two SYK systems l and r, N Majoranas each on N qubits (their eq. 2.1-2.2: block j of two qubits carries psi_l^{2j-1},
   psi_l^{2j}, psi_r^{2j-1}, psi_r^{2j} with the string (ZX)^{j-1}); H_l = -SUM_{i<j<k<l} J_{ijkl} psi psi psi psi (i^{q/2} = -1), H_r the
   same couplings (their 2.37 for q = 4); |I> annihilated by psi_l^j + i psi_r^j (2.35): the product over blocks of (|00> - |01> + |10> +
   |11>)/2; |tfd> = e^{-beta H/4} |I> normalized (2.36); V = (i/qN) SUM_j psi_l^j psi_r^j and U = e^{i mu V} (2.50), applied exactly as the
   product over j of cos(mu/2qN) + i sin(mu/2qN) W_j with W_j = 2i psi_l^j psi_r^j (W_j^2 = 1); the protocol inserts the qubit into the Dirac
   fermion chi_l^1 = (psi_l^1 + i psi_l^2)/sqrt2 at time -t' and reads chi_r^1 at +t; the reduced state of the reference R and the read
   qubit T is their (2.13) with (2.63)-(2.66):
       rho_TR = (1/2) [[rho11, 0, 0, rho14], [0, 1-rho11, 0, 0], [0, 0, 1-rho11, 0], [rho14*, 0, 0, rho11]],
       rho11 = (1 - <{A1, [A2, B12]}>)/2,   rho14 = <{A1, B1}>/2 - <{A2, A1 B2 A1}>,
   with A_i = psi_l^i(-t') = e^{-iHt'} psi_l^i e^{iHt'}, B_i = U^dag psi_r^i(t) U, B12 = B1 B2, H = H_l + H_r, all in |tfd>. The mutual
   information I_RT = S(R) + S(T) - S(RT) = 2 ln 2 + SUM lambda ln lambda over the four eigenvalues (rho11 +- |rho14|)/2, (1-rho11)/2 twice;
   K = <{A1, B1}> is printed beside it (their large-N variable; their Table 1: max I_RT/(2 ln 2) = 0.540 at q = 4).

   HOW, EXACTLY. Amplitudes are complex Q.40 (whq.h). Majoranas are applied as bare Pauli strings (sqrt2 psi, so each is an involution) and
   the 1/sqrt2 factors are restored at the end as exact shifts (two Majoranas per correlator term: 1/2; four: 1/4). H v is accumulated in
   signed 128-bit per amplitude over every term and rounded ONCE. e^{sH} v by Taylor in steps of dt = 1 (|H| dt < 1.3 for N <= 14, J = 1):
   term_{k+1} = s H term_k / (k+1), each rounded to nearest, the sum in 128-bit, stopped when a term is below one unit (2^-40). The random
   couplings J_{ijkl} ~ N(0, 6/N^3) come from xorshift64* with a declared seed (twelve 40-bit uniforms summed), so the ensemble is
   reproducible bit for bit. Self-tests printed: (psi_l^j + i psi_r^j)|I> = 0 exactly; {bare_a, bare_b} = 2 delta exactly; (H_l - H_r)|I> = 0
   exactly (the 128-bit sums are exact, so the rounded vectors must be identical); the norm of every real-time evolution (drift printed);
   mu = 0 must give rho11 = 1/2, rho14 = 0, I_RT = 0 up to rounding (the two sides' operators anticommute: an exact zero of the algebra).

   RUN:  whsyk N=8 model=syk seeds=1,2 beta=4 mu=-12,12 t=1,2,3 [tp=<list, default t>] [pair=1,2] [vnorm=N]
         whsyk N=8 model=m1 beta=4 mu=-12,12 t=2.8 [pair=1,2]     (the 2022 model: 7 Majoranas, V normalized with N = 10)
   Every value is printed with fmt (12 decimals) and folded into a chain hash; the last line is "WHSYK: done". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "whq.h"

#define MAXN 14
#define MAXDIM (1 << MAXN)
#define MAXTERMS 1001                      /* C(14, 4) */
#define NVEC 14

typedef struct { q re, im; } cq;
typedef struct { int64_t hi; uint64_t lo; } s128;
typedef struct { uint32_t x, z; int ph; } pstr;      /* i^ph X^x Z^z on N qubits; bit k of a basis index = qubit k */

static int N, DIM, NM;                                /* Majoranas a side (even), 2^N, Majoranas used (7 for model 1) */
static pstr ML[MAXN + 1], MR[MAXN + 1];               /* the bare Majoranas sqrt2 psi_l^j, sqrt2 psi_r^j, j = 1..N */
static pstr HL[MAXTERMS], HR[MAXTERMS]; static q HC[MAXTERMS]; static int NT;   /* the Hamiltonian's terms (bare Paulis) and coefficients */
static cq V[NVEC][MAXDIM];                            /* the working vectors */
static s128 AR[MAXDIM], AI[MAXDIM];                   /* the 128-bit accumulators */
static uint64_t chain = 0xcbf29ce484222325ull;
static int taylor_max = 0; static q drift_max = 0;

/* ---------------- 128-bit signed accumulation ---------------- */
static inline void s128_add(s128 *a, int64_t hi, uint64_t lo) { uint64_t l = a->lo + lo; a->hi += hi + (l < lo); a->lo = l; }
static inline void mulprod(q a, q b, int64_t *hi, uint64_t *lo)        /* the exact signed 128-bit product of two Q.40 values */
{
    int neg = (a < 0) ^ (b < 0); uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b, h, l;
    mul128(ua, ub, &h, &l);
    if (neg) { l = ~l + 1; h = ~h + (l == 0); }
    *hi = (int64_t)h; *lo = l;
}
static inline void s128_addq(s128 *a, q v) { s128_add(a, v >> (64 - FB) /* arithmetic */, (uint64_t)v << FB); }   /* v scaled to 2^80 */
static inline q s128_round(s128 a)                                      /* a / 2^40, rounded to nearest, half away from zero */
{
    int neg = a.hi < 0;
    if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t lo = a.lo, hi = (uint64_t)a.hi;
    uint64_t r = (lo >> FB) | (hi << (64 - FB));
    if ((lo & MASK) >= ((uint64_t)1 << (FB - 1))) r++;
    return neg ? -(q)r : (q)r;
}
static inline void acc_mac(s128 *are, s128 *aim, q c, q vre, q vim, int rot)   /* (are, aim) += i^rot * c * (vre + i vim) */
{
    int64_t h1, h2; uint64_t l1, l2;
    mulprod(c, vre, &h1, &l1); mulprod(c, vim, &h2, &l2);               /* p = c v = (h1:l1) + i (h2:l2) */
    switch (rot & 3) {
        case 0: s128_add(are, h1, l1); s128_add(aim, h2, l2); break;    /* +p */
        case 1: { uint64_t nl = ~l2 + 1; int64_t nh = ~h2 + (nl == 0); s128_add(are, nh, nl); s128_add(aim, h1, l1); break; }   /* i p = (-im, re) */
        case 2: { uint64_t n1 = ~l1 + 1; int64_t m1 = ~h1 + (n1 == 0); uint64_t n2 = ~l2 + 1; int64_t m2 = ~h2 + (n2 == 0); s128_add(are, m1, n1); s128_add(aim, m2, n2); break; }
        default: { uint64_t n1 = ~l1 + 1; int64_t m1 = ~h1 + (n1 == 0); s128_add(are, h2, l2); s128_add(aim, m1, n1); break; }   /* -i p = (im, -re) */
    }
}

/* ---------------- Pauli strings ---------------- */
static inline int parity32(uint32_t v) { v ^= v >> 16; v ^= v >> 8; v ^= v >> 4; v ^= v >> 2; v ^= v >> 1; return v & 1; }
static int pop32(uint32_t v) { int c = 0; while (v) { v &= v - 1; c++; } return c; }
static pstr pmul(pstr a, pstr b) { pstr r; r.x = a.x ^ b.x; r.z = a.z ^ b.z; r.ph = (a.ph + b.ph + 2 * pop32(a.z & b.x)) & 3; return r; }
static pstr PX(int k) { pstr p = { 1u << k, 0, 0 }; return p; }
static pstr PY(int k) { pstr p = { 1u << k, 1u << k, 1 }; return p; }     /* Y = i X Z */
static pstr PZ(int k) { pstr p = { 0, 1u << k, 0 }; return p; }
static pstr PID(void) { pstr p = { 0, 0, 0 }; return p; }

static void papply(const pstr *p, const cq *in, cq *out)                 /* out = P in (bare, exact: a permutation with phases) */
{
    for (int b = 0; b < DIM; b++) {
        int rot = (p->ph + 2 * parity32(p->z & (uint32_t)b)) & 3; const cq v = in[b]; cq w;
        switch (rot) { case 0: w = v; break; case 1: w.re = -v.im; w.im = v.re; break; case 2: w.re = -v.re; w.im = -v.im; break; default: w.re = v.im; w.im = -v.re; }
        out[b ^ p->x] = w;
    }
}
static void acc_clear(void) { memset(AR, 0, sizeof(s128) * DIM); memset(AI, 0, sizeof(s128) * DIM); }
static void acc_terms(const pstr *T, const q *C, int nt, const cq *in)   /* (AR, AI) += SUM_terms C_t P_t in */
{
    for (int t = 0; t < nt; t++) {
        const uint32_t x = T[t].x, z = T[t].z; const int ph = T[t].ph; const q c = C[t];
        for (int b = 0; b < DIM; b++) acc_mac(&AR[b ^ x], &AI[b ^ x], c, in[b].re, in[b].im, ph + 2 * parity32(z & (uint32_t)b));
    }
}
static void acc_round(cq *out) { for (int b = 0; b < DIM; b++) { out[b].re = s128_round(AR[b]); out[b].im = s128_round(AI[b]); } }
static void apply_H(const cq *in, cq *out) { acc_clear(); acc_terms(HL, HC, NT, in); acc_terms(HR, HC, NT, in); acc_round(out); }

/* ---------------- vectors ---------------- */
static q norm2(const cq *v)
{
    s128 a = { 0, 0 }; int64_t h; uint64_t l;
    for (int b = 0; b < DIM; b++) { mulprod(v[b].re, v[b].re, &h, &l); s128_add(&a, h, l); mulprod(v[b].im, v[b].im, &h, &l); s128_add(&a, h, l); }
    return s128_round(a);
}
static cq inner(const cq *a, const cq *b)                               /* <a|b> = SUM conj(a) b, exact sum, one rounding */
{
    s128 re = { 0, 0 }, im = { 0, 0 }; int64_t h; uint64_t l;
    for (int k = 0; k < DIM; k++) {
        mulprod(a[k].re, b[k].re, &h, &l); s128_add(&re, h, l); mulprod(a[k].im, b[k].im, &h, &l); s128_add(&re, h, l);
        mulprod(a[k].re, b[k].im, &h, &l); s128_add(&im, h, l); mulprod(-a[k].im, b[k].re, &h, &l); s128_add(&im, h, l);
    }
    cq r = { s128_round(re), s128_round(im) }; return r;
}
static q maxabs(const cq *v) { q m = 0; for (int b = 0; b < DIM; b++) { q a = v[b].re < 0 ? -v[b].re : v[b].re, c = v[b].im < 0 ? -v[b].im : v[b].im; if (a > m) m = a; if (c > m) m = c; } return m; }
static void normalize(cq *v) { q n = sqrtq(norm2(v)); for (int b = 0; b < DIM; b++) { v[b].re = divq(v[b].re, n); v[b].im = divq(v[b].im, n); } }

/* e^{s H} v in place, s = (sre, sim) a complex Q.40 multiplier for ONE Taylor step; term and u are work buffers */
static void taylor_step(cq *v, q sre, q sim, cq *term, cq *u)
{
    for (int b = 0; b < DIM; b++) { AR[b].hi = 0; AR[b].lo = 0; AI[b].hi = 0; AI[b].lo = 0; s128_addq(&AR[b], v[b].re); s128_addq(&AI[b], v[b].im); }
    memcpy(term, v, sizeof(cq) * DIM);
    int k, small = 0;
    for (k = 1; k <= 80; k++) {
        /* u = H term, using a private accumulator pass: we need AR/AI for the sum, so accumulate H term into the term buffer via a second set */
        static s128 BR[MAXDIM], BI[MAXDIM];
        memset(BR, 0, sizeof(s128) * DIM); memset(BI, 0, sizeof(s128) * DIM);
        for (int t = 0; t < NT; t++) {
            const uint32_t x = HL[t].x, z = HL[t].z; const int ph = HL[t].ph; const q c = HC[t];
            for (int b = 0; b < DIM; b++) acc_mac(&BR[b ^ x], &BI[b ^ x], c, term[b].re, term[b].im, ph + 2 * parity32(z & (uint32_t)b));
        }
        for (int t = 0; t < NT; t++) {
            const uint32_t x = HR[t].x, z = HR[t].z; const int ph = HR[t].ph; const q c = HC[t];
            for (int b = 0; b < DIM; b++) acc_mac(&BR[b ^ x], &BI[b ^ x], c, term[b].re, term[b].im, ph + 2 * parity32(z & (uint32_t)b));
        }
        for (int b = 0; b < DIM; b++) { u[b].re = s128_round(BR[b]); u[b].im = s128_round(BI[b]); }
        /* term = s u / k */
        for (int b = 0; b < DIM; b++) {
            q re = mulq(u[b].re, sre) - mulq(u[b].im, sim), im = mulq(u[b].re, sim) + mulq(u[b].im, sre);
            term[b].re = rdiv(re, k); term[b].im = rdiv(im, k);
            s128_addq(&AR[b], term[b].re); s128_addq(&AI[b], term[b].im);
        }
        if (maxabs(term) <= 1) { if (++small >= 2) break; } else small = 0;
    }
    if (k > taylor_max) taylor_max = k;
    for (int b = 0; b < DIM; b++) { v[b].re = s128_round(AR[b]); v[b].im = s128_round(AI[b]); }
}
/* e^{-i H t} v (sign = -1) or e^{+i H t} v (sign = +1), t >= 0 in Q.40, in steps of 1 then the remainder; the norm drift is tracked */
static void evolve(cq *v, q t, int sign, cq *term, cq *u)
{
    q n0 = norm2(v);
    while (t > 0) {
        q dt = t >= ONE ? ONE : t; t -= dt;
        taylor_step(v, 0, sign > 0 ? dt : -dt, term, u);                 /* s = +- i dt */
    }
    q n1 = norm2(v), d = n1 - n0; if (d < 0) d = -d; if (d > drift_max) drift_max = d;
}

/* ---------------- the model ---------------- */
static void build_majoranas(void)
{
    for (int j = 1; j <= N / 2; j++) {
        int f = 2 * (j - 1), s = f + 1; pstr str = PID();
        for (int jp = 1; jp < j; jp++) { str = pmul(str, PZ(2 * (jp - 1))); str = pmul(str, PX(2 * (jp - 1) + 1)); }   /* (ZX)^{j-1} */
        ML[2 * j - 1] = pmul(pmul(str, PX(f)), PX(s));                 /* XX */
        ML[2 * j]     = pmul(pmul(str, PY(f)), PX(s));                 /* YX */
        MR[2 * j - 1] = pmul(str, PY(s));                              /* IY */
        MR[2 * j]     = pmul(str, PZ(s));                              /* IZ */
    }
}
static void build_I(cq *v)                                              /* the product over blocks of (|00> - |01> + |10> + |11>)/2 */
{
    for (int b = 0; b < DIM; b++) {
        int sgn = 1;
        for (int j = 0; j < N / 2; j++) { int a = (b >> (2 * j)) & 1, c = (b >> (2 * j + 1)) & 1; if (a == 0 && c == 1) sgn = -sgn; }
        v[b].re = sgn * (ONE >> (N / 2)); v[b].im = 0;
    }
}
static uint64_t rng;
static uint64_t xs64(void) { rng ^= rng >> 12; rng ^= rng << 25; rng ^= rng >> 27; return rng * 0x2545F4914F6CDD1Dull; }   /* xorshift64* */
static q gauss40(void) { q s = 0; for (int i = 0; i < 12; i++) s += (q)(xs64() >> 24); return s - 6 * ONE; }   /* 12 uniforms: mean 0, variance 1 */

static void build_syk(uint64_t seed)
{
    rng = seed ? seed : 0x9E3779B97F4A7C15ull;
    q var = divq(6 * ONE, (q)N * N * N * ONE), sd = sqrtq(var);      /* <J^2> = 3!/N^3 (J = 1) */
    NT = 0;
    for (int i = 1; i <= N; i++) for (int j = i + 1; j <= N; j++) for (int k = j + 1; k <= N; k++) for (int l = k + 1; l <= N; l++) {
        q J = mulq(gauss40(), sd);
        HL[NT] = pmul(pmul(pmul(ML[i], ML[j]), ML[k]), ML[l]);         /* the bare 4-product; the (1/sqrt2)^4 = 1/4 and i^{q/2} = -1 go in the coefficient */
        HR[NT] = pmul(pmul(pmul(MR[i], MR[j]), MR[k]), MR[l]);
        HC[NT] = -rdiv(J, 4);
        NT++;
    }
    NM = N;
}
static void build_model1(void)                                          /* Kobrin-Schuster-Yao eq. (1): 7 Majoranas, five terms, {psi, psi} = delta */
{
    static const int idx[5][4] = { {1, 2, 4, 5}, {1, 3, 4, 7}, {1, 3, 5, 6}, {2, 3, 4, 6}, {2, 3, 5, 7} };
    static const double cf[5] = { -0.36, 0.19, -0.71, 0.22, 0.49 };
    NT = 0;
    for (int t = 0; t < 5; t++) {
        double c = cf[t] * (double)ONE / 4.0; q cq40 = (q)(c < 0 ? c - 0.5 : c + 0.5);   /* the coefficient over the four 1/sqrt2, declared */
        HL[NT] = pmul(pmul(pmul(ML[idx[t][0]], ML[idx[t][1]]), ML[idx[t][2]]), ML[idx[t][3]]);
        HR[NT] = pmul(pmul(pmul(MR[idx[t][0]], MR[idx[t][1]]), MR[idx[t][2]]), MR[idx[t][3]]);
        HC[NT] = cq40; NT++;
    }
    NM = 7;
}

/* U = PROD_j (c + i s W_j), W_j = 2i psi_l^j psi_r^j = i * (bare_l^j bare_r^j); sign = +1 for U, -1 for U^dag. v in place, w a work buffer */
static void apply_U(cq *v, q c, q s, int sign, cq *w, int vnormN)
{
    (void)vnormN;
    for (int j = 1; j <= NM; j++) {
        pstr W = pmul(ML[j], MR[j]); W.ph = (W.ph + 1) & 3;               /* times i */
        papply(&W, v, w);
        q ss = sign > 0 ? s : -s;
        for (int b = 0; b < DIM; b++) {                                  /* v = c v + i ss (W v) */
            q wr = w[b].re, wi = w[b].im;
            q re = mulq(v[b].re, c) - mulq(wi, ss), im = mulq(v[b].im, c) + mulq(wr, ss);
            v[b].re = re; v[b].im = im;
        }
    }
}

static void hashq(q v) { chain = fnv(chain, &v, 8); }
static q parseq(const char *s)                                           /* "-12", "2.8", "0.5" -> Q.40 (decimal part rounded) */
{
    int neg = 0; if (*s == '-') { neg = 1; s++; }
    q ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (*s++ - '0');
    q fp = 0, scale = 1; if (*s == '.') { s++; while (*s >= '0' && *s <= '9' && scale < 1000000000000ll) { fp = fp * 10 + (*s++ - '0'); scale *= 10; } }
    q v = ip * ONE + (fp ? divq(fp * ONE, scale * ONE) : 0);
    return neg ? -v : v;
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
    if (N < 4 || N > MAXN || (N & 1)) { printf("WHSYK: N must be even, 4..%d\n", MAXN); return 1; }
    if (!nmu) { mus[0] = -12 * ONE; mus[1] = 12 * ONE; nmu = 2; }
    if (!nt) { ts[0] = ONE; ts[1] = 2 * ONE; ts[2] = 3 * ONE; nt = 3; }
    if (!nseed) { seeds[0] = 1; nseed = 1; }
    if (model) { N = 8; nseed = 1; if (!vnormN) vnormN = 10; }
    if (!vnormN) vnormN = N;
    DIM = 1 << N;
    char A[32], B[32], C[32], D[32], E[32], F[32], G[32];
    fmt(beta, A); fmt(mus[0], B);
    printf("WHSYK: Q.%d (whq.h); %s; N = %d Majoranas a side on %d qubits (2^%d amplitudes); beta = %s; V normalized with N = %d; the Dirac fermion read = (psi^%d + i psi^%d)/sqrt2; Taylor dt = 1\n",
           FB, model ? "MODEL 1 of the 2022 experiment (Kobrin-Schuster-Yao eq. 1: 7 Majoranas, 5 terms)" : "the FULL q = 4 SYK, <J^2> = 6/N^3, J = 1", N, N, N, A, vnormN, pa, pb);
    build_majoranas();
    cq *I0 = V[0], *T = V[1], *w = V[2], *a1 = V[3], *a2 = V[4], *a21 = V[5], *a12 = V[6], *up = V[7], *b1 = V[8], *b2 = V[9], *bb = V[10], *s1 = V[11], *s2 = V[12], *s3 = V[13];

    /* self-test 1: (bare_l^j + i bare_r^j)|I> = 0 exactly, for every j */
    build_I(I0);
    q worst = 0;
    for (int j = 1; j <= N; j++) {
        papply(&ML[j], I0, s1); papply(&MR[j], I0, s2);
        for (int b = 0; b < DIM; b++) { q re = s1[b].re - s2[b].im, im = s1[b].im + s2[b].re; if (re < 0) re = -re; if (im < 0) im = -im; if (re > worst) worst = re; if (im > worst) worst = im; }
    }
    /* self-test 2: {bare_a, bare_b} = 2 delta_ab on |I> for a few pairs (l-l, r-r, l-r) */
    q worst2 = 0; int pairs[6][2] = { {1, 1}, {1, 2}, {2, 3}, {1, N}, {N, N}, {2, 2} };
    for (int p = 0; p < 6; p++) for (int side = 0; side < 3; side++) {
        pstr *Pa = side == 2 ? ML : (side ? MR : ML), *Pb = side == 2 ? MR : (side ? MR : ML); int a = pairs[p][0], c = pairs[p][1];
        papply(&Pa[a], I0, s1); papply(&Pb[c], s1, s2); papply(&Pb[c], I0, s1); papply(&Pa[a], s1, s3);
        q expect = (side != 2 && a == c) ? 2 * (ONE >> (N / 2)) : 0;
        for (int b = 0; b < DIM; b++) { q re = s2[b].re + s3[b].re - (I0[b].re > 0 ? expect : -expect) * (I0[b].re != 0), im = s2[b].im + s3[b].im; if (re < 0) re = -re; if (im < 0) im = -im; if (re > worst2) worst2 = re; if (im > worst2) worst2 = im; }
    }
    fmt(worst, A); fmt(worst2, B);
    printf("WHSYK selftest: max |(psi_l^j + i psi_r^j)|I>| = %s (must be 0.000000000000); max |{bare_a, bare_b} - 2 delta| on |I> = %s (must be 0)\n", A, B);

    for (int si = 0; si < nseed; si++) {
        if (model) build_model1(); else build_syk(seeds[si]);
        /* self-test 3: (H_l - H_r)|I> = 0 exactly */
        acc_clear(); acc_terms(HL, HC, NT, I0); acc_round(s1); acc_clear(); acc_terms(HR, HC, NT, I0); acc_round(s2);
        q worst3 = 0; for (int b = 0; b < DIM; b++) { q re = s1[b].re - s2[b].re, im = s1[b].im - s2[b].im; if (re < 0) re = -re; if (im < 0) im = -im; if (re > worst3) worst3 = re; if (im > worst3) worst3 = im; }
        q hnorm = norm2(s1); fmt(worst3, A); fmt(sqrtq(hnorm), B); fmt(HC[0], C);
        printf("WHSYK seed %llu: %d terms, first coefficient %s; max |(H_l - H_r)|I>| = %s (must be 0); |H_l I| = %s\n", (unsigned long long)seeds[si], NT, C, A, B);
        hashq(HC[0]); hashq(HC[NT - 1]);
        /* the thermofield double: e^{-beta H/4} |I>, normalized */
        memcpy(T, I0, sizeof(cq) * DIM);
        q bq = rdiv(beta, 4); drift_max = 0;
        while (bq > 0) { q d = bq >= ONE ? ONE : bq; bq -= d; taylor_step(T, -d, 0, s1, s2); }
        q nT = norm2(T); normalize(T); fmt(nT, A); fmt(norm2(T), B);
        printf("WHSYK seed %llu: |e^{-beta H/4} I|^2 = %s (Z_beta / Z_inf); normalized to %s; Taylor terms used so far at most %d\n", (unsigned long long)seeds[si], A, B, taylor_max);
        hashq(nT);
        for (int mi = 0; mi < nmu; mi++) {
            q arg = divq(mus[mi], (q)(2 * 4 * vnormN) * ONE), c = cosq(arg), s = sinq(arg);
            for (int ti = 0; ti < nt; ti++) {
                q t = ts[ti], tp = ntp ? tps[ti < ntp ? ti : ntp - 1] : t; drift_max = 0;
                /* the A side: w = e^{iHt'} T; a_i = e^{-iHt'} m_l^i w; a21 = A2 a1; a12 = A1 a2 */
                memcpy(w, T, sizeof(cq) * DIM); evolve(w, tp, +1, s1, s2);
                papply(&ML[pa], w, a1); evolve(a1, tp, -1, s1, s2);
                papply(&ML[pb], w, a2); evolve(a2, tp, -1, s1, s2);
                memcpy(s3, a1, sizeof(cq) * DIM); evolve(s3, tp, +1, s1, s2); papply(&ML[pb], s3, a21); evolve(a21, tp, -1, s1, s2);
                memcpy(s3, a2, sizeof(cq) * DIM); evolve(s3, tp, +1, s1, s2); papply(&ML[pa], s3, a12); evolve(a12, tp, -1, s1, s2);
                /* the B side on T: up = e^{-iHt} U T; b_i = U^dag e^{iHt} m_r^i up; bb = U^dag e^{iHt} m_r^a m_r^b up */
                memcpy(up, T, sizeof(cq) * DIM); apply_U(up, c, s, +1, s1, vnormN); evolve(up, t, -1, s1, s2);
                papply(&MR[pa], up, b1); evolve(b1, t, +1, s1, s2); apply_U(b1, c, s, -1, s1, vnormN);
                papply(&MR[pb], up, b2); evolve(b2, t, +1, s1, s2); apply_U(b2, c, s, -1, s1, vnormN);
                papply(&MR[pb], up, s3); papply(&MR[pa], s3, bb); evolve(bb, t, +1, s1, s2); apply_U(bb, c, s, -1, s1, vnormN);
                cq Y1 = inner(a1, b1);                                   /* <a1|b1>: K = Re Y1 (the 1/2 of two Majoranas times 2 Re) */
                /* B12 a1 and B2 a1 share e^{-iHt} U a1 */
                memcpy(up, a1, sizeof(cq) * DIM); apply_U(up, c, s, +1, s1, vnormN); evolve(up, t, -1, s1, s2);
                papply(&MR[pb], up, s3); papply(&MR[pa], s3, w); evolve(w, t, +1, s1, s2); apply_U(w, c, s, -1, s1, vnormN);   /* w = B12 a1 */
                cq X2b = inner(a2, w);                                   /* <a2|B12 a1> = -conj<a1|B12 a2>: used as Re X2 = -Re<a2|B12 a1> */
                papply(&MR[pb], up, s3); evolve(s3, t, +1, s1, s2); apply_U(s3, c, s, -1, s1, vnormN);                        /* s3 = B2 a1 */
                cq Y2 = inner(a12, s3);                                  /* <a12|B2 a1> */
                /* B12 a2 */
                memcpy(up, a2, sizeof(cq) * DIM); apply_U(up, c, s, +1, s1, vnormN); evolve(up, t, -1, s1, s2);
                papply(&MR[pb], up, s3); papply(&MR[pa], s3, w); evolve(w, t, +1, s1, s2); apply_U(w, c, s, -1, s1, vnormN);   /* w = B12 a2 */
                cq X2 = inner(a1, w);                                    /* <a1|B12 a2> */
                /* B2 a12 */
                memcpy(up, a12, sizeof(cq) * DIM); apply_U(up, c, s, +1, s1, vnormN); evolve(up, t, -1, s1, s2);
                papply(&MR[pb], up, s3); evolve(s3, t, +1, s1, s2); apply_U(s3, c, s, -1, s1, vnormN);                        /* s3 = B2 a12 */
                cq Y3 = inner(a1, s3);                                   /* <a1|B2 a12> */
                cq X1 = inner(a21, bb);                                  /* <a21|B12 T> */
                (void)X2b;
                /* <{A1,[A2,B12]}> = (1/4) (2 Re X1 - 2 Re X2); rho11 = (1 - that)/2 = 1/2 - (Re X1 - Re X2)/4 */
                q rho11 = (ONE >> 1) - rdiv(X1.re - X2.re, 4);
                /* rho14 = (1/2)<{A1,B1}> - <{A2, A1 B2 A1}> = (1/2)(Re Y1) - (1/4)(Y2 + Y3) */
                q K = Y1.re;
                q r14re = rdiv(Y1.re, 2) - rdiv(Y2.re + Y3.re, 4), r14im = -rdiv(Y2.im + Y3.im, 4);
                q m14 = sqrtq(mulq(r14re, r14re) + mulq(r14im, r14im));
                q l1 = rdiv(rho11 + m14, 2), l2 = rdiv(rho11 - m14, 2), l3 = rdiv(ONE - rho11, 2);
                q S = 0; int clamp = 0; q lam[4] = { l1, l2, l3, l3 };
                for (int e = 0; e < 4; e++) { if (lam[e] <= 0) { clamp = 1; continue; } S -= mulq(lam[e], lnq(lam[e])); }
                q I = 2 * LN2 - S, frac = divq(I, 2 * LN2);
                fmt(t, A); fmt(mus[mi], B); fmt(rho11, C); fmt(r14re, D); fmt(r14im, E); fmt(K, F); fmt(frac, G);
                char HB[32], HT[32], HD[32]; fmt(beta, HB); fmt(tp, HT); fmt(drift_max, HD);   /* v2 fix: one buffer per field (v1 shared one and printed beta thrice) */
                printf("WHSYK N=%d seed=%llu beta=%s mu=%s t=%s tp=%s: rho11 %s | rho14 %s %s i | K %s | I_RT/(2 ln 2) %s%s | norm drift %s | Taylor terms <= %d\n",
                       N, (unsigned long long)seeds[si], HB, B, A, HT, C, D, E, F, G, clamp ? " (an eigenvalue clamped at 0)" : "", HD, taylor_max);
                hashq(rho11); hashq(r14re); hashq(r14im); hashq(K); hashq(frac);
            }
        }
    }
    printf("WHSYK: chain %016llx\n", (unsigned long long)chain);
    printf("WHSYK: done\n");
    return 0;
}
