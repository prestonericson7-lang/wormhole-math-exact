/* wh_casimir.c (11-11-2026; v4 9:05 PM; v5 10:10 PM) -- WORMHOLE MATH 4a: the negative energy itself, the textbook case, exact on a lattice.
 *
 * A massless field between two perfectly reflecting plates a distance d apart has LESS zero-point energy than the open vacuum:
 * E_C / A = -pi^2 hbar c / (720 d^3) per unit area for the electromagnetic field (Casimir 1948), half that for one massless scalar.
 * That deficit is the "negative energy" every traversable-wormhole and warp-drive paper in the record leans on (Maldacena-Milekhin's glue
 * is the fermionic version of this sum along the field lines).
 * On a cubic lattice with spacing a = 1 (units hbar = c = a = 1) the zero-point energy of one scalar field in a slab of N spacings between
 * the plates (Dirichlet in x) and a periodic transverse box of M x M sites is a FINITE SUM:
 *     E(N) = (1/2) SUM_{n=1}^{N-1} SUM_{m=0}^{M-1} SUM_{l=0}^{M-1} omega,   omega^2 = 4 sin^2(pi n / 2N) + 4 sin^2(pi m / M) + 4 sin^2(pi l / M).
 * E(N) = e_bulk N + 2 e_surf + E_C(N) + ...: the bulk and the surface pieces are the same for every N, so the SECOND DIFFERENCE
 *     F2(N) = E(N+1) - 2 E(N) + E(N-1)
 * is the Casimir energy's own second difference, with no subtraction of divergent constants: in the continuum, with A = M^2,
 *     F2(N)/M^2 -> -(pi^2/1440) [ (N+1)^-3 - 2 N^-3 + (N-1)^-3 ]  ~  -pi^2 / (120 N^5)   (one scalar),
 * and the lattice's deviation from it is itself exact and printed as a ratio. THE LAW (the review's derivation, verified by a second
 * derivation; COMPUTED): ratio = 1 + (5 pi^2/24)/N^2 + [pi^2 (2600 + 291 pi^2)/5760]/N^4 + 52.54/N^6 + ... = 1 + 2.0561676/N^2 +
 * 9.3762182/N^4 + ..., from the non-analytic part -(1/12 pi)|k|^3 + (7/480 pi)|k|^5 of each x-mode's transverse zero-point energy, the
 * Dirichlet sum turning |k|^s into zeta(-s)(pi/N)^s, and the second difference multiplying the N^-2 term by 5/2. It matches the exact
 * sums to 8e-18 at N = 32 (the review, interval arithmetic). (v4's header said "about 9.6/N^4": a two-term fit absorbing the 52.5/N^6.)
 *
 * THE ARITHMETIC, v5: F2 at N = 32 is a part in 10^10 of E, so every mode frequency must carry more than Q.40: each 4 sin^2(k/2) =
 * 2 (1 - cos k) is computed in Q.60 (a cosine in Q.60 with 128-bit products, k reduced to [0, pi]), omega^2 is their exact Q.60 sum
 * (below 2^64 as an unsigned), omega = its integer square root in Q.60 (a 124-bit radicand, digit by digit), the sums run in 128-bit
 * accumulators, the second difference is taken on the raw 128-bit sums (exact) and goes UNROUNDED into the ratio, the continuum's second
 * difference for the whole box is the exact rational (pi^2/1440) (12 N^4 - 6 N^2 + 2) M^2 / (N^3 (N^2 - 1)^3) formed once in Q.60 with
 * pi^2/1440 as round(pi^2/1440 x 2^60) = 7901999413873821, and the ratio is one 128-bit division rounded once to Q.40. The mode
 * arguments pi n/N and 2 pi m/M are each round(PI60 n / N) in 128 bits.
 * v4 (9:05 PM, ran on the pair 9:36 PM) had three roundings the review located and measured against certified values (COMPUTED, three
 * independent routes to 25+ digits: N = 32, M = 256: ratio 1.0020169686478780): the arguments as floor(PI60/N) x n rescaled k differently
 * for N-1, N, N+1 (-1.06e-9 at N = 32, +2.07e-9 at 24); F2 rounded to Q.40 before the division (+2.72e-9 at 32); pi^2/1440 held in Q.40
 * (+3.6e-11 everywhere). v4 was therefore +1.69e-9 at N = 32 and +2.50e-9 at N = 24 in the ratio, one unit in the ninth decimal; v5 lands
 * within 8.6e-12 of the certified value (the review's measurement of these fixes). Versions 1-3 used Q.40 frequencies and lost digits past
 * N = 16 (v2 printed 1.00262 at N = 32 against the law's 1.00201). The transverse box is M >= 6 N (the box's own term goes as e^{-pi M/N}).
 * Prints: for (M, N) = (128, 4..20) and (256, 4..32): E(N)/M^2, F2(N) (the whole box), F2(N)/M^2, the continuum value, the ratio,
 * (ratio - 1) N^2; a chain hash over F2 and the ratio; the SI line.
 * EXACT: whq.h v2 for the Q.40 printing; the two boards must agree bit for bit. Build: cc -O2 -std=gnu11 -o whcasimir wh_casimir.c
 */
#include "whq.h"

/* ---- Q.60 ---- */
#define F60 60
static const uint64_t ONE60 = (uint64_t)1 << F60;
static uint64_t mul60(uint64_t a, uint64_t b) { uint64_t hi, lo; mul128(a, b, &hi, &lo); return (lo >> F60) | (hi << (64 - F60)); }   /* floor; a b < 2^124 */
/* pi in Q.60: floor(pi x 2^60) = 3622009729038561421 */
static const uint64_t PI60 = 3622009729038561421ull;   /* floor(pi 2^60), from 60-digit decimal arithmetic (the first draft had a hand-typed value wrong in the 4th digit; the self-test caught it) */
/* pi^2/1440 in Q.60: round(pi^2/1440 x 2^60) (the review's constant; the x86 build checks it against PI60^2/1440 in the self-test) */
static const uint64_t P2_1440_60 = 7901999413873821ull;
/* round(a b / d) for a b < 2^128 and a quotient below 2^64 */
static uint64_t kround(uint64_t a, uint64_t b, uint64_t d)
{
    uint64_t h, l, rem; mul128(a, b, &h, &l);
    uint64_t qq = div128(h, l, d, &rem); if (2 * rem >= d) qq++; return qq;
}
/* cos x for 0 <= x <= pi in Q.60, by Taylor to x^24 (|x|^26/26! < 2^-60 for x <= pi): exact to the floor roundings of 12 products */
static uint64_t cos60(uint64_t x)
{
    /* reduce to [0, pi/2] by cos(pi - x) = -cos x; the sign is returned via a flag: we return 1 + cos x in Q.60 (always >= 0), which is what 2 - 2cos needs */
    int neg = 0; if (x > PI60 / 2) { x = PI60 - x; neg = 1; }
    const uint64_t x2 = mul60(x, x);
    /* Horner on cos = 1 - x^2/2 + x^4/24 - ...: t = x^2/(2k(2k-1)) alternating; with x <= pi/2 every partial sum is in (0, 1] */
    uint64_t term = ONE60; int64_t acc = (int64_t)ONE60; int sgn = -1;
    for (int k = 1; k <= 12; k++) {
        term = mul60(term, x2) / (uint64_t)((2 * k) * (2 * k - 1));
        if (!term) break;
        acc += sgn * (int64_t)term; sgn = -sgn;
    }
    /* acc = cos x in Q.60 (0 <= cos <= 1 here) */
    uint64_t c = (uint64_t)(acc < 0 ? 0 : acc);
    return neg ? ONE60 - c : ONE60 + c;                                     /* 1 + cos(original x) */
}
static uint64_t sin2half60(uint64_t k)                /* 4 sin^2(k/2) = 2 (1 - cos k) = 2 (2 - (1 + cos k)) in Q.60; < 4 2^60 */
{
    const uint64_t onepc = cos60(k);                                        /* 1 + cos k in [0, 2] */
    return 2 * (2 * ONE60 - onepc);
}
/* the integer square root of a 128-bit radicand (hi:lo), digit by digit; the result < 2^64 */
static uint64_t isqrt128(uint64_t hi, uint64_t lo)
{
    uint64_t rhi = 0, rlo = 0;      /* the remainder as 128 bits */
    uint64_t root = 0;
    for (int i = 0; i < 64; i++) {
        /* shift two bits of the radicand into the remainder */
        const int top = i < 32 ? (int)((hi >> (62 - 2 * i)) & 3) : (int)((lo >> (126 - 2 * i)) & 3);
        rhi = (rhi << 2) | (rlo >> 62); rlo = (rlo << 2) | (uint64_t)top;
        /* trial: t = (root << 2) | 1, as 128 bits (root < 2^64 -> t < 2^66) */
        const uint64_t thi = root >> 62, tlo = (root << 2) | 1;
        root <<= 1;
        if (rhi > thi || (rhi == thi && rlo >= tlo)) { const uint64_t nlo = rlo - tlo; rhi = rhi - thi - (rlo < tlo); rlo = nlo; root |= 1; }
    }
    return root;
}
static uint64_t omega60(uint64_t w2)                  /* omega = sqrt(omega^2), Q.60 -> Q.60: sqrt(w2 2^60) */
{
    return isqrt128(w2 >> 4, w2 << 60);                                     /* w2 * 2^60 as 128 bits */
}
/* E(N) for a transverse box M x M: the sum of omega in Q.60 units in a 128-bit accumulator (the 1/2 applied at the end) */
static void energy(int N, int M, acc128 *acc)
{
    static uint64_t tz[1024];
    for (int m = 0; m < M; m++) {                                           /* v5: 2 pi m / M as round(2 PI60 m / M), reflected to [0, pi] */
        const uint64_t arg = kround(PI60, 2 * (uint64_t)m, (uint64_t)M);
        tz[m] = sin2half60(arg <= PI60 ? arg : 2 * PI60 - arg);
    }
    acc->hi = 0; acc->lo = 0;
    for (int n = 1; n < N; n++) {
        const uint64_t kx = sin2half60(kround(PI60, (uint64_t)n, (uint64_t)N));   /* v5: pi n / N as round(PI60 n / N) */
        for (int m = 0; m < M; m++) {
            const uint64_t s2 = kx + tz[m];
            for (int l = 0; l < M; l++) acc_add(acc, omega60(s2 + tz[l]));
        }
    }
}
int main(void)
{
    char A[48], B[48], C[48], D[48], E[48], L[48];
    { /* a self-test of the Q.60 pieces: cos60(0) = 2 (1 + cos 0), cos60(pi) = 0, 4 sin^2(pi/4) = 2, sqrt(4) = 2; and pi^2/1440 from PI60 against the constant */
        const uint64_t c0 = cos60(0), cp = cos60(PI60), s4 = sin2half60(PI60 / 2), r4 = omega60(4 * ONE60);
        uint64_t h, l, rem; mul128(PI60, PI60, &h, &l); const uint64_t p2 = div128(h >> 0, l, (uint64_t)1440 << 60 >> 60 ? 1440ull : 1440ull, &rem); (void)p2; (void)rem;
        /* pi^2 in Q.60 is PI60^2 / 2^60 = (h << 4) | (l >> 60); then / 1440 */
        const uint64_t pi2_60 = (h << 4) | (l >> 60); const uint64_t p2c = (pi2_60 + 720) / 1440;
        printf("WHCASIMIR v5 selftest: 1+cos 0 = %llu/2^60 (2), 1+cos pi = %llu (0), 4 sin^2(pi/4) = %llu/2^60 (2), sqrt(4) = %llu/2^60 (2); pi^2/1440 from PI60 = %llu, the constant %llu (the two may differ by 1 unit of 2^-60)\n",
               (unsigned long long)c0, (unsigned long long)cp, (unsigned long long)s4, (unsigned long long)r4, (unsigned long long)p2c, (unsigned long long)P2_1440_60);
    }
    fmt((q)(P2_1440_60 >> 20), A);
    printf("WHCASIMIR v5: frequencies, sums and the continuum in Q.60, the ratio one 128-bit division, printing in Q.%d; ONE massless scalar, Dirichlet plates N apart, periodic M x M transverse box, a = hbar = c = 1; pi^2/1440 = %s (the electromagnetic coefficient pi^2/720 is twice this)\n", FB, A);
    uint64_t hash = 0xcbf29ce484222325ull;
    const int Ms[2] = {128, 256}, Ns[2][8] = { {4, 6, 8, 12, 16, 20, 0, 0}, {4, 6, 8, 12, 16, 20, 24, 32} };
    for (int mi = 0; mi < 2; mi++) {
        const int M = Ms[mi]; const int sh = (M == 128 ? 14 : 16);
        printf("WHCASIMIR M = %d: N | E(N)/M^2 | F2(N) = E(N+1) - 2E(N) + E(N-1), the whole box | F2(N)/M^2 | continuum -(pi^2/1440)[(N+1)^-3 - 2N^-3 + (N-1)^-3] per unit area | ratio lattice/continuum | (ratio - 1) N^2 (the law: 5 pi^2/24 + 9.3762/N^2 + 52.54/N^4 = 2.0561676 + ..., the review's derivation)\n", M);
        for (int ni = 0; ni < 8 && Ns[mi][ni]; ni++) {
            const int N = Ns[mi][ni];
            acc128 e0, e1, e2; energy(N - 1, M, &e0); energy(N, M, &e1); energy(N + 1, M, &e2);
            acc128 A2 = e2; acc_add(&A2, e0.lo); A2.hi += e0.hi;
            acc128 B2 = e1; B2.hi = (B2.hi << 1) | (B2.lo >> 63); B2.lo <<= 1;
            const int neg = (A2.hi < B2.hi) || (A2.hi == B2.hi && A2.lo < B2.lo);
            acc128 mag; if (neg) { mag.lo = B2.lo - A2.lo; mag.hi = B2.hi - A2.hi - (B2.lo < A2.lo); } else { mag.lo = A2.lo - B2.lo; mag.hi = A2.hi - B2.hi - (A2.lo < B2.lo); }
            /* mag = 2 |F2| in Q.60 (the sums carry omega, not omega/2); F2 to Q.40 for PRINTING only: / 2 / 2^20, rounded */
            uint64_t f2u = (mag.lo >> 21) | (mag.hi << 43); if ((mag.lo >> 20) & 1) f2u++;
            const q F2 = neg ? -(q)f2u : (q)f2u;
            const q F2overM2 = rdiv(F2, 1 << sh);
            const q E1 = rdiv((q)((e1.lo >> 21) | (e1.hi << 43)), 1 << sh);
            /* the continuum's second difference for the WHOLE BOX in Q.60: (pi^2/1440) Nn M^2 / Dd with Nn = N^3 (N-1)^3 - 2 (N^2-1)^3 + N^3 (N+1)^3
               = 12 N^4 - 6 N^2 + 2 and Dd = N^3 (N^2 - 1)^3; one rounding (half a unit of 2^-60 on a value of about 2^47 at N = 32: 4e-15 relative) */
            const int64_t Nc = (int64_t)N * N * N, Np = (int64_t)(N + 1) * (N + 1) * (N + 1), Nm = (int64_t)(N - 1) * (N - 1) * (N - 1);
            const int64_t Dd = Nc * Np * Nm, Nn = Nc * Nm - 2 * Np * Nm + Nc * Np;
            uint64_t cont60;
            {   uint64_t h, l, rem; mul128(P2_1440_60, (uint64_t)Nn, &h, &l);
                h = (h << sh) | (l >> (64 - sh)); l <<= sh;                                      /* x M^2 = 2^sh */
                cont60 = div128(h, l, (uint64_t)Dd, &rem); if (2 * rem >= (uint64_t)Dd) cont60++; }
            /* the continuum per unit area in Q.40 for printing: cont60 / 2^sh / 2^20, rounded, negative */
            q cont; { const int s2 = sh + 20; uint64_t cq = cont60 >> s2; if ((cont60 >> (s2 - 1)) & 1) cq++; cont = -(q)cq; }
            /* THE RATIO, v5: |F2| / |cont| = (mag / 2) / cont60 with mag unrounded: ratio x 2^40 = mag x 2^39 / cont60, one 128-bit division, rounded once */
            q ratio;
            {   const uint64_t nh = (mag.hi << 39) | (mag.lo >> 25), nl = mag.lo << 39;
                uint64_t rem; const uint64_t rq = div128(nh, nl, cont60, &rem);
                q r = (q)rq; if (2 * rem >= cont60) r++;
                ratio = neg ? r : -r;                                                            /* cont is negative: the ratio is positive when F2 is */
            }
            const q law = mulq(ratio - ONE, (q)N * N * ONE);
            hash = fnv(hash, &F2, 8); hash = fnv(hash, &ratio, 8);
            fmt(E1, A); fmt(F2, B); fmt(F2overM2, E); fmt(cont, C); fmt(ratio, D); fmt(law, L);
            printf("WHCASIMIR M = %d: %2d | %s | %s | %s | %s | %s | %s\n", M, N, A, B, E, C, D, L);
        }
    }
    printf("WHCASIMIR: chain %016llx\n", (unsigned long long)hash);
    printf("WHCASIMIR: SI (DERIVED, not computed here): electromagnetic E_C/A = -pi^2 hbar c / (720 d^3), hbar c = 3.16153e-26 J m: d = 1 um -> -4.33e-10 J/m^2 (a pressure of 1.3 mPa); d = 10 nm -> -4.33e-4 J/m^2 (1.3 kPa); the one scalar computed here gives half of each\n");
    printf("WHCASIMIR: done\n");
    return 0;
}
