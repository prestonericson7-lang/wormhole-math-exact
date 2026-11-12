/* wh_casimir.c (11-11-2026; v2 7:35 PM after the math buddy's review) -- WORMHOLE MATH 4a: the negative energy itself, the textbook case,
 * exact on a lattice. v2: the continuum coefficient for ONE scalar field is pi^2/1440 (pi^2/720 is the electromagnetic value, two
 * polarizations); the second difference is taken in 128-bit before any shift (v1 shifted each E(N) to Q.40 first and lost digits); the
 * transverse box is M >= 6 N (the periodic box's own finite-size term goes as e^{-pi M/N}; M = 64 at N = 32 was 13% off); the ratio
 * divides by a negative number, which v1's divq could not do (whq.h v2 can).
 *
 * A massless scalar field between two perfectly reflecting plates a distance d apart has LESS zero-point energy than the open vacuum:
 * E_C / A = -pi^2 hbar c / (720 d^3) per unit area for the electromagnetic field (Casimir 1948), half that for one massless scalar. That deficit is the "negative energy" every traversable-wormhole and
 * warp-drive paper in the record leans on (Maldacena-Milekhin's glue is the fermionic version of this sum along the field lines).
 * On a cubic lattice with spacing a = 1 (units hbar = c = a = 1) the zero-point energy of the field in a slab of N spacings between the
 * plates (Dirichlet in x) and a periodic transverse box of M x M sites is a FINITE SUM:
 *     E(N) = (1/2) SUM_{n=1}^{N-1} SUM_{m=0}^{M-1} SUM_{l=0}^{M-1} omega,   omega^2 = 4 sin^2(pi n / 2N) + 4 sin^2(pi m / M) + 4 sin^2(pi l / M).
 * E(N) = e_bulk N + 2 e_surf + E_C(N) + ...: the bulk and the surface pieces are the same for every N, so the SECOND DIFFERENCE
 *     F2(N) = E(N+1) - 2 E(N) + E(N-1)
 * is the Casimir energy's own second difference, with no subtraction of divergent constants: in the continuum, with A = M^2,
 *     F2(N)/M^2 -> -(pi^2/1440) [ (N+1)^-3 - 2 N^-3 + (N-1)^-3 ]  ~  -pi^2 / (120 N^5)   (one scalar),
 * and the lattice's deviation from it (order 1/N^2) is itself exact and printed. Every term is a square root of sines in Q.40; the sums
 * run in 128-bit accumulators; the two boards must agree bit for bit. Rounding: sqrtq is a floor, so every term carries a bias of about
 * -1/2 unit; the bias is the same per term for every N and cancels in the second difference; what remains is the random part, about
 * 0.29 units x sqrt(terms) = 700 units at N = 32, M = 256 against F2 of 3.5e8 units: five digits of F2 at N = 32, more at small N.
 * Prints: for (M, N) = (128, 4..20) and (256, 4..32): E(N)/M^2, F2(N) (the whole box), F2(N)/M^2, the continuum value, their ratio.
 * EXACT: whq.h (Q.40, 128-bit intermediates, round to nearest). Build: cc -O2 -std=gnu11 -o whcasimir wh_casimir.c
 */
#include "whq.h"

static q sin2half(q k)                                /* 4 sin^2(k/2) = 2 (1 - cos k), exact to rounding */
{
    return 2 * (ONE - cosq(k));
}
/* E(N) for a transverse box M x M: returns the sum in a 128-bit accumulator of Q.40 units (the 1/2 applied at the end) */
static void energy(int N, int M, acc128 *acc)
{
    static q tz[1024];                                /* 4 sin^2(pi m / M) for m = 0..M-1 */
    for (int m = 0; m < M; m++) tz[m] = sin2half(divq(mulq(2 * PI, (q)m * ONE), (q)M * ONE));
    acc->hi = 0; acc->lo = 0;
    for (int n = 1; n < N; n++) {
        const q kx = sin2half(divq(mulq(PI, (q)n * ONE), (q)N * ONE));
        for (int m = 0; m < M; m++) {
            const q s2 = kx + tz[m];
            for (int l = 0; l < M; l++) acc_add(acc, (uint64_t)sqrtq(s2 + tz[l]));
        }
    }
}
int main(void)
{
    char A[48], B[48], C[48], D[48];
    const q pi2_1440 = divq(mulq(PI, PI), 1440 * ONE);
    fmt(pi2_1440, A);
    printf("WHCASIMIR v2: Q.%d (whq.h v2); ONE massless scalar, Dirichlet plates N apart, periodic M x M transverse box, a = hbar = c = 1; pi^2/1440 = %s (the electromagnetic coefficient pi^2/720 is twice this)\n", FB, A);
    uint64_t hash = 0xcbf29ce484222325ull;
    const int Ms[2] = {128, 256}, Ns[2][8] = { {4, 6, 8, 12, 16, 20, 0, 0}, {4, 6, 8, 12, 16, 20, 24, 32} };   /* M >= 6 N throughout */
    for (int mi = 0; mi < 2; mi++) {
        const int M = Ms[mi]; const int sh = (M == 128 ? 14 : 16);   /* M^2 as a power of two */
        printf("WHCASIMIR M = %d: N | E(N)/M^2 | F2(N) = E(N+1) - 2E(N) + E(N-1), the whole box | F2(N)/M^2 | continuum -(pi^2/1440)[(N+1)^-3 - 2N^-3 + (N-1)^-3] per unit area | ratio lattice/continuum\n", M);
        for (int ni = 0; ni < 8 && Ns[mi][ni]; ni++) {
            const int N = Ns[mi][ni];
            acc128 e0, e1, e2; energy(N - 1, M, &e0); energy(N, M, &e1); energy(N + 1, M, &e2);
            /* the raw sums are SUM omega in Q.40 units (the 1/2 not yet applied): below 2^63 for M <= 256, N <= 33, so the second difference is exact in int64 */
            if (e0.hi || e1.hi || e2.hi) { printf("WHCASIMIR: a sum overflowed 64 bits at M = %d N = %d -- stop\n", M, N); return 1; }
            const int64_t F2raw = (int64_t)e2.lo - 2 * (int64_t)e1.lo + (int64_t)e0.lo;   /* 2 F2 in Q.40 units, exact */
            const q F2 = rdiv(F2raw, 2);                                      /* F2 for the whole box, Q.40, one rounding */
            const q F2overM2 = rdiv(F2, 1 << sh);                             /* F2/M^2, a second rounding, for the table only */
            const q E1 = rdiv((q)(e1.lo >> 1), 1 << sh);                      /* E(N)/M^2 for the table */
            const q c1 = divq(ONE, (q)(N + 1) * (N + 1) * (N + 1) * ONE), c0 = divq(ONE, (q)N * N * N * ONE), cm = divq(ONE, (q)(N - 1) * (N - 1) * (N - 1) * ONE);
            const q cont = -mulq(pi2_1440, c1 - 2 * c0 + cm);                /* per unit area */
            const q contbox = cont * ((q)1 << sh);                            /* times M^2, exact */
            const q ratio = divq(F2, contbox);                                /* whole-box numbers: full precision; cont < 0 is fine in whq.h v2 */
            hash = fnv(hash, &F2, 8); hash = fnv(hash, &ratio, 8);
            char E[48]; fmt(E1, A); fmt(F2, B); fmt(F2overM2, E); fmt(cont, C); fmt(ratio, D);
            printf("WHCASIMIR M = %d: %2d | %s | %s | %s | %s | %s\n", M, N, A, B, E, C, D);
        }
    }
    printf("WHCASIMIR: chain %016llx\n", (unsigned long long)hash);
    printf("WHCASIMIR: SI (DERIVED, not computed here): electromagnetic E_C/A = -pi^2 hbar c / (720 d^3), hbar c = 3.16153e-26 J m: d = 1 um -> -4.33e-10 J/m^2 (a pressure of 1.3 mPa); d = 10 nm -> -4.33e-4 J/m^2 (1.3 kPa); the one scalar computed here gives half of each\n");
    printf("WHCASIMIR: done\n");
    return 0;
}
