/* limb_tile.c (11-12-2026) -- C = op(A) B with the three-limb kernel, tiled as the main program's v5.1 mmul (wh_syk5.c:150-173): the output in
 * 32 x 32 tiles, the k index in blocks of 32, the accumulators of each output entry carried across the k blocks, and ONE finalize (one
 * rounding) per entry after the last k block. The sums are exact (mod 2^128), so neither the tiling, nor the order of the k blocks,
 * nor the order of k inside them, nor skipping exact zeros changes a bit: the result equals the untiled s128 product exactly
 * (limb_tile_test.c).
 *   B is split once, transposed and k-block-major (block kb holds, for every column j, its 32 terms k in [32 kb, 32 kb + 32)
 *   contiguous); op(A) is split one block of 32 rows at a time (each row once in all), k-block-major the same way. The kernel's two
 *   operand streams are unit-stride, so the D = 512 strided read of the s128 version does not arise, and the 16 column segments a
 *   half tile reuses (16 x 32 x 36 B = 18 KB) are one contiguous run: no L1 set conflicts (a column-major Bt with its 18 KB column
 *   stride would put 4 of the 16 segments on the same sets, all 4 ways of the 32 KB L1).
 *   Inside a k block the tile's columns are taken in halves of 16, so that the B sub-block the 32 rows reuse (16 x 32 entries x 36 B
 *   = 18 KB) stays in the A9's 32 KB L1 next to the current row segment (32 x 36 B) and the 16 accumulators (16 x 112 B).
 *   Zero skip by ranges (see limb_tile.h): the k order is perm (the parity order makes parity blocks contiguous); each row of op(A) and
 *   each column of B keeps [first nonzero, last nonzero + 1) in that order; a dot product runs over the intersection only. */
#include "limb_tile.h"

void lq_parity_perm(int D, int *perm)
{
    int n = 0;
    for (int p = 0; p < 2; p++) for (int k = 0; k < D; k++) if ((__builtin_popcount((unsigned)k) & 1) == p) perm[n++] = k;
}

static inline int nz(lq_c64 v) { return v.re != 0 || v.im != 0; }

long lq60_mmul_tiled(const lq_c64 *A, int mode, const lq_c64 *B, lq_c64 *C, int D, const int *perm,
                     lq60_left *Al, lq60_right *Bt, lq60_acc *acc, int *rng)
{
    const int T = D < LQ_TILE ? D : LQ_TILE;
    int *ras = rng, *rae = rng + LQ_TILE, *rbs = rng + 2 * LQ_TILE, *rbe = rng + 2 * LQ_TILE + D;
    long macs = 0;
    for (int j = 0; j < D; j++) {                                   /* B: column j, k in the chosen order, and its nonzero range */
        int s = D, e = 0;
        for (int kk = 0; kk < D; kk++) {
            const lq_c64 b = B[(long)(perm ? perm[kk] : kk) * D + j];
            lq60_split_right(b, &Bt[((long)(kk / T) * D + j) * T + kk % T]);     /* k-block-major: block kb holds every column's 32 terms */
            if (nz(b)) { if (kk < s) s = kk; e = kk + 1; }
        }
        rbs[j] = s; rbe[j] = e;
    }
    for (int i0 = 0; i0 < D; i0 += T) {
        const int ni = D - i0 < T ? D - i0 : T;
        for (int i = 0; i < ni; i++) {                              /* op(A): rows i0 .. i0+ni-1, k in the chosen order */
            int s = D, e = 0;
            for (int kk = 0; kk < D; kk++) {
                const int k = perm ? perm[kk] : kk; lq_c64 a;
                if (mode == 0) a = A[(long)(i0 + i) * D + k];
                else { a = A[(long)k * D + i0 + i]; a.im = (int64_t)((uint64_t)0 - (uint64_t)a.im); }   /* conj, as wh_syk4.c:158, wh_syk5.c:142 */
                lq60_split_left(a, &Al[((long)(kk / T) * T + i) * T + kk % T]);     /* k-block-major, as Bt */
                if (nz(a)) { if (kk < s) s = kk; e = kk + 1; }
            }
            ras[i] = s; rae[i] = e;
        }
        for (int j0 = 0; j0 < D; j0 += T) {
            const int nj = D - j0 < T ? D - j0 : T;
            lq60_zero(acc, T * T);
            for (int k0 = 0, kb = 0; k0 < D; k0 += T, kb++) {
                const int k1 = D - k0 < T ? D : k0 + T;
                for (int jh = 0; jh < nj; jh += LQ_JSUB) {
                    const int je = jh + LQ_JSUB < nj ? jh + LQ_JSUB : nj;
                    for (int i = 0; i < ni; i++) {
                        const int as = ras[i] > k0 ? ras[i] : k0, ae = rae[i] < k1 ? rae[i] : k1;
                        if (as >= ae) continue;
                        const lq60_left *ar = &Al[((long)kb * T + i) * T] - k0;
                        for (int j = jh; j < je; j++) {
                            const int lo = rbs[j0 + j] > as ? rbs[j0 + j] : as, hi = rbe[j0 + j] < ae ? rbe[j0 + j] : ae;
                            if (lo >= hi) continue;
                            const lq60_right *bc = &Bt[((long)kb * D + j0 + j) * T] - k0;
                            lq60_accumulate(&acc[i * T + j], ar + lo, bc + lo, hi - lo); macs += hi - lo;
                        }
                    }
                }
            }
            for (int i = 0; i < ni; i++) for (int j = 0; j < nj; j++) {
                lq_c64 *c = &C[(long)(i0 + i) * D + j0 + j]; lq60_finalize(&acc[i * T + j], 60, &c->re, &c->im);
            }
        }
    }
    return macs;
}
