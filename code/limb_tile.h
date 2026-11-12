/* limb_tile.h (11-12-2026) -- the tiled product with the three-limb kernel (an example of the integration; see limb_tile.c) */
#ifndef LIMB_TILE_H
#define LIMB_TILE_H
#include "limb_q60.h"
#define LQ_TILE 32          /* output tiles LQ_TILE x LQ_TILE, k in blocks of LQ_TILE */
#define LQ_JSUB 16          /* inside a k block the tile's columns go in two halves: the B sub-block 16 x 32 x 36 B = 18 KB stays in L1 */
/* the k order that makes parity blocks contiguous: perm[k'] = the k of rank k' in the order (popcount parity of k, k) */
void lq_parity_perm(int D, int *perm);
/* C = op(A) B (D x D, row-major, Q.60 in and out), op(A) = A (mode 0) or A^dagger (mode 1).
   perm: the order in which the k index is stored and summed (NULL = 0, 1, ..., D-1); any order gives the same exact sums. Exact zeros
   are skipped the way v5's zero skip does, but by RANGES: after the reordering, each row of op(A) and each column of B records the
   range [first, last + 1) of its nonzero k positions, and a dot product runs only over the intersection. With the parity order, a
   parity-block-diagonal op(A) costs D^3/2 MACs and two block-diagonal operands D^3/4.
   Work buffers (the caller's): Al LQ_TILE*D' entries (576 KB at D = 512), Bt D*D' entries (9 MB at D = 512), with D' = D rounded up to a
   multiple of LQ_TILE (D' = D for every power of two from 32 up; below 32 the tile is D itself), acc LQ_TILE*LQ_TILE (112 KB),
   rng 2*LQ_TILE + 2*D ints. Returns the number of complex MACs done. */
long lq60_mmul_tiled(const lq_c64 *A, int mode, const lq_c64 *B, lq_c64 *C, int D, const int *perm,
                     lq60_left *Al, lq60_right *Bt, lq60_acc *acc, int *rng);
#endif
