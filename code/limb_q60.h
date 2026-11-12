/* limb_q60.h (11-12-2026) -- THE THREE-LIMB KERNEL: exact complex dot products of Q.60 operands up to the Q.60 ceiling of 8 (any
 * int64 except INT64_MIN), by 24-bit limbs, with the three-multiplication complex product; int32 x int32 -> int64 products, int64
 * and uint32 accumulators only (no __int128).
 *
 * THE ARITHMETIC. For one output, over the k terms a_k = ar + i ai (row of op(A)) and b_k = br + i bi (column of B):
 *     P1 = SUM ar br,   P2 = SUM ai bi,   P3 = SUM (ar + ai)(br + bi);     re = P1 - P2,   im = P3 - P1 - P2
 *   (exact identities over the integers: P3 - P1 - P2 = SUM (ar bi + ai br)), each P a real dot product formed EXACTLY modulo 2^128
 *   from 24-bit limbs, and re, im rounded ONCE. 27 limb products a complex MAC (9 a real product: 24 SMLAL + 3 MLA) instead of 36.
 * SPLIT (once per matrix): v = h 2^48 + m 2^24 + l,  l = v mod 2^24, m = floor(v / 2^24) mod 2^24 (both in [0, 2^24)),
 *   h = floor(v / 2^48): in [-2^15, 2^15 - 1] for an int64 v, in [-2^16, 2^16 - 1] for the 65-bit sums s = re + im (split exactly,
 *   with carries, so s needs no 65-bit integer).
 *   left operand op(A)[i][k] -> lq60_left (36 bytes): limbs of re, im, s = re + im. Row i of op(A) contiguous in k.
 *   right operand B[k][j] -> lq60_right (36 bytes): limbs of re, im, s = re + im, stored TRANSPOSED: column j of B contiguous in k.
 * ACCUMULATE: acc(i,j) += SUM_k over a block of k terms. For each P (x y = ar br, ai bi, s_a s_b) five sums by weight 2^(24c):
 *     C0 = SUM x_l y_l,  C1 = SUM (x_l y_m + x_m y_l),  C2 = SUM (x_l y_h + x_m y_m + x_h y_l),  C3 = SUM (x_m y_h + x_h y_m),  C4 = SUM x_h y_h,
 *   C0..C3 in signed int64, C4 modulo 2^32 (uint32: weight 2^96, so only C4 mod 2^32 reaches a 128-bit sum). Four passes over the k
 *   block: P1 classes 0-2 and 4; P2 the same; P3 the same; classes 3 of P1, P2, P3. At most seven registers of sums live in each.
 * FINALIZE: T_P = C0 + C1 2^24 + C2 2^48 + C3 2^72 + C4 2^96 modulo 2^128 for each P; re = T_P1 - T_P2, im = T_P3 - T_P1 - T_P2 modulo
 *   2^128; each rounded ONCE with the main program's rule: lq_round_sh(., sh), sh = 60 for Q.60 x Q.60 -> Q.60 (any 1..63).
 *
 * CONTRACT (proof in the note): left entries |re|, |im| < 2^Ba, right entries |re|, |im| < 2^Bb, Ba, Bb <= 63 (no INT64_MIN), L = the
 * number of k terms between lq60_zero and lq60_finalize. C0..C3 of all three P stay inside int64 for EVERY such input if and only if
 *     2 L (2^24 - 1)^2 <= 2^63 - 1  (L <= 16384)   and two weaker conditions on C2 and C3 of P3 (L <= 32513 and L <= 4194304 at
 *     Ba = Bb = 63), i.e.  L <= lq60_lmax(Ba, Bb) = 16384 for every operand size.
 * Under it every T_P is exact modulo 2^128, so re and im are the s128 kernel's 128-bit sums exactly, and the result is BIT-EQUAL to
 * the main program's s128 kernel (wh_syk5.c:61-75 and mmul 155-173) for every input, including results too large for an int64 (same bits).
 *   ANY Q.60 operands (|raw| <= 2^63 - 1, values below 8), any L <= 16384 k terms per finalize: bit-equal. (D = 512 needs L = 512.)
 * INT64_MIN is excluded because the s128 kernel's own -a.im (wh_syk5.c:167) and conj (142, 147) negate it, which C leaves undefined.
 */
#ifndef LIMB_Q60_H
#define LIMB_Q60_H
#include "limb_common.h"

typedef struct { int32_t rl, rm, rh, il, im, ih, sl, sm, sh; } lq60_left;          /* 36 bytes an entry */
typedef struct { int32_t rl, rm, rh, il, im, ih, sl, sm, sh; } lq60_right;         /* 36 bytes an entry */
typedef struct { int64_t p1[4], p2[4], p3[4]; uint32_t e1, e2, e3, pad; } lq60_acc;  /* 112 bytes an output: C0..C3 of P1, P2, P3; C4 */

#define LQ60_LMAX 16384L              /* lq60_lmax(Ba, Bb) for every Ba, Bb <= 63 */

void lq60_split_left(lq_c64 a, lq60_left *o);
void lq60_split_right(lq_c64 b, lq60_right *o);
/* whole-matrix helpers (op = identity): out[i*nk + k] = split(A[i*lda + k]); out[j*nk + k] = split(B[k*ldb + j]); both return the
   largest lq_bits over the entries' re and im (64 means an INT64_MIN entry: outside the contract) */
int lq60_split_left_rows(const lq_c64 *A, int lda, int ni, int nk, lq60_left *out);
int lq60_split_right_cols(const lq_c64 *B, int ldb, int nk, int nj, lq60_right *out);
void lq60_zero(lq60_acc *acc, int n);
void lq60_accumulate(lq60_acc *acc, const lq60_left *a, const lq60_right *b, int n);
lq_s128 lq60_sum_re(const lq60_acc *acc);          /* the exact sums modulo 2^128 (no rounding) */
lq_s128 lq60_sum_im(const lq60_acc *acc);
void lq60_finalize(const lq60_acc *acc, int sh, int64_t *re, int64_t *im);
long lq60_lmax(int Ba, int Bb);                    /* the contract's largest L (0 if an operand bound exceeds 63) */

#endif
