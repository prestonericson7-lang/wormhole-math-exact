/* limb_q60.c (11-12-2026) -- the three-limb kernel with the three-multiplication complex product; contract and layout in limb_q60.h.
 * Portable C99 (gnu11). The hot loop is lq60_accumulate: four passes over the k block, so that no pass needs more of the A9's 14
 * usable registers than it has: (1)-(3) classes 0, 1, 2 and 4 of P1 = SUM ar br, P2 = SUM ai bi, P3 = SUM (ar+ai)(br+bi) (seven
 * registers of sums, the three limbs of y and one of x live), (4) class 3 of all three. The sums are signed int64 (see LQ_P in
 * limb_common.h), class 4 a uint32; the contract keeps every signed sum inside int64.
 * Measured on the A9 (qemu instruction counts, count_a9.out): 76.7 instructions a complex MAC at 64 k terms a call, 79.2 at 32;
 * the four-multiplication schedule (limb_m4.c, limb_variants.c v_a4s) needs 86.0 / 87.9, the s128 kernel 235.6 (dot form). */
#include "limb_q60.h"

#define M24 0xffffff

static inline void split3(int64_t v, int32_t *h, int32_t *m, int32_t *l)
{
    *l = (int32_t)(v & M24);
    *m = (int32_t)((v >> 24) & M24);
    *h = (int32_t)(v >> 48);                      /* floor(v / 2^48) in [-2^15, 2^15 - 1] */
}
/* the limbs of u + v (a 65-bit integer) without forming it: low and middle limbs normalized to [0, 2^24), top limb in [-2^16, 2^16) */
static inline void split_sum(int64_t u, int64_t v, int32_t *h, int32_t *m, int32_t *l)
{
    const int64_t s0 = (u & M24) + (v & M24);                                   /* [0, 2^25 - 2] */
    const int64_t s1 = ((u >> 24) & M24) + ((v >> 24) & M24) + (s0 >> 24);     /* [0, 2^25 - 1] */
    *l = (int32_t)(s0 & M24); *m = (int32_t)(s1 & M24);
    *h = (int32_t)((u >> 48) + (v >> 48) + (s1 >> 24));
}

void lq60_split_left(lq_c64 a, lq60_left *o)
{
    split3(a.re, &o->rh, &o->rm, &o->rl); split3(a.im, &o->ih, &o->im, &o->il); split_sum(a.re, a.im, &o->sh, &o->sm, &o->sl);
}

void lq60_split_right(lq_c64 b, lq60_right *o)
{
    split3(b.re, &o->rh, &o->rm, &o->rl); split3(b.im, &o->ih, &o->im, &o->il); split_sum(b.re, b.im, &o->sh, &o->sm, &o->sl);
}

int lq60_split_left_rows(const lq_c64 *A, int lda, int ni, int nk, lq60_left *out)
{
    int bmax = 0;
    for (int i = 0; i < ni; i++) for (int k = 0; k < nk; k++) {
        const lq_c64 a = A[(long)i * lda + k]; const int b1 = lq_bits(a.re), b2 = lq_bits(a.im);
        if (b1 > bmax) bmax = b1; if (b2 > bmax) bmax = b2;
        lq60_split_left(a, &out[(long)i * nk + k]);
    }
    return bmax;
}

int lq60_split_right_cols(const lq_c64 *B, int ldb, int nk, int nj, lq60_right *out)
{
    int bmax = 0;
    for (int k = 0; k < nk; k++) for (int j = 0; j < nj; j++) {
        const lq_c64 b = B[(long)k * ldb + j]; const int b1 = lq_bits(b.re), b2 = lq_bits(b.im);
        if (b1 > bmax) bmax = b1; if (b2 > bmax) bmax = b2;
        lq60_split_right(b, &out[(long)j * nk + k]);
    }
    return bmax;
}

void lq60_zero(lq60_acc *acc, int n)
{
    for (int i = 0; i < n; i++) {
        for (int c = 0; c < 4; c++) acc[i].p1[c] = acc[i].p2[c] = acc[i].p3[c] = 0;
        acc[i].e1 = acc[i].e2 = acc[i].e3 = acc[i].pad = 0;
    }
}

void lq60_accumulate(lq60_acc *acc, const lq60_left *a, const lq60_right *b, int n)
{
    {   /* pass 1: P1 = SUM ar br, classes 0, 1, 2 and 4 */
        int64_t c0 = acc->p1[0], c1 = acc->p1[1], c2 = acc->p1[2]; uint32_t e4 = acc->e1;
        for (int k = 0; k < n; k++) {
            const lq60_left *x = a + k; const lq60_right *y = b + k; int32_t t;
            const int32_t yl = y->rl, ym = y->rm, yh = y->rh;
            t = x->rl; c0 += LQ_P(t, yl); c1 += LQ_P(t, ym); c2 += LQ_P(t, yh);
            t = x->rm; c1 += LQ_P(t, yl); c2 += LQ_P(t, ym);
            t = x->rh; c2 += LQ_P(t, yl); e4 += (uint32_t)t * (uint32_t)yh;
        }
        acc->p1[0] = c0; acc->p1[1] = c1; acc->p1[2] = c2; acc->e1 = e4;
    }
    {   /* pass 2: P2 = SUM ai bi */
        int64_t c0 = acc->p2[0], c1 = acc->p2[1], c2 = acc->p2[2]; uint32_t e4 = acc->e2;
        for (int k = 0; k < n; k++) {
            const lq60_left *x = a + k; const lq60_right *y = b + k; int32_t t;
            const int32_t yl = y->il, ym = y->im, yh = y->ih;
            t = x->il; c0 += LQ_P(t, yl); c1 += LQ_P(t, ym); c2 += LQ_P(t, yh);
            t = x->im; c1 += LQ_P(t, yl); c2 += LQ_P(t, ym);
            t = x->ih; c2 += LQ_P(t, yl); e4 += (uint32_t)t * (uint32_t)yh;
        }
        acc->p2[0] = c0; acc->p2[1] = c1; acc->p2[2] = c2; acc->e2 = e4;
    }
    {   /* pass 3: P3 = SUM (ar + ai)(br + bi) */
        int64_t c0 = acc->p3[0], c1 = acc->p3[1], c2 = acc->p3[2]; uint32_t e4 = acc->e3;
        for (int k = 0; k < n; k++) {
            const lq60_left *x = a + k; const lq60_right *y = b + k; int32_t t;
            const int32_t yl = y->sl, ym = y->sm, yh = y->sh;
            t = x->sl; c0 += LQ_P(t, yl); c1 += LQ_P(t, ym); c2 += LQ_P(t, yh);
            t = x->sm; c1 += LQ_P(t, yl); c2 += LQ_P(t, ym);
            t = x->sh; c2 += LQ_P(t, yl); e4 += (uint32_t)t * (uint32_t)yh;
        }
        acc->p3[0] = c0; acc->p3[1] = c1; acc->p3[2] = c2; acc->e3 = e4;
    }
    {   /* pass 4: class 3 of P1, P2, P3 */
        int64_t q1 = acc->p1[3], q2 = acc->p2[3], q3 = acc->p3[3];
        for (int k = 0; k < n; k++) {
            const lq60_left *x = a + k; const lq60_right *y = b + k; int32_t t, ym, yh;
            ym = y->rm; yh = y->rh; t = x->rm; q1 += LQ_P(t, yh); t = x->rh; q1 += LQ_P(t, ym);
            ym = y->im; yh = y->ih; t = x->im; q2 += LQ_P(t, yh); t = x->ih; q2 += LQ_P(t, ym);
            ym = y->sm; yh = y->sh; t = x->sm; q3 += LQ_P(t, yh); t = x->sh; q3 += LQ_P(t, ym);
        }
        acc->p1[3] = q1; acc->p2[3] = q2; acc->p3[3] = q3;
    }
}

static inline uint64_t sar64(uint64_t v, int s) { return (v >> s) | (((uint64_t)0 - (v >> 63)) << (64 - s)); }   /* arithmetic >>, 1..63 */
static inline void add128(uint64_t *hi, uint64_t *lo, uint64_t h, uint64_t l) { const uint64_t t = *lo + l; *hi += h + (t < l); *lo = t; }

/* T = C0 + C1 2^24 + C2 2^48 + C3 2^72 + C4 2^96 (mod 2^128): C0, C1, C2 signed; C3 and C4 enter only as residues */
static void sum5(const int64_t *c, uint32_t c4, uint64_t *hi, uint64_t *lo)
{
    const uint64_t a0 = (uint64_t)c[0], a1 = (uint64_t)c[1], a2 = (uint64_t)c[2], a3 = (uint64_t)c[3];
    *hi = sar64(a0, 63); *lo = a0;
    add128(hi, lo, sar64(a1, 40), a1 << 24);
    add128(hi, lo, sar64(a2, 16), a2 << 48);
    *hi += (a3 << 8) + ((uint64_t)c4 << 32);
}
static void sub128(uint64_t *hi, uint64_t *lo, uint64_t h, uint64_t l) { add128(hi, lo, ~h + (l == 0), ~l + 1); }   /* (hi:lo) -= (h:l) */

lq_s128 lq60_sum_re(const lq60_acc *acc)
{
    uint64_t h1, l1, h2, l2; sum5(acc->p1, acc->e1, &h1, &l1); sum5(acc->p2, acc->e2, &h2, &l2);
    sub128(&h1, &l1, h2, l2);
    lq_s128 s; s.hi = (int64_t)h1; s.lo = l1; return s;
}
lq_s128 lq60_sum_im(const lq60_acc *acc)
{
    uint64_t h1, l1, h2, l2, h3, l3; sum5(acc->p1, acc->e1, &h1, &l1); sum5(acc->p2, acc->e2, &h2, &l2); sum5(acc->p3, acc->e3, &h3, &l3);
    sub128(&h3, &l3, h1, l1); sub128(&h3, &l3, h2, l2);
    lq_s128 s; s.hi = (int64_t)h3; s.lo = l3; return s;
}

void lq60_finalize(const lq60_acc *acc, int sh, int64_t *re, int64_t *im)      /* = lq_round_sh(lq60_sum_re / _im (acc), sh), each T_P formed once */
{
    uint64_t h1, l1, h2, l2, h3, l3; sum5(acc->p1, acc->e1, &h1, &l1); sum5(acc->p2, acc->e2, &h2, &l2); sum5(acc->p3, acc->e3, &h3, &l3);
    sub128(&h3, &l3, h1, l1); sub128(&h3, &l3, h2, l2); sub128(&h1, &l1, h2, l2);
    lq_s128 r, i; r.hi = (int64_t)h1; r.lo = l1; i.hi = (int64_t)h3; i.lo = l3;
    *re = lq_round_sh(r, sh); *im = lq_round_sh(i, sh);
}

/* the contract's largest L (limb_q60.h; bounds.py): P1, P2 (operand bounds Ba, Bb) and P3 (sums: top limbs twice as large), L products
   each: C1 <= 2 L M^2 <= 2^63 - 1 binds (L <= 16384); C0 <= L M^2; C2 of P3 from above L (M^2 + (2Ha + 2Hb - 2) M) <= 2^63 - 1; C2 and C3
   of P3 from below L (2Ha + 2Hb) M <= 2^63 */
long lq60_lmax(int Ba, int Bb)
{
    if (Ba < 48) Ba = 48; if (Bb < 48) Bb = 48;
    if (Ba > 63 || Bb > 63) return 0;
    const uint64_t Ha = (uint64_t)1 << (Ba - 48), Hb = (uint64_t)1 << (Bb - 48), M = M24;
    uint64_t L = 16384;
    const uint64_t l2 = (((uint64_t)1 << 63) - 1) / (M * M + (2 * Ha + 2 * Hb - 2) * M);
    const uint64_t l3 = ((uint64_t)1 << 63) / ((2 * Ha + 2 * Hb) * M);
    if (l2 < L) L = l2; if (l3 < L) L = l3;
    return (long)L;
}
