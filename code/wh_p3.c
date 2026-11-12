/* wh_p3.c (11-12-2026, 6:35 AM) -- WORMHOLE MATH 1b, the near-brane throat (P3) of Maldacena-Milekhin's fine-tuned brane as an integer
 * ODE shot, per the reviewer's verified specification (MATH_WORK_01.md, the 4:30 AM entry of 11-12). Three rescaled variables are
 * integrated inward from the brane in the proper depth t: alpha = A - A_b, s = sigma - sigma_b and bh = (sigma' - 1)/Z with
 * Z = Q2 e^{-4 sigma}; A' comes from the constraint's + root as -c/(B + sqrt(B^2 - Z c)) (no cancellation as Z -> 0), the way wh_flow.c
 * takes Lprime from its constraint. RK4 with h = 2^-12 and the once-rounded update with the carried remainder (wh_flow.c v3's rule), in
 * Q.56 (values below 128: the shot ell runs 5.3 to 14.5 and the RK4 sum reaches 83 at MM's charge). One shooting parameter, ell = S^2/Q2,
 * which is exactly ln(L_IR/R5) in the paper's eq. (3.24), bisected on the FIRST EVENT of each trial (the specification's rule E): UP
 * (sigma' < 0: ell too large); DOWN (s < s*), TPA (A' <= 0) or TP (the radicand negative at a stage): ell too small. Exact integer tests,
 * no tolerance, no fitted depth. The arithmetic is whq.h's at 56 bits (whq56.h: mulq and divq round to nearest, sqrtq the floor root,
 * expq, lnq) in the specification's order of operations, so the result should match the reviewer's integer mimic bit for bit (its raw
 * ells: 383814628494007828 at Q2 = 6 2^20, 492682328441613039 at 6 2^28, 1043449860927232495 at MM's charge) or lie within its stated
 * tolerance for another unbiased arithmetic (8 units of 2^-56 at the reduced ratios, 40 at MM). The first integral W = e^{2 alpha - 2 s}
 * (ah - bh) is carried against its integrated right side P as a check (W - P printed on the profile rows). Inputs per case (section 6
 * of the specification): LQ = ln Q2, ln Z* (Z* = 6 + 3w/2, the fixed point), the bracket [ell_LO, ell_LO + w]. kappa = 1 (a round S2),
 * L = R5 = 1. No floating point, no libm: the same bits on every machine.
 * Usage: wh_p3 [case=20|28|mm] [steps=200000] [profile=1] [every=2048]       (every: the profile's row spacing in steps; 2048 = 0.5 in t) */
#include "whq56.h"
#include <stdlib.h>

static const q TWO3 = 48038396025285291LL;         /* 2/3 in Q.56 */
static q LQ, LNZS, LO, HI, EXPECT;                 /* the case: ln Q2, ln Z*, the bracket, the mimic's raw */
static q LZ, LM, KL, SSTAR;                        /* per trial: LZ = LQ + 2 ln ell, ell - 1/2, kappa ell, s* = -(ln Z* + LZ)/4 */
static int STEPS_MAX = 200000, EVERY = 2048;
static uint64_t CHAIN = 0xcbf29ce484222325ull;
enum { EV_NONE = 0, EV_UP, EV_DOWN, EV_TPA, EV_TP };
static const char *EVNAME[] = { "none", "UP (sigma' < 0: ell too large)", "DOWN (s < s*: ell too small)", "TPA (A' <= 0: ell too small)", "TP (the radicand negative at a stage: ell too small)" };

typedef struct { q al, s, bh, P; } st;             /* alpha, s, bh, and P (the first integral's right side, integrated) */
typedef struct { q al, s, bh, P, Z, Zb, Za, ah; int tp; } rhs_t;

static void fev(const st *y, rhs_t *k)             /* the right sides, in the specification's order of operations */
{
    const q Z = expq(-4 * y->s - LZ);                                 /* Q2 e^{-4 sigma} */
    const q XZ = mulq(LM, expq(4 * y->s - 2 * y->al));               /* e^{-2A}/Z */
    const q YZ = mulq(KL, expq(2 * y->s));                            /* kappa e^{-2 sigma}/Z */
    const q Zb = mulq(Z, y->bh);
    const q c = 6 * y->bh + mulq(Zb, y->bh) + XZ - YZ + ONE / 2;
    const q B = 3 * ONE + 2 * Zb;
    const q rad = mulq(B, B) - mulq(Z, c);
    k->Z = Z; k->Zb = Zb; k->tp = rad < 0;
    if (rad < 0) { k->al = k->s = k->bh = k->P = k->Za = k->ah = 0; return; }
    const q sq = sqrtq(rad);                                           /* = A' + 2 sigma' exactly */
    const q ah = divq(-c, B + sq), Za = mulq(Z, ah);                   /* (A' - 1)/Z and A' - 1 */
    k->ah = ah; k->Za = Za;
    k->al = -(ONE + Za);                                               /* d alpha/dt = -A' */
    k->s = -(ONE + Zb);                                                /* d s/dt = -sigma' */
    k->bh = 2 * ah + 2 * y->bh + TWO3 - YZ + 2 * mulq(Za, y->bh) - 2 * mulq(Zb, y->bh);
    k->P = -(expq(2 * y->al - 2 * y->s) - mulq(KL, expq(2 * y->al)) - mulq(LM, expq(2 * y->s)));   /* dP/dt */
}
static q upd(q X, q *R)                            /* the increment h X/6 for h = 2^-12, rounded once (half up), the remainder carried: v3's rule */
{
    const q D = 6 << 12; const q T = X + *R + D / 2; q inc = T / D, rem = T % D; if (rem < 0) { inc -= 1; rem += D; }
    *R = rem - D / 2; return inc;
}
static void fmt56_18(q v, char *out)              /* eighteen decimals of a Q.56 value (truncated) */
{
    const uint64_t u = v < 0 ? -(uint64_t)v : (uint64_t)v, ip = u >> 56, fp = u & MASK; uint64_t h, l;
    mul128(fp, 1000000000ull, &h, &l); const uint64_t d1 = (h << 8) | (l >> 56), r1 = l & MASK;
    mul128(r1, 1000000000ull, &h, &l); const uint64_t d2 = (h << 8) | (l >> 56);
    sprintf(out, "%s%llu.%09llu%09llu", v < 0 ? "-" : "", (unsigned long long)ip, (unsigned long long)d1, (unsigned long long)d2);
}
static void stage(const st *y, const rhs_t *k, q hf, st *o)   /* o = y + hf k (hf the mulq factor for h/2 or h) */
{ o->al = y->al + mulq(hf, k->al); o->s = y->s + mulq(hf, k->s); o->bh = y->bh + mulq(hf, k->bh); o->P = y->P + mulq(hf, k->P); }

/* one trial at ell: integrates to the first event; returns the event and the step count; with profile, prints rows every EVERY steps */
static int trial(q ell, int *steps, int profile, q sigb, q Ab)
{
    LZ = LQ + 2 * lnq(ell); LM = ell - ONE / 2; KL = ell; SSTAR = -rdiv(LNZS + LZ, 4);
    st y = {0, 0, 0, 0}; q R[4] = {0, 0, 0, 0}; rhs_t k1, k2, k3, k4; st y2, y3, y4;
    const q H2 = (q)1 << 43, H1 = (q)1 << 44;                        /* h/2 and h as Q.56 factors */
    fev(&y, &k1); if (k1.tp) { *steps = 0; return EV_TP; }
    for (int n = 0; n < STEPS_MAX; n++) {
        if (profile && n % EVERY == 0) {
            char a[40], b[40], c[40], d[40], e[40], f[40], g[40]; const q W = mulq(expq(2 * y.al - 2 * y.s), k1.ah - y.bh);
            fmt56_18(Ab + y.al, a); fmt56_18(sigb + y.s, b); fmt56_18(ONE + k1.Za, c); fmt56_18(ONE + k1.Zb, d); fmt56_18(y.bh, e); fmt56_18(k1.ah, f); fmt56_18(W - y.P, g);
            printf("WHP3 row t = %d/4096 = %s: A = %s, sigma = %s, A' = %s, sigma' = %s, bh = %s, ah = %s, W - P = %s\n", n, (n % 4096 == 0) ? (n / 4096 == 0 ? "0" : "int") : "half", a, b, c, d, e, f, g);
            const q raws[6] = { y.al, y.s, y.bh, k1.Za, k1.Zb, k1.ah }; CHAIN = fnv(CHAIN, raws, sizeof raws);
        }
        if (ONE + k1.Zb < 0) { *steps = n; return EV_UP; }           /* the events, tested at the state after each step */
        if (y.s < SSTAR) { *steps = n; return EV_DOWN; }
        if (ONE + k1.Za <= 0) { *steps = n; return EV_TPA; }
        stage(&y, &k1, H2, &y2); fev(&y2, &k2); if (k2.tp) { *steps = n; return EV_TP; }
        stage(&y, &k2, H2, &y3); fev(&y3, &k3); if (k3.tp) { *steps = n; return EV_TP; }
        stage(&y, &k3, H1, &y4); fev(&y4, &k4); if (k4.tp) { *steps = n; return EV_TP; }
        y.al += upd(k1.al + 2 * k2.al + 2 * k3.al + k4.al, &R[0]);
        y.s  += upd(k1.s  + 2 * k2.s  + 2 * k3.s  + k4.s,  &R[1]);
        y.bh += upd(k1.bh + 2 * k2.bh + 2 * k3.bh + k4.bh, &R[2]);
        y.P  += upd(k1.P  + 2 * k2.P  + 2 * k3.P  + k4.P,  &R[3]);
        fev(&y, &k1); if (k1.tp) { *steps = n + 1; return EV_TP; }
    }
    *steps = STEPS_MAX; return EV_NONE;
}
int main(int argc, char **argv)
{
    const char *cs = "20"; int profile = 1;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "case=", 5)) cs = argv[i] + 5;
        else if (!strncmp(argv[i], "steps=", 6)) STEPS_MAX = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "profile=", 8)) profile = atoi(argv[i] + 8);
        else if (!strncmp(argv[i], "every=", 6)) EVERY = atoi(argv[i] + 6);
    }
    const char *cname;
    if (!strcmp(cs, "20")) { cname = "Q2 = 6 2^20 = 6291456 (lhat = 1/32, w = 9.77e-4)"; LQ = 1128040239353705893LL; LNZS = 129127468633249144LL; LO = 383788833474402286LL; HI = 383859210809038831LL; EXPECT = 383814628494007828LL; }
    else if (!strcmp(cs, "28")) { cname = "Q2 = 6 2^28 = 1610612736 (lhat = 1/128, w = 6.10e-5)"; LQ = 1527612384516288882LL; LNZS = 129110975958876184LL; LO = 492680742577044628LL; HI = 492685140657110293LL; EXPECT = 492682328441613039LL; }
    else if (!strcmp(cs, "mm")) { cname = "MM's charge, Q2 = 6.8112e21 (re = 1.5e7 m by the paper's chain, R5 = 5e-5 m; w = 2.97e-11)"; LQ = 3622541031543053145LL; LNZS = 129109876447783086LL; LO = 1043449860854426327LL; HI = 1043449860998541516LL; EXPECT = 1043449860927232495LL; }
    else { printf("WHP3: case must be 20, 28 or mm\n"); return 1; }
    char a[40], b[40], c[40], d[40];
    fmt56_18(LQ, a); fmt56_18(LNZS, b); fmt56_18(LO, c); fmt56_18(HI, d);
    printf("WHP3: the near-brane throat (P3) as an integer shot in Q.56 (unit 2^-56 = 1.39e-17): alpha, s, bh = (sigma' - 1)/Z inward from the brane, A' from the constraint's + root, RK4 h = 2^-12 with the carried remainder, ell = S^2/Q2 = ln(L_IR/R5) bisected on the first event (rule E); arithmetic whq.h's at 56 bits\n");
    printf("WHP3 case %s: %s; ln Q2 = %s (raw %lld), ln Z* = %s (raw %lld), the bracket [%s, %s] (raws %lld, %lld), the mimic's raw ell %lld\n", cs, cname, a, (long long)LQ, b, (long long)LNZS, c, d, (long long)LO, (long long)HI, (long long)EXPECT);
    int st0, st1; const int e0 = trial(LO, &st0, 0, 0, 0), e1 = trial(HI, &st1, 0, 0, 0);
    printf("WHP3 bracket ends: ell_LO -> %s at t = %d/4096; ell_LO + w -> %s at t = %d/4096\n", EVNAME[e0], st0, EVNAME[e1], st1);
    if (e1 != EV_UP || e0 == EV_UP || e0 == EV_NONE) { printf("WHP3: the bracket does not have DOWN at the bottom and UP at the top -- stop\n"); return 1; }
    q lo = LO, hi = HI; int halvings = 0, lev = 0, lst = 0;
    while (hi - lo > 1) {
        const q mid = lo + (hi - lo) / 2; int stp; const int ev = trial(mid, &stp, 0, 0, 0);
        if (ev == EV_NONE) { printf("WHP3: a trial reached %d steps without an event -- stop\n", STEPS_MAX); return 1; }
        if (ev == EV_UP) hi = mid; else lo = mid;
        halvings++; lev = ev; lst = stp;
    }
    const q ell = lo; CHAIN = fnv(CHAIN, &ell, 8);
    fmt56_18(ell, a);
    printf("WHP3 ell = %s (raw %lld; the lower end of the final 1-unit bracket after %d halvings; the mimic's raw %lld: %lld units apart); the last trial ended %s at t = %d/4096 = %d.%04d\n",
           a, (long long)ell, halvings, (long long)EXPECT, (long long)(ell - EXPECT), EVNAME[lev], lst, lst / 4096, (lst % 4096) * 10000 / 4096);
    /* the outputs of section 5 */
    const q lnell = lnq(ell), ln6 = lnq(6 * ONE);
    const q kappa = ell - rdiv(ln6 + LQ + 2 * lnell, 4);                            /* kappa_P3 = ell - ln(l_B/R5), ln(l_B/R5) = (1/4) ln(6 Q2 ell^2) */
    const q sigb = rdiv(LQ + lnell, 2), one_m = ONE - divq(ONE, 2 * ell);           /* sigma_b = (ln Q2 + ln ell)/2;  1 - 1/(2 ell) */
    const q Ab = sigb - rdiv(lnq(one_m), 2);                                         /* A_b = sigma_b - (1/2) ln(1 - 1/(2 ell)) */
    const q ras = divq(ONE, sqrtq(one_m)) - ONE, invell = divq(ONE, ell);            /* R_A/S - 1 and 1/ell = R4 S^2 */
    const q lnS = rdiv(LQ + lnell, 2), lbrat = sqrtq(sqrtq(divq(2 * ell, lnS)));      /* (2 ell / ln S)^(1/4) = l_B / the paper's (3.25) l_B */
    fmt56_18(kappa, a); fmt56_18(ras, b); fmt56_18(invell, c); fmt56_18(lbrat, d);
    printf("WHP3 kappa_P3 = ell - ln(l_B/R5) = %s (expected kappa_flat + c w: 0.12853534371656 at 6 2^20, 0.12823136939461 at 6 2^28, 0.12821097053423 at MM); R_A/S - 1 = %s; R4 S^2 = 1/ell = %s; l_B against the paper's (3.25): %s\n", a, b, c, d);
    fmt56_18(sigb, a); fmt56_18(Ab, b); fmt56_18(SSTAR, c);
    printf("WHP3 brane: sigma_b = %s, A_b = %s; s* = %s (of the last trial)\n", a, b, c);
    if (profile) { int stp; const int ev = trial(ell, &stp, 1, sigb, Ab); printf("WHP3 profile at ell: the first event %s at t = %d/4096\n", EVNAME[ev], stp); }
    printf("WHP3: chain %016llx\n", (unsigned long long)CHAIN);
    printf("WHP3: done\n");
    return 0;
}
