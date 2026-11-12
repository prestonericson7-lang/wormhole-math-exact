/* wh_p1data.c (11-12-2026, 7:45 AM) -- WORMHOLE MATH 1b, P1's boundary data in the cylinder frame, as the reviewer's verified data entry
 * (2:41 AM, 11-12) prescribes, in Q.60 fixed point with the circle coordinate in Q.56: the four-dimensional boundary metric of Maldacena-
 * Milekhin's SO(3) example (the throat (2.7), the mouths (2.4), the S3 exterior) as one conformal class, written as g_cyl = -e^{-2 Phi(eta)}
 * dt^2 + d eta^2 + dOmega^2 with Phi = ln rho_o, its three pieces elementary (Phi_T = ln ell - ln cosh(eta + ln 2 gamma), Phi_M = ln(4 re) +
 * 2 ln cosh(eta/2), Phi_S = ln a - ln cosh(eta - ln(2a/re))), blended in ln rho_o by the C-infinity step St(u) = 1/(1 + exp(1/u - 1/(1 - u)))
 * over windows of half-width 1 centred at eta_c = -(2/3) ln gamma and eta_d = (2/3) ln(a/re). From Phi: the local conformal invariant
 * w = 2 Phi'' - 2 Phi'^2 + 2 (0 in the throat and the exterior, 3 sech^2(eta/2) in the mouth), the Fefferman-Graham coefficient g2 = minus the
 * Schouten tensor (mixed diag((3 - w)/6, (3 - w)/6, (w - 6)/12, (w - 6)/12)), the flux log term (3/(2 lhat^4)) z^4 ln(z/lhat) diag(1, 1, -1, -1),
 * and the Dirichlet data at the cutoff z = eps = lhat/30: ln T = -ln eps - Phi + (1/2) ln(1 + eps^2 (3 - w)/6 + f), ln S = -ln eps + (1/2) ln(1 +
 * eps^2 (w - 6)/12 - f), ln h_etaeta = -2 ln eps + ln(1 + eps^2 (3 - w)/6 + f), f = (3/(2 lhat^4)) eps^4 ln(eps/lhat), h_zz = 1/eps^2, h_eta z = 0.
 * Every piece that passes 8 is formed in 128 bits or in Q.56 (eta, Phi) or Q.40 (the step's 1/u^3 factors); no floating point.
 * Checks (the entry's table (e), gamma = 1e3, a = ell, re = 1): Phi at the throat window edge (eta = -5.6052) = ln 266.898, at the throat-mouth
 * centre (-4.6052) = ln 100.874 with w = 0.206, at the mouth window edge (-3.6052) = ln 38.815 with w = 0.309, at the photon sphere (0) = ln 4
 * with w = 3; the throat centre Phi = ln 1000, w = 0.
 * Usage: wh_p1data [gamma=1000] [a=1000] [lhat=0.25] [Ge=76] [every=4]      (a in units of re; Ge intervals over the quarter from the throat
 * centre eta = -ln 2 gamma to the photon sphere eta = 0; the other quarter is the mirror when a = ell) */
#include "whq60.h"
#include <string.h>
#include <stdlib.h>

#define ONE56 ((int64_t)1 << 56)
#define ONE40 ((int64_t)1 << 40)
static inline void split60(s128 x, int64_t *n, int64_t *f) { *n = (int64_t)(((uint64_t)x.hi << 4) | (x.lo >> 60)); *f = (int64_t)(x.lo & (((uint64_t)1 << 60) - 1)); }
static inline int64_t exp_of128(s128 x) { int64_t n, f; split60(x, &n, &f); return n < -62 ? 0 : exp60n(n, f); }
static inline s128 s128_neg(s128 a) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); return a; }
static s128 ln_s128(int64_t x) { int k = 0; while (x < ONE60 / 2) { x <<= 1; k++; } s128 r = {0, 0}; s128_addv(&r, ln60(x)); s128_addp(&r, -k, ln2_60()); return r; }   /* ln x for any positive Q.60 x */
static int64_t ln56_of_int(int64_t m) { s128 r = {0, 0}; int k = 0; while (m >= 2) { if (m & 1) break; m >>= 1; k++; }   /* ln of a positive integer in Q.56: m = 2^k odd... general: */
    /* general positive integer: m = 2^k y with y in [1, 2) by shifting the Q.60 value */
    return 0; }
static int64_t lnint56(int64_t m)                             /* ln m for an integer m >= 1, in Q.56 */
{
    int k = 0; while ((m >> k) >= 2) k++;                    /* 2^k <= m < 2^(k+1) */
    const int64_t y = (int64_t)((((unsigned __int128)m) << 60) >> k);   /* m/2^k in [1, 2) as Q.60 */
    s128 r = {0, 0}; s128_addv(&r, ln60(y)); s128_addp(&r, k, ln2_60());
    return s128_round_sh(r, 4);
}
static int64_t parse56(const char *s) { int neg = 0; if (*s == '-') { neg = 1; s++; } int64_t ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (*s++ - '0'); int64_t fr = 0; if (*s == '.') { char buf[40] = "0."; strncat(buf, s + 1, 30); fr = (parse60(buf) + 8) >> 4; } const int64_t v = (ip << 56) + fr; return neg ? -v : v; }
static int64_t fmt56v(int64_t v, char *out) { const uint64_t u = v < 0 ? -(uint64_t)v : (uint64_t)v, ip = u >> 56, fr = u & ((uint64_t)ONE56 - 1); uint64_t h, l; mul128(fr, 1000000000000ull, &h, &l); uint64_t q = (h << 8) | (l >> 56); if (l & ((uint64_t)1 << 55)) q++; uint64_t ipr = ip; if (q >= 1000000000000ull) { q -= 1000000000000ull; ipr++; } sprintf(out, "%s%llu.%012llu", v < 0 ? "-" : "", (unsigned long long)ipr, (unsigned long long)q); return 0; }
/* tanh and sech^2 of a Q.56 argument of any size, as Q.60 values: through e^{-2|x|} */
static void tanh_sech2(int64_t x56, int64_t *th, int64_t *sech2)
{
    const int neg = x56 < 0; const int64_t ax = neg ? -x56 : x56;
    s128 e = {0, 0}; s128_addp(&e, -32, ax);                    /* -2|x| in Q.60 (16 x 2) */
    const int64_t em = exp_of128(e);                            /* e^{-2|x|} <= 1 */
    const int64_t den = ONE60 + em, t = div60(ONE60 - em, den);
    *th = neg ? -t : t;
    *sech2 = div60(em, mul60(den, den) / 4);                             /* sech^2 = 4 e^{-2|x|}/(1 + e^{-2|x|})^2 = e^{-2|x|}/((1 + e^{-2|x|})^2/4) */
}
static int64_t lncosh56(int64_t x56)                          /* ln cosh x = |x| - ln 2 + ln(1 + e^{-2|x|}), Q.56 */
{
    const int64_t ax = x56 < 0 ? -x56 : x56;
    s128 e = {0, 0}; s128_addp(&e, -32, ax); const int64_t em = exp_of128(e);
    s128 r = {0, 0}; s128_addp(&r, 16, ax); s128_addv(&r, -ln2_60()); s128_addv(&r, ln60(ONE60 + em));
    return s128_round_sh(r, 4);
}
/* the step St(u) and its first two u-derivatives for a Q.60 u in (0, 1): sigma = 1/(1 + e^a), a = 1/u - 1/(1 - u);
   St' = sigma (1 - sigma) (1/u^2 + 1/(1-u)^2);  St'' = sigma (1 - sigma) [(1 - 2 sigma)(1/u^2 + 1/(1-u)^2)^2 - 2/u^3 + 2/(1-u)^3]
   1/u and 1/(1-u) in Q.56 (below 128 for u in [1/100, 99/100]); the polynomial factors in Q.40; outside that range St is 0 or 1 to 2^-60. */
static void step3(int64_t u, int64_t *st, int64_t *st1, int64_t *st2)
{
    if (u <= ONE60 / 100) { *st = 0; *st1 = 0; *st2 = 0; return; }
    if (u >= ONE60 - ONE60 / 100) { *st = ONE60; *st1 = 0; *st2 = 0; return; }
    const int64_t v = ONE60 - u;
    uint64_t rem; const int64_t iu = (int64_t)div128((uint64_t)1 << (116 - 64), 0, (uint64_t)u, &rem), iv = (int64_t)div128((uint64_t)1 << (116 - 64), 0, (uint64_t)v, &rem);   /* 1/u, 1/(1-u) in Q.56 */
    const int64_t a56 = iu - iv;                                 /* a = 1/u - 1/(1-u), Q.56, |a| < 100 */
    int64_t sig;                                                 /* sigma = 1/(1 + e^a) */
    if (a56 >= 0) { s128 e = {0, 0}; s128_addp(&e, -16, a56); const int64_t em = exp_of128(e); sig = div60(em, ONE60 + em); }
    else { s128 e = {0, 0}; s128_addp(&e, 16, a56); const int64_t ep = exp_of128(e); sig = div60(ONE60, ONE60 + ep); }
    const int64_t s1s = mul60(sig, ONE60 - sig);                 /* sigma (1 - sigma) */
    /* 1/u^2 + 1/(1-u)^2 and the cubes, in Q.40: 1/u <= 100 -> 1/u^2 <= 1e4, 1/u^3 <= 1e6 (Q.40: 1e6 x 2^40 = 1.1e18 fits) */
    const int64_t iu40 = iu >> 16, iv40 = iv >> 16;
    const int64_t iu2 = (int64_t)((((unsigned __int128)iu40 * (unsigned __int128)iu40)) >> 40), iv2 = (int64_t)((((unsigned __int128)iv40 * (unsigned __int128)iv40)) >> 40);
    const int64_t iu3 = (int64_t)((((unsigned __int128)iu2 * (unsigned __int128)iu40)) >> 40), iv3 = (int64_t)((((unsigned __int128)iv2 * (unsigned __int128)iv40)) >> 40);
    const int64_t p = iu2 + iv2;                                 /* Q.40, <= 2e4 */
    const int64_t p2 = (int64_t)((((unsigned __int128)p * (unsigned __int128)p)) >> 40);   /* Q.40, <= 4e8 (4e8 x 2^40 = 4.4e20 overflows!) -> keep p2 in Q.20 instead */
    (void)p2;
    /* St' = s1s p : a Q.60 x Q.40 product rounded to Q.60 */
    { s128 r = {0, 0}; s128_addp(&r, s1s, p); *st1 = s128_round_sh(r, 40); }
    /* St'' = s1s [(1 - 2 sigma) p^2 - 2 iu3 + 2 iv3]: p^2 in Q.20 (<= 4e8 x 2^20 = 4e14), the bracket in Q.20 */
    const int64_t p20 = p >> 20; const int64_t p2_20 = (int64_t)((((__int128)p20 * (__int128)p20)) >> 20);
    const int64_t one2s = ONE60 - 2 * sig;                        /* 1 - 2 sigma, Q.60 */
    s128 br = {0, 0}; s128_addp(&br, one2s, p2_20); const int64_t t1 = s128_round_sh(br, 60);   /* (1 - 2 sigma) p^2 in Q.20 */
    const int64_t bracket20 = t1 - 2 * (iu3 >> 20) + 2 * (iv3 >> 20);
    { s128 r = {0, 0}; s128_addp(&r, s1s, bracket20); *st2 = s128_round_sh(r, 20); }
    *st = sig;
}
typedef struct { int64_t PHI, PHI1, PHI2, w, lnT56, lnS56, lnH56; } prow;
static int64_t ln2g, ln2a, lnell, lna, ln4, eta_c, eta_d, eps, eps2, f, lneps;
static prow row_at(int64_t eta)
{
        /* the three pieces and their eta-derivatives (Q.56 for Phi; Q.60 for Phi', Phi'') */
        int64_t thT, s2T, thM, s2M, thS, s2S;
        tanh_sech2(eta + ln2g, &thT, &s2T); tanh_sech2(eta / 2, &thM, &s2M); tanh_sech2(eta - ln2a, &thS, &s2S);
        const int64_t PT = lnell - lncosh56(eta + ln2g), PM = ln4 + 2 * lncosh56(eta / 2), PS = lna - lncosh56(eta - ln2a);
        const int64_t PT1 = -thT, PT2 = -s2T, PM1 = thM, PM2 = s2M / 2, PS1 = -thS, PS2 = -s2S;
        /* the blends: u = (eta - eta_c + 1)/2 as Q.60 (from Q.56: x 16 / 2 = x 8) */
        int64_t Xc, Xc1, Xc2, Xd, Xd1, Xd2;
        { const int64_t uc56 = (eta - eta_c + ONE56) / 2, ud56 = (eta - eta_d + ONE56) / 2;
          const int64_t uc = uc56 <= 0 ? 0 : (uc56 >= ONE56 ? ONE60 : uc56 << 4), ud = ud56 <= 0 ? 0 : (ud56 >= ONE56 ? ONE60 : ud56 << 4);
          step3(uc, &Xc, &Xc1, &Xc2); step3(ud, &Xd, &Xd1, &Xd2); Xc1 /= 2; Xc2 /= 4; Xd1 /= 2; Xd2 /= 4; }                 /* d/d eta = (1/2) d/du */
        /* Phi_MS = (1 - Xd) Phi_M + Xd Phi_S and its derivatives (the difference Phi_S - Phi_M is small: Q.56 -> Q.60 fits) */
        const int64_t dMS56 = PS - PM; const int64_t dMS = dMS56 << 4;                                                        /* Phi_S - Phi_M, Q.60 (below 8 inside the windows) */
        const int64_t PMS = PM + s128_round_sh((s128){ 0, 0 }, 1) + ((int64_t)(mul60(Xd, dMS) >> 4));
        const int64_t PMS1 = mul60(Xd1, dMS) + mul60(ONE60 - Xd, PM1) + mul60(Xd, PS1);
        const int64_t PMS2 = mul60(Xd2, dMS) + 2 * mul60(Xd1, PS1 - PM1) + mul60(ONE60 - Xd, PM2) + mul60(Xd, PS2);
        const int64_t dTM56 = PMS - PT; const int64_t dTM = dTM56 << 4;
        const int64_t PHI = PT + ((int64_t)(mul60(Xc, dTM) >> 4));
        const int64_t PHI1 = mul60(Xc1, dTM) + mul60(ONE60 - Xc, PT1) + mul60(Xc, PMS1);
        const int64_t PHI2 = mul60(Xc2, dTM) + 2 * mul60(Xc1, PMS1 - PT1) + mul60(ONE60 - Xc, PT2) + mul60(Xc, PMS2);
        const int64_t w = 2 * PHI2 - 2 * mul60(PHI1, PHI1) + 2 * ONE60;                                                      /* the local conformal invariant */
        /* the Dirichlet data at z = eps: ln T, ln S, ln h_etaeta (ln eps and Phi in 128 bits where they add) */
        const int64_t g2t = rdiv64(3 * ONE60 - w, 6), g2s = rdiv64(w - 6 * ONE60, 12);                                        /* (3 - w)/6, (w - 6)/12 */
        const int64_t ct = ln60(ONE60 + mul60(eps2, g2t) + f), cs = ln60(ONE60 + mul60(eps2, g2s) - f);                       /* ln(1 + eps^2 g2 +- f), tiny */
        s128 lT = {0, 0}; s128_addv(&lT, -lneps); s128_addp(&lT, -16, PHI); s128_addv(&lT, ct / 2);                            /* ln T = -ln eps - Phi + ct/2 */
        s128 lS = {0, 0}; s128_addv(&lS, -lneps); s128_addv(&lS, cs / 2);                                                       /* ln S = -ln eps + cs/2 */
        s128 lH = {0, 0}; s128_addp(&lH, -2, lneps); s128_addv(&lH, ct);   /* 2 ln eps = -9.57 passes the ceiling: the product in 128 bits */                                                        /* ln h_etaeta = -2 ln eps + ct */
        const int64_t lnT56 = s128_round_sh(lT, 4), lnS56 = s128_round_sh(lS, 4), lnH56 = s128_round_sh(lH, 4);
    prow r = { PHI, PHI1, PHI2, w, lnT56, lnS56, lnH56 }; return r;
}
int main(int argc, char **argv)
{
    int64_t gamma = 1000, aa = 1000; int64_t lhat = ONE60 / 4; int GE = 76, every = 4;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "gamma=", 6)) gamma = atoll(argv[i] + 6);
        else if (!strncmp(argv[i], "a=", 2)) aa = atoll(argv[i] + 2);
        else if (!strncmp(argv[i], "lhat=", 5)) lhat = parse60(argv[i] + 5);
        else if (!strncmp(argv[i], "Ge=", 3)) GE = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "every=", 6)) every = atoi(argv[i] + 6);
    }
    char b1[48], b2[48], b3[48], b4[48], b5[48], b6[48], b7[48];
    ln2g = lnint56(2 * gamma); ln2a = lnint56(2 * aa); lnell = lnint56(gamma); lna = lnint56(aa); ln4 = lnint56(4);               /* Q.56 (re = 1) */
    eta_c = -rdiv64(2 * lnint56(gamma), 3); eta_d = rdiv64(2 * lnint56(aa), 3);                                                  /* the blend centres, Q.56 */
    const int64_t eta0 = -ln2g, eta1 = 0, DE = rdiv64(eta1 - eta0, GE);                                                            /* the quarter: throat centre to the photon sphere */
    /* the cutoff and the flux log term: eps = lhat/30; f = (3/(2 lhat^4)) eps^4 ln(eps/lhat) = (3/2) ln(1/30) / 30^4 (independent of lhat!) */
    eps = rdiv64(lhat, 30); eps2 = mul60(eps, eps);
    const int64_t ln30 = ln60(30 * (ONE60 / 32)) + 5 * ln2_60();                                                                   /* ln 30 = ln(30/32) + 5 ln 2 */
    f = -rdiv64(mul60(3 * (ONE60 / 2), ln30), 810000);                                                                           /* -(3/2) ln 30 / 30^4 = -6.30e-6 */
    lneps = ln60(eps);                                                                                                            /* ln eps (eps = 1/120: -4.787) */
    fmt56v(eta0, b1); fmt56v(eta_c, b2); fmt56v(eta_d, b3); fmt60(lhat, b4); fmt60(eps, b5); fmt60(f, b6);
    printf("WHP1DATA: P1's boundary data in the cylinder frame (the reviewer's verified data entry): gamma = ell/re = %lld, a/re = %lld, lhat = %s, the cutoff eps = lhat/30 = %s; eta from the throat centre %s to the photon sphere 0 in %d steps (the other quarter is the mirror when a = ell); blend centres eta_c = %s, eta_d = %s, half-width 1; the flux log term f = (3/2) ln(1/30)/30^4 = %s on the t, eta entries (minus on the S2 ones)\n",
           (long long)gamma, (long long)aa, b4, b5, b1, GE, b2, b3, b6);
    const int64_t worst = whq60_selftest("WHP1DATA"); printf("WHP1DATA selftest: largest deviation %lld units of 2^-60\n", (long long)worst);
    uint64_t chain = 0xcbf29ce484222325ull; int64_t wmax = 0, phimin = ONE56 * 100;
    for (int i = 0; i <= GE; i++) {
        const int64_t eta = eta0 + (int64_t)i * DE;                                                                          /* Q.56 */
        const prow r = row_at(eta);
        if (r.w > wmax) wmax = r.w; if (r.PHI < phimin) phimin = r.PHI;
        const int64_t raws[8] = { eta, r.PHI, r.PHI1, r.PHI2, r.w, r.lnT56, r.lnS56, r.lnH56 }; chain = fnv(chain, raws, sizeof raws);
        if (i % every == 0 || i == GE) {
            fmt56v(eta, b1); fmt56v(r.PHI, b2); fmt60(r.PHI1, b3); fmt60(r.w, b4); fmt56v(r.lnT56, b5); fmt56v(r.lnS56, b6); fmt56v(r.lnH56, b7);
            printf("WHP1DATA row %d: eta = %s, Phi = ln rho_o = %s, Phi' = %s, w = %s, ln T = %s, ln S = %s, ln h_etaeta = %s\n", i, b1, b2, b3, b4, b5, b6, b7);
        }
    }
    /* the entry's table (e), gamma = 1e3, a = ell: (eta, rho_o, w) at the window edges and centres, at the EXACT centres eta_c -+ 1 (the entry
       prints eta to four decimals; 3e-5 in eta is 3e-5 in Phi at these slopes) */
    { const char *nm[5] = { "throat window edge (eta_c - 1)", "throat-mouth centre (eta_c)", "mouth window edge (eta_c + 1)", "photon sphere (0)", "throat centre (-ln 2 gamma)" };
      const int64_t etas[5] = { eta_c - ONE56, eta_c, eta_c + ONE56, 0, eta0 }; const char *rho_s[5] = { "266.898", "100.874", "38.815", "4", "1000" }; const char *w_s[5] = { "0", "0.206", "0.309", "3", "0" };
      for (int k = 0; k < 5; k++) { const prow r = row_at(etas[k]); fmt56v(etas[k], b3); fmt56v(r.PHI, b1); fmt60(r.w, b2);
        printf("WHP1DATA check %s, eta = %s: Phi = %s (the entry: ln %s), w = %s (the entry: %s)\n", nm[k], b3, b1, rho_s[k], b2, w_s[k]); } }
    fmt60(wmax, b1); fmt56v(phimin, b2);
    printf("WHP1DATA: the largest w on the quarter %s (3 at the photon sphere), the smallest Phi %s (ln 4 = 1.386294361120 there)\n", b1, b2);
    printf("WHP1DATA: chain %016llx\n", (unsigned long long)chain);
    printf("WHP1DATA: done\n");
    return 0;
}
