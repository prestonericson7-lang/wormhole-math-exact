/* wh_ode.c (11-12-2026, 2:20 AM) -- WORMHOLE MATH, THE EXACT ODE INTEGRATOR FOR A VECTOR STATE: wh_flow.c version 3's fixed-step RK4 in Q.40
 * (int64 with 128-bit intermediates, products and divisions rounded to nearest; each component's update h (k1 + 2 k2 + 2 k3 + k4)/6 rounded
 * once with its remainder carried to the next step) generalized from the flow's three variables to a state y[0..nv-1] with a right-hand side
 * chosen by model=. Built so that the near-brane throat of 1b (an ODE integrated from the brane and matched to the 1a flow) is a right-hand
 * side, a start and a stop rule here, and nothing else changes.
 *
 * ACCEPTANCE: model=flow is wh_flow.c v3 (the two Maldacena-Milekhin flows, (3.19) flat and (A.30) sphere) through this frame: the same
 * stage arithmetic in the same order, the state hashed after every step in component order, so its chains must be wh_flow.c's bit for bit:
 * flat 24d2f870a2abea59, sphere 85a3df31c75d8935. model=decay is y' = -y from y(0) = 1 to t = 8: the integrator's global error at h = 2^-12 is
 * below a unit of 2^-40, so y(8) must equal the header's expq(-8) to a few units (a check of the frame that owes nothing to the flow).
 *
 * Usage: whode model=flow|decay [h=12]      (the step is 2^-h; h = 12 is the flow's)
 * Build: cc -O2 -std=gnu11 -o whode wh_ode.c      (whq.h beside it) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "whq.h"

#define NVMAX 8
static int HS = 12;                                   /* the step is 2^-HS */
static q REM[NVMAX];                                  /* the carried remainders, zero at the start of each integration */
static q upd(q X, q *R)                               /* the increment h X / 6 for h = 2^-HS, rounded once (half up), remainder kept in *R */
{
    const q D = (q)6 << HS;
    const q T = *R + X + D / 2;                       /* exact: |X| < 2^47 */
    q inc = T / D;
    if (T % D < 0) inc--;                             /* floor, also for T < 0 */
    *R = T - D / 2 - inc * D;                         /* = R + X - inc D, in [-D/2, D/2) */
    return inc;
}
typedef void (*rhs_fn)(const q *y, q *dy, int which);
static void rk4v(q *y, int nv, q h, rhs_fn f, int which)   /* one step: the stages as wh_flow.c's rk4, component by component */
{
    q k1[NVMAX], k2[NVMAX], k3[NVMAX], k4[NVMAX], y2[NVMAX], y3[NVMAX], y4[NVMAX]; const q h2 = h / 2;
    f(y, k1, which);
    for (int i = 0; i < nv; i++) y2[i] = y[i] + mulq(h2, k1[i]);
    f(y2, k2, which);
    for (int i = 0; i < nv; i++) y3[i] = y[i] + mulq(h2, k2[i]);
    f(y3, k3, which);
    for (int i = 0; i < nv; i++) y4[i] = y[i] + mulq(h, k3[i]);
    f(y4, k4, which);
    for (int i = 0; i < nv; i++) y[i] += upd(k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i], &REM[i]);
}

/* ---- model=flow: wh_flow.c v3's two flows (sphere = 0 flat (3.19)/(A.28-29), sphere = 1 (A.30) with l_B-hat = 1) ---- */
static q Lprime(q s, q sp, int sphere)
{
    const q e4 = expq(-4 * s), e2 = sphere ? expq(-2 * s) : 0;
    const q under = 3 * mulq(sp, sp) + 6 * ONE - 3 * e4 + e2;
    return -2 * sp + sqrtq(under);
}
static q spp(q s, q sp, int sphere)
{
    const q e4 = expq(-4 * s), e2 = sphere ? expq(-2 * s) : 0;
    const q lp = Lprime(s, sp, sphere);
    return 4 * ONE + e2 - 4 * e4 - 2 * mulq(lp, sp) - 2 * mulq(sp, sp);
}
static void flow_rhs(const q *y, q *dy, int sphere)   /* y = (sigma, sigma', lambda) */
{
    dy[0] = y[1]; dy[1] = spp(y[0], y[1], sphere); dy[2] = Lprime(y[0], y[1], sphere);
}
static const q START[2][3] = { { 1072463084LL, 2813843829LL, -15236010697140LL }, { -67469346530LL, 3031326166LL, -15938931920661LL } };
static const q GAMMA[2][2] = { { -5894808496714LL, 4315299320483LL }, { -6103177591682LL, 4955530204813LL } };

/* ---- model=decay: y' = -y ---- */
static void decay_rhs(const q *y, q *dy, int which) { (void)which; dy[0] = -y[0]; }

int main(int argc, char **argv)
{
    const char *model = "flow";
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "model=", 6)) model = argv[i] + 6;
        else if (!strncmp(argv[i], "h=", 2)) HS = atoi(argv[i] + 2);
    }
    if (HS < 4 || HS > 30) { printf("WHODE: h= must be 4..30\n"); return 1; }
    const q h = ONE >> HS;
    char a[48], b[48], c[48], d[48], e[48], f[48];
    printf("WHODE: Q.%d fixed point (int64, 128-bit intermediates, round to nearest); RK4 with the step 2^-%d, each component's update rounded once with its remainder carried; model = %s\n", FB, HS, model);
    if (!strcmp(model, "flow")) {
        for (int sphere = 0; sphere < 2; sphere++) {
            uint64_t hash = 0xcbf29ce484222325ull;
            q y[3] = { START[sphere][0], START[sphere][1], START[sphere][2] }, rho = -8 * ONE;
            memset(REM, 0, sizeof REM);
            const char *name = sphere ? "sphere" : "flat";
            fmt(y[0], a); fmt(y[1], b); fmt(y[2], c);
            printf("WHODE %s: start on the unstable manifold at rho = -8: sigma = %s, sigma' = %s, lambda = %s\n", name, a, b, c);
            int steps = 0;
            while (rho < 12 * ONE) {
                rk4v(y, 3, h, flow_rhs, sphere);
                hash = fnv(hash, &y[0], 8); hash = fnv(hash, &y[1], 8); hash = fnv(hash, &y[2], 8);
                rho += h; steps++;
                if (y[0] > 16 * ONE) break;
            }
            fmt(rho, a); fmt(y[1], b); fmt(Lprime(y[0], y[1], sphere), c); fmt(y[0] - rho, d); fmt(y[2] - rho, e); fmt(y[2] - y[0], f);
            printf("WHODE %s: END at rho = %s after %d steps: sigma' = %s (UV: 1), lambda' = %s (UV: 1), sigma - rho = %s, lambda - rho = %s, lambda - sigma = %s; chain %016llx (wh_flow.c v3: %s)\n",
                   name, a, steps, b, c, d, e, f, (unsigned long long)hash, sphere ? "85a3df31c75d8935" : "24d2f870a2abea59");
            fmt(y[0] - rho + GAMMA[sphere][0], a); fmt(y[2] - rho + GAMMA[sphere][1], b);
            printf("WHODE %s: invariants K_sigma = %s (%s), K_lambda = %s (%s)\n", name, a, sphere ? "0.377362726342" : "0.407800114190", b, sphere ? "-0.297030895775" : "-0.271699935370");
        }
    } else if (!strcmp(model, "decay")) {
        uint64_t hash = 0xcbf29ce484222325ull;
        q y[1] = { ONE }, t = 0; memset(REM, 0, sizeof REM); int steps = 0;
        while (t < 8 * ONE) { rk4v(y, 1, h, decay_rhs, 0); hash = fnv(hash, &y[0], 8); t += h; steps++; }
        const q ex = expq(-8 * ONE); fmt(y[0], a); fmt(ex, b); fmt(t, c);
        printf("WHODE decay: y(%s) = %s after %d steps; expq(-8) = %s; difference %lld units of 2^-40 (the RK4 error at this step is below one unit: a few units are the roundings); chain %016llx\n",
               c, a, steps, b, (long long)(y[0] - ex), (unsigned long long)hash);
    } else { printf("WHODE: unknown model %s\n", model); return 1; }
    printf("WHODE: done\n");
    return 0;
}
