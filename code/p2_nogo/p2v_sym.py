# P2 adversarial verifier: independent sympy checks.
# Own curvature code (Christoffel -> Ricci), written fresh here; nothing imported from the worker's p2_checks.py,
# stage 1's g5lib.py or the stage-1 verifier's vlib.py.
# Conventions: signature (-+++), R^a_bcd = d_c Gam^a_db - d_d Gam^a_cb + Gam^a_ce Gam^e_db - Gam^a_de Gam^e_cb,
# R_bd = R^a_bad, G_ab = R_ab - g_ab R/2, Einstein: G_ab + Lam g_ab = 8 pi G T_ab.
import sympy as sp, time
t0 = time.time()

def christoffel(g, x):
    n = len(x); gi = sp.simplify(g.inv())
    Gm = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                s = sum(gi[a, d] * (sp.diff(g[d, c], x[b]) + sp.diff(g[d, b], x[c]) - sp.diff(g[b, c], x[d])) for d in range(n))
                s = sp.simplify(s / 2)
                Gm[a][b][c] = s; Gm[a][c][b] = s
    return Gm

def ricci(g, x):
    n = len(x); Gm = christoffel(g, x)
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(b, n):
            s = 0
            for a in range(n):
                s += sp.diff(Gm[a][d][b], x[a]) - sp.diff(Gm[a][a][b], x[d])
                for e in range(n):
                    s += Gm[a][a][e] * Gm[e][d][b] - Gm[a][d][e] * Gm[e][a][b]
            s = sp.simplify(s)
            Ric[b, d] = s; Ric[d, b] = s
    return Ric

def scalar(g, Ric):
    gi = g.inv(); n = g.shape[0]
    return sp.simplify(sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n)))

def out(*args):
    print(*args, flush=True)

X, th, ph, t, r_, chi = sp.symbols('X theta phi t r chi', real=True)
G, q, g4, re, a, Lam = sp.symbols('G q g4 r_e a Lambda', positive=True)
rf = sp.Function('r')(X); N = sp.Function('N')(X); be = sp.Function('beta')(X)
rp, rpp = sp.diff(rf, X), sp.diff(rf, X, 2)

# ---- sanity: round S2 of radius a has R = 2/a^2 (sign convention) ----
gS2 = sp.diag(a**2, a**2 * sp.sin(th)**2)
out('V0 sanity: R[S2 radius a] =', scalar(gS2, ricci(gS2, [th, ph])))

# ---- V1: 3R of dX^2 + r(X)^2 dOmega^2 ----
h = sp.diag(1, rf**2, rf**2 * sp.sin(th)**2)
R3 = scalar(h, ricci(h, [X, th, ph]))
out('V1 3R =', R3)
out('V1 3R - [2(1 - r\'^2)/r^2 - 4 r\'\'/r] =', sp.simplify(R3 - (2 * (1 - rp**2) / rf**2 - 4 * rpp / rf)))
out('V1 at a maximum (r\'=0, r\'\'=-k, k>=0): 3R - 2/r^2 =', sp.simplify(R3.subs(rpp, -sp.Symbol('k')).subs(rp, 0) - 2 / rf**2))

# ---- V2: static 4D metric -N^2 dt^2 + h ----
x4 = [t, X, th, ph]
g = sp.diag(-N**2, 1, rf**2, rf**2 * sp.sin(th)**2)
Ric = ricci(g, x4); R4 = scalar(g, Ric)
Gt = sp.simplify(Ric - g * R4 / 2)
gi = g.inv()
Gmix = sp.simplify(gi * Gt)                      # G^a_b
out('V2 2 G_nn - 3R (G_nn = G_tt/N^2) =', sp.simplify(2 * Gt[0, 0] / N**2 - R3))
D2N = sp.diff(N, X, 2) + 2 * rp * sp.diff(N, X) / rf    # Laplacian on the slice
out('V2 R_tt - N D^2 N =', sp.simplify(Ric[0, 0] - N * D2N))
birk = sp.simplify(Gmix[0, 0] - Gmix[1, 1])
out('V2 G^t_t - G^X_X =', birk)
out('V2 G^t_t - G^X_X - (2/r)(r\'\' - r\' N\'/N) =', sp.simplify(birk - 2 / rf * (rpp - rp * sp.diff(N, X) / N)))
# lapse equation from the trace-reversed Einstein equation with a general diagonal T (rho, p1, p2, p3) and Lambda
rho, p1, p2, p3 = sp.symbols('rho p1 p2 p3')
Tlow = sp.diag(rho * N**2, p1, p2 * rf**2, p3 * rf**2 * sp.sin(th)**2)
trT = sp.simplify(sum(gi[i, i] * Tlow[i, i] for i in range(4)))
Rtt_einstein = 8 * sp.pi * G * (Tlow[0, 0] - trT * g[0, 0] / 2) + Lam * g[0, 0]
out('V2 lapse eq: R_tt[Einstein]/N - N [4 pi G (rho + p1 + p2 + p3) - Lam] =',
    sp.simplify(Rtt_einstein / N - N * (4 * sp.pi * G * (rho + p1 + p2 + p3) - Lam)))

# ---- V3: Maxwell field F = (q/2) sin(th) dth^dph, L = -(1/(4 g4^2)) F^2 ----
def Fmat():
    F = sp.zeros(4, 4); F[2, 3] = q / 2 * sp.sin(th); F[3, 2] = -F[2, 3]; return F
F = Fmat()
F2 = sp.simplify(sum(F[m, n] * F[al, bb] * gi[m, al] * gi[n, bb] for m in range(4) for n in range(4) for al in range(4) for bb in range(4)))
TM = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        TM[m, n] = sp.simplify((sum(F[m, al] * F[n, bb] * gi[al, bb] for al in range(4) for bb in range(4)) - g[m, n] * F2 / 4) / g4**2)
TMmix = sp.simplify(gi * TM)
rhoM = sp.simplify(-TMmix[0, 0])
out('V3 F^2 =', F2, '; rho_M =', rhoM, '; p_X/rho =', sp.simplify(TMmix[1, 1] / rhoM),
    '; p_th/rho =', sp.simplify(TMmix[2, 2] / rhoM), '; p_ph/rho =', sp.simplify(TMmix[3, 3] / rhoM))
out('V3 rho + sum p over rho_M =', sp.simplify((-TMmix[0, 0] + TMmix[1, 1] + TMmix[2, 2] + TMmix[3, 3]) / rhoM))
re2 = sp.pi * G * q**2 / g4**2
out('V3 16 pi G rho_M - 2 re^2/r^4 with re^2 = pi G q^2/g4^2 :', sp.simplify(16 * sp.pi * G * rhoM - 2 * re2 / rf**4))

# independent route: vary sqrt(-g) L with respect to diagonal metric components (T^mm = (2/sqrt(-g)) d(sqrt(-g) L)/d g_mm)
c0, c1, c2_, c3 = sp.symbols('c0 c1 c2 c3')
Fth = q / 2 * sp.sin(th)
def Tdiag_from_action(Lfun):
    gg = [c0, c1, c2_, c3]
    sq = sp.sqrt(-c0 * c1 * c2_ * c3)
    s = 2 * Fth**2 / (c2_ * c3)             # F_ab F^ab for the diagonal metric
    dens = sq * Lfun(s)
    res = []
    for i in range(4):
        Tup = 2 / sq * sp.diff(dens, gg[i])
        Tmix = gg[i] * Tup
        res.append(Tmix)
    sub = {c0: -N**2, c1: 1, c2_: rf**2, c3: rf**2 * sp.sin(th)**2}
    return [sp.simplify(e.subs(sub)) for e in res]
TA = Tdiag_from_action(lambda s: -s / (4 * g4**2))
out('V3 from the action: T^a_a =', TA, '; -T^t_t - rho_M =', sp.simplify(-TA[0] - rhoM))

# ---- V4: extremal RN (2.4) with this normalization ----
x4r = [t, r_, th, ph]
for label, f in [('extremal', (1 - re / r_)**2)]:
    gr = sp.diag(-f, 1 / f, r_**2, r_**2 * sp.sin(th)**2)
    Ricr = ricci(gr, x4r); Rr = scalar(gr, Ricr); Gr = sp.simplify(Ricr - gr * Rr / 2)
    gri = gr.inv()
    F2r = sp.simplify(2 * F[2, 3]**2 * gri[2, 2] * gri[3, 3])
    Tr = sp.zeros(4, 4)
    for m in range(4):
        for n in range(4):
            Tr[m, n] = sp.simplify((sum(F[m, al] * F[n, bb] * gri[al, bb] for al in range(4) for bb in range(4)) - gr[m, n] * F2r / 4) / g4**2)
    E = sp.simplify((Gr - 8 * sp.pi * G * Tr).subs(q, re * g4 / sp.sqrt(sp.pi * G)))
    out('V4 %s RN: G_ab - 8 pi G T_ab with q = re g4/sqrt(pi G): all zero =' % label, E == sp.zeros(4, 4))
M_, Q2 = sp.symbols('M Q2', positive=True)
fg = 1 - 2 * M_ / r_ + Q2 / r_**2
gg_ = sp.diag(-fg, 1 / fg, r_**2, r_**2 * sp.sin(th)**2)
Ric_g = ricci(gg_, x4r); Gg = sp.simplify(Ric_g - gg_ * scalar(gg_, Ric_g) / 2)
Gtt_mix = sp.simplify((gg_.inv() * Gg)[0, 0])
sol = sp.solve(sp.Eq(Gtt_mix, -8 * sp.pi * G * q**2 / (8 * g4**2 * r_**4)), Q2)
out('V4 general RN: G^t_t = -8 pi G rho_M gives Q2 =', sol, '; extremal f = (1 - re/r)^2 needs M = re, Q2 = re^2, so re^2 =', sol[0])

# ---- V5: stationary flow along the circle (shift beta along X) ----
gsh = sp.Matrix([[-N**2 + be**2, be, 0, 0], [be, 1, 0, 0], [0, 0, rf**2, 0], [0, 0, 0, rf**2 * sp.sin(th)**2]])
Ricsh = ricci(gsh, x4); Rsh = scalar(gsh, Ricsh); Gsh = sp.simplify(Ricsh - gsh * Rsh / 2)
nlow = sp.Matrix([-N, 0, 0, 0]); nup = sp.simplify(gsh.inv() * nlow)
out('V5 n.n =', sp.simplify((nlow.T * nup)[0, 0]))
Gnn = sp.simplify((nup.T * Gsh * nup)[0, 0])
diffsh = sp.factor(sp.simplify(2 * Gnn - R3))
out('V5 2 G_nn - 3R =', diffsh)
out('V5 same at r\' = 0:', sp.simplify(diffsh.subs(rp, 0)))
# with a shift the static lapse equation R_tt = N D^2 N no longer holds: show the extra terms
xi = sp.Matrix([1, 0, 0, 0])
Rxin = sp.simplify((xi.T * Ricsh * nup)[0, 0])
resid = sp.simplify(Rxin - D2N)                  # static case: R_ab xi^a n^b = R_tt/N = D^2 N
out('V5 R_ab xi^a n^b - D^2 N with a shift =', resid)
out('V5 that residual + (1/r^2)(r^2 beta beta\'/N)\' (zero: a divergence, so only Komar survives) =',
    sp.simplify(resid + sp.diff(rf**2 * be * sp.diff(be, X) / N, X) / rf**2))
# the stationary slice identity: INT N [K_ij K^ij - K^2 + 4 pi G (rho + S) - Lam] dV = 0; here K_ij K^ij - K^2 is:
Ka, Kb = sp.diff(be, X) / N, rp * be / (rf * N)        # K^X_X and K^th_th = K^ph_ph for this metric
out('V5 K_ij K^ij - K^2 =', sp.expand(Ka**2 + 2 * Kb**2 - (Ka + 2 * Kb)**2), '(either sign: the lapse identity is a static-only result)')

# ---- V6: Misner-Sharp identity and the integrated deficit ----
m = rf / (2 * G) * (1 - rp**2)
out('V6 dm/dX - r^2 r\' 3R/(4G) =', sp.simplify(sp.diff(m, X) - rf**2 * rp * R3 / (4 * G)))
rs, rt, rr = sp.symbols('r_s r_t rr', positive=True)
Ix = (rs - rt) / (2 * G) - sp.integrate(4 * sp.pi * rr**2 * re**2 / (8 * sp.pi * G * rr**4), (rr, rt, rs))
out('V6 INT_{r_t}^{r*} 4 pi r^2 rho_x dr =', sp.factor(sp.simplify(Ix)),
    '; minus (r*-r_t)(1 - re^2/(r_t r*))/(2G) =', sp.simplify(Ix - (rs - rt) * (1 - re**2 / (rt * rs)) / (2 * G)),
    '; at r_t = re:', sp.factor(sp.simplify(Ix.subs(rt, re))))

# ---- V7: Einstein static universe (dust + Lambda) ----
x4c = [t, chi, th, ph]
ge = sp.diag(-1, a**2, a**2 * sp.sin(chi)**2, a**2 * sp.sin(chi)**2 * sp.sin(th)**2)
Rice = ricci(ge, x4c); Ge = sp.simplify(Rice - ge * scalar(ge, Rice) / 2)
rd = sp.Symbol('rho_d')
Te = sp.diag(rd, 0, 0, 0)
eqs = [sp.simplify(Ge[i, i] + Lam * ge[i, i] - 8 * sp.pi * G * Te[i, i]) for i in range(4)]
solE = sp.solve([eqs[0], eqs[1]], [rd, Lam], dict=True)[0]
out('V7 ESU: solution', solE, '; theta-theta and phi-phi residuals:', [sp.simplify(e.subs(solE)) for e in eqs[2:]])
re_slice = sp.simplify(R3.subs(rf, a * sp.sin(X / a)).doit())
out('V7 slice r = a sin(X/a): 3R =', re_slice, '; 16 pi G rho_d + 2 Lambda =', sp.simplify(16 * sp.pi * G * solE[rd] + 2 * solE[Lam]),
    '; supply over the floor 2/a^2 =', sp.simplify((16 * sp.pi * G * solE[rd] + 2 * solE[Lam]) / (2 / a**2)),
    '; lapse 4 pi G rho_d - Lambda =', sp.simplify(4 * sp.pi * G * solE[rd] - solE[Lam]))

# ---- V8: running coupling, L = -(1/4) s/g^2(s), 1/g^2 = (R5/g5^2) ln(L_IR/R5), L_IR = lam s^(-1/4) (L_IR ~ l_B ~ B^(-1/2)) ----
R5s, g5s, lam_ = sp.symbols('R5 g5 lam', positive=True)
def Lrun(s):
    return -s / 4 * (R5s / g5s**2) * (sp.log(lam_ / R5s) - sp.log(s) / 4)
TR = Tdiag_from_action(Lrun)
rhoR = sp.simplify(-TR[0])
lnIR = (sp.log(lam_ / R5s) - sp.log(2 * Fth**2 / (rf**4 * sp.sin(th)**2)) / 4)
ratio = sp.simplify((-TR[0] + TR[1] + TR[2] + TR[3]) / (2 * rhoR))
out('V8 running: p_X/rho =', sp.simplify(TR[1] / rhoR), '; (rho + sum p)/(2 rho) - (1 - 1/(2 ln(L_IR/R5))) =',
    sp.simplify(sp.expand_log(ratio - (1 - 1 / (2 * lnIR)), force=True)))

# ---- V9: Lambda alone with the field: the uniform tube R x S1 x S2 (Plebanski-Hacyan type) ----
r0 = sp.Symbol('r0', positive=True)
gt_ = sp.diag(-1, 1, r0**2, r0**2 * sp.sin(th)**2)
Rict = ricci(gt_, x4); Gtu = sp.simplify(Rict - gt_ * scalar(gt_, Rict) / 2)
gti = gt_.inv()
F2t = sp.simplify(2 * F[2, 3]**2 * gti[2, 2] * gti[3, 3])
Tt = sp.zeros(4, 4)
for mm in range(4):
    for nn in range(4):
        Tt[mm, nn] = sp.simplify((sum(F[mm, al] * F[nn, bb] * gti[al, bb] for al in range(4) for bb in range(4)) - gt_[mm, nn] * F2t / 4) / g4**2)
Et = [sp.simplify((Gtu[i, i] + Lam * gt_[i, i] - 8 * sp.pi * G * Tt[i, i]).subs(q, re * g4 / sp.sqrt(sp.pi * G))) for i in range(4)]
solT = sp.solve([Et[0], Et[2]], [r0, Lam], dict=True)
out('V9 tube: solutions', solT, '; all residuals:', [[sp.simplify(e.subs(s_)) for e in Et] for s_ in solT])

# ---- V10: where the trace anomaly enters. Generic L = -(1/4) s f(s), s = F^2 ----
fF = sp.Function('f')
TG = Tdiag_from_action(lambda s: -s * fF(s) / 4)
rhoG = sp.simplify(-TG[0]); trG = sp.simplify(sum(TG))
Xs = q**2 / (2 * rf**4)
out('V10 generic f: rho - X f(X)/4 =', sp.simplify(rhoG - Xs * fF(Xs) / 4), '; rho contains a derivative of f:', rhoG.has(sp.Derivative) or rhoG.has(sp.Subs))
out('V10 trace T - X^2 f\'(X) =', sp.simplify(trG - Xs**2 * sp.Subs(sp.Derivative(fF(sp.Symbol('u')), sp.Symbol('u')), sp.Symbol('u'), Xs).doit()))
lt = sp.Symbol('lambda_t', positive=True)       # ln(L_IR/R5) at the field strength considered
k5 = R5s / g5s**2
fX = k5 * lt; fpX = -k5 / (4 * Xs)               # log model: f = k5 (ln(lam/R5) - ln(X)/4), f' = -k5/(4X)
out('V10 log model: R4/(16 pi G rho_M) = -8 pi G T/(16 pi G rho) =', sp.simplify(-(Xs**2 * fpX) / (2 * Xs * fX / 4)))

# ---- V11: throat AdS2(R_A) x S2(S) with the log-running coupling, 4D effective theory ----
S, RA = sp.symbols('S R_A', positive=True)
Gth = Gmix.subs({N: sp.cosh(X / RA)}).subs(rf, S).doit()
Gth = sp.simplify(Gth)
rho_t = Xs.subs(rf, S) * fX / 4
pth_t = rho_t + (Xs.subs(rf, S))**2 * fpX.subs(rf, S) / 2
eqA = sp.Eq(Gth[0, 0], -8 * sp.pi * G * rho_t)
eqB = sp.Eq(Gth[2, 2], 8 * sp.pi * G * pth_t)
out('V11 G^t_t, G^X_X, G^th_th of AdS2 x S2 =', Gth[0, 0], ',', Gth[1, 1], ',', Gth[2, 2])
solS = sp.solve(eqA, S)
re_t2 = sp.pi * G * q**2 * k5 * lt               # re^2 = pi G q^2/g4^2 with 1/g4^2 = k5 lambda_t at the throat
out('V11 S^2 from G^t_t =', [sp.simplify(s_**2) for s_ in solS], '; minus re^2(throat coupling) =', [sp.simplify(s_**2 - re_t2) for s_ in solS])
solRA = sp.solve(eqB.subs(S, solS[0]), RA)
out('V11 (S/R_A)^2 =', [sp.simplify(solS[0]**2 / r_**2) for r_ in solRA], '(expect 1 - 1/(2 lambda_t))')
out('V11 G^X_X - 8 pi G p_X residual:', sp.simplify((Gth[1, 1] + 8 * sp.pi * G * rho_t).subs(S, solS[0])))

# ---- V12: the r* deficit with the running coupling: g(x) = x^2 - 1 - ln(x)/lambda, x = r*/re ----
xx, lam2 = sp.symbols('x lambda', positive=True)
gx = xx**2 - 1 - sp.log(xx) / lam2
out('V12 g(1) =', gx.subs(xx, 1), "; g'(x) =", sp.diff(gx, xx), "; g'(1) =", sp.diff(gx, xx).subs(xx, 1), "; g''(x) =", sp.diff(gx, xx, 2),
    '-> g > 0 for all x > 1 whenever lambda > 1/2')

# ---- V13: running coupling in the integrated form: extra field mass INT 4 pi r^2 (kappa - 1) re^2/(8 pi G r^4) dr ----
xs_ = sp.Symbol('x_s', positive=True)            # x* = r*/re >= 1
dI = sp.simplify(sp.integrate(4 * sp.pi * rr**2 * (sp.log(rr / re) / lam2) * re**2 / (8 * sp.pi * G * rr**4), (rr, re, xs_ * re)))
out('V13 running adds to the field\'s INT 4 pi r^2 rho dr:', dI, '; times 2 G lambda / re =', sp.simplify(dI * 2 * G * lam2 / re),
    '(= 1 - (1 + ln x*)/x*, in [0, 1) for x* >= 1, so the extra per half drops by less than re/(2 G lambda))')
out('V13 check 1 - (1 + ln x)/x at x = 1, 10, 1e8:', [sp.N(1 - (1 + sp.log(v)) / v, 6) for v in [1, 10, sp.Integer(10)**8]])

# ---- V14: the tt - XX equation. G^t_t - G^X_X = (2N/r)(r'/N)' = 8 pi G (T^t_t - T^X_X) = -8 pi G (rho + p_X) ----
out('V14 G^t_t - G^X_X - (2N/r)(r\'/N)\' =', sp.simplify(Gmix[0, 0] - Gmix[1, 1] - 2 * N / rf * sp.diff(rp / N, X)))
out('V14 rho + p_X: field const g4 =', sp.simplify(-TMmix[0, 0] + TMmix[1, 1]), '; field running =', sp.simplify(-TR[0] + TR[1]),
    '; generic f =', sp.simplify(-TG[0] + TG[1]), '; Lambda (T = -Lam g/(8 pi G)) = 0 identically')
out('V14 so with only the field and Lambda: (r\'/N)\' = 0, r\' = C N with N > 0: r monotonic or constant (the tube)')

print('done in %.1f s' % (time.time() - t0))
