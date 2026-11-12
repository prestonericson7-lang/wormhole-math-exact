# P2 no-go: exact symbolic checks (sympy). Own Christoffel/Ricci code, independent of g5lib.py and vlib.py.
# Conventions: signature (-,+,+,+); R^a_bcd = d_c Gam^a_bd - ...; R_bc = d_a Gam^a_bc - d_c Gam^a_ba + Gam^a_ad Gam^d_bc - Gam^a_cd Gam^d_ba.
# Einstein: G_mn + Lambda g_mn = 8 pi G T_mn. Maxwell action -(1/(4 g4^2)) INT sqrt(-g) F^2.
import sympy as sp
import time

t0 = time.time()
X, t, th, ph, chi = sp.symbols('X t theta phi chi', real=True)
G, q, g4, re, a, rhod, Lam = sp.symbols('G q g4 r_e a rho_d Lambda', positive=True)
r = sp.Function('r')(X)
N = sp.Function('N')(X)
beta = sp.Function('beta')(X)


def christoffel(g, x):
    n = len(x)
    gi = sp.simplify(g.inv())
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for A in range(n):
        for b in range(n):
            for c in range(n):
                Gam[A][b][c] = sp.simplify(sum(gi[A, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b])
                                                           - sp.diff(g[b, c], x[d])) for d in range(n)) / 2)
    return Gam


def ricci(g, x):
    n = len(x)
    Gam = christoffel(g, x)
    R = sp.zeros(n)
    for b in range(n):
        for c in range(n):
            R[b, c] = sp.simplify(sum(sp.diff(Gam[A][b][c], x[A]) - sp.diff(Gam[A][b][A], x[c])
                                      + sum(Gam[A][A][d] * Gam[d][b][c] - Gam[A][c][d] * Gam[d][b][A] for d in range(n))
                                      for A in range(n)))
    return R


def scalar(g, Ric):
    gi = g.inv()
    n = g.shape[0]
    return sp.simplify(sum(gi[i, j] * Ric[i, j] for i in range(n) for j in range(n)))


rp, rpp = sp.diff(r, X), sp.diff(r, X, 2)

# ---- (1) 3R of the slice dX^2 + r(X)^2 dOmega^2
h = sp.diag(1, r**2, r**2 * sp.sin(th)**2)
R3 = scalar(h, ricci(h, [X, th, ph]))
claim = 2 * (1 - rp**2) / r**2 - 4 * rpp / r
print('(1) 3R =', R3)
print('(1) 3R - [2(1 - r\'^2)/r^2 - 4 r\'\'/r] =', sp.simplify(R3 - claim))

# ---- (2) static 4D metric: Hamiltonian constraint 2 G_nn = 3R, and lapse equation R_tt = N D^2 N
g = sp.diag(-N**2, 1, r**2, r**2 * sp.sin(th)**2)
x4 = [t, X, th, ph]
Ric = ricci(g, x4)
Rs = scalar(g, Ric)
Gt = sp.simplify(Ric - g * Rs / 2)
Gnn = sp.simplify(Gt[0, 0] / N**2)
print('(2) static: 2 G_nn - 3R =', sp.simplify(2 * Gnn - R3))
D2N = sp.diff(r**2 * sp.diff(N, X), X) / r**2
print('(2) static: R_tt - N D^2 N =', sp.simplify(Ric[0, 0] - N * D2N))

# ---- (3) stationary SO(3) metric with a shift along the circle: -N^2 dt^2 + (dX + beta dt)^2 + r^2 dOmega^2
gs = sp.Matrix([[-N**2 + beta**2, beta, 0, 0], [beta, 1, 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2 * sp.sin(th)**2]])
Rics = ricci(gs, x4)
Rss = scalar(gs, Rics)
Gs = sp.simplify(Rics - gs * Rss / 2)
nvec = sp.Matrix([1 / N, -beta / N, 0, 0])
Gnn_s = sp.simplify((nvec.T * Gs * nvec)[0, 0])
diff_s = sp.simplify(2 * Gnn_s - R3)
print('(3) shift: 2 G_nn - 3R =', sp.factor(diff_s))
print('(3) shift: same at r\' = 0 :', sp.simplify(diff_s.subs(rp, 0)))
nn = sp.simplify((nvec.T * gs * nvec)[0, 0])
print('(3) shift: n.n =', nn)

# ---- (4) magnetic energy density: F = (q/2) sin(th) dth ^ dph on the static metric
F = sp.zeros(4)
F[2, 3] = q / 2 * sp.sin(th)
F[3, 2] = -F[2, 3]
gi = g.inv()
Fup = sp.simplify(gi * F * gi)          # F^{mn}
F2 = sp.simplify(sum(F[i, j] * Fup[i, j] for i in range(4) for j in range(4)))
T = sp.zeros(4)
for m in range(4):
    for n_ in range(4):
        T[m, n_] = sp.simplify((sum(F[m, A] * F[n_, B] * gi[A, B] for A in range(4) for B in range(4)) - g[m, n_] * F2 / 4) / g4**2)
rhoM = sp.simplify(T[0, 0] / N**2)
re2 = sp.pi * G * q**2 / g4**2
print('(4) F^2 =', F2, '; rho_M =', rhoM)
print('(4) 16 pi G rho_M - 2 re^2/r^4 (re^2 = pi G q^2/g4^2) =', sp.simplify(16 * sp.pi * G * rhoM - 2 * re2 / r**4))
print('(4) the field lines carry tension: T^X_X / rho_M =', sp.simplify(T[1, 1] / rhoM), '; T^th_th / rho_M =', sp.simplify(T[2, 2] / r**2 / rhoM))
# variation check of T_tt from the action, with the inverse metric as independent symbols (diagonal metric)
u0, u1, u2, u3 = sp.symbols('u0 u1 u2 u3')   # g^tt, g^XX, g^thth, g^phph ; u0 < 0
Fth = q / 2 * sp.sin(th)
X2 = 2 * Fth**2 * u2 * u3                     # F_mn F^mn for the magnetic field
sqrtg = 1 / sp.sqrt(-u0 * u1 * u2 * u3)
Lmax = -X2 / (4 * g4**2)
Ttt = sp.simplify(-2 / sqrtg * sp.diff(sqrtg * Lmax, u0))
rho_var = sp.simplify(Ttt * (-u0))            # rho = T_tt / N^2 with N^2 = -1/u0
print('(4) from the action: rho =', sp.simplify(rho_var.subs({u2: 1 / r**2, u3: 1 / (r**2 * sp.sin(th)**2)})), '(same as above:',
      sp.simplify(rho_var.subs({u2: 1 / r**2, u3: 1 / (r**2 * sp.sin(th)**2)}) - rhoM) == 0, ')')
# running coupling: L = -(1/4) f(X) X ; the energy density has no f' term
fX = sp.Function('f')
Xs = sp.Symbol('Xs', positive=True)
Lrun = -sp.Rational(1, 4) * fX(X2) * X2
Ttt_run = sp.simplify(-2 / sqrtg * sp.diff(sqrtg * Lrun, u0))
rho_run = sp.simplify((Ttt_run * (-u0)).subs({u2: 1 / r**2, u3: 1 / (r**2 * sp.sin(th)**2)}))
print('(4) running coupling L = -(1/4) f(F^2) F^2: rho =', rho_run, '; contains f\' :', rho_run.has(sp.Derivative))

# ---- (5) extremal RN (MM 2.4-2.5): G^m_n = 8 pi G T^m_n for all components iff re^2 = pi G q^2/g4^2
rr = sp.Symbol('rr', positive=True)
f = (1 - re / rr)**2
gR = sp.diag(-f, 1 / f, rr**2, rr**2 * sp.sin(th)**2)
xR = [t, rr, th, ph]
RicR = ricci(gR, xR)
RsR = scalar(gR, RicR)
GR = sp.simplify(RicR - gR * RsR / 2)
giR = gR.inv()
F2R = sp.simplify(2 * F[2, 3]**2 * giR[2, 2] * giR[3, 3])
TR = sp.zeros(4)
for m in range(4):
    for n_ in range(4):
        TR[m, n_] = sp.simplify((sum(F[m, A] * F[n_, B] * giR[A, B] for A in range(4) for B in range(4)) - gR[m, n_] * F2R / 4) / g4**2)
mixed = sp.simplify(giR * (GR - 8 * sp.pi * G * TR))
mixed = sp.simplify(mixed.subs(q, sp.sqrt(re**2 * g4**2 / (sp.pi * G))))
print('(5) extremal RN, G^m_n - 8 pi G T^m_n with q^2 = re^2 g4^2/(pi G):', [sp.simplify(mixed[i, i]) for i in range(4)],
      '; off-diagonal all zero:', all(sp.simplify(mixed[i, j]) == 0 for i in range(4) for j in range(4) if i != j))
# its static slice in proper distance: r' = sqrt(f) = 1 - re/r ; then 3R - 16 pi G rho_M = 0
rpR = 1 - re / rr
rppR = sp.diff(rpR, rr) * rpR
R3R = 2 * (1 - rpR**2) / rr**2 - 4 * rppR / rr
print('(5) RN slice: 3R - 2 re^2/r^4 =', sp.simplify(R3R - 2 * re**2 / rr**4))

# ---- (6) Einstein static universe with dust and Lambda
gE = sp.diag(-1, a**2, a**2 * sp.sin(chi)**2, a**2 * sp.sin(chi)**2 * sp.sin(th)**2)
xE = [t, chi, th, ph]
RicE = ricci(gE, xE)
RsE = scalar(gE, RicE)
GE = sp.simplify(RicE - gE * RsE / 2)
eqs = [sp.Eq(GE[0, 0] + Lam * gE[0, 0], 8 * sp.pi * G * rhod), sp.Eq(GE[1, 1] + Lam * gE[1, 1], 0)]
sol = sp.solve(eqs, [rhod, Lam], dict=True)[0]
print('(6) ESU: rho_dust =', sol[rhod], '; Lambda4 =', sol[Lam],
      '; thth, phph eqs hold:', sp.simplify(GE[2, 2] + sol[Lam] * gE[2, 2]) == 0 and sp.simplify(GE[3, 3] + sol[Lam] * gE[3, 3]) == 0)
rE = a * sp.sin(X / a)
R3E = sp.simplify(2 * (1 - sp.diff(rE, X)**2) / rE**2 - 4 * sp.diff(rE, X, 2) / rE)
print('(6) ESU slice r = a sin(X/a): 3R =', R3E, '; 16 pi G rho_dust + 2 Lambda4 =', sp.simplify(16 * sp.pi * G * sol[rhod] + 2 * sol[Lam]))
print('(6) ESU at the equator (r* = a): 3R - 2/r*^2 = 4|r\'\'|/r* =', sp.simplify(R3E - 2 / a**2),
      '; supply over the bound (16 pi G rho_dust + 2 Lambda4)/(2/a^2) =', sp.simplify((16 * sp.pi * G * sol[rhod] + 2 * sol[Lam]) / (2 / a**2)))
print('(6) ESU lapse equation: 4 pi G (rho + 3p) - Lambda4 =', sp.simplify(4 * sp.pi * G * sol[rhod] - sol[Lam]))

# ---- (7) integrated form: Misner-Sharp mass m = (r/2G)(1 - r'^2) has m' = (r^2 r'/(4G)) 3R
m = r / (2 * G) * (1 - rp**2)
print('(7) dm/dX - (r^2 r\'/(4G)) 3R =', sp.simplify(sp.diff(m, X) - r**2 * rp / (4 * G) * R3))
rt, rs = sp.symbols('r_t r_s', positive=True)
rr2 = sp.Symbol('rr2', positive=True)
Mmag = sp.integrate(4 * sp.pi * rr2**2 * re**2 / (8 * sp.pi * G * rr2**4), (rr2, rt, rs))
extra = sp.factor(sp.simplify((rs - rt) / (2 * G) - Mmag))
print('(7) INT_{r_t}^{r*} 4 pi r^2 rho_x dr (throat minimum r_t, maximum r* = r_s, both r\' = 0) =', extra)
print('(7) at r_t = re:', sp.factor(extra.subs(rt, re)))
# ---- (8) lapse identity source rho + sum p for the magnetic field, constant and running coupling (all T^m_n from the action)
Xv = sp.Symbol('Xv', positive=True)
sub_sph = {u0: -1 / N**2, u1: 1, u2: 1 / r**2, u3: 1 / (r**2 * sp.sin(th)**2)}
for lab, Lag in [('constant g4', -X2 / (4 * g4**2)), ('running f(F^2)', -sp.Rational(1, 4) * fX(X2) * X2)]:
    Tm = [sp.simplify(-2 / sqrtg * sp.diff(sqrtg * Lag, uu) * uu) for uu in (u0, u1, u2, u3)]   # mixed T^m_m = g^mm T_mm
    rho8 = sp.simplify(-Tm[0].subs(sub_sph))
    pX, pth, pph = [sp.simplify(Tm[i].subs(sub_sph)) for i in (1, 2, 3)]
    print('(8) %s: rho = %s ; p_X = %s ; p_th = %s ; p_ph - p_th = %s ; rho + sum p = %s'
          % (lab, rho8, pX, pth, sp.simplify(pph - pth), sp.simplify(rho8 + pX + pth + pph)))
# with MM's log running 1/g4^2 = (R5/g5^2) ln(L_IR/R5), L_IR ~ l_B ~ (F^2)^(-1/4): f = (R5/g5^2)(l0 - ln(F^2/X0)/4)
R5s, g5s, l0, X0 = sp.symbols('R5 g5 l0 X0', positive=True)
frun = R5s / g5s**2 * (l0 - sp.log(Xv / X0) / 4)
srcrun = sp.simplify(frun * Xv / 2 + Xv**2 * sp.diff(frun, Xv))
print('(8) log running: rho + sum p = f X/2 + X^2 f\' =', sp.factor(srcrun), '; over 2 rho = f X/2:', sp.simplify(srcrun / (frun * Xv / 2)))
print('done in %.1f s' % (time.time() - t0))
