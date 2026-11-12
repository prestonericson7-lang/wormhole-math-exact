# P2 adversarial verifier: numbers. mpmath, 30 digits. Written fresh; does not import or run the worker's p2_numbers.py.
# MM inputs [M]: re = 1.5e7 m (MM:536), R5 = 5e-5 m (MM:442-443), d_min = 0.5 ly (MM:707, :745),
# c2 = 4 pi R5 re /(G4 sqrt(ln(re/R5))), ell = (4/pi) re^2 sqrt(ln)/R5 (MM:541-545), E_Cas = -c2/(8 ell) (MM:188).
# Constants [E: CODATA/IAU values from memory].
from mpmath import mp, mpf, pi, log, sqrt, quad, diff, nstr
mp.dps = 30
Gn = mpf('6.67430e-11'); cl = mpf('299792458'); hbar = mpf('1.054571817e-34')
ly = mpf('9460730472580800'); Msun = mpf('1.98847e30')
re = mpf('1.5e7'); R5 = mpf('5e-5')
lp2 = hbar * Gn / cl**3                        # G4 in m^2 (hbar = c = 1)
kg_per_m = cl**2 / Gn                          # geometric length -> kg
s = lambda v, n=4: nstr(v, n)

lnr = log(re / R5)
G4c2 = 4 * pi * R5 * re / sqrt(lnr); c2 = G4c2 / lp2
ell = 4 / pi * re**2 / R5 * sqrt(lnr); gam = ell / re
lB = (3 * lnr)**mpf('0.25') * sqrt(re * R5)
print('N1 ln(re/R5) =', s(lnr, 7), '; G4 c2 =', s(G4c2, 6), 'm^2 ; c2 =', s(c2, 4), '; ell =', s(ell, 5), 'm =', s(ell / ly, 5), 'ly ; gamma =', s(gam, 4), '; l_B =', s(lB, 5), 'm')
print('N1 consistency: ell - 16 re^3/(G4 c2) =', s(ell - 16 * re**3 / G4c2, 3), '; pi ell =', s(pi * ell / ly, 5), 'ly')

# 2D CFT density ratio at r*: 16 pi G |rho_2D| / (2/r*^2) with |rho_2D| ~ (pi c2/6)/(4 pi r*^2 L^2)
rat = lambda L: 16 * pi * (pi * G4c2 / 6) / (4 * pi * L**2) / 2
print('N2 ratio formula check (pi/3) G4c2/L^2 vs direct at L = re:', s(pi / 3 * G4c2 / re**2, 6), s(rat(re), 6))
for lab, L in [('L = re', re), ('L = 0.5 ly', mpf('0.5') * ly), ('L = pi ell', pi * ell), ('L = l_B(re)', lB)]:
    print('N2 %-12s ratio = %s' % (lab, s(rat(L), 3)))
Lone = sqrt(pi * G4c2 / 3)
print('N2 length at which the ratio would reach 1:', s(Lone, 4), 'm (below l_B =', s(lB, 4), 'm, where the 2D description itself fails)')

# total Casimir energy of a circle of length L relative to one mouth's mass Me = re/G4
for lab, E in [('flat-circle -pi c2/(6L)', lambda L: pi * G4c2 / (6 * L)), ("MM's -c2/(8 ell) with L = pi ell", lambda L: pi * G4c2 / (8 * L))]:
    print('N3 |E_Cas|/Me at L = re using %s: %s' % (lab, s(E(re) / re, 3)))
Ecas_kg = c2 * hbar / (8 * ell * cl); Me_kg = re * kg_per_m
print('N3 MM config: |E_Cas| =', s(Ecas_kg, 4), 'kg ; Me =', s(Me_kg, 4), 'kg =', s(Me_kg / Msun, 4), 'Msun ; ratio =', s(Ecas_kg / Me_kg, 4), '; 2/gamma^2 =', s(2 / gam**2, 4), '; |E_bin| =', s(Ecas_kg / 2, 4), 'kg')

# 4D CFT relative size, kappa
d = mpf('0.5') * ly; rs = d / pi
print('N4 r* = d/pi =', s(rs, 5), 'm ; (R5/re)^2 =', s((R5 / re)**2, 3), '; (R5/r*)^2 =', s((R5 / rs)**2, 3))
for lab, lnIR in [('ln(l_B/R5)', log(lB / R5)), ('ln(re/R5)/2', lnr / 2)]:
    print('N5 kappa with %-12s = %s: r* = d/pi -> %s ; r* = pi ell -> %s ; rho+sum p factor 1 - 1/(2 ln) = %s'
          % (lab, s(lnIR, 5), s(1 + log(rs / re) / lnIR, 4), s(1 + log(pi * ell / re) / lnIR, 4), s(1 - 1 / (2 * lnIR), 4)))
x = mpf(1); lnmin = lnr / 2
# does kappa(r) re^2 < r^2 for all r > re? f(x) = x^2 - 1 - ln(x)/lnmin, f(1) = 0, f'(x) = 2x - 1/(lnmin x) > 0
print('N5 f\'(1) = 2 - 1/ln_IR =', s(2 - 1 / lnmin, 4), '> 0, so r*^2 > kappa re^2 for every r* > re')

# closing at d = 0.5 ly
rho_need = (rs**2 - re**2) / (8 * pi * rs**4) * kg_per_m
rho_mag = re**2 / (8 * pi * rs**4) * kg_per_m
m_side = (rs - re)**2 / (2 * rs) * kg_per_m
print('N6 at r* = d/pi: extra density >=', s(rho_need, 4), 'kg/m^3 ; field density', s(rho_mag, 4), 'kg/m^3 ; extra mass per side (throat to r*) =',
      s(m_side, 4), 'kg =', s(m_side / Msun, 3), 'Msun ; both sides =', s(2 * m_side, 4), 'kg =', s(2 * m_side / Msun, 3), 'Msun')
Om = sqrt(re / d**3)
print('N6 MM:707 check: Kepler Omega at d = 0.5 ly =', s(Om, 4), '/m vs 1/ell =', s(1 / ell, 4), '/m')

# worked profile, units re = 1: F = r'^2 = (1 - 1/r)^2 (r*^k - r^k)/(r*^k - 1)
def prof(RS, k):
    F = lambda r: (1 - 1 / r)**2 * (RS**k - r**k) / (RS**k - 1)
    R3 = lambda r: 2 * (1 - F(r)) / r**2 - 2 * diff(F, r) / r
    Dl = lambda r: R3(r) - 2 / r**4
    return F, R3, Dl
RS = mpf(10)
F, R3, Dl = prof(RS, 1)
sq = sqrt(RS - 1)
Xan = lambda r: 2 * sq * sqrt(RS - r) + log((sq + sqrt(RS - r)) / (sq - sqrt(RS - r)))   # X* - X, closed form (k = 1)
print('N7 profile k=1, r* = 10 re:  r/re | X*-X closed form | X*-X quad | 3R re^2 | field | extra (16 pi G rho_x re^2)')
for r in [mpf('1.001'), mpf('1.01'), mpf('1.5'), mpf(3), mpf(7), mpf('9.9')]:
    Xq = quad(lambda y: 1 / sqrt(F(y)), [r, RS])
    print('N7   %-6s %-10s %-10s %-9s %-10s %s' % (s(r), s(Xan(r), 6), s(Xq, 6), s(R3(r), 4), s(2 / r**4, 3), s(Dl(r), 3)))
print('N7   10 (r*): 3R =', s(R3(RS), 6), '= 2/r*^2 + 4|r\'\'|/r* =', s(2 / RS**2, 4), '+', s(-2 * diff(F, RS) / RS, 4), '; field', s(2 / RS**4, 3),
      '; extra', s(Dl(RS), 6), '; bound 2(r*^2-1)/r*^4 =', s(2 * (RS**2 - 1) / RS**4, 6), '; ratio', s(Dl(RS) / (2 * (RS**2 - 1) / RS**4), 5))
I = quad(lambda r: r**2 * Dl(r) / 2, [1, 2, 5, RS])           # = 2 G * INT 4 pi r^2 rho_x dr per side
print('N7 per side: 2G INT 4 pi r^2 rho_x dr =', s(I, 15), 'vs (r*-1)^2/r* =', s((RS - 1)**2 / RS, 15), '; extra mass per side =', s(I / 2, 6), 're/G')
Iprop = quad(lambda r: r**2 * Dl(r) / 2 / sqrt(F(r)), [1, 1.5, 3, 7, 9.9, RS])
print('N7 per side, proper-volume version 2G INT rho_x dV =', s(Iprop, 6), '(dr-weighted is the Misner-Sharp mass; proper is larger)')
print('N7 at MM re: rho_x(r*) =', s(Dl(RS) / (16 * pi * re**2) * kg_per_m, 4), 'kg/m^3 ; extra mass per side =', s(I / 2 * re * kg_per_m, 4),
      'kg ; whole circle =', s(I * re * kg_per_m, 4), 'kg')
Lc = 2 * Xan(mpf('1.001'))
share_floor = pi / 3 * G4c2 / (Lc * re)**2
share_need = share_floor * (2 / RS**2) / Dl(RS)
print('N7 circle length > 2 (X* - X(1.001)) =', s(Lc, 5), 're ; 2D CFT share of the floor 2/r*^2 <=', s(share_floor, 3), '; of the extra needed <=', s(share_need, 3))
mins = min(Dl(1 + (RS - 1) * mpf(j) / 2000) for j in range(1, 2001))
print('N7 min extra on a 2000-point grid over (re, r*]:', s(mins, 4))
for k in [mpf(2), mpf('0.25'), mpf('0.01')]:
    Fk, R3k, Dlk = prof(RS, k)
    print('N8 k = %s: extra at r* over the bound = %s ; r\'\'(r*) = %s ; per side 2G INT = %s'
          % (s(k, 3), s(Dlk(RS) / (2 * (RS**2 - 1) / RS**4), 5), s(diff(Fk, RS) / 2, 4), s(quad(lambda r: r**2 * Dlk(r) / 2, [1, 2, 5, RS]), 10)))
