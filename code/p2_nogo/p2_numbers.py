# P2 no-go: numbers at MM's parameters, and a worked circle profile. mpmath, 30 digits.
# MM inputs [M]: re = 1.5e7 m (MM:536, (3.26)); R5 = 5e-5 m (MM:442-443); smallest mouth separation d = 0.5 ly (MM:707);
# c2 and ell as functions of re, R5 (MM:541-545); E_Casimir = -c2/(8 ell) (MM:188-189).
# Physical constants [E: CODATA values from memory].
from mpmath import mp, mpf, pi, log, sqrt, quad, diff
mp.dps = 30

Gn = mpf('6.67430e-11'); c = mpf('299792458'); hbar = mpf('1.054571817e-34')
lp = sqrt(hbar * Gn / c**3)                  # Planck length, m
ly = mpf('9.4607304725808e15'); Msun = mpf('1.98892e30')
re = mpf('1.5e7'); R5 = mpf('5e-5')

lnr = log(re / R5)
G4c2 = 4 * pi * R5 * re / sqrt(lnr)          # G4 * c2 in m^2 (hbar = c = 1)
c2 = G4c2 / lp**2
ell = 4 / pi * re**2 / R5 * sqrt(lnr)
gam = ell / re
lB = (3 * lnr)**mpf('0.25') * sqrt(re * R5)  # MM (3.25)
print('ln(re/R5) = %s ; G4 c2 = %s m^2 ; c2 = %s (MM:558 says ~7e72) ; ell = %s m = %s ly ; gamma = %s ; l_B = %s m'
      % (mp.nstr(lnr, 8), mp.nstr(G4c2, 6), mp.nstr(c2, 4), mp.nstr(ell, 5), mp.nstr(ell / ly, 5), mp.nstr(gam, 4), mp.nstr(lB, 6)))

print('\n[2D CFT] ratio 16 pi G4 |rho_2D| / (2/r*^2) ~ (pi/3) G4 c2 / L^2 = (4 pi^2/3) R5 re /(L^2 sqrt ln(re/R5))')
for name, L in [('L = re (absurdly short circle)', re), ('L = 0.5 ly (MM:707 smallest d)', mpf('0.5') * ly),
                ('L = pi ell (time of flight, MM:233)', pi * ell)]:
    print('   %-38s ratio = %s ; R5 re / L^2 = %s' % (name, mp.nstr(pi / 3 * G4c2 / L**2, 3), mp.nstr(R5 * re / L**2, 3)))

dmin = mpf('0.5') * ly
rs = dmin / pi                               # r* of a round-S3-like exterior with antipodal mouths at d = pi r*
print('\n[4D CFT] (R5/r*)^2 at r* = re: %s ; at r* = d/pi = %s m (d = 0.5 ly): %s'
      % (mp.nstr((R5 / re)**2, 3), mp.nstr(rs, 4), mp.nstr((R5 / rs)**2, 3)))

for lnIR, lab in [(log(lB / R5), 'ln(l_B/R5)'), (lnr / 2, 'MM ln(re/R5)/2')]:
    kap = 1 + log(rs / re) / lnIR
    kap2 = 1 + log(pi * ell / re) / lnIR
    print('[running g4, (3.24)] kappa(r*) = 1 + ln(r*/re)/%s (= %s): at r* = d/pi %s ; at r* = pi ell %s'
          % (lab, mp.nstr(lnIR, 6), mp.nstr(kap, 4), mp.nstr(kap2, 4)))

Me = re * c**2 / Gn
Ecas = c2 * hbar / (8 * ell * c)
print('\n[masses] Me = re c^2/G = %s kg = %s Msun ; |E_Casimir| = c2 hbar/(8 ell c) = %s kg ; |E_bin| = %s kg (MM:560 ~5e9) ; |E_Casimir|/Me = %s ; 2/gamma^2 = %s'
      % (mp.nstr(Me, 4), mp.nstr(Me / Msun, 4), mp.nstr(Ecas, 4), mp.nstr(Ecas / 2, 4), mp.nstr(Ecas / Me, 4), mp.nstr(2 / gam**2, 4)))
need_rho = c**2 / (8 * pi * Gn) * (rs**2 - re**2) / rs**4
mag_rho = c**2 / (16 * pi * Gn) * 2 * re**2 / rs**4
need_M = (rs - re)**2 / rs * c**2 / (2 * Gn)
print('[closing at d = 0.5 ly, r* = d/pi] extra density at r* >= %s kg/m^3 (magnetic there: %s kg/m^3) ; extra mass (r*-re)^2 c^2/(2 G r*) = %s kg = %s Msun ; over |E_Casimir|: %s'
      % (mp.nstr(need_rho, 4), mp.nstr(mag_rho, 4), mp.nstr(need_M, 4), mp.nstr(need_M / Msun, 4), mp.nstr(need_M / Ecas, 4)))

# ---------------- worked profile, units re = 1 ----------------
# r'(r)^2 = F(r) = (1 - 1/r)^2 P(r), P = (rs^k - r^k)/(rs^k - 1): extremal-RN mouth (P -> 1) closing at rs.
print('\n[profile] r\'^2 = (1 - re/r)^2 (r*^k - r^k)/(r*^k - re^k), units re = 1')
for RS, k in [(mpf(10), mpf(1)), (mpf(10), mpf(2)), (mpf(10), mpf('0.25')), (mpf(100), mpf(1))]:
    Pf = lambda x: (RS**k - x**k) / (RS**k - 1)
    F = lambda x: (1 - 1 / x)**2 * Pf(x)
    Fp = lambda x: diff(F, x)
    R3 = lambda x: 2 * (1 - F(x)) / x**2 - 2 * Fp(x) / x       # r'' = F'/2
    Dl = lambda x: R3(x) - 2 / x**4                            # 16 pi G4 rho_x * re^2
    Xf = lambda x: quad(lambda y: 1 / sqrt(F(y)), [x, RS])     # proper distance to the maximum
    bound = 2 * (RS**2 - 1) / RS**4
    rpp = Fp(RS) / 2
    print(' r* = %s re, k = %s: at r*: 3R = %s = 2/r*^2 + 4|r\'\'|/r* = %s + %s ; magnetic 2re^2/r*^4 = %s ; required extra Delta = %s ; bound 2(r*^2-re^2)/r*^4 = %s ; Delta/bound = %s'
          % (mp.nstr(RS, 3), mp.nstr(k, 3), mp.nstr(R3(RS), 6), mp.nstr(2 / RS**2, 6), mp.nstr(-4 * rpp / RS, 6), mp.nstr(2 / RS**4, 4),
             mp.nstr(Dl(RS), 6), mp.nstr(bound, 6), mp.nstr(Dl(RS) / bound, 5)))
    I = quad(lambda x: x**2 * Dl(x) / 2, [1, 1 + (RS - 1) / 4, RS])
    mins = min(Dl(1 + (RS - 1) * mpf(j) / 400) for j in range(1, 401))
    print('    2 G4 * INT 4 pi r^2 rho_x dr = %s vs (r*-re)^2/r* = %s ; min Delta on the grid = %s (>= 0: extra density needed everywhere outside)'
          % (mp.nstr(I, 10), mp.nstr((RS - 1)**2 / RS, 10), mp.nstr(mins, 4)))
    if RS == 10 and k == 1:
        print('    r/re    X* - X (re)   3R re^2      2re^2/r^4    Delta re^2   (16 pi G4 rho_x re^2)')
        for x in [mpf('1.001'), mpf('1.01'), mpf('1.1'), mpf('1.5'), mpf(2), mpf(3), mpf(5), mpf(7), mpf(9), mpf('9.9'), RS]:
            Xd = Xf(x) if x < RS else mpf(0)
            print('    %-7s %-13s %-12s %-12s %-12s' % (mp.nstr(x, 4), mp.nstr(Xd, 6), mp.nstr(R3(x), 5), mp.nstr(2 / x**4, 5), mp.nstr(Dl(x), 5)))
        Lmin = 2 * Xf(mpf('1.001'))
        print('    circle length > 2 (X* - X(1.001 re)) = %s re ; 2D CFT share at r* <= (pi/3) G4 c2/L^2 = %s (MM re, R5)'
              % (mp.nstr(Lmin, 5), mp.nstr(pi / 3 * G4c2 / (Lmin * re)**2, 3)))
        print('    at MM\'s re = 1.5e7 m: rho_x(r*) = %s kg/m^3 ; extra mass = (r*-re)^2/(2 G r*) = %s Me = %s kg'
              % (mp.nstr(Dl(RS) / re**2 * c**2 / (16 * pi * Gn), 4), mp.nstr((RS - 1)**2 / RS / 2, 4), mp.nstr((RS - 1)**2 / RS / 2 * Me, 4)))
