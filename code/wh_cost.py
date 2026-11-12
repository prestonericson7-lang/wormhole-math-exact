"""wh_cost.py (11-12-2026, 7:35 AM) -- what Maldacena-Milekhin's wormhole costs at each mouth size, with the paper's human-traversability
condition (3.26) dropped: a portal nothing living goes through can be smaller. The construction itself (section 3.2) needs re >> R5 and
lhat = l_B/re << 1. Every row is COMPUTED on the PC from the paper's own formula chain with our exact constant kappa_flat in place of its
estimate (the reviewer's chain for the paper's eq. 3.27, our 1a and P3 results):
  ln_eff = (1/2) ln(re/R5) + (1/4) ln(6 ln_eff) + kappa_flat        (the running coupling's log, the paper's (3.24)-(3.25) made exact)
  Q2/R5^2 = (re/R5)^2 / ln_eff;  w = (3 + sqrt(9 + 96 Q2))/(4 Q2);  S* = R5/sqrt(w);  L3 = R5/sqrt(3 + w/4)
  c2 = 6 pi L3 S*^2/G5 with G5 = G4 R5/2 (3.23); l_B = sqrt(re R5) (6 ln_eff)^{1/4}; lhat = l_B/re
  ell = 16 re^3/(G4 c2) (2.12);  E_bin = -c2/(16 ell) (2.13);  gamma = ell/re (2.14)
  q (the integer magnetic charge, Dirac units) = re g4/(sqrt(pi) l_p) (2.5) with g4 = e (alpha = 1/137.036)
  B at the mouth: the extremal horizon field, B = c^2 sqrt(mu0/(4 pi G))/re (gauge-coupling independent; = (2.6) in SI)
  tidal acceleration on a 0.5 m body: a = c^2 (0.5 m)/re^2 (3.26); the proper traversal time pi re/c; the outside traversal time pi ell/c.
SI constants CODATA 2018; G4 = l_p^2 (hbar = c = 1), l_p = 1.616255e-35 m; R5 = 5e-5 m (the paper's choice, the experimental bound).
"""
import mpmath as mp
mp.mp.dps = 30
R5 = mp.mpf('5e-5'); lp = mp.mpf('1.616255e-35'); G4 = lp**2
c = mp.mpf('299792458'); hbar = mp.mpf('1.054571817e-34'); G = mp.mpf('6.67430e-11'); mu0 = mp.mpf('1.25663706212e-6')
alpha = 1/mp.mpf('137.035999084'); e4 = mp.sqrt(4*mp.pi*alpha)          # g4 = e in Heaviside-Lorentz natural units
KAPPA = mp.mpf('0.128210970524')
ly = mp.mpf('9.4607304725808e15'); yr = mp.mpf('31557600'); AU = mp.mpf('1.495978707e11'); gacc = mp.mpf('9.80665')
Bsch = mp.mpf('4.414e9')                                                  # the Schwinger field in tesla, a scale only
def row(re):
    lnS = mp.log(re/R5)
    L = mp.findroot(lambda x: x - (lnS/2 + mp.log(6*x)/4 + KAPPA), max(mp.mpf(1), lnS/2 + 1))
    Q2 = (re/R5)**2/L; w = (3 + mp.sqrt(9 + 96*Q2))/(4*Q2); Sst = R5/mp.sqrt(w); L3 = R5/mp.sqrt(3 + w/4)
    G5 = G4*R5/2; c2 = 6*mp.pi*L3*Sst**2/G5
    lB = mp.sqrt(re*R5)*(6*L)**mp.mpf('0.25'); lhat = lB/re
    ell = 16*re**3/(G4*c2); Eb = -c2/(16*ell)*hbar/c; gam = ell/re
    q = re*e4/(mp.sqrt(mp.pi)*lp)
    B = c**2*mp.sqrt(mu0/(4*mp.pi*G))/re
    tidal = c**2*mp.mpf('0.5')/re**2/gacc
    tprop = mp.pi*re/c; tout = mp.pi*ell/c
    return dict(re=re, L=L, lhat=lhat, Q2=Q2, c2=c2, lB=lB, ell=ell, Eb=Eb, gam=gam, q=q, B=B, tidal=tidal, tprop=tprop, tout=tout)
def fm(x, d=3): return mp.nstr(x, d)
def length(m):
    m = mp.mpf(m)
    if m >= ly: return fm(m/ly, 4) + " ly"
    if m >= AU: return fm(m/AU, 4) + " AU"
    if m >= 1000: return fm(m/1000, 4) + " km"
    if m >= 1: return fm(m, 4) + " m"
    if m >= mp.mpf('1e-3'): return fm(m*1000, 4) + " mm"
    return fm(m*1e6, 4) + " um"
def tm(s):
    s = mp.mpf(s)
    if s >= yr: return fm(s/yr, 4) + " yr"
    if s >= 3600: return fm(s/3600, 4) + " h"
    if s >= 1: return fm(s, 4) + " s"
    if s >= mp.mpf('1e-3'): return fm(s*1000, 4) + " ms"
    if s >= mp.mpf('1e-6'): return fm(s*1e6, 4) + " us"
    return fm(s*1e9, 4) + " ns"
sizes = [("100 R5, the construction's floor", mp.mpf('5e-3')), ("5 cm", mp.mpf('5e-2')), ("1 m", mp.mpf(1)), ("1 km", mp.mpf(1000)),
         ("1,000 km", mp.mpf('1e6')), ("the paper's: a human survives the tides", mp.mpf('1.5e7')), ("one astronomical unit", AU)]
print("| mouth radius re | ln_eff | lhat = l_B/re | charge q (Dirac units) | field at the mouth B | c2 | throat length ell | traversal time, outside (pi ell/c) | proper traversal time (pi re/c) | binding energy E_bin | tidal acceleration, 0.5 m body | what bites |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for name, re in sizes:
    r = row(re)
    flags = []
    if r['lhat'] > mp.mpf('0.1'): flags.append("lhat %s: the construction's small-lhat premise is marginal" % fm(r['lhat'], 2))
    if r['tidal'] > 20: flags.append("tides %s g: no human" % fm(r['tidal'], 2))
    if r['B'] > Bsch: flags.append("B is %s x the Schwinger field" % fm(r['B']/Bsch, 2))
    if not flags: flags.append("the paper's own point")
    print("| %s (%s) | %s | %s | %s | %s T | %s | %s | %s | %s | %s kg | %s g | %s |" % (length(re), name, fm(r['L'], 5), fm(r['lhat'], 3), fm(r['q'], 3), fm(r['B'], 3), fm(r['c2'], 4),
          length(r['ell']), tm(r['tout']), tm(r['tprop']), fm(r['Eb'], 4), fm(r['tidal'], 3), "; ".join(flags)))
print()
print("Scaling laws (lhat -> 0): q ~ re; B ~ 1/re; c2 ~ re/sqrt(ln_eff); ell ~ re^2 sqrt(ln_eff)/R5; the outside traversal time ~ re^2; |E_bin| ~ 1/(re ln_eff) (smaller mouths need MORE negative energy, not less); the tides ~ 1/re^2.")
r = row(mp.mpf('1.5e7'))
print("Check against the 1a page at the paper's size: ln_eff %s (page 14.4574779), l_B %s m (83.578), c2 %s (6.709e72), ell %s ly (3,257), gamma %s (2.054e12), E_bin %s kg (-4.79e9)" % (fm(r['L'], 9), fm(r['lB'], 6), fm(r['c2'], 5), fm(r['ell']/ly, 5), fm(r['gam'], 4), fm(r['Eb'], 4)))
