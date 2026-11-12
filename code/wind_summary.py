#!/usr/bin/env python3
"""wind_summary.py (11-12 12:55 AM) -- Kobrin-Schuster-Yao's size-winding summary numbers from v6's WINDING lines (the machines print the
exact per-size P(l) and q(l) in Q.60; this script does the fits the paper defines, in float64, on the PC):
  r_l = |q(l)| / P(l);  r = SUM |q(l)|;  W = SUM q(l) (real);  rbar = (r - W)/(1 - W);
  C(m) = |SUM_l q(l) e^{-i m l}|, C* = max_m C(m), m* = argmax (period pi, quoted in (-pi/2, pi/2]);
  M = 2|q(l1)| + 2|q(l2)| - r for the two largest |q|;  L = max(W, M);  chi = (C* - L)/(r - L)  (KSY eq. D1).
Also 2|C(mu)| at m = -mu/(4 vnorm) from the same q(l), to compare with the machine's own C(mu) line.
usage: python wind_summary.py <whsyk6 log> [vnorm]     (vnorm defaults to N; model 1 uses 10)"""
import sys, re, cmath, math
log = sys.argv[1]; vn_arg = int(sys.argv[2]) if len(sys.argv) > 2 else None
LINE = re.compile(r'^WHSYK6 N=(\d+) seed=(\d+) beta=\S+ t=([\d.]+) op=(\d+): WINDING(.*?) \|\| ')
SIZE = re.compile(r'\| l=(\d+) P (-?[\d.]+) q (-?[\d.]+) (-?[\d.]+) i')
CL = re.compile(r'^WHSYK6 N=(\d+) seed=(\d+) beta=\S+ t=([\d.]+) op=(\d+) mu=(-?[\d.]+): C\(mu\) from q\(l\) (-?[\d.]+) (-?[\d.]+) i, 2\|C\| (-?[\d.]+)')
cm = {}
for raw in open(log, encoding='utf-8', errors='replace'):
    m = CL.match(raw)
    if m: cm[(int(m.group(1)), int(m.group(2)), float(m.group(3)), int(m.group(4)), int(round(float(m.group(5)))))] = float(m.group(8))
print('| N | seed | t | op | W | r | rbar | C* | m* | N m* | chi | 2|C(-12)| / C* | 2|C(+12)| / C* | machine 2|C(-12)|, 2|C(+12)| |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for raw in open(log, encoding='utf-8', errors='replace'):
    m = LINE.match(raw)
    if not m: continue
    N, seed, t, op, body = int(m.group(1)), int(m.group(2)), float(m.group(3)), int(m.group(4)), m.group(5)
    vn = vn_arg or N
    q = {}; P = {}
    for sm in SIZE.finditer(body):
        l = int(sm.group(1)); P[l] = float(sm.group(2)); q[l] = complex(float(sm.group(3)), float(sm.group(4)))
    ls = sorted(q)
    W = sum(q[l] for l in ls).real; r = sum(abs(q[l]) for l in ls)
    rbar = (r - W) / (1 - W) if W != 1 else float('nan')
    def C(mm): return abs(sum(q[l] * cmath.exp(-1j * mm * l) for l in ls))
    # the argmax over one period pi, coarse then refined
    best = max((C(-math.pi / 2 + math.pi * k / 20000), -math.pi / 2 + math.pi * k / 20000) for k in range(20001))
    m0 = best[1]; step = math.pi / 20000
    for _ in range(40):
        step /= 2
        cand = max((C(m0 + d), m0 + d) for d in (-step, 0.0, step)); m0 = cand[1]
    Cstar, mstar = C(m0), m0
    if mstar > math.pi / 2: mstar -= math.pi
    if mstar <= -math.pi / 2: mstar += math.pi
    top = sorted(ls, key=lambda l: -abs(q[l]))[:2]
    M = 2 * abs(q[top[0]]) + 2 * abs(q[top[1]]) - r
    L = max(W, M); chi = (Cstar - L) / (r - L) if r != L else float('nan')
    c_neg = C(12.0 / (4 * vn)); c_pos = C(-12.0 / (4 * vn))     # m = -mu/(4 vnorm)
    mach = cm.get((N, seed, t, op, -12)), cm.get((N, seed, t, op, 12))
    print('| %d | %d | %.1f | %d | %.9f | %.6f | %.4f | %.6f | %+.4f | %+.3f | %.4f | %.4f | %.4f | %s |' % (
        N, seed, t, op, W, r, rbar, Cstar, mstar, N * mstar, chi, c_neg / Cstar, c_pos / Cstar,
        ('%.6f, %.6f' % mach) if None not in mach else '(no C(mu) lines)'))
    print('|   per size: ' + '  '.join('l=%d p %.5f |q| %.5f arg-arg1 %+.3f r_l %.4f' % (l, P[l], abs(q[l]), cmath.phase(q[l] / q[ls[0]]), abs(q[l]) / P[l] if P[l] else float('nan')) for l in ls) + ' |')
