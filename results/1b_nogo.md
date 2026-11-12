# 1b, a first result: the paper's circle of mouths has no static solution on its fine-tuned brane

**What this page is.** Maldacena and Milekhin's humanly traversable wormhole (arXiv:2008.06618) lives on a Randall-Sundrum brane whose
tension is tuned so that the four-dimensional cosmological constant vanishes; the two mouths are the ends of a tube of magnetic field lines,
and the configuration they draw closes the field lines in a circle through the exterior. Before any of that geometry is solved numerically
(the plan's 1b), its first stage of derivation asked whether the whole closed circle can be static at all in the paper's own theory. The
answer, derived by the reviewer (the independent AI session that checks every program and page here) and checked by two further
independent rebuilds with their own symbolic and high-precision code, is no: at the circle's widest point the slice's curvature exceeds what
the field's energy can supply, and the integrated lapse equation shows that extra positive energy alone cannot make up the difference,
nor a positive cosmological constant alone; the brane would need both, or matter with rho + 3p below zero. The paper's own rotationally
symmetric example (the Einstein static universe with antipodal mouths) is a four-dimensional solution only with dust and a positive
cosmological constant, which it has; the paper uses it as boundary data for the five-dimensional problem rather than as a solution of its
four-dimensional theory, and so does this repository's next computation (P1 on `PLAN.md`).

Nothing here was computed on the two machines: it is a derivation with exact identities and estimated magnitudes, and it is labelled as the
review labels it. The scripts that check every formula and print every number are in `code/p2_nogo/`, with their printed output beside
them. The text below is the review's, as delivered at 2:07 AM on 11-12-2026, unchanged.

---

## P2: the whole circle on the fine-tuned brane has no static solution

This was derived in the review's first stage of 1b and checked in the review. The review's independent verifier reproduced both key formulas with its own code. A second independent check rebuilt every formula and number on this page.

Labels: [M] MEASURED (a program's output or a published value), [C] COMPUTED (a derivation or program, with its inputs given), [E] ESTIMATE (judgement).

MM = Maldacena and Milekhin, "Humanly traversable wormholes", arXiv:2008.06618.

**Claim.**
- P2 is the whole magnetic field-line circle: MM's throat (their eq. 2.7), the two mouths (2.4), and an exterior that closes the circle.
- It is static and SO(3)-symmetric. It sits on MM's Randall-Sundrum II brane, whose tension is fine-tuned so that the 4D cosmological constant Lambda4 is zero (MM section 3.2).
- In MM's theory and regime it has no static solution.
- Extra positive energy with rho + 3p >= 0 is not enough by itself, and a positive Lambda4 (a tension above the fine-tuned value) is not enough by itself either. The brane would need extra matter that is not of the Lambda type, and also a positive Lambda4 or other matter with rho + 3p < 0. The Einstein static universe has both.

**The argument [C; every step is exact].**
1. A static slice has zero extrinsic curvature. So the time-time Einstein equation is the Hamiltonian constraint 3R = 16 pi G4 rho + 2 Lambda4. Here 3R is the slice's scalar curvature, rho is the energy density seen by static observers, and Lambda4 = 0.
2. With SO(3), the slice is dX^2 + r(X)^2 dOmega^2, where X is the proper distance along the field lines. Then 3R = 2(1 - r'^2)/r^2 - 4 r''/r.
3. The circle closes, so r(X) has a maximum r*, where r' = 0 and r'' <= 0. There 3R = 2/r*^2 + 4|r''|/r* >= 2/r*^2.
4. MM's field is F = (q/2) sin(th) dth^dph, with the Maxwell term -(1/(4 g4^2)) INT F^2 of their eq. 2.1.
   - Its energy density is rho_M = q^2/(8 g4^2 r^4).
   - So 16 pi G4 rho_M = 2 re^2/r^4, with re^2 = pi G4 q^2/g4^2 (their eq. 2.5).
   - That is less than 2/r^2 wherever r > re.
5. So at the maximum the slice needs an extra energy density rho_x(r*) >= (r*^2 - re^2)/(8 pi G4 r*^4) on top of the field's. Once r* >> re, that is about 1/(8 pi G4 r*^2).

**Integrated form [C, exact].**
- The Misner-Sharp mass m = (r/(2 G4))(1 - r'^2) obeys dm/dX = 4 pi r^2 r' rho.
- From the throat's minimum r_t to the maximum, this gives INT 4 pi r^2 rho_x dr = (r* - r_t)(1 - re^2/(r_t r*))/(2 G4). For r_t = re, that is (r* - re)^2/(2 G4 r*).
- So each side of a closed exterior, from the throat to r*, must hold about the mass of a black hole of its own radius. The whole exterior must hold about twice that.

**Why extra positive energy alone does not help [C for the exact identity; E for the size of the quantum terms].**
- The same equations give the lapse equation D^2 N = N [4 pi G4 (rho + p1 + p2 + p3) - Lambda4]. Here N is the redshift factor, which is nowhere zero (MM section 3.2).
- On a closed slice the left side integrates to zero. So INT N [4 pi G4 (rho + sum p) - Lambda4] dV = 0.
- The field lines are under tension: p_X = -rho_M and p_th = p_ph = +rho_M. So rho + sum p = 2 rho_M > 0.
- MM's running coupling (their eq. 3.24) lowers this by the factor 1 - 1/(2 ln(L_IR/R5)), about 0.96. It stays positive.
- Its N-weighted integral is about the two mouths' mass, 2 re/G4 [E]. Dust, or any matter with rho + 3p >= 0, only adds to it.
- With Lambda4 = 0, only the quantum sector's negative rho + sum p could balance it.
  - At MM's numbers, the whole Casimir energy is 5e-25 of one mouth's mass [C].
  - For any circle longer than re, it is at most about 4e-12 of it [E].
- This identity needs neither SO(3) nor a maximum of r. So it also rules out circles that never grow past re.
- A positive Lambda4 alone is not enough either [C].
  - The field and Lambda4 both have rho + p_X = 0, and then the static equations make r'/N constant. So with only those two, a closed circle has constant r. The uniform tube R x S1 x S2 with r = sqrt(2) re and Lambda4 = 1/(4 re^2) is such a solution, and it has no throat.
  - Where r is largest, the extra matter needs rho + p_X = |r''|/(4 pi G4 r*).

**What MM's sector can supply.**
- **The 4D Maxwell energy.** It is already counted in step 4.
  - With MM's running coupling, the density at r* is larger by kappa = 1 + ln(r*/re)/ln(L_IR/R5).
  - That is 2.3 to 2.4 at r* = 1.5e15 m, and 3.1 to 3.2 at r* = pi ell [C from (3.24)-(3.25); E in the choice of L_IR].
  - The bound becomes (r*^2 - kappa re^2)/(8 pi G4 r*^4). A logarithm cannot close a gap of (r*/re)^2.
- **The 2D CFT along the field lines.**
  - Its Casimir energy is negative in total: E = -c2/(8 ell) (MM eq. 2.11).
  - Whatever its local sign, its density at r* is of order (pi c2/6)/(A L^2). Here A = 4 pi r*^2 is the transverse area. L is the circle's length: strictly, N(r*) times its optical length, the integral of dX/N, which includes the throat's delay.
  - So 16 pi G4 |rho_2D| / (2/r*^2) ~ (pi/3) G4 c2/L^2 = (4 pi^2/3) R5 re/(L^2 sqrt(ln(re/R5))) [E]. This uses c2 = 4 pi R5 re/(G4 sqrt(ln(re/R5))) from MM section 3.3.
  - At MM's numbers, re = 1.5e7 m and R5 = 50 micrometres. That gives G4 c2 = 1.83e3 m^2 and c2 = 7.0e72, as in their eq. 3.27 [C].
  - The ratio is then [E; the arithmetic is exact]:
    - 8.6e-29 for L = 0.5 light year, the smallest separation MM allow (appendix B);
    - 2.2e-37 for L = pi ell = 9,800 light years, the time-of-flight length (ell = 3,113 light years);
    - only 8.5e-12 even for L = re.
  - Its local anomaly part is bounded the same way, with L replaced by the redshift length N/|N'|. That length is at least re in MM's throat and mouths [E].
- **The 4D CFT (the AdS5 bulk).** It adds curvature terms of relative size about (R5/r*)^2: at most 1.1e-23 (at r* = re), and 1.1e-39 at r* = 1.5e15 m [E].
- **In mass units [C].** Suppose the circle closes with the mouths 0.5 light year apart, so that r* = d/pi = 1.5e15 m for a round exterior.
  - It then needs at least 2.4e-5 kg/m^3 of extra density at r*. The field's own density there is 2.3e-21 kg/m^3.
  - It also needs about 1.0e42 kg (5e11 solar masses) of extra mass on each side of r*, so about 2.0e42 kg (1e12 solar masses) in the whole exterior.
  - For comparison, the whole Casimir energy is -1.0e10 kg, and each mouth is Me = re/G4 = 2.0e34 kg.

**MM's SO(3) example, as a check.**
- MM's example (at the end of their section 3.2) puts the mouths at antipodal points of the Einstein static universe R x S3.
- They say it is not a solution of their 4D gravity theory, and they use it as boundary data for the 5D problem.
- This repository's P1 does the same: the 4D metric is fixed at the AdS5 boundary, and there is no dynamical brane.
- As a 4D solution, the Einstein static universe of radius a needs dust of density rho_d = 1/(4 pi G4 a^2) and Lambda4 = 1/a^2 [C]. It passes both tests:
  - On the slice r = a sin(X/a), 3R = 6/a^2 everywhere. At the equator (r* = a, r'' = -1/a), the supply 16 pi G4 rho_d + 2 Lambda4 = 4/a^2 + 2/a^2 is three times the bound 2/a^2.
  - With N = 1 the lapse equation holds, because 4 pi G4 rho_d - Lambda4 = 0.
- The dust pays for the bending, and the positive Lambda4 is what lets a closed slice be static.

**What is exact, what is estimated, what is assumed.**
- Exact: the inequality 3R(r*) >= 2/r*^2; the formula for 3R; 16 pi G4 rho_M = 2 re^2/r^4 for MM's action 2.1; the static constraint; and the two integral identities.
- Estimated: the bounds on the 2D and 4D CFT terms, and kappa.
- Assumed:
  - Static. The inequality needs it only at r*, and there the constraint also holds with a stationary flow along the circle [C]. The lapse identity needs a static metric everywhere: with a flow it gains a term K_ij K^ij - K^2, which can have either sign [C].
  - SO(3), and the circle closes.
  - The effective 4D theory holds at r* >> R5, with the CFT in its vacuum-like state. A hot CFT (a horizon in the bulk) would be extra energy with rho + 3p > 0.
  - MM's regime, L >> re >> R5.
  - The deficit is order one only when r* - re >> re/(2 ln(L_IR/R5)), about 4% of re. Closer to re than that, the log-running and anomaly corrections, of relative size about 1/(2 ln(L_IR/R5)), are as large as the deficit and would have to be computed. These are the terms behind the exact brane relation R4 = Q2/S^4. In MM's configuration the exterior has r* >> re, so the no-go holds there [C; the review's stage-1 verifier, MATH_WORK_01.md stage-1 entry].
  - Lambda4 = 0, and no other matter on the brane.

**Checks.**
- The review's sympy script uses its own curvature code. It confirms exactly [C]:
  - the 3R of the warped slice;
  - rho_M from the stress tensor, and again by varying the Maxwell action;
  - that MM's extremal black hole (2.4) solves every Einstein-Maxwell component with this normalization exactly when re^2 = pi G4 q^2/g4^2, which is their (2.5);
  - 2 G_nn = 3R and R_tt = N D^2 N on the static slice;
  - that the stationary shift term is proportional to r';
  - the Einstein static universe's dust and Lambda4;
  - both integral identities;
  - rho + sum p for a constant and for a running coupling.
- The review's independent verifier printed the same 3R and 16 pi G4 rho_M - 2 re^2/r^4 = 0 with its own code [M].
- **A worked profile [C].** Take an extremal-black-hole mouth that turns over at r* = 10 re: r'^2 = (1 - re/r)^2 (r* - r)/(r* - re). The field carries the mouth almost alone, and every point beyond it needs extra density.

| r/re | distance to r* (re) | 3R re^2 | field 2 re^2/r^4, times re^2 | extra needed, 16 pi G4 rho_x re^2 |
|---|---|---|---|---|
| 1.01 | 26.18 | 1.922 | 1.922 | 6.4e-5 |
| 1.5 | 21.74 | 0.4390 | 0.3951 | 0.0439 |
| 3 | 18.64 | 0.1015 | 0.0247 | 0.0768 |
| 7 | 11.71 | 0.0508 | 0.00083 | 0.0500 |
| 10 (r*) | 0 | 0.0380 | 0.0002 | 0.0378 |

- At r*, 3R = 0.020 + 0.018: the 2/r*^2 floor plus the bending.
  - The field gives 0.0002, so the extra needed is 0.0378/re^2.
  - That is 1.9 times the bound 0.0198/re^2.
- On each side, the extra mass integrates to (r* - re)^2/(2 G4 r*) = 4.05 re/G4, to every printed digit.
- At MM's re, that is 4.5e9 kg/m^3 at r* and 8.2e34 kg on each side (1.6e35 kg for the whole circle). The 2D CFT could supply at most 1.4e-15 of the density needed at r* [E].
- With r'^2 proportional to (r*^k - r^k) instead, k = 2 gives 2.7 times the bound and k = 1/4 gives 1.5. The bound is reached only as the top of the profile flattens (r'' -> 0).

---

## Files

`code/p2_nogo/p2_checks.py` (the exact identities in sympy, with its output `p2_checks.out`), `code/p2_nogo/p2_numbers.py` (the numbers at the
paper's parameters and the worked profile, in 30-digit arithmetic, with `p2_numbers.out`), and the second builder's independent rebuild
`code/p2_nogo/p2v_sym.py` and `p2v_num.py` with their outputs. Run any of them with Python 3 and the sympy and mpmath packages.
