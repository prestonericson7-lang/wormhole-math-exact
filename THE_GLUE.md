# The glue: what holds a wormhole open in known physics, and why those wormholes are not time machines

Preston Ericson, 2026-11-11. This page is the reading that led to the three computations in `PLAN.md`. Every paper is labelled
"read whole" or "abstract only" so the reader knows how far each claim here was checked at the source. Quotations are short and
attributed; the papers are linked, not re-hosted.

## 1. The starting sentence

The U.S. Defense Intelligence Agency's 2010 reference document on traversable wormholes (Eric W. Davis, "Traversable Wormholes,
Stargates, and Negative Energy", DIA, 2010; released under the Freedom of Information Act in 2022 with the other AAWSAP reference
documents; read whole) states, as a property of the theory, that a traversable wormhole can be made into a time machine. That is
classical general relativity: Morris, Thorne and Yurtsever, Phys. Rev. Lett. 61, 1446 (1988). The 2010 desk could not know what came
after 2016. What came after is the subject of this page.

## 2. What holds a wormhole open in known physics: negative Casimir energy of quantum fields

- Gao, Jafferis and Wall, "Traversable Wormholes via a Double Trace Deformation", arXiv:1608.05687 (2016; abstract only): couple the
  two mouths of an Einstein-Rosen bridge by a quantum interaction and the quantum matter stress tensor acquires negative average null
  energy, whose gravitational back-reaction renders the bridge traversable. They read it as a form of quantum teleportation and state
  that it cannot be used to violate causality.
- Maldacena, Milekhin and Popov, "Traversable wormholes in four dimensions", arXiv:1807.04726 (2018; abstract only): a four-dimensional
  wormhole in Einstein-Maxwell theory plus charged massless fermions; the fermions give a negative Casimir-like energy that makes the
  wormhole possible; it is a pair of entangled near-extremal black holes; a long wormhole that does not lead to causality violations;
  embeddable in the Standard Model only by making its size small compared to the electroweak scale (microscopic).

## 3. The human-sized one: Maldacena and Milekhin, "Humanly traversable wormholes", arXiv:2008.06618, Phys. Rev. D 103, 066007 (2021); read whole

Allowed in the quantum theory, with one catch in their words: the time it takes to go through the wormhole should be longer than the
time it takes to travel between the two mouths on the outside. Their footnote: this also implies that they cannot be converted into
time machines. Their recipe, which they themselves call science fiction: a dark sector (a conformal field theory with a dark U(1)) in
the Randall-Sundrum II model with one warped extra dimension; two opposite-charged extremal magnetic black holes whose field lines run
through the throat; the two-dimensional Casimir energy along each field line as the glue.

Their requirements (their eq. 3.26-3.27): a mouth radius re > 1.5 x 10^7 m so that a 0.5 m body feels under 20 g (our arithmetic: a
mass of about 2 x 10^34 kg, ten thousand suns, each mouth; they say the mouths resemble intermediate-mass charged black holes from
outside); the extra dimension at its largest allowed size, R5 = 50 micrometres (the bound of the 2020 Eot-Wash gravity test, their ref.
18); then c2 ~ 7 x 10^72, a throat length ell ~ 3,000 light years, a boost gamma ~ 2 x 10^12, binding energy ~ -5 x 10^9 kg; proper
time through it about a second (pi re), outside time pi ell: tens of thousands of years for somebody looking from outside.

The practical problems in their own section 3.4: a CMB photon falling in is boosted by gamma and seen inside by gamma squared, so the
huge black hole would have to sit inside a refrigerator; the dark sector must sit below 10^-26 eV; producing the wormhole requires
topology change, which they call difficult. Their conclusion: no plausible mechanism for formation has been given; only that these are
configurations allowed by the equations. The free-fermion version needs more than 10^52 species, too many if the cutoff of the local
field theory is to be above the TeV scale, since they expect fewer than M_Planck^2 / TeV^2 ~ 10^32.

What they left uncomputed is the subject of computations 1a and 1b in `PLAN.md`: their section 3.2 says the solution exists and
that one could imagine finding it numerically in a more explicit form, and that they will not do it there.

## 4. The time machine, both doors

- The classical door (the DIA's sentence) is Morris, Thorne and Yurtsever 1988: given a traversable wormhole, relative motion or a
  gravitational potential difference between the mouths produces closed timelike curves.
- The quantum constructions above are all "long" wormholes, slower than the outside path, and so are not time machines by
  construction. Every one of them says so.
- Hawking, "Chronology protection conjecture", Phys. Rev. D 46, 603 (1992); read whole: as a spacetime is warped toward closed timelike
  curves, vacuum polarization makes the expectation value of the energy-momentum tensor large, and the back-reaction would prevent
  closed timelike curves from appearing; against the Planck-scale cutoff of Kim and Thorne he argues that a cutoff in the invariant
  distance would give a metric perturbation of order one, changing the spacetime and probably making the Cauchy horizon impossible
  to cross. The conjecture, in his words, is that the laws of physics prevent the appearance of closed timelike curves. Status: a
  conjecture, unproven; a full theory of quantum gravity is needed to settle it.

So: allowed by the classical equations, forbidden by every quantum construction found so far, and argued against by Hawking. The
question is open, and this page marks it open, neither proven nor dismissed.

## 5. The laboratory side since 2010

- Wilson et al., "Observation of the dynamical Casimir effect in a superconducting circuit", Nature 479, 376 (2011), arXiv:1105.4714
  (abstract only): the effect the 2010 DIA review expected "before 2012" was observed.
- Squeezed vacuum, the 2010 review's first "negative energy generator", is standard equipment in the gravitational-wave detectors
  (general knowledge; to be cited to the LIGO squeezing papers).
- Jafferis et al., "Traversable wormhole dynamics on a quantum processor", Nature 612, 51 (2022) (abstract and the institutions'
  summaries): a nine-qubit, 164-gate circuit on a superconducting processor running a machine-learned sparse version of the SYK
  model, the holographic dual of a traversable wormhole, with a qubit inserted on one side and read on the other under a
  negative-energy coupling. Kobrin, Schuster and Yao, arXiv:2302.07897 (2023; pulled, not yet read whole) argue the learned model
  shows no gravitational physics. It is a simulation of the dual quantum system, not a spacetime. This dispute is computation 2.

## 6. The warp side

Bobrick and Martire, "Introducing physical warp drives", Class. Quantum Grav. 38, 105009 (2021), arXiv:2102.06824 (read whole): any
warp drive is a shell of ordinary or exotic material moving inertially; a general positive-energy spherically symmetric subluminal
class; the Alcubierre negative energy reduced by flattening and by an optimal wall shape. Fuchs et al., "Constant velocity physical
warp drive solution", Class. Quantum Grav. 41, 095013 (2024), arXiv:2405.02709 (read whole): a matter shell of 2.4 Jupiter masses,
inner radius 10 m, outer 20 m, with a 0.02 c shift inside, satisfying all four energy conditions to their numerical floor, with no
acceleration mechanism; their toolkit Warp Factory is public. Lentz, arXiv:2006.07125 (2020; abstract only) claims positive-energy
superluminal solitons; Bobrick and Martire write that their conclusions do not support that claim. Checking the 2021 and 2024
claims exactly, and computing the shift limit the 2024 paper leaves as future work, is computation 3.

## 7. Where the energy scale meets the collider

The glue of section 3 lives at the TeV scale three times over: the TeV cutoff caps the fermion count (10^32 against the 10^52
needed), the five-dimensional Planck scale must be well above 1 TeV, and the Standard-Model-embeddable wormhole must be smaller than
the electroweak scale. The Large Hadron Collider is the instrument at that scale. ATLAS, arXiv:1707.04147 (2017; abstract only; 36.7
fb^-1 at 13 TeV): spin-2 resonances with a warped extra-dimension model as benchmark and large-extra-dimension signals, no
significant deviation from the Standard Model. CMS, arXiv:1705.01403 (2017; abstract only; 2.3 fb^-1): semiclassical black holes
and string balls excluded up to 9.5 TeV and quantum black holes up to 9.0 TeV in models with extra dimensions. Later runs tightened
these bounds (to be cited). As the record stands, the extra dimensions the glue needs have been looked for up to about 10 TeV and
not found.

## 8. What this repository computes

The three computations of `PLAN.md`: the explicit Maldacena-Milekhin solution; the full SYK traversable-wormhole protocol at sizes
past nine qubits, to test the disputed 2022 signature; and the energy conditions of the 2021 and 2024 warp metrics, exactly. Results
others can check. None of them builds a door.
