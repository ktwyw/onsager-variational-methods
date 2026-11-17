# Changelog

## 0.1.0 (2026-11-17)

First release: Chapters 1–14 (Parts I–III), each an executed Jupyter notebook with
three figures and four exercises, and executed sample solutions to all 56 exercises.

### Part I. The principle

- Chapter 1: the Rayleighian of a particle; minimum dissipation in Stokes flow; trial
  functions for flow in a square duct.
- Chapter 2: Onsager's kinetic equations, reciprocal relations and the energy balance;
  a settling rod; two spheres with hydrodynamic interaction.
- Chapter 3: catalogue of free energies and dissipation functions; functional
  derivatives; constraints and Lagrange multipliers; the squeeze film.
- Chapter 4: the diffusion equation and the Einstein relation; the Stokes, Darcy and
  Brinkman equations.
- Chapter 5: the approximation method on diffusion into a half-space; a dynamical
  shape parameter; the Onsager–Machlup integral; what to trust in a trial solution.

### Part II. Liquids and interfaces

- Chapter 6: capillary rise and imbibition; the Lucas–Washburn law and Jurin's height
  from one slow variable; the exact solution of the kinetic equation; the energy
  balance; where inertia makes the overdamped description fail (Bosanquet's
  equation); imbibition into porous materials.
- Chapter 7: spreading of a droplet; the unbalanced Young force and the friction of a
  liquid wedge with its logarithmic cut-off; Tanner's law and partial wetting; the
  thin-film equation from the same two functions; comparison with its numerical
  solution, including the coefficient of the logarithm that a cap-shaped trial
  profile cannot give.
- Chapter 8: coalescence and sintering of viscous drops. Frenkel's model derived from
  the Rayleighian, with three ways of tying its flow to its shapes; Hopper's exact
  solution for two cylinders as the check; the diagnosis that Frenkel's flow leaves
  the filling of the neck without friction; a consistent trial (exact shapes and a
  potential flow) whose error has a guaranteed sign.
- Chapter 9: a drying droplet. The flow driven by evaporation in a cap; the Man–Doi
  equation from the Rayleighian with contact-line friction; paths of the solute and
  the deposit; the closed form for fast evaporation; a check against the thin-film
  equation with evaporation and solute transport.

### Part III. Soft materials

- Chapter 10: diffusion and sedimentation of colloids. The free energy (osmotic
  pressure and weight) and the friction between particles and solvent, with
  incompressibility imposed by the pressure; one hindrance function for collective
  diffusion and hindered settling; sedimentation equilibrium and the gravitational
  length. Mason and Weaver's exact solution tests two one-parameter trials, a smooth
  exponential and a falling front, each right in its own regime. A concentrated
  suspension of hard spheres checked against Kynch's front velocity and the exact
  equilibrium profile.
- Chapter 11: swelling and squeezing of a gel. The free energy of the network with a
  swelling pressure and the friction with the solvent; force balance and Darcy's law,
  and the diffusion equation for the network where the volume flux vanishes. A sheet
  on a substrate, swelling or squeezed under a porous piston (exact series). A
  swelling sphere: Tanaka and Fillmore's solution and the shear modulus; the slowest
  rate lies between π²K/fa² and 15K/fa², set by the bulk modulus; trial functions give
  upper bounds, one term (a uniform dilation) gives 15K/fa² for any G, three terms the
  exact rate to four digits. Finite-element checks and the stresses inside a swelling
  sphere.
- Chapter 12: phase separation. The Cahn–Hilliard equation from the Flory–Huggins
  free energy with a gradient term and the friction of interdiffusion; phase diagram,
  interface profile and tension; Cahn's linear theory of spinodal decomposition;
  trial functions of a few Fourier modes, from growth to saturation; coarsening and
  the Lifshitz–Slyozov law from a scaling argument with the two functions, compared
  with a simulation in two dimensions.
- Chapter 13: polymer dynamics. The Rouse equation from the free energy of entropic
  springs and the friction of the beads; normal modes and the Rouse time; time
  correlations of fluctuations from Onsager's regression hypothesis, checked against
  Brownian dynamics; stress relaxation; the dumbbell as a one-parameter trial and its
  improvement; the coil–stretch transition in elongational flow, and with finitely
  extensible springs the dumbbell against the full chain.
- Chapter 14: nematic liquid crystals. Frank elasticity, a field and rotational
  friction in a twist cell; the Fréedericksz threshold and the exact steady twist
  (elliptic integrals); switch-on and switch-off times; a one-mode trial,
  da/dt = −a + h²J₁(2a), exact at the threshold but beyond π/2 at strong fields; a
  director in shear flow, with the coupling of rotation and shear in the dissipation
  function and the Parodi relation as Onsager reciprocity; flow alignment and
  tumbling; wall layers in a sheared cell with their exact profile.

### The `onsager` module

- The Rayleighian, the friction matrix and the kinetic equations with constraints and
  Lagrange multipliers, in SymPy; the friction matrix of a trial profile.
- Reference solutions and solvers for the checks: conservative finite-volume solvers
  for diffusion (`smoluchowski_fv`), the thin-film equation (`thin_film_fv`), a drying
  drop (`drying_film_fv`) and a settling suspension (`sedimentation_fv`); exact flow
  in a rectangular duct; Hopper's exact coalescence; Mason and Weaver's solution; the
  exact linear swelling of a gel sheet and sphere, and a finite-element solver
  (`gel_fe`); a pseudo-spectral Cahn–Hilliard solver; Brownian dynamics of Rouse
  chains; the steady state of a FENE chain in elongational flow; a finite-difference
  solver for a nematic director and the exact Fréedericksz twist.
- Comparison helpers and the plotting style of the figures.
- 36 tests, run on Python 3.10 and 3.12.

### The book

- Three figures per chapter: two beside the theory and one collecting the checks; a
  gallery and two animations in the README.
- A Jupyter Book 2 site with an introduction, a notation page and a list of suggested
  references.
