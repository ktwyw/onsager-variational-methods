# Changelog

## Unreleased (towards 0.1.0)

Part I, the principle: 5 executed notebooks.

- Chapter 1: the Rayleighian of a particle; minimum dissipation in Stokes flow; trial
  functions for flow in a square duct.
- Chapter 2: Onsager's kinetic equations, reciprocal relations and the energy balance;
  a settling rod; two spheres with hydrodynamic interaction.
- Chapter 3: catalogue of free energies and dissipation functions; functional
  derivatives; constraints and Lagrange multipliers; the squeeze film.
- Chapter 4: the diffusion equation and the Einstein relation; the Stokes, Darcy and
  Brinkman equations.
- Chapter 5: the approximation method on diffusion into a half-space; a dynamical
  shape parameter; the Onsager–Machlup integral.
- The `onsager` module (Rayleighian and kinetic equations with constraints in SymPy;
  friction matrix of a trial profile; a conservative finite-volume diffusion solver;
  flow in a rectangular duct; comparison helpers), with 16 tests.
- Three figures per chapter: two beside the theory (schematics, maps of the Rayleighian,
  fields and fluxes) and one collecting the checks; a gallery in the README.
- Four exercises at the end of each chapter, with executed sample solutions in
  `solutions/part1_principle/`; a notation page.

Still to come before 0.1.0: Part II (Chapters 6–9) and Chapters 10–11 of Part III.
