# Variational methods in soft matter

*The Onsager principle, taught through worked problems, with every variational result
checked against an exact or numerical solution.*

The slow motion of soft matter is decided by a balance between thermodynamic forces
and friction. The Onsager principle states that balance as a minimum: the rates of
change are those that minimise the Rayleighian, the sum of the dissipation function
$\Phi$ and the rate of change of the free energy $A$. Two scalar functions replace a
list of forces; the equations that result obey the second law and Onsager's
reciprocal relations by construction; and a minimum can be approached with trial
functions when the exact solution is out of reach.

Each chapter is an executed Jupyter notebook that follows one recipe: choose the
slow variables, write $A$, write $\Phi$, minimise $\Phi + \dot A$, then solve and
interpret. Every result is compared with the exact or fully numerical solution of
the same problem. The notebooks share a small Python module, `onsager`, written to
be read as much as to be used.

The chapters assume thermodynamics and fluid mechanics at the level of a senior
undergraduate. No calculus of variations is assumed; functional derivatives are
introduced from the beginning in Chapter 3. Each chapter ends with exercises,
and a sample solution to every exercise is given in the *Solutions* part of the book. Symbols are collected in [Notation and
conventions](notation.md), and books and articles for further study in
[Suggested references](references.md).

## Part I. The principle

| Chapter | Track | Topics |
|---|:-:|---|
| [1. Why a variational principle](../notebooks/part1_principle/01_why_a_variational_principle.ipynb) | core | the Rayleighian of a particle; minimum dissipation in Stokes flow; trial functions for a square duct |
| [2. Onsager's kinetic equations](../notebooks/part1_principle/02_kinetic_equations.ipynb) | core | the friction matrix; reciprocal relations; the energy balance; a settling rod; two interacting spheres |
| [3. Writing the two functions](../notebooks/part1_principle/03_writing_the_two_functions.ipynb) | core | free energies and dissipation functions; functional derivatives; Lagrange multipliers; the squeeze film |
| [4. From the principle to field equations](../notebooks/part1_principle/04_field_equations.ipynb) | core | diffusion and the Einstein relation; the Stokes, Darcy and Brinkman equations |
| [5. The approximation method](../notebooks/part1_principle/05_approximation_method.ipynb) | core | trial profiles for diffusion; a dynamical shape parameter; the Onsager–Machlup integral; what to trust |

## In preparation

Part II, liquids and interfaces: capillary rise and imbibition; spreading of a
droplet; coalescence and sintering of viscous drops; a drying droplet. Part III,
soft materials: diffusion and sedimentation of colloids; swelling and squeezing of a
gel; phase separation; polymer dynamics; nematic liquid crystals.

## Where to start

| If you are… | Read |
|---|---|
| meeting the principle for the first time | Chapters 1 and 2 |
| familiar with the principle and want to use it | Chapters 3 and 5 |
| interested in where the transport equations come from | Chapter 4 |

## Companion books

[*Perturbation methods for differential equations*](https://github.com/ktwyw/math-perturbation-methods),
[*Asymptotic methods: expansions, integrals and sums*](https://github.com/ktwyw/math-asymptotic-methods)
and [*Error analysis and propagation*](https://github.com/ktwyw/error-analysis).

## Licence and citation

Code: MIT. Text and figures: CC BY 4.0. Please cite the material using the
metadata in `CITATION.cff` in the repository.
