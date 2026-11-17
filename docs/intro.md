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

The level is that of a graduate course in physics or engineering, and the chapters
assume thermodynamics and fluid mechanics at the level of a senior undergraduate.
No calculus of variations is assumed; functional derivatives are introduced from the
beginning in Chapter 3. Each chapter ends with four exercises, and a sample solution
to every exercise is given in the *Solutions* part of the book. Symbols are collected
in [Notation and conventions](notation.md), and books and articles for further study
in [Suggested references](references.md).

## Part I. The principle

| Chapter | Track | Topics |
|---|:-:|---|
| [1. Why a variational principle](../notebooks/part1_principle/01_why_a_variational_principle.ipynb) | core | the Rayleighian of a particle; minimum dissipation in Stokes flow; trial functions for a square duct |
| [2. Onsager's kinetic equations](../notebooks/part1_principle/02_kinetic_equations.ipynb) | core | the friction matrix; reciprocal relations; the energy balance; a settling rod; two interacting spheres |
| [3. Writing the two functions](../notebooks/part1_principle/03_writing_the_two_functions.ipynb) | core | free energies and dissipation functions; functional derivatives; Lagrange multipliers; the squeeze film |
| [4. From the principle to field equations](../notebooks/part1_principle/04_field_equations.ipynb) | core | diffusion and the Einstein relation; the Stokes, Darcy and Brinkman equations |
| [5. The approximation method](../notebooks/part1_principle/05_approximation_method.ipynb) | core | trial profiles for diffusion; a dynamical shape parameter; the Onsager–Machlup integral; what to trust |

## Part II. Liquids and interfaces

| Chapter | Track | Topics |
|---|:-:|---|
| [6. Capillary rise and imbibition](../notebooks/part2_liquids_interfaces/06_capillary_rise.ipynb) | core | the Lucas–Washburn law; Jurin's height; inertia and oscillations; imbibition into porous materials |
| [7. Spreading of a droplet](../notebooks/part2_liquids_interfaces/07_spreading_droplet.ipynb) | core | the friction of a liquid wedge; Tanner's law; partial wetting; the thin-film equation; what a trial shape cannot give |
| [8. Coalescence and sintering of viscous drops](../notebooks/part2_liquids_interfaces/08_coalescence_sintering.ipynb) | core | Frenkel's model; trial shapes and trial flows that must fit; the minimum principle as a bound; Hopper's exact solution |
| [9. A drying droplet](../notebooks/part2_liquids_interfaces/09_drying_droplet.ipynb) | core | flow driven by evaporation; contact-line friction; coffee rings, volcanoes and mountains |

## Part III. Soft materials

| Chapter | Track | Topics |
|---|:-:|---|
| [10. Diffusion and sedimentation of colloids](../notebooks/part3_soft_materials/10_diffusion_sedimentation.ipynb) | core | osmotic pressure and friction; collective diffusion and hindered settling; sedimentation equilibrium; trial solutions that must have the right physics |
| [11. Swelling and squeezing of a gel](../notebooks/part3_soft_materials/11_gel_swelling.ipynb) | core | stress and diffusion in a gel; swelling and squeezing of a sheet; a swelling sphere and the shear modulus; trial functions as bounds |
| [12. Phase separation](../notebooks/part3_soft_materials/12_phase_separation.ipynb) | core | the Cahn–Hilliard equation from the Rayleighian; phase diagram and interface; spinodal decomposition; coarsening |
| [13. Polymer dynamics](../notebooks/part3_soft_materials/13_polymer_dynamics.ipynb) | core | the Rouse model; normal modes; fluctuations and Onsager's regression hypothesis; the dumbbell as a trial; coil–stretch transition |
| [14. Nematic liquid crystals](../notebooks/part3_soft_materials/14_nematic_liquid_crystals.ipynb) | core | Frank elasticity and rotational friction; the Fréedericksz transition; switching; flow alignment, tumbling and the Parodi relation |

## Still to come

Part IV, on dissipation functions that are not quadratic and on error estimates for
trial solutions, will be added after the corresponding research papers are published.

## Where to start

| If you are… | Read |
|---|---|
| meeting the principle for the first time | Chapters 1 and 2 |
| familiar with the principle and want to use it | Chapters 3 and 5 |
| interested in where the transport equations come from | Chapter 4 |
| looking for a complete application | Chapters 6 and 7 |
| interested in where trial solutions go wrong | Chapters 5, 7, 8 and 10 |
| working with colloids or suspensions | Chapters 4 and 10 |
| working with gels | Chapters 10 and 11 |
| working with mixtures and phase separation | Chapters 4 and 12 |
| working with polymers | Chapters 2 and 13 |
| working with liquid crystals | Chapters 2 and 14 |

## Companion books

[*Perturbation methods for differential equations*](https://github.com/ktwyw/math-perturbation-methods),
[*Asymptotic methods: expansions, integrals and sums*](https://github.com/ktwyw/math-asymptotic-methods)
and [*Error analysis and propagation*](https://github.com/ktwyw/error-analysis).

## Licence and citation

Code: MIT. Text and figures: CC BY 4.0. Please cite the material using the
metadata in `CITATION.cff` in the repository.
