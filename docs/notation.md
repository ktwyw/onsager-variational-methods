# Notation and conventions

| Symbol | Meaning |
|---|---|
| $x_i$, $\dot x_i$ | slow variables and their rates of change |
| $A$ | free energy (of the system and of any reservoir it exchanges with) |
| $\dot A$ | rate of change of the free energy, $\sum_i (\partial A/\partial x_i)\,\dot x_i$ |
| $\Phi$ | dissipation function: half the rate at which work is turned into heat |
| $R = \Phi + \dot A$ | Rayleighian; minimised with respect to the rates at fixed state |
| $\zeta_{ij}$ | friction matrix, $\Phi = \tfrac12\sum\zeta_{ij}\dot x_i\dot x_j$; symmetric and positive definite |
| $\mu_{ij}$ | mobility matrix, the inverse of $\zeta_{ij}$ |
| $\zeta$ | friction constant of one particle ($6\pi\eta a$ for a sphere) |
| $\eta$ | viscosity |
| $\xi$ | friction per unit volume between a liquid and a porous solid, $\eta/\kappa$ |
| $\kappa$ | permeability (Chapters 4 and 6); gradient-energy coefficient (Chapter 3) |
| $c$ | number density; $f(c)$ free energy density |
| $\mu$ | chemical potential, $\delta A/\delta c$ (not to be confused with $\mu_{ij}$) |
| $\Pi$ | osmotic pressure, $c f' - f$ |
| $D$ | diffusion coefficient; $k_BT/\zeta$ in a dilute solution |
| $J$, $v$ | flux and mean velocity of particles, $J = c\,v$ |
| $\delta A/\delta c$ | functional derivative |
| $p$ | pressure; the Lagrange multiplier of incompressibility |
| $a_\alpha$ | parameters of a trial solution (Chapter 5) |
| $\gamma$, $\theta$ | surface tension; contact angle (Part II); in Chapter 7 $\theta$ is the dynamic and $\theta_e$ the equilibrium angle |
| $h_J$ | Jurin height, $2\gamma\cos\theta/\rho g a$ (Chapter 6) |
| $\phi$ | porosity (Chapter 6) |
| $\Omega$ | volume of a drop (Chapter 7) |
| $S$ | spreading coefficient, $\gamma_{SV} - \gamma_{SL} - \gamma$ (Chapter 7) |
| $h_c$, $R$ | thickness at the centre and radius of the wetted circle of a drop (Chapter 7) |
| $b$ | microscopic cut-off at a contact line; thickness of a precursor film (Chapter 7) |
| $\ell$ | the logarithm $\ln(h_c/b) - 1$ in the friction of a liquid wedge (Chapter 7); elsewhere a length |
| $\xi_{\rm cl}$, $\xi_{\rm hyd}$ | friction constant of a contact line; its hydrodynamic part, $3\eta\ell/\theta$ (Chapter 7) |

Conventions:

* The dissipation function carries the factor $\tfrac12$: the rate of dissipation is
  $2\Phi$, and at the minimum of the Rayleighian $\dot A = -2\Phi$.
* Numerical work is in dimensionless units, stated in each chapter.
* A result is "checked" when it has been compared with an exact solution or with a
  numerical solution accurate enough for the comparison, and the size of the
  difference is printed in the notebook.
