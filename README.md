# Faunus Notebooks

This repository contains Jupyter Notebooks for running
Faunus, a framework for Metropolis Monte Carlo simulations.

### Description of the notebooks: ###

#### Metropolis Monte Carlo ####

- [pH titration of a protein - explicit, grand canonical salt](titration/)
  (incl. all-atom to amino acid coarse graining)
- [A salt solution in contact with a surface with titrating surface sites](surface-titrating/)
- [Melted salt with off-center charges](offcenter-hardsphere/)
  (incl. parallel tempering of the dielectric constant)
- [PMF of rigid fibril in contact with a charged surface](fibril-surface/)
  (incl. `sbatch` submit script, Debye-Hückel electrostatics)
- [Free energy between two rigid bodies in explicit salt and constant pH](twobody/)
- [Uniformly distribute charges on a sphere using hyperspherical geometry](charges-on-sphere/)

#### Statistical mechanical theory ####

- [Equation of state using generalized van der Waals theory (interactive plot)](gvdw-interactive/)
- [Electric multipole interactions: MC, Keesom etc.](multipole/)
