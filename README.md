# MATRIX_MultiphaseTissueModelling

## What
Finite element simulation of two-phase models for general tissue-like structures. 

## Context 
Group discussions on multiphase tissue modelling as part of the [MATRIX program on the Mathematics of Tissue Dynamics](https://www.matrix-inst.org.au/events/mathematics-of-tissue-dynamics)

Group participants:
  - Adriana Zanca
  - Ishraq Ahmed
  - Joshua Won
  - Alexander Browning
  - Jennifer Flegg
  - Ricardo Ruiz Baier

## Steps 
Rewriting the Eulerian 1D formulation from 

```
@rticle{Breward2002,
  author  = {Breward C and others},
  title   = {The role of cell-cell interactions in a two-phase model for avascular tumour growth},
  year    = {2002},
  journal = {Journal of Mathematical Biology},
  volume  = {45},
  pages   = {125--152},
  doi     = {10.1007/s002850200149},
}
```

into multidimensional Lagrangian form. 

Then writing a weak formulation and a finite element discretisation in terms of cell displacement, fluid pressure, solid volume fraction, and oxygen concentration. 
