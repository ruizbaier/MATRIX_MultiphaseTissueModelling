# MATRIX_MultiphaseTissueModelling

## What
Finite element simulation of two-phase models for general tissue-like structures. 

## Context 
Group discussions on multiphase tissue modelling as part of the [MATRIX program on the Mathematics of Tissue Dynamics](https://www.matrix-inst.org.au/events/mathematics-of-tissue-dynamics)

![](https://github.com/ruizbaier/MATRIX_MultiphaseTissueModelling/blob/main/tissue_sketch.png)
![](https://github.com/ruizbaier/MATRIX_MultiphaseTissueModelling/blob/main/glioma_test.png)



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
  title   = {The role of cell-cell interactions in a two-phase model for avascular tumour growth},
  author  = {Breward, CJW and Byrne, HM and Lewis, CE},
  number  = {2},
  year    = {2002},
  journal = {Journal of Mathematical Biology},
  volume  = {45},
  pages   = {125--152},
  doi     = {10.1007/s002850200149},
}
```

into multidimensional Lagrangian form. 

Then writing a weak formulation and a finite element discretisation in terms of cell displacement, fluid pressure, solid volume fraction, and oxygen concentration. 
