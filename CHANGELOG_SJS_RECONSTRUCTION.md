# SJS reconstruction changelog

## 2026-08-17

- Reframed the manuscript around high-dimensional profile inference rather than screening and block-coordinate optimisation as separate contributions.
- Defined an unsmoothed, ridge-profile population target.
- Replaced the legacy residualised-design score by the exact envelope-theorem score using the original fixed-effect design.
- Added the missing ridge-curvature term to the residualised representation of the Schur-complement Hessian.
- Replaced full-matrix precision inversion by coordinate-wise Dantzig rows in the primary specification.
- Replaced observation-scale variance notation by an explicit independent-cluster sandwich variance.
- Removed empirical coverage calibration from the simulation design.
- Replaced the GSE121239 illustration by a preregistered GSE65391 longitudinal SLE analysis with an explicit outcome audit.
- Added a full modular Monte Carlo protocol, benchmark matrix, output contract and reproducibility scripts.
- Preserved the pre-reconstruction repository state on `legacy-pre-sjs-reconstruction`.
