# SparseAttack (PlatEMO experiments)

Brief project README for the SparseAttack repository. This project contains a collection
of PlatEMO-based problem definitions, data loaders and experiment harnesses for
research into sparse attacks on optimization problems.
## Requirements

- MATLAB (R2018a or newer recommended). This repository uses plain `.m` files and
  PlatEMO code.

## Layout (important files & folders)
- `main.m` — experiment entry point. It prepares the save folder and drives runs.
- `utils/Origin.m` — computes and caches baseline (original) performance values.
- `utils/Perturbation.m` — central attacker factory (maps problem name to a problem-specific perturbation/evaluation helper).
- `utils/computeOriginObjSingle.m`, `utils/computeOriginObjMulti.m` — canonical origin-evaluator helpers (single-/multi-objective) used across the codebase.
- `utils/LoadOrgData.m` — data loader and per-problem loader helpers.
- `PlatEMO 4.0/Problems/...` — problem definitions (AttackTSP, AttackKP, AttackMOTSP, AttackMOKP, AttackPO, AttackCase118, etc.).
- `PlatEMO 4.0/Algorithms/A_attack...` — attack implementations and helpers (algorithms
  or operators that generate/simulate adversarial perturbations). Place attack-specific
  code here.
- `PlatEMO 4.0/Algorithms/A_attacked...` — attacked algorithm implementations (the
  algorithms under test).

## How to run a quick smoke test
- Run main directly in the MATLAB command window (ensure the project folder is on the MATLAB path or start MATLAB from the project root).
- Or invoke from MATLAB’s command line with the following example (adjust parameters as needed):

```matlab
main(1,'TSP','GA','SparseEA',2500,3000,30,30,0.05,3,0.1,10,1)
```
