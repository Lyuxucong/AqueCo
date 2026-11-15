# AqueCo performance development build

## RDA/CCA
- `select="none"` skips forward selection and fits all chosen environmental variables.
- `tests="global"` skips term and axis permutation tests.
- `mode="fast"` changes the default 999 permutations to 99, but explicit nondefault counts are preserved.
- For reproducible scientific inference use publication mode and the full intended tests.

## Null models
- `aq_null(..., cache=TRUE)` caches the **filtered** phylogenetic distance matrix.
- `aq_beta_nri`, `aq_beta_nti`, and `aq_rcbray` accept `workers`.
- Unix/macOS: multicore parallel chunk execution. Windows: warning and sequential fallback.
- Each iteration uses deterministic per-iteration seeds when seed is supplied.
- Random streams differ from previous AqueCo releases even with the same seed.
- Algorithm definitions and filtering criteria are unchanged.
- Caching increases memory usage; large unfiltered phylogenies may not fit.
- The RCbray randomization is an AqueCo-specific model, not validated for exact microeco equivalence.

## Verification
The source and test files are included, but no R runtime is available in this environment.
Run `devtools::test(filter="performance")` and `devtools::check()` locally.
