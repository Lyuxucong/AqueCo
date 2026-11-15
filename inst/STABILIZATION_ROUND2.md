# AqueCo documentation and compatibility pass (round 2)

Based on AqueCo stabilization round 1 and the supplied R CMD check log.

Changes:
- Synchronize Rd usage signatures with corresponding exported function definitions.
- Add missing Rd argument entries and descriptions for documented functions.
- Correct invalid beta decomposition Rd cross-reference and wrap network usage.
- Replace non-ASCII symbols in source strings with Unicode escape sequences.
- Qualify stats::ave and utils::head/tail.
- Declare data-mask ggplot aesthetics for static checking.
- Remove unused methods Imports and move development notes into inst/.
- Keep the LMM module excluded.

Known limitations:
- `linkET::aes_()` deprecation warning arises in linkET internals and is not fixed by these source changes.
- Rd placeholder descriptions should be enhanced with scientifically detailed parameter definitions before public release.
- R runtime unavailable here: run devtools::test() and devtools::check() locally.
