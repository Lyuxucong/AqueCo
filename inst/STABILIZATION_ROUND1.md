# AqueCo stabilization — round 1

Changes:
- Replaced the source-file path inspection in `test-beta-decomp.R` with a behavioral plotting test.
- Replaced `linewidth` with `size` for `linkET::geom_square` compatibility.
- Qualified the `Error` aesthetic in `aq_env_plot` with `.data$Error`.
- Kept the existing modules and omitted LMM hierarchical partitioning.

Verification:
- ZIP structure verified.
- R runtime is not available in the packaging environment; `devtools::test()` and `devtools::check()` still need to run locally.
- Remaining `linkET::aes_()` deprecation warnings originate in the upstream package and may remain.
- Rd warnings/notes need to be re-evaluated from the fresh `devtools::check()` output.
