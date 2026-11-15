# Statistical Interpretation and Reproducibility

[← Table of contents](README.md)

## Record the analysis context

For every analysis, record sample counts, experimental units, grouping variables, filtering criteria, taxonomic rank, abundance scale, transformation, distance metric, statistical test, permutation count, P-value adjustment, and software versions.

## Randomness

Use a fixed random seed where supported:

```r
set.seed(123)
capture.output(sessionInfo(), file = "sessionInfo.txt")
```

If a function exposes its own `seed` argument, specify it explicitly.

## Common interpretation pitfalls

- A significant P value does not measure effect magnitude.
- Association and network analyses do not demonstrate causation.
- Distance matrices create nonindependent pairwise observations.
- Model R-squared does not quantify the causal contribution of an ecological process.
- Shared variation fractions are not independent effects.
- Null-model results depend on the randomization algorithm and filtering.
- Six samples per treatment are not sufficient for all ecological models.

## Writing a Methods section

A reproducible Methods section should state:

1. The software and exact version.
2. Input data and preprocessing.
3. Model selection and transformations.
4. Thresholds, permutations, and random seeds.
5. Multiple-testing correction and uncertainty reporting.
6. Model assumptions and limitations.

## Figure and table provenance

Save both the fitted R object and the exported image. Include a README for result directories describing which script created each output.

## Suggested project structure

```text
project/
├── data/
├── scripts/
├── results/
├── figures/
└── sessionInfo.txt
```

## Software citation

Author: Xucong Lyu. Package: AqueCo. Cite the exact GitHub release or commit used in the analysis; add a repository URL and persistent identifier after publication.
