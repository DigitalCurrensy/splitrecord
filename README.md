# SPLITRECORD

SPLITRECORD scores the disagreement between two records of the same place. It does not average them into one map.

**Owner:** Digital Currensy Inc.
**License:** Apache-2.0. Our code only. Cited data and papers stay with their authors.

## What it decides

Split, or no split. The residual is the result.

## The rule

Align the two series on shared months. Drop gaps. Do not fill them. On the overlap, the residual is z(A) minus z(B). Mann-Kendall and Sen’s slope score that residual. Benjamini–Hochberg is applied across units. A unit splits only when the overlap is long enough, the corrected trend clears the bar, the residual rises, and the second series falls.

## Worked cases

Ten hydrologic units in the southern San Joaquin. Each case is a monthly pair stored in this repository. The pair proves the rule. It is not a file a customer sent, and it is not an InSAR or GRACE delivery.

## What it will not do

- Average the two records onto one map.
- Treat a short series as a trend.
- Fetch granules in order to print the score.
- Sign a permit or send an invoice.

## Run

```
PYTHONPATH=src python -m unittest tests.test_kernel
```

Notes under `docs/` are the build record. This page is the description.
