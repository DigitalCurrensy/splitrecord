# SPLITRECORD

For a hydrologist who already holds two official records of one basin and cannot honestly average them.

**Owner:** Digital Currensy Inc.
**License:** Apache-2.0. Our code only. Cited data and papers stay with their authors.
**Status:** Private until the owner publishes it.

## What it decides

Split, or no split. The residual is the result.

## The rule

Align the two series on shared months. Drop gaps. Do not fill them. On the overlap, the residual is z(A) minus z(B). Mann-Kendall and Sen’s slope score that residual. A short series is not a trend.

## Worked cases

The pairs in this repository are worked months. They prove the rule. They are not a file a customer sent, and they are not an InSAR or GRACE delivery.

## What it will not do

- Average the two records onto one map.
- Fetch granules in order to print the score.
- Sign a permit or send an invoice.

## Run

```
git clone <this repo>
cd splitrecord
PYTHONPATH=src python -m unittest tests.test_kernel
```

Python 3.12. No third-party packages. The test is the demo.
