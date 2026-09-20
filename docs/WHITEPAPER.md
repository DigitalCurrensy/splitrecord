# SPLITRECORD whitepaper (plain language)

## The job
Agencies publish two measurements that people treat as the same story. Often they are not. Vegetation looks healthy while the aquifer falls. Radar sees a structure change that the camera still paints green. A climate model and a satellite point opposite ways in the same county.

The usual fix is to average them. That is how bad decisions get a pretty map.

SPLITRECORD keeps both records. It asks four questions:

1. Are these two products even measuring a related thing?
2. After you put them on the same grid and time, do they move together?
3. If they do not, is the split large enough to survive a significance test?
4. What human decision pretends they agree?

## Method (no mystery)
- Ingest two named NASA or partner products with version IDs.
- Align time and place. Refuse to compare mismatched resolutions without saying so.
- Compute a residual: A minus a scaled B, or rank correlation, depending on the pair.
- Test whether the residual has a trend or a break. Mann-Kendall and Sen slope with an autocorrelation correction. False-discovery control on a map.
- Write an 80-word brief from structured fields only. No model is allowed to invent a cause.

## Output
A ranked list of places. Each card carries: product A, product B, versions, dates, residual, q-value, decision at risk, and "not enough data" when that is the truth.

## Limits
We do not own NASA data. We do not certify court evidence in v1. We do not resolve which sensor is "correct" without a third source.
