# SPLITRECORD wiring

```
Earthdata pair A, pair B
        |
   [Ingest seat]   named product + version + checksum
        |
   [Align seat]    same grid, same time window, or abort
        |
   [Score seat]    residual + MK/Sen + FDR
        |
   [Gate seat]     refuse if n too small or versions missing
        |
   [Brief seat]    80-word paragraph from fields only
        |
   [Audit seat]    log every number back to a granule ID
```

Store: Parquet + a tiny Postgres of product versions.
Compute: CPU. No rented cluster for v1.
Interface: one page, one basin, ranked cards.
