# Seats that pass work, not chat

Each seat is a process with a typed handoff. If the handoff fails, the next seat does not run.

1. Ingest — product_id, version, time, bbox, checksum
2. Align — grid_id, time_id, nodata_mask
3. Score — residual, slope, q
4. Gate — pass / fail / not_enough
5. Brief — 80 words, field-locked
6. Audit — who ran it, when, which files

No seat may call a language model to decide significance.
