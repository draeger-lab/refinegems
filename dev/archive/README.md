# Development Archive

This directory contains completed or paused development experiments that are
kept for reproducibility but are no longer active workspaces.

## Database identifier gap filling (archived 19.08.2026)

`database-gapfilling-identifiers.zip` archives the completed `substance2db`
identifier-mapping and three-column primary-key migration work. It contains:

- the scoped input and manually reviewed mapping tables;
- the ModelSEED and VMH matching scripts;
- the reusable migration script and focused tests;
- audit and unresolved-identifier reports;
- the VMH model used for local identifier matching; and
- the experiment README with the full workflow and results.

The archive is retained as a reproducibility record. The current source of
truth is `src/refinegems/data/database/data.db`, with the reproducible media
schema in `src/refinegems/data/database/media_db.sql`.

To inspect the archived files without modifying the archive:

```bash
unzip -l dev/archive/database-gapfilling-identifiers.zip
```

