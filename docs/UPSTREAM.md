# Upstream provenance

Code here that was ported from another repository records exactly where it came from, so later
upstream work can be told apart from what is already here.

## The record

[upstream/bikelab.json](../upstream/bikelab.json):

| Field | Meaning |
|---|---|
| `repository` | The upstream repository |
| `commit` | Full SHA of the upstream commit the port matches |
| `watch` | Upstream paths whose new or changed files should be reviewed |
| `ignore` | Patterns under `watch` that are irrelevant here (viewers, styles, vendored code) |
| `ported[].source` | An upstream file that was ported |
| `ported[].blob` | Its git blob hash at `commit`; proves the record names the content that was read |
| `ported[].what` | What was taken from it |
| `ported[].targets` | Files here that carry it |

## Automation

- **CI (`provenance` job):** `tools/upstream_provenance.py check` fails if the commit is missing, a blob
  does not match, or a target file was moved without updating the record.
- **Weekly (`upstream-drift.yml`):** compares the commit with upstream `main` and keeps one issue,
  "Upstream drift: BikeLab", up to date: ported files that changed (with their commits) and new or
  changed watched files that are not ported. The issue closes itself when there is no drift.

## After porting upstream changes

```bash
python tools/upstream_provenance.py drift --upstream-dir ../BikeLab --report drift.md   # what to port
# ...port and test...
python tools/upstream_provenance.py record --upstream-dir ../BikeLab --commit <full-sha>
python tools/upstream_provenance.py check --upstream-dir ../BikeLab
```

Commit the code and the updated record together, naming the upstream commit in the message.
