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
| `driftOwner` | Optional, default true; on the manifest or a `ported` entry (entry wins). Only owned entries raise drift; see below |
| `ported[].what` | What was taken from it |
| `ported[].targets` | Files here that carry it |

## Shared tooling and ownership

The script and the weekly workflow live in Sub0Math (`tools/provenance/`,
`.github/workflows/upstream-drift-reusable.yml`), pinned by commit in `ci.yml` and `upstream-drift.yml`;
`<Sub0Math>` above is a checkout of it (CPM places it in your build tree). `check` always covers every
entry. An upstream file ported into several repositories is owned by exactly one (`driftOwner`), so a
change to it opens one issue, there; a non-owned entry carries `"driftOwner": false` and `"ownedBy"`.
Here `driftOwner` is false at manifest level (unported watched files are reported by the other
repository) and true on the six files ported here, so their drift issue is raised in this repository. Two further entries (`footprint.js`, `wheel-limits.js`) record the generic geometry and ring solver taken from those files; they are `driftOwner: false`, because the station they belong to is ported elsewhere, and that repository raises their drift.

## Automation

- **CI (`provenance` job):** `upstream_provenance.py check` fails if the commit is missing, a blob
  does not match, or a target file was moved without updating the record.
- **Weekly (`upstream-drift.yml`):** compares the commit with upstream `main` and keeps one issue,
  "Upstream drift: BikeLab", up to date: ported files that changed (with their commits) and new or
  changed watched files that are not ported. The issue closes itself when there is no drift.

## After porting upstream changes

```bash
python <Sub0Math>/tools/provenance/upstream_provenance.py drift --upstream-dir ../BikeLab --report drift.md   # what to port
# ...port and test...
python <Sub0Math>/tools/provenance/upstream_provenance.py record --upstream-dir ../BikeLab --commit <full-sha>
python <Sub0Math>/tools/provenance/upstream_provenance.py check --upstream-dir ../BikeLab
```

Commit the code and the updated record together, naming the upstream commit in the message.
