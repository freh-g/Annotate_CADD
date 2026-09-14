# CADD Annotation Pipeline

## `annotate.py` — annotate with CADD

Matches each variant against CADD's `whole_genome_SNVs_inclAnno.tsv.gz` (available at https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz) via
tabix, on `chrom + pos + ALT`.

```bash
python 11_annotate_parallel.py \
    --input-dir  ../training_sets/ \
    --output-dir ../annotated/ \
    --pattern '*_training.tsv' \
    --unmatched
```

## Input column requirements

The script locates input columns **by position, not by name**. The header row
is read but only carried through to the output — it is never parsed to find
columns. So the labels in your header can be anything; what matters is that the
right value sits at the index each flag points to.

With the default flags (`--chrom-col 0`, `--pos-col 1`, `--alt-col 3`), the
layout must be:

| Index | Role | Notes |
|---|---|---|
| **0** | chrom | any `chr` prefix is stripped automatically |
| **1** | pos | integer position |
| **2** | *(ignored)* | a placeholder must exist so index 3 is reachable, but its content is never read — despite being `rsid` in the example, it can hold anything |
| **3** | effect allele | the ALT used for the CADD match |

Consequences:

- **Names are irrelevant.** `effect_allele` is just the example label; rename any
  column freely.
- **Column 2 is not required to be an rsid** — it only needs to exist so the
  effect allele lands at position 3.
- **Reordering is fine if you update the flags.** For example, if the effect
  allele sits at position 2, pass `--alt-col 2` (then you'd only need 3 columns).
- **A REF column is not read at all in the default `alt_only` mode**, wherever it
  is and whatever it's called. It only matters under `--match-mode ref_alt`, and
  then you must point at it with `--ref-col <index>`. Renaming alone does
  nothing; position is what's used.
- **Everything after the matched columns**, is never inspected — it's carried straight through to the output.


| flag | default | notes |
|---|---|---|
| `--match-mode` | `alt_only` | use `ref_alt` only if your input has a real REF column |
| `--unmatched` | off | also writes `*_unmatched.tsv` for variants with no CADD hit |

## Command-line arguments

| Flag | Default | Purpose |
|---|---|---|
| `--input-dir` | `../training_sets/` | directory of input files |
| `--pattern` | `*_training.tsv` | glob deciding which files get annotated |
| `--tabix` | (CADD path) | bgzipped, tabix-indexed CADD file |
| `--output-dir` | `../annotated/` | where annotated files go |
| `--tabix-header` | `None` | optional external header file for the CADD table |
| `--match-mode` | `alt_only` | `alt_only` or `ref_alt` |
| `--chrom-col` / `--pos-col` / `--alt-col` | `0` / `1` / `3` | input column indices (0-based); `--alt-col` is the effect allele |
| `--ref-col` | `-1` | REF column, only used with `ref_alt` |
| `--no-keep-ref-alt` | off | omit the matched CADD Ref/Alt from output |
| `--unmatched` | off | write a `*_unmatched.tsv` per input file |
| `--threads` | all CPUs | worker-pool size |


## Example usage

```bash
python 11_annotate_parallel.py \
    --input-dir  folder/to/the/files/to/annotate \
    --tabix /path/to/CADD/annotation/file \ 
    --output-dir ../annotated/ \
    --pattern '*_training.tsv' \
    --unmatched
```


**Output:** `<file>_annotated.tsv` — original columns + `CADD_hit_index` +
`CADD_Ref`/`CADD_Alt` + all CADD feature columns. Row count ≥ input row count.


