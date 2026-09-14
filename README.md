Annotates GWAS training-set variants with **CADD** features by doing fast,
indexed lookups into a bgzipped, tabix-indexed whole-genome CADD file. It runs
one worker per chromosome in parallel.

## Input / Output

| | Path | Columns |
|---|---|---|
| **Input** | `training_sets/<EFO>_training.tsv` | `chrom pos rsid effect_allele efo study locus beta p_value z beta_sign pip label pip_target` |
| **Output** | `annotated/<EFO>_training.tsv_annotated.tsv` | same columns **+** `CADD_Ref CADD_Alt` **+** every CADD feature column|
| **Reference** | `whole_genome_SNVs_inclAnno.tsv.gz` | tabix-indexed CADD table laid out as `Chrom Pos Ref Alt <features…>` available here: https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz |

Optionally also writes a `<EFO>_training.tsv_unmatched.tsv` holding the input
rows that found no CADD match (`--unmatched`).
 
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

## Matching modes

- **`alt_only`** (default) — matches on `chrom + pos + ALT`, where `ALT` is the
  training set's `effect_allele`. This is the correct mode here because the
  training set has **no REF column**.
- **`ref_alt`** — matches on `chrom + pos + REF + ALT`; requires a real REF
  column supplied via `--ref-col`.

## How it works

1. **Read CADD metadata once** (`main`): grabs the CADD header to learn the
   feature columns, and reads the set of contig names from the tabix index.
2. **Normalize chromosome names**: builds a `contig_map` so bare names like `1`
   or `X` map onto whatever the index actually uses (`1` vs `chr1`, `MT` vs
   `M`, etc.).
3. **Per input file** (`annotate_file`):
   - Reads all rows and collects **unique** variant keys per chromosome. The
     same variant appears once per study, so de-duplicating avoids repeated
     tabix lookups.
   - Builds one task per chromosome; warns and skips any chromosome absent from
     the index.
   - Fans the tasks out over a multiprocessing `Pool` (one chromosome per
     worker), with a `tqdm` progress bar.
4. **Per chromosome** (`process_chunk`): opens the tabix file, `fetch`es the CADD
   rows at each position, and keeps rows whose `Pos` and `Alt` (and `Ref` in
   `ref_alt` mode) match. If a variant has more than one CADD row it keeps the
   first and increments a multi-hit counter. Feature lists are padded/trimmed to
   a fixed width so the output stays rectangular even with ragged source rows.
5. **Write output** (`annotate_file`): re-joins the annotations back to the
   original lines **in input order**, appends `CADD_Ref`/`CADD_Alt` (unless
   `--no-keep-ref-alt`) plus the feature columns, and routes unmatched rows to
   the unmatched file.
6. **Summary** (`main`): prints per-file and grand-total matched counts,
   percentages, multi-hit counts, and wall-clock time.

## Key functions

- **`get_tabix_header`** — resolves the CADD column header, trying, in order: an
  explicit header file, the tabix embedded header (last header line), then the
  first line of the uncompressed file.
- **`tabix_contigs`** — returns the set of contigs in the tabix index.
- **`process_chunk`** — the parallel worker; annotates one chromosome's unique
  keys and returns `{key: [CADD_Ref, CADD_Alt, *features]}`.
- **`annotate_file`** — orchestrates de-duplication, task creation, the worker
  pool, and ordered output for a single file.
- **`main`** — globs inputs, loads CADD metadata, builds the contig map, loops
  over files, and reports totals.

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