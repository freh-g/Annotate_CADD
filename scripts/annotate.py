#!/usr/bin/env python3
"""
Annotate training-set variants with CADD features via tabix.

Input  : training_sets/<EFO>_training.tsv
         columns: chrom pos rsid effect_allele efo study locus
                  beta p_value z beta_sign pip label pip_target
Output : annotated/<EFO>_training.tsv_annotated.tsv
         same columns + CADD_hit_index + CADD_Ref CADD_Alt + all CADD
         feature columns

Matching
    alt_only (default)  chrom + pos + ALT == effect_allele
    ref_alt             chrom + pos + REF + ALT   (needs a real REF column)

The training set has no REF column, so alt_only is the correct mode.

MULTI-HIT VARIANTS -- KEEPS ALL ROWS (does not take-first)
    CADD reports one row per overlapping transcript, so ~32% of variants in
    this pipeline return more than one row for the same chrom/pos/ALT. The
    site-level features (conservation, chromatin, GERP, etc.) are identical
    across those rows; the transcript-specific ones (Consequence, GeneID,
    Exon, cDNApos, oAA/nAA, etc.) genuinely differ per transcript and are
    NOT interchangeable -- keeping only the first, as the previous version
    of this script did, silently discards whichever other transcripts also
    overlapped that position (e.g. a variant missense in one transcript but
    intronic in another).

    This version instead EXPANDS each matched input row into one output row
    PER matching CADD transcript hit, all carrying the same original data
    (chrom, pos, beta, etc.) but each with its own CADD_Ref/CADD_Alt/feature
    columns. A new `CADD_hit_index` column (0-based) is added so rows
    originating from the same input variant can be grouped back together
    downstream (e.g. `groupby(['chrom','pos','effect_allele','study'])`).

    This means: (a) matched/total in the summary now refers to INPUT rows,
    not output rows -- a single matched input row can produce several output
    rows; (b) any downstream script doing per-row analysis (recall@k,
    region-grouped CV, etc.) must decide how to handle the resulting
    duplicates -- e.g. collapse back to one row per variant (as the
    previous take-first behaviour did) at the point where it matters,
    or explicitly treat transcript-level annotation as its own axis.

Usage
    python 11_annotate_parallel.py \
        --input-dir  ../training_sets/ \
        --output-dir ../annotated/ \
        --pattern '*_training.tsv' \
        --unmatched
"""
import argparse, sys, os, time, glob
import pysam
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir",  default="../training_sets/")
    p.add_argument("--pattern",    default="*training*",
                   help="only files matching this are annotated")
    p.add_argument("--tabix",  required=True,help="path of the CADD tabix file")
    p.add_argument("--output-dir", default="../annotated/")
    p.add_argument("--tabix-header", default=None)
    p.add_argument("--match-mode", choices=["alt_only", "ref_alt"],
                   default="alt_only")
    p.add_argument("--chrom-col", type=int, default=0)
    p.add_argument("--pos-col",   type=int, default=1)
    p.add_argument("--alt-col",   type=int, default=3,
                   help="effect_allele column (0-based)")
    p.add_argument("--ref-col",   type=int, default=-1,
                   help="REF column; only used with --match-mode ref_alt")
    p.add_argument("--no-keep-ref-alt", action="store_true",
                   help="do not append the CADD Ref/Alt actually matched")
    p.add_argument("--unmatched", action="store_true",
                   help="write a *_unmatched.tsv per input file")
    p.add_argument("--threads",   type=int, default=cpu_count())
    return p.parse_args()


def get_tabix_header(tabix_path, tabix_header_path=None):
    if tabix_header_path:
        with open(tabix_header_path) as f:
            return f.readline().rstrip("\n").lstrip("#").split("\t")
    tb = pysam.TabixFile(tabix_path)
    if tb.header:
        line = list(tb.header)[-1]
        tb.close()
        return line.lstrip("#").split("\t")
    tb.close()
    plain = tabix_path[:-3] if tabix_path.endswith(".gz") else tabix_path
    with open(plain) as f:
        return f.readline().rstrip("\n").lstrip("#").split("\t")


def tabix_contigs(tabix_path):
    tb = pysam.TabixFile(tabix_path)
    c = set(tb.contigs)
    tb.close()
    return c


# CADD tsv layout is fixed:  Chrom  Pos  Ref  Alt  <features...>
CADD_CHROM, CADD_POS, CADD_REF, CADD_ALT = 0, 1, 2, 3


def match_hits(hits, pos, alt, ref, match_mode):
    """
    Pure matching logic, factored out of process_chunk so it can be unit
    tested without a real tabix file. Returns a list of matched CADD rows
    (each a list of fields split on tab) -- ALL matches, not just the first.
    """
    found = []
    for hit in hits:
        hf = hit.split("\t") if isinstance(hit, str) else hit
        if int(hf[CADD_POS]) != pos:
            continue
        if hf[CADD_ALT].upper() != alt:
            continue
        if match_mode == "ref_alt" and hf[CADD_REF].upper() != ref:
            continue
        found.append(hf)
    return found


def process_chunk(task):
    """
    Annotate one chromosome's worth of UNIQUE variant keys.
    Returns {key: [[CADD_Ref, CADD_Alt, *features], ...]}  -- a LIST of
    matches per key, not a single row.
    """
    (keys, chrom_tag, contig, tabix_path, tabix_header_path,
     match_mode) = task

    tb      = pysam.TabixFile(tabix_path)
    n_extra = len(get_tabix_header(tabix_path, tabix_header_path)) - 4

    ann   = {}
    multi = 0

    for key in keys:
        if match_mode == "alt_only":
            _, pos, alt = key
            ref = None
        else:
            _, pos, ref, alt = key

        try:
            hits = tb.fetch(contig, pos - 1, pos)
        except ValueError:
            continue

        found = match_hits(hits, pos, alt, ref, match_mode)
        if not found:
            continue
        if len(found) > 1:
            multi += 1

        rows = []
        for hf in found:
            feats = hf[4:]
            # pad/trim ragged rows so every output line has the same width
            if len(feats) != n_extra:
                feats = (feats + [""] * n_extra)[:n_extra]
            rows.append([hf[CADD_REF], hf[CADD_ALT]] + feats)
        ann[key] = rows

    tb.close()
    return ann, multi, chrom_tag


def annotate_file(in_file, out_file, args, contig_map, n_extra, extra_cols):
    cc, pc, ac, rc = (args.chrom_col, args.pos_col,
                      args.alt_col, args.ref_col)
    alt_only = args.match_mode == "alt_only"

    with open(in_file) as fin:
        header = fin.readline().rstrip("\n").split("\t")
        lines  = [ln.rstrip("\n") for ln in fin if ln.strip()]

    if not lines:
        print("  empty file, skipping", file=sys.stderr)
        return 0, 0, 0, 0

    # ── collect UNIQUE variant keys per chromosome ───────────────────────────
    # the same variant appears once per study, so this avoids repeat lookups
    per_chrom = defaultdict(set)
    row_keys  = []
    for ln in lines:
        f     = ln.split("\t")
        chrom = f[cc].replace("chr", "")
        pos   = int(f[pc])
        alt   = f[ac].upper()
        key   = (chrom, pos, alt) if alt_only else \
                (chrom, pos, f[rc].upper(), alt)
        row_keys.append(key)
        per_chrom[chrom].add(key)

    chroms = sorted(per_chrom, key=lambda c: (len(c), c))
    n_uniq = sum(len(v) for v in per_chrom.values())
    print(f"  {len(lines):,} rows -> {n_uniq:,} unique variants "
          f"across {len(chroms)} chromosomes", file=sys.stderr)

    tasks = []
    for chrom in chroms:
        contig = contig_map.get(chrom)
        if contig is None:
            print(f"  WARNING: contig '{chrom}' absent from the tabix index "
                  f"- {len(per_chrom[chrom]):,} variants cannot be annotated",
                  file=sys.stderr)
            continue
        tasks.append((sorted(per_chrom[chrom]), chrom, contig,
                      args.tabix, args.tabix_header, args.match_mode))

    ann         = {}
    total_multi = 0
    if tasks:
        ncpus = max(1, min(args.threads, len(tasks)))
        with Pool(processes=ncpus) as pool:
            for sub, multi, _ in tqdm(
                pool.imap_unordered(process_chunk, tasks),
                total=len(tasks), desc=f"  {os.path.basename(in_file)}",
                unit="chr", file=sys.stderr
            ):
                ann.update(sub)
                total_multi += multi

    # ── write annotated + unmatched, preserving input row order ──────────────
    # one output row PER matching CADD transcript hit -- CADD_hit_index lets
    # you regroup back to the original variant downstream.
    out_header = list(header) + ["CADD_hit_index"]
    if not args.no_keep_ref_alt:
        out_header += ["CADD_Ref", "CADD_Alt"]
    out_header += extra_cols

    matched_rows   = 0   # INPUT rows that matched at least one CADD hit
    output_rows    = 0   # total OUTPUT lines written (>= matched_rows)
    unm_path = out_file.replace("_annotated.tsv", "_unmatched.tsv")
    fu = open(unm_path, "w") if args.unmatched else None
    if fu:
        fu.write("\t".join(header) + "\n")

    with open(out_file, "w") as fo:
        fo.write("\t".join(out_header) + "\n")
        for ln, key in zip(lines, row_keys):
            hit_rows = ann.get(key)
            if not hit_rows:
                if fu:
                    fu.write(ln + "\n")
                continue
            matched_rows += 1
            base_fields = ln.split("\t")
            for hit_idx, rec in enumerate(hit_rows):
                fields = list(base_fields) + [str(hit_idx)]
                fields += rec if not args.no_keep_ref_alt else rec[2:]
                fo.write("\t".join(fields) + "\n")
                output_rows += 1

    if fu:
        fu.close()

    return matched_rows, len(lines), total_multi, output_rows


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if args.match_mode == "ref_alt" and args.ref_col < 0:
        sys.exit("ERROR: --match-mode ref_alt requires --ref-col")

    input_files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    input_files = [f for f in input_files if os.path.isfile(f)]
    if not input_files:
        sys.exit(f"No files matching '{args.pattern}' in {args.input_dir}")

    # tabix metadata, read once
    tabix_cols = get_tabix_header(args.tabix, args.tabix_header)
    extra_cols = tabix_cols[4:]
    n_extra    = len(extra_cols)

    # map bare chromosome names onto whatever the index actually uses
    contigs    = tabix_contigs(args.tabix)
    contig_map = {}
    for c in list(map(str, range(1, 23))) + ["X", "Y", "MT", "M"]:
        for cand in (c, f"chr{c}"):
            if cand in contigs:
                contig_map[c] = cand
                break

    print(f"Found {len(input_files)} file(s) matching '{args.pattern}'",
          file=sys.stderr)
    print(f"Match mode: {args.match_mode}   CADD feature columns: {n_extra}",
          file=sys.stderr)
    print(f"Tabix contigs look like: {sorted(contigs)[:3]} ...", file=sys.stderr)
    print("NOTE: this version keeps ALL matching CADD rows per variant "
          "(one output row per transcript hit), not just the first. "
          "See CADD_hit_index in the output.", file=sys.stderr)

    grand_matched = grand_total = grand_multi = grand_output = 0
    wall_start = time.time()

    for in_file in input_files:
        base     = os.path.basename(in_file)
        out_file = os.path.join(args.output_dir, base + "_annotated.tsv")
        print(f"\n> {base}", file=sys.stderr)
        t0 = time.time()

        matched, total, multi, out_rows = annotate_file(
            in_file, out_file, args, contig_map, n_extra, extra_cols)

        elapsed = int(time.time() - t0)
        pct = 100 * matched / total if total else 0
        print(f"  {matched:,}/{total:,} input rows matched ({pct:.1f}%) "
              f"-> {out_rows:,} output rows in {elapsed}s -> {out_file}",
              file=sys.stderr)
        if multi:
            print(f"  {multi:,} variants had >1 CADD row (ALL kept, "
                  f"expanded into separate output rows)", file=sys.stderr)
        if args.unmatched and matched < total:
            print(f"  unmatched -> "
                  f"{out_file.replace('_annotated.tsv', '_unmatched.tsv')}",
                  file=sys.stderr)

        grand_matched += matched
        grand_total   += total
        grand_multi   += multi
        grand_output  += out_rows

    elapsed = int(time.time() - wall_start)
    h, r = divmod(elapsed, 3600)
    m, s = divmod(r, 60)
    pct = 100 * grand_matched / grand_total if grand_total else 0
    print(f"\nAll done in {h}h {m}m {s}s. "
          f"{grand_matched:,}/{grand_total:,} input rows matched ({pct:.1f}%), "
          f"{grand_output:,} total output rows.", file=sys.stderr)
    if grand_multi:
        print(f"{grand_multi:,} multi-hit variants overall "
              f"(expanded to {grand_output - grand_matched:,} extra rows).",
              file=sys.stderr)


if __name__ == "__main__":
    main()