#!/usr/bin/env python3
"""
Annotate training-set variants with CADD features via tabix.

Input  : training_sets/<EFO>_training.tsv
         columns: chrom pos rsid effect_allele efo study locus
                  beta p_value z beta_sign pip label pip_target
Output : annotated/<EFO>_training.tsv_annotated.tsv
         same columns + CADD_Ref CADD_Alt + all CADD feature columns

Matching
    alt_only (default)  chrom + pos + ALT == effect_allele
    ref_alt             chrom + pos + REF + ALT   (needs a real REF column)

The training set has no REF column, so alt_only is the correct mode.

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
    p.add_argument("--pattern",    default="*_training.tsv",
                   help="only files matching this are annotated")
    p.add_argument("--tabix",  required = True)
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


def process_chunk(task):
    """
    Annotate one chromosome's worth of UNIQUE variant keys.
    Returns {key: [CADD_Ref, CADD_Alt, *features]}
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

        found = []
        for hit in hits:
            hf = hit.split("\t")
            if int(hf[CADD_POS]) != pos:
                continue
            if hf[CADD_ALT].upper() != alt:
                continue
            if match_mode == "ref_alt" and hf[CADD_REF].upper() != ref:
                continue
            found.append(hf)

        if not found:
            continue
        if len(found) > 1:
            multi += 1

        hf    = found[0]
        feats = hf[4:]
        # pad/trim ragged rows so every output line has the same width
        if len(feats) != n_extra:
            feats = (feats + [""] * n_extra)[:n_extra]
        ann[key] = [hf[CADD_REF], hf[CADD_ALT]] + feats

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
        return 0, 0, 0

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
    out_header = list(header)
    if not args.no_keep_ref_alt:
        out_header += ["CADD_Ref", "CADD_Alt"]
    out_header += extra_cols

    matched  = 0
    unm_path = out_file.replace("_annotated.tsv", "_unmatched.tsv")
    fu = open(unm_path, "w") if args.unmatched else None
    if fu:
        fu.write("\t".join(header) + "\n")

    with open(out_file, "w") as fo:
        fo.write("\t".join(out_header) + "\n")
        for ln, key in zip(lines, row_keys):
            rec = ann.get(key)
            if rec is None:
                if fu:
                    fu.write(ln + "\n")
                continue
            fields = ln.split("\t")
            fields += rec if not args.no_keep_ref_alt else rec[2:]
            fo.write("\t".join(fields) + "\n")
            matched += 1

    if fu:
        fu.close()

    return matched, len(lines), total_multi


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

    grand_matched = grand_total = grand_multi = 0
    wall_start = time.time()

    for in_file in input_files:
        base     = os.path.basename(in_file)
        out_file = os.path.join(args.output_dir, base + "_annotated.tsv")
        print(f"\n> {base}", file=sys.stderr)
        t0 = time.time()

        matched, total, multi = annotate_file(
            in_file, out_file, args, contig_map, n_extra, extra_cols)

        elapsed = int(time.time() - t0)
        pct = 100 * matched / total if total else 0
        print(f"  {matched:,}/{total:,} matched ({pct:.1f}%) in {elapsed}s "
              f"-> {out_file}", file=sys.stderr)
        if multi:
            print(f"  {multi:,} sites had >1 CADD row for that ALT (kept first)",
                  file=sys.stderr)
        if args.unmatched and matched < total:
            print(f"  unmatched -> "
                  f"{out_file.replace('_annotated.tsv', '_unmatched.tsv')}",
                  file=sys.stderr)

        grand_matched += matched
        grand_total   += total
        grand_multi   += multi

    elapsed = int(time.time() - wall_start)
    h, r = divmod(elapsed, 3600)
    m, s = divmod(r, 60)
    pct = 100 * grand_matched / grand_total if grand_total else 0
    print(f"\nAll done in {h}h {m}m {s}s. "
          f"{grand_matched:,}/{grand_total:,} matched ({pct:.1f}%).",
          file=sys.stderr)
    if grand_multi:
        print(f"{grand_multi:,} multi-hit sites overall.", file=sys.stderr)


if __name__ == "__main__":
    main()