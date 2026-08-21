#!/usr/bin/env python3
"""gene_calls.py — split a CNVkit genemetrics table into amplified / deleted
gene lists using ABSOLUTE copy number (cn), thresholded relative to a given
ploidy, instead of a fixed log2-ratio cutoff.

Used by Snakefile rule `genelist`. Requires the genemetrics table to have
been generated with segments from `cnvkit call` (so a `cn` column is
present) and disabled log2/probe pre-filtering (`genemetrics -t 0 -m 1`).

    amp_cutoff = round(ploidy) + amp_offset   -> amplified if cn >= amp_cutoff
    del_cutoff = max(round(ploidy) - del_offset, 0) -> deleted if cn <= del_cutoff
"""

import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genemetrics", required=True, help="input genemetrics TSV")
    p.add_argument("--ploidy", required=True, type=float)
    p.add_argument("--amp-offset", required=True, type=int)
    p.add_argument("--del-offset", required=True, type=int)
    p.add_argument("--min-probes", required=True, type=int)
    p.add_argument("--amp-out", required=True)
    p.add_argument("--del-out", required=True)
    p.add_argument("--amp-genes-out", required=True)
    p.add_argument("--del-genes-out", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.genemetrics, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    if "cn" not in fieldnames:
        sys.exit(
            f"ERROR: {args.genemetrics!r} has no 'cn' column. "
            "genemetrics must be run with segments from `cnvkit call` "
            "(rule genemetrics -s {sample}.call.cns), not a plain .cns."
        )

    probes_col = "probes" if "probes" in fieldnames else None

    base = round(args.ploidy)
    amp_cutoff = base + args.amp_offset
    del_cutoff = max(base - args.del_offset, 0)

    def probes_ok(row):
        if probes_col is None:
            return True
        try:
            return int(row[probes_col]) >= args.min_probes
        except (TypeError, ValueError):
            return True

    def cn_of(row):
        try:
            return int(float(row["cn"]))
        except (TypeError, ValueError):
            return None

    amplified, deleted = [], []
    for row in rows:
        if not probes_ok(row):
            continue
        cn = cn_of(row)
        if cn is None:
            continue
        if cn >= amp_cutoff:
            amplified.append(row)
        elif cn <= del_cutoff:
            deleted.append(row)

    amplified.sort(key=lambda r: abs(cn_of(r) - args.ploidy), reverse=True)
    deleted.sort(key=lambda r: abs(cn_of(r) - args.ploidy), reverse=True)

    def write_tsv(path, out_rows, kind, cutoff):
        with open(path, "w", newline="") as fh:
            fh.write(
                f"# ploidy={args.ploidy} {kind}_cutoff={cutoff} "
                f"min_probes={args.min_probes}\n"
            )
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(out_rows)

    def write_genes(path, out_rows):
        gene_col = fieldnames[0]
        with open(path, "w") as fh:
            for row in out_rows:
                fh.write(row[gene_col] + "\n")

    write_tsv(args.amp_out, amplified, "amp", amp_cutoff)
    write_tsv(args.del_out, deleted, "del", del_cutoff)
    write_genes(args.amp_genes_out, amplified)
    write_genes(args.del_genes_out, deleted)

    print(f"amplified (cn >= {amp_cutoff}, ploidy {args.ploidy}): {len(amplified)} genes")
    print(f"deleted   (cn <= {del_cutoff}, ploidy {args.ploidy}): {len(deleted)} genes")


if __name__ == "__main__":
    main()
