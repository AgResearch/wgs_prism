#!/usr/bin/env python3
# 2026 Benjamin J Perry
# MIT License
# Copyright (c) 2026 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

"""Convert splitBarcode's BarcodeStat.txt into MultiQC custom content.

MultiQC has no splitBarcode module, so without this the demultiplexing rates
are absent from the report entirely.

Sample naming needs care. splitBarcode prefixes every row of BarcodeStat.txt
with the literal string "barcode" (sheet ZU7598125_All becomes barcodeZU7598125_All),
while the fastq files it writes are named <RUN>_<LANE>_<Sample_ID>.fq.gz and so
reach MultiQC as sample <RUN>_<LANE>_<Sample_ID>. Both edits are needed for this
table to line up with the FastQC section rather than silently sitting beside it.
"""

import argparse
import csv
import sys

BARCODE_PREFIX = "barcode"
TOTAL_ROW = "Total"


def parse_barcodestat(path):
    """Return ([(sample_id, correct, corrected, total, percentage)], total_row)."""
    rows = []
    total = None
    with open(path, newline="") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field.strip() for field in line.split("\t")]
            if len(fields) < 5:
                continue
            name = fields[0]
            if name == TOTAL_ROW:
                total = fields
                continue
            if name.startswith(BARCODE_PREFIX):
                name = name[len(BARCODE_PREFIX) :]
            rows.append([name] + fields[1:5])
    return rows, total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--barcodestat", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--sample-prefix",
        default="",
        help="prepended to each ID to match FastQC sample names, e.g. FT150034703_L01_",
    )
    args = parser.parse_args()

    rows, total = parse_barcodestat(args.barcodestat)
    if not rows:
        sys.exit(f"error: no barcode rows parsed from {args.barcodestat}")

    total_pct = total[4] if total else "unknown"
    description = (
        f"Demultiplexing rates reported by MGI splitBarcode. "
        f"{len(rows)} samples; {total_pct}% of reads were assigned to a sample."
    )

    with open(args.out, "w", newline="") as handle:
        handle.write("# id: 'splitbarcode'\n")
        handle.write("# section_name: 'splitBarcode Demultiplexing'\n")
        handle.write(f"# description: \"{description}\"\n")
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#    id: 'splitbarcode_table'\n")
        handle.write("#    namespace: 'splitBarcode'\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["Sample", "Reads Exact", "Reads Corrected", "Reads Total", "Percent of Run"]
        )
        for name, correct, corrected, read_total, percentage in rows:
            writer.writerow(
                [f"{args.sample_prefix}{name}", correct, corrected, read_total, percentage]
            )

    print(
        f"wrote {len(rows)} samples to {args.out} (total demultiplexed {total_pct}%)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
