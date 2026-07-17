#!/usr/bin/env python3
# 2026 Benjamin J Perry
# MIT License
# Copyright (c) 2026 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

"""Translate an Illumina-format sample sheet into an MGI splitBarcode barcode
file plus the positional arguments splitBarcode needs.

The barcode geometry is measured from the data rather than declared:
RunInfoMetrics.csv reports Read1Length as 101 when the true insert is 100, and
deriving offsets from that field silently destroys the demultiplexing. The
insert length is therefore computed as (observed read length - barcode lengths),
and the resulting offsets are checked against real reads before anything is
written.
"""

import argparse
import csv
import gzip
import json
import sys
from collections import Counter

# splitBarcode accepts these characters in a barcode sequence.
VALID_BASES = set("ACGTN")

# MGI requires barcodes longer than 2 nt.
MIN_BARCODE_LEN = 3


class GeometryError(Exception):
    """Raised when the run geometry cannot be trusted to demultiplex."""


def parse_sample_sheet(path):
    """Return [(Sample_ID, index, index2)] from the [Data] section.

    The sheet is written by Illumina tooling and carries a UTF-8 BOM.
    """
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))

    data_start = None
    for i, row in enumerate(rows):
        if row and row[0].strip().lower() == "[data]":
            data_start = i + 1
            break
    if data_start is None:
        raise GeometryError(f"no [Data] section found in {path}")

    header = [field.strip() for field in rows[data_start]]
    try:
        id_col = header.index("Sample_ID")
        i7_col = header.index("index")
        i5_col = header.index("index2")
    except ValueError:
        raise GeometryError(
            f"[Data] header of {path} must contain Sample_ID, index and index2; got {header}"
        )

    samples = []
    for row in rows[data_start + 1 :]:
        # Trailing commas pad short rows; a row with an empty Sample_ID is padding.
        if not row or not row[0].strip():
            continue
        if row[0].strip().startswith("["):
            break
        samples.append(
            (row[id_col].strip(), row[i7_col].strip().upper(), row[i5_col].strip().upper())
        )

    if not samples:
        raise GeometryError(f"[Data] section of {path} contains no samples")
    return samples


def parse_run_info(path):
    """Return the Name,Value pairs from an MGI RunInfoMetrics.csv."""
    info = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2:
                info[row[0].strip()] = row[1].strip()
    return info


def read_fastq(path, limit=None):
    """Yield sequence lines from a gzipped fastq, at most `limit` of them."""
    with gzip.open(path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 != 1:
                continue
            yield line.rstrip("\n")
            if limit is not None:
                limit -= 1
                if limit <= 0:
                    return


def observed_read_length(fastq):
    """Length of the first read. RunInfoMetrics' Read1Length is not trustworthy."""
    for seq in read_fastq(fastq, limit=1):
        return len(seq)
    raise GeometryError(f"{fastq} contains no reads")


def resolve_geometry(samples, run_info, fastq):
    """Work out where the two barcodes sit in the read.

    Lengths come from the sample sheet (that is what we match against);
    RunInfoMetrics is used only to cross-check them, because its Barcode1/2
    numbering does not follow read order.
    """
    read2_len = int(run_info.get("Read2Length", "0") or 0)
    if read2_len:
        raise GeometryError(
            f"Read2Length is {read2_len}: this run is paired-end, which is not yet supported. "
            "Add splitBarcode -2 and extend the validation gate to search both reads."
        )

    i7_lengths = {len(i7) for _, i7, _ in samples}
    i5_lengths = {len(i5) for _, _, i5 in samples}
    if len(i7_lengths) != 1 or len(i5_lengths) != 1:
        raise GeometryError(
            f"mixed index lengths in the sample sheet (index: {sorted(i7_lengths)}, "
            f"index2: {sorted(i5_lengths)}); splitBarcode needs one geometry per run"
        )
    i7_len = i7_lengths.pop()
    i5_len = i5_lengths.pop()

    # MGI numbers the barcodes independently of read order, so compare as a set.
    declared = {
        int(run_info[key])
        for key in ("Barcode1Length", "Barcode2Length")
        if run_info.get(key, "").isdigit()
    }
    if declared and declared != {i7_len, i5_len}:
        raise GeometryError(
            f"sample sheet index lengths {{{i7_len}, {i5_len}}} do not match "
            f"RunInfoMetrics Barcode1/2Length {declared}"
        )

    total_len = observed_read_length(fastq)
    insert = total_len - i7_len - i5_len
    if insert <= 0:
        raise GeometryError(
            f"read length {total_len} is too short for a {i7_len} bp + {i5_len} bp barcode"
        )

    # Sequence Order Read1-Read2-Dualbarcode-Barcode puts i7 first, then i5.
    return {
        "total_read_len": total_len,
        "insert_len": insert,
        "i7": (insert, i7_len),
        "i5": (insert + i7_len, i5_len),
    }


def validate_offsets(fastq, geometry, expected, sample_reads, min_rate):
    """Check the computed offsets against real reads before trusting them.

    Returns the exact-match rate. This is the gate that catches an off-by-one in
    the insert length, and what makes future paired-end/multi-lane changes safe.
    """
    i7_off, i7_len = geometry["i7"]
    i5_off, i5_len = geometry["i5"]

    matched = 0
    total = 0
    observed = Counter()
    for seq in read_fastq(fastq, limit=sample_reads):
        total += 1
        barcode = seq[i7_off : i7_off + i7_len] + seq[i5_off : i5_off + i5_len]
        if barcode in expected:
            matched += 1
        else:
            observed[barcode] += 1

    if not total:
        raise GeometryError(f"{fastq} contains no reads to validate against")

    rate = matched / total
    if rate < min_rate:
        common = "\n".join(
            f"      {bc!r}  ({count} reads)" for bc, count in observed.most_common(5)
        )
        expected_examples = "\n".join(f"      {bc!r}" for bc in sorted(expected)[:5])
        raise GeometryError(
            f"barcode validation failed: only {rate:.2%} of {total} reads match a sample "
            f"sheet barcode at the computed offsets (need >= {min_rate:.0%}).\n"
            f"  read length {geometry['total_read_len']}, insert {geometry['insert_len']}, "
            f"i7 {i7_len} bp at offset {i7_off}, i5 {i5_len} bp at offset {i5_off}\n"
            f"    most common observed barcodes:\n{common}\n"
            f"    examples expected from the sample sheet:\n{expected_examples}\n"
            "  The offsets or the sample sheet are wrong; refusing to write a barcode file."
        )
    return rate


def check_mgi_constraints(records):
    """Enforce the constraints splitBarcode places on a barcode file."""
    ids = [sample_id for sample_id, _ in records]
    barcodes = [barcode for _, barcode in records]

    for sample_id, barcode in records:
        if not sample_id:
            raise GeometryError("sample sheet contains an empty Sample_ID")
        if len(barcode) <= MIN_BARCODE_LEN - 1:
            raise GeometryError(f"barcode for {sample_id} is too short: {barcode!r}")
        if set(barcode) - VALID_BASES:
            raise GeometryError(
                f"barcode for {sample_id} contains non-ACGTN characters: {barcode!r}"
            )

    for label, values in (("Sample_ID", ids), ("barcode", barcodes)):
        duplicates = [value for value, count in Counter(values).items() if count > 1]
        if duplicates:
            raise GeometryError(
                f"duplicate {label} in sample sheet: {sorted(duplicates)[:5]} "
                f"({len(duplicates)} total); splitBarcode requires both to be unique"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--run-info", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--out-barcodes", required=True)
    parser.add_argument("--out-params", required=True)
    parser.add_argument("--mismatch", type=int, default=1)
    parser.add_argument(
        "--sample-reads",
        type=int,
        default=20000,
        help="reads to stream when validating the offsets",
    )
    parser.add_argument(
        "--min-match-rate",
        type=float,
        default=0.70,
        help="fail below this fraction of reads matching a sample sheet barcode",
    )
    parser.add_argument(
        "--force-insert",
        type=int,
        default=None,
        help="override the measured insert length (for testing the validation gate)",
    )
    args = parser.parse_args()

    samples = parse_sample_sheet(args.sample_sheet)
    run_info = parse_run_info(args.run_info)
    geometry = resolve_geometry(samples, run_info, args.fastq)

    if args.force_insert is not None:
        i7_len = geometry["i7"][1]
        geometry["insert_len"] = args.force_insert
        geometry["i7"] = (args.force_insert, i7_len)
        geometry["i5"] = (args.force_insert + i7_len, geometry["i5"][1])

    records = [(sample_id, i7 + i5) for sample_id, i7, i5 in samples]
    check_mgi_constraints(records)

    rate = validate_offsets(
        args.fastq,
        geometry,
        {barcode for _, barcode in records},
        args.sample_reads,
        args.min_match_rate,
    )

    i7_off, i7_len = geometry["i7"]
    i5_off, i5_len = geometry["i5"]
    params = {
        "b": [[i7_off, i7_len, args.mismatch], [i5_off, i5_len, args.mismatch]],
        "reverse": False,
        "read_length": geometry["total_read_len"],
        "insert_length": geometry["insert_len"],
        "n_samples": len(records),
        "validation_match_rate": round(rate, 4),
        "sequencer_id": run_info.get("Sequencer ID", ""),
    }

    with open(args.out_barcodes, "w", newline="") as handle:
        for sample_id, barcode in records:
            handle.write(f"{sample_id}\t{barcode}\n")

    with open(args.out_params, "w") as handle:
        json.dump(params, handle, indent=2)
        handle.write("\n")

    print(
        f"{len(records)} samples; read {geometry['total_read_len']} bp = insert "
        f"{geometry['insert_len']} + i7 {i7_len} @ {i7_off} + i5 {i5_len} @ {i5_off}; "
        f"{rate:.2%} of {args.sample_reads} reads matched a sample sheet barcode",
        file=sys.stderr,
    )


if __name__ == "__main__":
    try:
        main()
    except GeometryError as exc:
        sys.exit(f"error: {exc}")
