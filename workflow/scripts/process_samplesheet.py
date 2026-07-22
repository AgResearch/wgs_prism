#!/usr/bin/env python3
# 2026 Benjamin J Perry
# MIT License
# Copyright (c) 2026 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz
"""Sample-sheet processing for MGI ``splitBarcode`` demultiplexing.

Used two ways:

* imported by ``workflow/mgi_splitbarcodes.smk`` at parse time, so that every
  derivable value (read files, ``-b`` arguments, barcode entries) is resolved
  in Python before the DAG is built;
* invoked as a CLI by the ``process_samplesheet`` rule to write the per-lane
  barcode and barcode-parameter files.

Standard library only, and no ``dataclasses`` -- the bare ``python3`` on the
eRI compute nodes is 3.6.8, and the 3.11 interpreter only exists inside the
snakemake module. ``typing.NamedTuple`` keeps this importable under either.

Design principle: **fail loudly**. A wrong barcode offset produces a run that
looks successful while sending ~100% of reads to ``undecoded``. Nothing here
guesses -- the geometry is read from the run's own ``BioInfo.csv``, every
derived value is cross-checked against the reads, and anything inconsistent
raises.

Where the ``-b`` arguments come from: the per-lane ``BioInfo.csv``. The two
platforms use different keys -- T1+ writes ``Read1Len,Read2Len,BarcodeLen,
DualBarcodeLen``, G99 writes ``Read1 Cycles,Read2 Cycles,Barcode,Dual
Barcode`` -- and they sequence the two barcode blocks in *opposite* order,
which G99 states explicitly as ``Sequence Order,Read1-Read2-Dualbarcode-
Barcode`` and T1+ leaves implicit (Barcode first). Read that order and one
rule covers both platforms:

    the sample sheet's ``index`` is the **first** barcode block on the
    machine, ``index2`` the **second**

The mapping therefore looks inverted if key names are compared alone.
Verified against the hand-made ground truth for both platforms: G99
``FT150034703`` -> ``-b 100 8 1 -b 108 10 1``, T1+ ``DL100018479`` ->
``-b 100 10 1 -b 110 8 1``.

Three things ``BioInfo.csv`` will not tell you, so they are taken from the reads:

* **the anchor.** G99 overstates its declared read length by exactly one cycle
  (``Read1 Cycles,101`` for a 100-cycle read, ``151`` for a 150-cycle one --
  it counts the boundary where T1+ reports a length), and trusting it sends the
  whole lane to ``undecoded``. The barcode blocks always sit at the *end* of
  the read, so the anchor is derived as ``read_length - total barcode cycles``
  (118-8-10 = 100, 120-10-10 = 100, 170-10-10 = 150) and the declared value is
  only reported as a cross-check.
* **the orientation.** T1+ paired-end runs carry both indices reverse-
  complemented in ``read_2``; single-end runs carry them forward. Nothing in
  ``BioInfo.csv`` distinguishes the two, so orientation is confirmed by
  counting matches at the (already known) offsets and refused if neither
  orientation matches.
* **the slot order.** Which of the sheet's two indices the instrument
  sequenced first. ``DL100018479`` is index-first, ``DL100018466``
  index2-first. Measured by counting complete sheet *pairs*, and refused
  outright when the plate's pairs are symmetric -- see ``measure_layout``,
  which explains why that case cannot be resolved from the reads at all.

Orientation and slot order are absorbed into the barcode file's *content* by
``oriented_barcode``, so the ``-b`` offsets always come straight from the
geometry.
"""

import argparse
import csv
import gzip
import os
import re
import sys
from typing import Dict, List, NamedTuple, Optional, Sequence, Set, Tuple

# Reads sampled to confirm orientation and read length. The offsets are known
# before the fastq is opened, so this is a yes/no confirmation at two fixed
# positions, not a search: ~0.2 s per lane.
DEFAULT_N_READS = 2000

# Minimum fraction of scanned reads that must match both indices exactly at the
# BioInfo-derived offsets. Observed in production: 0.748 (T1+ SE), 0.904 (G99
# SE) and 0.538-0.776 (T1+ PE, revcomp), so the default leaves a wide margin while
# still catching a sheet paired with the wrong run.
DEFAULT_MIN_HIT_RATE = 0.10

# The winning orientation must beat every other by this factor. A genuine
# orientation wins by orders of magnitude (the runner-up is typically 0.000);
# anything closer means the layout is not being identified reliably, so refuse.
DEFAULT_DISCRIMINATION_RATIO = 2.0

# splitBarcode's -b startCycle is 0-based: the production `-b 100 10 1` on a
# 120 bp read with a 100 bp insert lines up exactly with a Python string index,
# so derived offsets are passed through without adjustment.
#
# In paired-end mode the barcode blocks follow read 2, and splitBarcode counts
# -b across the *concatenated* read, not from the start of read 2. Measured on
# DL100018466 L01 (1M read pairs): -b 300/310 decodes 79.394%, -b 150/160
# decodes 0.0001%. Its -p flag (read 1 length) proved inert on 2.0.0-4 -- the
# same 79.394% with and without it -- but is still emitted to state the intent.
# Kept as a config choice so a future version can be re-measured, not re-coded.
PE_OFFSET_BASES = ("read2-relative", "concatenated")

_SECTION_RE = re.compile(r"^\[(?P<name>[^\]]+)\]")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_VALID_BASES = set("ACGTN")

FORWARD = "forward"
REVCOMP = "revcomp"
_ORIENTATIONS = (
    (FORWARD, FORWARD),
    (FORWARD, REVCOMP),
    (REVCOMP, FORWARD),
    (REVCOMP, REVCOMP),
)

# Which of the sheet's two indices the instrument sequenced *first*. Not
# derivable from BioInfo.csv, and not always derivable from the reads either --
# see the swap-pair degeneracy in `measure_layout`.
INDEX_FIRST = "index-first"
INDEX2_FIRST = "index2-first"
SLOT_ORDERS = (INDEX_FIRST, INDEX2_FIRST)
AUTO = "auto"


class SamplesheetError(Exception):
    """Raised for any sample sheet the pipeline refuses to process."""


class BarcodeLayoutError(Exception):
    """Raised when the barcode layout cannot be established."""


class BioInfoError(Exception):
    """Raised when a run's BioInfo.csv is missing, unreadable or incomplete."""


def revcomp(seq: str) -> str:
    """Reverse complement of an IUPAC-free (ACGTN) sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def _oriented(seq: str, orientation: str) -> str:
    return seq if orientation == FORWARD else revcomp(seq)


# --------------------------------------------------------------------------
# Sample sheet reading
# --------------------------------------------------------------------------


class Samplesheet(NamedTuple):
    """A parsed MGI sample sheet.

    ``sections`` keeps every section verbatim (``Header``, ``Reads``, ``Data``,
    ``GenerateKeyfile``, ...) so later work -- sample-name processing, sheet
    validation -- can reach sections this module does not interpret. Only
    ``[Data]`` is promoted to dictionaries.
    """

    path: str
    sections: Dict[str, List[List[str]]]
    data: List[Dict[str, str]]
    data_columns: List[str]

    @property
    def header(self) -> Dict[str, str]:
        """``[Header]`` as key -> value (e.g. ``Flowcell``, ``Instrument``)."""
        out: Dict[str, str] = {}
        for row in self.sections.get("Header", []):
            if row and row[0]:
                out[row[0]] = row[1] if len(row) > 1 else ""
        return out


def read_samplesheet(path: str) -> Samplesheet:
    """Parse an MGI sample sheet into sections plus typed ``[Data]`` rows.

    Handles both platforms' quirks: T1+ sheets are CRLF with stray leading
    spaces in values (``Flowcell, DL100018749``), G99 sheets are LF; the two
    carry different column sets and different sections. Every field is
    stripped, and ``[Data]`` columns are addressed by name, never by position.
    """
    if not os.path.exists(path):
        raise SamplesheetError(f"sample sheet does not exist: {path}")

    # sections is a plain dict, which preserves insertion order, so the order
    # the sheet declares its sections in survives without tracking it here.
    sections: Dict[str, List[List[str]]] = {}
    current: Optional[str] = None

    # newline="" lets csv handle the CRLF sheets without leaving stray \r.
    with open(path, "r", newline="", encoding="utf-8-sig") as handle:
        for raw in csv.reader(handle):
            row = [field.strip() for field in raw]
            if not any(row):
                continue
            match = _SECTION_RE.match(row[0])
            if match:
                current = match.group("name").strip()
                sections.setdefault(current, [])
                continue
            if current is None:
                raise SamplesheetError(
                    f"{path}: content before the first [Section] header: {row!r}"
                )
            sections[current].append(row)

    if "Data" not in sections:
        raise SamplesheetError(f"{path}: no [Data] section found")
    if not sections["Data"]:
        raise SamplesheetError(f"{path}: [Data] section is empty")

    columns = sections["Data"][0]
    data: List[Dict[str, str]] = []
    for row in sections["Data"][1:]:
        # Pad short rows rather than zipping them away, so a missing trailing
        # field reads as empty instead of shifting every column left.
        padded = list(row) + [""] * (len(columns) - len(row))
        data.append(dict(zip(columns, padded)))

    if not data:
        raise SamplesheetError(f"{path}: [Data] section has a header but no sample rows")

    return Samplesheet(path=path, sections=sections, data=data, data_columns=columns)


def _column(ss: Samplesheet, *candidates: str) -> str:
    """Resolve the first present column name, case-insensitively."""
    lowered = {name.lower(): name for name in ss.data_columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    raise SamplesheetError(
        f"{ss.path}: [Data] has none of the expected columns {list(candidates)}; found {ss.data_columns}"
    )


# --------------------------------------------------------------------------
# Lanes
# --------------------------------------------------------------------------


def parse_lanes(value: str) -> List[int]:
    """Expand a ``Lanes`` field, which is a **digit string**.

    ``"12"`` -> ``[1, 2]``, ``"34"`` -> ``[3, 4]``, ``"1"`` -> ``[1]``. This is
    membership, not equality: a row with ``Lanes=12`` belongs to L01 *and* L02.
    Anything other than non-zero digits is refused rather than guessed at:
    comma- and dash-separated sheets have not been validated against.
    """
    text = (value or "").strip()
    if not text:
        raise SamplesheetError("empty Lanes field")
    if not text.isdigit():
        raise SamplesheetError(
            f"Lanes field {value!r} is not a digit string; expected e.g. '1', '12', '34'"
        )
    lanes = sorted({int(char) for char in text})
    if any(lane == 0 for lane in lanes):
        raise SamplesheetError(f"Lanes field {value!r} contains lane 0")
    return lanes


def lane_label(lane: int) -> str:
    """``1`` -> ``"L01"``."""
    return f"L{lane:02d}"


def lane_number(lane) -> int:
    """Accept ``"L01"``, ``"l01"``, ``"1"`` or ``1`` and return ``1``."""
    if isinstance(lane, int):
        return lane
    text = str(lane).strip()
    if text[:1].upper() == "L":
        text = text[1:]
    if not text.isdigit():
        raise SamplesheetError(f"cannot interpret lane {lane!r}")
    return int(text)


def lanes_in_samplesheet(ss: Samplesheet) -> List[str]:
    """Sorted, deduplicated lane labels across every ``[Data]`` row."""
    column = _column(ss, "Lanes", "Lane")
    found: Set[int] = set()
    for row in ss.data:
        found.update(parse_lanes(row[column]))
    return [lane_label(lane) for lane in sorted(found)]


def samples_for_lane(ss: Samplesheet, lane) -> List[Dict[str, str]]:
    """Rows belonging to ``lane`` -- by digit membership, not equality."""
    wanted = lane_number(lane)
    column = _column(ss, "Lanes", "Lane")
    rows = [row for row in ss.data if wanted in parse_lanes(row[column])]
    if not rows:
        raise SamplesheetError(
            f"{ss.path}: no samples for lane {lane_label(wanted)}; "
            f"lanes present: {', '.join(lanes_in_samplesheet(ss))}"
        )
    return rows


# --------------------------------------------------------------------------
# Barcodes
# --------------------------------------------------------------------------


def _validate_index(seq: str, what: str, sample_id: str) -> str:
    upper = seq.strip().upper()
    bad = set(upper) - _VALID_BASES
    if bad:
        raise SamplesheetError(
            f"sample {sample_id}: {what} contains non-ACGTN characters {sorted(bad)}: {seq!r}"
        )
    return upper


def barcode_entries(
    ss: Samplesheet,
    lane,
    orientation: Tuple[str, str] = (FORWARD, FORWARD),
    slot_order: str = INDEX_FIRST,
) -> List[Tuple[str, str]]:
    """``[(sample_id, slot 1 + slot 2), ...]`` for one lane.

    The two indices concatenated in ``-b`` order, which is the format
    ``splitBarcode -B`` expects: ``Sample_ID<TAB>barcode``, headerless and
    tab-delimited. Single-index lanes (empty ``index2``) collapse to just
    ``index``.

    ``orientation`` is how the *reads* carry each index and ``slot_order``
    which index the instrument sequenced first, both as resolved by
    ``measure_layout``. Everything the run does differently from the sheet is
    absorbed here, in the barcode file's *content*, so the ``-b`` offsets stay
    as the geometry gives them. (splitBarcode's own ``-r`` flag is deliberately
    unused: it applies to the whole concatenated barcode and cannot express a
    per-index difference.)
    """
    if slot_order not in SLOT_ORDERS:
        raise SamplesheetError(
            f"slot_order must be one of {list(SLOT_ORDERS)}, got {slot_order!r}"
        )
    rows = samples_for_lane(ss, lane)
    id_column = _column(ss, "Sample_ID", "SampleID", "Sample_Name")
    index_column = _column(ss, "index", "Index", "I7_Index", "index1")
    try:
        index2_column: Optional[str] = _column(ss, "index2", "Index2", "I5_Index")
    except SamplesheetError:
        index2_column = None

    entries: List[Tuple[str, str]] = []
    lengths: Set[Tuple[int, int]] = set()
    seen: Set[str] = set()

    for row in rows:
        sample_id = row[id_column].strip()
        if not sample_id:
            raise SamplesheetError(f"{ss.path}: a row in lane {lane} has an empty Sample_ID")
        if sample_id in seen:
            raise SamplesheetError(
                f"{ss.path}: sample {sample_id} appears more than once in lane {lane}"
            )
        seen.add(sample_id)

        index = _validate_index(row[index_column], "index", sample_id)
        if not index:
            raise SamplesheetError(f"{ss.path}: sample {sample_id} has an empty index")
        index2 = ""
        if index2_column is not None:
            index2 = _validate_index(row.get(index2_column, ""), "index2", sample_id)

        lengths.add((len(index), len(index2)))
        entries.append((sample_id, oriented_barcode(index, index2, orientation, slot_order)))

    if len(lengths) > 1:
        raise SamplesheetError(
            f"{ss.path}: lane {lane_label(lane_number(lane))} mixes index lengths "
            f"{sorted(lengths)}; splitBarcode needs one -b layout for the whole lane"
        )
    return entries


def oriented_barcode(
    index: str, index2: str, orientation: Tuple[str, str], slot_order: str
) -> str:
    """The barcode string as the instrument writes it, slot 1 then slot 2.

    One place decides how a sheet row becomes a read's barcode, so the ``-B``
    file and the layout scoring in ``measure_layout`` cannot drift apart.
    """
    first = _oriented(index, orientation[0])
    second = _oriented(index2, orientation[1])
    if slot_order == INDEX2_FIRST:
        return second + first
    return first + second


def index_pairs(ss: Samplesheet, lane) -> List[Tuple[str, str]]:
    """``[(index, index2), ...]`` for one lane, in sheet order."""
    rows = samples_for_lane(ss, lane)
    index_column = _column(ss, "index", "Index", "I7_Index", "index1")
    try:
        index2_column: Optional[str] = _column(ss, "index2", "Index2", "I5_Index")
    except SamplesheetError:
        index2_column = None
    return [
        (
            row[index_column].strip().upper(),
            row.get(index2_column, "").strip().upper() if index2_column else "",
        )
        for row in rows
    ]


def write_barcode_file(entries: Sequence[Tuple[str, str]], path: str) -> None:
    """Write the headerless, tab-delimited ``-B`` file."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w") as handle:
        for sample_id, barcode in entries:
            handle.write(f"{sample_id}\t{barcode}\n")


# --------------------------------------------------------------------------
# Read files
# --------------------------------------------------------------------------


class ReadFiles(NamedTuple):
    """The per-lane fastq(s), and which one carries the barcode."""

    r1: str
    r2: Optional[str]

    @property
    def is_pe(self) -> bool:
        return self.r2 is not None

    @property
    def barcode_fastq(self) -> str:
        """MGI appends the barcode to the end of the *last* read."""
        return self.r2 if self.r2 is not None else self.r1

    @property
    def cli_args(self) -> str:
        """``-1 a.fq.gz`` or ``-1 a.fq.gz -2 b.fq.gz``."""
        if self.r2 is None:
            return f"-1 {self.r1}"
        return f"-1 {self.r1} -2 {self.r2}"

    @property
    def paths(self) -> List[str]:
        return [self.r1] if self.r2 is None else [self.r1, self.r2]


def read_files(run_dir: str, run: str, lane) -> ReadFiles:
    """Locate a lane's fastq(s) and return them in a typed tuple.

    Layout is ``{run_dir}/{LANE}/{RUN}_{LANE}_read.fq.gz`` for single-end and
    ``..._read_1.fq.gz`` / ``..._read_2.fq.gz`` for paired-end.
    """
    label = lane_label(lane_number(lane))
    lane_dir = os.path.join(run_dir, label)
    single = os.path.join(lane_dir, f"{run}_{label}_read.fq.gz")
    first = os.path.join(lane_dir, f"{run}_{label}_read_1.fq.gz")
    second = os.path.join(lane_dir, f"{run}_{label}_read_2.fq.gz")

    if os.path.exists(first) and os.path.exists(second):
        return ReadFiles(r1=first, r2=second)
    if os.path.exists(single):
        return ReadFiles(r1=single, r2=None)
    if os.path.exists(first) != os.path.exists(second):
        raise SamplesheetError(
            f"lane {label}: found only one half of the paired-end pair in {lane_dir} -- "
            "refusing to demultiplex an incomplete pair"
        )
    raise SamplesheetError(
        f"lane {label}: no reads found. Expected one of:\n  {single}\n  {first} (+ _2)"
    )


def _sequence_lines(fastq: str, n_reads: int) -> List[str]:
    """The first ``n_reads`` sequence lines of a gzipped fastq."""
    lines: List[str] = []
    with gzip.open(fastq, "rt") as handle:
        for number, line in enumerate(handle):
            if number % 4 != 1:
                continue
            lines.append(line.rstrip("\r\n"))
            if len(lines) >= n_reads:
                break
    if not lines:
        raise BarcodeLayoutError(f"no reads could be read from {fastq}")
    return lines


def read_length(fastq: str, n_reads: int = 20) -> int:
    """Modal sequence-line length of the first few reads."""
    counts: Dict[int, int] = {}
    for line in _sequence_lines(fastq, n_reads):
        counts[len(line)] = counts.get(len(line), 0) + 1
    return max(counts.items(), key=lambda item: (item[1], item[0]))[0]


# --------------------------------------------------------------------------
# BioInfo.csv -- the run's own description of its cycle structure
# --------------------------------------------------------------------------

# Both platforms sequence the insert first and the barcode cycles last
# (T1+ says so as BarcodePosition,BarcodePosEnd). They differ in which of the
# two barcode blocks comes first, so the order is read from the file where the
# file states it and defaults to the T1+ convention where it does not.
BARCODE = "Barcode"
DUAL_BARCODE = "DualBarcode"
_T1_BLOCK_ORDER = (BARCODE, DUAL_BARCODE)

# Key aliases, normalised (lowercased, non-alphanumerics dropped) so that
# "Read1Len", "Read1 Cycles" and "read1_len" all resolve.
_BIOINFO_KEYS = {
    "read1_len": ("read1len", "read1cycles"),
    "read2_len": ("read2len", "read2cycles"),
    "barcode_len": ("barcodelen", "barcode"),
    "dual_barcode_len": ("dualbarcodelen", "dualbarcode"),
    "sequence_order": ("sequenceorder",),
    "machine_id": ("machineid",),
}


def _normalise_key(key: str) -> str:
    return "".join(char for char in key.lower() if char.isalnum())


class BarcodeBlock(NamedTuple):
    """One sequenced barcode block: its machine name and its cycle count."""

    name: str
    length: int


class BioInfo(NamedTuple):
    """The cycle structure of one lane, as the instrument recorded it."""

    path: str
    platform: str
    read1_len: int
    read2_len: int
    blocks: List[BarcodeBlock]
    values: Dict[str, str]

    @property
    def is_pe(self) -> bool:
        return self.read2_len > 0

    @property
    def barcode_cycles(self) -> int:
        """Total cycles the barcode blocks occupy at the end of the read."""
        return sum(block.length for block in self.blocks)

    @property
    def declared_insert_len(self) -> int:
        """The read length the instrument claims -- read 2's in PE mode.

        Used only as a cross-check: G99 overstates ``Read1 Cycles`` by one.
        """
        return self.read2_len if self.is_pe else self.read1_len


def bioinfo_path(run_dir: str, lane) -> str:
    """``{run_dir}/{LANE}/BioInfo.csv``, falling back to the run root."""
    label = lane_label(lane_number(lane))
    per_lane = os.path.join(run_dir, label, "BioInfo.csv")
    if os.path.exists(per_lane):
        return per_lane
    at_root = os.path.join(run_dir, "BioInfo.csv")
    if os.path.exists(at_root):
        return at_root
    raise BioInfoError(
        f"lane {label}: no BioInfo.csv found. Expected one of:\n  {per_lane}\n  {at_root}"
    )


def _bioinfo_int(values: Dict[str, str], field: str, required: bool = True) -> int:
    for alias in _BIOINFO_KEYS[field]:
        if alias in values and values[alias] != "":
            text = values[alias]
            try:
                return int(text)
            except ValueError:
                # from None: the ValueError adds nothing the message lacks.
                raise BioInfoError(
                    f"{field} is not an integer: {text!r}"
                ) from None
    if required:
        raise BioInfoError(
            f"no {field} field; looked for {list(_BIOINFO_KEYS[field])}"
        )
    return 0


def _block_order(values: Dict[str, str]) -> Tuple[str, ...]:
    """Barcode block order, from ``Sequence Order`` when the run states it.

    G99 writes ``Read1-Read2-Dualbarcode-Barcode``; T1+ writes nothing and
    sequences Barcode first. Read1/Read2 tokens are dropped -- only the
    relative order of the two barcode blocks matters here.
    """
    raw = ""
    for alias in _BIOINFO_KEYS["sequence_order"]:
        if values.get(alias):
            raw = values[alias]
            break
    if not raw:
        return _T1_BLOCK_ORDER

    order: List[str] = []
    for token in raw.split("-"):
        name = _normalise_key(token)
        if name.startswith("read"):
            continue
        if name == "dualbarcode":
            order.append(DUAL_BARCODE)
        elif name == "barcode":
            order.append(BARCODE)
        else:
            raise BioInfoError(
                f"Sequence Order {raw!r} contains an unrecognised block {token!r}"
            )
    if not order:
        raise BioInfoError(f"Sequence Order {raw!r} names no barcode blocks")
    return tuple(order)


def read_bioinfo(path: str) -> BioInfo:
    """Parse a lane's ``BioInfo.csv`` into read lengths and barcode blocks."""
    if not os.path.exists(path):
        raise BioInfoError(f"BioInfo.csv does not exist: {path}")

    values: Dict[str, str] = {}
    with open(path, "r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if not row or not row[0].strip():
                continue
            values[_normalise_key(row[0])] = row[1].strip() if len(row) > 1 else ""

    # Machine ID is not dependable: FT150034546 and FT150034753 are G99 runs
    # that report a bare serial (R21007100260015). The *schema* is the reliable
    # signal -- only G99 writes "Read1 Cycles" -- and it is what everything
    # below actually keys off, so fall back to it for the label too.
    machine = values.get("machineid", "")
    if "G99" in machine.upper():
        platform = "G99"
    elif "T1" in machine.upper():
        platform = "T1+"
    elif "read1cycles" in values:
        platform = "G99"
    elif "read1len" in values:
        platform = "T1+"
    else:
        platform = machine or "unknown"

    lengths = {
        BARCODE: _bioinfo_int(values, "barcode_len"),
        DUAL_BARCODE: _bioinfo_int(values, "dual_barcode_len", required=False),
    }
    # A zero-length block was not sequenced. Dropping it here lets a
    # single-barcode run collapse to one -b without a special case.
    blocks = [
        BarcodeBlock(name=name, length=lengths[name])
        for name in _block_order(values)
        if lengths[name] > 0
    ]
    if not blocks:
        raise BioInfoError(f"{path}: no barcode cycles declared")

    return BioInfo(
        path=path,
        platform=platform,
        read1_len=_bioinfo_int(values, "read1_len"),
        read2_len=_bioinfo_int(values, "read2_len", required=False),
        blocks=blocks,
        values=values,
    )


# --------------------------------------------------------------------------
# Barcode layout -- geometry from BioInfo.csv, confirmed against the reads
# --------------------------------------------------------------------------


class BarcodeLayout(NamedTuple):
    """Where the barcode sits: derived from BioInfo, confirmed on the reads."""

    fastq: str
    platform: str
    offset: int
    offset2: Optional[int]
    blocks: List[BarcodeBlock]
    index_len: int
    index2_len: int
    orientation: Tuple[str, str]
    hit_rate: float
    runner_up_rate: float
    n_scanned: int
    read_len: int
    declared_insert_len: int
    is_pe: bool
    slot_order: str = INDEX_FIRST
    slot_order_declared: bool = False
    symmetric_pairs: int = 0
    read1_len: Optional[int] = None

    @property
    def slot_lengths(self) -> Tuple[int, int]:
        """``(slot 1, slot 2)`` index lengths, in the order the machine read them."""
        if self.slot_order == INDEX2_FIRST:
            return self.index2_len, self.index_len
        return self.index_len, self.index2_len

    @property
    def anchor_disagrees(self) -> bool:
        """True when BioInfo's declared read length is not the real one.

        Reported, not fatal -- the read length is the authority for the anchor.
        """
        return self.offset != self.declared_insert_len

    @property
    def anchor_is_g99_fencepost(self) -> bool:
        """True for G99's systematic one-cycle overstatement.

        G99 declares ``Read1 Cycles,101`` for a 100-cycle read 1 and
        ``151`` for a 150-cycle one, while the same file's ``Sequence Type``
        (``SE100+8+10``, ``PE150+10``) gives the true figure -- it counts the
        boundary rather than the length. Confirmed on FT150034703,
        FT150034546 and FT150034753; T1+ never does it. Distinguished from a
        genuine mismatch so that a routine G99 run does not log something that
        reads like a corrupt file.
        """
        return self.declared_insert_len - self.offset == 1

    def describe(self) -> str:
        blocks = " + ".join(f"{block.name} {block.length}" for block in self.blocks)
        slot1, slot2 = self.slot_lengths
        first, second = (
            ("index2", "index") if self.slot_order == INDEX2_FIRST else ("index", "index2")
        )
        mode = "PE (read 2)" if self.is_pe else "SE"
        declared = " (declared)" if self.slot_order_declared else ""
        text = (
            f"{self.platform} {mode}: read {self.read_len} bp = "
            f"insert {self.offset} + [{blocks}]; {first} {slot1} @ {self.offset}"
        )
        if slot2:
            text += f", {second} {slot2} @ {self.offset2}"
        text += (
            f"; orientation index={self.orientation[0]} index2={self.orientation[1]}"
            f"; slot order {self.slot_order}{declared}"
        )
        text += (
            f"; hit-rate {self.hit_rate:.3f} (runner-up {self.runner_up_rate:.3f})"
            f" over {self.n_scanned} reads"
        )
        if self.anchor_is_g99_fencepost:
            text += (
                f"\n  note: BioInfo.csv declares {self.declared_insert_len} cycles of "
                f"insert against {self.offset} measured - the usual G99 one-cycle "
                "overstatement; using the reads."
            )
        elif self.anchor_disagrees:
            text += (
                f"\n  WARNING: BioInfo.csv declares {self.declared_insert_len} cycles "
                f"of insert but the reads give {self.offset} - that is not the known "
                "G99 off-by-one; check the run before trusting this lane."
            )
        if self.slot_order_declared and self.symmetric_pairs:
            text += (
                f"\n  note: {self.symmetric_pairs} index pairs in this lane are "
                "symmetric, so the slot order cannot be checked against the reads - "
                "it rests on the declared BARCODE_SLOT_ORDER."
            )
        return text


def barcode_geometry(
    bio: BioInfo, read_len: int, index_len: int, index2_len: int
) -> Tuple[int, Optional[int]]:
    """``(offset, offset2)`` for one lane, in barcode-fastq coordinates.

    The barcode blocks sit at the **end** of the read, so the anchor is
    ``read_len - total barcode cycles`` rather than BioInfo's declared read
    length -- G99 overstates that by one cycle and anchoring on it sends the
    whole lane to ``undecoded``.

    ``index2`` starts at ``offset + len(block 1)``, **not** at
    ``offset + len(index)``: a block can be longer than the index the sheet
    puts in it (T1+ sequences 10 barcode cycles for an 8-base index), and those
    unused cycles sit between the two indices.
    """
    cycles = bio.barcode_cycles
    if read_len <= cycles:
        raise BarcodeLayoutError(
            f"reads are {read_len} bp but {bio.path} declares {cycles} barcode cycles - "
            "there is no insert left"
        )

    offset = read_len - cycles
    first_block = bio.blocks[0]
    if index_len > first_block.length:
        raise BarcodeLayoutError(
            f"the sheet's index is {index_len} bases but {bio.path} sequenced only "
            f"{first_block.length} cycles for the first barcode block "
            f"({first_block.name})"
        )

    if not index2_len:
        return offset, None

    if len(bio.blocks) < 2:
        raise BarcodeLayoutError(
            f"the sheet has a dual index but {bio.path} declares only one barcode block "
            f"({first_block.name} {first_block.length} cycles)"
        )
    second_block = bio.blocks[1]
    if index2_len > second_block.length:
        raise BarcodeLayoutError(
            f"the sheet's index2 is {index2_len} bases but {bio.path} sequenced only "
            f"{second_block.length} cycles for the second barcode block "
            f"({second_block.name})"
        )
    return offset, offset + first_block.length


def _symmetric_pairs(pairs: Sequence[Tuple[str, str]]) -> int:
    """How many ``(index, index2)`` pairs also appear swapped in the same lane.

    Applying the slot order backwards maps each of these onto a pair the lane
    already contains, so the reads cannot tell the two orders apart -- see
    ``measure_layout``.
    """
    unique = set(pairs)
    return sum(1 for first, second in unique if (second, first) in unique)


def measure_layout(
    fastq: str,
    bio: BioInfo,
    pairs: Sequence[Tuple[str, str]],
    n_reads: int = DEFAULT_N_READS,
    min_hit_rate: float = DEFAULT_MIN_HIT_RATE,
    discrimination_ratio: float = DEFAULT_DISCRIMINATION_RATIO,
    slot_order: str = AUTO,
) -> BarcodeLayout:
    """Derive the geometry from ``bio``, then resolve the rest against the reads.

    ``BioInfo.csv`` gives the offsets but not two things that vary run to run,
    because they are properties of how the library kit's sheet was written
    rather than of the sequencing:

    * **orientation** -- T1+ paired-end runs carry both indices reverse-
      complemented in ``read_2``; the single-end runs carry them forward;
    * **slot order** -- which of the sheet's two indices the instrument
      sequenced first. ``DL100018479`` (SE) sequences ``index`` first;
      ``DL100018466`` (PE) sequences ``index2`` first.

    Both are measured here by counting, at the two known offsets, how many
    reads carry a **complete sheet pair** -- not two indices matched
    independently, because that cannot tell the orderings apart.

    The swap-pair degeneracy
    ------------------------
    When a plate uses the same sequences in both index positions, every sample
    ``(A, B)`` tends to have a partner sample ``(B, A)``. Both slot orders then
    decode *exactly* the same number of reads -- they just hand each read to
    the other sample of the pair. ``DL100018466`` is such a plate: all 96 pairs
    in L01 are symmetric, and both orders decode 79.394096%. Nothing in the
    reads can break that tie, so this refuses to pick one and asks for
    ``slot_order`` to be declared. (For that run the tie was broken by the two
    negative controls: 7 and 6 reads under ``index2-first``, 7,922 and 11,243
    under ``index-first``.)

    Three gates, any of which aborts:

    1. the hit rate is below ``min_hit_rate`` -- the sheet and the run are not
       a matching pair, or the geometry is wrong;
    2. the winning orientation does not beat the runner-up by
       ``discrimination_ratio``;
    3. the slot order is undecidable and none was declared.
    """
    if not pairs:
        raise BarcodeLayoutError(f"no indices supplied to check {fastq}")
    if slot_order not in SLOT_ORDERS and slot_order != AUTO:
        raise BarcodeLayoutError(
            f"slot_order must be {list(SLOT_ORDERS)} or {AUTO!r}, got {slot_order!r}"
        )

    lengths = {(len(index), len(index2)) for index, index2 in pairs}
    if len(lengths) > 1:
        raise BarcodeLayoutError(
            f"lane mixes index lengths {sorted(lengths)}; splitBarcode needs one -b layout for "
            f"the whole lane.\n  fastq: {fastq}"
        )
    index_len, index2_len = lengths.pop()

    lines = _sequence_lines(fastq, n_reads)
    # Use only the modal read length: a handful of trimmed or truncated reads
    # must not move the anchor for everything else.
    counts: Dict[int, int] = {}
    for line in lines:
        counts[len(line)] = counts.get(len(line), 0) + 1
    read_len = max(counts.items(), key=lambda item: (item[1], item[0]))[0]
    usable = [line for line in lines if len(line) == read_len]
    n_scanned = len(usable)

    # Candidate orders. With equal-length indices the geometry is identical
    # either way, but the barcode *content* differs, so both are still scored.
    # A single-index lane has one slot to fill, so slot order cannot apply: a
    # declared order is ignored rather than refused, because BARCODE_SLOT_ORDER
    # is set per run and a run may hold both single- and dual-index lanes.
    # Without this, index2-first on a single-index lane maps the sole index to
    # slot 2 and emits a zero-length -b.
    declared = slot_order != AUTO and bool(index2_len)
    if not index2_len:
        orders = (INDEX_FIRST,)
    elif declared:
        orders = (slot_order,)
    else:
        orders = SLOT_ORDERS

    scored: Dict[Tuple[Tuple[str, str], str], int] = {}
    geometry: Dict[str, Tuple[int, Optional[int]]] = {}
    failures: List[str] = []

    for order in orders:
        slot1_len, slot2_len = (
            (index2_len, index_len) if order == INDEX2_FIRST else (index_len, index2_len)
        )
        try:
            offset, offset2 = barcode_geometry(bio, read_len, slot1_len, slot2_len)
        except BarcodeLayoutError as error:
            # This order does not fit the blocks the instrument sequenced;
            # that is an answer, not a failure, unless no order fits.
            failures.append(f"{order}: {error}")
            continue
        geometry[order] = (offset, offset2)

        for orientation in _ORIENTATIONS:
            if not index2_len and orientation[1] == REVCOMP:
                continue  # no index2 to flip; would double-count the same candidate
            expected = {
                oriented_barcode(index, index2, orientation, order) for index, index2 in pairs
            }
            span2 = offset2 + slot2_len if offset2 is not None else offset + slot1_len
            hits = 0
            for line in usable:
                observed = line[offset : offset + slot1_len]
                if offset2 is not None:
                    observed += line[offset2:span2]
                if observed in expected:
                    hits += 1
            scored[(orientation, order)] = hits

    if not scored:
        # Joined into a local first: an f-string expression cannot contain a
        # backslash before 3.12.
        detail = "\n  ".join(failures)
        raise BarcodeLayoutError(
            f"no slot order fits the barcode blocks {bio.path} declares:\n  {detail}"
        )

    ranked = sorted(scored.items(), key=lambda item: item[1], reverse=True)
    (best_orientation, best_order), best_hits = ranked[0]
    runner_hits = ranked[1][1] if len(ranked) > 1 else 0
    offset, offset2 = geometry[best_order]

    # The runner-up that matters for orientation is the best *other* orientation
    # under the winning order; a different slot order is a separate question.
    orientation_runner = max(
        [hits for (orientation, order), hits in scored.items()
         if order == best_order and orientation != best_orientation] or [0]
    )
    order_runner = max(
        [hits for (orientation, order), hits in scored.items()
         if orientation == best_orientation and order != best_order] or [0]
    )

    layout = BarcodeLayout(
        fastq=fastq,
        platform=bio.platform,
        offset=offset,
        offset2=offset2,
        blocks=list(bio.blocks),
        index_len=index_len,
        index2_len=index2_len,
        orientation=best_orientation,
        slot_order=best_order,
        slot_order_declared=declared,
        hit_rate=best_hits / float(n_scanned),
        runner_up_rate=runner_hits / float(n_scanned),
        n_scanned=n_scanned,
        read_len=read_len,
        declared_insert_len=bio.declared_insert_len,
        is_pe=bio.is_pe,
        symmetric_pairs=_symmetric_pairs(pairs),
    )

    if layout.hit_rate < min_hit_rate:
        raise BarcodeLayoutError(
            f"the sheet's indices are not where {bio.path} says they are.\n  {layout.describe()}\n"
            f"  hit-rate is below the {min_hit_rate:.3f} threshold - the sample sheet and the "
            f"run are probably not a matching pair.\n  fastq: {fastq}"
        )

    if orientation_runner and best_hits < discrimination_ratio * orientation_runner:
        raise BarcodeLayoutError(
            f"barcode orientation is ambiguous in {fastq}.\n  {layout.describe()}\n"
            "  the winning orientation does not beat the runner-up by the "
            f"required {discrimination_ratio}x - refusing to guess."
        )

    if not declared and order_runner and best_hits < discrimination_ratio * order_runner:
        raise BarcodeLayoutError(
            "cannot tell which index this run sequenced first, and guessing "
            f"would silently swap samples.\n  {layout.describe()}\n"
            f"  {layout.symmetric_pairs} of {len(set(pairs))} index pairs in this lane "
            f"also appear swapped, so both slot orders decode the same {best_hits} "
            "reads and hand each one to the other sample of the pair.\n"
            "  Resolve it with external evidence -- a negative control or an empty "
            "well should receive almost no reads under the correct order -- then set "
            f"BARCODE_SLOT_ORDER to '{INDEX_FIRST}' or '{INDEX2_FIRST}' in the config."
        )

    return layout


def barcode_params(
    layout: BarcodeLayout,
    mismatch: int = 1,
    pe_offset_base: str = "read2-relative",
) -> List[Tuple[int, int, int]]:
    """``[(start, length, mismatch), ...]`` -- one tuple per ``-b``.

    Each ``-b`` takes only the bases the sheet actually uses, so a block with
    unused trailing cycles (T1+ sequences 10, the sheet supplies 8) yields
    ``-b <start> 8 <mismatch>`` and the 2 junk cycles are simply never read.

    For PE the geometry is expressed in *read 2* coordinates. Whether
    splitBarcode counts ``-b`` the same way, or across the concatenated read,
    is what ``pe_offset_base`` selects -- so the experiment flips a config
    value rather than editing this function.
    """
    if pe_offset_base not in PE_OFFSET_BASES:
        raise ValueError(
            f"PE_OFFSET_BASE must be one of {list(PE_OFFSET_BASES)}, got {pe_offset_base!r}"
        )

    shift = 0
    if layout.is_pe and pe_offset_base == "concatenated":
        if layout.read1_len is None:
            raise BarcodeLayoutError(
                "PE_OFFSET_BASE=concatenated needs read 1's length, which was not measured"
            )
        shift = layout.read1_len

    slot1_len, slot2_len = layout.slot_lengths
    params = [(layout.offset + shift, slot1_len, mismatch)]
    if slot2_len:
        params.append((layout.offset2 + shift, slot2_len, mismatch))
    return params


def format_b_args(params: Sequence[Tuple[int, int, int]]) -> str:
    """``-b 100 10 1 -b 110 8 1``."""
    return " ".join(f"-b {start} {length} {mismatch}" for start, length, mismatch in params)


# --------------------------------------------------------------------------
# Per-lane plan
# --------------------------------------------------------------------------


class LanePlan(NamedTuple):
    """Everything the rules need for one lane, resolved at parse time."""

    lane: str
    entries: List[Tuple[str, str]]
    reads: ReadFiles
    bio: BioInfo
    layout: BarcodeLayout
    params: List[Tuple[int, int, int]]
    pe_offset_base: str

    @property
    def b_args(self) -> str:
        return format_b_args(self.params)

    @property
    def orientation_arg(self) -> str:
        """``ff``/``fr``/``rf``/``rr`` -- how the reads carry the two indices."""
        return "".join(part[0] for part in self.layout.orientation)

    @property
    def slot_order_arg(self) -> str:
        return self.layout.slot_order

    @property
    def extra_args(self) -> str:
        """``-p <read1 length>``, which splitBarcode needs only when ``-b`` is
        counted across the concatenated read."""
        if self.layout.is_pe and self.pe_offset_base == "concatenated":
            return f"-p {self.layout.read1_len}"
        return ""


def lane_plan(
    ss: Samplesheet,
    lane,
    run_dir: str,
    run: str,
    mismatch: int = 1,
    pe_offset_base: str = "read2-relative",
    n_reads: int = DEFAULT_N_READS,
    min_hit_rate: float = DEFAULT_MIN_HIT_RATE,
    slot_order: str = AUTO,
) -> LanePlan:
    """Resolve one lane end to end: reads, BioInfo geometry, orientation, ``-b``.

    Per lane, deliberately: two lanes of one run can carry different index
    lengths, so the sheet's index lengths -- which decide how many bases each
    ``-b`` takes -- are resolved from that lane's rows against that lane's
    ``BioInfo.csv``.
    """
    label = lane_label(lane_number(lane))
    reads = read_files(run_dir, run, label)
    bio = read_bioinfo(bioinfo_path(run_dir, label))

    if bio.is_pe != reads.is_pe:
        declares = "paired-end" if bio.is_pe else "single-end"
        holds = "a paired-end pair of" if reads.is_pe else "a single"
        raise BarcodeLayoutError(
            f"lane {label}: {bio.path} declares {declares} but the run directory "
            f"holds {holds} fastq(s)"
        )

    layout = measure_layout(
        reads.barcode_fastq,
        bio,
        index_pairs(ss, label),
        n_reads=n_reads,
        min_hit_rate=min_hit_rate,
        slot_order=slot_order,
    )
    if reads.is_pe:
        layout = layout._replace(read1_len=read_length(reads.r1))

    return LanePlan(
        lane=label,
        entries=barcode_entries(
            ss, label, orientation=layout.orientation, slot_order=layout.slot_order
        ),
        reads=reads,
        bio=bio,
        layout=layout,
        params=barcode_params(layout, mismatch, pe_offset_base),
        pe_offset_base=pe_offset_base,
    )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _write(path: Optional[str], text: str) -> None:
    if path in (None, "-"):
        sys.stdout.write(text)
        return
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w") as handle:
        handle.write(text)


def _cmd_lanes(args: argparse.Namespace) -> int:
    ss = read_samplesheet(args.sample_sheet)
    _write(args.out, "".join(lane + "\n" for lane in lanes_in_samplesheet(ss)))
    return 0


ORIENTATION_ARGS = {
    "ff": (FORWARD, FORWARD),
    "fr": (FORWARD, REVCOMP),
    "rf": (REVCOMP, FORWARD),
    "rr": (REVCOMP, REVCOMP),
}


def _cmd_barcodes(args: argparse.Namespace) -> int:
    ss = read_samplesheet(args.sample_sheet)
    entries = barcode_entries(
        ss,
        args.lane,
        orientation=ORIENTATION_ARGS[args.orientation],
        slot_order=args.slot_order,
    )
    _write(args.out, "".join(f"{sample}\t{code}\n" for sample, code in entries))
    return 0


def _cmd_params(args: argparse.Namespace) -> int:
    ss = read_samplesheet(args.sample_sheet)
    plan = lane_plan(
        ss,
        args.lane,
        args.run_dir,
        args.run or os.path.basename(os.path.normpath(args.run_dir)),
        mismatch=args.mismatch,
        pe_offset_base=args.pe_offset_base,
        n_reads=args.n_reads,
        min_hit_rate=args.min_hit_rate,
        slot_order=args.slot_order,
    )
    sys.stderr.write(f"{plan.lane}: {plan.layout.describe()}\n")
    sys.stderr.write(
        f"{plan.lane}: orientation arg {plan.orientation_arg}, slot order {plan.slot_order_arg}\n"
    )
    _write(args.out, (plan.b_args + " " + plan.extra_args).strip() + "\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Process an MGI sample sheet into splitBarcode inputs.",
    )
    sub = parser.add_subparsers(dest="command")
    sub.required = True  # argparse in 3.6 does not accept required= in add_subparsers

    lanes = sub.add_parser("lanes", help="list lanes present in the sheet, one per line")
    lanes.add_argument("--sample-sheet", required=True)
    lanes.add_argument("--out", default="-")
    lanes.set_defaults(func=_cmd_lanes)

    barcodes = sub.add_parser("barcodes", help="write the splitBarcode -B file for one lane")
    barcodes.add_argument("--sample-sheet", required=True)
    barcodes.add_argument("--lane", required=True)
    barcodes.add_argument("--out", default="-")
    barcodes.add_argument(
        "--orientation",
        default="ff",
        choices=sorted(ORIENTATION_ARGS),
        help="how the reads carry (index, index2): f=forward, r=reverse complement",
    )
    barcodes.add_argument(
        "--slot-order",
        default=INDEX_FIRST,
        choices=list(SLOT_ORDERS),
        help="which index the instrument sequenced first",
    )
    barcodes.set_defaults(func=_cmd_barcodes)

    params = sub.add_parser("params", help="measure the layout and emit the -b arguments")
    params.add_argument("--sample-sheet", required=True)
    params.add_argument("--lane", required=True)
    params.add_argument("--run-dir", required=True)
    params.add_argument("--run", default=None, help="defaults to the run-dir basename")
    params.add_argument("--out", default="-")
    params.add_argument("--mismatch", type=int, default=1)
    params.add_argument("--pe-offset-base", default="read2-relative", choices=list(PE_OFFSET_BASES))
    params.add_argument("--n-reads", type=int, default=DEFAULT_N_READS)
    params.add_argument("--min-hit-rate", type=float, default=DEFAULT_MIN_HIT_RATE)
    params.add_argument(
        "--slot-order",
        default=AUTO,
        choices=[AUTO] + list(SLOT_ORDERS),
        help="which index the instrument sequenced first; 'auto' measures it "
        "and refuses when the plate's index pairs are symmetric",
    )
    params.set_defaults(func=_cmd_params)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (SamplesheetError, BarcodeLayoutError, BioInfoError) as error:
        sys.stderr.write(f"ERROR: {error}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
