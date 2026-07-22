# 2026 Benjamin J Perry
# MIT License
# Copyright (c) 2026 Benjamin J Perry
"""Unit tests for the MGI sample-sheet / barcode-layout logic.

Run under the snakemake module's Python 3.11::

    module load snakemake/7.32.3-foss-2023a-Python-3.11.6
    python3 -m unittest discover tests

Fixtures are *generated* in ``setUp`` rather than committed, so the repository
carries no binary ``.fq.gz`` blobs and the shapes under test (read length,
index lengths, trailing junk bases) stay visible and editable as code.
"""

import gzip
import os
import random
import sys
import tempfile
import unittest

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "workflow", "scripts")
)

import process_samplesheet as pss  # noqa: E402


# --------------------------------------------------------------------------
# Fixture builders
# --------------------------------------------------------------------------

# Real barcodes from the validated T1+ run (DL100018749 sheet, L01/L03).
T1_L01 = [
    ("SQ5522", "GAACTGAGCG", "AGCGCTAG"),
    ("SQ5523", "AGGTCAGATA", "GATATCGA"),
    ("SQ5524", "CGTCTCATAT", "CGCAGACG"),
]
T1_L03 = [
    ("SQ5565", "GAACTGAGCG", "AGCGCTAG"),
    ("SQ5568", "ATTCCATAAG", "TATGAGTA"),
    ("SQ5569", "GACGAGATTA", "ACATAGCG"),
]
# G99 shape: 8 bp index, 10 bp index2, index2 shared across all samples.
G99_L01 = [
    ("SSWAG2212002", "ATCACGTT", "TGTCCGTAGG"),
    ("SSWAG2212101", "GTAGAAGC", "TGTCCGTAGG"),
    ("GTSGPC01", "GATCAGCG", "TGTCCGTAGG"),
]


def write_samplesheet(path, rows, lanes="12", crlf=False, with_index2=True):
    """Write a minimal but structurally faithful MGI sheet.

    ``rows`` items are ``(sample_id, index, index2)`` -- placed on ``lanes`` --
    or ``(lanes, sample_id, index, index2)`` to mix lane blocks in one sheet,
    as the real T1+ sheet does with its 12 and 34 blocks.
    """
    columns = ["Lanes", "Sample_ID", "I7_Index_ID", "index"]
    if with_index2:
        columns += ["I5_Index_ID", "index2"]
    columns += ["Sample_Project"]

    lines = [
        "[Header]",
        "Date,7/07/2026",
        "Instrument Type,T1+",
        "Flowcell, TESTRUN",
        "",
        "[Reads]",
        "100",
        "",
        "[Data]",
        ",".join(columns),
    ]
    for row in rows:
        row_lanes, sample_id, index, index2 = row if len(row) == 4 else (lanes,) + tuple(row)
        fields = [row_lanes, sample_id, "UDP0001", index]
        if with_index2:
            fields += ["Q5001", index2]
        fields += ["proj"]
        lines.append(",".join(fields))
    lines += ["", "[GenerateKeyfile]", "Sample_ID,plateid", "SQ0001,13751"]

    terminator = "\r\n" if crlf else "\n"
    with open(path, "w", newline="") as handle:
        handle.write(terminator.join(lines) + terminator)
    return path


def _random_bases(length, rng):
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_plain_fastq(path, length=150, n_reads=400, seed=7):
    """A fastq with no barcode anywhere -- read 1 of a paired-end MGI run."""
    rng = random.Random(seed)
    with gzip.open(path, "wt") as handle:
        for number in range(n_reads):
            sequence = _random_bases(length, rng)
            handle.write("@read{}\n{}\n+\n{}\n".format(number, sequence, "F" * length))
    return path


def write_fastq(
    path,
    rows,
    insert_len=100,
    block1_len=None,
    block2_len=None,
    n_reads=400,
    revcomp_index2=False,
    revcomp_index=False,
    scramble=False,
    seed=1,
):
    """Build a gzipped fastq whose reads carry the barcodes at the end.

    Reads are laid out exactly as the instrument writes them::

        insert | block 1 | block 2

    where a block is the sheet's index followed by however many cycles the
    machine sequenced but the sheet does not use. Those unused cycles are what
    ``block1_len``/``block2_len`` model: T1+ sequences 10 barcode cycles for
    ``DL100018466``'s 8-base index, which puts **2 junk bases between the two
    indices**, not just after them.
    """
    rng = random.Random(seed)

    with gzip.open(path, "wt") as handle:
        for number in range(n_reads):
            sample_id, index, index2 = rows[number % len(rows)]
            if scramble:
                # No sample's barcode appears anywhere in the read.
                index = _random_bases(len(index), rng)
                index2 = _random_bases(len(index2), rng)
            pad1 = _random_bases((block1_len or len(index)) - len(index), rng)
            pad2 = _random_bases((block2_len or len(index2)) - len(index2), rng)
            if revcomp_index:
                index = pss.revcomp(index)
            if revcomp_index2:
                index2 = pss.revcomp(index2)
            insert = _random_bases(insert_len, rng)
            sequence = insert + index + pad1 + index2 + pad2
            handle.write(
                "@read{}\n{}\n+\n{}\n".format(number, sequence, "F" * len(sequence))
            )
    return path


def write_bioinfo(path, platform="T1+", read1_len=100, read2_len=0, barcode_len=10, dual_len=8):
    """Write a BioInfo.csv in the schema the named platform uses.

    The two schemas differ in more than key names: G99 states its block order
    explicitly and overstates ``Read1 Cycles`` by one, which is exactly the
    behaviour the anchor derivation has to survive.
    """
    if platform == "G99":
        lines = [
            "Machine ID,INV-MGI-G99",
            "Sequence Type,SE{}+{}+{}-GTSEQ".format(read1_len, dual_len, barcode_len),
            "Sequence Order,Read1-Read2-Dualbarcode-Barcode",
            "Read1 Cycles,{}".format(read1_len + 1),  # G99 really does report 101 for 100
            "Read2 Cycles,{}".format(read2_len),
            "Barcode,{}".format(barcode_len),
            "Dual Barcode,{}".format(dual_len),
        ]
    else:
        lines = [
            "Machine ID,INV-MGI-T1",
            "Read1Len,{}".format(read1_len),
            "Read2Len,{}".format(read2_len),
            "BarcodeLen,{}".format(barcode_len),
            "DualBarcodeLen,{}".format(dual_len),
            "BarcodePosition,BarcodePosEnd",
        ]
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


class FixtureCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

    def lane_dir(self, run, lane):
        path = os.path.join(self.tmp, run, lane)
        os.makedirs(path, exist_ok=True)
        return path

    def se_run(self, run, lane, rows, platform="T1+", bioinfo=True, **kwargs):
        """A single-end run directory; returns the run dir."""
        directory = self.lane_dir(run, lane)
        write_fastq(
            os.path.join(directory, "{}_{}_read.fq.gz".format(run, lane)), rows, **kwargs
        )
        if bioinfo:
            write_bioinfo(
                os.path.join(directory, "BioInfo.csv"),
                platform=platform,
                read1_len=kwargs.get("insert_len", 100),
                read2_len=0,
                barcode_len=kwargs.get("block1_len") or len(rows[0][1]),
                dual_len=kwargs.get("block2_len") or len(rows[0][2]),
            )
        return os.path.join(self.tmp, run)

    def pe_run(self, run, lane, rows, read1_len=150, bioinfo=True, **kwargs):
        """A paired-end run directory: read 1 is barcode-free, read 2 carries it."""
        directory = self.lane_dir(run, lane)
        write_plain_fastq(
            os.path.join(directory, "{}_{}_read_1.fq.gz".format(run, lane)), length=read1_len
        )
        write_fastq(
            os.path.join(directory, "{}_{}_read_2.fq.gz".format(run, lane)), rows, **kwargs
        )
        if bioinfo:
            write_bioinfo(
                os.path.join(directory, "BioInfo.csv"),
                read1_len=read1_len,
                read2_len=kwargs.get("insert_len", 150),
                barcode_len=kwargs.get("block1_len") or len(rows[0][1]),
                dual_len=kwargs.get("block2_len") or len(rows[0][2]),
            )
        return os.path.join(self.tmp, run)

    def bio(self, **kwargs):
        """A parsed BioInfo for the geometry/orientation tests."""
        path = os.path.join(self.tmp, "BioInfo.{}.csv".format(len(os.listdir(self.tmp))))
        return pss.read_bioinfo(write_bioinfo(path, **kwargs))


# --------------------------------------------------------------------------
# Sheet parsing and lane expansion
# --------------------------------------------------------------------------


class TestParseLanes(unittest.TestCase):
    def test_digit_string_expands(self):
        self.assertEqual(pss.parse_lanes("1"), [1])
        self.assertEqual(pss.parse_lanes("12"), [1, 2])
        self.assertEqual(pss.parse_lanes("34"), [3, 4])
        self.assertEqual(pss.parse_lanes("1234"), [1, 2, 3, 4])

    def test_whitespace_tolerated(self):
        self.assertEqual(pss.parse_lanes(" 12 "), [1, 2])

    def test_non_digit_refused(self):
        for value in ("1,2", "1-2", "L01", "", "  ", "one"):
            with self.assertRaises(pss.SamplesheetError):
                pss.parse_lanes(value)

    def test_lane_zero_refused(self):
        with self.assertRaises(pss.SamplesheetError):
            pss.parse_lanes("012")

    def test_lane_number_accepts_labels(self):
        self.assertEqual(pss.lane_number("L01"), 1)
        self.assertEqual(pss.lane_number("L04"), 4)
        self.assertEqual(pss.lane_number(3), 3)
        self.assertEqual(pss.lane_label(1), "L01")


class TestSamplesheet(FixtureCase):
    def test_lanes_union_across_rows(self):
        """Two lane blocks (12 and 34), exactly as the real T1+ sheet is laid out."""
        rows = [("12",) + row for row in T1_L01] + [("34",) + row for row in T1_L03]
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), rows)
        ss = pss.read_samplesheet(path)
        self.assertEqual(pss.lanes_in_samplesheet(ss), ["L01", "L02", "L03", "L04"])
        self.assertEqual(
            [row["Sample_ID"] for row in pss.samples_for_lane(ss, "L03")],
            ["SQ5565", "SQ5568", "SQ5569"],
        )

    def test_membership_not_equality(self):
        """Lanes=12 must place a row in BOTH L01 and L02."""
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01, lanes="12")
        ss = pss.read_samplesheet(path)
        for lane in ("L01", "L02"):
            ids = [row["Sample_ID"] for row in pss.samples_for_lane(ss, lane)]
            self.assertEqual(ids, ["SQ5522", "SQ5523", "SQ5524"])
        with self.assertRaises(pss.SamplesheetError):
            pss.samples_for_lane(ss, "L03")

    def test_crlf_and_leading_spaces(self):
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01, crlf=True)
        ss = pss.read_samplesheet(path)
        self.assertEqual(ss.header["Flowcell"], "TESTRUN")
        self.assertEqual(pss.lanes_in_samplesheet(ss), ["L01", "L02"])
        self.assertEqual(pss.barcode_entries(ss, "L01")[0], ("SQ5522", "GAACTGAGCGAGCGCTAG"))

    def test_sections_preserved(self):
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01)
        ss = pss.read_samplesheet(path)
        self.assertIn("GenerateKeyfile", ss.sections)
        self.assertIn("Reads", ss.sections)

    def test_missing_file(self):
        with self.assertRaises(pss.SamplesheetError):
            pss.read_samplesheet(os.path.join(self.tmp, "nope.csv"))


class TestBarcodeEntries(FixtureCase):
    def test_entries_are_id_tab_index_plus_index2(self):
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01)
        ss = pss.read_samplesheet(path)
        self.assertEqual(
            pss.barcode_entries(ss, "L01"),
            [
                ("SQ5522", "GAACTGAGCGAGCGCTAG"),
                ("SQ5523", "AGGTCAGATAGATATCGA"),
                ("SQ5524", "CGTCTCATATCGCAGACG"),
            ],
        )

    def test_single_index_collapses(self):
        rows = [(s, i, "") for s, i, _ in T1_L01]
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), rows, with_index2=False)
        ss = pss.read_samplesheet(path)
        self.assertEqual(pss.barcode_entries(ss, "L01")[0], ("SQ5522", "GAACTGAGCG"))

    def test_mixed_index_lengths_refused(self):
        rows = list(T1_L01) + [("SQ9999", "ACGT", "AGCGCTAG")]
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), rows)
        ss = pss.read_samplesheet(path)
        with self.assertRaises(pss.SamplesheetError):
            pss.barcode_entries(ss, "L01")

    def test_duplicate_sample_refused(self):
        rows = list(T1_L01) + [T1_L01[0]]
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), rows)
        ss = pss.read_samplesheet(path)
        with self.assertRaises(pss.SamplesheetError):
            pss.barcode_entries(ss, "L01")

    def test_non_acgtn_index_refused(self):
        rows = [("SQ5522", "GAACTGAGC?", "AGCGCTAG")]
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), rows)
        ss = pss.read_samplesheet(path)
        with self.assertRaises(pss.SamplesheetError):
            pss.barcode_entries(ss, "L01")

    def test_write_barcode_file_round_trip(self):
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01)
        ss = pss.read_samplesheet(path)
        out = os.path.join(self.tmp, "nested", "L01.barcodes")
        pss.write_barcode_file(pss.barcode_entries(ss, "L01"), out)
        with open(out) as handle:
            self.assertEqual(handle.readline(), "SQ5522\tGAACTGAGCGAGCGCTAG\n")


# --------------------------------------------------------------------------
# Read discovery
# --------------------------------------------------------------------------


class TestReadFiles(FixtureCase):
    def test_single_end(self):
        run_dir = self.se_run("RUN1", "L01", T1_L01, block2_len=10)
        reads = pss.read_files(run_dir, "RUN1", "L01")
        self.assertFalse(reads.is_pe)
        self.assertEqual(reads.barcode_fastq, reads.r1)
        self.assertTrue(reads.cli_args.startswith("-1 "))
        self.assertNotIn("-2", reads.cli_args)

    def test_paired_end_prefers_pair_and_scans_read2(self):
        run_dir = self.pe_run("RUN2", "L01", T1_L01, insert_len=150)
        reads = pss.read_files(run_dir, "RUN2", "L01")
        self.assertTrue(reads.is_pe)
        self.assertEqual(reads.barcode_fastq, reads.r2)
        self.assertIn("-2", reads.cli_args)
        self.assertEqual(len(reads.paths), 2)

    def test_missing_reads_raise(self):
        self.lane_dir("RUN3", "L01")
        with self.assertRaises(pss.SamplesheetError):
            pss.read_files(os.path.join(self.tmp, "RUN3"), "RUN3", "L01")

    def test_half_a_pair_refused(self):
        directory = self.lane_dir("RUN4", "L01")
        write_fastq(os.path.join(directory, "RUN4_L01_read_1.fq.gz"), T1_L01)
        with self.assertRaises(pss.SamplesheetError):
            pss.read_files(os.path.join(self.tmp, "RUN4"), "RUN4", "L01")


# --------------------------------------------------------------------------
# BioInfo.csv
# --------------------------------------------------------------------------


class TestReadBioInfo(FixtureCase):
    def test_t1_schema(self):
        bio = self.bio(platform="T1+", read1_len=100, barcode_len=10, dual_len=10)
        self.assertEqual(bio.platform, "T1+")
        self.assertEqual(bio.read1_len, 100)
        self.assertFalse(bio.is_pe)
        self.assertEqual(
            [(block.name, block.length) for block in bio.blocks],
            [(pss.BARCODE, 10), (pss.DUAL_BARCODE, 10)],
        )
        self.assertEqual(bio.barcode_cycles, 20)

    def test_g99_schema_reverses_the_block_order(self):
        """``Sequence Order,Read1-Read2-Dualbarcode-Barcode`` is the whole
        'the platforms are inverted' story -- read it and one rule fits both."""
        bio = self.bio(platform="G99", read1_len=100, barcode_len=10, dual_len=8)
        self.assertEqual(bio.platform, "G99")
        self.assertEqual(
            [(block.name, block.length) for block in bio.blocks],
            [(pss.DUAL_BARCODE, 8), (pss.BARCODE, 10)],
        )
        self.assertEqual(bio.barcode_cycles, 18)

    def test_g99_declared_read1_is_one_too_long(self):
        bio = self.bio(platform="G99", read1_len=100, barcode_len=10, dual_len=8)
        self.assertEqual(bio.declared_insert_len, 101)

    def test_pe_is_recognised_from_read2_len(self):
        bio = self.bio(read1_len=150, read2_len=150)
        self.assertTrue(bio.is_pe)
        self.assertEqual(bio.declared_insert_len, 150)

    def test_zero_length_block_is_dropped(self):
        bio = self.bio(barcode_len=10, dual_len=0)
        self.assertEqual([block.length for block in bio.blocks], [10])

    def test_missing_file(self):
        with self.assertRaises(pss.BioInfoError):
            pss.read_bioinfo(os.path.join(self.tmp, "nope.csv"))

    def test_bioinfo_path_prefers_the_lane(self):
        directory = self.lane_dir("RUNB", "L01")
        write_bioinfo(os.path.join(directory, "BioInfo.csv"))
        self.assertEqual(
            pss.bioinfo_path(os.path.join(self.tmp, "RUNB"), "L01"),
            os.path.join(directory, "BioInfo.csv"),
        )

    def test_bioinfo_path_raises_when_absent(self):
        self.lane_dir("RUNC", "L01")
        with self.assertRaises(pss.BioInfoError):
            pss.bioinfo_path(os.path.join(self.tmp, "RUNC"), "L01")


# --------------------------------------------------------------------------
# Geometry: offsets from BioInfo, lengths from the sheet
# --------------------------------------------------------------------------


class TestBarcodeGeometry(FixtureCase):
    def test_t1_se_shape(self):
        """DL100018479: 120 bp read, blocks 10+10, sheet indices 10/8."""
        bio = self.bio(read1_len=100, barcode_len=10, dual_len=10)
        self.assertEqual(pss.barcode_geometry(bio, 120, 10, 8), (100, 110))

    def test_g99_se_shape(self):
        """FT150034703: 118 bp read, blocks 8+10, sheet indices 8/10.

        The anchor comes out at 100 even though BioInfo declares 101 cycles.
        """
        bio = self.bio(platform="G99", read1_len=100, barcode_len=10, dual_len=8)
        self.assertEqual(pss.barcode_geometry(bio, 118, 8, 10), (100, 108))

    def test_pe_shape_leaves_a_gap_between_the_indices(self):
        """DL100018466: blocks 10+10 but 8-base indices -> index2 at 160, not 158.

        This is the case that the old adjacency assumption got wrong: the two
        unused cycles sit *between* the indices, not after them.
        """
        bio = self.bio(read1_len=150, read2_len=150, barcode_len=10, dual_len=10)
        self.assertEqual(pss.barcode_geometry(bio, 170, 8, 8), (150, 160))

    def test_single_index_has_no_second_offset(self):
        bio = self.bio(barcode_len=10, dual_len=10)
        self.assertEqual(pss.barcode_geometry(bio, 120, 10, 0), (100, None))

    def test_index_longer_than_its_block_is_refused(self):
        bio = self.bio(barcode_len=8, dual_len=10)
        with self.assertRaises(pss.BarcodeLayoutError) as caught:
            pss.barcode_geometry(bio, 118, 10, 10)
        self.assertIn("first barcode block", str(caught.exception))

    def test_index2_longer_than_its_block_is_refused(self):
        bio = self.bio(barcode_len=10, dual_len=8)
        with self.assertRaises(pss.BarcodeLayoutError) as caught:
            pss.barcode_geometry(bio, 118, 10, 10)
        self.assertIn("second barcode block", str(caught.exception))

    def test_dual_index_sheet_against_single_barcode_run_is_refused(self):
        bio = self.bio(barcode_len=10, dual_len=0)
        with self.assertRaises(pss.BarcodeLayoutError):
            pss.barcode_geometry(bio, 110, 10, 8)

    def test_read_shorter_than_the_barcode_cycles_is_refused(self):
        bio = self.bio(barcode_len=10, dual_len=10)
        with self.assertRaises(pss.BarcodeLayoutError):
            pss.barcode_geometry(bio, 20, 10, 8)


# --------------------------------------------------------------------------
# Orientation: the one thing BioInfo cannot tell us
# --------------------------------------------------------------------------


class TestMeasureLayout(FixtureCase):
    def measure(self, path, rows, bio, **kwargs):
        return pss.measure_layout(path, bio, [(i, i2) for _, i, i2 in rows], **kwargs)

    def test_t1_se_shape(self):
        path = os.path.join(self.tmp, "t1.fq.gz")
        write_fastq(path, T1_L01, insert_len=100, block2_len=10)
        bio = self.bio(read1_len=100, barcode_len=10, dual_len=10)
        layout = self.measure(path, T1_L01, bio)
        self.assertEqual((layout.offset, layout.offset2), (100, 110))
        self.assertEqual((layout.index_len, layout.index2_len), (10, 8))
        self.assertEqual(layout.read_len, 120)
        self.assertEqual(layout.orientation, (pss.FORWARD, pss.FORWARD))
        self.assertGreater(layout.hit_rate, 0.9)
        self.assertFalse(layout.anchor_disagrees)

    def test_g99_se_shape_flags_the_declared_length(self):
        path = os.path.join(self.tmp, "g99.fq.gz")
        write_fastq(path, G99_L01, insert_len=100)
        bio = self.bio(platform="G99", read1_len=100, barcode_len=10, dual_len=8)
        layout = self.measure(path, G99_L01, bio)
        self.assertEqual((layout.offset, layout.offset2), (100, 108))
        self.assertEqual(layout.read_len, 118)
        self.assertTrue(layout.anchor_disagrees)
        self.assertTrue(layout.anchor_is_g99_fencepost)
        self.assertIn("one-cycle overstatement", layout.describe())
        self.assertNotIn("WARNING", layout.describe())

    def test_a_bigger_anchor_mismatch_is_warned_about_not_normalised(self):
        """Only the +1 is routine; anything else is a reason to stop and look."""
        path = os.path.join(self.tmp, "offanchor.fq.gz")
        write_fastq(path, G99_L01, insert_len=100)
        bio = self.bio(platform="G99", read1_len=104, barcode_len=10, dual_len=8)
        layout = self.measure(path, G99_L01, bio)
        self.assertEqual(layout.offset, 100)  # the reads still win
        self.assertFalse(layout.anchor_is_g99_fencepost)
        self.assertIn("WARNING", layout.describe())

    def test_pe_revcomp_shape_is_detected_not_refused(self):
        """The DL100018466 shape: gap between the indices, both revcomp."""
        path = os.path.join(self.tmp, "pe.fq.gz")
        rows = [(s, i[:8], i2[:8]) for s, i, i2 in T1_L01]
        write_fastq(
            path,
            rows,
            insert_len=150,
            block1_len=10,
            block2_len=10,
            revcomp_index=True,
            revcomp_index2=True,
        )
        bio = self.bio(read1_len=150, read2_len=150, barcode_len=10, dual_len=10)
        layout = self.measure(path, rows, bio)
        self.assertEqual((layout.offset, layout.offset2), (150, 160))
        self.assertEqual(layout.orientation, (pss.REVCOMP, pss.REVCOMP))
        self.assertGreater(layout.hit_rate, 0.9)

    def test_single_revcomp_index2_is_detected(self):
        path = os.path.join(self.tmp, "rc2.fq.gz")
        write_fastq(path, T1_L01, insert_len=100, block2_len=10, revcomp_index2=True)
        bio = self.bio(read1_len=100, barcode_len=10, dual_len=10)
        layout = self.measure(path, T1_L01, bio)
        self.assertEqual(layout.orientation, (pss.FORWARD, pss.REVCOMP))

    def test_single_index(self):
        rows = [(s, i, "") for s, i, _ in T1_L01]
        path = os.path.join(self.tmp, "single.fq.gz")
        write_fastq(path, rows, insert_len=100)
        bio = self.bio(read1_len=100, barcode_len=10, dual_len=0)
        layout = self.measure(path, rows, bio)
        self.assertEqual(layout.offset, 100)
        self.assertEqual(layout.index2_len, 0)
        self.assertIsNone(layout.offset2)
        self.assertEqual(pss.format_b_args(pss.barcode_params(layout)), "-b 100 10 1")

    def test_raises_when_no_barcodes_match(self):
        """A mismatched sheet/run pair must abort, not emit a plausible offset."""
        path = os.path.join(self.tmp, "scrambled.fq.gz")
        write_fastq(path, T1_L01, insert_len=100, block2_len=10, scramble=True)
        bio = self.bio(read1_len=100, barcode_len=10, dual_len=10)
        with self.assertRaises(pss.BarcodeLayoutError) as caught:
            self.measure(path, T1_L01, bio)
        self.assertIn("hit-rate", str(caught.exception))

    def test_index2_first_run_is_detected(self):
        """The instrument sequenced index2 into the first block.

        Detectable here because these three samples' pairs are *not* symmetric.
        """
        rows = [(s, i, i2) for s, i, i2 in G99_L01]
        path = os.path.join(self.tmp, "swapped.fq.gz")
        # Reads carry index2 (10 bp) first, then index (8 bp).
        write_fastq(path, [(s, i2, i) for s, i, i2 in rows], insert_len=100)
        bio = self.bio(read1_len=100, barcode_len=10, dual_len=8)
        layout = self.measure(path, rows, bio)
        self.assertEqual(layout.slot_order, pss.INDEX2_FIRST)
        self.assertEqual(layout.slot_lengths, (10, 8))
        self.assertGreater(layout.hit_rate, 0.9)

    def test_symmetric_pairs_refuse_to_guess_the_slot_order(self):
        """The DL100018466 degeneracy: (A,B) and (B,A) both present.

        Both orders decode identically, so picking one would silently swap the
        two samples. The scan must stop and ask instead.
        """
        rows = [("S1", "ACGTTTTT", "GGGGCCAT"), ("S2", "GGGGCCAT", "ACGTTTTT")]
        path = os.path.join(self.tmp, "symmetric.fq.gz")
        write_fastq(path, rows, insert_len=100)
        bio = self.bio(read1_len=100, barcode_len=8, dual_len=8)
        with self.assertRaises(pss.BarcodeLayoutError) as caught:
            self.measure(path, rows, bio)
        message = str(caught.exception)
        self.assertIn("silently swap samples", message)
        self.assertIn("BARCODE_SLOT_ORDER", message)

    def test_declared_slot_order_resolves_the_degeneracy(self):
        rows = [("S1", "ACGTTTTT", "GGGGCCAT"), ("S2", "GGGGCCAT", "ACGTTTTT")]
        path = os.path.join(self.tmp, "symmetric2.fq.gz")
        write_fastq(path, rows, insert_len=100)
        bio = self.bio(read1_len=100, barcode_len=8, dual_len=8)
        layout = self.measure(path, rows, bio, slot_order=pss.INDEX2_FIRST)
        self.assertEqual(layout.slot_order, pss.INDEX2_FIRST)
        self.assertTrue(layout.slot_order_declared)
        self.assertEqual(layout.symmetric_pairs, 2)
        self.assertIn("cannot be checked against the reads", layout.describe())

    def test_pairs_are_scored_as_pairs_not_independently(self):
        """A read whose two halves come from *different* samples is not a hit."""
        rows = [("S1", "ACGTACGT", "TTTTGGGG"), ("S2", "CCCCAAAA", "TTGGCCAA")]
        crossed = [("X", "ACGTACGT", "TTGGCCAA")]  # S1's index with S2's index2
        path = os.path.join(self.tmp, "crossed.fq.gz")
        write_fastq(path, crossed, insert_len=100)
        bio = self.bio(read1_len=100, barcode_len=8, dual_len=8)
        with self.assertRaises(pss.BarcodeLayoutError) as caught:
            self.measure(path, rows, bio)
        self.assertIn("hit-rate", str(caught.exception))

    def test_wrong_block_lengths_are_caught_by_the_hit_rate(self):
        """BioInfo describing a shape the reads do not have must abort.

        The geometry alone cannot detect this -- the offsets it produces are
        arithmetically fine, they just point at the wrong bases.
        """
        path = os.path.join(self.tmp, "shifted.fq.gz")
        write_fastq(path, T1_L01, insert_len=100, block2_len=10)
        bio = self.bio(read1_len=100, barcode_len=12, dual_len=8)
        with self.assertRaises(pss.BarcodeLayoutError):
            self.measure(path, T1_L01, bio)


# --------------------------------------------------------------------------
# -b argument construction
# --------------------------------------------------------------------------


class TestBarcodeParams(unittest.TestCase):
    def layout(self, **kwargs):
        base = dict(
            fastq="x.fq.gz",
            platform="T1+",
            offset=100,
            offset2=110,
            blocks=[pss.BarcodeBlock(pss.BARCODE, 10), pss.BarcodeBlock(pss.DUAL_BARCODE, 10)],
            index_len=10,
            index2_len=8,
            orientation=(pss.FORWARD, pss.FORWARD),
            hit_rate=0.9,
            runner_up_rate=0.0,
            n_scanned=2000,
            read_len=120,
            declared_insert_len=100,
            is_pe=False,
        )
        base.update(kwargs)
        return pss.BarcodeLayout(**base)

    def test_t1_production_arguments(self):
        self.assertEqual(
            pss.format_b_args(pss.barcode_params(self.layout())), "-b 100 10 1 -b 110 8 1"
        )

    def test_g99_production_arguments(self):
        layout = self.layout(index_len=8, index2_len=10, offset2=108, read_len=118)
        self.assertEqual(
            pss.format_b_args(pss.barcode_params(layout)), "-b 100 8 1 -b 108 10 1"
        )

    def test_mismatch_is_threaded_through(self):
        self.assertEqual(
            pss.format_b_args(pss.barcode_params(self.layout(), mismatch=2)),
            "-b 100 10 2 -b 110 8 2",
        )

    def test_b_length_comes_from_the_sheet_not_the_block(self):
        """The instrument sequenced 10 cycles; the sheet supplies 8, so -b takes 8.

        The unused cycles still shift index2's start -- that is the block's job,
        not the index's.
        """
        layout = self.layout(index_len=8, index2_len=8, offset2=160, offset=150)
        self.assertEqual(
            pss.format_b_args(pss.barcode_params(layout)), "-b 150 8 1 -b 160 8 1"
        )

    def test_pe_read2_relative_is_the_default(self):
        layout = self.layout(
            offset=150, offset2=160, index_len=8, index2_len=8, read1_len=150, read_len=170, is_pe=True
        )
        self.assertEqual(
            pss.format_b_args(pss.barcode_params(layout)), "-b 150 8 1 -b 160 8 1"
        )

    def test_pe_concatenated_shifts_by_read1_length(self):
        layout = self.layout(
            offset=150, offset2=160, index_len=8, index2_len=8, read1_len=150, read_len=170, is_pe=True
        )
        self.assertEqual(
            pss.format_b_args(pss.barcode_params(layout, pe_offset_base="concatenated")),
            "-b 300 8 1 -b 310 8 1",
        )

    def test_concatenated_does_not_shift_a_single_end_run(self):
        self.assertEqual(
            pss.format_b_args(pss.barcode_params(self.layout(), pe_offset_base="concatenated")),
            "-b 100 10 1 -b 110 8 1",
        )

    def test_unknown_pe_offset_base_refused(self):
        with self.assertRaises(ValueError):
            pss.barcode_params(self.layout(), pe_offset_base="whatever")


# --------------------------------------------------------------------------
# End-to-end lane plan
# --------------------------------------------------------------------------


class TestLanePlan(FixtureCase):
    def test_se_lane_plan(self):
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01, lanes="12")
        run_dir = self.se_run("RUNX", "L01", T1_L01, insert_len=100, block2_len=10)
        ss = pss.read_samplesheet(path)
        plan = pss.lane_plan(ss, "L01", run_dir, "RUNX")
        self.assertEqual(plan.lane, "L01")
        self.assertEqual(plan.b_args, "-b 100 10 1 -b 110 8 1")
        self.assertEqual(plan.orientation_arg, "ff")
        self.assertEqual(plan.extra_args, "")
        self.assertEqual(len(plan.entries), 3)
        self.assertFalse(plan.reads.is_pe)
        # Forward run: the barcode file is the sheet's own sequences.
        self.assertEqual(plan.entries[0], ("SQ5522", "GAACTGAGCGAGCGCTAG"))

    def test_pe_revcomp_lane_plan_writes_revcomp_into_the_barcode_file(self):
        """The DL100018466 shape end to end: gap, revcomp, 8-base -b lengths."""
        rows = [(s, i[:8], i2[:8]) for s, i, i2 in T1_L01]
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), rows, lanes="1")
        run_dir = self.pe_run(
            "RUNY",
            "L01",
            rows,
            read1_len=150,
            insert_len=150,
            block1_len=10,
            block2_len=10,
            revcomp_index=True,
            revcomp_index2=True,
        )
        ss = pss.read_samplesheet(path)
        plan = pss.lane_plan(ss, "L01", run_dir, "RUNY")
        self.assertTrue(plan.reads.is_pe)
        self.assertEqual(plan.layout.read1_len, 150)
        self.assertEqual(plan.b_args, "-b 150 8 1 -b 160 8 1")
        self.assertEqual(plan.orientation_arg, "rr")
        # The -b offsets stay forward; the revcomp lives in the -B file, so the
        # entry is revcomp(index) + revcomp(index2).
        sample, barcode = plan.entries[0]
        self.assertEqual(
            barcode, pss.revcomp(rows[0][1]) + pss.revcomp(rows[0][2])
        )

    def test_pe_extra_args_carry_read1_length_when_concatenated(self):
        rows = [(s, i[:8], i2[:8]) for s, i, i2 in T1_L01]
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), rows, lanes="1")
        run_dir = self.pe_run(
            "RUNP", "L01", rows, read1_len=150, insert_len=150, block1_len=10, block2_len=10
        )
        ss = pss.read_samplesheet(path)
        plan = pss.lane_plan(ss, "L01", run_dir, "RUNP", pe_offset_base="concatenated")
        self.assertEqual(plan.b_args, "-b 300 8 1 -b 310 8 1")
        self.assertEqual(plan.extra_args, "-p 150")

    def test_mismatched_sheet_and_run_aborts(self):
        """The fail-loud contract, end to end."""
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01, lanes="1")
        run_dir = self.se_run("RUNZ", "L01", T1_L01, insert_len=100, scramble=True)
        ss = pss.read_samplesheet(path)
        with self.assertRaises(pss.BarcodeLayoutError):
            pss.lane_plan(ss, "L01", run_dir, "RUNZ")

    def test_pe_sheet_against_se_run_aborts(self):
        """BioInfo and the run directory disagreeing is a stop, not a warning."""
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01, lanes="1")
        run_dir = self.se_run("RUNW", "L01", T1_L01, insert_len=100, block2_len=10)
        write_bioinfo(
            os.path.join(run_dir, "L01", "BioInfo.csv"),
            read1_len=100,
            read2_len=100,
            barcode_len=10,
            dual_len=10,
        )
        ss = pss.read_samplesheet(path)
        with self.assertRaises(pss.BarcodeLayoutError) as caught:
            pss.lane_plan(ss, "L01", run_dir, "RUNW")
        self.assertIn("paired-end", str(caught.exception))

    def test_missing_bioinfo_aborts(self):
        path = write_samplesheet(os.path.join(self.tmp, "ss.csv"), T1_L01, lanes="1")
        run_dir = self.se_run(
            "RUNV", "L01", T1_L01, insert_len=100, block2_len=10, bioinfo=False
        )
        ss = pss.read_samplesheet(path)
        with self.assertRaises(pss.BioInfoError):
            pss.lane_plan(ss, "L01", run_dir, "RUNV")


if __name__ == "__main__":
    unittest.main()
