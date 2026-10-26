"""Verifier for the FFPE tumour / contaminated-normal somatic analysis task.

Expected values are the parameters and somatic SNVs planted when the synthetic data were
generated (truth by construction), not output of the reference solution. Contract checks
run first, so a failure names the violated field before any accuracy comparison.
"""
import json
import math
import os
import re
import stat
from pathlib import Path

SNVS = Path("/root/results/somatic_snvs.tsv")
ESTIMATES = Path("/root/results/sample_estimates.json")
TRUTH_SNVS = Path("/tests/data/truth_somatic.tsv")
TRUTH_PARAMS = Path("/tests/data/truth_params.json")
REFERENCE = Path("/tests/data/reference.fa")

MAX_BYTES = 1024 * 1024
MAX_ROWS = 5000
HEADER = "chrom\tpos\tref\talt"
CHROMS = ("chr1", "chr2", "chr3", "chr4")
BASES = ("A", "C", "G", "T")
POS_RE = re.compile(r"[1-9][0-9]*")
KEYS = {"tumor_purity", "tumor_ploidy", "normal_tumor_fraction"}

MIN_PRECISION = 0.90
MIN_RECALL = 0.85
PURITY_REL_TOL = 0.10
PLOIDY_REL_TOL = 0.10
NORMAL_FRACTION_REL_TOL = 0.15


def read_text(path):
    """Read a regular UTF-8 file of at most 1 MiB without blocking on FIFOs."""
    assert path.exists(), f"missing output {path}"
    fd = os.open(path, os.O_RDONLY | os.O_NONBLOCK)
    try:
        info = os.fstat(fd)
        assert stat.S_ISREG(info.st_mode), f"{path} must be a regular file"
        assert info.st_size <= MAX_BYTES, f"{path} must not exceed 1 MiB"
        data = os.read(fd, MAX_BYTES + 1)
    finally:
        os.close(fd)
    assert len(data) <= MAX_BYTES, f"{path} must not exceed 1 MiB"
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise AssertionError(f"{path} is not valid UTF-8: {exc}") from exc


def reference():
    seqs, name = {}, None
    for line in REFERENCE.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def snv_lines():
    text = read_text(SNVS)
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    lines = [ln[:-1] if ln.endswith("\r") else ln for ln in lines]
    assert lines, "somatic_snvs.tsv is empty"
    return lines


def reported_snvs():
    lines = snv_lines()
    assert lines[0] == HEADER, f"header must be exactly {HEADER!r}, got {lines[0]!r}"
    ref = reference()
    calls = []
    seen = set()
    for n, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        assert len(fields) == 4, f"line {n}: expected 4 tab-separated fields, got {len(fields)}"
        chrom, pos, rbase, abase = fields
        assert chrom in CHROMS, f"line {n}: chrom {chrom!r} is not one of {CHROMS}"
        assert POS_RE.fullmatch(pos), f"line {n}: pos {pos!r} is not a plain positive decimal integer"
        p = int(pos)
        assert p <= len(ref[chrom]), f"line {n}: pos {p} beyond the end of {chrom}"
        assert rbase in BASES, f"line {n}: ref {rbase!r} is not one of {BASES}"
        assert rbase == ref[chrom][p - 1], f"line {n}: ref {rbase!r} differs from reference base {ref[chrom][p - 1]!r}"
        assert abase in BASES, f"line {n}: alt {abase!r} is not one of {BASES}"
        assert abase != rbase, f"line {n}: alt equals ref"
        assert (chrom, p) not in seen, f"line {n}: duplicate position {chrom}:{p}"
        seen.add((chrom, p))
        calls.append((chrom, p, rbase, abase))
    assert 1 <= len(calls) <= MAX_ROWS, f"expected between 1 and {MAX_ROWS} SNV rows, got {len(calls)}"
    return set(calls)


def reject_constant(name):
    raise ValueError(f"non-finite JSON constant {name}")


def estimates():
    text = read_text(ESTIMATES)
    try:
        obj = json.loads(text, parse_constant=reject_constant)
    except ValueError as exc:
        raise AssertionError(f"sample_estimates.json is not valid finite JSON: {exc}") from exc
    assert isinstance(obj, dict), "sample_estimates.json must contain a JSON object"
    assert set(obj) == KEYS, f"keys must be exactly {sorted(KEYS)}, got {sorted(obj)}"
    for key in sorted(KEYS):
        v = obj[key]
        assert isinstance(v, (int, float)) and not isinstance(v, bool), f"{key} must be a JSON number"
        assert math.isfinite(v), f"{key} must be finite"
    assert 0 < obj["tumor_purity"] <= 1, "tumor_purity must be in (0, 1]"
    assert 0 < obj["tumor_ploidy"] <= 10, "tumor_ploidy must be in (0, 10]"
    assert 0 <= obj["normal_tumor_fraction"] < 1, "normal_tumor_fraction must be in [0, 1)"
    return obj


def truth_snvs():
    rows = TRUTH_SNVS.read_text().splitlines()
    out = set()
    for line in rows[1:]:
        chrom, pos, rbase, abase = line.split("\t")[:4]
        out.add((chrom, int(pos), rbase, abase))
    return out


def truth_params():
    return json.loads(TRUTH_PARAMS.read_text())


def relative_error(reported, true):
    return abs(reported - true) / abs(true)


def test_snv_table_contract():
    """somatic_snvs.tsv: regular UTF-8 file <= 1 MiB, exact header, 4 valid fields per row,
    ref matching reference.fa, alt != ref, unique positions, 1-5000 rows."""
    reported_snvs()


def test_estimates_contract():
    """sample_estimates.json: a JSON object with exactly the three keys, finite numeric values
    in their stated ranges."""
    estimates()


def test_snv_precision():
    """At least 90 % of reported SNVs match a planted somatic SNV on chrom, pos, ref and alt."""
    calls = reported_snvs()
    truth = truth_snvs()
    tp = len(calls & truth)
    precision = tp / len(calls)
    print(f"reported={len(calls)} true_positives={tp} precision={precision:.4f}")
    assert precision >= MIN_PRECISION, f"SNV precision {precision:.4f} < {MIN_PRECISION}"


def test_snv_recall():
    """At least 85 % of planted somatic SNVs (clonal and subclonal) are reported exactly."""
    calls = reported_snvs()
    truth = truth_snvs()
    tp = len(calls & truth)
    recall = tp / len(truth)
    print(f"true_positives={tp} truth={len(truth)} recall={recall:.4f}")
    assert recall >= MIN_RECALL, f"SNV recall {recall:.4f} < {MIN_RECALL}"


def test_tumor_purity():
    """tumor_purity is within 10 % relative error of the simulated malignant-cell fraction."""
    err = relative_error(estimates()["tumor_purity"], truth_params()["tumor_purity"])
    print(f"tumor_purity relative error {err:.4f}")
    assert err <= PURITY_REL_TOL, f"tumor_purity relative error {err:.4f} > {PURITY_REL_TOL}"


def test_tumor_ploidy():
    """tumor_ploidy is within 10 % relative error of the simulated length-weighted mean copy number."""
    err = relative_error(estimates()["tumor_ploidy"], truth_params()["tumor_ploidy"])
    print(f"tumor_ploidy relative error {err:.4f}")
    assert err <= PLOIDY_REL_TOL, f"tumor_ploidy relative error {err:.4f} > {PLOIDY_REL_TOL}"


def test_normal_tumor_fraction():
    """normal_tumor_fraction is within 15 % relative error of the simulated malignant-cell
    fraction of the blood sample."""
    err = relative_error(estimates()["normal_tumor_fraction"], truth_params()["normal_tumor_fraction"])
    print(f"normal_tumor_fraction relative error {err:.4f}")
    assert err <= NORMAL_FRACTION_REL_TOL, f"normal_tumor_fraction relative error {err:.4f} > {NORMAL_FRACTION_REL_TOL}"
