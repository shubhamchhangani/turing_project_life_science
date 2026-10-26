"""Generate the ffpe-somatic-calling dataset (author tool, never shipped to the agent).

Usage (inside the dev image):
    python3 generate.py SEED OUT_PUBLIC OUT_PRIVATE

OUT_PUBLIC  receives reference.fa(.fai), tumor.cram(.crai), normal.cram(.crai)
OUT_PRIVATE receives truth_somatic.tsv, truth_params.json and author diagnostics.

Physical model (all choices documented in the task README):
  * diploid synthetic genome (4 chromosomes, CpG depletion); germline het/hom SNVs
  * tumour: whole-genome-doubled malignant clone with clonal allele-specific copy-number
    segments; purity P; founder clone + nested subclones A (CCF_A) and B (CCF_B, B in A)
  * normal (PBMC): fraction T of cells are circulating malignant cells, same genome
  * DNA fragments are drawn in proportion to DNA content: in a segment with tumour copy
    number CN a fragment comes from a malignant cell with probability c*CN/(c*CN+2(1-c))
  * clonal SNVs carry multiplicity m (copies of the mutated allele): pre-doubling
    mutations on a doubled haplotype have m = n_hap, later ones m = 1; subclonal m = 1
  * FFPE cytosine deamination on single original strands of tumour molecules: ref C>T only
    in read pairs derived from the top strand (F1R2), ref G>A only in bottom-strand pairs (F2R1)
  * PCR duplicate families (identical coordinates and orientation, shared damage,
    independent sequencing errors); no duplicate flags set
  * sequencing errors from binned Illumina qualities that degrade along the read
  * recurrent homopolymer read-through errors at fixed sites in BOTH libraries (per read)
"""
import json
import subprocess
import sys
from array import array
from pathlib import Path

import numpy as np
import pysam

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
COMP = np.array([3, 2, 1, 0], dtype=np.uint8)  # A<->T, C<->G on 0..3 codes
READ_LEN = 150

CONTIGS = [("chr1", 300_000), ("chr2", 250_000), ("chr3", 250_000), ("chr4", 200_000)]
# clonal tumour copy number per segment: (contig, start, end, copies of hap0, copies of hap1)
SEGMENTS = [
    ("chr1", 0, 300_000, 2, 2),
    ("chr2", 0, 250_000, 2, 1),
    ("chr3", 0, 120_000, 2, 0),
    ("chr3", 120_000, 250_000, 2, 2),
    ("chr4", 0, 200_000, 3, 2),
]

PARAMS = dict(
    purity=0.62,
    normal_tumor_fraction=0.20,
    ccf_a=0.35,
    ccf_b=0.22,
    n_clonal=230,
    n_sub_a=110,
    n_sub_b=60,
    pre_doubling_fraction=0.5,
    het_rate=1 / 1100,
    hom_rate=1 / 3500,
    n_recurrent=260,
    recurrent_molecule_fraction=0.6,
    tumor_mean_unique_depth=170,
    normal_mean_unique_depth=90,
    tumor_family_q=0.31,   # geometric PCR family parameter (mean family 1/(1-q))
    normal_family_q=0.10,
    tumor_frag=(210.0, 0.24),   # lognormal median, sigma
    normal_frag=(320.0, 0.20),
    damage_mean=0.010,
    damage_shape=0.45,
)

QBINS = np.array([2, 12, 23, 37], dtype=np.uint8)


def make_genome(rng, n):
    freqs = np.array([0.295, 0.205, 0.205, 0.295])
    trans = np.tile(freqs, (4, 1))
    trans[1, 2] *= 0.25  # CpG depletion
    trans /= trans.sum(axis=1, keepdims=True)
    cum = np.cumsum(trans, axis=1)
    u = rng.random(n)
    seq = np.empty(n, dtype=np.uint8)
    seq[0] = 0
    for i in range(1, n):  # one-off cost, ~1 s per Mb
        seq[i] = np.searchsorted(cum[seq[i - 1]], u[i])
    return np.minimum(seq, 3)


def transition(b):
    return {0: 2, 2: 0, 1: 3, 3: 1}[int(b)]


def pick_alt(rng, ref, ts_prob=0.67):
    if rng.random() < ts_prob:
        return transition(ref)
    others = [x for x in range(4) if x != ref and x != transition(ref)]
    return int(rng.choice(others))


def homopolymer_followers(seq, min_run=5):
    """Positions immediately after a homopolymer run of >= min_run (and the run base)."""
    out = []
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
        else:
            if run >= min_run:
                out.append((i, int(seq[i - 1])))
            run = 1
    return out


def main(seed, out_pub, out_priv):
    p = PARAMS
    rng = np.random.default_rng(seed)
    out_pub.mkdir(parents=True, exist_ok=True)
    out_priv.mkdir(parents=True, exist_ok=True)

    # global coordinates: contigs concatenated
    offsets = {}
    o = 0
    for name, ln in CONTIGS:
        offsets[name] = o
        o += ln
    L = o
    contig_of = np.zeros(L, dtype=np.int64)
    for i, (name, ln) in enumerate(CONTIGS):
        contig_of[offsets[name]:offsets[name] + ln] = i
    ref = make_genome(rng, L)

    seg_n = np.zeros((2, L), dtype=np.int64)  # tumour copies of hap0 / hap1 per position
    segs = []
    for name, s, e, a, b in SEGMENTS:
        gs, ge = offsets[name] + s, offsets[name] + e
        seg_n[0, gs:ge], seg_n[1, gs:ge] = a, b
        segs.append((gs, ge, a, b))
    seg_cn = seg_n.sum(0)
    ploidy = float(seg_cn.mean())

    margin = 1000
    used = np.zeros(L, dtype=bool)
    for name, ln in CONTIGS:
        used[offsets[name]:offsets[name] + margin] = True
        used[offsets[name] + ln - margin:offsets[name] + ln] = True
    for gs, ge, _, _ in segs:  # keep variants away from breakpoints
        used[max(0, gs - 600):gs + 600] = True

    def reserve(pos):
        used[max(0, pos - 3):pos + 4] = True

    # recurrent read-through artifact sites (both libraries)
    cands = [(i, b) for i, b in homopolymer_followers(ref) if not used[i] and ref[i] != b
             and contig_of[i] == contig_of[i - 6]]
    idx = rng.choice(len(cands), size=p["n_recurrent"], replace=False)
    recurrent = []
    for k in sorted(idx):
        pos, base = cands[k]
        kind = "molecule" if rng.random() < p["recurrent_molecule_fraction"] else "read"
        recurrent.append(dict(pos=int(pos), alt=int(base), kind=kind, rate=float(rng.uniform(0.03, 0.16)),
                              jit_t=float(rng.uniform(0.8, 1.25)), jit_n=float(rng.uniform(0.8, 1.25))))
        reserve(pos)

    # germline variants on two haplotypes
    hap = np.stack([ref.copy(), ref.copy()])
    germline = []
    free = np.flatnonzero(~used)
    for kind, rate in (("het", p["het_rate"]), ("hom", p["hom_rate"])):
        n = rng.binomial(len(free), rate)
        for pos in rng.choice(free, size=n, replace=False):
            if used[pos]:
                continue
            alt = pick_alt(rng, ref[pos])
            h = int(rng.integers(2))
            if kind == "hom":
                hap[:, pos] = alt
            else:
                hap[h, pos] = alt
            germline.append(dict(pos=int(pos), ref=int(ref[pos]), alt=alt, kind=kind, hap=h if kind == "het" else -1))
            reserve(pos)

    # somatic SNVs: CpG C>T (SBS1-like) enriched in subclones, otherwise flat over 6 classes
    is_cpg = np.zeros(L, dtype=bool)
    is_cpg[:-1] |= (ref[:-1] == 1) & (ref[1:] == 2)
    is_cpg[1:] |= (ref[1:] == 2) & (ref[:-1] == 1)
    cpg_sites = np.flatnonzero(is_cpg)
    somatic = []
    for clone, n, cpg_frac in (("clonal", p["n_clonal"], 0.30), ("A", p["n_sub_a"], 0.60), ("B", p["n_sub_b"], 0.60)):
        made = 0
        while made < n:
            if rng.random() < cpg_frac:
                pos = int(rng.choice(cpg_sites))
                alt = transition(ref[pos])  # C>T or G>A at CpG
            else:
                pos = int(rng.integers(0, L))
                cls = int(rng.integers(6))  # C>A C>G C>T T>A T>C T>G (pyrimidine frame)
                pyr_ref, pyr_alt = [(1, 0), (1, 2), (1, 3), (3, 0), (3, 1), (3, 2)][cls]
                if ref[pos] in (1, 3):
                    if ref[pos] != pyr_ref:
                        continue
                    alt = pyr_alt
                else:
                    if COMP[ref[pos]] != pyr_ref:
                        continue
                    alt = int(COMP[pyr_alt])
            if used[pos]:
                continue
            haps = [h for h in (0, 1) if seg_n[h, pos] > 0]
            h = int(rng.choice(haps))
            nh = int(seg_n[h, pos])
            m = nh if (clone == "clonal" and nh >= 2 and rng.random() < p["pre_doubling_fraction"]) else 1
            somatic.append(dict(pos=pos, ref=int(ref[pos]), alt=int(alt), clone=clone, hap=h, mult=m))
            reserve(pos)
            made += 1
    somatic.sort(key=lambda s: s["pos"])
    s_pos = np.array([s["pos"] for s in somatic])
    s_alt = np.array([s["alt"] for s in somatic], dtype=np.uint8)
    s_hap = np.array([s["hap"] for s in somatic])
    s_lvl = np.array([{"clonal": 0, "A": 1, "B": 2}[s["clone"]] for s in somatic])
    s_frac = np.array([s["mult"] / seg_n[s["hap"], s["pos"]] for s in somatic])

    # FFPE damage propensities at every C (top strand) and G (bottom strand) position
    scale = p["damage_mean"] / p["damage_shape"]
    lam = np.minimum(rng.gamma(p["damage_shape"], scale, size=L), 0.5)
    c_sites = np.flatnonzero(ref == 1)
    g_sites = np.flatnonzero(ref == 2)

    # per-read sequencing errors (differ between mates) versus library-preparation errors
    # fixed into the molecule (shared by both mates and by every PCR duplicate)
    rec_read = [r for r in recurrent if r["kind"] == "read"]
    rec_mol = [r for r in recurrent if r["kind"] == "molecule"]
    rec_pos = np.array([r["pos"] for r in rec_read])
    rec_alt = np.array([r["alt"] for r in rec_read], dtype=np.uint8)
    mol_pos = np.array([r["pos"] for r in rec_mol])
    mol_alt = np.array([r["alt"] for r in rec_mol], dtype=np.uint8)

    with open(out_pub / "reference.fa", "w") as fh:
        for name, ln in CONTIGS:
            fh.write(f">{name}\n")
            s = BASES[ref[offsets[name]:offsets[name] + ln]].tobytes().decode()
            for i in range(0, ln, 60):
                fh.write(s[i:i + 60] + "\n")
    pysam.faidx(str(out_pub / "reference.fa"))

    def qualities(n_reads):
        cyc = np.arange(READ_LEN) / READ_LEN
        p2 = 0.002 + 0.012 * cyc ** 3
        p12 = 0.01 + 0.04 * cyc ** 2
        p23 = 0.04 + 0.10 * cyc ** 2
        cum = np.stack([p2, p2 + p12, p2 + p12 + p23], axis=1)
        u = rng.random((n_reads, READ_LEN))
        q = (u[..., None] > cum[None]).sum(-1)
        return QBINS[q]

    ERR = {2: 0.5, 12: 10 ** -1.2, 23: 10 ** -2.3, 37: 10 ** -3.7}
    err_lut = np.zeros(64)
    for k, v in ERR.items():
        err_lut[k] = v

    def library(name, c, mean_depth, fam_q, frag, damage, rg_id, out_bam):
        # fragments per segment in proportion to DNA content c*CN + 2(1-c)
        content = np.array([(ge - gs) * (c * (a + b) + 2 * (1 - c)) for gs, ge, a, b in segs])
        n_total = int(L * mean_depth / (2 * READ_LEN))
        per_seg = np.round(n_total * content / content.sum()).astype(int)
        med, sig = frag
        starts, flens = [], []
        for (gs, ge, _, _), nf in zip(segs, per_seg):
            fl = np.clip(np.round(med * np.exp(sig * rng.standard_normal(nf))), 160, 600).astype(int)
            # fragments stay inside their contig; segment of a fragment = segment of its start
            cid = contig_of[gs]
            cend = offsets[CONTIGS[cid][0]] + CONTIGS[cid][1]
            st = rng.integers(gs, ge, nf)
            st = np.minimum(st, cend - fl)
            starts.append(st)
            flens.append(fl)
        start = np.concatenate(starts)
        flen = np.concatenate(flens)
        n_frag = len(start)
        cn = seg_cn[start]
        tumor_cell = rng.random(n_frag) < c * cn / (c * cn + 2 * (1 - c))
        h0 = seg_n[0, start] / np.maximum(cn, 1)
        hapi = np.where(tumor_cell, (rng.random(n_frag) >= h0).astype(int), rng.integers(0, 2, n_frag))
        top = rng.random(n_frag) < 0.5           # original strand: top -> F1R2, bottom -> F2R1
        v = rng.random(n_frag)
        level = np.where(v < p["ccf_b"], 2, np.where(v < p["ccf_a"], 1, 0))
        fam = rng.geometric(1 - fam_q, n_frag)
        header = {"HD": {"VN": "1.6", "SO": "unsorted"},
                  "SQ": [{"SN": nm, "LN": ln} for nm, ln in CONTIGS],
                  "RG": [{"ID": rg_id, "SM": name, "LB": f"{name}-lib1", "PL": "ILLUMINA", "PU": f"HXK7{seed}.{rg_id}"}]}
        bam = pysam.AlignmentFile(str(out_bam), "wb", header=header)
        qual_pool = qualities(200_000)
        qp = 0
        tile = 0
        order = rng.permutation(n_frag)
        for f in order:
            s0, ln = int(start[f]), int(flen[f])
            cid = int(contig_of[s0])
            coff = offsets[CONTIGS[cid][0]]
            mol = hap[hapi[f], s0:s0 + ln].copy()
            if tumor_cell[f]:
                a, b = np.searchsorted(s_pos, [s0, s0 + ln])
                for k in range(a, b):
                    if s_hap[k] == hapi[f] and s_lvl[k] <= level[f] and rng.random() < s_frac[k]:
                        mol[s_pos[k] - s0] = s_alt[k]
            if damage:
                sites = c_sites if top[f] else g_sites
                a, b = np.searchsorted(sites, [s0, s0 + ln])
                if b > a:
                    sp = sites[a:b]
                    hit = sp[rng.random(b - a) < lam[sp]]
                    # deamination only where the molecule still carries the reference base
                    hit = hit[mol[hit - s0] == ref[hit]]
                    mol[hit - s0] = 3 if top[f] else 0  # C->T on top / G->A (C->T on bottom)
            a, b = np.searchsorted(mol_pos, [s0, s0 + ln])
            for k in range(a, b):
                r = rec_mol[k]
                if rng.random() < r["rate"] * (r["jit_t"] if damage else r["jit_n"]):
                    mol[mol_pos[k] - s0] = mol_alt[k]
            r1_fwd = bool(top[f])
            for copy in range(int(fam[f])):
                tile += 1
                qname = f"HXK7{seed}:{rg_id}:{1101 + tile % 97}:{tile}"
                reads = []
                for mate in (0, 1):
                    fwd = (mate == 0) == r1_fwd
                    rs = 0 if fwd else ln - READ_LEN
                    seq = mol[rs:rs + READ_LEN].copy()
                    if qp + 1 >= len(qual_pool):
                        qual_pool = qualities(200_000)
                        qp = 0
                    q = qual_pool[qp].copy()
                    qp += 1
                    if not fwd:  # cycles run from the 3' end of the forward coordinate
                        q = q[::-1]
                    # recurrent read-through errors (sequence-direction dependent)
                    ra, rb = np.searchsorted(rec_pos, [s0 + rs, s0 + rs + READ_LEN])
                    for k in range(ra, rb):
                        r = rec_read[k]
                        rate = r["rate"] * (r["jit_t"] if damage else r["jit_n"]) * (1.0 if fwd else 0.6)
                        if rng.random() < rate:
                            seq[rec_pos[k] - s0 - rs] = rec_alt[k]
                    e = rng.random(READ_LEN) < err_lut[q]
                    if e.any():
                        seq[e] = (seq[e] + rng.integers(1, 4, e.sum())) % 4
                    reads.append((fwd, s0 + rs - coff, seq, q))
                for mate, (fwd, pos, seq, q) in enumerate(reads):
                    mfwd, mpos = reads[1 - mate][0], reads[1 - mate][1]
                    a = pysam.AlignedSegment(bam.header)
                    a.query_name = qname
                    a.flag = 1 | 2 | (64 if mate == 0 else 128) | (0 if fwd else 16) | (0 if mfwd else 32)
                    a.reference_id = cid
                    a.reference_start = pos
                    a.mapping_quality = 60
                    a.cigarstring = f"{READ_LEN}M"
                    a.next_reference_id = cid
                    a.next_reference_start = mpos
                    a.template_length = ln if fwd else -ln
                    a.query_sequence = BASES[seq].tobytes().decode()
                    a.query_qualities = array("B", q.tobytes())
                    a.set_tags([("MC", f"{READ_LEN}M"), ("MQ", 60), ("RG", rg_id)])
                    bam.write(a)
        bam.close()
        return dict(fragments=int(n_frag), reads=int(2 * fam.sum()))

    stats = {}
    for name, frac, depth, q, frag, dmg, rg in (
        ("tumor", p["purity"], p["tumor_mean_unique_depth"], p["tumor_family_q"], p["tumor_frag"], True, "T1"),
        ("normal", p["normal_tumor_fraction"], p["normal_mean_unique_depth"], p["normal_family_q"], p["normal_frag"], False, "N1"),
    ):
        raw = out_priv / f"{name}.unsorted.bam"
        stats[name] = library(name, frac, depth, q, frag, dmg, rg, raw)
        sorted_bam = out_priv / f"{name}.sorted.bam"
        subprocess.run(["samtools", "sort", "-@2", "-m", "600M", "-o", str(sorted_bam), str(raw)], check=True)
        subprocess.run(["samtools", "view", "-C", "-T", str(out_pub / "reference.fa"),
                        "--output-fmt-option", "version=3.0",
                        "-o", str(out_pub / f"{name}.cram"), str(sorted_bam)], check=True)
        subprocess.run(["samtools", "index", str(out_pub / f"{name}.cram")], check=True)
        raw.unlink()
        sorted_bam.unlink()

    b = "ACGT"
    names = [nm for nm, _ in CONTIGS]
    lvl_ccf = {"clonal": 1.0, "A": p["ccf_a"], "B": p["ccf_b"]}

    def local(gpos):
        cid = int(contig_of[gpos])
        return names[cid], gpos - offsets[names[cid]] + 1

    P_, T_ = p["purity"], p["normal_tumor_fraction"]
    with open(out_priv / "truth_somatic.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tclone\tccf\tmultiplicity\ttumor_cn\texpected_vaf_tumor\texpected_vaf_normal\n")
        for s in somatic:
            ccf = lvl_ccf[s["clone"]]
            cn = int(seg_cn[s["pos"]])
            ch, lp = local(s["pos"])
            fh.write(f"{ch}\t{lp}\t{b[s['ref']]}\t{b[s['alt']]}\t{s['clone']}\t{ccf}\t{s['mult']}\t{cn}\t"
                     f"{P_ * ccf * s['mult'] / (P_ * cn + 2 * (1 - P_)):.5f}\t"
                     f"{T_ * ccf * s['mult'] / (T_ * cn + 2 * (1 - T_)):.5f}\n")
    with open(out_priv / "truth_germline.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tkind\n")
        for g in sorted(germline, key=lambda g: g["pos"]):
            ch, lp = local(g["pos"])
            fh.write(f"{ch}\t{lp}\t{b[g['ref']]}\t{b[g['alt']]}\t{g['kind']}\n")
    with open(out_priv / "truth_recurrent.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tkind\trate\n")
        for r in recurrent:
            ch, lp = local(r["pos"])
            fh.write(f"{ch}\t{lp}\t{b[ref[r['pos']]]}\t{b[r['alt']]}\t{r['kind']}\t{r['rate']:.4f}\n")
    np.save(out_priv / "damage_lambda.npy", lam.astype(np.float32))
    (out_priv / "truth_params.json").write_text(json.dumps(dict(
        seed=seed, tumor_purity=P_, tumor_ploidy=ploidy, normal_tumor_fraction=T_,
        n_somatic=len(somatic), segments=SEGMENTS, library_stats=stats, params=p), indent=2) + "\n")


if __name__ == "__main__":
    main(int(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
