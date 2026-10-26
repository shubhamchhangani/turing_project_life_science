"""Reference analysis: FFPE tumour versus a tumour-contaminated normal.

Usage: python3 call_somatic.py DATA_DIR RESULTS_DIR

Pipeline
  1. Fragment-level counting: duplicate families collapsed by (contig, left end, right end,
     pair orientation); overlapping mates counted once; BQ >= 20, MAPQ >= 20. Counts are kept
     separately for F1R2 and F2R1 pairs, and mate-discordant overlap observations are recorded.
  2. Allele-specific copy number: tumour depth in 5 kb bins plus B-allele fractions of germline
     heterozygous sites, optimal-partition segmentation, then a grid fit of purity and ploidy
     with integer allele-specific copy numbers per segment (ASCAT-style), refined by requiring
     clonal SNV multiplicities to be integers.
  3. Damage-aware SNV evidence: ref C>T is damage-prone in F1R2 pairs and ref G>A in F2R1
     pairs; those orientations enter the likelihood marginalised over a per-site deamination
     rate whose prior is fitted (nonparametric EM) on background C/G sites.
  4. Tumour-in-normal fraction T: two-class mixture (somatic with copy-number-aware expected
     normal allele fraction vs library-independent recurrent error) maximised over T.
  5. Calls: somatic class posterior > 0.5, tumour evidence over noise, not germline, and no
     mate-discordant support (per-read errors are discordant, molecular changes concordant).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pysam
from scipy.special import gammaln

import os

MIN_BQ = 20
DEBUG = bool(os.environ.get("SOMATIC_DEBUG"))
MIN_MQ = 20
BIN = 5000
B2I = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    B2I[_b] = _i


def count_library(path, ref_path, offsets, length):
    """Return (counts[orientation, position, base], discordant[position, base]).

    Positions are global (contigs concatenated in reference order). counts are over
    deduplicated fragments; orientation 0 = F1R2 (read 1 forward), 1 = F2R1 (read 1 reverse).
    discordant counts each mate base at overlap positions where the two mates disagree.
    """
    counts = np.zeros((2, length, 4), dtype=np.int32)
    disc = np.zeros(length * 4, dtype=np.int32)
    seen = set()
    pending = {}
    buf_idx, disc_idx = [], []
    n_buf = 0

    def flush():
        nonlocal buf_idx, disc_idx, n_buf
        if buf_idx:
            counts.reshape(-1)[:] += np.bincount(np.concatenate(buf_idx), minlength=counts.size).astype(np.int32)
        if disc_idx:
            disc[:] += np.bincount(np.concatenate(disc_idx), minlength=disc.size).astype(np.int32)
        buf_idx, disc_idx, n_buf = [], [], 0

    def bases(read, off):
        seq = np.frombuffer(read.query_sequence.encode(), dtype=np.uint8)
        qual = np.asarray(read.query_qualities, dtype=np.int64)
        if read.cigartuples == [(0, read.query_length)]:
            qpos = np.arange(read.query_length)
            rpos = read.reference_start + qpos
        else:
            arr = np.array(read.get_aligned_pairs(matches_only=True), dtype=np.int64).reshape(-1, 2)
            qpos, rpos = arr[:, 0], arr[:, 1]
        b = B2I[seq[qpos]]
        ok = (qual[qpos] >= MIN_BQ) & (b < 4)
        return rpos[ok] + off, b[ok]

    with pysam.AlignmentFile(path, "rc", reference_filename=ref_path) as fh:
        names = list(fh.references)
        for read in fh.fetch(until_eof=True):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail
                    or read.is_duplicate or read.mapping_quality < MIN_MQ or not read.is_paired
                    or read.mate_is_unmapped or read.reference_id != read.next_reference_id):
                continue
            mate = pending.pop(read.query_name, None)
            if mate is None:
                pending[read.query_name] = read
                continue
            r1, r2 = (read, mate) if read.is_read1 else (mate, read)
            orient = 1 if r1.is_reverse else 0
            off = offsets[names[r1.reference_id]]
            key = (r1.reference_id, min(r1.reference_start, r2.reference_start),
                   max(r1.reference_end, r2.reference_end), orient)
            if key in seen:  # duplicate family member: keep the first fragment only
                continue
            seen.add(key)
            p1, b1 = bases(r1, off)
            p2, b2 = bases(r2, off)
            common, i1, i2 = np.intersect1d(p1, p2, assume_unique=True, return_indices=True)
            keep2 = np.ones(len(p2), dtype=bool)
            keep2[i2] = False
            disagree = b1[i1] != b2[i2]
            keep1 = np.ones(len(p1), dtype=bool)
            keep1[i1[disagree]] = False
            if disagree.any():
                cp = common[disagree]
                disc_idx.append(np.concatenate([cp * 4 + b1[i1[disagree]], cp * 4 + b2[i2[disagree]]]))
            pos = np.concatenate([p1[keep1], p2[keep2]])
            base = np.concatenate([b1[keep1], b2[keep2]])
            buf_idx.append((orient * length + pos) * 4 + base)
            n_buf += len(pos)
            if n_buf > 20_000_000:
                flush()
    flush()
    return counts, disc.reshape(length, 4)


def logbinom(k, n, p):
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + k * np.log(p) + (n - k) * np.log1p(-p)


def segment_contig(x, penalty):
    """Optimal partitioning of rows of x (bins x features) under squared error + penalty."""
    n = len(x)
    cs = np.vstack([np.zeros(x.shape[1]), np.cumsum(x, 0)])
    cs2 = np.concatenate([[0.0], np.cumsum((x ** 2).sum(1))])

    def cost(i, j):  # bins i..j-1
        s = cs[j] - cs[i]
        return cs2[j] - cs2[i] - (s ** 2).sum() / (j - i)

    best = np.full(n + 1, np.inf)
    best[0] = -penalty
    back = np.zeros(n + 1, dtype=int)
    for j in range(1, n + 1):
        for i in range(j):
            c = best[i] + cost(i, j) + penalty
            if c < best[j]:
                best[j], back[j] = c, i
    cuts = []
    j = n
    while j > 0:
        cuts.append((back[j], j))
        j = back[j]
    return cuts[::-1]


def minor_fraction(alt, dep):
    """MLE of the minor-allele fraction for unphased hets (symmetric binomial mixture)."""
    grid = np.linspace(0.0, 0.5, 501)
    ll = np.logaddexp(logbinom(alt[:, None], dep[:, None], grid + 1e-4),
                      logbinom(alt[:, None], dep[:, None], 1 - grid - 1e-4)).sum(0)
    return float(grid[ll.argmax()])


def main(data, results):
    ref_path = str(data / "reference.fa")
    fasta = pysam.FastaFile(ref_path)
    names = list(fasta.references)
    lens = [fasta.get_reference_length(c) for c in names]
    starts = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(int)
    offsets = dict(zip(names, starts))
    length = int(sum(lens))
    ref = np.concatenate([B2I[np.frombuffer(fasta.fetch(c).upper().encode(), dtype=np.uint8)] for c in names])

    def local(g):
        ci = int(np.searchsorted(starts, g, side="right") - 1)
        return names[ci], int(g - starts[ci])

    tum, tum_disc = count_library(str(data / "tumor.cram"), ref_path, offsets, length)
    nor, nor_disc = count_library(str(data / "normal.cram"), ref_path, offsets, length)
    t_all = tum.sum(0)
    n_all = nor.sum(0)
    t_dep = t_all.sum(1)
    n_dep = n_all.sum(1)
    valid = ref < 4
    pos_idx = np.arange(length)

    # ---------------- allele-specific copy number, purity and ploidy ----------------
    n_alt_counts = n_all.copy()
    n_alt_counts[pos_idx[valid], ref[valid]] = -1
    n_alt = n_alt_counts.argmax(1)
    n_vaf = n_all[pos_idx, n_alt] / np.maximum(n_dep, 1)
    het = np.flatnonzero(valid & (n_dep >= 30) & (n_vaf > 0.3) & (n_vaf < 0.7) & (t_dep >= 40))
    het_alt = t_all[het, n_alt[het]].astype(float)
    het_dep = t_dep[het].astype(float)

    bins = []  # (contig index, global start, global end)
    for ci, ln in enumerate(lens):
        for s in range(0, ln, BIN):
            bins.append((ci, starts[ci] + s, starts[ci] + min(ln, s + BIN)))
    bins = np.array(bins)
    depth_bin = np.array([t_dep[s:e].mean() for _, s, e in bins])
    het_bin = np.searchsorted(bins[:, 1], het, side="right") - 1
    baf_bin = np.array([np.mean(np.abs(het_alt[het_bin == k] / het_dep[het_bin == k] - 0.5))
                        if (het_bin == k).sum() >= 2 else np.nan for k in range(len(bins))])
    mean_depth = float(t_dep[valid].mean())
    lr_bin = depth_bin / mean_depth
    sd_lr = np.median(np.abs(np.diff(lr_bin))) / 0.954 + 1e-9
    fin = np.isfinite(baf_bin)
    sd_baf = np.median(np.abs(np.diff(baf_bin[fin]))) / 0.954 + 1e-9
    baf_fill = np.where(fin, baf_bin, np.nanmedian(baf_bin))
    feats = np.stack([lr_bin / sd_lr, baf_fill / sd_baf], 1)

    segments = []  # (global start, global end)
    for ci in range(len(names)):
        rows = np.flatnonzero(bins[:, 0] == ci)
        for i, j in segment_contig(feats[rows], penalty=40.0):
            segments.append((int(bins[rows[i], 1]), int(bins[rows[j - 1], 2])))
    seg_len = np.array([e - s for s, e in segments], dtype=float)
    seg_r = np.array([t_dep[s:e].mean() / mean_depth for s, e in segments])
    seg_b = []
    for s, e in segments:
        m = (het >= s) & (het < e)
        seg_b.append(minor_fraction(het_alt[m], het_dep[m]) if m.sum() >= 5 else np.nan)
    seg_b = np.array(seg_b)

    # integer allele-specific states (major, minor) up to total 8
    states = np.array([(a, b) for a in range(0, 9) for b in range(0, a + 1) if 0 < a + b <= 8])
    st_n = states.sum(1).astype(float)
    sig_r = 0.02
    sig_b = 0.015

    def fit(P, psi):
        R = (P * st_n + 2 * (1 - P)) / (P * psi + 2 * (1 - P))
        B = (P * states[:, 1] + (1 - P)) / (P * st_n + 2 * (1 - P))
        cost_r = ((seg_r[:, None] - R[None]) / sig_r) ** 2
        cost_b = np.where(np.isfinite(seg_b)[:, None], ((np.nan_to_num(seg_b)[:, None] - B[None]) / sig_b) ** 2, 0.0)
        c = cost_r + cost_b
        k = c.argmin(1)
        implied = float((seg_len * st_n[k]).sum() / seg_len.sum())
        total = float((seg_len / seg_len.sum() * c[np.arange(len(k)), k]).sum()) + ((implied - psi) / 0.02) ** 2
        return total, k, implied

    fits = []
    for P in np.arange(0.10, 1.0001, 0.01):
        for psi in np.arange(1.5, 6.0001, 0.02):
            total, k, implied = fit(P, psi)
            fits.append((total, P, psi, k, implied))
    fits.sort(key=lambda f: f[0])

    # ---------------- SNV candidates and damage-aware evidence ----------------
    alt_counts = t_all.copy()
    alt_counts[pos_idx[valid], ref[valid]] = -1
    alt = alt_counts.argmax(1)
    a_t = t_all[pos_idx, alt]
    cand = np.flatnonzero(valid & (a_t >= 3) & (a_t >= 0.015 * np.maximum(t_dep, 1)))
    r = ref[cand]
    al = alt[cand]
    prone = np.where((r == 1) & (al == 3), 0, np.where((r == 2) & (al == 0), 1, -1))
    clean = np.where(prone >= 0, 1 - prone, -1)
    at_alt = np.where(clean >= 0, tum[np.maximum(clean, 0), cand, al], a_t[cand])
    at_dep = np.where(clean >= 0, tum[np.maximum(clean, 0), cand].sum(1), t_dep[cand])
    an_alt = n_all[cand, al]
    an_dep = n_dep[cand]
    vn = an_alt / np.maximum(an_dep, 1)
    seg_starts = np.array([s for s, _ in segments])
    cand_seg = np.searchsorted(seg_starts, cand, side="right") - 1

    mask = valid.copy()
    mask[cand] = False
    nonref = n_dep[mask].sum() - n_all[np.flatnonzero(mask), ref[mask]].sum()
    err = float(max(nonref / max(n_dep[mask].sum(), 1) / 3, 1e-5))

    lam = np.concatenate([[0.0], np.geomspace(1e-3, 0.5, 47)])
    bg = []
    for base, alt_b, o in ((1, 3, 0), (2, 0, 1)):
        sites = np.flatnonzero(ref == base)
        ok = (tum[1 - o, sites, alt_b] == 0) & (n_all[sites, alt_b] == 0)
        sites = sites[ok]
        bg.append(np.stack([tum[o, sites, alt_b], tum[o, sites].sum(1)], 1))
    bg = np.concatenate(bg).astype(float)
    bg = bg[np.random.default_rng(1).permutation(len(bg))[:200_000]]
    lb = logbinom(bg[:, :1], bg[:, 1:], np.clip(lam[None] + err, 1e-9, 1))
    wl = np.full(len(lam), 1 / len(lam))
    for _ in range(300):
        z = lb + np.log(wl + 1e-300)
        z = np.exp(z - z.max(1, keepdims=True))
        z /= z.sum(1, keepdims=True)
        wl = z.mean(0)
    log_wl = np.log(wl + 1e-300)

    vg = np.linspace(0.004, 0.7, 175)
    log_pv = -np.log(len(vg))
    p_alt = np.where(prone >= 0, tum[np.maximum(prone, 0), cand, al], 0)
    p_dep = np.where(prone >= 0, tum[np.maximum(prone, 0), cand].sum(1), 0)
    lt = np.empty((len(cand), len(vg)))
    lt0 = np.empty(len(cand))
    for lo in range(0, len(cand), 2000):
        sl = slice(lo, lo + 2000)
        ac, dc = at_alt[sl, None].astype(float), at_dep[sl, None].astype(float)
        ap, dp = p_alt[sl, None].astype(float), p_dep[sl, None].astype(float)
        has_p = (prone[sl] >= 0)[:, None]

        def tumour_ll(v):
            lc = logbinom(ac, dc, v + err)
            pr = np.logaddexp.reduce(
                logbinom(ap[..., None], dp[..., None], np.clip(v[..., None] + lam + err, 1e-9, 1 - 1e-9)) + log_wl, axis=-1)
            return lc + np.where(has_p, pr, 0.0)

        lt[sl] = tumour_ll(np.broadcast_to(vg, (ac.shape[0], len(vg))))
        lt0[sl] = tumour_ll(np.zeros((ac.shape[0], 1)))[:, 0]
    an = an_alt[:, None].astype(float)
    dn = an_dep[:, None].astype(float)
    ll_0 = lt0 + logbinom(an_alt, an_dep, err)
    germline = (vn >= 0.3) & (an_alt >= 3)
    tumour_evidence = np.logaddexp.reduce(lt + log_pv, axis=1) - lt0

    # Distinct copy-number solutions (different integer state assignments) are compared by
    # how close clonal SNV multiplicities VAF * (P n + 2(1-P)) / P sit to integers; within the
    # chosen solution the purity/ploidy optimum of the copy-number fit is kept.
    strong = np.flatnonzero(~germline & (tumour_evidence > 12.0) & (at_dep >= 40))
    vt_strong = at_alt[strong] / at_dep[strong]
    families = {}
    for total, P, psi, k, implied in fits:
        if total > fits[0][0] + 25:
            break
        key = tuple(k)
        if key not in families:  # fits are sorted, so the first entry is the family optimum
            families[key] = (total, P, psi, k, implied)
    best_fit = None
    for total, P, psi, k, implied in families.values():
        n_site = st_n[k][cand_seg[strong]]
        mult = vt_strong * (P * n_site + 2 * (1 - P)) / P
        big = mult > 0.75  # clonal-like sites
        if big.sum() < 20:
            continue
        resid = float(np.mean(np.minimum(np.abs(mult[big] - np.round(mult[big])), 0.5) ** 2))
        score = total + 400 * resid
        if DEBUG:
            print(f"fit P={P:.2f} psi={psi:.2f} implied={implied:.3f} cn_cost={total:.3f} mult_resid={resid:.4f} score={score:.3f}")
        if best_fit is None or score < best_fit[0]:
            best_fit = (score, P, psi, k, implied)
    _, purity, _, seg_state, ploidy = best_fit
    seg_n = st_n[seg_state]

    # recurrent library-independent error: same rate in both libraries up to a lognormal factor
    jg = np.exp(np.linspace(-0.75, 0.75, 9))
    log_pj = -0.5 * (np.log(jg) / 0.3) ** 2
    log_pj -= np.logaddexp.reduce(log_pj)
    ll_r = np.empty(len(cand))
    for lo in range(0, len(cand), 2000):
        sl = slice(lo, lo + 2000)
        rec = lt[sl, :, None] + logbinom(an[sl, :, None], dn[sl, :, None], np.clip(vg[None, :, None] * jg + err, 1e-9, 0.99)) + log_pj
        ll_r[sl] = np.logaddexp.reduce(rec.reshape(rec.shape[0], -1), axis=1) + log_pv

    n_c = seg_n[cand_seg]

    def ll_somatic(T, rows=slice(None)):
        # tumour VAF v -> mutant copies per malignant cell x -> expected normal VAF
        x = vg[None] * (purity * n_c[rows, None] + 2 * (1 - purity)) / purity
        pn = T * x / (T * n_c[rows, None] + 2 * (1 - T))
        return np.logaddexp.reduce(lt[rows] + logbinom(an[rows], dn[rows], np.clip(pn, 0, 0.99) + err) + log_pv, axis=1)

    sel = np.flatnonzero(~germline & (tumour_evidence > 8.0))
    best = None
    for T in np.arange(0.0, 0.8001, 0.0025):
        ls = ll_somatic(T, sel)
        w_ = 0.5
        for _ in range(100):
            z = 1 / (1 + np.exp(np.clip(np.log(1 - w_) + ll_r[sel] - np.log(w_) - ls, -50, 50)))
            w_ = float(np.clip(z.mean(), 1e-3, 1 - 1e-3))
        tot = float(np.logaddexp(np.log(w_) + ls, np.log(1 - w_) + ll_r[sel]).sum())
        if best is None or tot > best[0]:
            best = (tot, T, w_)
    _, tin, w_s = best

    ll_s = ll_somatic(tin)
    # sequence context: allele equal to the base of an adjacent homopolymer run (>= 4) on either
    # side, a known systematic error mode. Its frequency in each class is learned by EM.
    def run_len(pos, step, base):
        n = 0
        p_ = pos + step
        while 0 <= p_ < length and ref[p_] == base and n < 8:
            n += 1
            p_ += step
        return n
    hp = np.array([max(run_len(c, -1, b), run_len(c, 1, b)) >= 4 for c, b in zip(cand, al)])
    ev = np.flatnonzero(~germline & (ll_s - ll_0 > 3.0))
    q_s, q_r = 0.05, 0.5
    for _ in range(200):
        ls_ = np.log(w_s) + ll_s[ev] + np.where(hp[ev], np.log(q_s), np.log1p(-q_s))
        lr_ = np.log(1 - w_s) + ll_r[ev] + np.where(hp[ev], np.log(q_r), np.log1p(-q_r))
        z = 1 / (1 + np.exp(np.clip(lr_ - ls_, -50, 50)))
        w_s = float(np.clip(z.mean(), 1e-3, 1 - 1e-3))
        q_s = float(np.clip((z * hp[ev]).sum() / z.sum(), 1e-3, 1 - 1e-3))
        q_r = float(np.clip(((1 - z) * hp[ev]).sum() / (1 - z).sum(), 1e-3, 1 - 1e-3))
    ctx_s = np.where(hp, np.log(q_s), np.log1p(-q_s))
    ctx_r = np.where(hp, np.log(q_r), np.log1p(-q_r))
    post = (np.log(w_s) + ll_s + ctx_s) - np.logaddexp(np.log(w_s) + ll_s + ctx_s, np.log(1 - w_s) + ll_r + ctx_r)
    d_alt = tum_disc[cand, al] + nor_disc[cand, al]
    d_exp = 2 * err * (tum_disc[cand].sum(1) + nor_disc[cand].sum(1) + 1)
    read_error_site = (d_alt >= 2) & (d_alt > 4 * d_exp)
    somatic = ~germline & ~read_error_site & (ll_s - ll_0 > 3.0) & (post > np.log(0.5))

    results.mkdir(parents=True, exist_ok=True)
    bases = "ACGT"
    with open(results / "somatic_snvs.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\n")
        for i in np.flatnonzero(somatic):
            c, lp = local(cand[i])
            fh.write(f"{c}\t{lp + 1}\t{bases[r[i]]}\t{bases[al[i]]}\n")
    (results / "sample_estimates.json").write_text(json.dumps(
        {"tumor_purity": round(float(purity), 4), "tumor_ploidy": round(float(ploidy), 4),
         "normal_tumor_fraction": round(float(tin), 4)}, indent=2) + "\n")
    desc = []
    for (s, e), k, rr, bb in zip(segments, seg_state, seg_r, seg_b):
        c, ls_ = local(s)
        desc.append(f"{c}:{ls_}-{ls_ + e - s}={states[k][0]}+{states[k][1]}(r={rr:.3f},b={bb:.3f})")
    print("segments:", ", ".join(desc))
    print(f"hets={len(het)} err={err:.2e} candidates={len(cand)} calls={int(somatic.sum())} "
          f"read_error_sites={int((read_error_site & ~germline & (ll_s - ll_0 > 3)).sum())} "
          f"purity={purity:.3f} ploidy={ploidy:.3f} tin={tin:.4f} w_s={w_s:.2f} q_s={q_s:.3f} q_r={q_r:.3f}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
