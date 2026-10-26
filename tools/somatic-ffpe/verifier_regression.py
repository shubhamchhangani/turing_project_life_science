"""Build valid/invalid artifact cases and run the packaged verifier image on each.

python3 verifier_regression.py ORACLE_RESULTS_DIR PRIVATE_DIR OUT_DIR
Expect: every 'valid_*' case reward 1, every 'bad_*' case reward 0.
"""
import json, os, shutil, subprocess, sys
from pathlib import Path
src, priv, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
shutil.rmtree(out, ignore_errors=True)
snv = (src / "somatic_snvs.tsv").read_text()
est = json.loads((src / "sample_estimates.json").read_text())
lines = snv.rstrip("\n").split("\n")
truth = [l.split("\t")[:4] for l in (priv / "truth_somatic.tsv").read_text().splitlines()[1:]]
germ = [l.split("\t")[:4] for l in (priv / "truth_germline.tsv").read_text().splitlines()[1:]]
rec = [l.split("\t")[:4] for l in (priv / "truth_recurrent.tsv").read_text().splitlines()[1:]]
cases = {}
def case(name, tsv=snv, js=json.dumps(est), special=None):
    cases[name] = (tsv, js, special)
case("valid_oracle")
case("valid_crlf", tsv=snv.replace("\n", "\r\n"))
case("valid_no_final_newline", tsv=snv.rstrip("\n"))
case("valid_trailing_blank_lines", tsv=snv + "\n\n")
case("valid_integer_like_estimates", js='{"tumor_purity": 0.62, "tumor_ploidy": 4, "normal_tumor_fraction": 0.2}')
case("bad_missing_files", special="empty")
case("bad_header", tsv=snv.replace("chrom\tpos\tref\talt", "chrom\tposition\tref\talt", 1))
case("bad_header_commas", tsv=snv.replace("\t", ","))
case("bad_extra_column", tsv="\n".join(l + "\tPASS" for l in lines) + "\n")
case("bad_lowercase_alt", tsv="\n".join([lines[0]] + [l[:-1] + l[-1].lower() for l in lines[1:]]) + "\n")
case("bad_ref_mismatch", tsv="\n".join([lines[0], "\t".join(lines[1].split("\t")[:2] + ["A" if lines[1].split("\t")[2] != "A" else "C", "G" if lines[1].split("\t")[2] != "G" else "T"])] + lines[2:]) + "\n")
case("bad_leading_zero_pos", tsv="\n".join([lines[0], lines[1].replace("\t", "\t0", 1)] + lines[2:]) + "\n")
case("bad_space_in_pos", tsv="\n".join([lines[0], lines[1].replace("\t", "\t ", 1)] + lines[2:]) + "\n")
case("bad_unicode_digit_pos", tsv="\n".join([lines[0], lines[1].replace("\t", "\t١", 1)] + lines[2:]) + "\n")
case("bad_duplicate_position", tsv=snv + lines[1] + "\n")
case("bad_zero_rows", tsv=lines[0] + "\n")
case("bad_too_many_rows", tsv=lines[0] + "\n" + "".join(f"chr1\t{i}\tA\tC\n" for i in range(1, 6000)))
case("bad_unknown_chrom", tsv="\n".join([lines[0], "chr5\t10\tA\tC"] + lines[1:]) + "\n")
case("bad_blank_line", tsv=snv + "\n" + lines[1] + "\n")
case("bad_not_utf8", tsv=snv.encode().decode() + "\xff")
case("bad_nan_estimate", js='{"tumor_purity": NaN, "tumor_ploidy": 3.7, "normal_tumor_fraction": 0.2}')
case("bad_string_estimate", js='{"tumor_purity": "0.62", "tumor_ploidy": 3.7, "normal_tumor_fraction": 0.2}')
case("bad_bool_estimate", js='{"tumor_purity": true, "tumor_ploidy": 3.7, "normal_tumor_fraction": 0.2}')
case("bad_extra_key", js=json.dumps({**est, "ploidy_note": 1}))
case("bad_json_array", js="[0.62, 3.71, 0.2]")
case("bad_purity_diploid_assumption", js=json.dumps({**est, "tumor_purity": 0.54}))
case("bad_ploidy_diploid", js=json.dumps({**est, "tumor_ploidy": 2.0}))
case("bad_ploidy_doubled_solution", js=json.dumps({"tumor_purity": 0.45, "tumor_ploidy": 7.4, "normal_tumor_fraction": 0.11}))
case("bad_tin_zero", js=json.dumps({**est, "normal_tumor_fraction": 0.0}))
case("bad_tin_naive_ratio", js=json.dumps({**est, "normal_tumor_fraction": 0.255}))
case("bad_truth_plus_recurrent", tsv="chrom\tpos\tref\talt\n" + "".join("\t".join(r) + "\n" for r in truth + rec))
case("bad_truth_plus_germline", tsv="chrom\tpos\tref\talt\n" + "".join("\t".join(r) + "\n" for r in truth + germ[:200]))
case("bad_clonal_only", tsv="chrom\tpos\tref\talt\n" + "".join("\t".join(r[:4]) + "\n" for r in [l.split("\t") for l in (priv / "truth_somatic.tsv").read_text().splitlines()[1:]] if r[4] == "clonal"))
case("bad_wrong_alt_allele", tsv="chrom\tpos\tref\talt\n" + "".join(f"{c}\t{p}\t{r}\t{'ACGT'[('ACGT'.index(a) + 1) % 4] if 'ACGT'[('ACGT'.index(a) + 1) % 4] != r else 'ACGT'[('ACGT'.index(a) + 2) % 4]}\n" for c, p, r, a in truth))
case("bad_fifo", special="fifo")
case("bad_symlink_to_truth", special="symlink_truth")
case("bad_symlink_results_dir", special="symlink_dir")
case("bad_oversize", tsv=snv + "#" * (1024 * 1024))
results = {}
for name, (tsv, js, special) in cases.items():
    d = out / name
    d.mkdir(parents=True)
    if special == "empty":
        pass
    elif special == "symlink_truth":
        os.symlink("/tests/data/truth_somatic.tsv", d / "somatic_snvs.tsv")
        os.symlink("/tests/data/truth_params.json", d / "sample_estimates.json")
    elif special == "symlink_dir":
        (d / "real").mkdir()
        (d / "real" / "somatic_snvs.tsv").write_text(tsv)
        (d / "real" / "sample_estimates.json").write_text(js)
    elif special == "fifo":
        os.mkfifo(d / "somatic_snvs.tsv")
        (d / "sample_estimates.json").write_text(js)
    else:
        (d / "somatic_snvs.tsv").write_bytes(tsv.encode("utf-8", "surrogateescape") if "\xff" not in tsv else tsv[:-1].encode() + b"\xff")
        (d / "sample_estimates.json").write_text(js)
    if special == "symlink_dir":
        cmd = ["docker", "run", "--rm", "--network", "none", "-v", f"{d.resolve()}:/mnt/sub:ro", "ffpe-tests", "bash", "-c",
               "rmdir /root/results && ln -s /mnt/sub/real /root/results && bash /tests/test.sh >/dev/null 2>&1; cat /logs/verifier/reward.txt"]
    else:
        cmd = ["docker", "run", "--rm", "--network", "none", "-v", f"{d.resolve()}:/root/results:ro", "ffpe-tests",
               "bash", "-c", "bash /tests/test.sh >/dev/null 2>&1; cat /logs/verifier/reward.txt"]
    r = subprocess.run(cmd,
                       capture_output=True, text=True, timeout=300)
    reward = r.stdout.strip()
    ok = (reward == "1") == name.startswith("valid_")
    results[name] = {"reward": reward, "as_expected": ok}
    print(f"{'OK ' if ok else 'BAD'} {name:34s} reward={reward}")
summary = {"cases": len(results), "as_expected": sum(v["as_expected"] for v in results.values()), "results": results}
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(f"{summary['as_expected']}/{summary['cases']} as expected")
