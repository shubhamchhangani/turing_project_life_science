"""Check known Turing package failures; does not predict portal review outcomes.

Usage: python3 tools/preflight_task.py TASK --report REPORT.json [--zip ZIP]
This reads an existing ZIP if supplied; it never creates one. Behavioral
verifier checks, scientific review, and completed model trials remain required.
"""
import argparse
import ast
import hashlib
import json
import re
import subprocess
import sys
import tomllib
import zipfile
from pathlib import Path, PurePosixPath


def audit(root, archive=None):
    errors = []
    checks = []

    def check(label, condition, detail=""):
        checks.append({"check": label, "passed": bool(condition), "detail": detail})
        if not condition:
            errors.append(label + (": " + detail if detail else ""))

    required = {"environment", "tests", "solution", "task.toml", "instruction.md"}
    entries = {p.name for p in root.iterdir()}
    check("root structure", required <= entries <= required | {"README.md"},
          f"missing={sorted(required - entries)}, extra={sorted(entries - required - {'README.md'})}")
    for name in ("environment", "tests", "solution"):
        check(f"{name} is a directory", (root / name).is_dir())
    files = sorted(p for p in root.rglob("*") if p.is_file())
    manifest = {}
    unwanted = []
    canaries = []
    symlinks = [str(p.relative_to(root)) for p in root.rglob("*") if p.is_symlink()]
    for p in files:
        relative = p.relative_to(root).as_posix()
        data = p.read_bytes()
        manifest[relative] = hashlib.sha256(data).hexdigest()
        if any(s.startswith(".") or s == "__pycache__" for s in p.relative_to(root).parts) or p.suffix in {".pyc", ".nbi", ".nbc"}:
            unwanted.append(relative)
        if b"harbor-canary" in data.lower():
            canaries.append(relative)
    check("no hidden/cache files", not unwanted, ", ".join(unwanted))
    check("no canary markers", not canaries, ", ".join(canaries))
    check("no symlinks in package", not symlinks, ", ".join(symlinks))
    try:
        config = tomllib.loads((root / "task.toml").read_text())
    except (OSError, ValueError) as exc:
        check("valid TOML", False, str(exc))
        config = {}
    check("schema 1.4", config.get("schema_version") == "1.4")
    check("top-level artifacts", isinstance(config.get("artifacts"), list) and
          bool(config.get("artifacts")) and all(isinstance(v, str) and v.startswith("/root/results/") for v in config.get("artifacts", [])))
    fields = {
        "task": {"name", "description", "version", "keywords"},
        "metadata": {"domain", "category", "subcategory", "reward_type", "source", "public_overlap"},
        "environment": {"network_mode", "cpus", "memory_mb", "storage_mb"},
        "agent": {"timeout_sec", "network_mode"},
        "verifier": {"timeout_sec", "network_mode", "environment_mode", "environment"},
    }
    for section, keys in fields.items():
        actual = config.get(section, {})
        check(f"{section} required fields", isinstance(actual, dict) and keys <= actual.keys())
    task, metadata = config.get("task", {}), config.get("metadata", {})
    environment, agent = config.get("environment", {}), config.get("agent", {})
    verifier = config.get("verifier", {})
    check("turing task namespace", bool(re.fullmatch(r"turing/[a-z0-9]+(?:-[a-z0-9]+)*", task.get("name", ""))))
    check("no legacy metadata", not ({"expert_time_estimate_hours", "field", "subfield", "tags"} & metadata.keys()))
    # The category pair requested for these biology packages; this is not an
    # invented validation of categories for unrelated Turing domains.
    check("life sciences biology category", metadata.get("domain") == "life-sciences" and metadata.get("category") == "biology")
    check("binary reward and overlap type", metadata.get("reward_type") == "binary" and isinstance(metadata.get("public_overlap"), bool))
    check("agent execution contract", agent.get("timeout_sec") == 18000 and agent.get("network_mode") == "public")
    check("agent environment network", environment.get("network_mode") == "public")
    check("isolated offline verifier", verifier.get("network_mode") == "no-network" and verifier.get("environment_mode") == "separate" and verifier.get("environment", {}).get("network_mode") == "no-network")
    for label, section, keys in [
        ("agent resources", environment, ["cpus", "memory_mb", "storage_mb"]),
        ("verifier resources", verifier.get("environment", {}), ["cpus", "memory_mb"]),
        ("verifier timeout", verifier, ["timeout_sec"]),
    ]:
        check(label, all(isinstance(section.get(k), (int, float)) and not isinstance(section.get(k), bool) and section[k] > 0 for k in keys))
    runner = root / "tests/test.sh"
    if runner.is_file():
        lines = [line.strip() for line in runner.read_text().splitlines() if line.strip() and not line.startswith("#")]
        check("test runner errexit", bool(lines) and lines[0] == "set -euo pipefail")
        result = subprocess.run(["bash", "-n", str(runner)], capture_output=True, text=True)
        check("test runner shell syntax", result.returncode == 0, result.stderr.strip())
    else:
        check("test runner exists", False)
    for p in files:
        if p.suffix == ".py":
            try:
                tree = ast.parse(p.read_text())
                undocumented = [n.name for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name.startswith("test_") and not ast.get_docstring(n)]
                check(f"Python syntax/docstrings: {p.relative_to(root)}", not undocumented, ", ".join(undocumented))
            except (SyntaxError, UnicodeError) as exc:
                check(f"Python syntax: {p.relative_to(root)}", False, str(exc))
    instruction = (root / "instruction.md").read_text() if (root / "instruction.md").is_file() else ""
    closing = "You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task."
    check("execution-contract closing", instruction.rstrip().endswith(closing))
    check("all artifacts stated", all(a in instruction for a in config.get("artifacts", [])))
    check("no scoring language in instructions", not re.search(r"\b(?:reward|grader|pass@\d|binary score)\b", instruction, re.I))
    if archive:
        with zipfile.ZipFile(archive) as z:
            names = [n for n in z.namelist() if not n.endswith("/")]
            unsafe = [n for n in names if PurePosixPath(n).is_absolute() or ".." in PurePosixPath(n).parts]
            check("ZIP safe paths", not unsafe, ", ".join(unsafe))
            check("ZIP unique members", len(names) == len(set(names)))
            check("ZIP direct root and exact file set", set(names) == set(manifest))
            check("ZIP matches source bytes", set(names) == set(manifest) and all(hashlib.sha256(z.read(n)).hexdigest() == manifest[n] for n in names))
    fingerprint = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    return {"task": str(root.resolve()), "static_checks_passed": not errors,
            "portal_acceptance_guaranteed": False, "package_fingerprint": fingerprint,
            "files": manifest, "checks": checks, "errors": errors,
            "required_additional_evidence": ["offline reference solution and verifier", "valid and invalid artifact regressions", "scientific coverage and independent truth review", "completed substantive trials on the required model roster"],
            "zip_created": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--zip", type=Path)
    args = parser.parse_args()
    if args.report.resolve().is_relative_to(args.task.resolve()):
        parser.error("Keep author reports outside the task package.")
    report = audit(args.task, args.zip)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"{sum(c['passed'] for c in report['checks'])}/{len(report['checks'])} static checks passed; report: {args.report}")
    for error in report["errors"]:
        print("FAIL:", error)
    sys.exit(0 if report["static_checks_passed"] else 1)
