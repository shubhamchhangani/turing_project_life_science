"""Summarize completed local trials without counting API/setup errors as wins.

This author-only tool never creates a ZIP. A solved candidate is held even if
another run fails. Local Codex runs cannot certify the platform's different
Terminus harness or replace a missing Opus evaluation.
"""
import argparse
import json
from datetime import datetime
from pathlib import Path

ROSTER = {("openai/gpt-6-sol", "max"), ("anthropic/claude-opus-5-5", "low")}


def classify(path):
    data = json.loads(path.read_text())
    config = data.get("config", {}).get("agent", {})
    model = config.get("model_name")
    effort = config.get("kwargs", {}).get("reasoning_effort")
    reward = (data.get("verifier_result") or {}).get("rewards", {}).get("reward")
    error = data.get("exception_info")
    commands = 0
    runtime_errors = []
    for log in path.parent.glob("agent/codex.txt"):
        for line in log.read_text().splitlines():
            try:
                event = json.loads(line)
            except ValueError:
                continue
            if event.get("type") == "item.completed" and event.get("item", {}).get("type") == "command_execution":
                commands += 1
            if event.get('type') in {'error', 'turn.failed'}:
                runtime_errors.append(str(event.get('error', event.get('message', ''))))
    trajectory = path.parent / 'agent/trajectory.json'
    if commands == 0 and trajectory.is_file():
        parsed = json.loads(trajectory.read_text())
        commands = sum(len(step.get('tool_calls') or []) for step in parsed.get('steps', []))
    for log in path.parent.glob('agent/claude-code.txt'):
        claude_calls = 0
        for line in log.read_text().splitlines():
            try:
                event = json.loads(line)
            except ValueError:
                continue
            if event.get('type') == 'assistant':
                claude_calls += sum(item.get('type') == 'tool_use' for item in event.get('message', {}).get('content', []))
            if event.get('type') == 'result' and event.get('is_error'):
                runtime_errors.append(str(event.get('result', event.get('errors', 'Claude runtime error'))))
        commands = max(commands, claude_calls)
    quota_interrupted = any(any(token in message.lower() for token in
                            ['usage limit', 'hit your limit', 'rate limit', 'quota', 'insufficient credit'])
                            for message in runtime_errors)
    execution = data.get('agent_execution') or {}
    try:
        elapsed = (datetime.fromisoformat(execution['finished_at']) - datetime.fromisoformat(execution['started_at'])).total_seconds()
    except (KeyError, TypeError, ValueError):
        elapsed = None
    full_budget_timeout = ((error or {}).get('exception_type') == 'AgentTimeoutError'
                           and elapsed is not None and elapsed >= 17999)
    finished = bool(data.get("finished_at"))
    if finished and reward == 1:
        status = "solved_do_not_package"
    elif not finished:
        status = "incomplete"
    elif quota_interrupted or runtime_errors or (error and not full_budget_timeout) or reward not in (0, 1) or not commands:
        status = "inconclusive_infrastructure_or_no_attempt"
    else:
        status = "substantive_attempt_failed"
    return {"result": str(path.resolve()), "model": model, "effort": effort,
            "status": status, "reward": reward, "completed_commands": commands,
            "error_type": (error or {}).get("exception_type"),
            "quota_interrupted": quota_interrupted,
            "agent_seconds": elapsed,
            "full_budget_timeout": full_budget_timeout,
            "counts_as_model_failure": status == "substantive_attempt_failed",
            "task_checksum": data.get("task_checksum"),
            "agent": (data.get("agent_info") or {}).get("name")}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("jobs", nargs="+", type=Path)
    p.add_argument("--report", required=True, type=Path)
    args = p.parse_args()
    results = []
    for job in args.jobs:
        for result in sorted(job.glob("*/result.json")):
            results.append(classify(result))
    failed = {(r["model"], r["effort"]) for r in results if r["counts_as_model_failure"]}
    solved = any(r["status"] == "solved_do_not_package" for r in results)
    # Deliberately do not authorize packaging from mixed checksums, incomplete
    # trials, or this script alone: reference/regression and platform evidence
    # must apply to the exact same package and required harness.
    report = {"decision": "DO_NOT_PACKAGE_SOLVED" if solved else "HOLD_FOR_COMPLETE_EVIDENCE",
              "zip_allowed": False, "trials": results,
              "required_roster_without_local_failure": [dict(model=m, effort=e) for m,e in sorted(ROSTER - failed)],
              "limitation": "Single local trials do not guarantee failure on future runs or the Turing harness."}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(report["decision"])
    for result in results:
        print(result["model"], result["effort"], result["status"])
