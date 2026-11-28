#!/usr/bin/env python3
"""Résumés de la CI et de tools/check.sh (Python seul, sans dépendance).

    summary.py job --name NOM --build-log F --junit F [--build-dir D] [--json F]
        Résumé d'un job : bibliothèques testées, références restaurées ou calculées, ccache,
        tableau des tests (rapport JUnit de CTest) ; annotations GitHub des erreurs de
        compilation, des assert et des empreintes différentes. Écrit dans
        $GITHUB_STEP_SUMMARY (sinon sur la sortie), et en JSON pour la Synthèse.

    summary.py synthesis --summaries DOSSIER
        Synthèse d'un run : tableau bloc × compilateur (API de GitHub), échecs en clair ;
        code de sortie non nul si un job requis a échoué ou a été annulé (variable RESULTS :
        toJSON(needs)).
"""

import argparse
import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

WORKSPACE = Path(os.environ.get("GITHUB_WORKSPACE", Path(__file__).resolve().parents[2]))
IN_CI = "GITHUB_ACTIONS" in os.environ

DIAGNOSTIC = re.compile(r"^(?P<file>[^\s:][^:]*):(?P<line>\d+):(?:(?P<col>\d+):)? (?P<kind>fatal error|error|warning): (?P<msg>.*)$")
ASSERT_APPLE = re.compile(r"Assertion failed: \((?P<expr>.*)\), function .*, file (?P<file>.+), line (?P<line>\d+)\.")
ASSERT_GLIBC = re.compile(r"^[^:]+: (?P<file>[^:]+):(?P<line>\d+): .*Assertion `(?P<expr>.*)' failed\.")
FINGERPRINT = re.compile(r"^(?P<label>\S.*?)\s+[0-9a-f]{16} DIFF[ÉE]RENTE")
REFERENCES = re.compile(r"Références (?P<name>\S+) : (?P<what>restaurées|calculées en (?P<seconds>\d+) s)")


def relative(path):
    try:
        return str(Path(path).resolve().relative_to(WORKSPACE.resolve()))
    except (ValueError, OSError):
        return path


def annotate(kind, file, line, message):
    """Annotation GitHub (rien hors CI)."""
    if not IN_CI:
        return
    level = "warning" if kind == "warning" else "error"
    message = message.replace("%", "%25").replace("\r", "").replace("\n", "%0A")
    print(f"::{level} file={relative(file)},line={line}::{message}")


def output(markdown):
    target = os.environ.get("GITHUB_STEP_SUMMARY")
    if target:
        with open(target, "a", encoding="utf-8") as f:
            f.write(markdown + "\n")
    else:
        print(markdown)


def duration(seconds):
    seconds = float(seconds)
    if seconds < 10:
        return f"{seconds:.1f} s".replace(".", ",")
    if seconds < 60:
        return f"{seconds:.0f} s"
    return f"{int(seconds // 60)} min {int(seconds % 60):02d}"


def read_tests(junit):
    tests = []
    if not junit or not Path(junit).exists():
        return tests
    for case in ET.parse(junit).getroot().iter("testcase"):
        out = case.findtext("system-out") or ""
        status = "ok"
        if case.find("failure") is not None or case.get("status") == "fail":
            status = "fail"
        elif case.find("skipped") is not None or case.get("status") in ("notrun", "disabled"):
            status = "skip"
        tests.append({"name": case.get("name"), "time": float(case.get("time") or 0), "status": status, "output": out})
    return tests


def labels_of(build_dir):
    if not build_dir:
        return {}
    try:
        data = json.loads(subprocess.run(["ctest", "--test-dir", build_dir, "--show-only=json-v1"],
                                         capture_output=True, text=True, check=True).stdout)
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError):
        return {}
    result = {}
    for test in data.get("tests", []):
        for prop in test.get("properties", []):
            if prop.get("name") == "LABELS":
                result[test["name"]] = " ".join(prop.get("value", []))
    return result


def failure_detail(test):
    """Première ligne utile de la sortie d'un test en échec, avec son annotation."""
    source = WORKSPACE / "tests" / (test["name"].replace(".", "/") + ".cpp")
    for line in test["output"].splitlines():
        for pattern in (ASSERT_APPLE, ASSERT_GLIBC):
            m = pattern.search(line)
            if m:
                annotate("error", m["file"], m["line"], f"{test['name']} : assert ({m['expr']}) échoue")
                return f"assert ({m['expr']}), {relative(m['file'])}:{m['line']}"
        m = FINGERPRINT.search(line)
        if m:
            where = 1
            if source.exists():
                for number, text in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
                    if f'"{m["label"]}"' in text:
                        where = number
                        break
            annotate("error", source, where, f"{test['name']} : empreinte « {m['label']} » différente")
            return line.strip()
    lines = [l for l in test["output"].splitlines() if l.strip()]
    return lines[-1].strip() if lines else "échec"


def ccache_rate():
    try:
        stats = dict(line.split("\t")[:2] for line in subprocess.run(
            ["ccache", "--print-stats"], capture_output=True, text=True, check=True).stdout.splitlines() if "\t" in line)
    except (OSError, subprocess.CalledProcessError):
        return None
    hits = int(stats.get("direct_cache_hit", 0)) + int(stats.get("preprocessed_cache_hit", 0))
    total = hits + int(stats.get("cache_miss", 0))
    return f"{100 * hits / total:.0f} % ({hits} sur {total})" if total else None


def job(args):
    log = Path(args.build_log).read_text(encoding="utf-8", errors="replace") if args.build_log and Path(args.build_log).exists() else ""
    diagnostics = []
    for line in log.splitlines():
        m = DIAGNOSTIC.match(line)
        if m and (m["file"], m["line"], m["msg"]) not in {(d[0], d[1], d[3]) for d in diagnostics}:
            diagnostics.append((m["file"], m["line"], m["kind"], m["msg"]))
    for file, line, kind, msg in diagnostics[:20]:
        annotate(kind, file, line, msg)
    restored = [m["name"] for m in REFERENCES.finditer(log) if m["what"] == "restaurées"]
    computed = [f"{m['name']} ({duration(m['seconds'])})" for m in REFERENCES.finditer(log) if m["seconds"]]
    tests = read_tests(args.junit)
    labels = labels_of(args.build_dir)
    failed = [t for t in tests if t["status"] == "fail"]
    errors = [d for d in diagnostics if d[2] != "warning"]
    ok = bool(tests) and not failed and not errors if not args.allow_empty else not failed and not errors
    passed = sum(t["status"] == "ok" for t in tests)
    head = f"### {args.name} — {'✓' if ok else '✗'} {passed} tests sur {len(tests)}"
    lines = [head]
    tested = Path(args.build_dir or "", "tested_libraries.txt")
    if args.build_dir and tested.exists():
        lines.append(f"Bibliothèques testées : {tested.read_text().strip().replace(';', ', ')}")
    if restored or computed:
        lines.append(f"Références : {len(restored)} familles restaurées"
                     + (f", calculées : {', '.join(computed)}" if computed else ""))
    rate = ccache_rate()
    if rate:
        lines.append(f"ccache : {rate}")
    if diagnostics:
        warnings = len(diagnostics) - len(errors)
        lines.append(f"Diagnostics du compilateur : {len(errors)} erreur{'s' * (len(errors) > 1)}, "
                     f"{warnings} avertissement{'s' * (warnings > 1)}")
    details = [{"test": t["name"], "detail": failure_detail(t)} for t in failed]
    details += [{"test": "compilation", "detail": f"{relative(d[0])}:{d[1]}: {d[3]}"} for d in errors[:5]]
    if details:
        lines.append("")
        lines += [f"- ✗ **{d['test']}** : {d['detail']}" for d in details]
    if tests:
        lines += ["", "| Test | Étiquettes | Durée | |", "|---|---|---|---|"]
        mark = {"ok": "✓", "fail": "✗", "skip": "—"}
        for t in sorted(tests, key=lambda t: (t["status"] != "fail", t["name"])):
            lines.append(f"| {t['name']} | {labels.get(t['name'], '')} | {duration(t['time'])} | {mark[t['status']]} |")
    if args.brief:
        brief = [f"{args.name} : {'✓' if ok else '✗'} {passed} tests sur {len(tests)}"]
        brief += [f"  ✗ {d['test']} : {d['detail']}" for d in details]
        print("\n".join(brief))
    else:
        output("\n".join(lines))
    if args.json:
        Path(args.json).write_text(json.dumps({"name": args.name, "ok": ok, "passed": passed, "total": len(tests),
                                               "failures": details}, ensure_ascii=False), encoding="utf-8")
    return 0 if ok else 1


def github(path):
    result = subprocess.run(["gh", "api", "--paginate", path], capture_output=True, text=True)
    return json.loads(result.stdout) if result.returncode == 0 and result.stdout.strip() else None


def synthesis(args):
    results = json.loads(os.environ.get("RESULTS", "{}"))
    repo, run, attempt = os.environ.get("GITHUB_REPOSITORY"), os.environ.get("GITHUB_RUN_ID"), os.environ.get("GITHUB_RUN_ATTEMPT", "1")
    bad = {name: r["result"] for name, r in results.items() if r.get("result") in ("failure", "cancelled")}
    mode = results.get("changes", {}).get("outputs", {}).get("mode", "")
    lines = [f"## CI — {'✓' if not bad else '✗'}"]
    if mode == "light":
        sha = os.environ.get("GITHUB_SHA", "")
        pulls = github(f"repos/{repo}/commits/{sha}/pulls") or []
        if pulls:
            lines.append(f"Fusion de la PR #{pulls[0]['number']} ({pulls[0]['title']}) : contenu testé par sa CI.")
        else:
            lines.append("Run léger après fusion : contenu testé par la CI de sa PR.")
    jobs = (github(f"repos/{repo}/actions/runs/{run}/attempts/{attempt}/jobs?per_page=100") or {}).get("jobs", [])
    from datetime import datetime
    cells, compilers = {}, []
    mark = {"success": "✓", "failure": "✗", "cancelled": "⊘", "skipped": "—"}
    for j in jobs:
        if " · " not in j["name"]:
            continue
        block, compiler = j["name"].split(" · ", 1)
        if compiler not in compilers:
            compilers.append(compiler)
        took = ""
        if j.get("started_at") and j.get("completed_at"):
            delta = datetime.fromisoformat(j["completed_at"].replace("Z", "+00:00")) - datetime.fromisoformat(j["started_at"].replace("Z", "+00:00"))
            took = f" {int(delta.total_seconds() // 60)}:{int(delta.total_seconds() % 60):02d}"
        cells[(block, compiler)] = mark.get(j.get("conclusion") or "", "…") + took
    if cells:
        blocks = [b for b in ("Release", "Sanitizers", "FMA") if any(k[0] == b for k in cells)]
        lines += ["", "| | " + " | ".join(compilers) + " |", "|---|" + "---|" * len(compilers)]
        for b in blocks:
            lines.append(f"| {b} | " + " | ".join(cells.get((b, c), "") for c in compilers) + " |")
    for name in ("package", "references"):
        if name in results and results[name]["result"] != "skipped":
            label = {"package": "Paquet CMake", "references": "Références"}[name]
            lines.append(f"\n{label} : {mark.get(results[name]['result'], results[name]['result'])}")
    if all(r.get("result") == "skipped" for n, r in results.items() if n != "changes") and mode == "tests":
        lines.append("\nRien à tester : documentation seule.")
    for summary in sorted(Path(args.summaries).glob("*/summary.json")) if Path(args.summaries).exists() else []:
        data = json.loads(summary.read_text(encoding="utf-8"))
        for failure in data.get("failures", []):
            lines.append(f"- ✗ {data['name']} : **{failure['test']}** — {failure['detail']}")
    if bad:
        lines.append("\nÉchecs : " + ", ".join(f"{n} ({r})" for n, r in bad.items()))
    output("\n".join(lines))
    return 1 if bad else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    j = sub.add_parser("job")
    j.add_argument("--name", required=True)
    j.add_argument("--build-log")
    j.add_argument("--junit")
    j.add_argument("--build-dir")
    j.add_argument("--json")
    j.add_argument("--allow-empty", action="store_true", help="aucun test n'est pas un échec")
    j.add_argument("--brief", action="store_true", help="résumé court sur la sortie (tools/check.sh)")
    s = sub.add_parser("synthesis")
    s.add_argument("--summaries", default="summaries")
    args = parser.parse_args()
    return job(args) if args.command == "job" else synthesis(args)


if __name__ == "__main__":
    sys.exit(main())
