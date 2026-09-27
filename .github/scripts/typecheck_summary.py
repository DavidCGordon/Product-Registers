"""Summarize pyright errors by file, for the advisory CI job and for local use.

The codebase carries a standing backlog of type errors, so the useful output is
not the list of findings but where they are concentrated: a per-file count tells
you which module to open if you want to reduce the number, and makes a sudden
jump in one file obvious.

Usage::

    python .github/scripts/typecheck_summary.py src test
    python .github/scripts/typecheck_summary.py src --limit 25

Always exits 0 -- this is advisory and must never fail a build.
"""
import argparse
import collections
import json
import os
import subprocess
import sys


def pyright_diagnostics(path: str) -> list[dict]:
    """Run pyright over `path` and return its error-level diagnostics.

    :param path: File or directory to check.
    :type path: str
    :return: The diagnostics pyright reported, empty if it could not be run.
    :rtype: list[dict]
    """
    try:
        completed = subprocess.run(
            [sys.executable, "-m", "pyright", "--level", "error", "--outputjson", path],
            capture_output=True, text=True,
        )
    except OSError as exc:                      # pyright missing entirely
        print(f"could not run pyright on {path}: {exc}", file=sys.stderr)
        return []

    # pyright exits nonzero whenever it reports anything, so the exit code says
    # nothing about whether the run itself succeeded -- parse and see.
    try:
        return json.loads(completed.stdout).get("generalDiagnostics", [])
    except json.JSONDecodeError:
        print(f"could not parse pyright output for {path}:\n{completed.stderr}",
              file=sys.stderr)
        return []


def counts_by_file(paths: list[str]) -> collections.Counter:
    """Error counts keyed by repository-relative file path."""
    root = os.getcwd()
    counts = collections.Counter()
    for path in paths:
        for diagnostic in pyright_diagnostics(path):
            name = diagnostic.get("file", "<unknown>")
            try:
                name = os.path.relpath(name, root)
            except ValueError:                  # different drive on Windows
                pass
            counts[name.replace(os.sep, "/")] += 1
    return counts


def render(counts: collections.Counter, limit: int) -> str:
    """Render the counts as a GitHub-flavoured markdown table."""
    total = sum(counts.values())
    if not total:
        return "## Type check (pyright, basic)\n\nNo errors.\n"

    lines = [
        "## Type check (pyright, basic)",
        "",
        f"**{total} errors across {len(counts)} files.** "
        f"Showing the {min(limit, len(counts))} with the most.",
        "",
        "| file | errors |",
        "| --- | ---: |",
    ]
    for name, n in counts.most_common(limit):
        lines.append(f"| `{name}` | {n} |")
    if len(counts) > limit:
        remaining = total - sum(n for _, n in counts.most_common(limit))
        lines.append(f"| _{len(counts) - limit} other files_ | {remaining} |")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="files or directories to check")
    parser.add_argument("--limit", type=int, default=15,
                        help="how many files to list (default: 15)")
    args = parser.parse_args()

    report = render(counts_by_file(args.paths), args.limit)
    print(report)

    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as handle:
            handle.write(report + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
