"""Command-line interface for evaluating synthetic organ-toxicity vignettes."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from organ_tox.evaluators import evaluate_patient_vignette

DEFAULT_DEMO_PATH = (
    Path(__file__).resolve().parent.parent / "data" / "synthetic_vignettes.jsonl"
)


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            stripped = raw.strip()
            if not stripped:
                continue
            try:
                records.append(json.loads(stripped))
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON on line {line_no} of {path}: {exc}") from exc
    return records


def run_batch(input_path: Path) -> dict[str, Any]:
    records = load_jsonl(input_path)
    results = [evaluate_patient_vignette(record) for record in records]
    critical_cases = sum(1 for item in results if item["max_severity_tier"] == "CRITICAL")
    moderate_cases = sum(1 for item in results if item["max_severity_tier"] == "MODERATE")
    return {
        "input_path": str(input_path),
        "vignettes_evaluated": len(results),
        "critical_cases": critical_cases,
        "moderate_cases": moderate_cases,
        "per_case_results": results,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="organ-tox",
        description="Evaluate zero-PHI synthetic vignettes against FDA/KDIGO/CTCAE/ASTCT toxicity criteria.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--demo",
        action="store_true",
        help="Run bundled synthetic toxicity vignettes (data/synthetic_vignettes.jsonl).",
    )
    group.add_argument(
        "--input",
        type=Path,
        help="Path to a synthetic JSONL file matching the organ-tox schema.",
    )
    args = parser.parse_args(argv)

    target_path = DEFAULT_DEMO_PATH if args.demo else args.input
    if not target_path.exists():
        print(f"Error: input file not found: {target_path}", file=sys.stderr)
        return 1

    summary = run_batch(target_path)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
