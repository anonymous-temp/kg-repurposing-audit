"""Command line tools with a database-free synthetic demonstration."""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .comparison import paired_interval
from .evidence import assess_strategy


def demo(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    frame = pd.DataFrame(
        {
            "compound": np.repeat(["toy:A", "toy:B", "toy:C", "toy:D"], 8),
            "disease": np.tile(["toy:X", "toy:Y", "toy:Z", "toy:W", "toy:U", "toy:V", "toy:R", "toy:S"], 4),
            "label": np.tile([1, 0, 0, 0, 0, 0, 0, 0], 4),
        }
    )
    frame["model"] = rng.normal(size=len(frame)) + frame.label * 0.8
    frame["reference"] = rng.normal(size=len(frame))
    result = paired_interval(frame.label, frame.model, frame.reference, frame.compound, frame.disease, n_bootstrap=300)
    record = {
        "id": "toy:strategy",
        "drug": "toy:drug",
        "indication": "toy:disease",
        "population": None,
        "comparator": None,
        "regimen": None,
        "endpoint": None,
        "required_exposure": None,
        "evidence": [],
        "constraints": [],
        "failures": [],
    }
    result["synthetic_example"] = True
    result["documentation_gate"] = assess_strategy(record, "2026-01-01")
    frame.to_csv(output / "synthetic_scores.csv", index=False)
    (output / "demo_report.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(prog="kg-audit", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    d = commands.add_parser("demo")
    d.add_argument("--output", default="outputs/demo")
    c = commands.add_parser("compare")
    c.add_argument("--scores", required=True)
    c.add_argument("--model", nargs="+", required=True)
    c.add_argument("--reference", nargs="+", required=True)
    c.add_argument("--clusters", choices=["disease", "compound", "two_way"], default="disease")
    c.add_argument("--bootstrap", type=int, default=5000)
    c.add_argument("--seed", type=int, default=20261003)
    c.add_argument("--output", required=True)
    e = commands.add_parser("evidence")
    e.add_argument("--record", required=True)
    e.add_argument("--as-of", required=True)
    e.add_argument("--output", required=True)
    a = parser.parse_args()
    if a.command == "demo":
        demo(a.output)
        return
    if a.command == "compare":
        path = Path(a.scores)
        sep = "\t" if ".tsv" in path.name else ","
        frame = pd.read_csv(path, sep=sep)
        needed = ["compound", "disease", "label"] + a.model + a.reference
        if any(x not in frame for x in needed):
            parser.error("Input score columns are missing; see docs/input-contracts.md")
        if frame.groupby(["compound", "disease"]).label.nunique().max() > 1:
            parser.error("A pair has conflicting evaluation labels")
        result = paired_interval(
            frame.label,
            frame[a.model].to_numpy().T,
            frame[a.reference].to_numpy().T,
            frame.compound,
            frame.disease,
            a.clusters,
            a.bootstrap,
            a.seed,
        )
        result["input_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    else:
        result = assess_strategy(json.loads(Path(a.record).read_text()), a.as_of)
    target = Path(a.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
