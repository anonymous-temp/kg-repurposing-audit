"""Export benchmark predictions as explicitly incomplete evidence-review records."""

import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from .evidence import assess_strategy, DOMAINS


def export_records(scores_path, score_column, output, graph, task, sampler, seed, as_of):
    path = Path(scores_path)
    frame = pd.read_csv(path, sep="\t" if ".tsv" in path.name else ",")
    keys = ["compound", "disease"]
    if frame.duplicated(keys).any():
        raise ValueError("Export requires unique candidate pairs; aggregate repeated sampling draws first")
    if score_column not in frame:
        raise ValueError("Requested score column is absent")
    frame["rank"] = frame.groupby("disease")[score_column].rank(ascending=False, method="min").astype(int)
    frame["candidate_count"] = frame.groupby("disease").compound.transform("count")
    source_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    counts = {}
    n = 0
    with out.open("w") as handle:
        for row in frame.to_dict("records"):
            pair = json.dumps([graph, task, sampler, seed, row["compound"], row["disease"]], separators=(",", ":"))
            record = {
                "id": "prediction:" + hashlib.sha256(pair.encode()).hexdigest()[:24],
                "drug": row["compound"],
                "indication": row["disease"],
                "assessment_date": as_of,
                "prediction": {
                    "model": score_column,
                    "graph": graph,
                    "task": task,
                    "candidate_sampler": sampler,
                    "seed": seed,
                    "score": float(row[score_column]),
                    "rank": int(row["rank"]),
                    "rank_scope": "same-disease candidates in the supplied evaluation manifest",
                    "candidate_count": int(row["candidate_count"]),
                    "source_sha256": source_hash,
                    "score_is_probability": False,
                },
                "evidence": [],
                "constraints": [{"domain": d, "status": "unknown", "evidence_ids": []} for d in DOMAINS],
                "failures": [],
            }
            status = assess_strategy(record, as_of)["status"]
            counts[status] = counts.get(status, 0) + 1
            handle.write(json.dumps(record, separators=(",", ":")) + "\n")
            n += 1
    report = {
        "records": n,
        "documentation_status": counts,
        "source_sha256": source_hash,
        "rank_domain": "Per-disease observed evaluation candidate pool; not the full therapeutic search space",
        "biomedical_fields_inferred": 0,
    }
    out.with_suffix(".summary.json").write_text(json.dumps(report, indent=2))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ["scores", "score-column", "output", "graph", "task", "sampler", "as-of"]:
        p.add_argument("--" + k, required=True)
    p.add_argument("--seed", type=int, required=True)
    a = p.parse_args()
    print(
        json.dumps(
            export_records(a.scores, a.score_column, a.output, a.graph, a.task, a.sampler, a.seed, a.as_of), indent=2
        )
    )


if __name__ == "__main__":
    main()
