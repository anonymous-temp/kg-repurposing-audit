"""Frozen public registry retrieval and descriptive reason-text audits.

Registry status and lexical cues are metadata. They are not efficacy labels.
Raw responses remain in the user's data directory and are never required in Git.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request

BASE = "https://clinicaltrials.gov/api/v2/studies"
PATTERNS = {
    "efficacy": r"\bfutil\w*\b|lack of (?:efficacy|effectiveness)|insufficient efficacy|inefficacy|lack of benefit|no (?:clinical )?(?:benefit|efficacy)",
    "safety": r"\bsafety\b|\btoxic\w*\b|\badverse event\w*\b|\btolerab\w*\b",
    "operational": r"recruit\w*|enrol\w*|accrual|funding|financial|business|commercial|sponsor decision|strategic|budget|administrative|staffing|investigator|pandemic|covid",
}


def reason_cues(text):
    return sorted(k for k, pattern in PATTERNS.items() if re.search(pattern, text or "", re.I))


def extract(study):
    p = study["protocolSection"]
    status = p["statusModule"]
    ident = p["identificationModule"]
    interventions = p.get("armsInterventionsModule", {}).get("interventions", [])
    return {
        "nct_id": ident["nctId"],
        "title": ident.get("briefTitle", ""),
        "registry_status": status.get("overallStatus"),
        "why_stopped": status.get("whyStopped", ""),
        "last_update": status.get("lastUpdatePostDateStruct", {}).get("date"),
        "has_results": bool(study.get("hasResults", False)),
        "drug_names": [x["name"] for x in interventions if x.get("type") == "DRUG"],
        "conditions": p.get("conditionsModule", {}).get("conditions", []),
        "source": f"https://clinicaltrials.gov/study/{ident['nctId']}",
        "reason_cues": reason_cues(status.get("whyStopped", "")),
        "writeback_action": "defer",
    }


def fetch_sample(outdir, cutoff, per_status=100):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    records = []
    manifest = []
    for status in ["TERMINATED", "WITHDRAWN", "SUSPENDED"]:
        params = {
            "query.term": f"AREA[StudyType]INTERVENTIONAL AND AREA[InterventionType]DRUG AND AREA[StudyFirstPostDate]RANGE[2015-01-01,2025-12-31] AND AREA[LastUpdatePostDate]RANGE[MIN,{cutoff}]",
            "filter.overallStatus": status,
            "pageSize": per_status,
            "sort": "StudyFirstPostDate:asc",
            "countTotal": "true",
            "format": "json",
        }
        url = BASE + "?" + urllib.parse.urlencode(params)
        path = outdir / f"{status.lower()}.json"
        request_path = outdir / f"{status.lower()}_request.json"
        if path.exists():
            if not request_path.exists() or json.loads(request_path.read_text()).get("url") != url:
                raise ValueError("Cached snapshot has a different or unrecorded request; use a new output directory")
            raw = path.read_bytes()
        else:
            request = urllib.request.Request(url, headers={"User-Agent": "kg-repurposing-audit/0.2 scholarly research"})
            with urllib.request.urlopen(request, timeout=60) as response:
                raw = response.read()
            path.write_bytes(raw)
            request_path.write_text(json.dumps({"url": url}, indent=2))
            time.sleep(0.4)
        data = json.loads(raw)
        studies = data.get("studies", [])
        if not studies:
            raise ValueError(f"No studies returned for {status}")
        chunk = [extract(s) for s in studies]
        if any(not x["drug_names"] or x["registry_status"] != status for x in chunk):
            raise ValueError("Registry response does not satisfy requested intervention/status filter")
        if any(x["last_update"] is None or x["last_update"] > cutoff for x in chunk):
            raise ValueError("Registry update exceeds specified cutoff or is missing")
        records.extend(chunk)
        manifest.append(
            {
                "status": status,
                "url": url,
                "total_matching": data.get("totalCount"),
                "sample_size": len(chunk),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
    if len({r["nct_id"] for r in records}) != len(records):
        raise ValueError("Duplicate registry identifier")
    summary = {}
    for status in ["TERMINATED", "WITHDRAWN", "SUSPENDED"]:
        rows = [r for r in records if r["registry_status"] == status]
        summary[status] = {
            "n": len(rows),
            "reason_present": sum(bool(x["why_stopped"].strip()) for x in rows),
            "results_posted": sum(x["has_results"] for x in rows),
            "efficacy_cue": sum("efficacy" in x["reason_cues"] for x in rows),
            "safety_cue": sum("safety" in x["reason_cues"] for x in rows),
            "operational_cue": sum("operational" in x["reason_cues"] for x in rows),
            "no_recognised_cue": sum(not x["reason_cues"] for x in rows),
            "multiple_cues": sum(len(x["reason_cues"]) > 1 for x in rows),
            "global_negative_edges_added": 0,
        }
    result = {
        "design": "First N records returned in ascending StudyFirstPostDate order within each discontinued-status stratum; interventional drug studies first posted during 2015-2025; last update at or before cutoff. Convenience sample, not a population-prevalence estimate.",
        "cutoff": cutoff,
        "per_status": per_status,
        "summary": summary,
        "interpretation": "Lexical descriptors, not expert-adjudicated failure reasons or efficacy labels. A posted result does not establish a negative outcome.",
        "sources": manifest,
    }
    (outdir / "records.json").write_text(json.dumps(records, indent=2))
    (outdir / "summary.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True)
    p.add_argument("--cutoff", required=True)
    p.add_argument("--per-status", type=int, default=100)
    args = p.parse_args()
    if not 1 <= args.per_status <= 1000:
        p.error("per-status must be between 1 and 1000")
    print(json.dumps(fetch_sample(args.output, args.cutoff, args.per_status)["summary"], indent=2))


if __name__ == "__main__":
    main()
