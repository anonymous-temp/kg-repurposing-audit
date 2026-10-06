"""Fetch arm-level design of stopped trials from the ClinicalTrials.gov API v2 (batches of 500 ids).
Idempotent: each batch is written to its own file and skipped if present."""
import json, os, sys, time, urllib.parse, urllib.request
ID_FILE = "work/data/ctgov/stopped_nct_ids.txt"
if not os.path.exists(ID_FILE):     # NCT identifiers of all stopped trials in Open Targets 26.09
    import pyarrow.parquet as pq
    t = pq.read_table("work/data/ot/clinical_report.parquet", columns=["id", "origin", "trialOverallStatus"]).to_pandas()
    s = t[(t.origin == "CLINICAL_TRIAL") & t.trialOverallStatus.isin(["TERMINATED", "WITHDRAWN", "SUSPENDED"])]
    os.makedirs(os.path.dirname(ID_FILE), exist_ok=True)
    open(ID_FILE, "w").write("\n".join(sorted(set(s.id.str.upper()))) + "\n")
IDS = [l.strip() for l in open(ID_FILE) if l.strip()]
OUT = "work/data/ctgov/batches"; os.makedirs(OUT, exist_ok=True)
FIELDS = "NCTId,ArmGroup,Intervention,Condition,Phase,OverallStatus,WhyStopped,StartDate,DesignAllocation,HasResults,LeadSponsorName"
B = 500
for k in range(0, len(IDS), B):
    f = f"{OUT}/b{k // B:03d}.json"
    if os.path.exists(f):
        continue
    ids = IDS[k:k + B]
    q = urllib.parse.urlencode({"filter.ids": ",".join(ids), "fields": FIELDS, "pageSize": 1000})
    for attempt in range(5):
        try:
            with urllib.request.urlopen("https://clinicaltrials.gov/api/v2/studies?" + q, timeout=120) as r:
                d = json.load(r)
            break
        except Exception as e:
            print("retry", k, e, flush=True); time.sleep(5 * (attempt + 1))
    else:
        sys.exit(f"failed batch {k}")
    assert not d.get("nextPageToken"), "unexpected pagination"
    json.dump(d["studies"], open(f + ".tmp", "w")); os.replace(f + ".tmp", f)
    print(k, len(ids), len(d["studies"]), flush=True); time.sleep(1)
print("done")
