#!/usr/bin/env python
"""Case review (Additional file 1, Section S21): select recorded or approved indications whose cleanest failure record
(same-concept efficacy or safety stop of the investigational agent, trial started before 2015) the strictest rule would
negate, fetch their trials from the ClinicalTrials.gov API v2 and write one dossier per case.

  python 20_case_review_dossiers.py hetionet|primekg   (DATA and HET as for the other pipelines; OUT=work/cases)
Then: python 20_case_review_dossiers.py fetch ; python 20_case_review_dossiers.py dossier
The codes themselves (paper_results/selection/case_taxonomy.tsv) follow the codebook in Section S21.
"""
import glob, json, os, random, sys, time, urllib.parse, urllib.request
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
OUT = os.environ.get("OUT", "work/cases"); os.makedirs(OUT, exist_ok=True)
stage = sys.argv[1]

if stage in ("hetionet", "primekg"):
    import pandas as pd
    from kg_audit import writeback as L
    if stage == "hetionet":
        comps, dises, names, ctd, cpd = L.load_graph(); ctd = set(ctd)
    else:
        G = json.load(open(f"{L.DATA}/graph.json")); ctd = {tuple(p) for p in G["indication"]}; names = G["names"]
    ev = L.load_evidence(); ind = ctd | ev["approved"]
    st = ev["stopped"]; st = st[st.start < L.CUTOFF_WRITEBACK].copy()
    roles = pd.read_csv(f"{L.DATA}/trial_roles.tsv", sep="\t"); dm = pd.read_csv(f"{L.DATA}/disease_map.tsv", sep="\t")
    gap = dict(zip(zip(dm.ot_disease, dm.disease), dm.depth_gap))
    st["gap"] = [gap.get((o, d), -1) for o, d in zip(st.ot_disease, st.disease)]
    st = st.merge(roles, on=["report_id", "compound"], how="left")
    st["sci"] = st.stop_categories.fillna("").apply(lambda s: bool(set(s.split("|")) & L.SCIENTIFIC))
    strict = st[st.sci & (st.role == "investigational") & (st.gap == 0)]
    pairs = sorted(p for p in set(strict.pair) if p in ind)
    random.seed(20261007)
    if len(pairs) > 40:
        pairs = sorted(random.sample(pairs, 40))
    npre = ev["trials"][ev["trials"].start < L.CUTOFF_WRITEBACK].groupby("pair").report_id.nunique()
    out = [{"graph": stage, "compound": p[0], "drug": names.get(p[0], p[0]), "disease_id": p[1], "disease": names.get(p[1], p[1]),
            "label": "recorded treatment" if p in ctd else "approval evidence only", "n_trials_before_2015": int(npre.get(p, 0)),
            "trials": [{"nct": r.report_id.upper(), "phase": r.phase, "stop_categories": r.stop_categories, "why_stopped": r.why_stopped,
                        "ot_condition": r.ot_disease, "start": str(r.start_date)} for r in strict[strict.pair == p].drop_duplicates("report_id").itertuples()]}
           for p in pairs]
    json.dump(out, open(f"{OUT}/pairs_{stage}.json", "w"), indent=1); print(stage, len(out), "cases")

elif stage == "fetch":
    ids = sorted({t["nct"] for f in glob.glob(f"{OUT}/pairs_*.json") for p in json.load(open(f)) for t in p["trials"]})
    fields = ("protocolSection.identificationModule,protocolSection.statusModule,protocolSection.conditionsModule,protocolSection.designModule,"
              "protocolSection.armsInterventionsModule,protocolSection.outcomesModule.primaryOutcomes,protocolSection.descriptionModule.briefSummary,"
              "protocolSection.eligibilityModule,hasResults")
    got = {}
    for i in range(0, len(ids), 50):
        q = urllib.parse.urlencode({"filter.ids": ",".join(ids[i:i + 50]), "fields": fields, "pageSize": 100, "format": "json"})
        req = urllib.request.Request("https://clinicaltrials.gov/api/v2/studies?" + q, headers={"User-Agent": "kg-repurposing-audit"})
        for s in json.load(urllib.request.urlopen(req, timeout=60)).get("studies", []):
            got[s["protocolSection"]["identificationModule"]["nctId"]] = s
        time.sleep(1)
    json.dump(got, open(f"{OUT}/ctgov_cases.json", "w")); print("fetched", len(got), "of", len(ids))

elif stage == "dossier":
    raw = json.load(open(f"{OUT}/ctgov_cases.json")); lines = []
    for g in ("hetionet", "primekg"):
        for k, p in enumerate(json.load(open(f"{OUT}/pairs_{g}.json")), 1):
            lines.append(f"### CASE {'H' if g == 'hetionet' else 'P'}{k:02d} | graph={g} | drug={p['drug']} | disease={p['disease']} ({p['disease_id']}) | "
                         f"label={p['label']} | registered trials started before 2015 for this pair: {p['n_trials_before_2015']}")
            for t in p["trials"]:
                ps = raw.get(t["nct"], {}).get("protocolSection", {}); dm = ps.get("designModule", {}); am = ps.get("armsInterventionsModule", {})
                el = ps.get("eligibilityModule", {}); sm = ps.get("statusModule", {}); idm = ps.get("identificationModule", {})
                lines.append(f"- {t['nct']} | phase {','.join(dm.get('phases', []))} | allocation {dm.get('designInfo', {}).get('allocation', '')} | "
                             f"enrollment {dm.get('enrollmentInfo', {}).get('count', '')} | status {sm.get('overallStatus', '')} | has results: {raw.get(t['nct'], {}).get('hasResults')}")
                lines.append(f"  Title: {idm.get('officialTitle') or idm.get('briefTitle', '')}")
                lines.append(f"  Conditions: {'; '.join(ps.get('conditionsModule', {}).get('conditions', []))}")
                lines.append(f"  Why stopped: {sm.get('whyStopped', t.get('why_stopped'))} | Open Targets category: {t['stop_categories']}")
                for a in am.get("armGroups", []):
                    lines.append(f"  Arm [{a.get('type', '')}] {a.get('label', '')}: {', '.join(a.get('interventionNames', []))} — {(a.get('description') or '')[:250]}")
                po = ps.get("outcomesModule", {}).get("primaryOutcomes", [])
                if po:
                    lines.append(f"  Primary outcome: {po[0].get('measure', '')[:200]}")
                lines.append(f"  Ages {el.get('minimumAge', '')}–{el.get('maximumAge', '')}; eligibility excerpt: {(el.get('eligibilityCriteria') or '')[:400].replace(chr(10), ' ')}")
            lines.append("")
    open(f"{OUT}/dossiers.md", "w").write("\n".join(lines)); print("dossiers written")
