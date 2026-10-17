#!/usr/bin/env python
"""Literature survey: how often do drug-repurposing studies use stopped trials as negative labels?

Two Europe PMC searches (full text where available):
  Q1  articles that mention repoDB, whose failed pairs are terminated, withdrawn or suspended trials;
  Q2  drug-repurposing articles that mention terminated, withdrawn or failed trials together with
      negative samples, labels or sets, without mentioning repoDB.
For each open-access hit with a PMC full text, the paragraphs that mention repoDB, stopped or failed
trials or negative examples are extracted for screening (screening codes are assigned separately and
stored in literature_screening.tsv).

Usage: 24_literature_survey.py OUT_DIR
Outputs: records.tsv (all hits), passages.jsonl (extracted paragraphs per open-access article),
         query_log.json (queries, dates and hit counts). Idempotent: full texts are cached in OUT_DIR/xml.
"""
import datetime, json, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
import pandas as pd

OUT = sys.argv[1] if len(sys.argv) > 1 else "work/literature"
os.makedirs(f"{OUT}/xml", exist_ok=True)
API = "https://www.ebi.ac.uk/europepmc/webservices/rest"
Q1 = '"repoDB"'
Q2 = ('("drug repurposing" OR "drug repositioning" OR "drug-disease association" OR "drug-disease associations" OR '
      '"drug indication") AND (terminated OR withdrawn OR "failed clinical trials" OR "failed trials" OR '
      '"clinical trial failures" OR "failed drug trials") AND ("negative samples" OR "negative examples" OR '
      '"negative pairs" OR "negative labels" OR "true negatives" OR "negative set" OR "negative data") NOT "repoDB"')
KEY = re.compile(r"repodb|terminat|withdrawn|suspended|failed (?:clinical )?trial|trial failure|negative (?:sample|example|pair|label|set|data|control)|true negative",
                 re.I)


def get(url, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=90).read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def search(q):
    rows, cursor = [], "*"
    while True:
        u = f"{API}/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 1000, "cursorMark": cursor,
                                                       "resultType": "lite"})
        d = json.loads(get(u))
        rows += d["resultList"]["result"]
        if d.get("nextCursorMark") in (None, cursor) or not d["resultList"]["result"]:
            return rows, d["hitCount"]
        cursor = d["nextCursorMark"]


def paragraphs(xml_bytes):
    root = ET.fromstring(xml_bytes)
    body = root.find(".//body")
    out = []
    for sec in (body.iter() if body is not None else []):
        if sec.tag in ("p", "td", "caption"):
            txt = re.sub(r"\s+", " ", "".join(sec.itertext())).strip()
            if txt:
                out.append(txt)
    return out


log = {"date": datetime.date.today().isoformat(), "queries": {}}
recs = []
for name, q in [("Q1", Q1), ("Q2", Q2)]:
    rows, n = search(q)
    log["queries"][name] = {"query": q, "hits": n}
    for r in rows:
        recs.append({"query": name, "id": r.get("id"), "source": r.get("source"), "pmid": r.get("pmid"), "pmcid": r.get("pmcid"),
                     "doi": r.get("doi"), "title": r.get("title"), "journal": r.get("journalTitle"),
                     "year": r.get("pubYear"), "pub_type": r.get("pubType"), "open_access": r.get("isOpenAccess"),
                     "in_pmc": r.get("inPMC")})
rec = pd.DataFrame(recs)
rec.to_csv(f"{OUT}/records.tsv", sep="\t", index=False)
with open(f"{OUT}/passages.jsonl", "w") as fh:
    for r in rec.itertuples(index=False):
        if not isinstance(r.pmcid, str) or r.open_access != "Y":
            continue
        path = f"{OUT}/xml/{r.pmcid}.xml"
        if not os.path.exists(path):
            try:
                open(path, "wb").write(get(f"{API}/{r.pmcid}/fullTextXML"))
            except Exception as e:  # no full text deposited
                open(path, "wb").write(b"")
                print("no full text", r.pmcid, e)
            time.sleep(0.3)
        raw = open(path, "rb").read()
        if not raw:
            continue
        try:
            paras = paragraphs(raw)
        except ET.ParseError:
            continue
        hits = [p for p in paras if KEY.search(p)]
        fh.write(json.dumps({"query": r.query, "pmcid": r.pmcid, "title": r.title, "year": r.year, "journal": r.journal,
                             "n_paragraphs": len(paras), "passages": [p[:1500] for p in hits][:25]}) + "\n")
log["records"] = int(len(rec))
log["open_access_pmc"] = int(((rec.open_access == "Y") & rec.pmcid.notna()).sum())
json.dump(log, open(f"{OUT}/query_log.json", "w"), indent=1)
print(json.dumps(log, indent=1))
