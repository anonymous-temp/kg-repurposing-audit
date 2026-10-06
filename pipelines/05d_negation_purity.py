"""Size and purity of the negative set of each write-back policy (share of negated pairs that are
recorded treatments or approved indications), in Hetionet and PrimeKG."""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
out = {}
for g, data in [("hetionet", "work/mapped"), ("primekg", "work/mapped_primekg")]:
    os.environ["DATA"] = data
    import importlib
    from kg_audit import writeback as L
    L = importlib.reload(L)
    ev = L.load_evidence()
    if g == "hetionet":
        ctd = set(L.load_graph()[3])
    else:
        ctd = {tuple(p) for p in json.load(open(f"{data}/graph.json"))["indication"]}
    good = ctd | ev["approved"]
    out[g] = {}
    for name in ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped"]:
        neg = {p for p, a in L.policy_sets(ev, L.NAMED[name]).items() if a == "negate"}
        out[g][name] = {"negated_pairs": len(neg), "recorded": len(neg & ctd), "recorded_or_approved": len(neg & good),
                        "pct": round(100 * len(neg & good) / len(neg), 1)}
    print(g, json.dumps(out[g]))
json.dump(out, open("work/results_role/negation_purity.json", "w"), indent=1)
