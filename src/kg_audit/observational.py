"""Descriptive overlap-weighted associations for a locally supplied cohort.

This module does not establish incident use, exchangeability, valid covariate
timing or a causal effect. It writes aggregate estimates only, never patient rows.
"""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

FULL_CONT = [
    "age",
    "sofa_resp",
    "sofa_coag",
    "sofa_liver",
    "sofa_cardio",
    "sofa_cns",
    "sofa_renal",
    "lactate",
    "comorb_count",
]
FULL_BIN = [
    "male",
    "emergency_adm",
    "vent",
    "strong_vaso",
    "other_vaso",
    "diabetes",
    "ckd",
    "chf",
    "copd",
    "liver_dis",
    "malignancy",
    "htn",
    "cad",
]


def landmark_records(frame):
    """Classify orders using only information through 48 hours.

    Exact hospital death times take precedence. Date-only deaths are intervals;
    records crossing an eligibility/outcome boundary are flagged, not assigned an
    invented death hour. Follow-up ends 28 days after ICU entry, not 28 after entry
    into the landmark risk set.
    """
    out = frame.copy()
    for c in ["t0", "pre_t0_statin", "first_post_statin", "deathtime", "dod", "dischtime"]:
        out[c] = pd.to_datetime(out[c], errors="coerce")
    landmark = out.t0 + pd.Timedelta(hours=48)
    horizon = out.t0 + pd.Timedelta(days=28)
    exact = out.deathtime.notna()
    dated = (~exact) & out.dod.notna()
    lower = out.deathtime.where(exact, out.dod)
    upper = lower.where(exact, lower + pd.Timedelta(days=1))
    out["treat"] = ((out.first_post_statin >= out.t0) & (out.first_post_statin <= landmark)).astype(int)
    out["prevalent_user"] = out.pre_t0_statin.notna().astype(int)
    out["ambiguous_landmark"] = (dated & (lower <= landmark) & (upper > landmark)).astype(int)
    out["ambiguous_horizon"] = (dated & (lower <= horizon) & (upper > horizon)).astype(int)
    out["alive_48h"] = (lower.isna() | (lower > landmark)).astype(int)
    out["hospitalized_48h"] = (out.dischtime >= landmark).astype(int)
    out["event28"] = ((lower > landmark) & (upper <= horizon)).astype(int)
    return out


def estimate(frame, mode):
    cont = ["age"] if mode == "demographic" else FULL_CONT
    binary = ["male", "emergency_adm"] if mode == "demographic" else FULL_BIN
    columns = []
    for name in cont:
        x = pd.to_numeric(frame[name], errors="coerce")
        median = x.median()
        if pd.isna(median):
            raise ValueError(f"Covariate {name} has no observed values")
        columns.append(x.fillna(median).to_numpy())
    for name in binary:
        columns.append(pd.to_numeric(frame[name], errors="raise").fillna(0).to_numpy())
    for group in ["Black", "Hispanic", "Asian", "Other"]:
        columns.append((frame.race_grp == group).to_numpy(dtype=float))
    age = frame.age.to_numpy(dtype=float)
    columns.extend([age**2, age**3])
    x = np.column_stack(columns).astype(float)
    x = StandardScaler().fit_transform(x)
    t = frame.treat.to_numpy(dtype=int)
    y = frame.event28.to_numpy(dtype=int)
    if not np.isin(t, [0, 1]).all() or not np.isin(y, [0, 1]).all():
        raise ValueError("Treatment/outcome indicators must be binary")
    fit = LogisticRegression(C=0.5, solver="liblinear", intercept_scaling=100, max_iter=5000).fit(x, t)
    logits = np.sum(x * fit.coef_[0], axis=1) + fit.intercept_[0]
    ps = 1 / (1 + np.exp(-np.clip(logits, -30, 30)))
    ps = np.clip(ps, 0.001, 0.999)
    w = np.where(t == 1, 1 - ps, ps)
    r1 = np.sum(w[t == 1] * y[t == 1]) / np.sum(w[t == 1])
    r0 = np.sum(w[t == 0] * y[t == 0]) / np.sum(w[t == 0])
    return {
        "risk_order": float(r1),
        "risk_comparison": float(r0),
        "risk_ratio": float(r1 / r0),
        "risk_difference": float(r1 - r0),
        "ess_order": float(w[t == 1].sum() ** 2 / np.sum(w[t == 1] ** 2)),
        "ess_comparison": float(w[t == 0].sum() ** 2 / np.sum(w[t == 0] ** 2)),
    }


def run(cohort_path, code_status_path, output, bootstrap=1000, seed=20261005, database=None):
    full = pd.read_csv(cohort_path)
    frame = full.copy()
    if code_status_path is not None:
        status = pd.read_csv(code_status_path)
        if status.stay_id.duplicated().any():
            raise ValueError("Duplicate stay in code-status input")
        frame = full.merge(status, on="stay_id", how="left", validate="one_to_one")
    if database is None:
        raise ValueError("A local cohort database is required to rebuild exposure and outcome timing")
    import duckdb

    con = duckdb.connect(str(database), read_only=True)
    timing = con.execute(
        "SELECT s.stay_id,s.t0,s.pre_t0_statin,s.first_post_statin,c.deathtime,c.dod,c.dischtime,c.anchor_age,c.anchor_year FROM statin_exposure s JOIN cohort c USING(stay_id)"
    ).df()
    con.close()
    timing = landmark_records(timing)
    columns = ["treat", "prevalent_user", "alive_48h", "event28"]
    frame = frame.drop(columns=columns).merge(
        timing[
            ["stay_id"]
            + columns
            + ["ambiguous_landmark", "ambiguous_horizon", "hospitalized_48h", "anchor_age", "anchor_year"]
        ],
        on="stay_id",
        validate="one_to_one",
    )
    frame["age"] = frame.anchor_age + pd.to_datetime(frame.t0).dt.year - frame.anchor_year
    base = (frame.sepsis3 == 1) & (frame.age >= 18)
    stages = {"sepsis_flag_and_adult_in_existing_cohort": int(base.sum())}
    eligible = base & (frame.prevalent_user == 0)
    stages["without_observed_pre_icu_order"] = int(eligible.sum())
    eligible &= (frame.hospitalized_48h == 1) & (frame.alive_48h == 1)
    stages["hospitalized_and_recorded_alive_at_landmark"] = int(eligible.sum())
    eligible &= (frame.ambiguous_landmark == 0) & (frame.ambiguous_horizon == 0)
    stages["unambiguous_available_death_timing"] = int(eligible.sum())
    if code_status_path is not None:
        eligible &= frame.comfort_care_48h.fillna(0) == 0
        stages["optional_recorded_comfort_care_exclusion"] = int(eligible.sum())
    frame = frame.loc[eligible].reset_index(drop=True)
    if frame.subject_id.duplicated().any():
        raise ValueError("Patient bootstrap expects one eligible stay per patient")
    if frame.age.isna().any():
        raise ValueError("Missing age requires explicit handling")
    result = {
        "n": len(frame),
        "n_order": int((frame.treat == 1).sum()),
        "n_comparison": int((frame.treat == 0).sum()),
        "events_order": int(frame.loc[frame.treat == 1, "event28"].sum()),
        "events_comparison": int(frame.loc[frame.treat == 0, "event28"].sum()),
        "bootstrap": bootstrap,
        "seed": seed,
        "cohort_attrition": stages,
        "estimand": "Descriptive overlap-population contrast of orders by 48 hours versus no order by 48 hours, among eligible patients still hospitalized and recorded alive at 48 hours. Subsequent initiation does not change the comparator assignment. Recorded deaths are counted through day 28 from ICU entry.",
        "limitations": [
            "No pre-admission washout established",
            "Order is not administration",
            "Landmark selection changes population",
            "Full adjustment includes measurements potentially after an early order",
            "Demographic adjustment leaves substantial confounding",
        ],
        "input_sha256": hashlib.sha256(Path(cohort_path).read_bytes()).hexdigest(),
        "models": {},
    }
    for mode in ["demographic", "full"]:
        point = estimate(frame, mode)
        rng = np.random.default_rng(seed)
        samples = []
        for i in range(bootstrap):
            selected = rng.integers(len(frame), size=len(frame))
            value = estimate(frame.iloc[selected], mode)
            samples.append([value["risk_ratio"], value["risk_difference"]])
            if (i + 1) % 200 == 0:
                print(f"{mode}: {i + 1}/{bootstrap} resamples", flush=True)
        samples = np.array(samples)
        point["risk_ratio_CI95"] = np.quantile(samples[:, 0], [0.025, 0.975]).tolist()
        point["risk_difference_CI95"] = np.quantile(samples[:, 1], [0.025, 0.975]).tolist()
        result["models"][mode] = point
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort", required=True)
    p.add_argument("--code-status", help="Optional sensitivity exclusion; not used in the primary descriptive rerun")
    p.add_argument("--database", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--bootstrap", type=int, default=1000)
    a = p.parse_args()
    run(a.cohort, a.code_status, a.output, a.bootstrap, database=a.database)


if __name__ == "__main__":
    main()
