#!/usr/bin/env python
"""Local MIMIC-IV extraction supporting a descriptive, order-defined illustration.
Credentialed data are required. Legacy field names do not establish incident use
or causal identification; final timing is rebuilt by kg_audit.observational.
"""

import duckdb, time

ROOT = __import__("os").environ["KG_AUDIT_WORKDIR"]
M = f"{ROOT}/data/mimic-iv-3.1"
OUT = f"{ROOT}/outputs/mimic"
con = duckdb.connect(f"{OUT}/mimic.duckdb")
con.execute("PRAGMA threads=6")


def R(t, s):
    return f"read_csv_auto('{M}/{s}/{t}.csv.gz', types={{'icd_code':'VARCHAR'}}, ignore_errors=false)"


t0 = time.time()

# ---------- comorbidities from diagnoses_icd (ICD-9 + ICD-10 prefixes) ----------
con.execute(f"""
CREATE OR REPLACE TABLE comorb AS
WITH dx AS (SELECT hadm_id, icd_version, icd_code FROM {R("diagnoses_icd", "hosp")})
SELECT hadm_id,
 max(CASE WHEN (icd_version=9 AND starts_with(icd_code,'250')) OR (icd_version=10 AND icd_code SIMILAR TO 'E1[0-4].*') THEN 1 ELSE 0 END) AS diabetes,
 max(CASE WHEN (icd_version=9 AND starts_with(icd_code,'585')) OR (icd_version=10 AND starts_with(icd_code,'N18')) THEN 1 ELSE 0 END) AS ckd,
 max(CASE WHEN (icd_version=9 AND starts_with(icd_code,'428')) OR (icd_version=10 AND starts_with(icd_code,'I50')) THEN 1 ELSE 0 END) AS chf,
 max(CASE WHEN (icd_version=9 AND (starts_with(icd_code,'490') OR starts_with(icd_code,'491') OR starts_with(icd_code,'492') OR starts_with(icd_code,'496'))) OR (icd_version=10 AND starts_with(icd_code,'J44')) THEN 1 ELSE 0 END) AS copd,
 max(CASE WHEN (icd_version=9 AND starts_with(icd_code,'571')) OR (icd_version=10 AND icd_code SIMILAR TO 'K7[0-4].*') THEN 1 ELSE 0 END) AS liver_dis,
 max(CASE WHEN (icd_version=9 AND (icd_code >= '140' AND icd_code < '209')) OR (icd_version=10 AND icd_code SIMILAR TO 'C[0-8].*') THEN 1 ELSE 0 END) AS malignancy,
 max(CASE WHEN (icd_version=9 AND starts_with(icd_code,'401')) OR (icd_version=10 AND starts_with(icd_code,'I10')) THEN 1 ELSE 0 END) AS htn,
 max(CASE WHEN (icd_version=9 AND (starts_with(icd_code,'410') OR starts_with(icd_code,'412') OR starts_with(icd_code,'414'))) OR (icd_version=10 AND (starts_with(icd_code,'I21') OR starts_with(icd_code,'I25'))) THEN 1 ELSE 0 END) AS cad
FROM dx GROUP BY hadm_id
""")
print(f"[comorb] rows={con.sql('SELECT count(*) FROM comorb').fetchone()[0]} ({time.time() - t0:.0f}s)", flush=True)

# ---------- assemble analytic cohort ----------
con.execute("""
CREATE OR REPLACE TABLE analytic AS
SELECT c.subject_id,c.stay_id,c.hadm_id,c.t0,
  c.anchor_age AS age, CASE WHEN c.gender='M' THEN 1 ELSE 0 END AS male,
  CASE WHEN lower(c.race) LIKE '%white%' THEN 'White'
       WHEN lower(c.race) LIKE '%black%' THEN 'Black'
       WHEN lower(c.race) LIKE '%hispan%' THEN 'Hispanic'
       WHEN lower(c.race) LIKE '%asian%' THEN 'Asian' ELSE 'Other' END AS race_grp,
  CASE WHEN c.admission_type LIKE '%EMER%' OR c.admission_type LIKE '%URGENT%' THEN 1 ELSE 0 END AS emergency_adm,
  c.anchor_year AS cal_year,
  s.sofa_total, s.sofa_resp,s.sofa_coag,s.sofa_liver,s.sofa_cardio,s.sofa_cns,s.sofa_renal,
  s.strong_vaso, s.other_vaso, s.vent, s.lactate, s.creat, s.bili, s.platelets,
  COALESCE(cm.diabetes,0) diabetes, COALESCE(cm.ckd,0) ckd, COALESCE(cm.chf,0) chf,
  COALESCE(cm.copd,0) copd, COALESCE(cm.liver_dis,0) liver_dis, COALESCE(cm.malignancy,0) malignancy,
  COALESCE(cm.htn,0) htn, COALESCE(cm.cad,0) cad,
  (COALESCE(cm.diabetes,0)+COALESCE(cm.ckd,0)+COALESCE(cm.chf,0)+COALESCE(cm.copd,0)+COALESCE(cm.liver_dis,0)+COALESCE(cm.malignancy,0)+COALESCE(cm.htn,0)+COALESCE(cm.cad,0)) AS comorb_count,
  se.prevalent_user, se.statin_initiator, se.never_statin,
  COALESCE(si.suspected_infection,0) AS suspected_infection,
  COALESCE(sic.sepsis_icd,0) AS sepsis_icd,
  o.days_to_death, o.dead28, o.death_dt
FROM cohort c
LEFT JOIN sofa s USING(stay_id)
LEFT JOIN comorb cm ON cm.hadm_id=c.hadm_id
LEFT JOIN statin_exposure se USING(stay_id)
LEFT JOIN suspicion si USING(stay_id)
LEFT JOIN sepsis_icd sic ON sic.hadm_id=c.hadm_id
LEFT JOIN outcome2 o USING(stay_id)
""")

# treatment arm + landmark(48h) survival
con.execute("""
CREATE OR REPLACE TABLE analytic2 AS
SELECT *,
  CASE WHEN statin_initiator=1 THEN 1 WHEN never_statin=1 THEN 0 ELSE NULL END AS treat,
  -- Sepsis-3 (primary) and ICD (sensitivity)
  CASE WHEN suspected_infection=1 AND sofa_total>=2 THEN 1 ELSE 0 END AS sepsis3,
  -- alive at 48h landmark?
  CASE WHEN days_to_death IS NULL OR days_to_death>2 THEN 1 ELSE 0 END AS alive_48h,
  -- event within 28d (from t0); survival time from 48h landmark (cap 28d)
  CASE WHEN days_to_death IS NOT NULL AND days_to_death<=28 AND days_to_death>2 THEN 1 ELSE 0 END AS event28,
  CASE WHEN days_to_death IS NOT NULL AND days_to_death<=28 AND days_to_death>2 THEN days_to_death-2
       ELSE 26 END AS surv_days
FROM analytic
""")
con.sql(f"COPY analytic2 TO '{OUT}/analytic_full.csv' (HEADER, DELIMITER ',')")


# report cohort sizes
def rep(name, where):
    q = f"""SELECT count(*) n, sum(treat) treated, sum(CASE WHEN treat=0 THEN 1 ELSE 0 END) control,
          avg(CASE WHEN treat=1 THEN event28 END) m_t, avg(CASE WHEN treat=0 THEN event28 END) m_c
          FROM analytic2 WHERE {where}"""
    r = con.sql(q).df().iloc[0]
    print(
        f"[{name}] n={r.n:.0f} treated={r.treated:.0f} control={r.control:.0f} | "
        f"28d mort treated={r.m_t:.3f} control={r.m_c:.3f}"
    )


print("\n=== cohort sizes (treat=initiator vs never; alive at 48h) ===")
rep("Sepsis-3 primary", "sepsis3=1 AND treat IS NOT NULL AND alive_48h=1 AND prevalent_user=0")
rep("Sepsis ICD sens.", "sepsis_icd=1 AND treat IS NOT NULL AND alive_48h=1 AND prevalent_user=0")
print("\n=== NAIVE (any statin incl. prevalent vs never; no landmark/adjust) ===")
naive = (
    con.sql("""SELECT
  avg(CASE WHEN (statin_initiator=1 OR prevalent_user=1) THEN dead28 END) m_statin,
  avg(CASE WHEN never_statin=1 THEN dead28 END) m_never,
  sum(CASE WHEN (statin_initiator=1 OR prevalent_user=1) THEN 1 ELSE 0 END) n_statin,
  sum(CASE WHEN never_statin=1 THEN 1 ELSE 0 END) n_never
  FROM analytic2 WHERE sepsis3=1""")
    .df()
    .iloc[0]
)
print(
    f"Sepsis-3 NAIVE: statin 28d mort={naive.m_statin:.3f} (n={naive.n_statin:.0f}) vs never={naive.m_never:.3f} (n={naive.n_never:.0f})"
)
print(f"DONE assemble in {time.time() - t0:.0f}s")
