#!/usr/bin/env python
"""Local MIMIC-IV extraction supporting a descriptive, order-defined illustration.
Credentialed data are required. Legacy field names do not establish incident use
or causal identification; final timing is rebuilt by kg_audit.observational.
"""
import duckdb, os, time
ROOT=__import__("os").environ["KG_AUDIT_WORKDIR"]; M=f"{ROOT}/data/mimic-iv-3.1"; OUT=f"{ROOT}/outputs/mimic"
os.makedirs(OUT,exist_ok=True)
con=duckdb.connect(f"{OUT}/mimic.duckdb")
con.execute("PRAGMA threads=6")
t0=time.time()

def R(tbl, sub="hosp"):
    types = ", types={'icd_code':'VARCHAR'}" if tbl == "diagnoses_icd" else ""
    return f"read_csv_auto('{M}/{sub}/{tbl}.csv.gz'{types}, ignore_errors=false)"

# ---------- Stage 1: base cohort (first ICU stay, adults) ----------
con.execute(f"""
CREATE OR REPLACE TABLE icu0 AS
WITH ic AS (
  SELECT subject_id,hadm_id,stay_id,first_careunit,intime,outtime,los,
         row_number() OVER (PARTITION BY subject_id ORDER BY intime) rn
  FROM {R('icustays','icu')}
)
SELECT * FROM ic WHERE rn=1
""")
con.execute(f"""
CREATE OR REPLACE TABLE cohort AS
SELECT c.subject_id, c.hadm_id, c.stay_id, c.first_careunit,
       c.intime AS t0, c.outtime, c.los,
       p.gender, p.anchor_age, p.anchor_year, p.dod,
       a.admittime, a.dischtime, a.deathtime, a.admission_type, a.race, a.hospital_expire_flag,
       date_diff('year', a.admittime::DATE, c.intime::DATE) AS yrs_from_adm
FROM icu0 c
JOIN {R('patients')} p USING(subject_id)
JOIN {R('admissions')} a USING(subject_id,hadm_id)
WHERE p.anchor_age >= 18
""")
n=con.sql("SELECT count(*) FROM cohort").fetchone()[0]
print(f"[1] base cohort (first ICU, adults): {n}", flush=True)

# ---------- 28-day mortality + survival time from t0 ----------
# death datetime: prefer deathtime (in-hospital); else dod (date). days from t0.
con.execute("""
CREATE OR REPLACE TABLE outcome AS
SELECT stay_id, subject_id, hadm_id, t0,
  COALESCE(deathtime, CASE WHEN dod IS NOT NULL THEN dod::TIMESTAMP END) AS death_dt
FROM cohort
""")
con.execute("""
CREATE OR REPLACE TABLE outcome2 AS
SELECT *,
  CASE WHEN death_dt IS NOT NULL THEN date_diff('day', t0, death_dt) END AS days_to_death,
  CASE WHEN death_dt IS NOT NULL AND date_diff('day', t0, death_dt) BETWEEN 0 AND 28 THEN 1 ELSE 0 END AS dead28
FROM outcome
""")
d=con.sql("SELECT count(*) n, sum(dead28) d, avg(dead28) FROM outcome2").df()
print(f"[1] 28-day deaths: {d.iloc[0,1]:.0f}/{d.iloc[0,0]:.0f} ({d.iloc[0,2]:.3f})", flush=True)

# ---------- Stage 2: statin exposure (new-user proxy) ----------
STATIN = r"(atorvastatin|simvastatin|pravastatin|rosuvastatin|lovastatin|fluvastatin|pitavastatin|livalo|lipitor|crestor|zocor|pravachol|mevacor|lescol|altoprev)"
con.execute(f"""
CREATE OR REPLACE TABLE statin_orders AS
SELECT subject_id,hadm_id, lower(drug) drug, starttime, stoptime, route
FROM {R('prescriptions')}
WHERE regexp_matches(lower(drug),'{STATIN}')
  AND NOT regexp_matches(lower(drug),'(nystatin|cilastatin|pentostatin|sandostatin|somatostatin|mycostatin)')
""")
# join to cohort, classify relative to t0 with grace window 48h
con.execute("""
CREATE OR REPLACE TABLE statin_exposure AS
WITH s AS (
  SELECT c.stay_id, c.t0,
    min(CASE WHEN so.starttime < c.t0 THEN so.starttime END) AS pre_t0_statin,
    min(CASE WHEN so.starttime >= c.t0 THEN so.starttime END) AS first_post_statin
  FROM cohort c LEFT JOIN statin_orders so
    ON so.hadm_id=c.hadm_id
  GROUP BY c.stay_id, c.t0
)
SELECT stay_id, t0, pre_t0_statin, first_post_statin,
  CASE WHEN pre_t0_statin IS NOT NULL THEN 1 ELSE 0 END AS prevalent_user,
  CASE WHEN pre_t0_statin IS NULL AND first_post_statin IS NOT NULL
            AND date_diff('hour', t0, first_post_statin) BETWEEN 0 AND 48 THEN 1 ELSE 0 END AS statin_initiator,
  CASE WHEN pre_t0_statin IS NULL AND first_post_statin IS NULL THEN 1 ELSE 0 END AS never_statin
FROM s
""")
e=con.sql("""SELECT sum(prevalent_user) prevalent, sum(statin_initiator) initiators,
             sum(never_statin) never, count(*) tot FROM statin_exposure""").df()
print(f"[2] statin: prevalent={e.prevalent[0]:.0f} initiators(0-48h)={e.initiators[0]:.0f} never={e.never[0]:.0f} tot={e.tot[0]:.0f}", flush=True)

# ---------- Stage 3a: simplified sepsis phenotype (ICD) ----------
con.execute(f"""
CREATE OR REPLACE TABLE sepsis_icd AS
SELECT DISTINCT hadm_id, 1 AS sepsis_icd FROM {R('diagnoses_icd')}
WHERE (icd_version=9 AND icd_code IN ('99591','99592','78552'))
   OR (icd_version=10 AND (starts_with(icd_code,'A41') OR icd_code IN ('R6520','R6521')))
""")
ns=con.sql("SELECT count(*) FROM sepsis_icd").fetchone()[0]
print(f"[3a] hadm with sepsis ICD: {ns}", flush=True)

# ---------- Stage 3b: suspicion of infection (antibiotic + culture near t0) ----------
ABX=r"(vancomycin|piperacillin|tazobactam|cefepime|ceftriaxone|ceftazidime|cefazolin|meropenem|imipenem|ertapenem|metronidazole|ciprofloxacin|levofloxacin|moxifloxacin|azithromycin|gentamicin|tobramycin|amikacin|ampicillin|amoxicillin|aztreonam|clindamycin|linezolid|daptomycin|doxycycline|nafcillin|oxacillin|penicillin|sulfamethoxazole|trimethoprim|ceftaroline|tigecycline|colistin)"
con.execute(f"""
CREATE OR REPLACE TABLE abx AS
SELECT c.stay_id, min(pr.starttime) AS first_abx
FROM cohort c JOIN {R('prescriptions')} pr ON pr.hadm_id=c.hadm_id
WHERE regexp_matches(lower(pr.drug),'{ABX}')
  AND pr.starttime BETWEEN c.t0 - INTERVAL 24 HOUR AND c.t0 + INTERVAL 72 HOUR
GROUP BY c.stay_id
""")
con.execute(f"""
CREATE OR REPLACE TABLE cult AS
SELECT c.stay_id, min(mb.charttime) AS first_cult
FROM cohort c JOIN {R('microbiologyevents')} mb ON mb.hadm_id=c.hadm_id
WHERE COALESCE(mb.charttime, mb.chartdate::TIMESTAMP) BETWEEN c.t0 - INTERVAL 24 HOUR AND c.t0 + INTERVAL 24 HOUR
GROUP BY c.stay_id
""")
con.execute("""
CREATE OR REPLACE TABLE suspicion AS
SELECT c.stay_id,
  CASE WHEN a.first_abx IS NOT NULL AND cu.first_cult IS NOT NULL THEN 1 ELSE 0 END AS suspected_infection
FROM cohort c LEFT JOIN abx a USING(stay_id) LEFT JOIN cult cu USING(stay_id)
""")
si=con.sql("SELECT sum(suspected_infection) s, count(*) t FROM suspicion").df()
print(f"[3b] suspected infection near t0: {si.s[0]:.0f}/{si.t[0]:.0f}", flush=True)
print(f"DONE stage1 in {time.time()-t0:.0f}s")
