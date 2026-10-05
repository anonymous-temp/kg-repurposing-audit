#!/usr/bin/env python
"""Local MIMIC-IV extraction supporting a descriptive, order-defined illustration.
Credentialed data are required. Legacy field names do not establish incident use
or causal identification; final timing is rebuilt by kg_audit.observational.
"""
import duckdb, time
ROOT=__import__("os").environ["KG_AUDIT_WORKDIR"]; M=f"{ROOT}/data/mimic-iv-3.1"; OUT=f"{ROOT}/outputs/mimic"
con=duckdb.connect(f"{OUT}/mimic.duckdb"); con.execute("PRAGMA threads=6")
def R(t,s): return f"read_csv_auto('{M}/{s}/{t}.csv.gz', ignore_errors=false)"
t0=time.time()

# windows table
con.execute("CREATE OR REPLACE TABLE win AS SELECT stay_id,hadm_id,t0, t0 - INTERVAL 6 HOUR AS w0, t0 + INTERVAL 24 HOUR AS w1 FROM cohort")

# ---------- LABS (labevents keyed by hadm_id) ----------
print("scanning labevents...",flush=True)
con.execute(f"""
CREATE OR REPLACE TABLE labs AS
SELECT w.stay_id,
  max(CASE WHEN le.itemid=50912 THEN le.valuenum END) AS creat,
  max(CASE WHEN le.itemid=50885 THEN le.valuenum END) AS bili,
  min(CASE WHEN le.itemid=51265 THEN le.valuenum END) AS platelets,
  max(CASE WHEN le.itemid=50813 THEN le.valuenum END) AS lactate,
  min(CASE WHEN le.itemid=50821 THEN le.valuenum END) AS pao2
FROM win w JOIN {R('labevents','hosp')} le
  ON le.hadm_id=w.hadm_id AND le.itemid IN (50912,50885,51265,50813,50821)
  AND le.charttime BETWEEN w.w0 AND w.w1 AND le.valuenum IS NOT NULL
GROUP BY w.stay_id
""")
print(f"  labs rows={con.sql('SELECT count(*) FROM labs').fetchone()[0]} ({time.time()-t0:.0f}s)",flush=True)

# ---------- CHART (chartevents keyed by stay_id) ----------
print("scanning chartevents (large)...",flush=True)
con.execute(f"""
CREATE OR REPLACE TABLE chart_raw AS
SELECT ce.stay_id, ce.charttime, ce.itemid, ce.valuenum
FROM win w JOIN {R('chartevents','icu')} ce
  ON ce.stay_id=w.stay_id AND ce.itemid IN (220052,220181,220739,223900,223901,223835,220277,223849)
  AND ce.charttime BETWEEN w.w0 AND w.w1 AND ce.valuenum IS NOT NULL
""")
# GCS total per charttime then min per stay
con.execute("""
CREATE OR REPLACE TABLE gcs AS
WITH g AS (
  SELECT stay_id, charttime,
    sum(CASE WHEN itemid IN (220739,223900,223901) THEN valuenum END) AS gcs_t,
    count(*) FILTER (WHERE itemid IN (220739,223900,223901)) AS ncomp
  FROM chart_raw GROUP BY stay_id,charttime)
SELECT stay_id, min(CASE WHEN ncomp=3 THEN gcs_t END) AS gcs_min FROM g GROUP BY stay_id
""")
con.execute("""
CREATE OR REPLACE TABLE chartagg AS
SELECT stay_id,
  min(CASE WHEN itemid IN (220052,220181) THEN valuenum END) AS map_min,
  max(CASE WHEN itemid=223835 THEN valuenum END) AS fio2_max,
  min(CASE WHEN itemid=220277 THEN valuenum END) AS spo2_min,
  max(CASE WHEN itemid=223849 THEN 1 ELSE 0 END) AS vent
FROM chart_raw GROUP BY stay_id
""")
print(f"  chart rows={con.sql('SELECT count(*) FROM chart_raw').fetchone()[0]} ({time.time()-t0:.0f}s)",flush=True)

# ---------- VASOPRESSORS (inputevents keyed by stay_id) ----------
print("scanning inputevents...",flush=True)
con.execute(f"""
CREATE OR REPLACE TABLE vaso AS
SELECT w.stay_id,
  max(CASE WHEN ie.itemid IN (221906,221289) THEN 1 ELSE 0 END) AS strong_vaso,   -- norepi/epi
  max(CASE WHEN ie.itemid IN (221662,221653,222315,221749) THEN 1 ELSE 0 END) AS other_vaso
FROM win w JOIN {R('inputevents','icu')} ie
  ON ie.stay_id=w.stay_id AND ie.itemid IN (221906,221289,221662,221653,222315,221749)
  AND ie.starttime BETWEEN w.w0 AND w.w1
GROUP BY w.stay_id
""")
print(f"  vaso rows={con.sql('SELECT count(*) FROM vaso').fetchone()[0]} ({time.time()-t0:.0f}s)",flush=True)

# ---------- SOFA scoring ----------
con.execute("""
CREATE OR REPLACE TABLE sofa AS
WITH j AS (
  SELECT c.stay_id, l.creat,l.bili,l.platelets,l.lactate,l.pao2,
         ca.map_min, ca.fio2_max, ca.spo2_min, ca.vent, g.gcs_min,
         v.strong_vaso, v.other_vaso
  FROM cohort c
  LEFT JOIN labs l USING(stay_id)
  LEFT JOIN chartagg ca USING(stay_id)
  LEFT JOIN gcs g USING(stay_id)
  LEFT JOIN vaso v USING(stay_id)
), s AS (
  SELECT *,
    -- respiration via PaO2/FiO2 (FiO2 % -> fraction); missing -> 0
    CASE WHEN pao2 IS NOT NULL AND fio2_max IS NOT NULL AND fio2_max>0 THEN
      (CASE WHEN pao2/(CASE WHEN fio2_max>1 THEN fio2_max/100.0 ELSE fio2_max END) < 100 AND vent=1 THEN 4
            WHEN pao2/(CASE WHEN fio2_max>1 THEN fio2_max/100.0 ELSE fio2_max END) < 200 AND vent=1 THEN 3
            WHEN pao2/(CASE WHEN fio2_max>1 THEN fio2_max/100.0 ELSE fio2_max END) < 300 THEN 2
            WHEN pao2/(CASE WHEN fio2_max>1 THEN fio2_max/100.0 ELSE fio2_max END) < 400 THEN 1 ELSE 0 END)
      ELSE 0 END AS sofa_resp,
    CASE WHEN platelets IS NULL THEN 0 WHEN platelets<20 THEN 4 WHEN platelets<50 THEN 3
         WHEN platelets<100 THEN 2 WHEN platelets<150 THEN 1 ELSE 0 END AS sofa_coag,
    CASE WHEN bili IS NULL THEN 0 WHEN bili>=12 THEN 4 WHEN bili>=6 THEN 3
         WHEN bili>=2 THEN 2 WHEN bili>=1.2 THEN 1 ELSE 0 END AS sofa_liver,
    CASE WHEN strong_vaso=1 THEN 3 WHEN other_vaso=1 THEN 2
         WHEN map_min IS NOT NULL AND map_min<70 THEN 1 ELSE 0 END AS sofa_cardio,
    CASE WHEN gcs_min IS NULL THEN 0 WHEN gcs_min<6 THEN 4 WHEN gcs_min<10 THEN 3
         WHEN gcs_min<13 THEN 2 WHEN gcs_min<15 THEN 1 ELSE 0 END AS sofa_cns,
    CASE WHEN creat IS NULL THEN 0 WHEN creat>=5 THEN 4 WHEN creat>=3.5 THEN 3
         WHEN creat>=2 THEN 2 WHEN creat>=1.2 THEN 1 ELSE 0 END AS sofa_renal
  FROM j
)
SELECT *, (sofa_resp+sofa_coag+sofa_liver+sofa_cardio+sofa_cns+sofa_renal) AS sofa_total
FROM s
""")
d=con.sql("SELECT avg(sofa_total) m, median(sofa_total) md, sum(CASE WHEN sofa_total>=2 THEN 1 ELSE 0 END) ge2, count(*) n FROM sofa").df()
print(f"[SOFA] mean={d.m[0]:.2f} median={d.md[0]:.1f} SOFA>=2: {d.ge2[0]:.0f}/{d.n[0]:.0f}",flush=True)
print(con.sql("SELECT round(avg(sofa_resp),2) resp,round(avg(sofa_coag),2) coag,round(avg(sofa_liver),2) liver,round(avg(sofa_cardio),2) cardio,round(avg(sofa_cns),2) cns,round(avg(sofa_renal),2) renal FROM sofa").df().to_string(index=False))
print(f"DONE sofa in {time.time()-t0:.0f}s")
