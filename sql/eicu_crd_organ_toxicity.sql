-- =============================================================================
-- eICU-CRD (BigQuery Standard SQL) Schema-Only Organ Toxicity Cohort Template
-- Author: Goktug Ozkan, MD (ORCID: 0000-0002-5022-9124)
-- Scope: Public PhysioNet BigQuery schema template (`physionet-data.eicu_crd`).
--        Contains ZERO patient data or exported query results. Credentialed
--        PhysioNet Data Use Agreement (DUA) is required to run against live tables.
-- =============================================================================

WITH icu_stays AS (
  SELECT
    hospitalid,
    patientunitstayid,
    gender,
    SAFE_CAST(REGEXP_REPLACE(age, r'[^0-9]', '') AS INT64) AS age_years
  FROM `physionet-data.eicu_crd.patient`
  WHERE SAFE_CAST(REGEXP_REPLACE(age, r'[^0-9]', '') AS INT64) >= 18
),

icu_labs AS (
  SELECT
    stay.hospitalid,
    stay.patientunitstayid,
    lab.labresultoffset,
    CASE
      WHEN LOWER(lab.labname) = 'creatinine' THEN 'scr'
      WHEN LOWER(lab.labname) IN ('alt (sgpt)', 'alt') THEN 'alt'
      WHEN LOWER(lab.labname) IN ('ast (sgot)', 'ast') THEN 'ast'
      WHEN LOWER(lab.labname) IN ('alkaline phos.', 'alkaline phosphatase') THEN 'alp'
      WHEN LOWER(lab.labname) IN ('total bilirubin', 'bilirubin, total') THEN 'tbl'
    END AS analyte,
    lab.labresult AS valuenum
  FROM icu_stays AS stay
  INNER JOIN `physionet-data.eicu_crd.lab` AS lab
    ON stay.patientunitstayid = lab.patientunitstayid
  WHERE lab.labresult IS NOT NULL
    AND lab.labresult > 0
    AND lab.labresultoffset >= -1440
    AND LOWER(lab.labname) IN (
      'creatinine',
      'alt (sgpt)', 'alt',
      'ast (sgot)', 'ast',
      'alkaline phos.', 'alkaline phosphatase',
      'total bilirubin', 'bilirubin, total'
    )
),

stay_metrics AS (
  SELECT
    hospitalid,
    patientunitstayid,
    ARRAY_AGG(
      IF(analyte = 'scr', valuenum, NULL) IGNORE NULLS
      ORDER BY labresultoffset ASC
      LIMIT 1
    )[SAFE_OFFSET(0)] AS baseline_scr_mg_dl,
    MAX(
      IF(analyte = 'scr' AND labresultoffset <= 2880, valuenum, NULL)
    ) AS peak_scr_48h_mg_dl,
    MAX(
      IF(analyte = 'scr' AND labresultoffset <= 10080, valuenum, NULL)
    ) AS peak_scr_7d_mg_dl,
    MAX(IF(analyte = 'alt', valuenum, NULL)) AS peak_alt_u_l,
    MAX(IF(analyte = 'ast', valuenum, NULL)) AS peak_ast_u_l,
    MAX(IF(analyte = 'alp', valuenum, NULL)) AS peak_alp_u_l,
    MAX(IF(analyte = 'tbl', valuenum, NULL)) AS peak_tbl_mg_dl
  FROM icu_labs
  GROUP BY hospitalid, patientunitstayid
)

SELECT
  stay.hospitalid,
  stay.patientunitstayid,
  stay.gender,
  stay.age_years,
  m.baseline_scr_mg_dl,
  m.peak_scr_48h_mg_dl,
  m.peak_scr_7d_mg_dl,
  m.peak_alt_u_l,
  m.peak_ast_u_l,
  m.peak_alp_u_l,
  m.peak_tbl_mg_dl,
  (
    GREATEST(COALESCE(m.peak_alt_u_l / 40.0, 0), COALESCE(m.peak_ast_u_l / 40.0, 0)) >= 3.0
    AND COALESCE(m.peak_tbl_mg_dl / 1.2, 0) >= 2.0
    AND COALESCE(m.peak_alp_u_l / 120.0, 999.0) < 2.0
  ) AS potential_hys_law_hepatocellular_screen
FROM icu_stays AS stay
INNER JOIN stay_metrics AS m
  ON stay.hospitalid = m.hospitalid
 AND stay.patientunitstayid = m.patientunitstayid;
