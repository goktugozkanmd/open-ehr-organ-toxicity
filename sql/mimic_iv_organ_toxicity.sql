-- =============================================================================
-- MIMIC-IV (BigQuery Standard SQL) Schema-Only Organ Toxicity Cohort Template
-- Author: Goktug Ozkan, MD (ORCID: 0000-0002-5022-9124)
-- Scope: Public PhysioNet BigQuery schema template (`physionet-data.mimiciv_hosp`).
--        Contains ZERO patient data or exported query results. Credentialed
--        PhysioNet Data Use Agreement (DUA) is required to run against live tables.
-- =============================================================================

WITH target_labitems AS (
  SELECT
    itemid,
    CASE
      WHEN itemid = 50912 THEN 'scr'   -- Creatinine (mg/dL)
      WHEN itemid = 50861 THEN 'alt'   -- Alanine Aminotransferase (IU/L)
      WHEN itemid = 50878 THEN 'ast'   -- Aspartate Aminotransferase (IU/L)
      WHEN itemid = 50863 THEN 'alp'   -- Alkaline Phosphatase (IU/L)
      WHEN itemid = 50885 THEN 'tbl'   -- Bilirubin, Total (mg/dL)
    END AS analyte
  FROM `physionet-data.mimiciv_hosp.d_labitems`
  WHERE itemid IN (50912, 50861, 50878, 50863, 50885)
),

adult_admissions AS (
  SELECT
    adm.subject_id,
    adm.hadm_id,
    adm.admittime,
    adm.dischtime,
    pat.gender,
    pat.anchor_age
  FROM `physionet-data.mimiciv_hosp.admissions` AS adm
  INNER JOIN `physionet-data.mimiciv_hosp.patients` AS pat
    ON adm.subject_id = pat.subject_id
  WHERE pat.anchor_age >= 18
),

admission_labs AS (
  SELECT
    adm.subject_id,
    adm.hadm_id,
    adm.admittime,
    tl.analyte,
    le.charttime,
    le.valuenum
  FROM adult_admissions AS adm
  INNER JOIN `physionet-data.mimiciv_hosp.labevents` AS le
    ON adm.hadm_id = le.hadm_id
  INNER JOIN target_labitems AS tl
    ON le.itemid = tl.itemid
  WHERE le.valuenum IS NOT NULL
    AND le.valuenum > 0
    AND le.charttime >= adm.admittime
    AND le.charttime <= adm.dischtime
),

renal_kdigo_summary AS (
  SELECT
    subject_id,
    hadm_id,
    ARRAY_AGG(
      IF(analyte = 'scr', valuenum, NULL) IGNORE NULLS
      ORDER BY charttime ASC
      LIMIT 1
    )[SAFE_OFFSET(0)] AS baseline_scr_mg_dl,
    MAX(
      IF(
        analyte = 'scr'
        AND charttime <= TIMESTAMP_ADD(admittime, INTERVAL 48 HOUR),
        valuenum,
        NULL
      )
    ) AS peak_scr_48h_mg_dl,
    MAX(
      IF(
        analyte = 'scr'
        AND charttime <= TIMESTAMP_ADD(admittime, INTERVAL 7 DAY),
        valuenum,
        NULL
      )
    ) AS peak_scr_7d_mg_dl
  FROM admission_labs
  GROUP BY subject_id, hadm_id
),

hepatic_edish_summary AS (
  SELECT
    subject_id,
    hadm_id,
    MAX(IF(analyte = 'alt', valuenum, NULL)) AS peak_alt_u_l,
    MAX(IF(analyte = 'ast', valuenum, NULL)) AS peak_ast_u_l,
    MAX(IF(analyte = 'alp', valuenum, NULL)) AS peak_alp_u_l,
    MAX(IF(analyte = 'tbl', valuenum, NULL)) AS peak_tbl_mg_dl
  FROM admission_labs
  GROUP BY subject_id, hadm_id
)

SELECT
  adm.subject_id,
  adm.hadm_id,
  adm.gender,
  adm.anchor_age AS age_years,
  ren.baseline_scr_mg_dl,
  ren.peak_scr_48h_mg_dl,
  ren.peak_scr_7d_mg_dl,
  CASE
    WHEN ren.baseline_scr_mg_dl IS NULL OR ren.baseline_scr_mg_dl <= 0 THEN NULL
    WHEN SAFE_DIVIDE(COALESCE(ren.peak_scr_7d_mg_dl, ren.peak_scr_48h_mg_dl), ren.baseline_scr_mg_dl) >= 3.0
      OR (
        COALESCE(ren.peak_scr_7d_mg_dl, ren.peak_scr_48h_mg_dl) >= 4.0
        AND (
          COALESCE(ren.peak_scr_48h_mg_dl, 0) - ren.baseline_scr_mg_dl >= 0.3
          OR SAFE_DIVIDE(COALESCE(ren.peak_scr_7d_mg_dl, ren.peak_scr_48h_mg_dl), ren.baseline_scr_mg_dl) >= 1.5
        )
      ) THEN 3
    WHEN SAFE_DIVIDE(COALESCE(ren.peak_scr_7d_mg_dl, ren.peak_scr_48h_mg_dl), ren.baseline_scr_mg_dl) >= 2.0 THEN 2
    WHEN SAFE_DIVIDE(COALESCE(ren.peak_scr_7d_mg_dl, ren.peak_scr_48h_mg_dl), ren.baseline_scr_mg_dl) >= 1.5
      OR (COALESCE(ren.peak_scr_48h_mg_dl, 0) - ren.baseline_scr_mg_dl >= 0.3) THEN 1
    ELSE 0
  END AS kdigo_aki_stage_lab,
  hep.peak_alt_u_l,
  hep.peak_ast_u_l,
  hep.peak_alp_u_l,
  hep.peak_tbl_mg_dl,
  (
    GREATEST(COALESCE(hep.peak_alt_u_l / 40.0, 0), COALESCE(hep.peak_ast_u_l / 40.0, 0)) >= 3.0
    AND COALESCE(hep.peak_tbl_mg_dl / 1.2, 0) >= 2.0
    AND COALESCE(hep.peak_alp_u_l / 120.0, 999.0) < 2.0
  ) AS potential_hys_law_hepatocellular_screen
FROM adult_admissions AS adm
LEFT JOIN renal_kdigo_summary AS ren
  ON adm.subject_id = ren.subject_id AND adm.hadm_id = ren.hadm_id
LEFT JOIN hepatic_edish_summary AS hep
  ON adm.subject_id = hep.subject_id AND adm.hadm_id = hep.hadm_id;
