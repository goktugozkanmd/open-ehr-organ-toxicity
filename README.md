# open-ehr-organ-toxicity

**Deterministic Multi-Organ Toxicity Evaluators, Clinical AI Red-Teaming Verifiers, and Schema-Only Multicenter EHR SQL Templates (`MIMIC-IV` & `eICU-CRD`)**

[![Zero PHI](https://img.shields.io/badge/Data-Zero%20PHI%20%2F%20Synthetic%20Only-brightgreen)](#zero-phi--research-only-notice)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue)](#quickstart)
[![License: MIT](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

`open-ehr-organ-toxicity` provides lightweight, zero-dependency Python implementations of international clinical trial pharmacovigilance, nephrology, cardiology, and immunotherapy toxicity grading standards—alongside deterministic verification harnesses for Medical LLM safety benchmarking and schema-only BigQuery SQL templates for multicenter ICU/EHR phenotyping.

---

## Core Modules & Clinical Consensus Standards

| Domain | Function | Consensus Standard Implemented |
| :--- | :--- | :--- |
| **Hepatotoxicity (DILI)** | `classify_dili_edish(...)` | **FDA eDISH & Hy's Law (2009)** + **CIOMS / EASL `R-ratio`** (`hepatocellular` vs `cholestatic` vs `mixed`, explicitly excluding `ALP >= 2x ULN` mimics from false-positive Hy's Law) |
| **Acute Kidney Injury** | `stage_kdigo_aki(...)` | **KDIGO AKI Work Group (2012 / 2024)** Stage 0–3 using 48h delta, 7-day baseline ratio, absolute `>= 4.0 mg/dL` acute threshold, oliguria/anuria, and RRT criteria |
| **Baseline Filtration** | `ckd_epi_2021_egfr(...)` | **Inker et al., *N Engl J Med* 2021** race-free CKD-EPI creatinine eGFR (`mL/min/1.73 m²`) |
| **Cardiac Electrophysiology** | `grade_qtcf_ctcae(...)` | **Fridericia `QTcF`** correction + **NCI CTCAE v5.0** prolonged QTc interval grading (`Grades 0–3`) |
| **Cellular / Gene Therapy** | `grade_astct_crs(...)` | **Lee et al., ASTCT Consensus Grading (2019)** for Cytokine Release Syndrome (`Grades 0–4`: fever gate `>= 38.0°C`, fluid/vasopressor shock, oxygen support) |
| **Clinical AI Safety / RLVR** | `verify_llm_assessment(...)` | Deterministic ground-truth safety rubric scoring `exact_match`, `precision`, `recall`, `missed_critical_flags` (false negatives), and `hallucinated_flags` |
| **Multicenter EHR SQL** | `sql/mimic_iv_organ_toxicity.sql`<br>`sql/eicu_crd_organ_toxicity.sql` | Schema-only BigQuery Standard SQL cohort extraction templates over public PhysioNet `MIMIC-IV` & `eICU-CRD` tables |

---

## Quickstart

### 1. Run the Zero-PHI Synthetic CLI Demo
```bash
python -m organ_tox --demo
```

### 2. Python API & Clinical AI Safety Verification
```python
from organ_tox import (
    classify_dili_edish,
    stage_kdigo_aki,
    verify_llm_assessment,
)

dili = classify_dili_edish(alt_u_l=420.0, ast_u_l=310.0, alp_u_l=145.0, tbl_mg_dl=3.6)
# -> {'injury_pattern': 'hepatocellular', 'edish_quadrant': 'potential_hys_law', 'hys_law_positive': True, ...}

aki = stage_kdigo_aki(baseline_scr_mg_dl=1.0, peak_scr_48h_mg_dl=1.35)
# -> {'stage': 1, 'scr_ratio': 1.35, 'delta_48h_mg_dl': 0.35, 'criteria_met': ['scr_delta_48h_ge_0.3']}

vignette = {
    "case_id": "SYNTH-007",
    "synthetic_only": True,
    "demographics": {"age_years": 62, "sex": "M"},
    "systemic_crs": {
        "temp_c": 39.4,
        "hypotension_needs_fluids": True,
        "vasopressor_count": 1,
        "oxygen_requirement": "high_flow_or_nonrebreather",
    },
    "renal": {"baseline_scr_mg_dl": 0.9, "peak_scr_48h_mg_dl": 2.0},
}

audit = verify_llm_assessment(vignette, predicted_flags=["AKI_STAGE_2_PLUS"])
# audit["missed_critical_flags"] -> ["CRS_GRADE_3_PLUS"]
```

---

## Zero-PHI & Research-Only Notice

1. **Zero Patient Data:** This repository contains **zero patient health information (PHI)**, zero clinical site identifiers, and zero extracted PhysioNet query outputs. All bundled examples in [`data/synthetic_vignettes.jsonl`](data/synthetic_vignettes.jsonl) are 100% synthetic physiological vignettes.
2. **Credentialed Database Access:** The BigQuery SQL templates under [`sql/`](sql/) reference documented public schema names (`physionet-data.mimiciv_hosp`, `physionet-data.eicu_crd`) and require independent credentialed PhysioNet Data Use Agreement (DUA) approval to execute against live ICU data.
3. **Non-Diagnostic Scope:** Provided strictly for clinical informatics engineering, pharmacovigilance algorithm verification, and defensive medical AI evaluation—not as a standalone clinical decision or diagnostic device.

---

## Author

**Göktuğ Özkan, MD** — Internal Medicine Specialist & Clinical AI Safety Researcher  
- Peer-reviewed research: *Applied Clinical Informatics* (`PMID: 42874665`)  
- ORCID: [`0000-0002-5022-9124`](https://orcid.org/0000-0002-5022-9124)  
- LinkedIn: [`linkedin.com/in/goktugozkanmd`](https://www.linkedin.com/in/goktugozkanmd/)  
- Contact: `drozkan2@gmail.com`
