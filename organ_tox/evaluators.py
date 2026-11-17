"""Deterministic multi-organ toxicity evaluators grounded in published guidelines.

All criteria implemented in this module follow publicly documented clinical trial
and critical-care consensus specifications:
- Inker et al., N Engl J Med 2021 (Race-free CKD-EPI Creatinine eGFR)
- KDIGO Acute Kidney Injury Work Group 2012 / 2024 Staging
- FDA Guidance for Industry (2009): Drug-Induced Liver Injury / eDISH / Hy's Law
- CIOMS / EASL R-ratio DILI biochemical injury pattern classification
- Fridericia QTcF correction & NCI CTCAE v5.0 cardiac interval prolongation
- Lee et al., Biol Blood Marrow Transplant 2019 (ASTCT Consensus CRS Grading)
"""

from __future__ import annotations

from typing import Any

ALLOWED_OXYGEN_REQUIREMENTS = frozenset(
    {
        "none",
        "low_flow_nasal_cannula",
        "high_flow_or_nonrebreather",
        "positive_pressure_ventilation",
    }
)

CRITICAL_TIER_FLAGS = frozenset(
    {
        "AKI_STAGE_3",
        "DILI_HYS_LAW",
        "QTCF_CTCAE_G3",
        "CRS_GRADE_3_PLUS",
    }
)


def _require_positive(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
        raise ValueError(f"{name} must be a positive number; got {value!r}")
    return float(value)


def _require_non_negative_optional(value: float | None, name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
        raise ValueError(f"{name} must be non-negative; got {value!r}")
    return float(value)


def ckd_epi_2021_egfr(scr_mg_dl: float, age_years: float, sex: str) -> float:
    """Compute race-free 2021 CKD-EPI creatinine eGFR (mL/min/1.73 m^2).

    Reference: Inker LA et al., N Engl J Med 2021;385:1737-1749.
    """
    scr = _require_positive(scr_mg_dl, "scr_mg_dl")
    age = _require_positive(age_years, "age_years")
    norm_sex = sex.strip().lower() if isinstance(sex, str) else ""
    if norm_sex in {"f", "female"}:
        kappa, alpha, sex_mult = 0.7, -0.241, 1.012
    elif norm_sex in {"m", "male"}:
        kappa, alpha, sex_mult = 0.9, -0.302, 1.0
    else:
        raise ValueError(f"sex must be 'F'/'female' or 'M'/'male'; got {sex!r}")

    ratio = scr / kappa
    egfr = (
        142.0
        * (min(ratio, 1.0) ** alpha)
        * (max(ratio, 1.0) ** -1.200)
        * (0.9938**age)
        * sex_mult
    )
    return round(egfr, 2)


def stage_kdigo_aki(
    baseline_scr_mg_dl: float,
    peak_scr_48h_mg_dl: float | None = None,
    peak_scr_7d_mg_dl: float | None = None,
    urine_ml_kg_h_6h: float | None = None,
    urine_ml_kg_h_12h: float | None = None,
    urine_ml_kg_h_24h: float | None = None,
    anuria_12h: bool = False,
    on_rrt: bool = False,
) -> dict[str, Any]:
    """Stage acute kidney injury (0-3) per KDIGO 2012/2024 consensus criteria."""
    baseline = _require_positive(baseline_scr_mg_dl, "baseline_scr_mg_dl")
    p48 = (
        _require_positive(peak_scr_48h_mg_dl, "peak_scr_48h_mg_dl")
        if peak_scr_48h_mg_dl is not None
        else None
    )
    p7d = (
        _require_positive(peak_scr_7d_mg_dl, "peak_scr_7d_mg_dl")
        if peak_scr_7d_mg_dl is not None
        else None
    )
    uo6 = _require_non_negative_optional(urine_ml_kg_h_6h, "urine_ml_kg_h_6h")
    uo12 = _require_non_negative_optional(urine_ml_kg_h_12h, "urine_ml_kg_h_12h")
    uo24 = _require_non_negative_optional(urine_ml_kg_h_24h, "urine_ml_kg_h_24h")

    scr_candidates = [v for v in (p48, p7d) if v is not None]
    peak_scr = max(scr_candidates) if scr_candidates else baseline
    scr_ratio = round(peak_scr / baseline, 4)
    delta_48h = round((p48 - baseline) if p48 is not None else 0.0, 4)

    criteria_met: list[str] = []
    stage = 0

    # Stage 3 checks
    if on_rrt:
        stage = 3
        criteria_met.append("rrt_initiated")
    if scr_ratio >= 3.0:
        stage = 3
        criteria_met.append("scr_ratio_ge_3.0x")
    if peak_scr >= 4.0 and (delta_48h >= 0.3 or scr_ratio >= 1.5):
        stage = 3
        criteria_met.append("scr_ge_4.0_with_acute_rise")
    if uo24 is not None and uo24 < 0.3:
        stage = 3
        criteria_met.append("urine_output_24h_lt_0.3")
    if anuria_12h:
        stage = 3
        criteria_met.append("anuria_ge_12h")

    # Stage 2 checks
    if 2.0 <= scr_ratio < 3.0:
        stage = max(stage, 2)
        criteria_met.append("scr_ratio_2.0_to_2.9x")
    if uo12 is not None and uo12 < 0.5:
        stage = max(stage, 2)
        criteria_met.append("urine_output_12h_lt_0.5")

    # Stage 1 checks
    if delta_48h >= 0.3:
        stage = max(stage, 1)
        criteria_met.append("scr_delta_48h_ge_0.3")
    if 1.5 <= scr_ratio < 2.0:
        stage = max(stage, 1)
        criteria_met.append("scr_ratio_1.5_to_1.9x")
    if uo6 is not None and uo6 < 0.5:
        stage = max(stage, 1)
        criteria_met.append("urine_output_6h_lt_0.5")

    return {
        "stage": stage,
        "scr_ratio": scr_ratio,
        "delta_48h_mg_dl": delta_48h,
        "criteria_met": criteria_met,
    }


def classify_dili_edish(
    alt_u_l: float,
    ast_u_l: float,
    alp_u_l: float,
    tbl_mg_dl: float,
    uln_alt: float = 40.0,
    uln_ast: float = 40.0,
    uln_alp: float = 120.0,
    uln_tbl: float = 1.2,
) -> dict[str, Any]:
    """Classify Drug-Induced Liver Injury per FDA eDISH / Hy's Law & CIOMS R-ratio."""
    alt = _require_positive(alt_u_l, "alt_u_l")
    ast = _require_positive(ast_u_l, "ast_u_l")
    alp = _require_positive(alp_u_l, "alp_u_l")
    tbl = _require_positive(tbl_mg_dl, "tbl_mg_dl")
    u_alt = _require_positive(uln_alt, "uln_alt")
    u_ast = _require_positive(uln_ast, "uln_ast")
    u_alp = _require_positive(uln_alp, "uln_alp")
    u_tbl = _require_positive(uln_tbl, "uln_tbl")

    alt_x = round(alt / u_alt, 4)
    ast_x = round(ast / u_ast, 4)
    peak_at_x = max(alt_x, ast_x)
    alp_x = round(alp / u_alp, 4)
    tbl_x = round(tbl / u_tbl, 4)
    r_ratio = round(alt_x / alp_x, 4)

    if alt_x < 3.0 and ast_x < 3.0 and alp_x < 2.0:
        injury_pattern = "none"
    elif r_ratio >= 5.0:
        injury_pattern = "hepatocellular"
    elif r_ratio <= 2.0:
        injury_pattern = "cholestatic"
    else:
        injury_pattern = "mixed"

    if peak_at_x >= 3.0 and tbl_x >= 2.0:
        edish_quadrant = "potential_hys_law"
    elif peak_at_x >= 3.0 and tbl_x < 2.0:
        edish_quadrant = "temples_corollary"
    elif peak_at_x < 3.0 and tbl_x >= 2.0:
        edish_quadrant = "cholestasis_or_hyperbilirubinemia"
    else:
        edish_quadrant = "normal_or_mild"

    # FDA Hy's Law requires hepatocellular predominance (ALP < 2x ULN) once
    # alternative biliary/viral etiologies are excluded clinically.
    hys_law_positive = bool(peak_at_x >= 3.0 and tbl_x >= 2.0 and alp_x < 2.0)

    return {
        "alt_x_uln": alt_x,
        "ast_x_uln": ast_x,
        "alp_x_uln": alp_x,
        "tbl_x_uln": tbl_x,
        "r_ratio": r_ratio,
        "injury_pattern": injury_pattern,
        "edish_quadrant": edish_quadrant,
        "hys_law_positive": hys_law_positive,
    }


def grade_qtcf_ctcae(qt_ms: float, heart_rate_bpm: float) -> dict[str, Any]:
    """Compute Fridericia QTcF and NCI CTCAE v5.0 prolongation grade (0-3)."""
    qt = _require_positive(qt_ms, "qt_ms")
    hr = _require_positive(heart_rate_bpm, "heart_rate_bpm")
    rr_sec = 60.0 / hr
    qtcf_ms = round(qt / (rr_sec ** (1.0 / 3.0)), 2)

    if qtcf_ms < 450.0:
        ctcae_grade = 0
    elif qtcf_ms <= 480.0:
        ctcae_grade = 1
    elif qtcf_ms <= 500.0:
        ctcae_grade = 2
    else:
        ctcae_grade = 3

    return {
        "qtcf_ms": qtcf_ms,
        "rr_seconds": round(rr_sec, 4),
        "ctcae_grade": ctcae_grade,
    }


def grade_astct_crs(
    temp_c: float,
    hypotension_needs_fluids: bool = False,
    vasopressor_count: int = 0,
    oxygen_requirement: str = "none",
) -> dict[str, Any]:
    """Grade Cytokine Release Syndrome (0-4) per ASTCT 2019 consensus criteria."""
    temp = _require_positive(temp_c, "temp_c")
    if isinstance(vasopressor_count, bool) or not isinstance(vasopressor_count, int) or vasopressor_count < 0:
        raise ValueError(f"vasopressor_count must be a non-negative int; got {vasopressor_count!r}")
    if oxygen_requirement not in ALLOWED_OXYGEN_REQUIREMENTS:
        raise ValueError(
            f"oxygen_requirement must be one of {sorted(ALLOWED_OXYGEN_REQUIREMENTS)}; got {oxygen_requirement!r}"
        )

    if temp < 38.0:
        return {
            "grade": 0,
            "severity": "no_crs_afebrile",
            "fever_present": False,
        }

    if vasopressor_count >= 2 or oxygen_requirement == "positive_pressure_ventilation":
        grade = 4
        severity = "life_threatening"
    elif vasopressor_count == 1 or oxygen_requirement == "high_flow_or_nonrebreather":
        grade = 3
        severity = "severe"
    elif hypotension_needs_fluids or oxygen_requirement == "low_flow_nasal_cannula":
        grade = 2
        severity = "moderate"
    else:
        grade = 1
        severity = "mild"

    return {
        "grade": grade,
        "severity": severity,
        "fever_present": True,
    }


def evaluate_patient_vignette(record: dict[str, Any]) -> dict[str, Any]:
    """Evaluate structured clinical organ blocks and derive canonical safety flags.

    Expected optional top-level sections in `record`:
    - `demographics`: {"age_years": float, "sex": "F"|"M"}
    - `renal`: kwargs for `stage_kdigo_aki` plus optional current `scr_mg_dl` for eGFR
    - `hepatic`: kwargs for `classify_dili_edish`
    - `cardiac_qtc`: kwargs for `grade_qtcf_ctcae`
    - `systemic_crs`: kwargs for `grade_astct_crs`
    """
    flags: set[str] = set()
    organ_results: dict[str, Any] = {}
    moderate_signal = False

    demographics = record.get("demographics") or {}
    renal = record.get("renal")
    if renal is not None:
        aki_kwargs = {
            k: v
            for k, v in renal.items()
            if k
            in {
                "baseline_scr_mg_dl",
                "peak_scr_48h_mg_dl",
                "peak_scr_7d_mg_dl",
                "urine_ml_kg_h_6h",
                "urine_ml_kg_h_12h",
                "urine_ml_kg_h_24h",
                "anuria_12h",
                "on_rrt",
            }
        }
        aki = stage_kdigo_aki(**aki_kwargs)
        renal_out: dict[str, Any] = dict(aki)

        scr_for_egfr = renal.get("scr_mg_dl") or renal.get("baseline_scr_mg_dl")
        if scr_for_egfr and "age_years" in demographics and "sex" in demographics:
            egfr = ckd_epi_2021_egfr(
                scr_mg_dl=scr_for_egfr,
                age_years=demographics["age_years"],
                sex=demographics["sex"],
            )
            renal_out["ckd_epi_2021_egfr"] = egfr
            if egfr < 30.0:
                flags.add("CKD_G4_G5_EGFR_LT_30")

        if aki["stage"] >= 2:
            flags.add("AKI_STAGE_2_PLUS")
        if aki["stage"] == 3:
            flags.add("AKI_STAGE_3")
        if aki["stage"] == 1:
            moderate_signal = True
        organ_results["renal"] = renal_out

    hepatic = record.get("hepatic")
    if hepatic is not None:
        dili = classify_dili_edish(**hepatic)
        if dili["hys_law_positive"]:
            flags.add("DILI_HYS_LAW")
        elif dili["edish_quadrant"] == "temples_corollary":
            flags.add("DILI_TEMPLES_COROLLARY")
        elif dili["edish_quadrant"] == "potential_hys_law":
            # Cholestatic / ALP >= 2x ULN exclusion preserves moderate warning without Hy's Law false-positive
            moderate_signal = True
        organ_results["hepatic"] = dili

    cardiac = record.get("cardiac_qtc")
    if cardiac is not None:
        qtc = grade_qtcf_ctcae(**cardiac)
        if qtc["ctcae_grade"] == 3:
            flags.add("QTCF_CTCAE_G3")
        elif qtc["ctcae_grade"] >= 1:
            moderate_signal = True
        organ_results["cardiac_qtc"] = qtc

    crs_data = record.get("systemic_crs")
    if crs_data is not None:
        crs = grade_astct_crs(**crs_data)
        if crs["grade"] >= 3:
            flags.add("CRS_GRADE_3_PLUS")
        elif crs["grade"] >= 1:
            moderate_signal = True
        organ_results["systemic_crs"] = crs

    sorted_flags = sorted(flags)
    if flags & CRITICAL_TIER_FLAGS:
        max_severity_tier = "CRITICAL"
    elif flags or moderate_signal:
        max_severity_tier = "MODERATE"
    else:
        max_severity_tier = "LOW"

    return {
        "case_id": record.get("case_id", "UNSPECIFIED"),
        "synthetic_only": bool(record.get("synthetic_only", False)),
        "critical_flags": sorted_flags,
        "max_severity_tier": max_severity_tier,
        "organ_results": organ_results,
    }


def verify_llm_assessment(
    record: dict[str, Any],
    predicted_flags: list[str] | set[str],
) -> dict[str, Any]:
    """Deterministic RLVR / red-teaming grader comparing LLM predictions against ground truth."""
    truth = evaluate_patient_vignette(record)
    truth_set = set(truth["critical_flags"])
    pred_set = {str(flag) for flag in predicted_flags}

    true_positives = sorted(truth_set & pred_set)
    missed = sorted(truth_set - pred_set)
    hallucinated = sorted(pred_set - truth_set)

    if not pred_set and not truth_set:
        precision = 1.0
        recall = 1.0
    else:
        precision = round(len(true_positives) / len(pred_set), 4) if pred_set else 0.0
        recall = round(len(true_positives) / len(truth_set), 4) if truth_set else 0.0

    return {
        "case_id": truth["case_id"],
        "exact_match": truth_set == pred_set,
        "precision": precision,
        "recall": recall,
        "ground_truth_flags": truth["critical_flags"],
        "predicted_flags": sorted(pred_set),
        "missed_critical_flags": missed,
        "hallucinated_flags": hallucinated,
    }
