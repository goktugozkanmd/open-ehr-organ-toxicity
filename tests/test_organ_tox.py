"""Unit, CLI, and privacy-guardrail tests for open-ehr-organ-toxicity."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from organ_tox import (
    ckd_epi_2021_egfr,
    classify_dili_edish,
    evaluate_patient_vignette,
    grade_astct_crs,
    grade_qtcf_ctcae,
    stage_kdigo_aki,
    verify_llm_assessment,
)
from organ_tox.__main__ import DEFAULT_DEMO_PATH, load_jsonl, main, run_batch

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_ckd_epi_2021_egfr_female_and_male_reference() -> None:
    # Female 55y, SCr 1.0 mg/dL: kappa=0.7, (1.0/0.7)^-1.200 * 0.9938^55 * 1.012 * 142
    egfr_f = ckd_epi_2021_egfr(scr_mg_dl=1.0, age_years=55, sex="F")
    assert 66.0 <= egfr_f <= 68.5

    # Male 68y, SCr 4.2 mg/dL: severe CKD/AKI-range filtration (< 20)
    egfr_m = ckd_epi_2021_egfr(scr_mg_dl=4.2, age_years=68, sex="male")
    assert 13.0 <= egfr_m <= 16.5

    with pytest.raises(ValueError):
        ckd_epi_2021_egfr(scr_mg_dl=0.0, age_years=50, sex="M")
    with pytest.raises(ValueError):
        ckd_epi_2021_egfr(scr_mg_dl=1.0, age_years=50, sex="unknown")


def test_stage_kdigo_aki_stages_0_to_3() -> None:
    s0 = stage_kdigo_aki(baseline_scr_mg_dl=0.9, peak_scr_48h_mg_dl=1.0)
    assert s0["stage"] == 0
    assert s0["criteria_met"] == []

    # Stage 1 via >= 0.3 mg/dL delta within 48h despite <1.5x ratio
    s1 = stage_kdigo_aki(baseline_scr_mg_dl=1.0, peak_scr_48h_mg_dl=1.35)
    assert s1["stage"] == 1
    assert "scr_delta_48h_ge_0.3" in s1["criteria_met"]

    # Stage 2 via 2.0-2.9x baseline over 7 days
    s2 = stage_kdigo_aki(baseline_scr_mg_dl=1.0, peak_scr_7d_mg_dl=2.2)
    assert s2["stage"] == 2
    assert "scr_ratio_2.0_to_2.9x" in s2["criteria_met"]

    # Stage 3 via RRT and anuria even with moderate SCr rise
    s3_rrt = stage_kdigo_aki(
        baseline_scr_mg_dl=1.1, peak_scr_48h_mg_dl=1.8, anuria_12h=True, on_rrt=True
    )
    assert s3_rrt["stage"] == 3
    assert "rrt_initiated" in s3_rrt["criteria_met"]
    assert "anuria_ge_12h" in s3_rrt["criteria_met"]


def test_classify_dili_edish_hys_law_vs_cholestatic_exclusion() -> None:
    # True FDA Hy's Law hepatocellular injury: ALT/AST >= 3x, TBL >= 2x, ALP < 2x ULN
    hys = classify_dili_edish(alt_u_l=420.0, ast_u_l=310.0, alp_u_l=145.0, tbl_mg_dl=3.6)
    assert hys["injury_pattern"] == "hepatocellular"
    assert hys["edish_quadrant"] == "potential_hys_law"
    assert hys["hys_law_positive"] is True

    # Prominent cholestasis (ALP = 4x ULN): falls in potential quadrant plot, NOT Hy's Law
    chol = classify_dili_edish(alt_u_l=140.0, ast_u_l=128.0, alp_u_l=480.0, tbl_mg_dl=4.2)
    assert chol["injury_pattern"] == "cholestatic"
    assert chol["edish_quadrant"] == "potential_hys_law"
    assert chol["hys_law_positive"] is False

    # Temple's Corollary: isolated AT >= 3x without bilirubin >= 2x
    temple = classify_dili_edish(alt_u_l=240.0, ast_u_l=180.0, alp_u_l=110.0, tbl_mg_dl=1.0)
    assert temple["edish_quadrant"] == "temples_corollary"
    assert temple["hys_law_positive"] is False

    # Mixed R-ratio (2 < R < 5): ALT=4x, ALP=1.5x -> R=2.667
    mixed = classify_dili_edish(alt_u_l=160.0, ast_u_l=140.0, alp_u_l=180.0, tbl_mg_dl=1.0)
    assert mixed["injury_pattern"] == "mixed"


def test_grade_qtcf_ctcae_boundaries() -> None:
    g0 = grade_qtcf_ctcae(qt_ms=420.0, heart_rate_bpm=60.0)
    assert g0["qtcf_ms"] == 420.0
    assert g0["ctcae_grade"] == 0

    g1 = grade_qtcf_ctcae(qt_ms=465.0, heart_rate_bpm=60.0)
    assert g1["ctcae_grade"] == 1

    g2 = grade_qtcf_ctcae(qt_ms=495.0, heart_rate_bpm=60.0)
    assert g2["ctcae_grade"] == 2

    g3 = grade_qtcf_ctcae(qt_ms=515.0, heart_rate_bpm=60.0)
    assert g3["qtcf_ms"] == 515.0
    assert g3["ctcae_grade"] == 3


def test_grade_astct_crs_consensus_grades() -> None:
    assert grade_astct_crs(temp_c=37.4, hypotension_needs_fluids=True)["grade"] == 0
    assert grade_astct_crs(temp_c=38.5)["grade"] == 1
    assert (
        grade_astct_crs(temp_c=38.9, oxygen_requirement="low_flow_nasal_cannula")[
            "grade"
        ]
        == 2
    )
    assert (
        grade_astct_crs(
            temp_c=39.4,
            hypotension_needs_fluids=True,
            vasopressor_count=1,
            oxygen_requirement="high_flow_or_nonrebreather",
        )["grade"]
        == 3
    )
    assert (
        grade_astct_crs(
            temp_c=39.8,
            vasopressor_count=2,
            oxygen_requirement="positive_pressure_ventilation",
        )["grade"]
        == 4
    )

    with pytest.raises(ValueError):
        grade_astct_crs(temp_c=39.0, oxygen_requirement="invalid_mode")


def test_verify_llm_assessment_detects_missed_and_hallucinated_flags() -> None:
    vignettes = {row["case_id"]: row for row in load_jsonl(DEFAULT_DEMO_PATH)}

    # Model misses CRS_GRADE_3_PLUS in dual-toxicity SYNTH-007 and hallucinates DILI_HYS_LAW
    audit = verify_llm_assessment(
        vignettes["SYNTH-007"],
        predicted_flags=["AKI_STAGE_2_PLUS", "DILI_HYS_LAW"],
    )
    assert audit["exact_match"] is False
    assert audit["missed_critical_flags"] == ["CRS_GRADE_3_PLUS"]
    assert audit["hallucinated_flags"] == ["DILI_HYS_LAW"]
    assert audit["precision"] == 0.5
    assert audit["recall"] == 0.5

    # Exact match on clean Hy's Law case SYNTH-001
    perfect = verify_llm_assessment(vignettes["SYNTH-001"], predicted_flags=["DILI_HYS_LAW"])
    assert perfect["exact_match"] is True
    assert perfect["precision"] == 1.0
    assert perfect["recall"] == 1.0


def test_cli_demo_execution_and_vignette_counts(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = main(["--demo"])
    assert exit_code == 0
    captured = json.loads(capsys.readouterr().out)
    assert captured["vignettes_evaluated"] == 8
    assert captured["critical_cases"] == 4
    assert captured["moderate_cases"] == 3

    by_id = {item["case_id"]: item for item in captured["per_case_results"]}
    assert by_id["SYNTH-001"]["critical_flags"] == ["DILI_HYS_LAW"]
    assert by_id["SYNTH-002"]["critical_flags"] == []
    assert by_id["SYNTH-002"]["max_severity_tier"] == "MODERATE"
    assert by_id["SYNTH-004"]["critical_flags"] == [
        "AKI_STAGE_2_PLUS",
        "AKI_STAGE_3",
        "CKD_G4_G5_EGFR_LT_30",
    ]
    assert by_id["SYNTH-006"]["critical_flags"] == ["QTCF_CTCAE_G3"]
    assert by_id["SYNTH-008"]["max_severity_tier"] == "LOW"


def test_ckd_epi_2021_egfr_low_scr_branch_and_validation() -> None:
    # Exercises min(ratio, 1.0) ** alpha branch (SCr <= kappa) for female and male
    egfr_f_low = ckd_epi_2021_egfr(scr_mg_dl=0.6, age_years=35, sex="female")
    egfr_m_low = ckd_epi_2021_egfr(scr_mg_dl=0.7, age_years=35, sex="M")
    assert egfr_f_low > 105.0
    assert egfr_m_low > 105.0

    with pytest.raises(ValueError):
        ckd_epi_2021_egfr(scr_mg_dl=0.9, age_years=-5, sex="M")


def test_stage_kdigo_aki_urine_output_and_absolute_scr_ge_4() -> None:
    s1_uo = stage_kdigo_aki(baseline_scr_mg_dl=1.0, urine_ml_kg_h_6h=0.4)
    assert s1_uo["stage"] == 1
    assert "urine_output_6h_lt_0.5" in s1_uo["criteria_met"]

    s2_uo = stage_kdigo_aki(baseline_scr_mg_dl=1.0, urine_ml_kg_h_12h=0.35)
    assert s2_uo["stage"] == 2
    assert "urine_output_12h_lt_0.5" in s2_uo["criteria_met"]

    s3_uo = stage_kdigo_aki(baseline_scr_mg_dl=1.0, urine_ml_kg_h_24h=0.2)
    assert s3_uo["stage"] == 3
    assert "urine_output_24h_lt_0.3" in s3_uo["criteria_met"]

    # Pre-existing elevated baseline 3.6 -> 4.1 (+0.5 delta in 48h, ratio 1.139x) still hits Stage 3
    s3_abs = stage_kdigo_aki(baseline_scr_mg_dl=3.6, peak_scr_48h_mg_dl=4.1)
    assert s3_abs["stage"] == 3
    assert "scr_ge_4.0_with_acute_rise" in s3_abs["criteria_met"]

    with pytest.raises(ValueError):
        stage_kdigo_aki(baseline_scr_mg_dl=1.0, urine_ml_kg_h_6h=-0.1)


def test_classify_dili_edish_normal_hyperbili_and_validation() -> None:
    normal = classify_dili_edish(alt_u_l=24.0, ast_u_l=21.0, alp_u_l=78.0, tbl_mg_dl=0.7)
    assert normal["injury_pattern"] == "none"
    assert normal["edish_quadrant"] == "normal_or_mild"
    assert normal["hys_law_positive"] is False

    hyperbili = classify_dili_edish(alt_u_l=60.0, ast_u_l=55.0, alp_u_l=100.0, tbl_mg_dl=3.0)
    assert hyperbili["edish_quadrant"] == "cholestasis_or_hyperbilirubinemia"
    assert hyperbili["hys_law_positive"] is False

    with pytest.raises(ValueError):
        classify_dili_edish(alt_u_l=40.0, ast_u_l=40.0, alp_u_l=100.0, tbl_mg_dl=1.0, uln_alt=0.0)


def test_verify_llm_empty_baseline_and_cli_missing_file(tmp_path: Path) -> None:
    vignettes = {row["case_id"]: row for row in load_jsonl(DEFAULT_DEMO_PATH)}
    clean = verify_llm_assessment(vignettes["SYNTH-008"], predicted_flags=[])
    assert clean["exact_match"] is True
    assert clean["precision"] == 1.0
    assert clean["recall"] == 1.0

    missing_exit = main(["--input", str(tmp_path / "does_not_exist.jsonl")])
    assert missing_exit == 1


def test_privacy_and_unpublished_manuscript_zero_leak_guardrail() -> None:
    forbidden_patterns = [
        re.compile(r"\b0?538\b"),
        re.compile(r"\b98\s*30\b"),
        re.compile(r"\+90"),
        re.compile(r"LABMED-"),
        re.compile(r"BCF-26"),
        re.compile(r"TIH-26"),
        re.compile(r"JAPHA-D"),
        re.compile(r"AJHSP-"),
    ]
    ignored_dirs = {".git", "__pycache__", ".pytest_cache"}
    for path in REPO_ROOT.rglob("*"):
        if any(part in ignored_dirs for part in path.parts):
            continue
        if path.is_file() and path.name != "test_organ_tox.py":
            text = path.read_text(encoding="utf-8", errors="ignore")
            for pattern in forbidden_patterns:
                assert not pattern.search(text), f"Forbidden leak {pattern.pattern} in {path}"
