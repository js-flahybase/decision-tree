"""
Condition Evaluator Engine
==========================

Reads patient lab values (from a CSV) and evaluates a set of clinical
"Significant Pattern" / "Early Pattern" rules, one function per
condition, following the template:

    def evaluate_<condition>(labs, patient):
        # 2. Pull values from labs / patient
        some_value = labs.get("some_value")
        # 3. Check if values are abnormal
        some_flag = is_elevated(some_value, SOME_THRESHOLD)
        # 4. Decide category
        if <strong combination met>:
            category = "Significant Pattern"
        elif <lighter combination met>:
            category = "Early Pattern"
        else:
            category = ""
        return [{
            "Condition": "<Condition Name>",
            "Category": category
        }]

OUTPUT FILTERING
-----------------
Only conditions whose Category is "Significant Pattern" or
"Early Pattern" are written to the output CSV. Conditions that
evaluate to "" (no pattern met) are silently dropped from the output
so the report only lists conditions actually flagged for the patient.

CSV INPUT FORMAT
-----------------
The input CSV is expected to have one row per patient, with lab test
names as columns (snake_case, see LAB_KEYS below) plus optional patient
columns: patient_id, sex ("M"/"F"), age.

Example columns:
    patient_id, sex, age, fasting_glucose, hba1c, tsh, free_t4, ...

Run:
    python condition_evaluators.py input.csv output.csv
"""

import csv
import sys
from typing import Optional, Dict, Any, List

# ---------------------------------------------------------------------------
# 1. Generic helpers & Domain mapping
# ---------------------------------------------------------------------------

CONDITION_DOMAINS = {
    "Type 2 Diabetes": "Endocrine Health",
    "Hypothyroidism": "Endocrine Health",
    "Maturity-Onset Diabetes of the Young": "Endocrine Health",
    "Hyperthyroidism": "Endocrine Health",
    "Coronary Artery Disease": "Cardiac Health",
    "Asthma": "Respiratory Health",
    "Rhinitis": "Respiratory Health",
    "NAFLD": "Digestive Health",
    "Inflammatory Bowel Disease": "Digestive Health",
    "Hereditary Hemochromatosis": "Digestive Health",
    "Atopic Dermatitis/Eczema": "Skin Health",
}


def _map_category_message(category: str) -> str:
    if category == "Significant Pattern":
        return "To discuss with General Practitioner"
    elif category in ("Early Pattern", "Elevated susceptibility"):
        return "Worth acting on for prevention"
    else:
        return "Typical - nothing to act on"
    
    
def _to_float(value) -> Optional[float]:
    """Safely convert a value to float, returning None if not possible/blank."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    value = str(value).strip()
    if value == "" or value.lower() in ("na", "n/a", "none", "null"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def is_elevated(value, threshold) -> bool:
    """True if value >= threshold (High/Abnormal-high check)."""
    v = _to_float(value)
    t = _to_float(threshold)
    if v is None or t is None:
        return False
    return v >= t


def is_above(value, threshold) -> bool:
    """True if value > threshold (strict)."""
    v = _to_float(value)
    t = _to_float(threshold)
    if v is None or t is None:
        return False
    return v > t


def is_low(value, threshold) -> bool:
    """True if value <= threshold (Low check)."""
    v = _to_float(value)
    t = _to_float(threshold)
    if v is None or t is None:
        return False
    return v <= t


def is_below(value, threshold) -> bool:
    """True if value < threshold (strict low)."""
    v = _to_float(value)
    t = _to_float(threshold)
    if v is None or t is None:
        return False
    return v < t


def in_range(value, low, high) -> bool:
    """True if low <= value <= high."""
    v = _to_float(value)
    lo = _to_float(low)
    hi = _to_float(high)
    if v is None or lo is None or hi is None:
        return False
    return lo <= v <= hi


def is_male(patient: Dict[str, Any]) -> bool:
    sex = str(patient.get("sex", "")).strip().lower()
    return sex in ("m", "male")


def is_female(patient: Dict[str, Any]) -> bool:
    sex = str(patient.get("sex", "")).strip().lower()
    return sex in ("f", "female")


def get_age(patient: Dict[str, Any]) -> Optional[float]:
    return _to_float(patient.get("age"))


# def fib4_score(age, ast, platelets, alt) -> Optional[float]:
#     """FIB-4 = (age * AST) / (platelets * sqrt(ALT))"""
#     age = _to_float(age)
#     ast = _to_float(ast)
#     platelets = _to_float(platelets)
#     alt = _to_float(alt)
#     if None in (age, ast, platelets, alt) or platelets == 0 or alt <= 0:
#         return None
#     try:
#         return (age * ast) / (platelets * math.sqrt(alt))
#     except (ValueError, ZeroDivisionError):
#         return None


def _fmt_values(**kwargs) -> str:
    parts = []
    for k, v in kwargs.items():
        if v is None or str(v).strip() == "":
            continue
        parts.append(f"{k} = {v}")
    return "; ".join(parts) if parts else "no data"

def _track(triggered_list: list, name: str, value: Any, threshold: Any, op: str, condition_met: bool):
    if condition_met:
        v = _to_float(value)
        if v is not None:
            triggered_list.append(f"{name}={v} ({op}{threshold})")


def _track_range(triggered_list: list, name: str, value: Any, low: Any, high: Any, condition_met: bool):
    if condition_met:
        v = _to_float(value)
        if v is not None:
            triggered_list.append(f"{name}={v} ({low}-{high})")
# ---------------------------------------------------------------------------
# 2. Thresholds (threshold1 - "Significant Pattern" tier)
# ---------------------------------------------------------------------------

THRESHOLDS = {
    # Endocrine Health
    "fasting_plasma_glucose_high": 126,          # mg/dL
    "hba1c_high": 6.5,                           # %
    "eag_high": 150,                              # mg/dL
    "fasting_insulin_high": 25.0,                # µIU/mL
    "tsh_high_hypothyroid": 10.0,                # µIU/mL
    "free_t4_low_hypothyroid": 0.8,              # ng/dL
    "tsh_low_hyperthyroid": 0.1,                 # µIU/mL
    "free_t4_high_hyperthyroid": 1.8,            # ng/dL
    "free_t3_high_hyperthyroid": 4.4,              # pg/mL
    "age_low": 25,
    # Cardiac Health
    "ldl_high": 160,                             # mg/dL
    "crp_high_cardiac":  3,                     # mg/L
    "lpa_high": 50,                              # mg/dL

    # Respiratory Health
    "eosinophils_high_resp":0.5,                # x10^9/L
    "eosinophils_high_rhinitis": 2.0,               # x10^9/L
    "ige_high_resp": 114,                  # IU/mL
    "neutrophils_high": 7.5,                            # x10^9/L
    "crp_high_resp": 3,                       # mg/L

    # Digestive Health - NAFLD
    "alt_high_men": 30, "alt_high_women": 20,    # U/L
    "ggt_high": 64,                               # U/L
    "triglycerides_high": 150,                   # mg/dL
    "hba1c_high_nafld": 6.5,
    "hdl_low_men": 40, "hdl_low_women": 50,      # mg/dL
    "non_hdl_high": 160,                          # mg/dL
    "fib4_high_35_65": 1.30,
    "fib4_high_65_plus": 2.00,
    "ast_alt_ratio_low": 1,
    "age_35_65": 35,
    "age_65_plus": 65,

    # Digestive Health - IBD
    "crp_high_ibd": 10,                     # mg/L 
    "esr_high_ibd": 30,                      #mm/hr 
    "platelets_high_ibd": 450,                 #x10^ 3 / µl 
    "hemoglobin_low_men": 13, "hemoglobin_low_women": 12,     # gm/dL
    "nlr_high": 2.5,
    "plr_high": 150,
    "lmr_high": 2.5,
    "rdw_high": 14.5,                       # %
    "mcv_low": 80,                           #fL
    "wbc_high": 11000,                          #cells/cu.mm
    "albumin_low": 3.5,                          #gm/dL 
    "globulin_high": 3.6,                          #gm/dL 

    # Digestive Health - Hemochromatosis
    "tsat_high_men": 45, "tsat_high_women": 55,   # %
    "ferritin_high_men": 300,                     # µg/L 
    "ferritin_high_premeno_women": 200,           # µg/L

    # Skin Health - Eczema
    "eosinophils_high_eczema": 1.5,               # 10^9/L
    "ige_high_eczema": 2000,                # IU/mL
    "crp_high_skin": 3,
    "esr_high_eczema": 40,

    # # Skin Health - Atopic Dermatitis
    # "eosinophils_high_ad": 0.5,                   # 10^9/L
    # "ige_high_ad": 300,                     # IU/mL
    # "esr_high_ad": 15,

    # # Skin Health - Psoriasis (Category always "" per table)
    # "hb_low_psoriasis_men": 13, "hb_low_psoriasis_women": 12,
    # "tlc_high_psoriasis": 11000,
    # "neutrophil_pct_high": 80,
    # "lymphocyte_pct_high": 40,
    # "nlr_high_psoriasis": 2.5,
    # "plr_high_psoriasis": 150,
    # "esr_high_psoriasis": 20,
    # "crp_high_psoriasis": 5,
}


# ---------------------------------------------------------------------------
# 2b. Threshold2 — "Early Pattern" tier
# ---------------------------------------------------------------------------

THRESHOLDS2 = {
    # Endocrine Health
    "fasting_plasma_glucose_borderline_low": 110,
    "fasting_plasma_glucose_borderline_high": 135,
    "hba1c_borderline": 6.5,
    "fasting_plasma_glucose_borderline_high_MODY": 99,
    "hba1c_borderline_MODY": 5.6,
    "tsh_borderline_hypothyroid": (4.0 , 10.0),
    "tsh_borderline_hyperthyroid": (0.1 , 0.4),
    "free_t4_range": (0.8, 1.8),
    "free_t3_range": (1.4, 4.4),

    # Cardiac Health
    "lpa_borderline_high": 50,


    # Digestive Health - NAFLD
    "alt_high_men": 30, "alt_high_women": 20,
    "fib4_high_35_65": 1.30,
    "fib4_high_65_plus": 2.00,
    "ast_alt_ratio_high": 1,
    "age_35_65": 35,
    "age_65_plus": 65,
    "triglycerides_borderline": 150,
    "hba1c_borderline_nafld": 6.5,
    "hdl_borderline_men": 40,
    "hdl_borderline_women": 50,

    # Digestive Health - IBD
    "crp_borderline_low": 5, "crp_borderline_high": 10,
    "esr_borderline": 20,

}


# ---------------------------------------------------------------------------
# 3. Evaluators — one per condition
# ---------------------------------------------------------------------------

def evaluate_type2_diabetes(labs, patient):
    # 2. Pull values from labs / patient
    fpg = labs.get("fasting_glucose")
    hba1c = labs.get("hba1c")
    eag = labs.get("estimated_average_glucose_(eag)")
    fasting_insulin = labs.get("fasting_insulin")

    # 3. Check if values are abnormal
    fpg_flag = is_elevated(fpg, THRESHOLDS["fasting_plasma_glucose_high"])
    hba1c_flag = is_elevated(hba1c, THRESHOLDS["hba1c_high"])
    eag_flag = is_elevated(eag, THRESHOLDS["eag_high"])
    insulin_flag = is_elevated(fasting_insulin, THRESHOLDS["fasting_insulin_high"])

    fpg_borderline_flag = in_range(fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_low"], THRESHOLDS2["fasting_plasma_glucose_borderline_high"])
    hba1c_borderline_flag = is_elevated(hba1c, THRESHOLDS2["hba1c_borderline"])

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if fpg_flag and hba1c_flag and eag_flag and insulin_flag:
        category = "Significant Pattern"
        _track(triggered, "fasting_glucose", fpg, THRESHOLDS["fasting_plasma_glucose_high"], ">=", fpg_flag)
        _track(triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high"], ">=", hba1c_flag)
        _track(triggered, "eag", eag, THRESHOLDS["eag_high"], ">=", eag_flag)
        _track(triggered, "fasting_insulin", fasting_insulin, THRESHOLDS["fasting_insulin_high"], ">=", insulin_flag)
    elif fpg_borderline_flag and hba1c_borderline_flag:
        category = "Early Pattern"
        _track(triggered, "fasting_glucose", fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_low"], ">=", fpg_borderline_flag)
        _track(triggered, "hba1c", hba1c, THRESHOLDS2["hba1c_borderline"], ">=", hba1c_borderline_flag)
    else:
        category = "Typical"
        if fpg_flag:
            _track(partial_triggered, "fasting_glucose", fpg, THRESHOLDS["fasting_plasma_glucose_high"], ">=", fpg_flag)
        elif fpg_borderline_flag:
            _track(partial_triggered, "fasting_glucose", fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_low"], ">=", fpg_borderline_flag)
        if hba1c_flag:
            _track(partial_triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high"], ">=", hba1c_flag)
        elif hba1c_borderline_flag:
            _track(partial_triggered, "hba1c", hba1c, THRESHOLDS2["hba1c_borderline"], ">=", hba1c_borderline_flag)
        _track(partial_triggered, "eag", eag, THRESHOLDS["eag_high"], ">=", eag_flag)
        _track(partial_triggered, "fasting_insulin", fasting_insulin, THRESHOLDS["fasting_insulin_high"], ">=", insulin_flag)

    return [{
        "Condition": "Type 2 Diabetes",
        "Category": category,
        "Values": _fmt_values(fasting_glucose=fpg, hba1c=hba1c, eag=eag, fasting_insulin=fasting_insulin),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]
def evaluate_maturity_onset_diabetes_of_the_young(labs, patient):
    # 2. Pull values from labs / patient
    fpg = labs.get("fasting_glucose")
    hba1c = labs.get("hba1c")
    age = labs.get("age")

    # 3. Check if values are abnormal
    fpg_flag = is_elevated(fpg, THRESHOLDS["fasting_plasma_glucose_high"])
    hba1c_flag = is_elevated(hba1c, THRESHOLDS["hba1c_high"])
    age_flag = is_below(age, THRESHOLDS["age_low"])

    fpg_borderline_flag = is_elevated(fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_high_MODY"])
    hba1c_borderline_flag = is_elevated(hba1c, THRESHOLDS2["hba1c_borderline_MODY"])

    # 4. Decide category
    # Note: definitive diagnosis requires genetic testing; this is a lab-pattern flag only.
    triggered = []
    partial_triggered = []
    if age_flag and fpg_flag and hba1c_flag:
        category = "Significant Pattern"
        _track(triggered, "age", age, THRESHOLDS["age_low"], "<", age_flag)
        _track(triggered, "fasting_glucose", fpg, THRESHOLDS["fasting_plasma_glucose_high"], ">=", fpg_flag)
        _track(triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high"], ">=", hba1c_flag)
    elif age_flag and fpg_borderline_flag and hba1c_borderline_flag:
        category = "Early Pattern"
        _track(triggered, "age", age, THRESHOLDS["age_low"], "<", age_flag)
        _track(triggered, "fasting_glucose", fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_high_MODY"], ">=", fpg_borderline_flag)
        _track(triggered, "hba1c", hba1c, THRESHOLDS2["hba1c_borderline_MODY"], ">=", hba1c_borderline_flag)
    else:
        category = "Typical"
        _track(partial_triggered, "age", age, THRESHOLDS["age_low"], "<", age_flag)
        if fpg_flag:
            _track(partial_triggered, "fasting_glucose", fpg, THRESHOLDS["fasting_plasma_glucose_high"], ">=", fpg_flag)
        elif fpg_borderline_flag:
            _track(partial_triggered, "fasting_glucose", fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_high_MODY"], ">=", fpg_borderline_flag)
        if hba1c_flag:
            _track(partial_triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high"], ">=", hba1c_flag)
        elif hba1c_borderline_flag:
            _track(partial_triggered, "hba1c", hba1c, THRESHOLDS2["hba1c_borderline_MODY"], ">=", hba1c_borderline_flag)

    return [{
        "Condition": "Maturity-Onset Diabetes of the Young",
        "Category": category,
        "Values": _fmt_values(fasting_glucose=fpg, hba1c=hba1c, age=age),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]

def evaluate_hypothyroidism(labs, patient):
    # 2. Pull values from labs / patient
    tsh = labs.get("tsh")
    free_t4 = labs.get("free_t4")

    # 3. Check if values are abnormal
    tsh_flag = is_elevated(tsh, THRESHOLDS["tsh_high_hypothyroid"])
    t4_flag = is_below(free_t4, THRESHOLDS["free_t4_low_hypothyroid"])

    tsh_borderline_flag = in_range(
        tsh, THRESHOLDS2["tsh_borderline_hypothyroid"][0],
        THRESHOLDS2["tsh_borderline_hypothyroid"][1]
    )
    t4_normal_flag = in_range(
        free_t4, THRESHOLDS2["free_t4_range"][0],
        THRESHOLDS2["free_t4_range"][1]
    )

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if tsh_flag and t4_flag:
        category = "Significant Pattern"
        _track(triggered, "tsh", tsh, THRESHOLDS["tsh_high_hypothyroid"], ">=", tsh_flag)
        _track(triggered, "free_t4", free_t4, THRESHOLDS["free_t4_low_hypothyroid"], "<", t4_flag)
    elif tsh_borderline_flag and t4_normal_flag:
        category = "Early Pattern"
        _track_range(triggered, "tsh", tsh, THRESHOLDS2["tsh_borderline_hypothyroid"][0], THRESHOLDS2["tsh_borderline_hypothyroid"][1], tsh_borderline_flag)
        _track_range(triggered, "free_t4", free_t4, THRESHOLDS2["free_t4_range"][0], THRESHOLDS2["free_t4_range"][1], t4_normal_flag)
    else:
        category = "Typical"
        if tsh_flag:
            _track(partial_triggered, "tsh", tsh, THRESHOLDS["tsh_high_hypothyroid"], ">=", tsh_flag)
        elif tsh_borderline_flag:
            _track_range(partial_triggered, "tsh", tsh, THRESHOLDS2["tsh_borderline_hypothyroid"][0], THRESHOLDS2["tsh_borderline_hypothyroid"][1], tsh_borderline_flag)
        if t4_flag:
            _track(partial_triggered, "free_t4", free_t4, THRESHOLDS["free_t4_low_hypothyroid"], "<", t4_flag)
        elif t4_normal_flag:
            _track_range(partial_triggered, "free_t4", free_t4, THRESHOLDS2["free_t4_range"][0], THRESHOLDS2["free_t4_range"][1], t4_normal_flag)


    return [{
        "Condition": "Hypothyroidism",
        "Category": category,
        "Values": _fmt_values(tsh=tsh, free_t4=free_t4),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]




def evaluate_hyperthyroidism(labs, patient):
    # 2. Pull values from labs / patient
    tsh = labs.get("tsh")
    free_t4 = labs.get("free_t4")
    free_t3 = labs.get("free_t3")

    # 3. Check if values are abnormal
    tsh_flag = is_below(tsh, THRESHOLDS["tsh_low_hyperthyroid"])
    t4_flag = is_above(free_t4, THRESHOLDS["free_t4_high_hyperthyroid"])
    t3_flag = is_above(free_t3, THRESHOLDS["free_t3_high_hyperthyroid"])

    tsh_borderline_flag = in_range(
        tsh, THRESHOLDS2["tsh_borderline_hyperthyroid"][0],
        THRESHOLDS2["tsh_borderline_hyperthyroid"][1]
    )
    t4_normal_flag = in_range(
        free_t4, THRESHOLDS2["free_t4_range"][0],
        THRESHOLDS2["free_t4_range"][1]
    )
    t3_normal_flag = in_range(
        free_t3, THRESHOLDS2["free_t3_range"][0],
        THRESHOLDS2["free_t3_range"][1]
    )

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if tsh_flag and t4_flag and t3_flag:
        category = "Significant Pattern"
        _track(triggered, "tsh", tsh, THRESHOLDS["tsh_low_hyperthyroid"], "<", tsh_flag)
        _track(triggered, "free_t4", free_t4, THRESHOLDS["free_t4_high_hyperthyroid"], ">", t4_flag)
        _track(triggered, "free_t3", free_t3, THRESHOLDS["free_t3_high_hyperthyroid"], ">", t3_flag)
    elif (tsh_borderline_flag or tsh_flag) and t4_normal_flag and t3_normal_flag:
        category = "Early Pattern"
        _track_range(triggered, "tsh", tsh, THRESHOLDS2["tsh_borderline_hyperthyroid"][0], THRESHOLDS2["tsh_borderline_hyperthyroid"][1], tsh_borderline_flag)
        _track_range(triggered, "free_t4", free_t4, THRESHOLDS2["free_t4_range"][0], THRESHOLDS2["free_t4_range"][1], t4_normal_flag)
        _track_range(triggered, "free_t3", free_t3, THRESHOLDS2["free_t3_range"][0], THRESHOLDS2["free_t3_range"][1], t3_normal_flag)
    else:
        category = "Typical"
        if tsh_flag:
            _track(partial_triggered, "tsh", tsh, THRESHOLDS["tsh_low_hyperthyroid"], "<", tsh_flag)
        elif tsh_borderline_flag:
            _track_range(partial_triggered, "tsh", tsh, THRESHOLDS2["tsh_borderline_hyperthyroid"][0], THRESHOLDS2["tsh_borderline_hyperthyroid"][1], tsh_borderline_flag)
        if t4_flag:
            _track(partial_triggered, "free_t4", free_t4, THRESHOLDS["free_t4_high_hyperthyroid"], ">", t4_flag)
        elif t4_normal_flag:
            _track_range(partial_triggered, "free_t4", free_t4, THRESHOLDS2["free_t4_range"][0], THRESHOLDS2["free_t4_range"][1], t4_normal_flag)
        if t3_flag:
            _track(partial_triggered, "free_t3", free_t3, THRESHOLDS["free_t3_high_hyperthyroid"], ">", t3_flag)
        elif t3_normal_flag:
            _track_range(partial_triggered, "free_t3", free_t3, THRESHOLDS2["free_t3_range"][0], THRESHOLDS2["free_t3_range"][1], t3_normal_flag)

    return [{
        "Condition": "Hyperthyroidism",
        "Category": category,
        "Values": _fmt_values(tsh=tsh, free_t4=free_t4, free_t3=free_t3),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]


def evaluate_coronary_artery_disease(labs, patient):
    # 2. Pull values from labs / patient
    ldl_c = labs.get("ldl_cholesterol")
    crp = labs.get("crp")
    lpa = labs.get("lpa")

    # 3. Check if values are abnormal
    ldl_flag = is_elevated(ldl_c, THRESHOLDS["ldl_high"])
    crp_flag = is_elevated(crp, THRESHOLDS["crp_high_cardiac"])
    lpa_flag = is_elevated(lpa, THRESHOLDS["lpa_high"])

    lpa_borderline_flag = is_elevated(lpa, THRESHOLDS2["lpa_borderline_high"])

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if ldl_flag and crp_flag and lpa_flag:
        category = "Significant Pattern"
        _track(triggered, "ldl_cholesterol", ldl_c, THRESHOLDS["ldl_high"], ">=", ldl_flag)
        _track(triggered, "crp", crp, THRESHOLDS["crp_high_cardiac"], ">=", crp_flag)
        _track(triggered, "lpa", lpa, THRESHOLDS["lpa_high"], ">=", lpa_flag)
    elif lpa_borderline_flag:
        category = "Elevated susceptibility"
        _track(triggered, "lpa", lpa, THRESHOLDS2["lpa_borderline_high"], ">=", lpa_borderline_flag)
    else:
        category = "Typical"
        _track(partial_triggered, "ldl_cholesterol", ldl_c, THRESHOLDS["ldl_high"], ">=", ldl_flag)
        _track(partial_triggered, "crp", crp, THRESHOLDS["crp_high_cardiac"], ">=", crp_flag)
        if lpa_flag:
            _track(partial_triggered, "lpa", lpa, THRESHOLDS["lpa_high"], ">=", lpa_flag)
        elif lpa_borderline_flag:
            _track(partial_triggered, "lpa", lpa, THRESHOLDS2["lpa_borderline_high"], ">=", lpa_borderline_flag)

    return [{
        "Condition": "Coronary Artery Disease",
        "Category": category,
        "Values": _fmt_values(ldl_cholesterol=ldl_c, crp=crp, lpa=lpa),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]


def evaluate_asthma(labs, patient):
    # 2. Pull values from labs / patient
    eosinophils = labs.get("eosinophils")
    ige = labs.get("ige")
    neutrophils = labs.get("neutrophils")
    crp = labs.get("crp")

    # 3. Check if values are abnormal
    eos_flag = is_elevated(eosinophils, THRESHOLDS["eosinophils_high_resp"])
    ige_flag = is_above(ige, THRESHOLDS["ige_high_resp"])
    neutrophils_flag = is_elevated(neutrophils, THRESHOLDS["neutrophils_high"])
    crp_flag = is_above(crp, THRESHOLDS["crp_high_resp"])

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if eos_flag and ige_flag and neutrophils_flag and crp_flag:
        category = "Elevated susceptibility"
        _track(triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_resp"], ">=", eos_flag)
        _track(triggered, "ige", ige, THRESHOLDS["ige_high_resp"], ">", ige_flag)
        _track(triggered, "neutrophils", neutrophils, THRESHOLDS["neutrophils_high"], ">=", neutrophils_flag)
        _track(triggered, "crp", crp, THRESHOLDS["crp_high_resp"], ">", crp_flag)
    else:
        category = "Typical"
        _track(partial_triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_resp"], ">=", eos_flag)
        _track(partial_triggered, "ige", ige, THRESHOLDS["ige_high_resp"], ">", ige_flag)
        _track(partial_triggered, "neutrophils", neutrophils, THRESHOLDS["neutrophils_high"], ">=", neutrophils_flag)
        _track(partial_triggered, "crp", crp, THRESHOLDS["crp_high_resp"], ">", crp_flag)

    return [{
        "Condition": "Asthma",
        "Category": category,
        "Values": _fmt_values(eosinophils=eosinophils, ige=ige, neutrophils=neutrophils, crp=crp),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]


def evaluate_allergic_rhinitis(labs, patient):
    # 2. Pull values from labs / patient
    eosinophils = labs.get("eosinophils")
    ige = labs.get("ige")

    # 3. Check if values are abnormal
    eos_flag = in_range(eosinophils, THRESHOLDS["eosinophils_high_resp"] , THRESHOLDS["eosinophils_high_rhinitis"])
    ige_flag = is_above(ige, THRESHOLDS["ige_high_resp"])

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if eos_flag and ige_flag:
        category = "Elevated susceptibility"
        _track(triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_resp"], ">=", eos_flag)
        _track(triggered, "ige", ige, THRESHOLDS["ige_high_resp"], ">", ige_flag)
    else:
        category = "Typical"
        _track(partial_triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_resp"], ">=", eos_flag)
        _track(partial_triggered, "ige", ige, THRESHOLDS["ige_high_resp"], ">", ige_flag)
        
    return [{
        "Condition": "Rhinitis",
        "Category": category,
        "Values": _fmt_values(eosinophils=eosinophils, ige=ige),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]


def evaluate_nafld(labs, patient):
    # 2. Pull values from labs / patient
    alt = labs.get("alt")
    ast = labs.get("ast")
    ggt = labs.get("ggt")
    triglycerides = labs.get("triglycerides")
    hba1c = labs.get("hba1c")
    hdl_c = labs.get("hdl_c")
    # platelets = labs.get("platelets")
    non_hdl_c = labs.get("non_hdl_c")
    fib4 = labs.get("fib4")
    age = get_age(patient)
    ast_alt_ratio = labs.get("ast/alt")

    # 3. Check if values are abnormal
    alt_threshold = THRESHOLDS["alt_high_men"] if is_male(patient) else THRESHOLDS["alt_high_women"]
    alt_flag = is_above(alt, alt_threshold)

    # ast_alt_ratio = None
    # if _to_float(ast) is not None and _to_float(alt) not in (None, 0):
    #     ast_alt_ratio = _to_float(ast) / _to_float(alt)
    ast_alt_ratio_flag = is_above(ast_alt_ratio, THRESHOLDS["ast_alt_ratio_low"])

    ggt_flag = is_above(ggt, THRESHOLDS["ggt_high"])
    tg_flag = is_elevated(triglycerides, THRESHOLDS["triglycerides_high"])
    hba1c_flag = is_elevated(hba1c, THRESHOLDS["hba1c_high_nafld"])

    hdl_threshold = THRESHOLDS["hdl_low_men"] if is_male(patient) else THRESHOLDS["hdl_low_women"]
    hdl_c_flag = is_below(hdl_c, hdl_threshold)

    non_hdl_c_flag = is_above(non_hdl_c, THRESHOLDS["non_hdl_high"])

     # FIB-4 score: >= 1.30 (age 35-65), or >= 2.0 if age >= 65
    # fib4 = fib4_score(age, ast, platelets, alt)
    fib4_threshold = THRESHOLDS["fib4_high_65_plus"] if (age is not None and age >= 65) else THRESHOLDS["fib4_high_35_65"]
    # fib4_flag = fib4 is not None and fib4 >= fib4_threshold
    fib4_flag = is_elevated(fib4, fib4_threshold)

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if (alt_flag and ast_alt_ratio_flag and ggt_flag and tg_flag
            and hba1c_flag and hdl_c_flag and non_hdl_c_flag and fib4_flag):
        category = "Significant Pattern"
        _track(triggered, "alt", alt, alt_threshold, ">", alt_flag)
        _track(triggered, "ggt", ggt, THRESHOLDS["ggt_high"], ">", ggt_flag)
        _track(triggered, "triglycerides", triglycerides, THRESHOLDS["triglycerides_high"], ">=", tg_flag)
        _track(triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high_nafld"], ">=", hba1c_flag)
        _track(triggered, "hdl_c", hdl_c, hdl_threshold, "<=", hdl_c_flag)
        _track(triggered, "non_hdl_c", non_hdl_c, THRESHOLDS["non_hdl_high"], ">", non_hdl_c_flag)
        _track(triggered, "fib4", fib4, fib4_threshold, ">=", fib4_flag)
        _track(triggered, "ast/alt ratio", ast_alt_ratio, THRESHOLDS["ast_alt_ratio_low"], ">", ast_alt_ratio_flag)
    elif alt_flag and ast_alt_ratio_flag and fib4_flag and tg_flag and hdl_c_flag and hba1c_flag:
        category = "Early Pattern"
        _track(triggered, "alt", alt, alt_threshold, ">", alt_flag)
        _track(triggered, "ast/alt ratio", ast_alt_ratio, THRESHOLDS["ast_alt_ratio_low"], ">", ast_alt_ratio_flag)
        _track(triggered, "fib4", fib4, fib4_threshold, ">=", fib4_flag)
        _track(triggered, "triglycerides", triglycerides, THRESHOLDS["triglycerides_high"], ">=", tg_flag)
        _track(triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high_nafld"], ">=", hba1c_flag)
        _track(triggered, "hdl_c", hdl_c, hdl_threshold, "<=", hdl_c_flag)
    else:
        category = "Typical"
        _track(partial_triggered, "alt", alt, alt_threshold, ">", alt_flag)
        _track(partial_triggered, "ggt", ggt, THRESHOLDS["ggt_high"], ">", ggt_flag)
        _track(partial_triggered, "triglycerides", triglycerides, THRESHOLDS["triglycerides_high"], ">=", tg_flag)
        _track(partial_triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high_nafld"], ">=", hba1c_flag)
        _track(partial_triggered, "hdl_c", hdl_c, hdl_threshold, "<=", hdl_c_flag)
        _track(partial_triggered, "non_hdl_c", non_hdl_c, THRESHOLDS["non_hdl_high"], ">", non_hdl_c_flag)
        _track(partial_triggered, "fib4", fib4, fib4_threshold, ">=", fib4_flag)
        _track(partial_triggered, "ast/alt ratio", ast_alt_ratio, THRESHOLDS["ast_alt_ratio_low"], ">", ast_alt_ratio_flag)

    return [{
        "Condition": "NAFLD",
        "Category": category,
        "Values": _fmt_values(alt=alt, ast=ast, ast_alt_ratio=ast_alt_ratio, ggt=ggt, triglycerides=triglycerides, hba1c=hba1c, hdl_c=hdl_c, non_hdl_c=non_hdl_c, fib4=fib4),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]

def evaluate_inflammatory_bowel_disease(labs, patient):
    # 2. Pull values from labs / patient
    crp = labs.get("crp")
    esr = labs.get("esr")
    platelets = labs.get("platelets")
    hemoglobin = labs.get("hemoglobin")
    albumin = labs.get("albumin")
    wbc = labs.get("wbc")
    mcv = labs.get("mcv")
    rdw = labs.get("rdw")
    Globulin = labs.get("globulin")
    nlr = labs.get("nlr")
    plr = labs.get("platelet/lymphocyte_ratio")
    lmr = labs.get("lymphocyte/monocyte_ratio")

    # 3. Check if values are abnormal
    crp_flag = is_above(crp, THRESHOLDS["crp_high_ibd"])
    esr_flag = is_above(esr, THRESHOLDS["esr_high_ibd"])
    platelets_flag = is_above(platelets, THRESHOLDS["platelets_high_ibd"])
    hemoglobin_threshold = THRESHOLDS["hemoglobin_low_men"] if is_male(patient) else THRESHOLDS["hemoglobin_low_women"]
    hemoglobin_flag = is_below(hemoglobin, hemoglobin_threshold)
    albumin_flag = is_below(albumin, THRESHOLDS["albumin_low"])
    wbc_flag = is_above(wbc, THRESHOLDS["wbc_high"])
    mcv_flag = is_below(mcv, THRESHOLDS["mcv_low"])
    rdw_flag = is_above(rdw, THRESHOLDS["rdw_high"])
    mcv_rdw_abnormal_flag = mcv_flag and rdw_flag
    Globulin_flag = is_above(Globulin, THRESHOLDS["globulin_high"])
    nlr_flag = is_above(nlr, THRESHOLDS["nlr_high"])
    plr_flag = is_above(plr, THRESHOLDS["plr_high"])
    lmr_flag = is_above(lmr, THRESHOLDS["lmr_high"])
    nlr_plr_lmr_flag = nlr_flag and plr_flag and lmr_flag

    crp_borderline_flag = in_range(crp, THRESHOLDS2["crp_borderline_low"], THRESHOLDS2["crp_borderline_high"])
    esr_borderline_flag = is_above(esr, THRESHOLDS2["esr_borderline"])

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if (crp_flag and esr_flag and platelets_flag
            and hemoglobin_flag and nlr_plr_lmr_flag and albumin_flag and Globulin_flag
            and wbc_flag and mcv_rdw_abnormal_flag):
        category = "Significant Pattern"
        _track(triggered, "crp", crp, THRESHOLDS["crp_high_ibd"], ">", crp_flag)
        _track(triggered, "esr", esr, THRESHOLDS["esr_high_ibd"], ">", esr_flag)
        _track(triggered, "platelets", platelets, THRESHOLDS["platelets_high_ibd"], ">", platelets_flag)
        _track(triggered, "hemoglobin", hemoglobin, hemoglobin_threshold, "<", hemoglobin_flag)
        _track(triggered, "albumin", albumin, THRESHOLDS["albumin_low"], "<", albumin_flag)
        _track(triggered, "wbc", wbc, THRESHOLDS["wbc_high"], ">", wbc_flag)
        _track(triggered, "mcv", mcv, THRESHOLDS["mcv_low"], "<", mcv_flag)
        _track(triggered, "rdw", rdw, THRESHOLDS["rdw_high"], ">", rdw_flag)
        _track(triggered, "nlr", nlr, THRESHOLDS["nlr_high"], ">", nlr_flag)
        _track(triggered, "plr", plr, THRESHOLDS["plr_high"], ">", plr_flag)
        _track(triggered, "lmr", lmr, THRESHOLDS["lmr_high"], ">", lmr_flag)
        _track(triggered, "Globulin", Globulin, THRESHOLDS["globulin_high"], ">", Globulin_flag)
        
    elif (crp_borderline_flag and esr_borderline_flag and platelets_flag
            and hemoglobin_flag and albumin_flag and wbc_flag and mcv_rdw_abnormal_flag):
        category = "Early Pattern"
        _track_range(triggered, "crp", crp, THRESHOLDS2["crp_borderline_low"], THRESHOLDS2["crp_borderline_high"], crp_borderline_flag)
        _track(triggered, "esr", esr, THRESHOLDS2["esr_borderline"], ">", esr_borderline_flag)
        _track(triggered, "platelets", platelets, THRESHOLDS["platelets_high_ibd"], ">", platelets_flag)
        _track(triggered, "hemoglobin", hemoglobin, hemoglobin_threshold, "<", hemoglobin_flag)
        _track(triggered, "albumin", albumin, THRESHOLDS["albumin_low"], "<", albumin_flag)
        _track(triggered, "wbc", wbc, THRESHOLDS["wbc_high"], ">", wbc_flag)
        _track(triggered, "mcv", mcv, THRESHOLDS["mcv_low"], "<", mcv_flag)
        _track(triggered, "rdw", rdw, THRESHOLDS["rdw_high"], ">", rdw_flag)
        
    else:
        category = "Typical"
        if crp_flag:
            _track(partial_triggered, "crp", crp, THRESHOLDS["crp_high_ibd"], ">", crp_flag)
        elif crp_borderline_flag:
            _track_range(partial_triggered, "crp", crp, THRESHOLDS2["crp_borderline_low"], THRESHOLDS2["crp_borderline_high"], crp_borderline_flag)
        if esr_flag:
            _track(partial_triggered, "esr", esr, THRESHOLDS["esr_high_ibd"], ">", esr_flag)
        elif esr_borderline_flag:
            _track(partial_triggered, "esr", esr, THRESHOLDS2["esr_borderline"], ">", esr_borderline_flag)
        _track(partial_triggered, "platelets", platelets, THRESHOLDS["platelets_high_ibd"], ">", platelets_flag)
        _track(partial_triggered, "hemoglobin", hemoglobin, hemoglobin_threshold, "<", hemoglobin_flag)
        _track(partial_triggered, "albumin", albumin, THRESHOLDS["albumin_low"], "<", albumin_flag)
        _track(partial_triggered, "wbc", wbc, THRESHOLDS["wbc_high"], ">", wbc_flag)
        _track(partial_triggered, "mcv", mcv, THRESHOLDS["mcv_low"], "<", mcv_flag)
        _track(partial_triggered, "rdw", rdw, THRESHOLDS["rdw_high"], ">", rdw_flag)
        _track(partial_triggered, "nlr", nlr, THRESHOLDS["nlr_high"], ">", nlr_flag)
        _track(partial_triggered, "plr", plr, THRESHOLDS["plr_high"], ">", plr_flag)
        _track(partial_triggered, "lmr", lmr, THRESHOLDS["lmr_high"], ">", lmr_flag)
        _track(partial_triggered, "Globulin", Globulin, THRESHOLDS["globulin_high"], ">", Globulin_flag)
 
    return [{
        "Condition": "Inflammatory Bowel Disease",
        "Category": category,
        "Values": _fmt_values(crp=crp, esr=esr, platelets=platelets, hemoglobin=hemoglobin, albumin=albumin, wbc=wbc, mcv=mcv, rdw=rdw, nlr=nlr, plr=plr, lmr=lmr),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]


def evaluate_hereditary_hemochromatosis(labs, patient):
    # 2. Pull values from labs / patient
    tsat = labs.get("transferrin_saturation")
    ferritin = labs.get("ferritin")

    # 3. Check if values are abnormal
    tsat_threshold = THRESHOLDS["tsat_high_men"] if is_male(patient) else THRESHOLDS["tsat_high_women"]
    tsat_flag = is_above(tsat, tsat_threshold)

    ferritin_threshold = THRESHOLDS["ferritin_high_men"] if is_male(patient) else THRESHOLDS["ferritin_high_premeno_women"]
    ferritin_flag = is_above(ferritin, ferritin_threshold)

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if tsat_flag and ferritin_flag:
        category = "Significant Pattern"
        _track(triggered, "transferrin_saturation", tsat, tsat_threshold, ">", tsat_flag)
        _track(triggered, "ferritin", ferritin, ferritin_threshold, ">", ferritin_flag)
    else:
        category = "Typical"
        _track(partial_triggered, "transferrin_saturation", tsat, tsat_threshold, ">", tsat_flag)
        _track(partial_triggered, "ferritin", ferritin, ferritin_threshold, ">", ferritin_flag)

    return [{
        "Condition": "Hereditary Hemochromatosis",
        "Category": category,
        "Values": _fmt_values(transferrin_saturation=tsat, ferritin=ferritin),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]


def evaluate_eczema(labs, patient):
    # 2. Pull values from labs / patient
    eosinophils = labs.get("eosinophils")
    ige = labs.get("ige")
    crp = labs.get("crp")
    esr = labs.get("esr")

    # 3. Check if values are abnormal
    eos_flag = is_elevated(eosinophils, THRESHOLDS["eosinophils_high_eczema"])
    ige_flag = is_above(ige, THRESHOLDS["ige_high_eczema"])
    crp_flag = is_above(crp, THRESHOLDS["crp_high_skin"])
    esr_flag = is_above(esr, THRESHOLDS["esr_high_eczema"])

    # 4. Decide category
    triggered = []
    partial_triggered = []
    if eos_flag and ige_flag and crp_flag and esr_flag:
        category = "Significant Pattern"
        _track(triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_eczema"], ">=", eos_flag)
        _track(triggered, "ige", ige, THRESHOLDS["ige_high_eczema"], ">", ige_flag)
        _track(triggered, "crp", crp, THRESHOLDS["crp_high_skin"], ">", crp_flag)
        _track(triggered, "esr", esr, THRESHOLDS["esr_high_eczema"], ">", esr_flag)
    else:
        category = "Typical"
        _track(partial_triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_eczema"], ">=", eos_flag)
        _track(partial_triggered, "ige", ige, THRESHOLDS["ige_high_eczema"], ">", ige_flag)
        _track(partial_triggered, "crp", crp, THRESHOLDS["crp_high_skin"], ">", crp_flag)
        _track(partial_triggered, "esr", esr, THRESHOLDS["esr_high_eczema"], ">", esr_flag)
        
    return [{
        "Condition": "Atopic Dermatitis/Eczema",
        "Category": category,
        "Values": _fmt_values(eosinophils=eosinophils, ige=ige, crp=crp, esr=esr),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]


# def evaluate_atopic_dermatitis(labs, patient):
#     # 2. Pull values from labs / patient
#     eosinophils = labs.get("eosinophils")
#     ige = labs.get("ige")
#     crp = labs.get("crp")
#     esr = labs.get("esr")

#     # 3. Check if values are abnormal
#     eos_flag = is_elevated(eosinophils, THRESHOLDS["eosinophils_high_eczema"])
#     ige_flag = is_above(ige, THRESHOLDS["ige_high_eczema"])
#     crp_flag = is_above(crp, THRESHOLDS["crp_high_skin"])
#     esr_flag = is_above(esr, THRESHOLDS["esr_high_eczema"])

#     # 4. Decide category
#     triggered = []
#     partial_triggered = []
#     if eos_flag and ige_flag and crp_flag and esr_flag:
#         category = "Significant Pattern"
#         _track(triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_eczema"], ">=", eos_flag)
#         _track(triggered, "ige", ige, THRESHOLDS["ige_high_eczema"], ">", ige_flag)
#         _track(triggered, "crp", crp, THRESHOLDS["crp_high_skin"], ">", crp_flag)
#         _track(triggered, "esr", esr, THRESHOLDS["esr_high_eczema"], ">", esr_flag)
#     else:
#         category = "Typical"
#         _track(partial_triggered, "eosinophils", eosinophils, THRESHOLDS["eosinophils_high_eczema"], ">=", eos_flag)
#         _track(partial_triggered, "ige", ige, THRESHOLDS["ige_high_eczema"], ">", ige_flag)
#         _track(partial_triggered, "crp", crp, THRESHOLDS["crp_high_skin"], ">", crp_flag)
#         _track(partial_triggered, "esr", esr, THRESHOLDS["esr_high_eczema"], ">", esr_flag)

#     return [{
#         "Condition": "Atopic Dermatitis",
#         "Category": category,
#         "Values": _fmt_values(eosinophils=eosinophils, ige=ige, crp=crp, esr=esr),
#         "TriggeredValues": ", ".join(triggered),
#         "PartialTriggered": ", ".join(partial_triggered)
#     }]


# ---------------------------------------------------------------------------
# 4. Registry of all evaluators
# ---------------------------------------------------------------------------

EVALUATORS = [
    evaluate_type2_diabetes,
    evaluate_hypothyroidism,
    evaluate_maturity_onset_diabetes_of_the_young,
    evaluate_hyperthyroidism,
    evaluate_coronary_artery_disease,
    evaluate_asthma,
    evaluate_allergic_rhinitis,
    evaluate_nafld,
    evaluate_inflammatory_bowel_disease,
    evaluate_hereditary_hemochromatosis,
    evaluate_eczema,
    # evaluate_atopic_dermatitis,
]


def evaluate_all(labs: Dict[str, Any], patient: Dict[str, Any]) -> List[Dict[str, str]]:
    results = []
    for fn in EVALUATORS:
        results.extend(fn(labs, patient))
    return results


# ---------------------------------------------------------------------------
# 5. CSV I/O
# ---------------------------------------------------------------------------

PATIENT_COLUMNS = {"patient_id", "sex", "age"}

# Only these categories are written to the output. Anything else
# (i.e. "" / no pattern met) is filtered out of the report.
OUTPUT_FIELDNAMES = ["Domain", "Condition", "Category", "DNA Marker(s)", "Blood Marker(s)", "Triggering PRS", "Snapshot Category", "All Blood Marker(s)"]
# VISIBLE_CATEGORIES = {"Significant Pattern", "Early Pattern", "Elevated susceptibility", "Typical"}
VISIBLE_CATEGORIES = {"Significant Pattern", "Typical"}
# Columns that indicate a "long format" CSV (one row per lab parameter,
# rather than one column per lab). If both are present we pivot the
# rows into one wide record per patient before evaluating.
LONG_FORMAT_PARAM_COL_CANDIDATES = {"parameter", "test", "test_name", "lab", "lab_name"}
LONG_FORMAT_VALUE_COL_CANDIDATES = {"value", "result", "result_value", "lab_value"}


def _normalize_key(name: str) -> str:
    return str(name).strip().lower().replace(" ", "_")


def _load_rows(input_path: str) -> List[Dict[str, str]]:
    with open(input_path, newline="", encoding="utf-8-sig") as f_in:
        reader = csv.DictReader(f_in)
        return [
            {_normalize_key(k): (v.strip() if isinstance(v, str) else v) for k, v in raw_row.items()}
            for raw_row in reader
        ]


def _pivot_long_format(rows: List[Dict[str, str]], param_col: str, value_col: str) -> List[Dict[str, Any]]:
    """
    Turn long-format rows (one row per parameter/value pair) into one
    wide record per patient. If there's no patient_id column, every row
    is assumed to belong to a single patient (id left blank).
    """
    patients: Dict[str, Dict[str, Any]] = {}
    order: List[str] = []

    for row in rows:
        pid = row.get("patient_id", "") or ""
        if pid not in patients:
            patients[pid] = {"patient_id": pid}
            order.append(pid)
            for pc in PATIENT_COLUMNS:
                if pc in row and row.get(pc):
                    patients[pid][pc] = row[pc]

        param = row.get(param_col)
        value = row.get(value_col)
        if not param:
            continue
        key = _normalize_key(param)
        patients[pid][key] = value

        # If the "parameter" row itself is sex/age (common in long-format
        # exports), also promote it into the patient-level fields so
        # is_male()/is_female()/get_age() can see it.
        if key in PATIENT_COLUMNS:
            patients[pid][key] = value

        # carry forward patient-level columns (sex/age) whenever present
        # as their own CSV columns
        for pc in PATIENT_COLUMNS:
            if pc in row and row.get(pc):
                patients[pid][pc] = row[pc]

    return [patients[pid] for pid in order]


def process_csv(input_path: str, output_path: str) -> None:
    rows = _load_rows(input_path)
    if not rows:
        print("No rows found in input CSV.")
        with open(output_path, "w", newline="", encoding="utf-8") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=OUTPUT_FIELDNAMES)
            writer.writeheader()
        return

    all_cols = set(rows[0].keys())
    param_col = next((c for c in LONG_FORMAT_PARAM_COL_CANDIDATES if c in all_cols), None)
    value_col = next((c for c in LONG_FORMAT_VALUE_COL_CANDIDATES if c in all_cols), None)

    if param_col and value_col:
        print(f"[FORMAT] Detected long-format CSV (parameter column='{param_col}', value column='{value_col}'). Pivoting to one record per patient.")
        wide_rows = _pivot_long_format(rows, param_col, value_col)
    else:
        print("[FORMAT] Detected wide-format CSV (one column per lab). Using rows as-is.")
        wide_rows = rows

    rows_out = []
    for row in wide_rows:
        patient = {k: row.get(k) for k in PATIENT_COLUMNS if k in row}
        labs = {k: v for k, v in row.items() if k not in PATIENT_COLUMNS}

        # --- DEBUG: print exactly what was fetched from the CSV for
        # this patient, BEFORE running any evaluator. Use this to
        # confirm column names/values are being read correctly and
        # that a rule not firing is a threshold issue, not a
        # missing/misnamed-column issue.
        print(f"\n[FETCHED] patient_id={patient.get('patient_id', '')}")
        print(f"  patient: {patient}")
        print(f"  labs:")
        for k, v in labs.items():
            print(f"    {k} = {v!r}")

        results = evaluate_all(labs, patient)

        # --- DEBUG: print every evaluator's result (even "" / no
        # pattern), so you can confirm whether a rule is firing based
        # on the fetched values above.
        print(f"  evaluations:")
        for r in results:
            cat_display = r["Category"] if r["Category"] else "(no pattern)"
            print(f"    {r['Condition']}: {cat_display}  |  {r.get('Values', '')}")

        for r in results:
            # Only show the condition if it hit "Significant Pattern"
            # or "Early Pattern" — skip blank/no-pattern results.
            if r["Category"] not in VISIBLE_CATEGORIES:
                continue
            cat = r["Category"]
            domain = CONDITION_DOMAINS.get(r["Condition"], "General Health")
            mapped_category = _map_category_message(cat)
            
            blood_marker = r.get("TriggeredValues", "") if cat in {"Significant Pattern", "Early Pattern", "Elevated susceptibility"} else ""
            all_markers = r.get("TriggeredValues") or r.get("PartialTriggered") or ""

            rows_out.append({
                "Domain": domain,
                "Condition": r["Condition"],
                "Category": cat,
                "DNA Marker(s)": "",  # Placeholder for future integration with genetic data
                "Blood Marker(s)": blood_marker,
                "Triggering PRS": "",
                "Snapshot Category": mapped_category,
                "All Blood Marker(s)": all_markers,
            })

    with open(output_path, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=OUTPUT_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows_out)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python condition_evaluators.py <input_labs.csv> <output_results.csv>")
        sys.exit(1)

    process_csv(sys.argv[1], sys.argv[2])
    print(f"Done. Results written to {sys.argv[2]}")