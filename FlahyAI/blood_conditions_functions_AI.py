# dead weight CLI approach
# import argparse
import functools
import inspect
import json
# dead weight CLI approach
# import csv
import os

# example run command:
# python script.py monogenic.json prs.json apoe.json blood.csv results.csv --thresholds thresholds_lab_a.json --sex male --age 35
#
# This module is also directly importable (e.g. to register individual
# evaluate_* functions as model tool-calls) without running as a script -
# THRESHOLDS is loaded with a sensible default at import time below, not
# only inside the CLI block, and the CLI-only pipeline run (reading
# args.*, writing output_csv) is confined to the bottom
# `if __name__ == "__main__":` block so importing never requires `args`
# or tries to run the full batch pipeline.
#
# family_history/past_history/symptoms are passed in directly by the
# caller (e.g. a model's tool-call, which selects matching terms from
# blood_conditions_user_context.json's vocab itself, enforced by that
# tool's schema) - this script does no vocab matching of its own.
#
# USER_CONTEXT below is read only to populate each evaluate_*_schema's
# enum lists at import time, straight from the JSON - never copied into
# the schemas as frozen literals, so editing the JSON (or CONDITION_GENES
# further below) changes what every schema offers without touching this file.

# checks these files in same directory as script by default
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_THRESHOLDS_PATH = os.path.join(SCRIPT_DIR, "metropolis_thresholds.json")

with open(os.path.join(SCRIPT_DIR, "blood_conditions_user_context.json")) as f:
    USER_CONTEXT = json.load(f)

def load_thresholds(path=DEFAULT_THRESHOLDS_PATH):
    """(Re)loads the module-level THRESHOLDS dict. Called once below with
    the default lab's thresholds so this module's functions are directly
    callable without running as __main__; --thresholds overrides this
    when run as a CLI script."""
    global THRESHOLDS
    with open(path) as f:
        THRESHOLDS = json.load(f)


load_thresholds()  # default load, so the _t()-based constants further below always have a value


def _t(key):  # added
    """Looks up key in the loaded lab thresholds JSON, returning a tuple if the value is a list."""
    value = THRESHOLDS[key]
    return tuple(value) if isinstance(value, list) else value

PATHOGENIC_VALUES = {"pathogenic", "likely pathogenic", "likely_pathogenic", "pathogenic/likely pathogenic", "pathogenic/likely_pathogenic"}
ELEVATED_CATEGORIES = {"elevated", "moderately_elevated", "moderately elevated"}

# dead weight CLI approach
# # genentics part from result jsons
# def build_genetics_from_jsons(monogenic_json_path, prs_json_path, apoe_json_path):
#     """
#     Reads the three JSON files (monogenic pathogenic variants, polygenic PRS
#     scores, APOE status) and returns a single genetics dict:
#         {
#             "flagged_genes": [...],              # monogenic pathogenic/likely pathogenic hits (+ APOE if e4 present)
#             "prs_elevated_conditions": {...},     # set of PRS json's condition keys where category is elevated/moderately_elevated
#             "apoe_status": "e3/e4",               # raw APOE genotype string, for reference/future use
#         }
#     """
#     # --- 1. monogenic pathogenic/likely pathogenic gene hits ---
#     flagged_genes = []
#     acmg_genes = []  # subset of flagged_genes
#     monogenic_lookup = {} # (gene, condition_lower) -> (rsid, pathogenicity)
#     with open(monogenic_json_path) as f:
#         monogenic_data = json.load(f)
#     for entry in monogenic_data:
#         sig = (entry.get("ClinVar Significance") or "").strip().lower()
#         if sig in PATHOGENIC_VALUES:
#             gene = (entry.get("Gene Symbol") or "").strip()
#             cond = (entry.get("condition / ACMG Diseases") or "").strip()
#             if gene:
#                 flagged_genes.append(gene)
#                 if entry.get("in_acmg") is True:
#                     acmg_genes.append(gene)
#                 if cond:
#                     monogenic_lookup.setdefault((gene, cond.lower()), []).append((
#                         entry.get("preferredVarTitle", ""),
#                         entry.get("ClinVar Significance", ""),
#                     ))
# 
#     # --- 2. Polygenic PRS elevated conditions ---
#     prs_elevated_conditions = set()
#     prs_categories = {}
#     with open(prs_json_path) as f:
#         prs_data = json.load(f)
#     for condition_key, info in prs_data.items():
#         category = (info.get("category") or "").strip().lower()
#         prs_categories[condition_key] = category
#         if category in ELEVATED_CATEGORIES:
#             prs_elevated_conditions.add(condition_key)
# 
#     # --- 3. APOE status ---
#     apoe_status = ""
#     with open(apoe_json_path) as f:
#         apoe_data = json.load(f)
#     apoe_status = (apoe_data.get("APOE_Status") or "").strip()
#     if "4" in apoe_status:  # e.g. "e3/e4" or "e4/e4" -> ε4 allele present
#         flagged_genes.append("APOE")
#         # acmg_genes.append("APOE") # apoe now not assigned to acmg list by default.
# 
#     return {
#         "flagged_genes": flagged_genes,
#         "acmg_genes": acmg_genes,
#         "prs_elevated_conditions": prs_elevated_conditions,
#         "prs_categories": prs_categories,
#         "apoe_status": apoe_status,
#         "monogenic_lookup": monogenic_lookup,
#     }

# match prs json syntax with condition's names used here
PRS_CONDITION_KEY_MAP = {
    "alzheimers": "Alzheimer's Disease",
    "cad": "Coronary Artery Disease",
    # "asthma": "Asthma",
    "atopicdermatitis": "Atopic Dermatitis",
    "eczema": "Eczema",
    "copd": "COPD",
    "hyperthyroidism": "Hyperthyroidism",
    "hypothyroidism": "Hypothyroidism",
    "ibd": "Inflammatory Bowel Disease",
    "nafld": "NAFLD",
    # "osteoarthritis": "Osteoarthritis",
    # "parkinson": "Parkinson's Disease",
    "psoriasis": "Psoriasis",
    "rheumatoidarthritis": "Rheumatoid Arthritis",
    # "rhinitis": "Rhinitis",
    "type2diabetes": "Type 2 Diabetes",
}

# dead weight CLI approach
# # combining genentic (gene, prs), blood (labs), and patient info (sex, age)
# def load_patient_data(monogenic_json, prs_json, apoe_json, blood_csv_path, sex, age):
#     """
#     Patient-reported context (family_history, symptoms, past_history) is 
#     NOT included here. sex/age are read from the blood csv by default. 
#     the args based --sex/--age on command line overrides csv values.
#     """
#     # -----------------------------
#     # Blood report
#     # -----------------------------
#     labs = {}
#     csv_sex = None
#     csv_age = None
#     with open(blood_csv_path) as f:
#         reader = csv.DictReader(f)
#         for row in reader:
#             key = row["parameter"].strip().lower()
#             # value = row["value"].strip()
#             value = row["value"].strip().replace(",", "")  #replace commas in numbers (e.g. 1,234.56) to avoid float conversion errors
# 
#             if key == "sex":
#                 csv_sex = value.strip().lower() or None
#                 continue
#             if key == "age":
#                 csv_age = int(float(value)) if value else None
#                 continue
# 
#             labs[key] = float(value) if value else None
# 
#     # Approach changed for hereditary hemochromatosis, iron profile is now category assigner and not just eligibility gate.
#     # Derive iron_profile_done: True only if BOTH ferritin and
#     # transferrin_saturation were present (non-null) in the blood CSV (for evaluate_hemochromatosis).
#     # labs["iron_profile_done"] = (
#     #     labs.get("ferritin") is not None and labs.get("transferrin_saturation") is not None
#     # )
#     
#     genetics = build_genetics_from_jsons(monogenic_json, prs_json, apoe_json)
#     
#     final_sex = sex or csv_sex
#     final_age = age if age is not None else csv_age
# 
#     if final_sex not in SEX_VALUES:
#         raise ValueError(f"Sex must be 'male' or 'female', got: {final_sex!r} (from --sex or csv)")
#     if final_age is None:
#         raise ValueError("Age not found in blood CSV and not provided via --age")
# 
#     patient = {"sex": final_sex, "age": final_age}
# 
#     # just to check
#     print(f"Patient: {patient}")
#     print(f"Genetics: {genetics}")
#     print(f"Labs: {labs}")
# 
#     return {
#         "patient": patient,
#         "genetics": genetics,
#         "labs": labs,
#     }

# =============================================================================
# Shared lab reference ranges / thresholds
# =============================================================================

# dead weight CLI approach
# SEX_VALUES = ["male", "female"]

# GENETIC_RISK_CLASSIFICATIONS = ["Elevated", "Moderately Elevated", "Typical", "Not Reported"]
# ELEVATED_CLASSIFICATIONS = {"Elevated", "Moderately Elevated"}

# All values below are now read from THRESHOLDS (loaded from --thresholds JSON, or
# thresholds_default.json) instead of being hardcoded, so a different lab's reference
# ranges/thresholds can be swapped in via the CLI arg without touching this file or any
# of the functions that use these names. # added

# --- Lipid Profile ---
TOTAL_CHOLESTEROL_AT_RISK = _t("TOTAL_CHOLESTEROL_AT_RISK")       # mg/dL
TOTAL_CHOLESTEROL_ELEVATED = _t("TOTAL_CHOLESTEROL_ELEVATED")     # mg/dL
LDL_C_NORMAL_CAD = _t("LDL_C_NORMAL_CAD")                         # mg/dL, preferred ceiling if CAD susceptibility present
LDL_C_AT_RISK = _t("LDL_C_AT_RISK")                               # mg/dL
LDL_C_ELEVATED = _t("LDL_C_ELEVATED")                             # mg/dL
LDL_C_SEVERE = _t("LDL_C_SEVERE")                                 # mg/dL
HDL_C_MALE_MIN = _t("HDL_C_MALE_MIN")                             # mg/dL
HDL_C_FEMALE_MIN = _t("HDL_C_FEMALE_MIN")                         # mg/dL
HDL_C_MALE_SEVERE = _t("HDL_C_MALE_SEVERE")                       # mg/dL
HDL_C_FEMALE_SEVERE = _t("HDL_C_FEMALE_SEVERE")                   # mg/dL
NON_HDL_C_NORMAL = _t("NON_HDL_C_NORMAL")                         # mg/dL
NON_HDL_C_ELEVATED = _t("NON_HDL_C_ELEVATED")                     # mg/dL
NON_HDL_C_SEVERE = _t("NON_HDL_C_SEVERE")                         # mg/dL
VLDL_UPPER_NORMAL = _t("VLDL_UPPER_NORMAL")                       # mg/dL
TRIGLYCERIDES_AT_RISK = _t("TRIGLYCERIDES_AT_RISK")               # mg/dL
TRIGLYCERIDES_ELEVATED = _t("TRIGLYCERIDES_ELEVATED")             # mg/dL
TRIGLYCERIDES_SEVERE = _t("TRIGLYCERIDES_SEVERE")                 # mg/dL
TRIGLYCERIDES_VERY_SEVERE = _t("TRIGLYCERIDES_VERY_SEVERE")       # mg/dL
# TG_HDL_RATIO_ELEVATED = 2.5
APOB_AT_RISK = _t("APOB_AT_RISK")                                 # mg/dL
APOB_ELEVATED = _t("APOB_ELEVATED")                               # mg/dL
LPA_AT_RISK = _t("LPA_AT_RISK")                                   # mg/dL
LPA_ELEVATED = _t("LPA_ELEVATED")                                 # mg/dL

# --- Glucose / Metabolic ---
FASTING_GLUCOSE_AT_RISK = _t("FASTING_GLUCOSE_AT_RISK")           # mg/dL
FASTING_GLUCOSE_ELEVATED = _t("FASTING_GLUCOSE_ELEVATED")         # mg/dL
HBA1C_AT_RISK = _t("HBA1C_AT_RISK")                               # %
HBA1C_ELEVATED = _t("HBA1C_ELEVATED")                             # %
FASTING_INSULIN_RANGE = _t("FASTING_INSULIN_RANGE")               # uIU/mL

# --- Inflammatory Markers ---
CRP_AT_RISK = _t("CRP_AT_RISK")                                   # mg/L
CRP_AT_RISK_COPD = _t("CRP_AT_RISK_COPD")                         # mg/L
CRP_ELEVATED = _t("CRP_ELEVATED")                                 # mg/L
IGE_UPPER_NORMAL = _t("IGE_UPPER_NORMAL")                         # IU/mL
FIBRINOGEN_RANGE = _t("FIBRINOGEN_RANGE")                         # mg/dL (updated from 200-400)
# FECAL_CALPROTECTIN_AT_RISK = 50       # ug/g
# FECAL_CALPROTECTIN_ELEVATED = 200     # ug/g

def esr_upper_normal(age):
    """Returns ESR upper limit of normal (mm/hr) by age."""
    if age is None:
        return None

    # added: age cutoff/limits now come from THRESHOLDS so they can vary per lab too
    if age <= THRESHOLDS["ESR_AGE_CUTOFF"]:
        return THRESHOLDS["ESR_UPPER_NORMAL_YOUNG"]
    else:
        return THRESHOLDS["ESR_UPPER_NORMAL_OLDER"]

# --- CBC ---
EOSINOPHILS_UPPER_NORMAL = _t("EOSINOPHILS_UPPER_NORMAL")         # x10^9/L
EOSINOPHILS_BORDERLINE = _t("EOSINOPHILS_BORDERLINE")             # x10^9/L
EOSINOPHILS_LIKELY = _t("EOSINOPHILS_LIKELY")                     # x10^9/L
EOSINOPHILS_ELEVATED = _t("EOSINOPHILS_ELEVATED")                 # x10^9/L
NEUTROPHILS_RANGE = _t("NEUTROPHILS_RANGE")                       # x10^9/L
# PLATELETS_RANGE = (150, 400)          # x10^3/μL
NLR_RANGE = _t("NLR_RANGE")                                       # ratio
HEMOGLOBIN_MALE_RANGE = _t("HEMOGLOBIN_MALE_RANGE")               # g/dL
HEMOGLOBIN_FEMALE_RANGE = _t("HEMOGLOBIN_FEMALE_RANGE")           # g/dL

# --- Liver Function / Fibrosis ---
AST_RANGE = _t("AST_RANGE")                                       # U/L
ALT_RANGE = _t("ALT_RANGE")                                       # U/L
GGT_RANGE = _t("GGT_RANGE")                                       # U/L
AST_ALT_RATIO_RANGE = _t("AST_ALT_RATIO_RANGE")
FIB4_AT_RISK = _t("FIB4_AT_RISK")
FIB4_ELEVATED = _t("FIB4_ELEVATED")
ALBUMIN_RANGE = _t("ALBUMIN_RANGE")                               # g/dL
# FIBROSIS_SCORE_ELEVATED_THRESHOLD = None  # "fibrosis scores" referenced generically (TM6SF2 rule), no defined metric/cutoff

# --- Thyroid ---
TSH_RANGE = _t("TSH_RANGE")                                       # mIU/L
TSH_SUPPRESSED = _t("TSH_SUPPRESSED")                             # mIU/L
TSH_ELEVATED = _t("TSH_ELEVATED")                                 # mIU/L
TSH_ELEVATED_RANGE = _t("TSH_ELEVATED_RANGE")                     # mIU/L
FREE_T4_RANGE = _t("FREE_T4_RANGE")                               # ng/dL
FREE_T3_RANGE = _t("FREE_T3_RANGE")                               # pg/mL
TPOAB_UPPER_NORMAL = _t("TPOAB_UPPER_NORMAL")                     # IU/mL

# --- Iron studies (HFE) ---
# FERRITIN_ELEVATED_THRESHOLD = None  # prompt says "ferritin abnormal", no numeric cutoff given
# TRANSFERRIN_SATURATION_ELEVATED_THRESHOLD = None  # same as above, no numeric cutoff given
# SERUM_IRON_ELEVATED_THRESHOLD = None  # not referenced with a number anywhere in prompt, included for completeness if needed

CALCIUM_UPPER_NORMAL = _t("CALCIUM_UPPER_NORMAL")                 # mg/dL
PTH_UPPER_NORMAL = _t("PTH_UPPER_NORMAL")                         # pg/mL
CK_UPPER_NORMAL = _t("CK_UPPER_NORMAL")                           # U/L

TOTAL_BILIRUBIN_RANGE = _t("TOTAL_BILIRUBIN_RANGE")               # mg/dL (Gilbert)
DIRECT_BILIRUBIN_LOWER = _t("DIRECT_BILIRUBIN_LOWER")             # mg/dL (Gilbert)
INDIRECT_BILIRUBIN_UPPER = _t("INDIRECT_BILIRUBIN_UPPER")         # mg/dL (Gilbert)
AST_UPPER = _t("AST_UPPER")                                       # U/L (Gilbert)
ALT_UPPER = _t("ALT_UPPER")                                       # U/L (Gilbert)
ALP_RANGE = _t("ALP_RANGE")                                       # U/L (Gilbert)
GGT_LFT_RANGE = _t("GGT_LFT_RANGE")                               # U/L (Gilbert)
FERRITIN_ELEVATED_MALE = _t("FERRITIN_ELEVATED_MALE")             # ng/mL (Hemochromatosis)
FERRITIN_ELEVATED_FEMALE = _t("FERRITIN_ELEVATED_FEMALE")         # ng/mL (Hemochromatosis)
TRANSFERRIN_SATURATION_ELEVATED = _t("TRANSFERRIN_SATURATION_ELEVATED")  # % (Hemochromatosis)
HEMOGLOBIN_UPPER_HEMOCHROMATOSIS = _t("HEMOGLOBIN_UPPER_HEMOCHROMATOSIS")  # g/dL
# MCV_UPPER_HEMOCHROMATOSIS = _t("MCV_UPPER_HEMOCHROMATOSIS")       # fL
# MCH_UPPER_HEMOCHROMATOSIS = _t("MCH_UPPER_HEMOCHROMATOSIS")       # pg
# VITAMIN_D_DEFICIENT = 20                       # ng/mL (Celiac)
# MCV_LOWER_CELIAC = 82                          # fL
# MCH_LOWER_CELIAC = 27                          # pg
# MCHC_LOWER_CELIAC = 31.5                       # g/dL
TTG_RANGE = _t("TTG_RANGE")                      # AU/mL (Celiac)
IGE_EARLY = _t("IGE_EARLY")                      # IU/mL (Atopic/Eczema)
IGE_LIKELY = _t("IGE_LIKELY")                    # IU/mL (Atopic/Eczema)

# =============================================================================
# Shared helper functions
# =============================================================================

def sex_based_threshold(sex, male_value, female_value):
    """Returns the sex-specific threshold, or None if sex is unknown —
    letting downstream is_elevated/is_below/is_outside_range checks
    safely resolve to False rather than silently defaulting to female threshold if sex != male"""
    if sex == "male":
        return male_value
    elif sex == "female":
        return female_value
    return None

def is_elevated(value, threshold):
    """Returns True if value is present and >= threshold."""
    return value is not None and threshold is not None and value >= threshold

def is_above(value, threshold):
    """Returns True if value is present and strictly greater than threshold."""
    return value is not None and threshold is not None and value > threshold

def is_below(value, threshold):
    """Returns True if value is present and < threshold."""
    return value is not None and threshold is not None and value < threshold

def is_low(value, threshold):
    """Returns True if value is present and <= threshold."""
    return value is not None and threshold is not None and value <= threshold

def is_outside_range(value, low, high):
    """Returns True if value is present and outside [low, high]."""
    return (
        value is not None
        and low is not None
        and high is not None
        and (value < low or value > high)
    )

# def is_classification_elevated(genetic_risk_classification):
#     """Returns True if the genetic risk classification counts as elevated
#     (Elevated or Moderately Elevated)"""
#     return genetic_risk_classification in ELEVATED_CLASSIFICATIONS

# added
def note(triggered, label, value, condition_bool, threshold_desc=""):
    """Appends 'label=value (threshold_desc)' to triggered if condition_bool is True. Returns
    condition_bool unchanged, so this can be dropped in-place around any
    is_elevated/is_below/is_above/is_outside_range call. threshold_desc should describe the
    threshold that was crossed, e.g. ">=130" or "<40" or "outside [150,400]"."""
    if condition_bool:
        entry = f"{label}={value}"
        if threshold_desc:
            entry += f" ({threshold_desc})"
        triggered.append(entry)
    return condition_bool

# for conditions to be run be default without any gene/prs check 
NO_GENE_GATE = "NO_GENE_GATE"

# Gene requirements per condition, from the slide deck
CONDITION_GENES = {
    "Alzheimer's Disease": ["APOE"],
    # "Asthma": ["IL13"],
    "Atopic Dermatitis": ["C11orf32/LRRC32", "TSLP", "IL13", "IL4R"],
    "Eczema": [], #no gene listed, only PRS gate
    "COPD": [],  # no gene listed, only PRS gate
    "Hyperthyroidism": ["CTLA4", "PTPN22"],
    "Hypothyroidism": ["CTLA4"],
    "Inflammatory Bowel Disease": ["NOD2", "ATG16L1", "IL23R", "MYC/POU5F1B", "IL27"],
    "NAFLD": ["PNPLA3", "TM6SF2"],
    # "Osteoarthritis": ["COL11A1", "SMAD3"],
    # "Parkinson's Disease": ["LRRK2", "GBA"],
    "Psoriasis": [],  # no gene listed, only PRS gate
    "Rheumatoid Arthritis": ["PTPN22", "STAT4", "TNFAIP3", "CTLA4"],
    # "Rhinitis": ["TSLP"],
    "Type 2 Diabetes": ["TCF7L2", "SLC30A8"],#------------------------------------------------------------------
    "Maturity-Onset Diabetes of the Young": ["HNF1A"],
    # "Multiple Endocrine Neoplasia Type 4": ["CDKN1B"],
    "Muscular Dystrophy": ["DMD", "FKRP", "LMNA", "SGCA", "SGCB"],
    "Gilbert Syndrome": ["UGT1A1"],
    "Hereditary Hemochromatosis": ["HFE"],
    "Celiac Disease": ["HLA-DQA1/HLA-DQB1"],
    "Familial Hypercholesterolemia": ["LDLR", "APOB", "PCSK9"],
    "Cardiomyopathy": [
        "ACTC1", "BAG3", "DES", "FLNC", "LMNA", "MYBPC3", "MYH7", "MYL2",
        "MYL3", "PRKAG2", "RBM20", "TNNC1", "TNNI3", "TNNT2", "TPM1", "TTN",
        "PLN", "DSC2", "DSG2", "DSP", "PKP2", "TMEM43",
    ],
    "Coronary Artery Disease": [
        "ABCA1", "ABCG5", "ABCG8", "APOA1", "APOA5", "APOB", "APOC2", "APOC3",
        "CREB3L3", "GPIHBP1", "LDLR", "LDLRAP1", "LMF1", "LPL", "PCSK9",
    ],
    "Familial Hypertriglyceridemia": ["CREB3L3", "LMF1", "LPL"],
    "HDL Deficiency": ["ABCA1", "LCAT", "APOA1"],
    # "Elevated Apolipoprotein B": ["APOB"],
    "Li-Fraumeni Syndrome": ["TP53"],
    "PTEN Hamartoma Tumor Syndrome": ["PTEN"],
    "Prostate Cancer": [
        "HOXB13", "ELAC2", "CDH1", "CHEK2", "BRCA2", "PALB2", "ATM", "PTEN",
        "TP53", "MSH6", "XPC", "AR", "BRCA1"
    ],
}

# which report "Domain" each condition is grouped under - kept live (unlike the rest of
# the CLI-only tail below) so the dispatcher can tag a condition's result with its
# domain the same way the CLI's results assembly did, without needing any of that
# dead weight.
CONDITION_DOMAINS = {
    "Alzheimer's Disease": "Brain Health",
    # "Asthma": "Respiratory Health",
    "Atopic Dermatitis": "Skin Health",
    "Eczema": "Skin Health",
    "COPD": "Respiratory Health",
    "Hyperthyroidism": "Endocrine Health",
    "Hypothyroidism": "Endocrine Health",
    "Inflammatory Bowel Disease": "Digestive Health",
    "NAFLD": "Digestive Health",
    # "Osteoarthritis": "Joints Health",
    # "Parkinson's Disease": "Brain Health",
    "Psoriasis": "Skin Health",
    "Rheumatoid Arthritis": "Immune Health",
    # "Rhinitis": "Respiratory Health",
    "Type 2 Diabetes": "Endocrine Health",
    "Maturity-Onset Diabetes of the Young": "Endocrine Health",
    # "Multiple Endocrine Neoplasia Type 4": "Endocrine Health",
    "Muscular Dystrophy": "Muscle Health",
    "Gilbert Syndrome": "Kidney & Metabolic Health",
    "Hereditary Hemochromatosis": "Digestive Health",
    "Celiac Disease": "Immune Health",
    "Familial Hypercholesterolemia": "Cardiac Health",
    "Cardiomyopathy": "Cardiac Health",
    "Coronary Artery Disease": "Cardiac Health",
    "Familial Hypertriglyceridemia": "Cardiac Health",
    "HDL Deficiency": "Cardiac Health",
    "Li-Fraumeni Syndrome": "Cancer",
    "PTEN Hamartoma Tumor Syndrome": "Cancer",
    "Prostate Cancer": "Cancer",
}

# only conditions filtered through this will have their functions run on patient data block
def has_flagged_gene(genetics, condition_name):
    """Returns True if:
    - the condition is exempt from the gene gate (NO_GENE_GATE), OR
    - at least one of its required genes is flagged (monogenic hit), OR
    - the condition has a polygenic PRS entry AND that PRS category is
      elevated/moderately elevated.
    Returns False (blocked) otherwise. genetics=None (not just an empty
    dict) is treated as "not available" everywhere below, not as a tested
    negative - callers wanting that distinction should check
    `genetics is None` before calling this."""
    genetics = genetics or {}
    required_genes = CONDITION_GENES.get(condition_name, [])

    if required_genes == NO_GENE_GATE:
        return True

    prs_key = None
    for key, mapped_name in PRS_CONDITION_KEY_MAP.items():
        if mapped_name == condition_name:
            prs_key = key
            break
    if prs_key is not None and prs_key in genetics.get("prs_elevated_conditions", set()):
        return True


    if not required_genes:
        return False
    flagged = genetics.get("flagged_genes", [])
    return any(gene in flagged for gene in required_genes)

# added
def get_gene_trigger(genetics, condition_name):
    """Returns a string describing what passed the gene/PRS gate for this condition:
    matched flagged gene(s) and/or 'PRS:<key>' if the condition's PRS entry is elevated.
    Returns '' for NO_GENE_GATE conditions or when nothing matched (shouldn't happen if
    has_flagged_gene already returned True)."""
    genetics = genetics or {}
    required_genes = CONDITION_GENES.get(condition_name, [])
    parts = []

    if required_genes != NO_GENE_GATE and required_genes:
        flagged = genetics.get("flagged_genes", [])
        parts.extend(gene for gene in required_genes if gene in flagged)

    # Triggering Genetics renamed as DNA Marker(s), now onlt includes gene flag, PRS classification in Triggering PRS, elevation info already being added in report.
    # prs_key = None
    # for key, mapped_name in PRS_CONDITION_KEY_MAP.items():
    #     if mapped_name == condition_name:
    #         prs_key = key
    #         break
    # if prs_key is not None and prs_key in genetics.get("prs_elevated_conditions", set()):
    #     parts.append(f"PRS:{prs_key}")

    return ", ".join(parts)

def get_monogenic_variant(genetics, condition_name):
    genetics = genetics or {}
    lookup = genetics.get("monogenic_lookup", {})
    required_genes = CONDITION_GENES.get(condition_name, [])
    if required_genes == NO_GENE_GATE or not required_genes:
        return "", ""
    key_lc = condition_name.lower()
    rsids, patho = [], []
    for gene in required_genes:
        for hit in lookup.get((gene, key_lc), []):
            rsids.append(hit[0])
            patho.append(hit[1])
    return ", ".join(rsids), ", ".join(patho)

def get_apoe_status(genetics, condition_name):
    genetics = genetics or {}
    required_genes = CONDITION_GENES.get(condition_name, [])
    if required_genes == NO_GENE_GATE or "APOE" not in required_genes:
        return ""
    if "APOE" not in genetics.get("flagged_genes", []):
        return ""
    return genetics.get("apoe_status", "")

def get_hla_allele(genetics, condition_name):
    genetics = genetics or {}
    required_genes = CONDITION_GENES.get(condition_name, [])
    if required_genes == NO_GENE_GATE or "HLA-DQA1/HLA-DQB1" not in required_genes:
        return ""
    if "HLA-DQA1/HLA-DQB1" not in genetics.get("flagged_genes", []):
        return ""
    return genetics.get("hla_allele", "")

def is_gene_flagged(genetics, condition_name):
    genetics = genetics or {}
    required_genes = CONDITION_GENES.get(condition_name, [])
    if not required_genes or required_genes == NO_GENE_GATE:
        return False
    flagged = genetics.get("flagged_genes", [])
    return any(gene in flagged for gene in required_genes)

def is_gene_acmg(genetics, condition_name):
    genetics = genetics or {}
    required_genes = CONDITION_GENES.get(condition_name, [])
    if not required_genes or required_genes == NO_GENE_GATE:
        return False
    acmg = genetics.get("acmg_genes", [])
    return any(gene in acmg for gene in required_genes)

def get_prs_category(genetics, condition_name):
    genetics = genetics or {}
    prs_key = None
    for key, mapped_name in PRS_CONDITION_KEY_MAP.items():
        if mapped_name == condition_name:
            prs_key = key
            break
    if prs_key is None:
        return None  # no PRS entry for this condition at all
    return genetics.get("prs_categories", {}).get(prs_key, "")

# =============================================================================
# Shared snapshot evaluation logic

PATTERN_CATEGORIES_INDICATING_TRIGGER = {
    "Early Pattern", "Significant Pattern"
}

ELEVATED_PRS = {"elevated", "moderately_elevated", "moderately elevated"}
REDUCED_TYPICAL_PRS = {"typical", "reduced", "moderately_reduced", "moderately reduced"}

case1 = "To discuss with General Practitioner"
case2 = "Worth acting on for prevention"
case3 = "Typical - nothing to act on"

CASE_PRIORITY = {case1: 0, case2: 1, case3: 2} #for case where both prs and gene are there, a lower priority doesn't suppress an upper one in final csv.

def get_snapshot_category(condition_name, result_category, genetics):
    if result_category == INSUFFICIENT_DATA:
        return "Not enough data to assess - testing recommended"

    gene_flagged = is_gene_flagged(genetics, condition_name)
    gene_acmg = is_gene_acmg(genetics, condition_name)
    gene_ran_pattern = result_category in PATTERN_CATEGORIES_INDICATING_TRIGGER

    if gene_flagged:
        gene_case = case1 if (gene_acmg or gene_ran_pattern) else case2
    else:
        gene_case = None

    prs_category = get_prs_category(genetics, condition_name)
    if prs_category in ELEVATED_PRS:
        prs_ran_pattern = gene_ran_pattern #same check for result_category
        prs_case = case1 if prs_ran_pattern else case2
    elif prs_category in REDUCED_TYPICAL_PRS:
        prs_case = case3
    else:
        prs_case = None  # no PRS entry, or unrecognized category

    # --- combine per priority rule: case1 > case2 > case3 ---
    candidates = [c for c in (gene_case, prs_case) if c is not None]
    if candidates:
        return min(candidates, key=lambda c: CASE_PRIORITY[c])

    # no gene, no PRS entry, the gate blocked entirely
    if result_category == "Typical":
        return case3

    return "" # fallback for error check (no clean answer)

PRS_DISPLAY_VALUE = {
    "elevated": "Polygenic Risk Classification",
    "moderately_elevated": "Polygenic Risk Classification",
    "moderately elevated": "Polygenic Risk Classification",
}

def get_triggering_prs(condition_name, genetics):
    """Returns the display value for Triggering PRS if this condition's PRS
    category is elevated/moderately elevated, else ''."""
    prs_category = get_prs_category(genetics, condition_name)
    return PRS_DISPLAY_VALUE.get(prs_category, "")


# =============================================================================
# Per-result enrichment (Domain/Triggering PRS/Snapshot Category/RSID_m/
# Pathogenicity_m/APOE_Status/Elevation_p/Threshold_p/All Blood Marker(s)) -
# this used to only happen in the CLI's now-dead post-processing loop. Kept
# here, applied via with_enrichment() below, so every evaluate_*'s OWN
# return value is already the fully enriched result - no separate step
# needed wherever evaluate_* gets called from (dispatcher or otherwise).
# =============================================================================

def enrich_result(entry, genetics):
    condition = entry["Condition"]
    entry["Domain"] = CONDITION_DOMAINS.get(condition, "")
    entry["Triggering PRS"] = get_triggering_prs(condition, genetics)
    entry["Snapshot Category"] = get_snapshot_category(condition, entry["Category"], genetics)
    entry["RSID_m"], entry["Pathogenicity_m"] = get_monogenic_variant(genetics, condition)
    entry["APOE_Status"] = get_apoe_status(genetics, condition)
    entry["HLA_Allele"] = get_hla_allele(genetics, condition)

    elevation = get_prs_category(genetics, condition) or ""
    entry["Elevation_p"] = elevation
    entry["Threshold_p"] = (
        ">99 percentile" if elevation == "elevated"
        else ">97 percentile" if elevation in {"moderately elevated", "moderately_elevated"}
        else ""
    )

    entry["All Blood Marker(s)"] = entry.get("Blood Marker(s)", "")
    if entry.get("Category") not in PATTERN_CATEGORIES_INDICATING_TRIGGER:
        entry["Blood Marker(s)"] = ""
    return entry


def with_enrichment(fn):
    """Wraps an evaluate_* function so its return value is already enriched.
    Binds against fn's real signature (not just kwargs.get("genetics")) so
    this works whether a caller passes genetics positionally or by keyword -
    callers (the dispatcher, tests, CLI-style calls) can use either, and
    nothing calling evaluate_* needs to know enrich_result exists at all."""
    sig = inspect.signature(fn)

    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        bound = sig.bind(*args, **kwargs)
        genetics = bound.arguments.get("genetics")
        return [enrich_result(entry, genetics) for entry in fn(*args, **kwargs)]
    return wrapped


# def resolve_typical_label(condition_name, category, genetics):
#     """If the function landed on 'Typical' (gate passed but nothing triggered clinically),
#     relabel it to show what actually passed the gate: PRS elevation level, or gene flag."""
#     if category != "Typical":
#         return category

#     prs_category = get_prs_category(genetics, condition_name)
#     if prs_category == "elevated":
#         return "Elevated (PRS)"
#     if prs_category in {"moderately_elevated", "moderately elevated"}:
#         return "Moderately Elevated (PRS)"
#     if is_gene_flagged(genetics, condition_name):
#         return "Gene Flagged"
#     return category
# 
# =============================================================================

# Family history and symptoms of respective conditions, are mentioned in the deck for ref

GENE_NOT_FOUND = "Typical"
INSUFFICIENT_DATA = "Insufficient Data"  # genetics unavailable (not just gate-negative), or none of this condition's labs were provided

# Alzheimer's Disease
def evaluate_alzheimers(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Alzheimer's Disease", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Alzheimer's Disease"):
        return [{"Condition": "Alzheimer's Disease", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    # hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)
    hdl_min = sex_based_threshold(sex, 60, 70) #absolute assignment as per the deck

    total_c = labs.get("total_cholesterol")
    ldl_c = labs.get("ldl_c")
    triglycerides = labs.get("triglycerides")
    hdl_c = labs.get("hdl_c")
    crp = labs.get("crp")

    triggered = []
    # --- LDL ---
    ldl_likely = is_elevated(ldl_c, LDL_C_ELEVATED)
    ldl_early = is_elevated(ldl_c, LDL_C_AT_RISK)
    if ldl_likely:
        note(triggered, "ldl_c", ldl_c, True, f">={LDL_C_ELEVATED}")
    elif ldl_early:
        note(triggered, "ldl_c", ldl_c, True, f">={LDL_C_AT_RISK}")

    # --- Total cholesterol ---
    total_c_likely = is_elevated(total_c, TOTAL_CHOLESTEROL_ELEVATED)
    total_c_early = is_elevated(total_c, TOTAL_CHOLESTEROL_AT_RISK)
    if total_c_likely:
        note(triggered, "total_cholesterol", total_c, True, f">={TOTAL_CHOLESTEROL_ELEVATED}")
    elif total_c_early:
        note(triggered, "total_cholesterol", total_c, True, f">={TOTAL_CHOLESTEROL_AT_RISK}")

    # --- Triglycerides ---
    tg_likely = is_elevated(triglycerides, TRIGLYCERIDES_ELEVATED)
    tg_early = is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK)
    if tg_likely:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_ELEVATED}")
    elif tg_early:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_AT_RISK}")

    # --- HDL (single threshold, no early/likely tiers) ---
    hdl_low = note(triggered, "hdl_c", hdl_c, is_below(hdl_c, hdl_min), f"<{hdl_min}")

    dyslipidemia_flag_early = ldl_early or total_c_early or tg_early or hdl_low
    dyslipidemia_flag_likely = ldl_likely or total_c_likely or tg_likely or hdl_low

    crp_early = is_elevated(crp, CRP_AT_RISK) #and not is_above(crp, CRP_ELEVATED)
    crp_likely = is_above(crp, CRP_ELEVATED)

    if crp_likely:
        note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
    elif crp_early:
        note(triggered, "crp", crp, True, f">={CRP_AT_RISK}")

    if (
        (dyslipidemia_flag_likely and crp_likely) or
        (dyslipidemia_flag_likely and crp_likely and (family_history or past_history or symptoms))
    ):
        category = "Significant Pattern"
    elif dyslipidemia_flag_early and crp_early:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Alzheimer's Disease", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Alzheimer's Disease"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_alzheimers_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_alzheimers',
        'description': "Evaluate Alzheimer's Disease susceptibility from lipid panel and crp, gated on APOE or PRS alzheimers genetics.",
        'parameters': {
            'type': 'object',
            'properties': {
                'total_cholesterol': {
                    'type': ['number', 'null'],
                    'description': 'total_cholesterol, mg/dL. Null if not tested.',
                },
                'ldl_c': {
                    'type': ['number', 'null'],
                    'description': 'ldl_c, mg/dL. Null if not tested.',
                },
                'triglycerides': {
                    'type': ['number', 'null'],
                    'description': 'triglycerides, mg/dL. Null if not tested.',
                },
                'hdl_c': {
                    'type': ['number', 'null'],
                    'description': 'hdl_c, mg/dL. Null if not tested.',
                },
                'crp': {
                    'type': ['number', 'null'],
                    'description': 'crp, mg/L. Null if not tested.',
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES["Alzheimer's Disease"],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'alzheimers_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'alzheimers' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get("Alzheimer's Disease", {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get("Alzheimer's Disease", {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get("Alzheimer's Disease", {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['total_cholesterol', 'ldl_c', 'triglycerides', 'hdl_c', 'crp', 'sex', 'flagged_genes', 'alzheimers_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Asthma
# def evaluate_asthma(labs, patient, genetics, family_history, past_history, symptoms=False):
#     if genetics is None:
#         return [{"Condition": "Asthma", "Category": INSUFFICIENT_DATA}]
#     if not has_flagged_gene(genetics, "Asthma"):
#         return [{"Condition": "Asthma", "Category": GENE_NOT_FOUND}]

#     eosinophils = labs.get("eosinophils")
#     neutrophils = labs.get("neutrophils")
#     triggered = []

#     eosinophil_likely = is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
#     eosinophil_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
#     if eosinophil_likely:
#         note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_UPPER_NORMAL}")
#     elif eosinophil_early:
#         note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")

#     neutrophil_flag = note(triggered, "neutrophils", neutrophils, is_above(neutrophils, NEUTROPHILS_RANGE[0]), f">{NEUTROPHILS_RANGE[0]}")

#     if (eosinophil_likely and neutrophil_flag) or (eosinophil_likely and neutrophil_flag and (symptoms or family_history or past_history)):
#         category = "Significant Pattern"
#     elif eosinophil_early:
#         category = "Early Pattern"
#     else:
#         category = "Elevated susceptibility"

#     return [{
#         "Condition": "Asthma", "Category": category,
#         "DNA Marker(s)": get_gene_trigger(genetics, "Asthma"),
#         "Blood Marker(s)": ", ".join(triggered),
#     }]
# evaluate_asthma_schema = {
#     'type': 'function',
#     'function': {
#         'name': 'evaluate_asthma',
#         'description': 'Evaluate Asthma susceptibility from eosinophil/neutrophil blood counts, gated on IL13 or PRS asthma genetics.',
#         'parameters': {
#             'type': 'object',
#             'properties': {
#                 'eosinophils': {
#                     'type': ['number', 'null'],
#                     'description': 'eosinophils, x10^9/L. Null if not tested.',
#                 },
#                 'neutrophils': {
#                     'type': ['number', 'null'],
#                     'description': 'neutrophils, x10^9/L. Null if not tested.',
#                 },
#                 'flagged_genes': {
#                     'type': ['array', 'null'],
#                     'items': {
#                         'type': 'string',
#                         'enum': CONDITION_GENES['Asthma'],
#                     },
#                     'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
#                 },
#                 'asthma_prs_elevated': {
#                     'type': ['boolean', 'null'],
#                     'description': "True if this patient's PRS 'asthma' category is elevated or moderately_elevated; null if PRS unavailable.",
#                 },
#                 'family_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Asthma', {}).get('family_history', []),
#                     },
#                     'description': 'Which of these the patient reports in family history.',
#                 },
#                 'past_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Asthma', {}).get('past_history', []),
#                     },
#                     'description': "Which of these are in the patient's own past medical history.",
#                 },
#                 'symptoms': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Asthma', {}).get('symptoms', []),
#                     },
#                     'description': 'Which of these symptoms the patient reports.',
#                 },
#             },
#             'required': ['eosinophils', 'neutrophils', 'flagged_genes', 'asthma_prs_elevated', 'family_history', 'past_history', 'symptoms'],
#             'additionalProperties': False,
#         },
#         'strict': True,
#     },
# }

# Atopic Dermatitis
def evaluate_atopic_dermatitis(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Atopic Dermatitis", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Atopic Dermatitis"):
        return [{"Condition": "Atopic Dermatitis", "Category": GENE_NOT_FOUND}]

    eosinophils = labs.get("eosinophils")
    ige = labs.get("ige")
    triggered = []

    eosinophil_likely = is_above(eosinophils, EOSINOPHILS_LIKELY)
    eosinophil_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_LIKELY)
    ige_likely = is_elevated(ige, IGE_LIKELY)
    ige_early = is_elevated(ige, IGE_EARLY)
    
    if eosinophil_likely:
        note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_LIKELY}")
    elif eosinophil_early:
        note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")

    if ige_likely:
        note(triggered, "ige", ige, True, f">={IGE_LIKELY}")
    elif ige_early:
        note(triggered, "ige", ige, True, f">={IGE_EARLY}")

    if (eosinophil_likely and ige_likely) or ((eosinophil_likely and ige_likely) and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif eosinophil_early and ige_early:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Atopic Dermatitis", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Atopic Dermatitis"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_atopic_dermatitis_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_atopic_dermatitis',
        'description': 'Evaluate Atopic Dermatitis susceptibility from eosinophil count and Total IgE, gated on C11orf32/LRRC32, TSLP, IL13, IL4R, or PRS atopicdermatitis genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'eosinophils': {
                    'type': ['number', 'null'],
                    'description': 'eosinophils, x10^9/L. Null if not tested.',
                },
                'ige': {
                    'type': ['number', 'null'],
                    'description': 'Total IgE, IU/mL. Null if not tested.',
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Atopic Dermatitis'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'atopicdermatitis_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'atopicdermatitis' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Atopic Dermatitis', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Atopic Dermatitis', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Atopic Dermatitis', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['eosinophils', 'ige', 'flagged_genes', 'atopicdermatitis_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Eczema
def evaluate_eczema(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Eczema", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Eczema"):
        return [{"Condition": "Eczema", "Category": GENE_NOT_FOUND}]

    eosinophils = labs.get("eosinophils")
    ige = labs.get("ige")
    triggered = []

    eosinophil_likely = is_above(eosinophils, EOSINOPHILS_LIKELY)
    eosinophil_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_LIKELY)
    ige_likely = is_elevated(ige, IGE_LIKELY)
    ige_early = is_elevated(ige, IGE_EARLY)
    
    if eosinophil_likely:
        note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_LIKELY}")
    elif eosinophil_early:
        note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")
    
    if ige_likely:
        note(triggered, "ige", ige, True, f">={IGE_LIKELY}")
    elif ige_early:
        note(triggered, "ige", ige, True, f">={IGE_EARLY}")

    if (eosinophil_likely and ige_likely) or ((eosinophil_likely and ige_likely) and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif eosinophil_early and ige_early:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Eczema", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Eczema"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_eczema_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_eczema',
        'description': 'Evaluate Eczema susceptibility from eosinophil count and Total IgE, gated on PRS eczema only (no required gene).',
        'parameters': {
            'type': 'object',
            'properties': {
                'eosinophils': {
                    'type': ['number', 'null'],
                    'description': 'eosinophils, x10^9/L. Null if not tested.',
                },
                'ige': {
                    'type': ['number', 'null'],
                    'description': 'Total IgE, IU/mL. Null if not tested.',
                },
                'eczema_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'eczema' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Eczema', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Eczema', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Eczema', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['eosinophils', 'ige', 'eczema_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Chronic Obstruction of Pulmonary Disorder​
def evaluate_copd(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "COPD", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "COPD"):
        return [{"Condition": "COPD", "Category": GENE_NOT_FOUND}]

    crp = labs.get("crp")
    fibrinogen = labs.get("fibrinogen")
    eosinophils = labs.get("eosinophils")
    neutrophils = labs.get("neutrophils")

    triggered = []
    crp_elevated_flag = is_above(crp, CRP_ELEVATED)
    crp_at_risk_flag = is_elevated(crp, CRP_AT_RISK_COPD) #and not is_above(crp, CRP_ELEVATED)
    if crp_elevated_flag:
        note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
    elif crp_at_risk_flag:
        note(triggered, "crp", crp, True, f">={CRP_AT_RISK_COPD}")

    eosinophil_flag_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_LIKELY)
    eosinophil_flag_likely = is_above(eosinophils, EOSINOPHILS_LIKELY)
    if eosinophil_flag_likely:
        note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_LIKELY}")
    elif eosinophil_flag_early:
        note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")

    fibrinogen_flag = note(triggered, "fibrinogen", fibrinogen, is_above(fibrinogen, FIBRINOGEN_RANGE[0]), f">={FIBRINOGEN_RANGE[0]}")
    neutrophil_flag = note(triggered, "neutrophils", neutrophils, is_above(neutrophils, NEUTROPHILS_RANGE[0]), f">{NEUTROPHILS_RANGE[0]}")

    if (eosinophil_flag_likely and neutrophil_flag and fibrinogen_flag and crp_elevated_flag) or (eosinophil_flag_likely and neutrophil_flag and fibrinogen_flag and crp_elevated_flag and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif (crp_at_risk_flag and eosinophil_flag_early):
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "COPD", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "COPD"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_copd_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_copd',
        'description': 'Evaluate COPD susceptibility from crp, fibrinogen, eosinophil and neutrophil counts, gated on PRS copd only (no specific gene).',
        'parameters': {
            'type': 'object',
            'properties': {
                'crp': {
                    'type': ['number', 'null'],
                    'description': 'crp, mg/L. Null if not tested.',
                },
                'fibrinogen': {
                    'type': ['number', 'null'],
                    'description': 'fibrinogen, mg/dL. Null if not tested.',
                },
                'eosinophils': {
                    'type': ['number', 'null'],
                    'description': 'eosinophils, x10^9/L. Null if not tested.',
                },
                'neutrophils': {
                    'type': ['number', 'null'],
                    'description': 'neutrophils, x10^9/L. Null if not tested.',
                },
                'copd_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'copd' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('COPD', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('COPD', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('COPD', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['crp', 'fibrinogen', 'eosinophils', 'neutrophils', 'copd_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Hyperthyroidism
def evaluate_hyperthyroidism(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Hyperthyroidism", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Hyperthyroidism"):
        return [{"Condition": "Hyperthyroidism", "Category": GENE_NOT_FOUND}]

    tsh = labs.get("tsh")
    free_t3 = labs.get("free_t3")
    free_t4 = labs.get("free_t4")
    tpoab = labs.get("tpoab")

    triggered = []
    # tsh_borederline_suppressed = is_below(tsh, TSH_BORDERLINE)
    tsh_suppressed = is_below(tsh, TSH_SUPPRESSED)
    if tsh_suppressed:
        note(triggered, "tsh", tsh, True, f"<{TSH_SUPPRESSED}")
    # elif tsh_borederline_suppressed:
    #     note(triggered, "tsh", tsh, True, f"<{TSH_BORDERLINE}")

    free_t3_flag = is_above(free_t3, FREE_T3_RANGE[1])
    free_t3_normal = free_t3 is not None and FREE_T3_RANGE[0] <= free_t3 <= FREE_T3_RANGE[1]
    free_t4_flag = is_above(free_t4, FREE_T4_RANGE[1])
    free_t4_normal = free_t4 is not None and FREE_T4_RANGE[0] <= free_t4 <= FREE_T4_RANGE[1]

    hormone_flag = free_t3_flag or free_t4_flag
    hormone_normal = free_t3_normal and free_t4_normal

    tpoab_flag = note(triggered, "tpoab", tpoab, is_elevated(tpoab, TPOAB_UPPER_NORMAL), f">={TPOAB_UPPER_NORMAL}")

    # only elevated T3/T4 reported, normal ones only used in logic and not reported.
    if free_t3_flag:
        note(triggered, "free_t3", free_t3, True, f">{FREE_T3_RANGE[1]}")
    # elif free_t3_normal and (tsh_suppressed or tpoab_flag):
    #     note(triggered, "free_t3", free_t3, True, f"({FREE_T3_RANGE[0]}-{FREE_T3_RANGE[1]})")

    if free_t4_flag:
        note(triggered, "free_t4", free_t4, True, f">{FREE_T4_RANGE[1]}")
    # elif free_t4_normal and (tsh_suppressed or tpoab_flag):
    #     note(triggered, "free_t4", free_t4, True, f"({FREE_T4_RANGE[0]}-{FREE_T4_RANGE[1]})")

    if (tsh_suppressed and hormone_flag and tpoab_flag) or (tsh_suppressed and hormone_flag and tpoab_flag and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif tsh_suppressed and (hormone_normal or tpoab_flag):
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Hyperthyroidism", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Hyperthyroidism"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_hyperthyroidism_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_hyperthyroidism',
        'description': 'Evaluate Hyperthyroidism susceptibility from tsh, free_t3, free_t4 and tpoab, gated on CTLA4, PTPN22, or PRS hyperthyroidism genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'tsh': {
                    'type': ['number', 'null'],
                    'description': 'tsh, uIU/mL. Null if not tested.',
                },
                'free_t3': {
                    'type': ['number', 'null'],
                    'description': 'free_t3, pg/mL. Null if not tested.',
                },
                'free_t4': {
                    'type': ['number', 'null'],
                    'description': 'free_t4, ng/dL. Null if not tested.',
                },
                'tpoab': {
                    'type': ['number', 'null'],
                    'description': 'tpoab, IU/mL. Null if not tested.',
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Hyperthyroidism'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'hyperthyroidism_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'hyperthyroidism' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hyperthyroidism', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hyperthyroidism', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hyperthyroidism', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['tsh', 'free_t3', 'free_t4', 'tpoab', 'flagged_genes', 'hyperthyroidism_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Hypothyroidism
def evaluate_hypothyroidism(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Hypothyroidism", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Hypothyroidism"):
        return [{"Condition": "Hypothyroidism", "Category": GENE_NOT_FOUND}]

    tsh = labs.get("tsh")
    # free_t3 = labs.get("free_t3")
    free_t4 = labs.get("free_t4")
    tpoab = labs.get("tpoab")

    triggered = []
    tsh_high = is_above(tsh, TSH_ELEVATED_RANGE[1])
    # tsh_mild = tsh is not None and TSH_ELEVATED_RANGE[0] <= tsh <= TSH_ELEVATED_RANGE[1] #this had an edge case
    tsh_mild = is_elevated(tsh, TSH_ELEVATED_RANGE[0])
    if tsh_high:
        note(triggered, "tsh", tsh, True, f">{TSH_ELEVATED_RANGE[1]}")
    elif tsh_mild:
        note(triggered, "tsh", tsh, True, f"[{TSH_ELEVATED_RANGE[0]}-{TSH_ELEVATED_RANGE[1]}]")
    
    # free_t3_flag = note(triggered, "free_t3", free_t3, is_below(free_t3, FREE_T3_ELEVATED), f"<{FREE_T3_ELEVATED}") #unused in triggers
    free_t4_flag = is_below(free_t4, FREE_T4_RANGE[0])
    free_t4_normal = free_t4 is not None and FREE_T4_RANGE[0] <= free_t4 <= FREE_T4_RANGE[1]
    tpoab_flag = note(triggered, "tpoab", tpoab, is_elevated(tpoab, TPOAB_UPPER_NORMAL), f">={TPOAB_UPPER_NORMAL}")

    # Only low T4 reported, normal T4 only used in logic, not reported 
    if free_t4_flag:
        note(triggered, "free_t4", free_t4, True, f"<{FREE_T4_RANGE[0]}")
    # elif free_t4_normal and (tsh_high or tsh_mild or tpoab_flag):
    #     note(triggered, "free_t4", free_t4, True, f"({FREE_T4_RANGE[0]}-{FREE_T4_RANGE[1]})")

    if (tsh_high and free_t4_flag and tpoab_flag) or (tsh_high and free_t4_flag and tpoab_flag and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    # tsh_mild and tpoab_flag case covers early pattern w/o user context
    elif tsh_mild and (tpoab_flag or free_t4_normal):
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Hypothyroidism", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Hypothyroidism"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_hypothyroidism_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_hypothyroidism',
        'description': 'Evaluate Hypothyroidism susceptibility from tsh, free_t4 and tpoab, gated on CTLA4 or PRS hypothyroidism genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'tsh': {
                    'type': ['number', 'null'],
                    'description': 'tsh, uIU/mL. Null if not tested.',
                },
                'free_t4': {
                    'type': ['number', 'null'],
                    'description': 'free_t4, ng/dL. Null if not tested.',
                },
                'tpoab': {
                    'type': ['number', 'null'],
                    'description': 'tpoab, IU/mL. Null if not tested.',
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Hypothyroidism'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'hypothyroidism_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'hypothyroidism' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hypothyroidism', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hypothyroidism', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hypothyroidism', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['tsh', 'free_t4', 'tpoab', 'flagged_genes', 'hypothyroidism_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Inflammatory Bowel Disease
def evaluate_ibd(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Inflammatory Bowel Disease", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Inflammatory Bowel Disease"):
        return [{"Condition": "Inflammatory Bowel Disease", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    hb_min = sex_based_threshold(sex, HEMOGLOBIN_MALE_RANGE[0], HEMOGLOBIN_FEMALE_RANGE[0])

    crp = labs.get("crp")
    hemoglobin = labs.get("hemoglobin")
    albumin = labs.get("albumin")

    triggered = []
    crp_elevated_flag = is_above(crp, CRP_ELEVATED)
    crp_at_risk_flag = is_elevated(crp, CRP_AT_RISK_COPD) #and not is_above(crp, CRP_ELEVATED)
    if crp_elevated_flag:
        note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
    elif crp_at_risk_flag:
        note(triggered, "crp", crp, True, f">={CRP_AT_RISK_COPD}")
    anemia_flag = note(triggered, "hemoglobin", hemoglobin, is_below(hemoglobin, hb_min), f"<{hb_min}")
    albumin_flag = note(triggered, "albumin", albumin, is_below(albumin, ALBUMIN_RANGE[0]), f"<{ALBUMIN_RANGE[0]}")
    if (crp_elevated_flag and anemia_flag and albumin_flag) or (crp_elevated_flag and anemia_flag and albumin_flag and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif crp_at_risk_flag:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Inflammatory Bowel Disease", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Inflammatory Bowel Disease"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_ibd_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_ibd',
        'description': 'Evaluate Inflammatory Bowel Disease susceptibility from crp, hemoglobin and albumin, gated on NOD2, ATG16L1, IL23R, MYC/POU5F1B, IL27, or PRS ibd genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'crp': {
                    'type': ['number', 'null'],
                    'description': 'crp, mg/L. Null if not tested.',
                },
                'hemoglobin': {
                    'type': ['number', 'null'],
                    'description': 'hemoglobin, g/dL. Null if not tested.',
                },
                'albumin': {
                    'type': ['number', 'null'],
                    'description': 'albumin, g/dL. Null if not tested.',
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Inflammatory Bowel Disease'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'ibd_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'ibd' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Inflammatory Bowel Disease', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Inflammatory Bowel Disease', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Inflammatory Bowel Disease', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['crp', 'hemoglobin', 'albumin', 'sex', 'flagged_genes', 'ibd_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Non alcoholic fatty liver disease
def evaluate_nafld(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "NAFLD", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "NAFLD"):
        return [{"Condition": "NAFLD", "Category": GENE_NOT_FOUND}]
    
    age = patient.get("age")
    fib4 = labs.get("fib4")
    ast = labs.get("ast")
    alt = labs.get("alt")
    ggt = labs.get("ggt")
    platelets = labs.get("platelets")
    ast_alt_ratio = labs.get("ast/alt")

    triggered = []
    fib4_elevated = is_elevated(fib4, FIB4_ELEVATED)
    fib4_at_risk = is_elevated(fib4, FIB4_AT_RISK[0]) if age < 65 else is_elevated(fib4, FIB4_AT_RISK[1])
    
    if fib4_elevated:
        note(triggered, "fib4", fib4, True, f">={FIB4_ELEVATED}")
    elif fib4_at_risk:
        threshold = FIB4_AT_RISK[0] if age < 65 else FIB4_AT_RISK[1]
        note(triggered, "fib4", fib4, True, f">={threshold}")
        
    ast_flag = note(triggered, "ast", ast, is_above(ast, AST_RANGE[1]), f">{AST_RANGE[1]}")
    alt_flag = note(triggered, "alt", alt, is_above(alt, ALT_RANGE[1]), f">{ALT_RANGE[1]}")
    ggt_flag = note(triggered, "ggt", ggt, is_above(ggt, GGT_RANGE[1]), f">{GGT_RANGE[1]}")
    ratio_flag = note(triggered, "ast/alt", ast_alt_ratio, is_above(ast_alt_ratio, AST_ALT_RATIO_RANGE[1]), f">{AST_ALT_RATIO_RANGE[1]}")
    lft_flag = (
        ast_flag
        or alt_flag
        or ggt_flag
        or ratio_flag
    )

    if (lft_flag and fib4_elevated) or (lft_flag and fib4_elevated and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif fib4_at_risk:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "NAFLD", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "NAFLD"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_nafld_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_nafld',
        'description': 'Evaluate NAFLD susceptibility from fib4, ast, alt, ggt and ast/alt ratio, gated on PNPLA3, TM6SF2, or PRS nafld genetics. Age affects which fib4 threshold tier applies.',
        'parameters': {
            'type': 'object',
            'properties': {
                'fib4': {
                    'type': ['number', 'null'],
                    'description': 'fib4, score. Null if not tested.',
                },
                'ast': {
                    'type': ['number', 'null'],
                    'description': 'ast, U/L. Null if not tested.',
                },
                'alt': {
                    'type': ['number', 'null'],
                    'description': 'alt, U/L. Null if not tested.',
                },
                'ggt': {
                    'type': ['number', 'null'],
                    'description': 'ggt, U/L. Null if not tested.',
                },
                'ast_alt_ratio': {
                    'type': ['number', 'null'],
                    'description': "ast/alt ratio. Null if not tested. NOTE: the labs dict key this function reads is the literal string \"ast/alt\" (not \"ast_alt_ratio\")",
                },
                'platelets': {
                    'type': ['number', 'null'],
                    'description': 'platelets, x10^3/uL. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['NAFLD'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'nafld_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'nafld' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('NAFLD', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('NAFLD', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('NAFLD', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['fib4', 'ast', 'alt', 'ggt', 'ast_alt_ratio', 'platelets', 'age', 'flagged_genes', 'nafld_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Osteoarthritis
# only based on user context, will assign "Elevated susceptibility" only based on PRS/gene and significant if user context present.
# def evaluate_osteoarthritis(labs, patient, genetics, family_history, past_history, symptoms=False):
#     if genetics is None:
#         return [{"Condition": "Osteoarthritis", "Category": INSUFFICIENT_DATA}]
#     if not has_flagged_gene(genetics, "Osteoarthritis"):
#         return [{"Condition": "Osteoarthritis", "Category": GENE_NOT_FOUND}]

#     if symptoms or family_history or past_history:
#         category = "Significant Pattern"
#     else:
#         category = "Elevated susceptibility"
#     # else:
#     #     category = "Typical" #no typical, blood parameters removed and only PRS and user context is used.

#     return [{
#         "Condition": "Osteoarthritis", "Category": category,
#         "DNA Marker(s)": get_gene_trigger(genetics, "Osteoarthritis"),
#         "Blood Marker(s)": "",
#     }]
# evaluate_osteoarthritis_schema = {
#     'type': 'function',
#     'function': {
#         'name': 'evaluate_osteoarthritis',
#         'description': 'Evaluate Osteoarthritis susceptibility - no lab dependency, gated on COL11A1, SMAD3, or PRS osteoarthritis genetics. Once gated, category is driven entirely by symptoms/family_history/past_history.',
#         'parameters': {
#             'type': 'object',
#             'properties': {
#                 'flagged_genes': {
#                     'type': ['array', 'null'],
#                     'items': {
#                         'type': 'string',
#                         'enum': CONDITION_GENES['Osteoarthritis'],
#                     },
#                     'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
#                 },
#                 'osteoarthritis_prs_elevated': {
#                     'type': ['boolean', 'null'],
#                     'description': "True if this patient's PRS 'osteoarthritis' category is elevated or moderately_elevated; null if PRS unavailable.",
#                 },
#                 'family_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Osteoarthritis', {}).get('family_history', []),
#                     },
#                     'description': 'Which of these the patient reports in family history.',
#                 },
#                 'past_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Osteoarthritis', {}).get('past_history', []),
#                     },
#                     'description': "Which of these are in the patient's own past medical history.",
#                 },
#                 'symptoms': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Osteoarthritis', {}).get('symptoms', []),
#                     },
#                     'description': 'Which of these symptoms the patient reports.',
#                 },
#             },
#             'required': ['flagged_genes', 'osteoarthritis_prs_elevated', 'family_history', 'past_history', 'symptoms'],
#             'additionalProperties': False,
#         },
#         'strict': True,
#     },
# }

# Parkinson's Disease
# def evaluate_parkinsons(labs, patient, genetics, family_history, past_history, symptoms=False):
#     if genetics is None:
#         return [{"Condition": "Parkinson's Disease", "Category": INSUFFICIENT_DATA}]
#     if not has_flagged_gene(genetics, "Parkinson's Disease"):
#         return [{"Condition": "Parkinson's Disease", "Category": GENE_NOT_FOUND}]

#     sex = patient.get("sex")
#     hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

#     total_c = labs.get("total_cholesterol")
#     ldl_c = labs.get("ldl_c")
#     triglycerides = labs.get("triglycerides")
#     hdl_c = labs.get("hdl_c")
#     crp = labs.get("crp")

#     triggered = []
#     # --- LDL ---
#     ldl_likely = is_elevated(ldl_c, LDL_C_ELEVATED)
#     ldl_early = is_elevated(ldl_c, LDL_C_AT_RISK)
#     if ldl_likely:
#         note(triggered, "ldl_c", ldl_c, True, f">={LDL_C_ELEVATED}")
#     elif ldl_early:
#         note(triggered, "ldl_c", ldl_c, True, f">={LDL_C_AT_RISK}")

#     # --- Total cholesterol ---
#     total_c_likely = is_elevated(total_c, TOTAL_CHOLESTEROL_ELEVATED)
#     total_c_early = is_elevated(total_c, TOTAL_CHOLESTEROL_AT_RISK)
#     if total_c_likely:
#         note(triggered, "total_cholesterol", total_c, True, f">={TOTAL_CHOLESTEROL_ELEVATED}")
#     elif total_c_early:
#         note(triggered, "total_cholesterol", total_c, True, f">={TOTAL_CHOLESTEROL_AT_RISK}")

#     # --- Triglycerides ---
#     tg_likely = is_elevated(triglycerides, TRIGLYCERIDES_ELEVATED)
#     tg_early = is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK)
#     if tg_likely:
#         note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_ELEVATED}")
#     elif tg_early:
#         note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_AT_RISK}")

#     # --- HDL (single threshold, no early/likely tiers) ---
#     hdl_low = note(triggered, "hdl_c", hdl_c, is_below(hdl_c, hdl_min), f"<{hdl_min}")

#     # dyslipidemia_flag_early = ldl_early or total_c_early or tg_early or hdl_low
#     # dyslipidemia_flag_likely = ldl_likely or total_c_likely or tg_likely or hdl_low
    
#     # manually put down for now, will just report Elevated susceptibility untill updated
#     dyslipidemia_flag_early = False
#     dyslipidemia_flag_likely = False

#     crp_early = is_elevated(crp, CRP_AT_RISK) #and not is_above(crp, CRP_ELEVATED)
#     crp_likely = is_above(crp, CRP_ELEVATED)

#     if crp_likely:
#         note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
#     elif crp_early:
#         note(triggered, "crp", crp, True, f">={CRP_AT_RISK}")

#     if (dyslipidemia_flag_likely and crp_likely) or (dyslipidemia_flag_likely and crp_likely and (family_history or symptoms  or past_history)):
#         category = "Significant Pattern"
#     elif dyslipidemia_flag_early and crp_early:
#         category = "Early Pattern"
#     else:
#         category = "Elevated susceptibility"

#     return [{
#         "Condition": "Parkinson's Disease", "Category": category,
#         "DNA Marker(s)": get_gene_trigger(genetics, "Parkinson's Disease"),
#         "Blood Marker(s)": ", ".join(triggered),
#     }]
# evaluate_parkinsons_schema = {
#     'type': 'function',
#     'function': {
#         'name': 'evaluate_parkinsons',
#         'description': "Evaluate Parkinson's Disease susceptibility, gated on LRRK2, GBA, or PRS parkinson genetics.",
#         'parameters': {
#             'type': 'object',
#             'properties': {
#                 'total_cholesterol': {
#                     'type': ['number', 'null'],
#                     'description': 'total_cholesterol, mg/dL. Null if not tested.',
#                 },
#                 'ldl_c': {
#                     'type': ['number', 'null'],
#                     'description': 'ldl_c, mg/dL. Null if not tested.',
#                 },
#                 'triglycerides': {
#                     'type': ['number', 'null'],
#                     'description': 'triglycerides, mg/dL. Null if not tested.',
#                 },
#                 'hdl_c': {
#                     'type': ['number', 'null'],
#                     'description': 'hdl_c, mg/dL. Null if not tested.',
#                 },
#                 'crp': {
#                     'type': ['number', 'null'],
#                     'description': 'crp, mg/L. Null if not tested.',
#                 },
#                 'sex': {
#                     'type': 'string',
#                     'enum': ['male', 'female'],
#                     'description': "Patient's biological sex.",
#                 },
#                 'flagged_genes': {
#                     'type': ['array', 'null'],
#                     'items': {
#                         'type': 'string',
#                         'enum': CONDITION_GENES["Parkinson's Disease"],
#                     },
#                     'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
#                 },
#                 'parkinson_prs_elevated': {
#                     'type': ['boolean', 'null'],
#                     'description': "True if this patient's PRS 'parkinson' category is elevated or moderately_elevated; null if PRS unavailable.",
#                 },
#                 'family_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get("Parkinson's Disease", {}).get('family_history', []),
#                     },
#                     'description': 'Which of these the patient reports in family history.',
#                 },
#                 'past_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get("Parkinson's Disease", {}).get('past_history', []),
#                     },
#                     'description': "Which of these are in the patient's own past medical history.",
#                 },
#                 'symptoms': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get("Parkinson's Disease", {}).get('symptoms', []),
#                     },
#                     'description': 'Which of these symptoms the patient reports.',
#                 },
#             },
#             'required': ['total_cholesterol', 'ldl_c', 'triglycerides', 'hdl_c', 'crp', 'sex', 'flagged_genes', 'parkinson_prs_elevated', 'family_history', 'past_history', 'symptoms'],
#             'additionalProperties': False,
#         },
#         'strict': True,
#     },
# }

# Psoriasis
def evaluate_psoriasis(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Psoriasis", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Psoriasis"):
        return [{"Condition": "Psoriasis", "Category": GENE_NOT_FOUND}]

    # sex = patient.get("sex")
    age = patient.get("age")
    esr_max = esr_upper_normal(age) # updated to only age dependent, sex based threshold removed (not in metropolis)

    crp = labs.get("crp")
    esr = labs.get("esr")
    nlr = labs.get("nlr")

    triggered = []
    crp_early = is_elevated(crp, CRP_AT_RISK) #and not is_above(crp, CRP_ELEVATED)
    crp_likely = is_above(crp, CRP_ELEVATED)
    esr_flag = note(triggered, "esr", esr, is_above(esr, esr_max), f">{esr_max}")
    nlr_early = is_elevated(nlr, NLR_RANGE[0]) #and not is_above(nlr, NLR_RANGE[1])
    nlr_likely = is_above(nlr, NLR_RANGE[1])

    if crp_likely:
        note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
    elif crp_early:
        note(triggered, "crp", crp, True, f">={CRP_AT_RISK}")

    if nlr_likely:
        note(triggered, "nlr", nlr, True, f">{NLR_RANGE[1]}")
    elif nlr_early:
        note(triggered, "nlr", nlr, True, f">={NLR_RANGE[0]}")

    if (crp_likely and esr_flag and nlr_likely) or (crp_likely and esr_flag and nlr_likely and (family_history or symptoms or past_history)):
        category = "Significant Pattern"
    elif crp_early and esr_flag and nlr_early:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Psoriasis", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Psoriasis"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_psoriasis_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_psoriasis',
        'description': 'Evaluate Psoriasis susceptibility from crp, esr and nlr, gated on PRS psoriasis only (no required gene). Age affects which esr threshold applies.',
        'parameters': {
            'type': 'object',
            'properties': {
                'crp': {
                    'type': ['number', 'null'],
                    'description': 'crp, mg/L. Null if not tested.',
                },
                'esr': {
                    'type': ['number', 'null'],
                    'description': 'esr, mm/hr. Null if not tested.',
                },
                'nlr': {
                    'type': ['number', 'null'],
                    'description': 'nlr, ratio. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'psoriasis_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'psoriasis' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Psoriasis', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Psoriasis', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Psoriasis', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['crp', 'esr', 'nlr', 'age', 'psoriasis_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

#  Rheumatoid Arthritis​
def evaluate_rheumatoid_arthritis(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Rheumatoid Arthritis", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Rheumatoid Arthritis"):
        return [{"Condition": "Rheumatoid Arthritis", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    age = patient.get("age")
    esr_max = esr_upper_normal(age) # updated to only age dependent, sex based threshold removed (not in metropolis)
    hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

    crp = labs.get("crp")
    esr = labs.get("esr")
    total_c = labs.get("total_cholesterol")
    ldl_c = labs.get("ldl_c")
    hdl_c = labs.get("hdl_c")
    triglycerides = labs.get("triglycerides")

    triggered = []
    # --- LDL ---
    ldl_likely = is_elevated(ldl_c, LDL_C_ELEVATED)
    ldl_early = is_elevated(ldl_c, LDL_C_AT_RISK)
    if ldl_likely:
        note(triggered, "ldl_c", ldl_c, True, f">={LDL_C_ELEVATED}")
    elif ldl_early:
        note(triggered, "ldl_c", ldl_c, True, f">={LDL_C_AT_RISK}")

    # --- Total cholesterol ---
    total_c_likely = is_elevated(total_c, TOTAL_CHOLESTEROL_ELEVATED)
    total_c_early = is_elevated(total_c, TOTAL_CHOLESTEROL_AT_RISK)
    if total_c_likely:
        note(triggered, "total_cholesterol", total_c, True, f">={TOTAL_CHOLESTEROL_ELEVATED}")
    elif total_c_early:
        note(triggered, "total_cholesterol", total_c, True, f">={TOTAL_CHOLESTEROL_AT_RISK}")

    # --- Triglycerides ---
    tg_likely = is_elevated(triglycerides, TRIGLYCERIDES_ELEVATED)
    tg_early = is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK)
    if tg_likely:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_ELEVATED}")
    elif tg_early:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_AT_RISK}")

    # --- HDL (single threshold, no early/likely tiers) ---
    hdl_low = note(triggered, "hdl_c", hdl_c, is_below(hdl_c, hdl_min), f"<{hdl_min}")

    dyslipidemia_flag_early = ldl_early or total_c_early or tg_early or hdl_low
    dyslipidemia_flag_likely = ldl_likely or total_c_likely or tg_likely or hdl_low

    esr_flag = note(triggered, "esr", esr, is_above(esr, esr_max), f">{esr_max}")
    crp_early = is_elevated(crp, CRP_AT_RISK) #and not is_above(crp, CRP_ELEVATED)
    crp_likely = is_above(crp, CRP_ELEVATED)
    if crp_likely:
        note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
    elif crp_early:
        note(triggered, "crp", crp, True, f">={CRP_AT_RISK}")

    if (dyslipidemia_flag_likely and crp_likely and esr_flag) or (dyslipidemia_flag_likely and crp_likely and esr_flag and (family_history or symptoms or past_history)):
        category = "Significant Pattern"
    elif dyslipidemia_flag_early and crp_early and esr_flag:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Rheumatoid Arthritis", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Rheumatoid Arthritis"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_rheumatoid_arthritis_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_rheumatoid_arthritis',
        'description': 'Evaluate Rheumatoid Arthritis susceptibility from crp, esr and lipid panel, gated on PTPN22, STAT4, TNFAIP3, CTLA4, or PRS rheumatoidarthritis genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'crp': {
                    'type': ['number', 'null'],
                    'description': 'crp, mg/L. Null if not tested.',
                },
                'esr': {
                    'type': ['number', 'null'],
                    'description': 'esr, mm/hr. Null if not tested.',
                },
                'total_cholesterol': {
                    'type': ['number', 'null'],
                    'description': 'total_cholesterol, mg/dL. Null if not tested.',
                },
                'ldl_c': {
                    'type': ['number', 'null'],
                    'description': 'ldl_c, mg/dL. Null if not tested.',
                },
                'hdl_c': {
                    'type': ['number', 'null'],
                    'description': 'hdl_c, mg/dL. Null if not tested.',
                },
                'triglycerides': {
                    'type': ['number', 'null'],
                    'description': 'triglycerides, mg/dL. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Rheumatoid Arthritis'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'rheumatoidarthritis_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'rheumatoidarthritis' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Rheumatoid Arthritis', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Rheumatoid Arthritis', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Rheumatoid Arthritis', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['crp', 'esr', 'total_cholesterol', 'ldl_c', 'hdl_c', 'triglycerides', 'age', 'sex', 'flagged_genes', 'rheumatoidarthritis_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Rhinitis 
# def evaluate_rhinitis(labs, patient, genetics, family_history, past_history, symptoms=False):
#     if genetics is None:
#         return [{"Condition": "Rhinitis", "Category": INSUFFICIENT_DATA}]
#     if not has_flagged_gene(genetics, "Rhinitis"):
#         return [{"Condition": "Rhinitis", "Category": GENE_NOT_FOUND}]

#     eosinophils = labs.get("eosinophils")
#     ige = labs.get("ige")

#     triggered = []
#     eosinophil_flag_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
#     eosinophil_flag_likely = is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
#     if eosinophil_flag_likely:
#         note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_UPPER_NORMAL}")
#     elif eosinophil_flag_early:
#         note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")
#     ige_flag = note(triggered, "ige", ige, is_above(ige, IGE_UPPER_NORMAL), f">{IGE_UPPER_NORMAL}")

#     if (eosinophil_flag_likely and ige_flag) or (eosinophil_flag_likely and ige_flag and (symptoms or family_history or past_history)):
#         category = "Significant Pattern"
#     elif eosinophil_flag_early or ige_flag:
#         category = "Early Pattern"
#     else:
#         category = "Elevated susceptibility"

#     return [{
#         "Condition": "Rhinitis", "Category": category,
#         "DNA Marker(s)": get_gene_trigger(genetics, "Rhinitis"),
#         "Blood Marker(s)": ", ".join(triggered),
#     }]
# evaluate_rhinitis_schema = {
#     'type': 'function',
#     'function': {
#         'name': 'evaluate_rhinitis',
#         'description': 'Evaluate Rhinitis susceptibility from eosinophil count and ige, gated on TSLP or PRS rhinitis genetics.',
#         'parameters': {
#             'type': 'object',
#             'properties': {
#                 'eosinophils': {
#                     'type': ['number', 'null'],
#                     'description': 'eosinophils, x10^9/L. Null if not tested.',
#                 },
#                 'ige': {
#                     'type': ['number', 'null'],
#                     'description': 'ige, IU/mL. Null if not tested.',
#                 },
#                 'flagged_genes': {
#                     'type': ['array', 'null'],
#                     'items': {
#                         'type': 'string',
#                         'enum': CONDITION_GENES['Rhinitis'],
#                     },
#                     'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
#                 },
#                 'rhinitis_prs_elevated': {
#                     'type': ['boolean', 'null'],
#                     'description': "True if this patient's PRS 'rhinitis' category is elevated or moderately_elevated; null if PRS unavailable.",
#                 },
#                 'family_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Rhinitis', {}).get('family_history', []),
#                     },
#                     'description': 'Which of these the patient reports in family history.',
#                 },
#                 'past_history': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Rhinitis', {}).get('past_history', []),
#                     },
#                     'description': "Which of these are in the patient's own past medical history.",
#                 },
#                 'symptoms': {
#                     'type': 'array',
#                     'items': {
#                         'type': 'string',
#                         'enum': USER_CONTEXT.get('Rhinitis', {}).get('symptoms', []),
#                     },
#                     'description': 'Which of these symptoms the patient reports.',
#                 },
#             },
#             'required': ['eosinophils', 'ige', 'flagged_genes', 'rhinitis_prs_elevated', 'family_history', 'past_history', 'symptoms'],
#             'additionalProperties': False,
#         },
#         'strict': True,
#     },
# }

# Type 2 Diabetes
def evaluate_t2d(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Type 2 Diabetes", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Type 2 Diabetes"):
        return [{"Condition": "Type 2 Diabetes", "Category": GENE_NOT_FOUND}]

    fasting_glucose = labs.get("fasting_glucose")
    fasting_insulin = labs.get("fasting_insulin")
    hba1c = labs.get("hba1c")

    triggered = []
    hba1c_elevated = is_elevated(hba1c, HBA1C_ELEVATED)
    hba1c_at_risk = is_elevated(hba1c, HBA1C_AT_RISK)
    if hba1c_elevated:
        note(triggered, "hba1c", hba1c, True, f">={HBA1C_ELEVATED}")
    elif hba1c_at_risk:
        note(triggered, "hba1c", hba1c, True, f">={HBA1C_AT_RISK}")

    fg_elevated = is_elevated(fasting_glucose, FASTING_GLUCOSE_ELEVATED)
    fg_at_risk = is_elevated(fasting_glucose, FASTING_GLUCOSE_AT_RISK)
    if fg_elevated:
        note(triggered, "fasting_glucose", fasting_glucose, True, f">={FASTING_GLUCOSE_ELEVATED}")
    elif fg_at_risk:
        note(triggered, "fasting_glucose", fasting_glucose, True, f">={FASTING_GLUCOSE_AT_RISK}")

    hyperglycemia = hba1c_elevated or fg_elevated
    hyperglycemia_early = hba1c_at_risk or fg_at_risk
    insulin_resistance_elevated = is_above(fasting_insulin, FASTING_INSULIN_RANGE[1])
    insulin_resistance_early = is_above(fasting_insulin, FASTING_INSULIN_RANGE[0])

    if insulin_resistance_elevated:
        note(triggered, "fasting_insulin", fasting_insulin, True, f">{FASTING_INSULIN_RANGE[1]}")
    elif insulin_resistance_early:
        note(triggered, "fasting_insulin", fasting_insulin, True, f">{FASTING_INSULIN_RANGE[0]}")

    if (hyperglycemia and insulin_resistance_elevated) or (hyperglycemia and insulin_resistance_elevated and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif hyperglycemia_early and insulin_resistance_early:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Type 2 Diabetes", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Type 2 Diabetes"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_t2d_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_t2d',
        'description': 'Evaluate Type 2 Diabetes susceptibility from fasting glucose, fasting insulin and hba1c, gated on TCF7L2, SLC30A8, or PRS type2diabetes genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'fasting_glucose': {
                    'type': ['number', 'null'],
                    'description': 'fasting_glucose, mg/dL. Null if not tested.',
                },
                'fasting_insulin': {
                    'type': ['number', 'null'],
                    'description': 'fasting_insulin, uIU/mL. Null if not tested.',
                },
                'hba1c': {
                    'type': ['number', 'null'],
                    'description': 'hba1c, %. Null if not tested.',
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Type 2 Diabetes'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'type2diabetes_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'type2diabetes' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Type 2 Diabetes', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Type 2 Diabetes', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Type 2 Diabetes', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['fasting_glucose', 'fasting_insulin', 'hba1c', 'flagged_genes', 'type2diabetes_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# =============================================================================
# (14)
# =============================================================================

# Maturity-Onset Diabetes of the Young"
def evaluate_mody(labs, patient, genetics, family_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Maturity-Onset Diabetes of the Young", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Maturity-Onset Diabetes of the Young"):
        return [{"Condition": "Maturity-Onset Diabetes of the Young", "Category": GENE_NOT_FOUND}]

    hba1c = labs.get("hba1c")
    fasting_glucose = labs.get("fasting_glucose")
    age = patient.get("age")

    # ref = user_context["MODY"]
    # fhx_terms = ref["family_history"]
    # family_history = family_history or []
    # family_history = any(term in family_history for term in fhx_terms)

    triggered = []
    hba1c_elevated = is_elevated(hba1c, HBA1C_ELEVATED)
    hba1c_at_risk = is_elevated(hba1c, HBA1C_AT_RISK)
    if hba1c_elevated:
        note(triggered, "hba1c", hba1c, True, f">={HBA1C_ELEVATED}")
    elif hba1c_at_risk:
        note(triggered, "hba1c", hba1c, True, f">={HBA1C_AT_RISK}")

    fg_elevated = is_elevated(fasting_glucose, FASTING_GLUCOSE_ELEVATED)
    fg_at_risk = is_elevated(fasting_glucose, FASTING_GLUCOSE_AT_RISK)
    if fg_elevated:
        note(triggered, "fasting_glucose", fasting_glucose, True, f">={FASTING_GLUCOSE_ELEVATED}")
    elif fg_at_risk:
        note(triggered, "fasting_glucose", fasting_glucose, True, f">={FASTING_GLUCOSE_AT_RISK}")

    hyperglycemia = hba1c_elevated or fg_elevated
    hyperglycemia_early = hba1c_at_risk or fg_at_risk

    if (hyperglycemia and age is not None and age < 25) or (hyperglycemia and (family_history or symptoms)):
        category = "Significant Pattern"
    elif hyperglycemia_early or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Maturity-Onset Diabetes of the Young", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Maturity-Onset Diabetes of the Young"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_mody_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_mody',
        'description': 'Evaluate Maturity-Onset Diabetes of the Young (MODY) susceptibility from hba1c/fasting glucose, patient age, and HNF1A genetics. No past_history dependency.',
        'parameters': {
            'type': 'object',
            'properties': {
                'hba1c': {
                    'type': ['number', 'null'],
                    'description': 'hba1c, %. Null if not tested.',
                },
                'fasting_glucose': {
                    'type': ['number', 'null'],
                    'description': 'fasting_glucose, mg/dL. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Maturity-Onset Diabetes of the Young'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Maturity-Onset Diabetes of the Young', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Maturity-Onset Diabetes of the Young', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['hba1c', 'fasting_glucose', 'age', 'flagged_genes', 'family_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# pth not covered in blood panel

## Multiple Endocrine Neoplasia Type 4
# def evaluate_men4(labs, patient, genetics, family_history, endocrine_conditions=None):
#     if not has_flagged_gene(genetics, "Multiple Endocrine Neoplasia Type 4"):
#         return [{"Condition": "Multiple Endocrine Neoplasia Type 4", "Category": GENE_NOT_FOUND}]

#     calcium = labs.get("calcium")
#     pth = labs.get("pth")

#     triggered = []
#     hypercalcemia = note(triggered, "calcium", calcium, is_above(calcium, CALCIUM_UPPER_NORMAL), f">{CALCIUM_UPPER_NORMAL}")
#     elevated_pth = note(triggered, "pth", pth, is_above(pth, PTH_UPPER_NORMAL), f">{PTH_UPPER_NORMAL}")
#     phpt_pattern = hypercalcemia and elevated_pth

#     if phpt_pattern and (endocrine_conditions or family_history):
#         category = "Significant Pattern"
#     elif phpt_pattern or (endocrine_conditions and family_history):
#         category = "Early Pattern"
#     else:
#         category = "Elevated susceptibility"

#     return [{
#         "Condition": "Multiple Endocrine Neoplasia Type 4", "Category": category,
#         "DNA Marker(s)": get_gene_trigger(genetics, "Multiple Endocrine Neoplasia Type 4"),
#         "Blood Marker(s)": ", ".join(triggered),
#     }]


# Muscular Dystrophy
def evaluate_muscular_dystrophy(labs, patient, genetics, family_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Muscular Dystrophy", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Muscular Dystrophy"):
        return [{"Condition": "Muscular Dystrophy", "Category": GENE_NOT_FOUND}]

    ck = labs.get("ck")

    ck_significant_multiple = CK_UPPER_NORMAL * 10
    ck_moderate_multiple = CK_UPPER_NORMAL * 5

    triggered = []
    marked_ck_elevation = is_elevated(ck, ck_significant_multiple)
    moderate_ck_elevation = is_elevated(ck, ck_moderate_multiple)
    if marked_ck_elevation:
        note(triggered, "ck", ck, True, f">={ck_significant_multiple}")
    elif moderate_ck_elevation:
        note(triggered, "ck", ck, True, f">={ck_moderate_multiple}")

    if marked_ck_elevation and (symptoms or family_history):
        category = "Significant Pattern"
    elif marked_ck_elevation or (moderate_ck_elevation and (symptoms or family_history)):
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Muscular Dystrophy", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Muscular Dystrophy"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_muscular_dystrophy_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_muscular_dystrophy',
        'description': 'Evaluate Muscular Dystrophy susceptibility from creatine kinase (ck), gated on DMD, FKRP, LMNA, SGCA, SGCB genetics. No past_history dependency.',
        'parameters': {
            'type': 'object',
            'properties': {
                'ck': {
                    'type': ['number', 'null'],
                    'description': 'ck, U/L. Null if not tested.',
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Muscular Dystrophy'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Muscular Dystrophy', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Muscular Dystrophy', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['ck', 'flagged_genes', 'family_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Familial Hypercholesterolemia
def evaluate_fh(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Familial Hypercholesterolemia", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Familial Hypercholesterolemia"):
        return [{"Condition": "Familial Hypercholesterolemia", "Category": GENE_NOT_FOUND}]

    ldl = labs.get("ldl_c")
    non_hdl = labs.get("non_hdl_c")
    lpa = labs.get("lpa")
    age = patient.get("age")

    triggered = []
    elevated_early = (
        is_above(ldl, LDL_C_AT_RISK) or 
        is_above(non_hdl, NON_HDL_C_ELEVATED)
        # is_elevated(lpa, LPA_AT_RISK) #is only getting noted but not responsible for any trigger (will only be printed if combined with another trigger) 
    )
    elevated_likely = (
        is_above(ldl, LDL_C_ELEVATED) or 
        is_above(non_hdl, NON_HDL_C_SEVERE)
        # is_elevated(lpa, LPA_ELEVATED) #is only getting noted but not responsible for any trigger (will only be printed if combined with another trigger) 
    )

    # changed for noting mild crossed values as well if significant is already triggered by another marker (for flahyAI and GP reporting)
    ldl_likely_flag = is_above(ldl, LDL_C_ELEVATED)
    ldl_early_flag = is_above(ldl, LDL_C_AT_RISK)
    if ldl_likely_flag:
        note(triggered, "ldl_c", ldl, True, f">{LDL_C_ELEVATED}")
    elif ldl_early_flag:
        note(triggered, "ldl_c", ldl, True, f">{LDL_C_AT_RISK}")

    non_hdl_likely_flag = is_above(non_hdl, NON_HDL_C_SEVERE)
    non_hdl_early_flag = is_above(non_hdl, NON_HDL_C_ELEVATED)
    if non_hdl_likely_flag:
        note(triggered, "non_hdl_c", non_hdl, True, f">{NON_HDL_C_SEVERE}")
    elif non_hdl_early_flag:
        note(triggered, "non_hdl_c", non_hdl, True, f">{NON_HDL_C_ELEVATED}")

    lpa_likely_flag = is_elevated(lpa, LPA_ELEVATED)
    lpa_early_flag = is_elevated(lpa, LPA_AT_RISK)
    if lpa_likely_flag:
        note(triggered, "lpa", lpa, True, f">={LPA_ELEVATED}")
    elif lpa_early_flag:
        note(triggered, "lpa", lpa, True, f">={LPA_AT_RISK}")

    if elevated_likely or (elevated_early and age is not None and age < 40) or (elevated_early and (family_history or symptoms or past_history)):
        category = "Significant Pattern"
    elif elevated_early or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Familial Hypercholesterolemia", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Familial Hypercholesterolemia"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_fh_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_fh',
        'description': 'Evaluate Familial Hypercholesterolemia susceptibility from ldl_c, non_hdl_c and lpa, gated on LDLR, APOB, PCSK9 genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'ldl_c': {
                    'type': ['number', 'null'],
                    'description': 'ldl_c, mg/dL. Null if not tested.',
                },
                'non_hdl_c': {
                    'type': ['number', 'null'],
                    'description': 'non_hdl_c, mg/dL. Null if not tested.',
                },
                'lpa': {
                    'type': ['number', 'null'],
                    'description': 'lpa, mg/dL. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Familial Hypercholesterolemia'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Familial Hypercholesterolemia', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Familial Hypercholesterolemia', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Familial Hypercholesterolemia', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['ldl_c', 'non_hdl_c', 'lpa', 'age', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Cardiomyopathy
def evaluate_cardiomyopathy(labs, patient, genetics, family_history, past_history, symptoms=False, imaging_performed=False, imaging_findings=None,):
    if genetics is None:
        return [{"Condition": "Cardiomyopathy", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Cardiomyopathy"):
        return [{"Condition": "Cardiomyopathy", "Category": GENE_NOT_FOUND}]

    hemoglobin = labs.get("hemoglobin")
    eosinophils = labs.get("eosinophils")
    sex = patient.get("sex")

    HIGH_SPECIFICITY_FINDINGS = USER_CONTEXT.get("Cardiomyopathy", {}).get("high_specificity_findings", [])  # subset of "imaging_findings" in blood_conditions_user_context.json

    imaging_findings = imaging_findings or []
    real_findings = [f for f in imaging_findings if f != "Normal"]

    imaging_abnormal = imaging_performed and len(real_findings) > 0
    imaging_multiple = imaging_performed and len(real_findings) >= 2
    imaging_high_spec = imaging_performed and any(
        f in HIGH_SPECIFICITY_FINDINGS for f in real_findings
    )

    triggered = []
    hb_normal = sex_based_threshold(sex, 12.5, 11.5)
    hb_severely_low = hb_normal is not None and note(triggered, "hemoglobin", hemoglobin, is_below(hemoglobin, hb_normal - 2), f"<{hb_normal - 2}") #because python evaluates arithmetic part before helper function's "not None" check.
    eos_high = note(triggered, "eosinophils", eosinophils, is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL), f">{EOSINOPHILS_UPPER_NORMAL}")
    labs_severe = hb_severely_low or eos_high

    gene_trigger = get_gene_trigger(genetics, "Cardiomyopathy")

    if hemoglobin is None and eosinophils is None and not imaging_performed:
        return [{
            "Condition": "Cardiomyopathy",
            "Category": "Elevated susceptibility",
            "DNA Marker(s)": gene_trigger,
            "Blood Marker(s)": ", ".join(triggered),
        }]

    # --- Case 1: imaging performed and abnormal -> imaging drives the call ---
    if imaging_performed and imaging_abnormal:
        if imaging_high_spec or imaging_multiple or symptoms or family_history or past_history:
            category = "Significant Pattern"
        else:
            category = "Early Pattern"

    # --- Case 2: imaging performed and normal ---
    elif imaging_performed and not imaging_abnormal:
        risk_factor_count = sum([family_history, past_history, symptoms, labs_severe])
        if risk_factor_count >= 2:
            category = "Early Pattern" #Re-imaging recommended (Echocardiogram/Cardiac MRI)
        else:
            category = "Elevated susceptibility"

    # --- Case 3: no imaging performed -> route to imaging ---
    else:
        if labs_severe or family_history or symptoms or past_history:
            category = "Early Pattern" #Advanced cardiac imaging recommended (Echocardiogram/Cardiac MRI)
        else:
            category = "Elevated susceptibility"

    return [{
        "Condition": "Cardiomyopathy", "Category": category,
        "DNA Marker(s)": gene_trigger,
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_cardiomyopathy_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_cardiomyopathy',
        'description': "Evaluate Cardiomyopathy susceptibility, gated on any of a 22-gene cardiomyopathy panel. Imaging result (if performed) drives the category most directly; hemoglobin/eosinophils and symptoms/family_history/past_history are the fallback risk factors when imaging hasn't been done or was normal.",
        'parameters': {
            'type': 'object',
            'properties': {
                'hemoglobin': {
                    'type': ['number', 'null'],
                    'description': 'hemoglobin, g/dL. Null if not tested.',
                },
                'eosinophils': {
                    'type': ['number', 'null'],
                    'description': 'eosinophils, x10^9/L. Null if not tested.',
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'imaging_performed': {
                    'type': 'boolean',
                    'description': 'True if cardiac imaging (echocardiogram/cardiac MRI) was performed.',
                },
                'imaging_findings': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Cardiomyopathy', {}).get('imaging_findings', []),
                    },
                    'description': "Imaging findings if performed - one or more of these, e.g. ['Normal'] or specific abnormalities. Null if imaging not performed.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Cardiomyopathy'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Cardiomyopathy', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Cardiomyopathy', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Cardiomyopathy', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['hemoglobin', 'eosinophils', 'sex', 'imaging_performed', 'imaging_findings', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Coronary Artery Disease
def evaluate_cad(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Coronary Artery Disease", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Coronary Artery Disease"):
        return [{"Condition": "Coronary Artery Disease", "Category": GENE_NOT_FOUND}]

    ldl = labs.get("ldl_c")
    non_hdl = labs.get("non_hdl_c")
    hdl = labs.get("hdl_c")
    triglycerides = labs.get("triglycerides")
    hba1c = labs.get("hba1c")
    apob = labs.get("apob")
    lpa = labs.get("lpa")
    age = patient.get("age")
    sex = patient.get("sex")

    hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

    triggered = []
    ldl_moderate = (
        is_elevated(ldl, LDL_C_AT_RISK) or 
        is_elevated(non_hdl, NON_HDL_C_ELEVATED) or
        is_elevated(apob, APOB_AT_RISK)
        # is_elevated(lpa, LPA_AT_RISK) #is only getting noted but not responsible for any trigger (will only be printed if combined with another trigger) 
    )
    ldl_high = (
        is_elevated(ldl, LDL_C_ELEVATED) or 
        is_elevated(non_hdl, NON_HDL_C_SEVERE) or
        is_elevated(apob, APOB_ELEVATED)
        # is_elevated(lpa, LPA_ELEVATED) #is only getting noted but not responsible for any trigger (will only be printed if combined with another trigger)
    )

    # changed for noting mild crossed values as well if significant is already triggered by another marker (for flahyAI and GP reporting)
    ldl_likely_flag = is_elevated(ldl, LDL_C_ELEVATED)
    ldl_early_flag = is_elevated(ldl, LDL_C_AT_RISK)
    if ldl_likely_flag:
        note(triggered, "ldl_c", ldl, True, f">={LDL_C_ELEVATED}")
    elif ldl_early_flag:
        note(triggered, "ldl_c", ldl, True, f">={LDL_C_AT_RISK}")

    non_hdl_likely_flag = is_elevated(non_hdl, NON_HDL_C_SEVERE)
    non_hdl_early_flag = is_elevated(non_hdl, NON_HDL_C_ELEVATED)
    if non_hdl_likely_flag:
        note(triggered, "non_hdl_c", non_hdl, True, f">={NON_HDL_C_SEVERE}")
    elif non_hdl_early_flag:
        note(triggered, "non_hdl_c", non_hdl, True, f">={NON_HDL_C_ELEVATED}")  

    apob_likely_flag = is_elevated(apob, APOB_ELEVATED)
    apob_early_flag = is_elevated(apob, APOB_AT_RISK)
    if apob_likely_flag:
        note(triggered, "apob", apob, True, f">={APOB_ELEVATED}")
    elif apob_early_flag:
        note(triggered, "apob", apob, True, f">={APOB_AT_RISK}")

    lpa_likely_flag = is_elevated(lpa, LPA_ELEVATED)
    lpa_early_flag = is_elevated(lpa, LPA_AT_RISK)
    if lpa_likely_flag:
        note(triggered, "lpa", lpa, True, f">={LPA_ELEVATED}")
    elif lpa_early_flag:
        note(triggered, "lpa", lpa, True, f">={LPA_AT_RISK}")

    low_hdl = note(triggered, "hdl_c", hdl, is_below(hdl, hdl_min), f"<{hdl_min}")
    tg_moderate = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")

    prediabetes = is_elevated(hba1c, HBA1C_AT_RISK)
    diabetes = is_elevated(hba1c, HBA1C_ELEVATED)
    if diabetes:
        note(triggered, "hba1c", hba1c, True, f">={HBA1C_ELEVATED}")
    elif prediabetes:
        note(triggered, "hba1c", hba1c, True, f">={HBA1C_AT_RISK}")

    any_moderate = ldl_moderate or low_hdl or tg_moderate or prediabetes
    combined_atherogenic = ldl_moderate and low_hdl and tg_moderate

    if (
        ldl_high or
        diabetes or
        combined_atherogenic or
        (any_moderate and symptoms and age is not None and age < 40) or
        (any_moderate and symptoms and (family_history or past_history))
    ):
        category = "Significant Pattern"
    # (any_moderate + age) added for early pattern w/o user context (labs only)
    elif family_history or (any_moderate and (symptoms or past_history)) or (any_moderate and age is not None and age < 40):
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Coronary Artery Disease", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Coronary Artery Disease"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_cad_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_cad',
        'description': 'Evaluate Coronary Artery Disease susceptibility from lipid panel, hba1c and apob, gated on a 15-gene CAD panel or PRS cad genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'ldl_c': {
                    'type': ['number', 'null'],
                    'description': 'ldl_c, mg/dL. Null if not tested.',
                },
                'non_hdl_c': {
                    'type': ['number', 'null'],
                    'description': 'non_hdl_c, mg/dL. Null if not tested.',
                },
                'hdl_c': {
                    'type': ['number', 'null'],
                    'description': 'hdl_c, mg/dL. Null if not tested.',
                },
                'triglycerides': {
                    'type': ['number', 'null'],
                    'description': 'triglycerides, mg/dL. Null if not tested.',
                },
                'hba1c': {
                    'type': ['number', 'null'],
                    'description': 'hba1c, %. Null if not tested.',
                },
                'apob': {
                    'type': ['number', 'null'],
                    'description': 'apob, mg/dL. Null if not tested.',
                },
                'lpa': {
                    'type': ['number', 'null'],
                    'description': 'lpa, mg/dL. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Coronary Artery Disease'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'cad_prs_elevated': {
                    'type': ['boolean', 'null'],
                    'description': "True if this patient's PRS 'cad' category is elevated or moderately_elevated; null if PRS unavailable.",
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Coronary Artery Disease', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Coronary Artery Disease', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Coronary Artery Disease', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['ldl_c', 'non_hdl_c', 'hdl_c', 'triglycerides', 'hba1c', 'apob', 'lpa', 'age', 'sex', 'flagged_genes', 'cad_prs_elevated', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Familial Hypertriglyceridemia
def evaluate_hypertriglyceridemia(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Familial Hypertriglyceridemia", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Familial Hypertriglyceridemia"):
        return [{"Condition": "Familial Hypertriglyceridemia", "Category": GENE_NOT_FOUND}]

    triglycerides = labs.get("triglycerides")
    hdl = labs.get("hdl_c")
    vldl = labs.get("vldl")
    sex = patient.get("sex")
    age = patient.get("age")

    triggered = []
    hdl_normal = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)
    hdl_severe_threshold = sex_based_threshold(sex, HDL_C_MALE_SEVERE, HDL_C_FEMALE_SEVERE)
    hdl_low = is_below(hdl, hdl_normal)
    hdl_severely_low = is_below(hdl, hdl_severe_threshold)
    if hdl_severely_low:
        note(triggered, "hdl_c", hdl, True, f"<{hdl_severe_threshold}")
    elif hdl_low:
        note(triggered, "hdl_c", hdl, True, f"<{hdl_normal}")

    tg_borderline = is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK)
    tg_moderate = is_elevated(triglycerides, TRIGLYCERIDES_ELEVATED)
    tg_severe = is_elevated(triglycerides, TRIGLYCERIDES_SEVERE)
    tg_very_severe = is_elevated(triglycerides, TRIGLYCERIDES_VERY_SEVERE)
    if tg_very_severe:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_VERY_SEVERE}")
    elif tg_severe:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_SEVERE}")
    elif tg_moderate:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_ELEVATED}")
    elif tg_borderline:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_AT_RISK}")

    vldl_elevated = note(triggered, "vldl", vldl, is_above(vldl, VLDL_UPPER_NORMAL), f">{VLDL_UPPER_NORMAL}")
    any_moderate_abnormal = tg_moderate or hdl_low or vldl_elevated

    if (
        tg_very_severe or
        (tg_severe and age is not None and age < 40) or
        (tg_severe and (family_history or past_history)) or
        (tg_severe and hdl_severely_low)
    ):
        category = "Significant Pattern"
    elif (
        any_moderate_abnormal or
        (tg_borderline and hdl_low) or
        family_history
    ):
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Familial Hypertriglyceridemia", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Familial Hypertriglyceridemia"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_hypertriglyceridemia_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_hypertriglyceridemia',
        'description': 'Evaluate Familial Hypertriglyceridemia susceptibility from triglycerides, hdl_c and vldl, gated on CREB3L3, LMF1, LPL genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'triglycerides': {
                    'type': ['number', 'null'],
                    'description': 'triglycerides, mg/dL. Null if not tested.',
                },
                'hdl_c': {
                    'type': ['number', 'null'],
                    'description': 'hdl_c, mg/dL. Null if not tested.',
                },
                'vldl': {
                    'type': ['number', 'null'],
                    'description': 'vldl, mg/dL. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Familial Hypertriglyceridemia'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Familial Hypertriglyceridemia', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Familial Hypertriglyceridemia', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Familial Hypertriglyceridemia', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['triglycerides', 'hdl_c', 'vldl', 'age', 'sex', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# HDL Deficiency
def evaluate_hdl_deficiency(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "HDL Deficiency", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "HDL Deficiency"):
        return [{"Condition": "HDL Deficiency", "Category": GENE_NOT_FOUND}]

    triglycerides = labs.get("triglycerides")
    hdl = labs.get("hdl_c")
    ldl = labs.get("ldl_c")
    sex = patient.get("sex")
    age = patient.get("age")

    triggered = []
    hdl_normal = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)
    hdl_severe_threshold = sex_based_threshold(sex, 25, 30)
    hdl_low = is_below(hdl, hdl_normal)
    hdl_severe = is_below(hdl, hdl_severe_threshold)
    if hdl_severe:
        note(triggered, "hdl_c", hdl, True, f"<{hdl_severe_threshold}")
    elif hdl_low:
        note(triggered, "hdl_c", hdl, True, f"<{hdl_normal}")

    tg_elevated = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")
    ldl_elevated = note(triggered, "ldl_c", ldl, is_elevated(ldl, LDL_C_AT_RISK), f">={LDL_C_AT_RISK}")

    isolated_significant = hdl_severe and not tg_elevated and not ldl_elevated
    isolated_moderate = hdl_low and not tg_elevated and not ldl_elevated
    mixed_moderate = hdl_low and (tg_elevated or ldl_elevated)

    if (
        isolated_significant or
        (isolated_moderate and ((age is not None and age < 40) or family_history or past_history))
    ):
        category = "Significant Pattern"
    elif isolated_moderate or mixed_moderate or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "HDL Deficiency", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "HDL Deficiency"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_hdl_deficiency_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_hdl_deficiency',
        'description': 'Evaluate HDL Deficiency susceptibility from hdl_c, triglycerides and ldl_c, gated on ABCA1, LCAT, APOA1 genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'triglycerides': {
                    'type': ['number', 'null'],
                    'description': 'triglycerides, mg/dL. Null if not tested.',
                },
                'hdl_c': {
                    'type': ['number', 'null'],
                    'description': 'hdl_c, mg/dL. Null if not tested.',
                },
                'ldl_c': {
                    'type': ['number', 'null'],
                    'description': 'ldl_c, mg/dL. Null if not tested.',
                },
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['HDL Deficiency'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('HDL Deficiency', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('HDL Deficiency', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('HDL Deficiency', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['triglycerides', 'hdl_c', 'ldl_c', 'age', 'sex', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# all terminology are kept old, not updated for this, to have reference.
# Elevated Apolipoprotein B
# def evaluate_apob(labs, patient, genetics, family_history, symptoms=False):
#     if not has_flagged_gene(genetics, "Elevated Apolipoprotein B"):
#         return [{"Condition": "Elevated Apolipoprotein B", "Category": GENE_NOT_FOUND}]

#     apob = labs.get("apob")
#     age = patient.get("age")

#     triggered = []
#     elevated_100 = is_elevated(apob, 100)
#     elevated_130 = is_elevated(apob, 130)
#     if elevated_130:
#         note(triggered, "apob", apob, True, f">={130}")
#     elif elevated_100:
#         note(triggered, "apob", apob, True, f">={100}")

#     if elevated_130 or (elevated_100 and age is not None and age < 40) or (elevated_100 and family_history):
#         category = "Likely Disease Onset"
#     elif elevated_100 or family_history:
#         category = "Early Pattern"
#     else:
#         category = "Typical"

#     return [{
#         "Condition": "Elevated Apolipoprotein B", "Category": category,
#         "Triggering Genetics": get_gene_trigger(genetics, "Elevated Apolipoprotein B"),
#         "Triggering Parameters": ", ".join(triggered),
#     }]


# Li-Fraumeni Syndrome
def evaluate_li_fraumeni_syndrome(labs, patient, genetics, family_history, symptoms=False, conditions=None):
    if genetics is None:
        return [{"Condition": "Li-Fraumeni Syndrome", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Li-Fraumeni Syndrome"):
        return [{"Condition": "Li-Fraumeni Syndrome", "Category": GENE_NOT_FOUND}]

    age = patient.get("age")
    conditions = conditions or []

    pathognomonic_tumors = USER_CONTEXT.get("Li-Fraumeni Syndrome", {}).get("pathognomonic_tumors", [])  #subset of conditions
    lfs_spectrum_tumors = USER_CONTEXT.get("Li-Fraumeni Syndrome", {}).get("lfs_spectrum_tumors", [])  #subset of conditions

    has_pathognomonic_tumor = any(t in conditions for t in pathognomonic_tumors)
    proband_lfs_tumors = [t for t in lfs_spectrum_tumors if t in conditions]
    has_lfs_tumor = len(proband_lfs_tumors) > 0
    multiple_primaries = len(proband_lfs_tumors) >= 2

    early_onset = age is not None and age < 46

    very_early_breast_cancer = (
        "Breast cancer" in conditions and age is not None and age < 31
    )

    if (
        has_pathognomonic_tumor or
        very_early_breast_cancer or
        family_history or
        multiple_primaries or
        (has_lfs_tumor and early_onset and family_history)
    ):
        category = "Significant Pattern"
    elif has_lfs_tumor or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Li-Fraumeni Syndrome", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Li-Fraumeni Syndrome"),
        "Blood Marker(s)": "",
    }]
evaluate_li_fraumeni_syndrome_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_li_fraumeni_syndrome',
        'description': 'Evaluate Li-Fraumeni Syndrome susceptibility - no lab dependency, gated on TP53 genetics. No past_history dependency.',
        'parameters': {
            'type': 'object',
            'properties': {
                'age': {
                    'type': 'integer',
                    'description': "Patient's age in years.",
                },
                'conditions': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Li-Fraumeni Syndrome', {}).get('lfs_spectrum_tumors', []),
                    },
                    'description': "Prior cancer diagnoses this patient has had, from this condition's LFS-spectrum tumor list (a subset of these are also pathognomonic on their own).",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Li-Fraumeni Syndrome'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Li-Fraumeni Syndrome', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Li-Fraumeni Syndrome', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['age', 'conditions', 'flagged_genes', 'family_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# PTEN Hamartoma Tumor Syndrome
def evaluate_phts(labs, patient, genetics, family_history, symptoms=False, conditions=None):
    if genetics is None:
        return [{"Condition": "PTEN Hamartoma Tumor Syndrome", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "PTEN Hamartoma Tumor Syndrome"):
        return [{"Condition": "PTEN Hamartoma Tumor Syndrome", "Category": GENE_NOT_FOUND}]

    conditions = conditions or []
    major_criteria = USER_CONTEXT.get("PTEN Hamartoma Tumor Syndrome", {}).get("major_criteria", [])   # subset of conditions
    minor_criteria = USER_CONTEXT.get("PTEN Hamartoma Tumor Syndrome", {}).get("minor_criteria", [])   # subset of conditions
    macrocephaly_plus_criteria = USER_CONTEXT.get("PTEN Hamartoma Tumor Syndrome", {}).get("macrocephaly_plus_criteria", [])  # subset of conditions

    major_count = sum(c in conditions for c in major_criteria)
    minor_count = sum(c in conditions for c in minor_criteria)

    macrocephaly = "Macrocephaly" in conditions
    gi_hamartoma = "Multiple gastrointestinal hamartomas" in conditions

    macrocephaly_plus = (
        macrocephaly and
        any(c in conditions for c in macrocephaly_plus_criteria)
    )

    if (
        macrocephaly_plus or
        (major_count >= 3 and (macrocephaly or gi_hamartoma)) or
        (major_count >= 2 and minor_count >= 3)
    ):
        category = "Significant Pattern"
    elif major_count >= 1 or minor_count >= 2 or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "PTEN Hamartoma Tumor Syndrome", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "PTEN Hamartoma Tumor Syndrome"),
        "Blood Marker(s)": "",
    }]
evaluate_phts_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_phts',
        'description': 'Evaluate PTEN Hamartoma Tumor Syndrome susceptibility - no lab dependency, gated on PTEN genetics. No past_history dependency.',
        'parameters': {
            'type': 'object',
            'properties': {
                'conditions': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('PTEN Hamartoma Tumor Syndrome', {}).get('major_criteria', []) + USER_CONTEXT.get('PTEN Hamartoma Tumor Syndrome', {}).get('minor_criteria', []),
                    },
                    'description': "Prior relevant diagnoses/features this patient has had, from this condition's major and minor criteria lists (macrocephaly_plus_criteria is a subset of minor_criteria, already covered).",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['PTEN Hamartoma Tumor Syndrome'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('PTEN Hamartoma Tumor Syndrome', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('PTEN Hamartoma Tumor Syndrome', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['conditions', 'flagged_genes', 'family_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Gilbert Syndrome
def evaluate_gilbert_syndrome(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Gilbert Syndrome", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Gilbert Syndrome"):
        return [{"Condition": "Gilbert Syndrome", "Category": GENE_NOT_FOUND}]

    total_bilirubin = labs.get("total_bilirubin")
    direct_bilirubin = labs.get("direct_bilirubin")
    indirect_bilirubin = labs.get("indirect_bilirubin")
    alt = labs.get("alt")
    ast = labs.get("ast")
    alp = labs.get("alp")
    ggt = labs.get("ggt")

    triggered = []

    total_bilirubin_early = is_elevated(total_bilirubin, TOTAL_BILIRUBIN_RANGE[0]) and is_below(total_bilirubin, TOTAL_BILIRUBIN_RANGE[1])
    total_bilirubin_elevated = is_elevated(total_bilirubin, TOTAL_BILIRUBIN_RANGE[1]) and is_low(total_bilirubin, TOTAL_BILIRUBIN_RANGE[2])
    if total_bilirubin_elevated:
        note(triggered, "total_bilirubin", total_bilirubin, True, f"{TOTAL_BILIRUBIN_RANGE[1]}-{TOTAL_BILIRUBIN_RANGE[2]}")
    elif total_bilirubin_early:
        note(triggered, "total_bilirubin", total_bilirubin, True, f"{TOTAL_BILIRUBIN_RANGE[0]}-{TOTAL_BILIRUBIN_RANGE[1]}")

    direct_bilirubin_flag = note(triggered, "direct_bilirubin", direct_bilirubin, is_below(direct_bilirubin, DIRECT_BILIRUBIN_LOWER), f"<{DIRECT_BILIRUBIN_LOWER}")
    indirect_bilirubin_flag = note(triggered, "indirect_bilirubin", indirect_bilirubin, is_above(indirect_bilirubin, INDIRECT_BILIRUBIN_UPPER), f">{INDIRECT_BILIRUBIN_UPPER}")
    
    alt_normal = is_low(alt, ALT_UPPER)
    ast_normal = is_low(ast, AST_UPPER)
    alp_normal = is_elevated(alp, ALP_RANGE[0]) and is_low(alp, ALP_RANGE[1])
    ggt_normal = is_elevated(ggt, GGT_LFT_RANGE[0]) and is_low(ggt, GGT_LFT_RANGE[1])

    # normal enzymes only used in logic and not reported
    # if total_bilirubin_early or total_bilirubin_elevated or direct_bilirubin_flag or indirect_bilirubin_flag:
    #     note(triggered, "alt", alt, alt_normal, f"<= {ALT_UPPER}")
    #     note(triggered, "ast", ast, ast_normal, f"<= {AST_UPPER}")
    #     note(triggered, "alp", alp, alp_normal, f"{ALP_RANGE[0]}-{ALP_RANGE[1]}")
    #     note(triggered, "ggt", ggt, ggt_normal, f"{GGT_LFT_RANGE[0]}-{GGT_LFT_RANGE[1]}")

    # noting down logic updated, 
    # trigger logic unchanged — all 7 must be True together
    parameter_elevated = (
        total_bilirubin_elevated and
        direct_bilirubin_flag and
        indirect_bilirubin_flag and
        alt_normal and
        ast_normal and
        alp_normal and
        ggt_normal
    )
    
    parameter_early = (
        total_bilirubin_early and
        direct_bilirubin_flag and
        indirect_bilirubin_flag and
        alt_normal and
        ast_normal and
        alp_normal and
        ggt_normal
    )

    if (
        parameter_elevated or
        (parameter_early and (family_history or symptoms or past_history))
    ):
        category = "Significant Pattern"
    elif (
        family_history or
        parameter_early 
    ):
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Gilbert Syndrome", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Gilbert Syndrome"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_gilbert_syndrome_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_gilbert_syndrome',
        'description': 'Evaluate Gilbert Syndrome susceptibility from bilirubin fractions and liver enzymes, gated on UGT1A1 genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'total_bilirubin': {
                    'type': ['number', 'null'],
                    'description': 'total_bilirubin, mg/dL. Null if not tested.',
                },
                'direct_bilirubin': {
                    'type': ['number', 'null'],
                    'description': 'direct_bilirubin, mg/dL. Null if not tested.',
                },
                'indirect_bilirubin': {
                    'type': ['number', 'null'],
                    'description': 'indirect_bilirubin, mg/dL. Null if not tested.',
                },
                'alt': {
                    'type': ['number', 'null'],
                    'description': 'alt, U/L. Null if not tested.',
                },
                'ast': {
                    'type': ['number', 'null'],
                    'description': 'ast, U/L. Null if not tested.',
                },
                'alp': {
                    'type': ['number', 'null'],
                    'description': 'alp, U/L. Null if not tested.',
                },
                'ggt': {
                    'type': ['number', 'null'],
                    'description': 'ggt, U/L. Null if not tested.',
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Gilbert Syndrome'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Gilbert Syndrome', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Gilbert Syndrome', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Gilbert Syndrome', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['total_bilirubin', 'direct_bilirubin', 'indirect_bilirubin', 'alt', 'ast', 'alp', 'ggt', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Hereditary Hemochromatosis
def evaluate_hemochromatosis(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Hereditary Hemochromatosis", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Hereditary Hemochromatosis"):
        return [{"Condition": "Hereditary Hemochromatosis", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    ferritin = labs.get("ferritin")
    transferrin_saturation = labs.get("transferrin_saturation")

    triggered = []
    ferritin_threshold = sex_based_threshold(sex, FERRITIN_ELEVATED_MALE, FERRITIN_ELEVATED_FEMALE)

    ferritin_elevated = note(triggered, "ferritin", ferritin, is_above(ferritin, ferritin_threshold), f">{ferritin_threshold}")
    tsat_elevated = note(triggered, "transferrin_saturation", transferrin_saturation, is_elevated(transferrin_saturation, TRANSFERRIN_SATURATION_ELEVATED), f">={TRANSFERRIN_SATURATION_ELEVATED}")

    if (ferritin_elevated and tsat_elevated) or ((ferritin_elevated and tsat_elevated) and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
    elif tsat_elevated:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Hereditary Hemochromatosis", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Hereditary Hemochromatosis"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_hemochromatosis_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_hemochromatosis',
        'description': 'Evaluate Hereditary Hemochromatosis susceptibility from ferritin and transferrin saturation, gated on HFE genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'ferritin': {
                    'type': ['number', 'null'],
                    'description': 'ferritin, ng/mL. Null if not tested.',
                },
                'transferrin_saturation': {
                    'type': ['number', 'null'],
                    'description': 'transferrin_saturation, %. Null if not tested.',
                },
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Hereditary Hemochromatosis'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hereditary Hemochromatosis', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hereditary Hemochromatosis', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Hereditary Hemochromatosis', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['ferritin', 'transferrin_saturation', 'sex', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Prostate Cancer
def evaluate_prostate(labs, patient, genetics, family_history, past_history, symptoms=False):
    sex = patient.get("sex")
    if sex != "male":
        return []

    if genetics is None:
        return [{"Condition": "Prostate Cancer", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Prostate Cancer"):
        return [{"Condition": "Prostate Cancer", "Category": GENE_NOT_FOUND}]

    return [{
        "Condition": "Prostate Cancer", "Category": "Elevated susceptibility", #hard-coded, update when evaluating E.P. and S.P.
        "DNA Marker(s)": get_gene_trigger(genetics, "Prostate Cancer"),
        "Blood Marker(s)": "", #PSA not evaluated, so hard-coded empty for now, update when evaluating based on reference update from lab
    }]
evaluate_prostate_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_prostate',
        'description': "Evaluate Prostate Cancer susceptibility - no lab dependency, gated on a 13-gene prostate cancer panel. Only evaluated for male patients (returns no result for others)",
        'parameters': {
            'type': 'object',
            'properties': {
                'sex': {
                    'type': 'string',
                    'enum': ['male', 'female'],
                    'description': "Patient's biological sex.",
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Prostate Cancer'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Prostate Cancer', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Prostate Cancer', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Prostate Cancer', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['sex', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# ttg-iga added back to blood panel

# # Celiac Disease
def evaluate_celiac_disease(labs, patient, genetics, family_history, past_history, symptoms=False):
    if genetics is None:
        return [{"Condition": "Celiac Disease", "Category": INSUFFICIENT_DATA}]
    if not has_flagged_gene(genetics, "Celiac Disease"):
        return [{"Condition": "Celiac Disease", "Category": GENE_NOT_FOUND}]

    ttg = labs.get("ttg")

    TTG_Early = is_elevated(ttg, TTG_RANGE[0])
    TTG_Likely = is_above(ttg, TTG_RANGE[1])

    triggered = []
    if TTG_Likely:
            note(triggered, "ttg", ttg, True, f">{TTG_RANGE[1]}")
    elif TTG_Early:
        note(triggered, "ttg", ttg, True, f">={TTG_RANGE[0]}")
    

    if TTG_Likely or (TTG_Early and (family_history or symptoms or past_history)):
        category = "Significant Pattern"
    elif TTG_Early:
        category = "Early Pattern"
    else:
        category = "Elevated susceptibility"

    return [{
        "Condition": "Celiac Disease", "Category": category,
        "DNA Marker(s)": get_gene_trigger(genetics, "Celiac Disease"),
        "Blood Marker(s)": ", ".join(triggered),
    }]
evaluate_celiac_disease_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_celiac_disease',
        'description': 'Evaluate Celiac Disease susceptibility from tissue transglutaminase IgA (ttg-iga) levels, gated on HLA-DQ2/DQ8 genetics.',
        'parameters': {
            'type': 'object',
            'properties': {
                'ttg': {
                    'type': ['number', 'null'],
                    'description': 'ttg, AU/mL. Null if not tested.',
                },
                'flagged_genes': {
                    'type': ['array', 'null'],
                    'items': {
                        'type': 'string',
                        'enum': CONDITION_GENES['Celiac Disease'],
                    },
                    'description': 'Which of these genes are flagged pathogenic/likely-pathogenic for this patient. Null if genetics is entirely unavailable (not just none of these genes).',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Celiac Disease', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Celiac Disease', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Celiac Disease', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['ttg', 'flagged_genes', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}

# Apply with_enrichment to every live evaluate_* in one place - this is the
# single choke point that guarantees Domain/Triggering PRS/Snapshot Category/
# RSID_m/Pathogenicity_m/APOE_Status/Elevation_p/Threshold_p/All Blood
# Marker(s) are always present on every evaluate_* call's own return value,
# without touching any of the 23 functions' bodies above. Add a new
# condition's function name here when it's enabled.
for _evaluator_name in [
    "evaluate_alzheimers", "evaluate_atopic_dermatitis", "evaluate_eczema", "evaluate_copd",
    "evaluate_hyperthyroidism", "evaluate_hypothyroidism", "evaluate_ibd", "evaluate_nafld",
    "evaluate_psoriasis", "evaluate_rheumatoid_arthritis", "evaluate_t2d", "evaluate_mody",
    "evaluate_muscular_dystrophy", "evaluate_fh", "evaluate_cardiomyopathy", "evaluate_cad",
    "evaluate_hypertriglyceridemia", "evaluate_hdl_deficiency", "evaluate_li_fraumeni_syndrome",
    "evaluate_phts", "evaluate_gilbert_syndrome", "evaluate_hemochromatosis", "evaluate_prostate", "evaluate_celiac_disease"
]:
    globals()[_evaluator_name] = with_enrichment(globals()[_evaluator_name])
del _evaluator_name



# dead weight CLI approach
# calling: CLI/batch entry point only - importing this module for its
# individual evaluate_* functions (e.g. as model tool-calls) never runs
# any of this.

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Run genetic risk evaluation across all conditions")
#     parser.add_argument("monogenic_json", help="Path to the monogenic pathogenic/likely-pathogenic variants JSON")
#     parser.add_argument("prs_json", help="Path to the polygenic PRS scores JSON")
#     parser.add_argument("apoe_json", help="Path to the APOE status JSON")
#     parser.add_argument("blood_csv", help="Path to the blood report CSV")
#     parser.add_argument("output_csv", help="Path to write the results CSV")
#     parser.add_argument("--sex", required=False, choices=["male", "female"], default=None, help="Patient's biological sex (optional override on csv extracted value)") # making it optional override
#     parser.add_argument("--age", required=False, type=int, default=None, help="Patient's age in years (optional override on csv extracted value)") #making it optional override
#     parser.add_argument("--thresholds", required=True, help="Path to the lab-specific reference-ranges/thresholds JSON") # added
#     args = parser.parse_args()
# 
#     load_thresholds(args.thresholds)  # override the default thresholds loaded at import time
# 
#     data = load_patient_data(args.monogenic_json, args.prs_json, args.apoe_json, args.blood_csv, sex=args.sex, age=args.age) #remains same, these args remain override call.
# 
#     results = [
#         # --- 14: labs/genetics-driven conditions ---
#         evaluate_alzheimers(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         # evaluate_asthma(
#         #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#         #     family_history=False,
#         #     past_history=False,
#         # ),
#         evaluate_atopic_dermatitis(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_eczema(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_copd(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_hyperthyroidism(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_hypothyroidism(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_ibd(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_nafld(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         # evaluate_osteoarthritis(
#         #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#         #     family_history=False,
#         #     past_history=False,
#         # ),
#         # evaluate_parkinsons(
#         #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#         #     family_history=False,
#         #     past_history=False,
#         # ),
#         evaluate_psoriasis(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_rheumatoid_arthritis(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         # evaluate_rhinitis(
#         #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#         #     family_history=False,
#         #     past_history=False,
#         # ),
#         evaluate_t2d(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
# 
#         # --- 14: second batch ---
#         evaluate_mody(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#         ),
#         # evaluate_men4(
#         #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#         #     family_history=False, endocrine_conditions=None,
#         # ),
#         evaluate_muscular_dystrophy(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             symptoms=False,
#         ),
#         evaluate_gilbert_syndrome(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#             symptoms=False,
#         ),
#         evaluate_hemochromatosis(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#             symptoms=False,
#         ),
#         # evaluate_celiac_disease(
#         #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#         #     family_history=False, past_history=False, symptoms=False,
#         # ),
#         evaluate_fh(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_cardiomyopathy(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             symptoms=False,
#             imaging_performed=False, imaging_findings=None,
#             past_history=False,
#         ),
#         evaluate_cad(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             symptoms=False,
#             past_history=False,
#         ),
#         evaluate_hypertriglyceridemia(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         evaluate_hdl_deficiency(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             past_history=False,
#         ),
#         # evaluate_apob(
#         #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#         #     family_history=False,
#         # ),
#         evaluate_li_fraumeni_syndrome(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             conditions=None,
#         ),
#         evaluate_phts(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             conditions=None,
#         ),
#         evaluate_prostate(
#             labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
#             family_history=False,
#             symptoms=False,
#             past_history=False,
#         ),
#     ]
# 
#     CONDITION_DOMAINS = {
#         "Alzheimer's Disease": "Brain Health",
#         "Asthma": "Respiratory Health",
#         "Atopic Dermatitis": "Skin Health",
#         "Eczema": "Skin Health",
#         "COPD": "Respiratory Health",
#         "Hyperthyroidism": "Endocrine Health",
#         "Hypothyroidism": "Endocrine Health",
#         "Inflammatory Bowel Disease": "Digestive Health",
#         "NAFLD": "Digestive Health",
#         "Osteoarthritis": "Joints Health",
#         "Parkinson's Disease": "Brain Health",
#         "Psoriasis": "Skin Health",
#         "Rheumatoid Arthritis": "Immune Health",
#         "Rhinitis": "Respiratory Health",
#         "Type 2 Diabetes": "Endocrine Health",
#         "Maturity-Onset Diabetes of the Young": "Endocrine Health",
#         # "Multiple Endocrine Neoplasia Type 4": "Endocrine Health",
#         "Muscular Dystrophy": "Muscle Health",
#         "Gilbert Syndrome": "Kidney & Metabolic Health",
#         "Hereditary Hemochromatosis": "Digestive Health",
#         # "Celiac Disease": "Immune Health",
#         "Familial Hypercholesterolemia": "Cardiac Health",
#         "Cardiomyopathy": "Cardiac Health",
#         "Coronary Artery Disease": "Cardiac Health",
#         "Familial Hypertriglyceridemia": "Cardiac Health",
#         "HDL Deficiency": "Cardiac Health",
#         "Li-Fraumeni Syndrome": "Cancer",
#         "PTEN Hamartoma Tumor Syndrome": "Cancer",
#         "Prostate Cancer": "Cancer"
#     }
# 
#     active_findings = [r[0] for r in results if r]
# 
#     for entry in active_findings:
#         entry["Domain"] = CONDITION_DOMAINS.get(entry["Condition"], "")
#         entry["Triggering PRS"] = get_triggering_prs(entry["Condition"], data["genetics"])
#         entry["Snapshot Category"] = get_snapshot_category(entry["Condition"], entry["Category"], data["genetics"])
#         rsid_m, patho_m = get_monogenic_variant(data["genetics"], entry["Condition"])
#         entry["RSID_m"] = rsid_m
#         entry["Pathogenicity_m"] = patho_m
#         entry["APOE_Status"] = get_apoe_status(data["genetics"], entry["Condition"])
# 
#         elevation = get_prs_category(data["genetics"], entry["Condition"]) or ""
#         entry["Elevation_p"] = elevation
#         entry["Threshold_p"] = (
#             ">99 percentile" if elevation == "elevated"
#             else ">97 percentile" if elevation in {"moderately elevated", "moderately_elevated"}
#             else ""
#         )
#     
#         entry["All Blood Marker(s)"] = entry.get("Blood Marker(s)", "")
#         if entry.get("Category") not in PATTERN_CATEGORIES_INDICATING_TRIGGER:
#             entry["Blood Marker(s)"] = ""
# 
# 
#     with open(args.output_csv, "w", newline="") as f:
#         writer = csv.DictWriter(f, fieldnames=["Domain", "Condition", "Category", "DNA Marker(s)", "Blood Marker(s)", "Triggering PRS", "Snapshot Category", "All Blood Marker(s)", "RSID_m", "Pathogenicity_m", "Elevation_p", "Threshold_p", "APOE_Status"])
#         writer.writeheader()
#         writer.writerows(active_findings)
# 
#     print(f"Results written to {args.output_csv}")
