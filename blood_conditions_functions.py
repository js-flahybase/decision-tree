import argparse
import json
import csv

# example run command:
# python script.py monogenic.json prs.json apoe.json blood.csv results.csv --sex male --age 35

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run genetic risk evaluation across all conditions")
    parser.add_argument("monogenic_json", help="Path to the monogenic pathogenic/likely-pathogenic variants JSON")
    parser.add_argument("prs_json", help="Path to the polygenic PRS scores JSON")
    parser.add_argument("apoe_json", help="Path to the APOE status JSON")
    parser.add_argument("blood_csv", help="Path to the blood report CSV")
    parser.add_argument("output_csv", help="Path to write the results CSV")
    parser.add_argument("--sex", required=True, choices=["male", "female"], help="Patient's biological sex") ##
    parser.add_argument("--age", required=True, type=int, help="Patient's age in years") ##
    args = parser.parse_args()

# loading the reference user-context json (can also be done using args (but since this doesn't change with samples, so hard-coded using))
# with open("/home/azureuser/decision-tree/blood_conditions_user_context.json") as f:
#     user_context = json.load(f)

# genentics part from result jsons 
def build_genetics_from_jsons(monogenic_json_path, prs_json_path, apoe_json_path):
    """
    Reads the three JSON files (monogenic pathogenic variants, polygenic PRS
    scores, APOE status) and returns a single genetics dict:
        {
            "flagged_genes": [...],              # monogenic pathogenic/likely pathogenic hits (+ APOE if e4 present)
            "prs_elevated_conditions": {...},     # set of PRS json's condition keys where category is elevated/moderately_elevated
            "apoe_status": "e3/e4",               # raw APOE genotype string, for reference/future use
        }
    """
    PATHOGENIC_VALUES = {"pathogenic", "likely pathogenic", "likely_pathogenic", "pathogenic/likely pathogenic", "pathogenic/likely_pathogenic"}
    ELEVATED_CATEGORIES = {"elevated", "moderately_elevated", "moderately elevated"}

    # --- 1. monogenic pathogenic/likely pathogenic gene hits ---
    flagged_genes = []
    acmg_genes = []  # subset of flagged_genes
    with open(monogenic_json_path) as f:
        monogenic_data = json.load(f)
    for entry in monogenic_data:
        sig = (entry.get("ClinVar Significance") or "").strip().lower()
        if sig in PATHOGENIC_VALUES:
            gene = (entry.get("Gene Symbol") or "").strip()
            if gene:
                flagged_genes.append(gene)
                if entry.get("in_acmg") is True:
                    acmg_genes.append(gene)

    # --- 2. Polygenic PRS elevated conditions ---
    prs_elevated_conditions = set()
    prs_categories = {}
    with open(prs_json_path) as f:
        prs_data = json.load(f)
    for condition_key, info in prs_data.items():
        category = (info.get("category") or "").strip().lower()
        prs_categories[condition_key] = category
        if category in ELEVATED_CATEGORIES:
            prs_elevated_conditions.add(condition_key)

    # --- 3. APOE status ---
    apoe_status = ""
    with open(apoe_json_path) as f:
        apoe_data = json.load(f)
    apoe_status = (apoe_data.get("APOE_Status") or "").strip()
    if "4" in apoe_status:  # e.g. "e3/e4" or "e4/e4" -> ε4 allele present
        flagged_genes.append("APOE")
        acmg_genes.append("APOE") # since apoe is coming from a different json, e4 presence currently make it to acmg list (change approach if needed)

    return {
        "flagged_genes": flagged_genes,
        "acmg_genes": acmg_genes,
        "prs_elevated_conditions": prs_elevated_conditions,
        "prs_categories": prs_categories,
        "apoe_status": apoe_status,
    }

# match prs json syntax with condition's names used here
PRS_CONDITION_KEY_MAP = {
    "alzheimers": "Alzheimer's Disease",
    "cad": "Coronary Artery Disease",
    "asthma": "Asthma",
    "atopicdermatitis": "Atopic Dermatitis/Eczema",
    "eczema": "Atopic Dermatitis/Eczema",
    "copd": "COPD",
    "hyperthyroidism": "Hyperthyroidism",
    "hypothyroidism": "Hypothyroidism",
    "ibd": "Inflammatory Bowel Disease",
    "nafld": "NAFLD",
    "osteoarthritis": "Osteoarthritis",
    "parkinson": "Parkinson's Disease",
    "psoriasis": "Psoriasis",
    "rheumatoidarthritis": "Rheumatoid Arthritis",
    "rhinitis": "Rhinitis",
    "type2diabetes": "Type 2 Diabetes",
}

# combining genentic (gene, prs), blood (labs), and patient info (sex, age)
def load_patient_data(monogenic_json, prs_json, apoe_json, blood_csv_path, sex, age):
    """
    Patient-reported context (family_history, symptoms, past_history) is 
    NOT included here
    """
    # -----------------------------
    # Blood report
    # -----------------------------
    labs = {}
    with open(blood_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row["parameter"].strip().lower()
            value = row["value"].strip()
            labs[key] = float(value) if value else None

    # Derive iron_profile_done: True only if BOTH ferritin and
    # transferrin_saturation were present (non-null) in the blood CSV (for evaluate_hemochromatosis).
    labs["iron_profile_done"] = (
        labs.get("ferritin") is not None and labs.get("transferrin_saturation") is not None
    )
    
    genetics = build_genetics_from_jsons(monogenic_json, prs_json, apoe_json)
    
    patient = {"sex": sex, "age": age}

    # just to check
    print(f"Patient: {patient}")
    print(f"Genetics: {genetics}")
    print(f"Labs: {labs}") #would also add an "iron_profile_done" entry

    return {
        "patient": patient,
        "genetics": genetics,
        "labs": labs,
    }

# =============================================================================
# Shared lab reference ranges / thresholds
# =============================================================================

SEX_VALUES = ["male", "female"]

# GENETIC_RISK_CLASSIFICATIONS = ["Elevated", "Moderately Elevated", "Typical", "Not Reported"]
# ELEVATED_CLASSIFICATIONS = {"Elevated", "Moderately Elevated"}

# --- Lipid Profile ---
TOTAL_CHOLESTEROL_AT_RISK = 200       # mg/dL
TOTAL_CHOLESTEROL_ELEVATED = 240      # mg/dL
LDL_C_NORMAL_CAD = 70                 # mg/dL, preferred ceiling if CAD susceptibility present
LDL_C_AT_RISK = 130                   # mg/dL
LDL_C_ELEVATED = 160                  # mg/dL
LDL_C_SEVERE = 190                    # mg/dL
HDL_C_MALE_MIN = 40                   # mg/dL
HDL_C_FEMALE_MIN = 50                 # mg/dL
HDL_C_MALE_SEVERE = 30                # mg/dL
HDL_C_FEMALE_SEVERE = 35              # mg/dL
NON_HDL_C_NORMAL = 130                # mg/dL
NON_HDL_C_ELEVATED = 160              # mg/dL
NON_HDL_C_SEVERE = 190                # mg/dL
VLDL_UPPER_NORMAL = 30                # mg/dL
TRIGLYCERIDES_AT_RISK = 150           # mg/dL
TRIGLYCERIDES_ELEVATED = 200          # mg/dL
TRIGLYCERIDES_SEVERE = 500            # mg/dL
TRIGLYCERIDES_VERY_SEVERE = 1000      # mg/dL
# TG_HDL_RATIO_ELEVATED = 2.5
APOB_AT_RISK = 100                    # mg/dL
APOB_ELEVATED = 130                   # mg/dL
LPA_AT_RISK = 50                      # mg/dL
LPA_ELEVATED = 100                    # mg/dL

# --- Glucose / Metabolic ---
FASTING_GLUCOSE_AT_RISK = 100         # mg/dL
FASTING_GLUCOSE_ELEVATED = 126        # mg/dL
HBA1C_AT_RISK = 5.7                   # %
HBA1C_ELEVATED = 6.5                  # %
FASTING_INSULIN_RANGE = (2, 25)       # uIU/mL

# --- Inflammatory Markers ---
CRP_AT_RISK = 0.8                     # mg/L
CRP_AT_RISK_COPD = 1.0                # mg/L
CRP_ELEVATED = 3.0                    # mg/L
IGE_UPPER_NORMAL = 100                # IU/mL
FIBRINOGEN_RANGE = (200, 400)         # mg/dL
# FECAL_CALPROTECTIN_AT_RISK = 50       # ug/g
# FECAL_CALPROTECTIN_ELEVATED = 200     # ug/g

def esr_upper_normal(sex, age):
    """Returns ESR upper limit of normal (mm/hr) by sex and age."""
    if sex not in SEX_VALUES or age is None:
        return None
    
    if sex == "male":
        return 16 if age < 50 else 21
    else:
        return 21 if age < 50 else 31

# --- CBC ---
EOSINOPHILS_UPPER_NORMAL = 0.4        # x10^9/L
EOSINOPHILS_BORDERLINE = 0.38         # x10^9/L
EOSINOPHILS_LIKELY = 0.42             # x10^9/L
EOSINOPHILS_ELEVATED = 0.5            # x10^9/L
NEUTROPHILS_RANGE = (1.8, 7.5)        # x10^9/L
PLATELETS_RANGE = (150, 400)          # x10^9/L
NLR_RANGE = (3, 6)                    # ratio
HEMOGLOBIN_MALE_RANGE = (13, 18)      # g/dL
HEMOGLOBIN_FEMALE_RANGE = (12, 16)    # g/dL

# --- Liver Function / Fibrosis ---
AST_RANGE = (10, 40)                  # U/L
ALT_RANGE = (10, 40)                  # U/L
GGT_RANGE = (9, 50)                   # U/L
AST_ALT_RATIO_RANGE = (0.7, 1.2)
FIB4_AT_RISK = 1.3
FIB4_ELEVATED = 2.67
ALBUMIN_RANGE = (3.5, 5.5)            # g/dL
# FIBROSIS_SCORE_ELEVATED_THRESHOLD = None  # "fibrosis scores" referenced generically (TM6SF2 rule), no defined metric/cutoff

# --- Thyroid ---
TSH_RANGE = (0.4, 4.0)                # mIU/L
TSH_BORDERLINE = 0.45                 # mIU/L
TSH_ELEVATED_RANGE = (4.1, 10.0)      # mIU/L
FREE_T4_RANGE = (0.8, 1.8)            # ng/dL
FREE_T3_RANGE = (2.3, 4.4)            # pg/mL
TPOAB_UPPER_NORMAL = 35               # IU/mL

# --- Iron studies (HFE) ---
# FERRITIN_ELEVATED_THRESHOLD = None  # prompt says "ferritin abnormal", no numeric cutoff given
# TRANSFERRIN_SATURATION_ELEVATED_THRESHOLD = None  # same as above, no numeric cutoff given
# SERUM_IRON_ELEVATED_THRESHOLD = None  # not referenced with a number anywhere in prompt, included for completeness if needed

CALCIUM_UPPER_NORMAL = 10.2                    # mg/dL
PTH_UPPER_NORMAL = 65                          # pg/mL
CK_UPPER_NORMAL = 200                          # U/L

TOTAL_BILIRUBIN_RANGE = (1.2, 4)               # mg/dL (Gilbert)
DIRECT_BILIRUBIN_LOWER = 0.3                   # mg/dL (Gilbert)
INDIRECT_BILIRUBIN_UPPER = 5                   # mg/dL (Gilbert)
AST_UPPER_STRICT = 41                          # U/L (Gilbert, Celiac, Hemochromatosis)
ALP_RANGE = (40, 129)                          # U/L (Gilbert, Celiac)
GGT_LFT_RANGE = (10, 47)                       # U/L (Gilbert, Celiac)
FERRITIN_ELEVATED = 300                        # ng/mL (Hemochromatosis)
TRANSFERRIN_SATURATION_ELEVATED = 45           # % (Hemochromatosis)
HEMOGLOBIN_UPPER_HEMOCHROMATOSIS = 18          # g/dL
MCV_UPPER_HEMOCHROMATOSIS = 100                # fL
MCH_UPPER_HEMOCHROMATOSIS = 34                 # pg
# VITAMIN_D_DEFICIENT = 20                       # ng/mL (Celiac)
# MCV_LOWER_CELIAC = 82                          # fL
# MCH_LOWER_CELIAC = 27                          # pg
# MCHC_LOWER_CELIAC = 31.5                       # g/dL
TTG_RANGE = (20,30)                              # U/mL (Celiac)

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
    "Asthma": ["IL13"],
    "Atopic Dermatitis/Eczema": ["C11orf32/LRRC32", "TSLP", "IL13", "IL4R"],
    "COPD": [],  # no gene listed, only PRS gate
    "Hyperthyroidism": ["CTLA4", "PTPN22"],
    "Hypothyroidism": ["CTLA4"],
    "Inflammatory Bowel Disease": ["NOD2", "ATG16L1", "IL23R", "MYC/POU5F1B", "IL27"],
    "NAFLD": ["PNPLA3", "TM6SF2"],
    "Osteoarthritis": ["COL11A1", "SMAD3"],
    "Parkinson's Disease": ["LRRK2", "GBA"],
    "Psoriasis": [],  # no gene listed, only PRS gate
    "Rheumatoid Arthritis": ["PTPN22", "STAT4", "TNFAIP3", "CTLA4"],
    "Rhinitis": ["TSLP"],
    "Type 2 Diabetes": ["TCF7L2", "SLC30A8"],#------------------------------------------------------------------
    "Maturity-Onset Diabetes of the Young": ["HNF1A"],
    "Multiple Endocrine Neoplasia Type 4": ["CDKN1B"],
    "Muscular Dystrophy": ["DMD", "FKRP", "LMNA", "SGCA", "SGCB"],
    "Gilbert Syndrome": ["UGT1A1"],
    "Hereditary Hemochromatosis": ["HFE"],
    "Celiac Disease": ["HLA-DQ2", "HLA-DQ8"],
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
    "Elevated Apolipoprotein B": ["APOB"],
    "Li-Fraumeni Syndrome": ["TP53"],
    "PTEN Hamartoma Tumor Syndrome": ["PTEN"],
}

# only conditions filtered through this will have their functions run on patient data block
def has_flagged_gene(genetics, condition_name):
    """Returns True if:
    - the condition is exempt from the gene gate (NO_GENE_GATE), OR
    - at least one of its required genes is flagged (monogenic hit), OR
    - the condition has a polygenic PRS entry AND that PRS category is
      elevated/moderately elevated.
    Returns False (blocked) otherwise."""
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
    required_genes = CONDITION_GENES.get(condition_name, [])
    parts = []

    if required_genes != NO_GENE_GATE and required_genes:
        flagged = genetics.get("flagged_genes", [])
        parts.extend(gene for gene in required_genes if gene in flagged)

    prs_key = None
    for key, mapped_name in PRS_CONDITION_KEY_MAP.items():
        if mapped_name == condition_name:
            prs_key = key
            break
    if prs_key is not None and prs_key in genetics.get("prs_elevated_conditions", set()):
        parts.append(f"PRS:{prs_key}")

    return ", ".join(parts)

def is_gene_flagged(genetics, condition_name):
    required_genes = CONDITION_GENES.get(condition_name, [])
    if not required_genes or required_genes == NO_GENE_GATE:
        return False
    flagged = genetics.get("flagged_genes", [])
    return any(gene in flagged for gene in required_genes)

def is_gene_acmg(genetics, condition_name):
    required_genes = CONDITION_GENES.get(condition_name, [])
    if not required_genes or required_genes == NO_GENE_GATE:
        return False
    acmg = genetics.get("acmg_genes", [])
    return any(gene in acmg for gene in required_genes)

def get_prs_category(genetics, condition_name):
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
    "elevated": "97.5",
    "moderately_elevated": "95",
    "moderately elevated": "95",
}

def get_triggering_prs(condition_name, genetics):
    """Returns the display value for Triggering PRS if this condition's PRS
    category is elevated/moderately elevated, else ''."""
    prs_category = get_prs_category(genetics, condition_name)
    return PRS_DISPLAY_VALUE.get(prs_category, "")

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
# =============================================================================
# (14)
# =============================================================================

# Family history and symptoms of respective conditions, are mentioned in the deck for ref

GENE_NOT_FOUND = "Typical"

# Alzheimer's Disease
def evaluate_alzheimers(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Alzheimer's Disease"):
        return [{"Condition": "Alzheimer's Disease", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

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
        note(triggered, "total_c", total_c, True, f">={TOTAL_CHOLESTEROL_ELEVATED}")
    elif total_c_early:
        note(triggered, "total_c", total_c, True, f">={TOTAL_CHOLESTEROL_AT_RISK}")

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
        (dyslipidemia_flag_likely and crp_likely and family_history) or
        (dyslipidemia_flag_likely and crp_likely and symptoms) or
        (dyslipidemia_flag_likely and crp_likely and family_history and symptoms)
    ):
        category = "Significant Pattern"
    elif dyslipidemia_flag_early and crp_early:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Alzheimer's Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Alzheimer's Disease"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Asthma
def evaluate_asthma(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Asthma"):
        return [{"Condition": "Asthma", "Category": GENE_NOT_FOUND}]

    eosinophils = labs.get("eosinophils")
    neutrophils = labs.get("neutrophils")
    triggered = []

    eosinophil_likely = is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
    eosinophil_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
    if eosinophil_likely:
        note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_UPPER_NORMAL}")
    elif eosinophil_early:
        note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")

    neutrophil_flag = note(triggered, "neutrophils", neutrophils, is_above(neutrophils, NEUTROPHILS_RANGE[0]), f">{NEUTROPHILS_RANGE[0]}")

    if (eosinophil_likely and neutrophil_flag) or (eosinophil_likely and neutrophil_flag and (symptoms or family_history)):
        category = "Significant Pattern"
    elif eosinophil_early:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Asthma", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Asthma"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Atopic Dermatitis/Eczema
def evaluate_atopic_dermatitis(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Atopic Dermatitis/Eczema"):
        return [{"Condition": "Atopic Dermatitis/Eczema", "Category": GENE_NOT_FOUND}]

    eosinophils = labs.get("eosinophils")
    triggered = []

    eosinophil_likely = is_above(eosinophils, EOSINOPHILS_LIKELY)
    eosinophil_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_LIKELY)
    if eosinophil_likely:
        note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_LIKELY}")
    elif eosinophil_early:
        note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")

    if eosinophil_likely or (eosinophil_likely and (symptoms or family_history)):
        category = "Significant Pattern"
    elif eosinophil_early:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Atopic Dermatitis/Eczema", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Atopic Dermatitis/Eczema"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Chronic Obstruction of Pulmonary Disorder​
def evaluate_copd(labs, patient, genetics, family_history, symptoms=False):
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

    if (eosinophil_flag_likely and neutrophil_flag and fibrinogen_flag and crp_elevated_flag) or (eosinophil_flag_likely and neutrophil_flag and fibrinogen_flag and crp_elevated_flag and (symptoms or family_history)):
        category = "Significant Pattern"
    elif (crp_at_risk_flag and eosinophil_flag_early):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "COPD", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "COPD"),
        "Triggering Parameters": ", ".join(triggered),
    }]

# to be updated with thresholds once confirmed
# Hyperthyroidism
def evaluate_hyperthyroidism(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Hyperthyroidism"):
        return [{"Condition": "Hyperthyroidism", "Category": GENE_NOT_FOUND}]

    tsh = labs.get("tsh")
    free_t3 = labs.get("free_t3")
    free_t4 = labs.get("free_t4")
    tpoab = labs.get("tpoab")

    triggered = []
    tsh_borederline_suppressed = is_below(tsh, TSH_BORDERLINE)
    tsh_suppressed = is_below(tsh, TSH_RANGE[0])
    if tsh_suppressed:
        note(triggered, "tsh", tsh, True, f"<{TSH_RANGE[0]}")
    elif tsh_borederline_suppressed:
        note(triggered, "tsh", tsh, True, f"<{TSH_BORDERLINE}")

    hormone_flag = (
        note(triggered, "free_t3", free_t3, is_above(free_t3, FREE_T3_RANGE[0]), f">{FREE_T3_RANGE[0]}")
        or note(triggered, "free_t4", free_t4, is_above(free_t4, FREE_T4_RANGE[0]), f">{FREE_T4_RANGE[0]}")
    )
    tpoab_flag = note(triggered, "tpoab", tpoab, is_elevated(tpoab, TPOAB_UPPER_NORMAL), f">={TPOAB_UPPER_NORMAL}")

    if (tsh_suppressed and hormone_flag and tpoab_flag) or (tsh_suppressed and hormone_flag and tpoab_flag and (symptoms or family_history)):
        category = "Significant Pattern"
    elif tsh_borederline_suppressed and (hormone_flag or tpoab_flag):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Hyperthyroidism", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Hyperthyroidism"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Hypothyroidism
def evaluate_hypothyroidism(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Hypothyroidism"):
        return [{"Condition": "Hypothyroidism", "Category": GENE_NOT_FOUND}]

    tsh = labs.get("tsh")
    free_t4 = labs.get("free_t4")
    tpoab = labs.get("tpoab")

    triggered = []
    tsh_high = is_above(tsh, TSH_ELEVATED_RANGE[1])
    tsh_mild = tsh is not None and TSH_ELEVATED_RANGE[0] <= tsh <= TSH_ELEVATED_RANGE[1]
    if tsh_high:
        note(triggered, "tsh", tsh, True, f">{TSH_ELEVATED_RANGE[1]}")
    elif tsh_mild:
        note(triggered, "tsh", tsh, True, f"[{TSH_ELEVATED_RANGE[0]}-{TSH_ELEVATED_RANGE[1]}]")
    free_t4_low = note(triggered, "free_t4", free_t4, is_below(free_t4, FREE_T4_RANGE[0]), f"<{FREE_T4_RANGE[0]}")
    tpoab_flag = note(triggered, "tpoab", tpoab, is_elevated(tpoab, TPOAB_UPPER_NORMAL), f">={TPOAB_UPPER_NORMAL}")

    if (tsh_high and free_t4_low and tpoab_flag):
        category = "Significant Pattern"
    # tsh_mild and tpoab_flag case covers early pattern w/o user context
    elif tsh_mild and (tpoab_flag or family_history):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Hypothyroidism", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Hypothyroidism"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Inflammatory Bowel Disease
def evaluate_ibd(labs, patient, genetics, family_history, symptoms=False):
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
    if (crp_elevated_flag and anemia_flag and albumin_flag) or (crp_elevated_flag and anemia_flag and albumin_flag and (symptoms or family_history)):
        category = "Significant Pattern"
    elif crp_at_risk_flag:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Inflammatory Bowel Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Inflammatory Bowel Disease"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Non alcoholic fatty liver disease
def evaluate_nafld(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "NAFLD"):
        return [{"Condition": "NAFLD", "Category": GENE_NOT_FOUND}]

    ast = labs.get("ast")
    alt = labs.get("alt")
    ggt = labs.get("ggt")
    fib4 = labs.get("fib4")
    platelets = labs.get("platelets")
    ast_alt_ratio = labs.get("ast/alt")

    triggered = []
    lft_flag = (
        note(triggered, "ast", ast, is_above(ast, AST_RANGE[0]), f">{AST_RANGE[0]}")
        or note(triggered, "alt", alt, is_above(alt, ALT_RANGE[0]), f">{ALT_RANGE[0]}")
        or note(triggered, "ggt", ggt, is_above(ggt, GGT_RANGE[0]), f">{GGT_RANGE[0]}")
        or note(triggered, "ast/alt", ast_alt_ratio, is_above(ast_alt_ratio, AST_ALT_RATIO_RANGE[0]), f">{AST_ALT_RATIO_RANGE[0]}")
    )
    fib4_elevated = note(triggered, "fib4", fib4, is_elevated(fib4, FIB4_ELEVATED), f">={FIB4_ELEVATED}")
    fib4_at_risk = note(triggered, "fib4", fib4, fib4 is not None and FIB4_AT_RISK <= fib4 < FIB4_ELEVATED, f"[{FIB4_AT_RISK}-{FIB4_ELEVATED})")
    platelets_low = note(triggered, "platelets", platelets, is_below(platelets, PLATELETS_RANGE[0]), f"<{PLATELETS_RANGE[0]}")

    if (lft_flag and fib4_elevated and platelets_low) or (lft_flag and fib4_elevated and platelets_low and (symptoms or family_history)):
        category = "Significant Pattern"
    elif (lft_flag or fib4_at_risk):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "NAFLD", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "NAFLD"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Osteoarthritis
# only based on user context, will assign "early pattern" only based on PRS/gene and significant if user context present.
def evaluate_osteoarthritis(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Osteoarthritis"):
        return [{"Condition": "Osteoarthritis", "Category": GENE_NOT_FOUND}]

    if symptoms or family_history:
        category = "Significant Pattern"
    else:
        category = "Elevated Susceptibility"
    # else:
    #     category = "Typical" #no typical, blood parameters removed and only PRS and user context is used.

    return [{
        "Condition": "Osteoarthritis", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Osteoarthritis"),
        "Triggering Parameters": "",
    }]


# Parkinson's Disease
def evaluate_parkinsons(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Parkinson's Disease"):
        return [{"Condition": "Parkinson's Disease", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

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
        note(triggered, "total_c", total_c, True, f">={TOTAL_CHOLESTEROL_ELEVATED}")
    elif total_c_early:
        note(triggered, "total_c", total_c, True, f">={TOTAL_CHOLESTEROL_AT_RISK}")

    # --- Triglycerides ---
    tg_likely = is_elevated(triglycerides, TRIGLYCERIDES_ELEVATED)
    tg_early = is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK)
    if tg_likely:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_ELEVATED}")
    elif tg_early:
        note(triggered, "triglycerides", triglycerides, True, f">={TRIGLYCERIDES_AT_RISK}")

    # --- HDL (single threshold, no early/likely tiers) ---
    hdl_low = note(triggered, "hdl_c", hdl_c, is_below(hdl_c, hdl_min), f"<{hdl_min}")

    # dyslipidemia_flag_early = ldl_early or total_c_early or tg_early or hdl_low
    # dyslipidemia_flag_likely = ldl_likely or total_c_likely or tg_likely or hdl_low
    
    # manually put down for now, will just report elevated susceptibility untill updated
    dyslipidemia_flag_early = False
    dyslipidemia_flag_likely = False

    crp_early = is_elevated(crp, CRP_AT_RISK) #and not is_above(crp, CRP_ELEVATED)
    crp_likely = is_above(crp, CRP_ELEVATED)

    if crp_likely:
        note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
    elif crp_early:
        note(triggered, "crp", crp, True, f">={CRP_AT_RISK}")

    if (dyslipidemia_flag_likely and crp_likely) or (dyslipidemia_flag_likely and crp_likely and (family_history or symptoms)):
        category = "Significant Pattern"
    elif dyslipidemia_flag_early and crp_early:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Parkinson's Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Parkinson's Disease"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Psoriasis
def evaluate_psoriasis(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Psoriasis"):
        return [{"Condition": "Psoriasis", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    age = patient.get("age")
    esr_max = esr_upper_normal(sex, age) #as per the prompt (children/newborn age case not defined)

    crp = labs.get("crp")
    esr = labs.get("esr")
    nlr = labs.get("nlr")

    triggered = []
    crp_early = is_elevated(crp, CRP_AT_RISK) #and not is_above(crp, CRP_ELEVATED)
    crp_likely = is_above(crp, CRP_ELEVATED)
    esr_flag = note(triggered, "esr", esr, is_elevated(esr, esr_max), f">={esr_max}")
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

    if (crp_likely and esr_flag and nlr_likely) or (crp_likely and esr_flag and nlr_likely and (family_history or symptoms)):
        category = "Significant Pattern"
    elif crp_early and esr_flag and nlr_early:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Psoriasis", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Psoriasis"),
        "Triggering Parameters": ", ".join(triggered),
    }]


#  Rheumatoid Arthritis​
def evaluate_rheumatoid_arthritis(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Rheumatoid Arthritis"):
        return [{"Condition": "Rheumatoid Arthritis", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    age = patient.get("age")
    esr_max = esr_upper_normal(sex, age) #as per the prompt (children/newborn age case not defined)
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
        note(triggered, "total_c", total_c, True, f">={TOTAL_CHOLESTEROL_ELEVATED}")
    elif total_c_early:
        note(triggered, "total_c", total_c, True, f">={TOTAL_CHOLESTEROL_AT_RISK}")

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

    esr_flag = note(triggered, "esr", esr, is_elevated(esr, esr_max), f">={esr_max}")
    crp_early = is_elevated(crp, CRP_AT_RISK) #and not is_above(crp, CRP_ELEVATED)
    crp_likely = is_above(crp, CRP_ELEVATED)
    if crp_likely:
        note(triggered, "crp", crp, True, f">{CRP_ELEVATED}")
    elif crp_early:
        note(triggered, "crp", crp, True, f">={CRP_AT_RISK}")

    if (dyslipidemia_flag_likely and crp_likely and esr_flag) or (dyslipidemia_flag_likely and crp_likely and esr_flag and (family_history or symptoms)):
        category = "Significant Pattern"
    # added for early pattern w/o user context (labs only)
    elif dyslipidemia_flag_early and crp_early and esr_flag:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Rheumatoid Arthritis", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Rheumatoid Arthritis"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Rhinitis 
def evaluate_rhinitis(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Rhinitis"):
        return [{"Condition": "Rhinitis", "Category": GENE_NOT_FOUND}]

    eosinophils = labs.get("eosinophils")
    ige = labs.get("ige")

    triggered = []
    eosinophil_flag_early = is_elevated(eosinophils, EOSINOPHILS_BORDERLINE) #and not is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
    eosinophil_flag_likely = is_above(eosinophils, EOSINOPHILS_UPPER_NORMAL)
    if eosinophil_flag_likely:
        note(triggered, "eosinophils", eosinophils, True, f">{EOSINOPHILS_UPPER_NORMAL}")
    elif eosinophil_flag_early:
        note(triggered, "eosinophils", eosinophils, True, f">={EOSINOPHILS_BORDERLINE}")
    ige_flag = note(triggered, "ige", ige, is_elevated(ige, IGE_UPPER_NORMAL), f">={IGE_UPPER_NORMAL}")

    if (eosinophil_flag_likely and ige_flag) or (eosinophil_flag_likely and ige_flag and (symptoms or family_history)):
        category = "Significant Pattern"
    elif eosinophil_flag_early or ige_flag:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Rhinitis", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Rhinitis"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Type 2 Diabetes
def evaluate_t2d(labs, patient, genetics, family_history, symptoms=False):
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
    insulin_resistance_flag = note(triggered, "fasting_insulin", fasting_insulin, is_above(fasting_insulin, FASTING_INSULIN_RANGE[0]), f">={FASTING_INSULIN_RANGE[0]}")

    if (hyperglycemia) or (hyperglycemia and (symptoms or family_history)):
        category = "Significant Pattern"
    elif hyperglycemia_early and insulin_resistance_flag:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Type 2 Diabetes", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Type 2 Diabetes"),
        "Triggering Parameters": ", ".join(triggered),
    }]

# =============================================================================
# (14)
# =============================================================================

# Maturity-Onset Diabetes of the Young"
def evaluate_mody(labs, patient, genetics, family_history, symptoms=False):
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

    if (hyperglycemia and age is not None and age < 25) or (hyperglycemia and family_history):
        category = "Significant Pattern"
    elif hyperglycemia_early or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Maturity-Onset Diabetes of the Young", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Maturity-Onset Diabetes of the Young"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Multiple Endocrine Neoplasia Type 4
def evaluate_men4(labs, patient, genetics, family_history, endocrine_conditions=None):
    if not has_flagged_gene(genetics, "Multiple Endocrine Neoplasia Type 4"):
        return [{"Condition": "Multiple Endocrine Neoplasia Type 4", "Category": GENE_NOT_FOUND}]

    calcium = labs.get("calcium")
    pth = labs.get("pth")

    triggered = []
    hypercalcemia = note(triggered, "calcium", calcium, is_above(calcium, CALCIUM_UPPER_NORMAL), f">{CALCIUM_UPPER_NORMAL}")
    elevated_pth = note(triggered, "pth", pth, is_above(pth, PTH_UPPER_NORMAL), f">{PTH_UPPER_NORMAL}")
    phpt_pattern = hypercalcemia and elevated_pth

    if phpt_pattern and (endocrine_conditions or family_history):
        category = "Significant Pattern"
    elif phpt_pattern or (endocrine_conditions and family_history):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Multiple Endocrine Neoplasia Type 4", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Multiple Endocrine Neoplasia Type 4"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Muscular Dystrophy
def evaluate_muscular_dystrophy(labs, patient, genetics, family_history, symptoms=False):
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
    # the just marked case added for early pattern w/o user context
    elif marked_ck_elevation or (moderate_ck_elevation and (symptoms or family_history)):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Muscular Dystrophy", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Muscular Dystrophy"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Familial Hypercholesterolemia
def evaluate_fh(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Familial Hypercholesterolemia"):
        return [{"Condition": "Familial Hypercholesterolemia", "Category": GENE_NOT_FOUND}]

    ldl = labs.get("ldl_c")
    non_hdl = labs.get("non_hdl_c")
    lpa = labs.get("lpa")
    age = patient.get("age")

    triggered = []
    elevated_early = (
        is_above(ldl, LDL_C_AT_RISK) or 
        is_above(non_hdl, NON_HDL_C_ELEVATED) or
        is_elevated(lpa, LPA_AT_RISK)
    )
    elevated_likely = (
        is_above(ldl, LDL_C_ELEVATED) or 
        is_above(non_hdl, NON_HDL_C_SEVERE) or
        is_elevated(lpa, LPA_ELEVATED)
    )

    if elevated_likely:
        note(triggered, "ldl_c", ldl, is_above(ldl, LDL_C_ELEVATED), f">{LDL_C_ELEVATED}") \
            or note(triggered, "non_hdl_c", non_hdl, is_above(non_hdl, NON_HDL_C_SEVERE), f">{NON_HDL_C_SEVERE}") \
            or note(triggered, "lpa", lpa, is_elevated(lpa, LPA_ELEVATED), f">={LPA_ELEVATED}")
    elif elevated_early:
        note(triggered, "ldl_c", ldl, is_above(ldl, LDL_C_AT_RISK), f">{LDL_C_AT_RISK}") \
            or note(triggered, "non_hdl_c", non_hdl, is_above(non_hdl, NON_HDL_C_ELEVATED), f">{NON_HDL_C_ELEVATED}") \
            or note(triggered, "lpa", lpa, is_elevated(lpa, LPA_AT_RISK), f">={LPA_AT_RISK}")

    if elevated_likely or (elevated_early and age is not None and age < 40) or (elevated_early and family_history):
        category = "Significant Pattern"
    elif elevated_early or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Familial Hypercholesterolemia", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Familial Hypercholesterolemia"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Cardiomyopathy
def evaluate_cardiomyopathy(labs, patient, genetics, family_history, symptoms=False, imaging_performed=False, imaging_findings=None,):
    if not has_flagged_gene(genetics, "Cardiomyopathy"):
        return [{"Condition": "Cardiomyopathy", "Category": GENE_NOT_FOUND}]

    hemoglobin = labs.get("hemoglobin")
    eosinophils = labs.get("eosinophils")
    sex = patient.get("sex")

    HIGH_SPECIFICITY_FINDINGS = [] # a subset of imaging findings
    
    imaging_findings = imaging_findings or []
    real_findings = [f for f in imaging_findings if f != "Normal"]

    imaging_abnormal = imaging_performed and len(real_findings) > 0
    imaging_multiple = imaging_performed and len(real_findings) >= 2
    imaging_high_spec = imaging_performed and any(
        f in HIGH_SPECIFICITY_FINDINGS for f in real_findings
    )

    triggered = []
    hb_normal = sex_based_threshold(sex, 13, 11.5)
    hb_severely_low = hb_normal is not None and note(triggered, "hemoglobin", hemoglobin, is_below(hemoglobin, hb_normal - 2), f"<{hb_normal - 2}") #because python evaluates arithmetic part before helper function's "not None" check.
    eos_high = note(triggered, "eosinophils", eosinophils, is_above(eosinophils, EOSINOPHILS_ELEVATED), f">={EOSINOPHILS_ELEVATED}")
    labs_severe = hb_severely_low or eos_high

    gene_trigger = get_gene_trigger(genetics, "Cardiomyopathy")

    if hemoglobin is None and eosinophils is None and not imaging_performed:
        return [{
            "Condition": "Cardiomyopathy",
            "Category": "Elevated Susceptibility",
            "Triggering Genetics": gene_trigger,
            "Triggering Parameters": ", ".join(triggered),
        }]

    # --- Case 1: imaging performed and abnormal -> imaging drives the call ---
    if imaging_performed and imaging_abnormal:
        if imaging_high_spec or imaging_multiple or symptoms or family_history:
            category = "Significant Pattern"
        else:
            category = "Early Pattern"

    # --- Case 2: imaging performed and normal ---
    elif imaging_performed and not imaging_abnormal:
        risk_factor_count = sum([family_history, symptoms, labs_severe])
        if risk_factor_count >= 2:
            category = "Early Pattern" #Re-imaging recommended (Echocardiogram/Cardiac MRI)
        else:
            category = "Elevated Susceptibility"

    # --- Case 3: no imaging performed -> route to imaging ---
    else:
        if labs_severe or family_history or symptoms:
            category = "Early Pattern" #Advanced cardiac imaging recommended (Echocardiogram/Cardiac MRI)
        else:
            category = "Elevated Susceptibility"

    return [{
        "Condition": "Cardiomyopathy", "Category": category,
        "Triggering Genetics": gene_trigger,
        "Triggering Parameters": ", ".join(triggered),
    }]


# Coronary Artery Disease
def evaluate_cad(labs, patient, genetics, family_history, symptoms=False):
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
        is_elevated(apob, APOB_AT_RISK) or
        is_elevated(lpa, LPA_AT_RISK)
    )
    ldl_high = (
        is_elevated(ldl, LDL_C_ELEVATED) or 
        is_elevated(non_hdl, NON_HDL_C_SEVERE) or
        is_elevated(apob, APOB_ELEVATED) or
        is_elevated(lpa, LPA_ELEVATED)
    )
    if ldl_high:
        note(triggered, "ldl_c", ldl, is_elevated(ldl, LDL_C_ELEVATED), f">={LDL_C_ELEVATED}") \
            or note(triggered, "non_hdl_c", non_hdl, is_elevated(non_hdl, NON_HDL_C_SEVERE), f">={NON_HDL_C_SEVERE}") \
            or note(triggered, "apob", apob, is_elevated(apob, APOB_ELEVATED), f">={APOB_ELEVATED}") \
            or note(triggered, "lpa", lpa, is_elevated(lpa, LPA_ELEVATED), f">={LPA_ELEVATED}")
    elif ldl_moderate:
        note(triggered, "ldl_c", ldl, is_elevated(ldl, LDL_C_AT_RISK), f">={LDL_C_AT_RISK}") \
            or note(triggered, "non_hdl_c", non_hdl, is_elevated(non_hdl, NON_HDL_C_ELEVATED), f">={NON_HDL_C_ELEVATED}") \
            or note(triggered, "apob", apob, is_elevated(apob, APOB_AT_RISK), f">={APOB_AT_RISK}") \
            or note(triggered, "lpa", lpa, is_elevated(lpa, LPA_AT_RISK), f">={LPA_AT_RISK}")

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
        (any_moderate and symptoms and family_history)
    ):
        category = "Significant Pattern"
    # (any_moderate + age) added for early pattern w/o user context (labs only)
    elif family_history or (any_moderate and symptoms) or (any_moderate and age is not None and age < 40):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Coronary Artery Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Coronary Artery Disease"),
        "Triggering Parameters": ", ".join(triggered),
    }]

# Familial Hypertriglyceridemia
def evaluate_hypertriglyceridemia(labs, patient, genetics, family_history, symptoms=False):
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

    vldl_elevated = note(triggered, "vldl", vldl, is_elevated(vldl, VLDL_UPPER_NORMAL), f">={VLDL_UPPER_NORMAL}")
    any_moderate_abnormal = tg_moderate or hdl_low or vldl_elevated

    if (
        tg_very_severe or
        (tg_severe and age is not None and age < 40) or
        (tg_severe and family_history) or
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
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Familial Hypertriglyceridemia", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Familial Hypertriglyceridemia"),
        "Triggering Parameters": ", ".join(triggered),
    }]

# HDL Deficiency
def evaluate_hdl_deficiency(labs, patient, genetics, family_history, symptoms=False):
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
        (isolated_moderate and ((age is not None and age < 40) or family_history))
    ):
        category = "Significant Pattern"
    elif isolated_moderate or mixed_moderate or family_history:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "HDL Deficiency", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "HDL Deficiency"),
        "Triggering Parameters": ", ".join(triggered),
    }]


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
    if not has_flagged_gene(genetics, "Li-Fraumeni Syndrome"):
        return [{"Condition": "Li-Fraumeni Syndrome", "Category": GENE_NOT_FOUND}]

    age = patient.get("age")
    conditions = conditions or []

    pathognomonic_tumors = [] #subset of conditions
    lfs_spectrum_tumors = [] #subset of conditions

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
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Li-Fraumeni Syndrome", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Li-Fraumeni Syndrome"),
        "Triggering Parameters": "",
    }]


# PTEN Hamartoma Tumor Syndrome
def evaluate_phts(labs, patient, genetics, family_history, symptoms=False, conditions=None):
    if not has_flagged_gene(genetics, "PTEN Hamartoma Tumor Syndrome"):
        return [{"Condition": "PTEN Hamartoma Tumor Syndrome", "Category": GENE_NOT_FOUND}]

    conditions = conditions or []
    major_criteria = []   # subset of conditions
    minor_criteria = []   # subset of conditions
    macrocephaly_plus_criteria = []  # subset of conditions

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
        category = "Elevated Susceptibility"

    return [{
        "Condition": "PTEN Hamartoma Tumor Syndrome", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "PTEN Hamartoma Tumor Syndrome"),
        "Triggering Parameters": "",
    }]


# Gilbert Syndrome
def evaluate_gilbert_syndrome(labs, patient, genetics, family_history, past_history, symptoms=False):
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
    parameter = (
        total_bilirubin is not None and TOTAL_BILIRUBIN_RANGE[0] < total_bilirubin < TOTAL_BILIRUBIN_RANGE[1] and
        direct_bilirubin is not None and direct_bilirubin > DIRECT_BILIRUBIN_LOWER and
        indirect_bilirubin is not None and indirect_bilirubin < INDIRECT_BILIRUBIN_UPPER and
        alt is not None and alt <= ALT_RANGE[1] and
        ast is not None and ast <= AST_UPPER_STRICT and
        alp is not None and ALP_RANGE[0] < alp <= ALP_RANGE[1] and
        ggt is not None and GGT_LFT_RANGE[0] < ggt <= GGT_LFT_RANGE[1]
    )
    # all 7 sub-checks are AND'd into a single pattern, so only record them once the
    # whole pattern is confirmed true — otherwise a false prefix would leave misleading
    # partial entries in `triggered` for a pattern that didn't actually fire
    if parameter:
        triggered.append(f"total_bilirubin={total_bilirubin} ({TOTAL_BILIRUBIN_RANGE[0]}-{TOTAL_BILIRUBIN_RANGE[1]})")
        triggered.append(f"direct_bilirubin={direct_bilirubin} (>{DIRECT_BILIRUBIN_LOWER})")
        triggered.append(f"indirect_bilirubin={indirect_bilirubin} (<{INDIRECT_BILIRUBIN_UPPER})")
        triggered.append(f"alt={alt} (<={ALT_RANGE[1]})")
        triggered.append(f"ast={ast} (<={AST_UPPER_STRICT})")
        triggered.append(f"alp={alp} ({ALP_RANGE[0]}-{ALP_RANGE[1]}]")
        triggered.append(f"ggt={ggt} ({GGT_LFT_RANGE[0]}-{GGT_LFT_RANGE[1]}]")

    if (
        (parameter and family_history and symptoms and past_history) or
        (parameter and family_history and past_history) or
        (parameter and family_history) or
        (family_history and symptoms) or
        (family_history and past_history)
    ):
        category = "Significant Pattern"
    elif (
        family_history or
        (parameter and symptoms and past_history) or
        (parameter and past_history) or
        parameter #added for early pattern w/o user context
    ):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Gilbert Syndrome", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Gilbert Syndrome"),
        "Triggering Parameters": ", ".join(triggered),
    }]

# Hereditary Hemochromatosis
def evaluate_hemochromatosis(labs, patient, genetics, family_history, past_history, symptoms=False):
    if not has_flagged_gene(genetics, "Hereditary Hemochromatosis"):
        return [{"Condition": "Hereditary Hemochromatosis", "Category": GENE_NOT_FOUND}]

    iron_profile_done = labs.get("iron_profile_done") #####
    ferritin = labs.get("ferritin")
    transferrin_saturation = labs.get("transferrin_saturation")
    alt = labs.get("alt")
    ast = labs.get("ast")
    hemoglobin = labs.get("hemoglobin")
    mcv = labs.get("mcv")
    mch = labs.get("mch")

    triggered = []
    abnormal_labs = (
        is_above(alt, ALT_RANGE[1]) and
        is_above(ast, AST_UPPER_STRICT) and
        is_above(hemoglobin, HEMOGLOBIN_UPPER_HEMOCHROMATOSIS) and
        is_above(mcv, MCV_UPPER_HEMOCHROMATOSIS) and
        is_above(mch, MCH_UPPER_HEMOCHROMATOSIS)
    )
    # all 5 sub-checks are AND'd into a single pattern, so only record them once the
    # whole pattern is confirmed true — otherwise a false prefix would leave misleading
    # partial entries in `triggered` for a pattern that didn't actually fire
    if abnormal_labs:
        triggered.append(f"alt={alt} (>{ALT_RANGE[1]})")
        triggered.append(f"ast={ast} (>{AST_UPPER_STRICT})")
        triggered.append(f"hemoglobin={hemoglobin} (>{HEMOGLOBIN_UPPER_HEMOCHROMATOSIS})")
        triggered.append(f"mcv={mcv} (>{MCV_UPPER_HEMOCHROMATOSIS})")
        triggered.append(f"mch={mch} (>{MCH_UPPER_HEMOCHROMATOSIS})")

    gene_trigger = get_gene_trigger(genetics, "Hereditary Hemochromatosis")

    if iron_profile_done is True:
        tsat_low = note(triggered, "transferrin_saturation", transferrin_saturation, is_below(transferrin_saturation, TRANSFERRIN_SATURATION_ELEVATED), f"<{TRANSFERRIN_SATURATION_ELEVATED}")
        ferritin_not_elevated = note(triggered, "ferritin", ferritin, ferritin is not None and ferritin <= FERRITIN_ELEVATED, f"<={FERRITIN_ELEVATED}")

        if tsat_low and ferritin_not_elevated:
            return [{
                "Condition": "Hereditary Hemochromatosis", "Category": "Elevated Susceptibility",
                "Triggering Genetics": gene_trigger,
                "Triggering Parameters": ", ".join(triggered),
            }]

        tsat_elevated = note(triggered, "transferrin_saturation", transferrin_saturation, is_elevated(transferrin_saturation, TRANSFERRIN_SATURATION_ELEVATED), f">={TRANSFERRIN_SATURATION_ELEVATED}")
        ferritin_elevated = note(triggered, "ferritin", ferritin, is_above(ferritin, FERRITIN_ELEVATED), f">{FERRITIN_ELEVATED}")

        if tsat_elevated and ferritin_elevated:
            score = sum([family_history, symptoms, past_history, abnormal_labs])
            if score >= 3:
                category = "Significant Pattern"
            elif score >= 1:
                category = "Early Pattern" #the score=1 covers early pattern w/o user context (abnormal labs only)
            else:
                category = "Elevated Susceptibility"
            return [{
                "Condition": "Hereditary Hemochromatosis", "Category": category,
                "Triggering Genetics": gene_trigger,
                "Triggering Parameters": ", ".join(triggered),
            }]

        return [{
            "Condition": "Hereditary Hemochromatosis", "Category": "Elevated Susceptibility",
            "Triggering Genetics": gene_trigger,
            "Triggering Parameters": ", ".join(triggered),
        }]

    if iron_profile_done is False:
        score = sum([family_history, symptoms, past_history, abnormal_labs])
        if score >= 3:
            category = "Significant Pattern"
        elif score >= 1:
            category = "Early Pattern" #the score=1 covers early pattern w/o user context (abnormal labs only)
        else:
            category = "Elevated Susceptibility"
        return [{
            "Condition": "Hereditary Hemochromatosis", "Category": category,
            "Triggering Genetics": gene_trigger,
            "Triggering Parameters": ", ".join(triggered),
        }]
    return [{
        "Condition": "Hereditary Hemochromatosis",
        "Category": "Elevated Susceptibility",
        "Triggering Genetics": gene_trigger,
        "Triggering Parameters": ", ".join(triggered),
    }]


# Celiac Disease
def evaluate_celiac_disease(labs, patient, genetics, family_history, past_history, symptoms=False):
    if not has_flagged_gene(genetics, "Celiac Disease"):
        return [{"Condition": "Celiac Disease", "Category": GENE_NOT_FOUND}]

    ttg = labs.get("ttg-iga") #or labs.get("Tissue Transglutaminase (tTG) Antibody, IgA")

    TTG_Early = is_elevated(ttg, TTG_RANGE[0]) #and not is_above(ttg, TTG_RANGE[1])
    TTG_Likely = is_above(ttg, TTG_RANGE[1])

    triggered = []
    if TTG_Likely:
            note(triggered, "tTG-IgA", ttg, True, f">{TTG_RANGE[1]}")
    elif TTG_Early:
        note(triggered, "tTG-IgA", ttg, True, f">={TTG_RANGE[0]}")
    

    if (
        TTG_Likely or
        (TTG_Early and (family_history or symptoms or past_history))
    ):
        category = "Significant Pattern"
    elif (
        TTG_Early
    ):
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "Celiac Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Celiac Disease"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# calling:

data = load_patient_data(args.monogenic_json, args.prs_json, args.apoe_json, args.blood_csv, sex=args.sex, age=args.age)

results = [
    # --- 14: labs/genetics-driven conditions ---
    evaluate_alzheimers(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_asthma(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_atopic_dermatitis(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_copd(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_hyperthyroidism(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_hypothyroidism(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_ibd(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_nafld(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_osteoarthritis(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_parkinsons(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_psoriasis(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_rheumatoid_arthritis(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_rhinitis(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_t2d(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),

    # --- 14: second batch ---
    evaluate_mody(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    # make family_history=[] as default instead of False for when including user context json.
    evaluate_men4(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, endocrine_conditions=None,
    ),
    evaluate_muscular_dystrophy(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, symptoms=False,
    ),
    evaluate_gilbert_syndrome(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, past_history=False, symptoms=False,
    ),
    evaluate_hemochromatosis(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, past_history=False, symptoms=False,
    ),
    evaluate_celiac_disease(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, past_history=False, symptoms=False,
    ),
    evaluate_fh(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_cardiomyopathy(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, symptoms=False, imaging_performed=False, imaging_findings=None,
    ),
    evaluate_cad(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, symptoms=False,
    ),
    evaluate_hypertriglyceridemia(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_hdl_deficiency(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    # evaluate_apob(
    #     labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
    #     family_history=False,
    # ),
    evaluate_li_fraumeni_syndrome(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, conditions=None,
    ),
    evaluate_phts(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, conditions=None,
    ),
]

CONDITION_DOMAINS = {
    "Alzheimer's Disease": "Brain Health",
    "Asthma": "Respiratory Health",
    "Atopic Dermatitis/Eczema": "Skin Health",
    "COPD": "Respiratory Health",
    "Hyperthyroidism": "Endocrine Health",
    "Hypothyroidism": "Endocrine Health",
    "Inflammatory Bowel Disease": "Digestive Health",
    "NAFLD": "Digestive Health",
    "Osteoarthritis": "Joints Health",
    "Parkinson's Disease": "Brain Health",
    "Psoriasis": "Skin Health",
    "Rheumatoid Arthritis": "Immune Health",
    "Rhinitis": "Respiratory Health",
    "Type 2 Diabetes": "Endocrine Health",
    "Maturity-Onset Diabetes of the Young": "Endocrine Health",
    "Multiple Endocrine Neoplasia Type 4": "Endocrine Health",
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
}

for r in results:
    entry = r[0]
    entry["Domain"] = CONDITION_DOMAINS.get(entry["Condition"], "")
    entry["Triggering PRS"] = get_triggering_prs(entry["Condition"], data["genetics"])
    entry["Snapshot Category"] = get_snapshot_category(entry["Condition"], entry["Category"], data["genetics"])
    # entry["Category"] = resolve_typical_label(entry["Condition"], entry["Category"], data["genetics"])

active_findings = [
    r[0]
    for r in results
]

with open(args.output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Domain", "Condition", "Category", "Triggering Genetics", "Triggering Parameters", "Triggering PRS", "Snapshot Category"])
    writer.writeheader()
    writer.writerows(active_findings)

print(f"Results written to {args.output_csv}")
