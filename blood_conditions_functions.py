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
    with open(monogenic_json_path) as f:
        monogenic_data = json.load(f)
    for entry in monogenic_data:
        sig = (entry.get("ClinVar Significance") or "").strip().lower()
        if sig in PATHOGENIC_VALUES:
            gene = (entry.get("Gene Symbol") or "").strip()
            if gene:
                flagged_genes.append(gene)

    # --- 2. Polygenic PRS elevated conditions ---
    prs_elevated_conditions = set()
    with open(prs_json_path) as f:
        prs_data = json.load(f)
    for condition_key, info in prs_data.items():
        category = (info.get("category") or "").strip().lower()
        if category in ELEVATED_CATEGORIES:
            prs_elevated_conditions.add(condition_key)

    # --- 3. APOE status ---
    apoe_status = ""
    with open(apoe_json_path) as f:
        apoe_data = json.load(f)
    apoe_status = (apoe_data.get("APOE_Status") or "").strip()
    if "4" in apoe_status:  # e.g. "e3/e4" or "e4/e4" -> ε4 allele present
        flagged_genes.append("APOE")

    return {
        "flagged_genes": flagged_genes,
        "prs_elevated_conditions": prs_elevated_conditions,
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
TG_HDL_RATIO_ELEVATED = 2.5
APOB_MALE_RANGE = (60, 130)           # mg/dL
APOB_FEMALE_RANGE = (55, 125)         # mg/dL

# --- Glucose / Metabolic ---
FASTING_GLUCOSE_AT_RISK = 100         # mg/dL
FASTING_GLUCOSE_ELEVATED = 126        # mg/dL
HBA1C_AT_RISK = 5.7                   # %
HBA1C_ELEVATED = 6.5                  # %
FASTING_INSULIN_RANGE = (2, 25)       # uIU/mL

# --- Inflammatory Markers ---
CRP_AT_RISK = 1.0                     # mg/dL
CRP_ELEVATED = 3.0                    # mg/dL
IGE_UPPER_NORMAL = 100                # IU/mL
FIBRINOGEN_RANGE = (200, 400)         # mg/dL
FECAL_CALPROTECTIN_AT_RISK = 50       # ug/g
FECAL_CALPROTECTIN_ELEVATED = 200     # ug/g

def esr_upper_normal(sex, age):
    """Returns ESR upper limit of normal (mm/hr) by sex and age."""
    if sex not in SEX_VALUES or age is None:
        return None
    
    if sex == "male":
        return 15 if age < 50 else 20
    else:
        return 20 if age < 50 else 30

# --- CBC ---
EOSINOPHILS_UPPER_NORMAL = 0.4        # x10^9/L
NEUTROPHILS_RANGE = (1.8, 7.5)        # x10^9/L
PLATELETS_RANGE = (150, 400)          # x10^9/L
NLR_RANGE = (0.8, 3.5)
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
FREE_T4_RANGE = (0.8, 1.8)            # ng/dL
FREE_T3_RANGE = (2.3, 4.4)            # pg/mL
TPOAB_UPPER_NORMAL = 35                # IU/mL

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
VITAMIN_D_DEFICIENT = 20                       # ng/mL (Celiac)
MCV_LOWER_CELIAC = 82                          # fL
MCH_LOWER_CELIAC = 27                          # pg
MCHC_LOWER_CELIAC = 31.5                       # g/dL

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
    "COPD": NO_GENE_GATE,  # no gene listed
    "Hyperthyroidism": ["CTLA4", "PTPN22"],
    "Hypothyroidism": ["CTLA4"],
    "Inflammatory Bowel Disease": ["NOD2", "ATG16L1", "IL23R", "MYC/POU5F1B", "IL27"],
    "NAFLD": ["PNPLA3", "TM6SF2"],
    "Osteoarthritis": ["COL11A1", "SMAD3"],
    "Parkinson's Disease": ["LRRK2", "GBA"],
    "Psoriasis": NO_GENE_GATE,  # no gene listed
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

# =============================================================================
# (14)
# =============================================================================

# Family history and symptoms of respective conditions, are mentioned in the deck for ref

GENE_NOT_FOUND = "Related mutation in gene(s) not found"

# Alzheimer's Disease
def evaluate_alzheimers(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Alzheimer's Disease"):
        return [{"Condition": "Alzheimer's Disease", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

    ldl_c = labs.get("ldl_c")
    non_hdl_c = labs.get("non_hdl_c")
    triglycerides = labs.get("triglycerides")
    hdl_c = labs.get("hdl_c")
    fasting_glucose = labs.get("fasting_glucose")
    hba1c = labs.get("hba1c")
    crp = labs.get("crp")

    triggered = []
    dyslipidemia_flag = (
        note(triggered, "ldl_c", ldl_c, is_elevated(ldl_c, LDL_C_AT_RISK), f">={LDL_C_AT_RISK}")
        or note(triggered, "non_hdl_c", non_hdl_c, is_elevated(non_hdl_c, NON_HDL_C_NORMAL), f">={NON_HDL_C_NORMAL}")
        or note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")
        or note(triggered, "hdl_c", hdl_c, is_below(hdl_c, hdl_min), f"<{hdl_min}")
    )
    glycemic_flag = (
        note(triggered, "fasting_glucose", fasting_glucose, is_elevated(fasting_glucose, FASTING_GLUCOSE_AT_RISK), f">={FASTING_GLUCOSE_AT_RISK}")
        or note(triggered, "hba1c", hba1c, is_elevated(hba1c, HBA1C_AT_RISK), f">={HBA1C_AT_RISK}")
    )
    crp_flag = note(triggered, "crp", crp, is_elevated(crp, CRP_AT_RISK), f">={CRP_AT_RISK}")
    prediabetic_flag = (
        note(triggered, "fasting_glucose", fasting_glucose, fasting_glucose is not None and 100 <= fasting_glucose <= 125, "[100-125]")
        or note(triggered, "hba1c", hba1c, hba1c is not None and 5.7 <= hba1c <= 6.4, "[5.7-6.4]")
    )

    if dyslipidemia_flag and glycemic_flag and crp_flag and family_history:
        category = "Likely Disease Onset"
    elif ((dyslipidemia_flag or crp_flag) and prediabetic_flag and family_history) or (dyslipidemia_flag and glycemic_flag and crp_flag):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    # neutrophils = labs.get("neutrophils") #just mentioned in parameters, but not in triggers for cases
    triggered = []
    eosinophil_flag = note(triggered, "eosinophils", eosinophils, is_elevated(eosinophils, EOSINOPHILS_UPPER_NORMAL), f">={EOSINOPHILS_UPPER_NORMAL}")
    # neutrophil_flag = is_elevated(neutrophils, NEUTROPHILS_RANGE[1])

    if eosinophil_flag and symptoms and family_history:
        category = "Likely Disease Onset"
    # just eosinophil_flag updated for early pattern w/o user context (labs only)
    elif eosinophil_flag:
        category = "Early Pattern"
    else:
        category = "Typical"

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
    eosinophil_flag = note(triggered, "eosinophils", eosinophils, is_elevated(eosinophils, EOSINOPHILS_UPPER_NORMAL), f">={EOSINOPHILS_UPPER_NORMAL}")

    if eosinophil_flag and symptoms and family_history:
        category = "Likely Disease Onset"
    # just eosinophil_flag for early pattern w/o user context (labs only)
    elif eosinophil_flag:
        category = "Early Pattern"
    else:
        category = "Typical"

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
    crp_elevated_flag = note(triggered, "crp", crp, is_elevated(crp, CRP_ELEVATED), f">={CRP_ELEVATED}")
    crp_at_risk_flag = note(triggered, "crp", crp, is_elevated(crp, CRP_AT_RISK), f">={CRP_AT_RISK}")
    eosinophil_flag = note(triggered, "eosinophils", eosinophils, is_elevated(eosinophils, EOSINOPHILS_UPPER_NORMAL), f">={EOSINOPHILS_UPPER_NORMAL}")
    fibrinogen_flag = note(triggered, "fibrinogen", fibrinogen, is_elevated(fibrinogen, FIBRINOGEN_RANGE[1]), f">={FIBRINOGEN_RANGE[1]}")
    neutrophil_flag = note(triggered, "neutrophils", neutrophils, is_elevated(neutrophils, NEUTROPHILS_RANGE[1]), f">={NEUTROPHILS_RANGE[1]}")

    if (crp_elevated_flag or fibrinogen_flag) and (eosinophil_flag or neutrophil_flag) and symptoms and family_history:
        category = "Likely Disease Onset"
    # just eosinophil_flag or crp_at_risk_flag added for early pattern w/o user context (labs only)
    elif (crp_at_risk_flag or eosinophil_flag):
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "COPD", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "COPD"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Eczema (merged with Atopic dermatitis (as done in prompt))
# def evaluate_eczema(labs, patient, genetics, family_history, symptoms=False):
#     if not has_flagged_gene(genetics, "Eczema"):
#         return [{"Condition": "Eczema", "Category": GENE_NOT_FOUND}]

#     eosinophils = labs.get("eosinophils")
#     eosinophil_flag = is_elevated(eosinophils, EOSINOPHILS_UPPER_NORMAL)

#     if eosinophil_flag and symptoms and family_history:
#         category = "Likely Disease Onset"
#     elif eosinophil_flag and (symptoms or family_history):
#         category = "Early Pattern"
#     else:
#         category = "Typical"

#     return [{"Condition": "Eczema", "Category": category}]


# Hyperthyroidism
def evaluate_hyperthyroidism(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Hyperthyroidism"):
        return [{"Condition": "Hyperthyroidism", "Category": GENE_NOT_FOUND}]

    tsh = labs.get("tsh")
    free_t3 = labs.get("free_t3")
    free_t4 = labs.get("free_t4")
    tpoab = labs.get("tpoab")

    triggered = []
    tsh_suppressed = note(triggered, "tsh", tsh, is_below(tsh, TSH_RANGE[0]), f"<{TSH_RANGE[0]}")
    hormone_flag = (
        note(triggered, "free_t3", free_t3, free_t3 is not None and free_t3 > FREE_T3_RANGE[1], f">{FREE_T3_RANGE[1]}")
        or note(triggered, "free_t4", free_t4, free_t4 is not None and free_t4 > FREE_T4_RANGE[1], f">{FREE_T4_RANGE[1]}")
    )
    tpoab_flag = note(triggered, "tpoab", tpoab, is_elevated(tpoab, TPOAB_UPPER_NORMAL), f">={TPOAB_UPPER_NORMAL}")

    if (tsh_suppressed and hormone_flag and tpoab_flag):
        category = "Likely Disease Onset"
    elif tsh_suppressed and (hormone_flag or tpoab_flag):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    tsh_high = note(triggered, "tsh", tsh, is_elevated(tsh, 10.0), f">={10.0}")
    tsh_mild = note(triggered, "tsh", tsh, tsh is not None and 4.1 <= tsh <= 10.0, "[4.1-10.0]")
    free_t4_low = note(triggered, "free_t4", free_t4, is_below(free_t4, FREE_T4_RANGE[0]), f"<{FREE_T4_RANGE[0]}")
    tpoab_flag = note(triggered, "tpoab", tpoab, is_elevated(tpoab, TPOAB_UPPER_NORMAL), f">={TPOAB_UPPER_NORMAL}")

    if (tsh_high and free_t4_low and tpoab_flag):
        category = "Likely Disease Onset"
    # tsh_mild and tpoab_flag case covers early pattern w/o user context
    elif tsh_mild and (tpoab_flag or family_history):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    crp_flag = note(triggered, "crp", crp, is_elevated(crp, CRP_ELEVATED), f">={CRP_ELEVATED}")
    anemia_flag = note(triggered, "hemoglobin", hemoglobin, is_below(hemoglobin, hb_min), f"<{hb_min}")
    albumin_flag = note(triggered, "albumin", albumin, is_below(albumin, ALBUMIN_RANGE[0]), f"<{ALBUMIN_RANGE[0]}")
    # (crp_flag and anemia_flag and albumin_flag) added for early pattern w/o user context (labs only)
    if ((crp_flag or anemia_flag or albumin_flag) and symptoms and family_history) or (crp_flag and anemia_flag and albumin_flag):
        category = "Likely Disease Onset"
    # just crp_flag for early pattern w/o user context (labs only)
    elif crp_flag:
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "Inflammatory Bowel Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Inflammatory Bowel Disease"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Non alcoholic fatty liver disease
def evaluate_nafld(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "NAFLD"):
        return [{"Condition": "Non Alcoholic Fatty Liver Disease", "Category": GENE_NOT_FOUND}]

    ast = labs.get("ast")
    alt = labs.get("alt")
    ggt = labs.get("ggt")
    fib4 = labs.get("fib4")
    platelets = labs.get("platelets")

    triggered = []
    lft_flag = (
        note(triggered, "ast", ast, is_elevated(ast, AST_RANGE[1]), f">={AST_RANGE[1]}")
        or note(triggered, "alt", alt, is_elevated(alt, ALT_RANGE[1]), f">={ALT_RANGE[1]}")
        or note(triggered, "ggt", ggt, is_elevated(ggt, GGT_RANGE[1]), f">={GGT_RANGE[1]}")
    )
    fib4_elevated = note(triggered, "fib4", fib4, is_elevated(fib4, FIB4_ELEVATED), f">={FIB4_ELEVATED}")
    fib4_at_risk = note(triggered, "fib4", fib4, fib4 is not None and FIB4_AT_RISK <= fib4 < FIB4_ELEVATED, f"[{FIB4_AT_RISK}-{FIB4_ELEVATED})")
    platelets_low = note(triggered, "platelets", platelets, is_below(platelets, PLATELETS_RANGE[0]), f"<{PLATELETS_RANGE[0]}")
    # (lft_flag and fib4_elevated and platelets_low) added for early pattern w/o user context (labs only)
    if (lft_flag and fib4_elevated and platelets_low):
        category = "Likely Disease Onset"
    # lft_flag or fib4_at_risk for early pattern w/o user context (labs only)
    elif (lft_flag or fib4_at_risk):
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "Non Alcoholic Fatty Liver Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "NAFLD"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Osteoarthritis
def evaluate_osteoarthritis(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Osteoarthritis"):
        return [{"Condition": "Osteoarthritis", "Category": GENE_NOT_FOUND}]

    crp = labs.get("crp")
    triggered = []
    crp_flag = note(triggered, "crp", crp, is_elevated(crp, CRP_AT_RISK), f">={CRP_AT_RISK}")

    if crp_flag and symptoms and family_history:
        category = "Likely Disease Onset"
    # just crp_flag for early pattern w/o user context (labs only)
    elif crp_flag:
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "Osteoarthritis", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Osteoarthritis"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Parkinson's Disease
def evaluate_parkinsons(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Parkinson's Disease"):
        return [{"Condition": "Parkinson's Disease", "Category": GENE_NOT_FOUND}]

    sex = patient.get("sex")
    hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

    ldl_c = labs.get("ldl_c")
    triglycerides = labs.get("triglycerides")
    hdl_c = labs.get("hdl_c")
    fasting_glucose = labs.get("fasting_glucose")
    hba1c = labs.get("hba1c")
    crp = labs.get("crp")

    triggered = []
    metabolic_flag = (
        note(triggered, "ldl_c", ldl_c, is_elevated(ldl_c, LDL_C_AT_RISK), f">={LDL_C_AT_RISK}")
        or note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")
        or note(triggered, "hdl_c", hdl_c, is_below(hdl_c, hdl_min), f"<{hdl_min}")
        or note(triggered, "fasting_glucose", fasting_glucose, is_elevated(fasting_glucose, FASTING_GLUCOSE_AT_RISK), f">={FASTING_GLUCOSE_AT_RISK}")
        or note(triggered, "hba1c", hba1c, is_elevated(hba1c, HBA1C_AT_RISK), f">={HBA1C_AT_RISK}")
    )
    crp_flag = note(triggered, "crp", crp, is_elevated(crp, CRP_AT_RISK), f">={CRP_AT_RISK}")

    if metabolic_flag and crp_flag and symptoms and family_history:
        category = "Likely Disease Onset"
    # added for early pattern w/o user context (labs only)
    elif (metabolic_flag or crp_flag):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    crp_flag = note(triggered, "crp", crp, is_elevated(crp, CRP_ELEVATED), f">={CRP_ELEVATED}")
    esr_flag = note(triggered, "esr", esr, is_elevated(esr, esr_max), f">={esr_max}")
    nlr_flag = note(triggered, "nlr", nlr, is_elevated(nlr, NLR_RANGE[1]), f">={NLR_RANGE[1]}")

    if crp_flag and esr_flag and nlr_flag and symptoms:
        category = "Likely Disease Onset"
    # added for early pattern w/o user context (labs only)
    elif (crp_flag or esr_flag or nlr_flag):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    hdl_c = labs.get("hdl_c")
    triglycerides = labs.get("triglycerides")

    triggered = []
    inflammatory_flag = (
        note(triggered, "crp", crp, is_elevated(crp, CRP_ELEVATED), f">={CRP_ELEVATED}")
        or note(triggered, "esr", esr, is_elevated(esr, esr_max), f">={esr_max}")
    )
    lipid_flag = (
        note(triggered, "hdl_c", hdl_c, is_below(hdl_c, hdl_min), f"<{hdl_min}")
        or note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")
    )

    if inflammatory_flag and lipid_flag and symptoms and family_history:
        category = "Likely Disease Onset"
    # added for early pattern w/o user context (labs only)
    elif (inflammatory_flag and lipid_flag):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    eosinophil_flag = note(triggered, "eosinophils", eosinophils, is_elevated(eosinophils, EOSINOPHILS_UPPER_NORMAL), f">={EOSINOPHILS_UPPER_NORMAL}")
    ige_flag = note(triggered, "ige", ige, is_elevated(ige, IGE_UPPER_NORMAL), f">={IGE_UPPER_NORMAL}")

    if eosinophil_flag and ige_flag and symptoms and family_history:
        category = "Likely Disease Onset"
    # added for early pattern w/o user context (labs only)
    elif (eosinophil_flag or ige_flag):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    elevated_flag = (
        note(triggered, "fasting_glucose", fasting_glucose, is_elevated(fasting_glucose, FASTING_GLUCOSE_ELEVATED), f">={FASTING_GLUCOSE_ELEVATED}")
        or note(triggered, "hba1c", hba1c, is_elevated(hba1c, HBA1C_ELEVATED), f">={HBA1C_ELEVATED}")
    )
    prediabetic_flag = (
        note(triggered, "fasting_glucose", fasting_glucose, fasting_glucose is not None and 100 <= fasting_glucose <= 125, "[100-125]")
        or note(triggered, "hba1c", hba1c, hba1c is not None and 5.7 <= hba1c <= 6.4, "[5.7-6.4]")
    )
    insulin_resistance_flag = note(triggered, "fasting_insulin", fasting_insulin, is_elevated(fasting_insulin, FASTING_INSULIN_RANGE[1]), f">={FASTING_INSULIN_RANGE[1]}")

    if elevated_flag and symptoms:
        category = "Likely Disease Onset"
    elif prediabetic_flag and insulin_resistance_flag:
        category = "Early Pattern"
    else:
        category = "Typical"

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

    triggered = []
    hyperglycemia = (
        note(triggered, "hba1c", hba1c, is_elevated(hba1c, HBA1C_ELEVATED), f">={HBA1C_ELEVATED}") or
        note(triggered, "fasting_glucose", fasting_glucose, is_elevated(fasting_glucose, FASTING_GLUCOSE_ELEVATED), f">={FASTING_GLUCOSE_ELEVATED}")
    )

    if (hyperglycemia and age is not None and age < 25) or (hyperglycemia and family_history):
        category = "Likely Disease Onset"
    elif hyperglycemia or family_history:
        category = "Early Pattern"
    else:
        category = "Typical"

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
        category = "Likely Disease Onset"
    elif phpt_pattern or (endocrine_conditions and family_history):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    marked_ck_elevation = note(triggered, "ck", ck, is_elevated(ck, ck_significant_multiple), f">={ck_significant_multiple}")
    moderate_ck_elevation = note(triggered, "ck", ck, is_elevated(ck, ck_moderate_multiple), f">={ck_moderate_multiple}")

    if marked_ck_elevation and (symptoms or family_history):
        category = "Likely Disease Onset"
    # the just marked case added for early pattern w/o user context
    elif marked_ck_elevation or (moderate_ck_elevation and (symptoms or family_history)):
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "Muscular Dystrophy", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Muscular Dystrophy"),
        "Triggering Parameters": ", ".join(triggered),
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
        category = "Likely Disease Onset"
    elif (
        family_history or
        (parameter and symptoms and past_history) or
        (parameter and past_history) or
        parameter #added for early pattern w/o user context
    ):
        category = "Early Pattern"
    else:
        category = "Typical"

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
                "Condition": "Hereditary Hemochromatosis", "Category": "Unlikely",
                "Triggering Genetics": gene_trigger,
                "Triggering Parameters": ", ".join(triggered),
            }]

        tsat_elevated = note(triggered, "transferrin_saturation", transferrin_saturation, is_elevated(transferrin_saturation, TRANSFERRIN_SATURATION_ELEVATED), f">={TRANSFERRIN_SATURATION_ELEVATED}")
        ferritin_elevated = note(triggered, "ferritin", ferritin, is_above(ferritin, FERRITIN_ELEVATED), f">{FERRITIN_ELEVATED}")

        if tsat_elevated and ferritin_elevated:
            score = sum([family_history, symptoms, past_history, abnormal_labs])
            if score >= 3:
                category = "Likely Disease Onset"
            elif score >= 1:
                category = "Early Pattern" #the score=1 covers early pattern w/o user context (abnormal labs only)
            else:
                category = "Typical"
            return [{
                "Condition": "Hereditary Hemochromatosis", "Category": category,
                "Triggering Genetics": gene_trigger,
                "Triggering Parameters": ", ".join(triggered),
            }]

        return [{
            "Condition": "Hereditary Hemochromatosis", "Category": "Typical",
            "Triggering Genetics": gene_trigger,
            "Triggering Parameters": ", ".join(triggered),
        }]

    if iron_profile_done is False:
        score = sum([family_history, symptoms, past_history, abnormal_labs])
        if score >= 3:
            category = "Likely Disease Onset"
        elif score >= 1:
            category = "Early Pattern" #the score=1 covers early pattern w/o user context (abnormal labs only)
        else:
            category = "Typical"
        return [{
            "Condition": "Hereditary Hemochromatosis", "Category": category,
            "Triggering Genetics": gene_trigger,
            "Triggering Parameters": ", ".join(triggered),
        }]
    return [{
        "Condition": "Hereditary Hemochromatosis",
        "Category": "Relevant lab parameters not found",
        "Triggering Genetics": gene_trigger,
        "Triggering Parameters": ", ".join(triggered),
    }]


# Celiac Disease
def evaluate_celiac_disease(labs, patient, genetics, family_history, past_history, symptoms=False):
    if not has_flagged_gene(genetics, "Celiac Disease"):
        return [{"Condition": "Celiac Disease", "Category": GENE_NOT_FOUND}]

    vitamin_d = labs.get("vitamin_d")
    alt = labs.get("alt")
    ast = labs.get("ast")
    alp = labs.get("alp")
    ggt = labs.get("ggt")
    mcv = labs.get("mcv")
    mch = labs.get("mch")
    mchc = labs.get("mchc")

    triggered = []
    parameter = (
        is_below(vitamin_d, VITAMIN_D_DEFICIENT) and
        is_above(alt, ALT_RANGE[1]) and
        is_above(ast, AST_UPPER_STRICT) and
        is_above(ggt, GGT_LFT_RANGE[1]) and
        is_above(alp, ALP_RANGE[1]) and
        is_below(mcv, MCV_LOWER_CELIAC) and
        is_below(mch, MCH_LOWER_CELIAC) and
        is_below(mchc, MCHC_LOWER_CELIAC)
    )
    # all 8 sub-checks are AND'd into a single pattern, so only record them once the
    # whole pattern is confirmed true — otherwise a false prefix would leave misleading
    # partial entries in `triggered` for a pattern that didn't actually fire
    if parameter:
        triggered.append(f"vitamin_d={vitamin_d} (<{VITAMIN_D_DEFICIENT})")
        triggered.append(f"alt={alt} (>{ALT_RANGE[1]})")
        triggered.append(f"ast={ast} (>{AST_UPPER_STRICT})")
        triggered.append(f"ggt={ggt} (>{GGT_LFT_RANGE[1]})")
        triggered.append(f"alp={alp} (>{ALP_RANGE[1]})")
        triggered.append(f"mcv={mcv} (<{MCV_LOWER_CELIAC})")
        triggered.append(f"mch={mch} (<{MCH_LOWER_CELIAC})")
        triggered.append(f"mchc={mchc} (<{MCHC_LOWER_CELIAC})")

    if (
        (parameter and family_history and symptoms and past_history) or
        (parameter and family_history and past_history) or
        (parameter and family_history) or
        (family_history and symptoms) or
        (family_history and past_history)
    ):
        category = "Likely Disease Onset"
    elif (
        family_history or
        (parameter and symptoms and past_history) or
        (parameter and past_history) or
        parameter #added for early pattern w/o user context
    ):
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "Celiac Disease", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Celiac Disease"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Familial Hypercholesterolemia
def evaluate_fh(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Familial Hypercholesterolemia"):
        return [{"Condition": "Familial Hypercholesterolemia", "Category": GENE_NOT_FOUND}]

    ldl = labs.get("ldl_c")
    non_hdl = labs.get("non_hdl_c")
    age = patient.get("age")

    triggered = []
    elevated_160 = (
        note(triggered, "ldl_c", ldl, is_above(ldl, LDL_C_ELEVATED), f">{LDL_C_ELEVATED}")
        or note(triggered, "non_hdl_c", non_hdl, is_above(non_hdl, NON_HDL_C_ELEVATED), f">{NON_HDL_C_ELEVATED}")
    )
    elevated_190 = (
        note(triggered, "ldl_c", ldl, is_above(ldl, LDL_C_SEVERE), f">{LDL_C_SEVERE}")
        or note(triggered, "non_hdl_c", non_hdl, is_above(non_hdl, NON_HDL_C_SEVERE), f">{NON_HDL_C_SEVERE}")
    )

    if elevated_190 or (elevated_160 and age is not None and age < 40) or (elevated_160 and family_history):
        category = "Likely Disease Onset"
    elif elevated_160 or family_history:
        category = "Early Pattern"
    else:
        category = "Typical"

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
    eosinophils_absolute = labs.get("eosinophils_absolute")
    sex = patient.get("sex")

    HIGH_SPECIFICITY_FINDINGS = [] #to be defined as a subset of imaging findings
    
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
    eos_high = note(triggered, "eosinophils_absolute", eosinophils_absolute, is_elevated(eosinophils_absolute, 1500), f">={1500}")
    labs_severe = hb_severely_low or eos_high

    gene_trigger = get_gene_trigger(genetics, "Cardiomyopathy")

    if hemoglobin is None and eosinophils_absolute is None and not imaging_performed:
        return [{
            "Condition": "Cardiomyopathy",
            "Category": "Relevant lab/imaging parameters not found",
            "Triggering Genetics": gene_trigger,
            "Triggering Parameters": ", ".join(triggered),
        }]

    # --- Case 1: imaging performed and abnormal -> imaging drives the call ---
    if imaging_performed and imaging_abnormal:
        if imaging_high_spec or imaging_multiple or symptoms or family_history:
            category = "Likely Disease Onset"
        else:
            category = "Early Pattern"

    # --- Case 2: imaging performed and normal ---
    elif imaging_performed and not imaging_abnormal:
        risk_factor_count = sum([family_history, symptoms, labs_severe])
        if risk_factor_count >= 2:
            category = "Re-imaging recommended (Echocardiogram/Cardiac MRI)"
        else:
            category = "Typical"

    # --- Case 3: no imaging performed -> route to imaging ---
    else:
        if labs_severe or family_history or symptoms:
            category = "Advanced cardiac imaging recommended (Echocardiogram/Cardiac MRI)"
        else:
            category = "Typical"

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
    age = patient.get("age")
    sex = patient.get("sex")

    hdl_min = sex_based_threshold(sex, HDL_C_MALE_MIN, HDL_C_FEMALE_MIN)

    triggered = []
    ldl_moderate = (
        note(triggered, "ldl_c", ldl, is_elevated(ldl, LDL_C_AT_RISK), f">={LDL_C_AT_RISK}")
        or note(triggered, "non_hdl_c", non_hdl, is_elevated(non_hdl, NON_HDL_C_ELEVATED), f">={NON_HDL_C_ELEVATED}")
    )
    low_hdl = note(triggered, "hdl_c", hdl, is_below(hdl, hdl_min), f"<{hdl_min}")
    tg_moderate = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")
    prediabetes = note(triggered, "hba1c", hba1c, is_elevated(hba1c, HBA1C_AT_RISK), f">={HBA1C_AT_RISK}")
    any_moderate = ldl_moderate or low_hdl or tg_moderate or prediabetes

    ldl_high = (
        note(triggered, "ldl_c", ldl, is_elevated(ldl, LDL_C_ELEVATED), f">={LDL_C_ELEVATED}")
        or note(triggered, "non_hdl_c", non_hdl, is_elevated(non_hdl, NON_HDL_C_SEVERE), f">={NON_HDL_C_SEVERE}")
    )
    diabetes = note(triggered, "hba1c", hba1c, is_elevated(hba1c, HBA1C_ELEVATED), f">={HBA1C_ELEVATED}")
    combined_atherogenic = ldl_moderate and low_hdl and tg_moderate

    if (
        ldl_high or
        diabetes or
        combined_atherogenic or
        (any_moderate and symptoms and age is not None and age < 40) or
        (any_moderate and symptoms and family_history)
    ):
        category = "Likely Disease Onset"
    # (any_moderate + age) added for early pattern w/o user context (labs only)
    elif family_history or (any_moderate and symptoms) or (any_moderate and age is not None and age < 40):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    hdl_low = note(triggered, "hdl_c", hdl, is_below(hdl, hdl_normal), f"<{hdl_normal}")
    hdl_severely_low = note(triggered, "hdl_c", hdl, is_below(hdl, sex_based_threshold(sex, HDL_C_MALE_SEVERE, HDL_C_FEMALE_SEVERE)), f"<{sex_based_threshold(sex, HDL_C_MALE_SEVERE, HDL_C_FEMALE_SEVERE)}")

    tg_borderline = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")
    tg_moderate = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_ELEVATED), f">={TRIGLYCERIDES_ELEVATED}")
    tg_severe = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_SEVERE), f">={TRIGLYCERIDES_SEVERE}")
    tg_very_severe = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_VERY_SEVERE), f">={TRIGLYCERIDES_VERY_SEVERE}")

    vldl_elevated = note(triggered, "vldl", vldl, is_elevated(vldl, VLDL_UPPER_NORMAL), f">={VLDL_UPPER_NORMAL}")
    any_moderate_abnormal = tg_moderate or hdl_low or vldl_elevated

    if (
        tg_very_severe or
        (tg_severe and age is not None and age < 40) or
        (tg_severe and family_history) or
        (tg_severe and hdl_severely_low)
    ):
        category = "Likely Disease Onset"
    elif (
        tg_severe or
        any_moderate_abnormal or
        (tg_borderline and hdl_low) or
        family_history
    ):
        category = "Early Pattern"
    else:
        category = "Typical"

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
    hdl_low = note(triggered, "hdl_c", hdl, is_below(hdl, hdl_normal), f"<{hdl_normal}")
    hdl_severe = note(triggered, "hdl_c", hdl, is_below(hdl, sex_based_threshold(sex, 25, 30)), f"<{sex_based_threshold(sex, 25, 30)}")

    tg_elevated = note(triggered, "triglycerides", triglycerides, is_elevated(triglycerides, TRIGLYCERIDES_AT_RISK), f">={TRIGLYCERIDES_AT_RISK}")
    ldl_elevated = note(triggered, "ldl_c", ldl, is_elevated(ldl, LDL_C_ELEVATED), f">={LDL_C_ELEVATED}")

    isolated_significant = hdl_severe and not tg_elevated and not ldl_elevated
    isolated_moderate = hdl_low and not tg_elevated and not ldl_elevated
    mixed_moderate = hdl_low and (tg_elevated or ldl_elevated)

    if (
        isolated_significant or
        (isolated_moderate and ((age is not None and age < 40) or family_history))
    ):
        category = "Likely Disease Onset"
    elif isolated_moderate or mixed_moderate or family_history:
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "HDL Deficiency", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "HDL Deficiency"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Elevated Apolipoprotein B
def evaluate_apob(labs, patient, genetics, family_history, symptoms=False):
    if not has_flagged_gene(genetics, "Elevated Apolipoprotein B"):
        return [{"Condition": "Elevated Apolipoprotein B", "Category": GENE_NOT_FOUND}]

    apob = labs.get("apob")
    age = patient.get("age")

    triggered = []
    elevated_100 = note(triggered, "apob", apob, is_elevated(apob, 100), f">={100}")
    elevated_130 = note(triggered, "apob", apob, is_elevated(apob, 130), f">={130}")

    if elevated_130 or (elevated_100 and age is not None and age < 40) or (elevated_100 and family_history):
        category = "Likely Disease Onset"
    elif elevated_100 or family_history:
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "Elevated Apolipoprotein B", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "Elevated Apolipoprotein B"),
        "Triggering Parameters": ", ".join(triggered),
    }]


# Li-Fraumeni Syndrome
def evaluate_li_fraumeni_syndrome(labs, patient, genetics, family_history, symptoms=False, conditions=None):
    if not has_flagged_gene(genetics, "Li-Fraumeni Syndrome"):
        return [{"Condition": "Li-Fraumeni Syndrome", "Category": GENE_NOT_FOUND}]

    age = patient.get("age")
    conditions = conditions or []

    pathognomonic_tumors = [
        "Adrenocortical carcinoma", "Choroid plexus carcinoma",
        "Embryonal rhabdomyosarcoma, anaplastic subtype"
    ]
    lfs_spectrum_tumors = [
        "Sarcoma", "Osteosarcoma", "Breast cancer", "Brain tumor",
        "Adrenocortical carcinoma", "Leukemia", "Choroid plexus carcinoma",
        "Embryonal rhabdomyosarcoma, anaplastic subtype"
    ]

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
        category = "Likely Disease Onset"
    elif has_lfs_tumor or family_history:
        category = "Early Pattern"
    else:
        category = "Typical"

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

    major_count = sum([
        "Macrocephaly" in conditions,
        "Multiple trichilemmomas" in conditions,
        "Multiple oral papillomas" in conditions,
        "Multiple gastrointestinal hamartomas" in conditions,
    ])
    minor_count = sum([
        "Autism spectrum disorder" in conditions,
        "Developmental delay" in conditions,
        "Vascular anomaly" in conditions,
        "Gastrointestinal polyps" in conditions,
        "Intellectual disability" in conditions,
        "Thyroid nodule" in conditions,
        "Lipoma" in conditions,
    ])

    macrocephaly = "Macrocephaly" in conditions
    gi_hamartoma = "Multiple gastrointestinal hamartomas" in conditions

    macrocephaly_plus = (
        macrocephaly and
        any([
            "Autism spectrum disorder" in conditions,
            "Developmental delay" in conditions,
            "Vascular anomaly" in conditions,
            "Gastrointestinal polyps" in conditions,
        ])
    )

    if (
        macrocephaly_plus or
        (major_count >= 3 and (macrocephaly or gi_hamartoma)) or
        (major_count >= 2 and minor_count >= 3)
    ):
        category = "Likely Disease Onset"
    elif major_count >= 1 or minor_count >= 2 or family_history:
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "PTEN Hamartoma Tumor Syndrome", "Category": category,
        "Triggering Genetics": get_gene_trigger(genetics, "PTEN Hamartoma Tumor Syndrome"),
        "Triggering Parameters": "",
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
    evaluate_men4(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, endocrine_conditions=None,
    ),
    evaluate_muscular_dystrophy(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_gilbert_syndrome(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, past_history=False,
    ),
    evaluate_hemochromatosis(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, past_history=False,
    ),
    evaluate_celiac_disease(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, past_history=False,
    ),
    evaluate_fh(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_cardiomyopathy(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, imaging_performed=False, imaging_findings=None,
    ),
    evaluate_cad(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_hypertriglyceridemia(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_hdl_deficiency(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_apob(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False,
    ),
    evaluate_li_fraumeni_syndrome(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, conditions=None,
    ),
    evaluate_phts(
        labs=data["labs"], patient=data["patient"], genetics=data["genetics"],
        family_history=False, conditions=None,
    ),
]

active_findings = [
    r[0]
    for r in results
    if r[0]["Category"] not in {"Related mutation in gene(s) not found"}
]

with open(args.output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Condition", "Category", "Triggering Genetics", "Triggering Parameters"])
    writer.writeheader()
    writer.writerows(active_findings)

print(f"Results written to {args.output_csv}")
