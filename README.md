# Genetic Risk Evaluation Functions

## What is this?

A set of Python functions — one per condition — that evaluate a patient's risk for genetic condition onset.

Each function looks at three things:

1. **Was the relevant gene flagged, or is the polygenic risk score (PRS) elevated, for this condition?**
2. **Are the relevant blood test values abnormal?**
3. **Does the patient have relevant symptoms, family history, or past history?**

Based on these inputs, each function returns one simple result: a **category indicating how strong the signal is**.

---

## Standard Function Structure

Every function follows the same general pattern:

```python
def evaluate_<condition>(labs, patient, genetics, family_history, symptoms=False, ...):
    # 1. Check gene/PRS first — nothing else matters if neither is present
    if not has_flagged_gene(genetics, "<Condition Name>"):
        return [{
            "Condition": "<Condition Name>",
            "Category": GENE_NOT_FOUND
        }]

    # 2. Pull values from labs / patient
    some_value = labs.get("some_value")

    # 3. Check if values are abnormal
    some_flag = is_elevated(some_value, SOME_THRESHOLD)

    # 4. Decide category
    if <strong combination met>:
        category = "Significant Pattern"
    elif <partial combination met>:
        category = "Early Pattern"
    else:
        category = "Elevated Susceptibility"

    return [{
        "Condition": "<Condition Name>",
        "Category": category
    }]
```

`GENE_NOT_FOUND` is set to `"Typical"` — meaning a condition whose gene/PRS gate never opened is labeled the same as a true "nothing to act on" result, since neither lab nor context data was ever checked for it.

---

## Gene Gate: Monogenic Flag OR Polygenic (PRS) Elevation

`has_flagged_gene(genetics, condition_name)` decides whether a condition even gets evaluated, and it now passes if **either** of these is true:

1. **Monogenic hit** — a specific gene tied to this condition (via `CONDITION_GENES`) was flagged as Pathogenic/Likely Pathogenic.
2. **Polygenic (PRS) elevation** — this condition has a polygenic risk score entry, and its category is `elevated` or `moderately_elevated`.

If neither a gene nor a PRS signal is present, the function returns immediately with `GENE_NOT_FOUND` (Typical) — no lab or symptom data is checked at all.

---

## Output Columns

Each result includes these columns:

- **Domain** — The broad health category this condition falls under (e.g. Cardiac Health, Endocrine Health, Digestive Health)
- **Condition** — The condition name
- **Category** — `Typical`, `Elevated Susceptibility`, `Early Pattern`, `Significant Pattern`
- **DNA Marker(s)** — which gene(s) entry actually passed the gate for that condition (e.g. `HFE`)
- **Blood Marker(s)** — which specific lab values crossed a threshold to support `Category` assignment, and by how much (e.g. `crp=3.2 (>=3.0)`)
- **Triggering PRS** — Shown when the condition's PRS category is `elevated` or `moderately_elevated` and contributed to passing the gate
- **Snapshot Category** — A plain-language recommendation: `To discuss with General Practitioner`, `Worth acting on for prevention`, or `Typical - nothing to act on`

This makes it possible to see *why* a category was reached, not just the final result.

---

## Five Types of Input

| Input | Where it comes from | Type | Example |
|---|---|---|---|
| **Blood test values (`labs`)** | Blood report CSV | Numbers | `hba1c: 6.8`, `ldl_c: 175` |
| **Monogenic gene flags** | Monogenic (ClinVar-style) JSON — Pathogenic/Likely Pathogenic entries | List of gene names | `["HFE", "TP53"]` |
| **Polygenic (PRS) scores** | Polygenic risk score JSON | Per-condition category (`typical` / `elevated` / `moderately_elevated`) | `"cad": {"category": "elevated"}` |
| **APOE status** | APOE-specific JSON | Genotype string | `"e3/e4"` — contributes an `APOE` gene flag if an ε4 allele is present |
| **Patient-reported context** (`family_history`, `symptoms`, `past_history`) | Asked to / reported by the patient | `True` / `False` | `family_history of Alzheimer's: True` |

`genetics` (passed into every function) is a single dict combining all three genetic sources:
```python
{
    "flagged_genes": [...],            # from the monogenic JSON + APOE (if ε4 present)
    "acmg_genes": [...],                 # subset of flagged_genes also listed as ACMG-actionable
    "prs_elevated_conditions": {...},  # set of condition keys from the PRS JSON where category is elevated
    "prs_categories": {...},             # raw PRS category string per condition key
    "apoe_status": "e3/e4",            # raw APOE genotype string
}
```
This dict is built once per patient by `build_genetics_from_jsons(...)`, not assembled inside each function.

---
 
## Patient Sex and Age
 
`sex` and `age` are read from the blood report CSV by default, using rows like:
```
parameter,value
sex,male
age,35
```
The `--sex`/`--age` command-line flags act as **overrides**: if provided, they take priority over whatever is in the CSV. If neither the CSV nor the command line provides a value, the script raises an error rather than guessing.

---

## Why Family History, Symptoms, and Past History Are Boolean

The patient does not provide a laboratory value for something like **"Do you have a family history of X?"** They either report it or they do not.

By the time these values reach the function, they are therefore simplified into a single boolean:

- `True` → The patient reported at least one relevant item.
  - Example: They reported a family history of diabetes or sudden cardiac death in the family.
- `False` → They did not report any relevant item, or they did not answer.

The function itself does not need to know which specific item was reported. It only needs to know whether that category applies to the condition.

---

## Helper Functions

The following helper functions are used throughout the evaluation functions:

| Helper | What it checks |
|---|---|
| `is_elevated(value, threshold)` | `value >= threshold` |
| `is_above(value, threshold)` | `value > threshold` (strict) |
| `is_below(value, threshold)` | `value < threshold` |
| `is_outside_range(value, low, high)` | `value < low` or `value > high` |

All four helpers automatically return `False` when the value is missing (`None`).

Therefore, a missing laboratory value **never causes the function to crash**; it simply does not count as an abnormal result.

---

## Running the Script

```bash
python script.py <monogenic.json> <prs.json> <apoe.json> <blood.csv> <output.csv>
```
Optional overrides (only needed if the blood CSV doesn't include `sex`/`age` rows, or we want to override them):
```bash
--sex {male, female} --age int(value)
```

This runs all 25 condition functions against the supplied data and writes result to `output.csv`.

---

## Summary

The evaluation framework follows a consistent pattern:

**Genetic flag (monogenic OR polygenic) → Relevant clinical/laboratory data → Clinical context → Risk category**

This allows each condition-specific function to use the same overall structure while applying condition-specific thresholds and combinations of clinical evidence.