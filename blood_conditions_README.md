# Genetic Risk Evaluation Functions

## What is this?

A set of Python functions — one per condition — that evaluate a patient's risk for genetic condition onset.

Each function looks at three things:

1. **Was the relevant gene flagged in the patient's genetic report?**
2. **Are the relevant blood test values abnormal?**
3. **Does the patient have relevant symptoms, family history, or past history?**

Based on these inputs, each function returns one simple result: a **category indicating how strong the signal is**.

---

## Standard Function Structure

Every function follows the same general pattern:

```python
def evaluate_<condition>(labs, patient, genetics, family_history, symptoms=False, ...):
    # 1. Check gene first — nothing else matters if it's not present
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
        category = "Likely Disease Onset"
    elif <partial combination met>:
        category = "Early Pattern"
    else:
        category = "Typical"

    return [{
        "Condition": "<Condition Name>",
        "Category": category
    }]
```

---

## Three Types of Input

| Input | Where it comes from | Type | Example |
|---|---|---|---|
| **Blood test values (`labs`)** | Lab report PDF | Numbers | `hba1c: 6.8`, `ldl_c: 175` |
| **Gene flags (`genetics`)** | Genetic report PDF | List of gene names | `["HFE", "TP53"]` |
| **Patient-reported context** (`family_history`, `symptoms`, `past_history`) | Asked to / reported by the patient | `True` / `False` | `family_history of Alzheimers: True` |

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

## Summary

The evaluation framework follows a consistent pattern:

**Genetic flag → Relevant clinical/laboratory data → Clinical context → Risk category**

This allows each condition-specific function to use the same overall structure while applying condition-specific thresholds and combinations of clinical evidence.