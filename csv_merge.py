import csv
import sys

# Preference order (lower index = higher priority) ------------------------------------

# "Typical" from csv1 (gene + blood) is treated as "Typical_gene",
# "Typical" from csv2 (blood) is treated as "Typical_blood" - only for ranking purposes.

PREFERENCE = ["Significant Pattern", "Early Pattern", "Elevated susceptibility", "Typical_blood", "Typical_gene"]

BLOOD_MARKER_COLS = ["Blood Marker(s)", "All Blood Marker(s)"]

def rank(category, source):
    """Returns priority rank for category, and tags 'Typical' by source csv"""
    if category == "Typical":
        category = f"Typical_{source}" #source is "gene" or "blood"
    if category not in PREFERENCE:
        raise ValueError(f"Unrecognized category '{category}' from {source} CSV - not in PREFERENCE list: {PREFERENCE}")
    return PREFERENCE.index(category)

def load_csv(path):
    """Returns a dictionary: Condition : row dict."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return {row["Condition"]: row for row in reader}, reader.fieldnames

def _is_empty(val):
    return val is None or val.strip() == "" or val.strip().lower() == "none"

def apply_cad_override(condition, chosen_row, other_row):
    """
    Special case: Coronary Artery Disease + Elevated susceptibility ->
    for each column in BLOOD_MARKER_COLS, if the chosen row's value is empty/None,
    fill it from the other row's value for that same column.
    """
    if condition != "Coronary Artery Disease":
        return chosen_row
    if chosen_row.get("Category") != "Elevated susceptibility":
        return chosen_row
    if other_row is None:
        return chosen_row

    for col in BLOOD_MARKER_COLS:
        if _is_empty(chosen_row.get(col, "")):
            other_val = other_row.get(col, "").strip()
            if not _is_empty(other_val):
                chosen_row[col] = other_val

    return chosen_row

def main(gene_csv_path, blood_csv_path, output_csv_path):
    gene_rows, fieldnames = load_csv(gene_csv_path)
    blood_rows, _ = load_csv(blood_csv_path) #(fine since both have same header, if header changes, this needs to change)

    final_rows = []

    for condition, gene_row in gene_rows.items():
        blood_row = blood_rows.get(condition)

        if blood_row is None:
            #condition not present in blood_csv -> use from gene_csv as-is
            chosen_row = apply_cad_override(condition, gene_row, None)
            final_rows.append(chosen_row)
            continue

        gene_rank = rank(gene_row["Category"], "gene")
        blood_rank = rank(blood_row["Category"], "blood")

        # lower rank means higher preference
        if gene_rank <= blood_rank:
            chosen_row, other_row = gene_row, blood_row
        else:
            chosen_row, other_row = blood_row, gene_row

        chosen_row = apply_cad_override(condition, chosen_row, other_row)
        final_rows.append(chosen_row)

    with open(output_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_rows)

    print(f"Final results written to {output_csv_path}")

if __name__=="__main__":
    if len(sys.argv) != 4:
        print("Usage: python csv_merge.py <gene_csv> <blood_csv> <output_csv>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])