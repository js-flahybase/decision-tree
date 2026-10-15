import csv
import sys

# Preference order (lower index = higher priority) ------------------------------------

# "Typical" from csv1 (gene + blood) is treated as "Typical_gene",
# "Typical" from csv2 (blood) is treated as "Typical_blood" - only for ranking purposes.

# PREFERENCE = ["Significant Pattern", "Early Pattern", "Elevated susceptibility", "Typical_blood", "Typical_gene"]

SIG, EARLY, ELEV, TYP = "significant pattern", "early pattern", "elevated susceptibility", "typical"
VALID = {SIG, EARLY, ELEV, TYP}

BLOOD_MARKER_COLS = ["Blood Marker(s)", "All Blood Marker(s)"]

# def rank(category, source):
#     """Returns priority rank for category, and tags 'Typical' by source csv"""
#     if category == "Typical":
#         category = f"Typical_{source}" #source is "gene" or "blood"
#     if category not in PREFERENCE:
#         raise ValueError(f"Unrecognized category '{category}' from {source} CSV - not in PREFERENCE list: {PREFERENCE}")
#     return PREFERENCE.index(category)

# to normalize category strings for comparison
def _norm(cat):
    return (cat or "").strip().lower()

def load_csv(path):
    """Returns a dictionary: Condition : row dict."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return {row["Condition"]: row for row in reader}, reader.fieldnames

def _is_empty(val):
    return val is None or val.strip() == "" or val.strip().lower() == "none"

GENETIC_COLS = ["RSID_m", "Pathogenicity_m", "Elevation_p", "Threshold_p", "APOE_Status"]

def merge_rows(condition, gene_row, blood_row):
    """Returns the row to keep. Raises ValueError for impossible combinations."""
    g, b = _norm(gene_row["Category"]), _norm(blood_row["Category"])
    for cat, src in ((g, "gene"), (b, "blood")):
        if cat not in VALID:
            raise ValueError(f"{condition}: unrecognized category '{cat}' in {src} CSV")

    if g == TYP:
        for col in GENETIC_COLS:
            blood_row[col] = gene_row.get(col,"")
        return blood_row #"Typical" in gene csv is always replaced (it was gene gate filtered out, so no info)
    if g == SIG:
        return gene_row #"Significant Pattern" from gene csv is at highest priority
    if g == EARLY:
        if b == SIG:
            raise ValueError(f"{condition}: gene Early Pattern but blood Significant Pattern") #error case, should not happen, still a fallback check to prevent slipping out any
        return gene_row
    # only g == ELEV cases left now (anything else handled in unrecognized category error or abouve cases)
    if b in (SIG, EARLY):
        raise ValueError(f"{condition}: gene Elevated susceptibility and blood {b}") #error case, should not happen, still a fallback check to prevent slipping out any
    if b == ELEV:
        for col in BLOOD_MARKER_COLS:
            if _is_empty(gene_row.get(col)):
                gene_row[col] = blood_row.get(col, "") #if both Elevated susceptibility and cols are empty, replace "Blood Marker(s)" and "All Blood Marker(s)" columns with blood csv's cols
    elif b == TYP and not _is_empty(blood_row.get("All Blood Marker(s)")):
        if _is_empty(gene_row.get("All Blood Marker(s)")):
            gene_row["All Blood Marker(s)"] = blood_row["All Blood Marker(s)"] #only add "All Blood Markers" from blood's Typical if it is non empty and gene's entry is empty
    return gene_row

def main(gene_csv_path, blood_csv_path, output_csv_path):
    gene_rows, fieldnames = load_csv(gene_csv_path)
    blood_rows, _ = load_csv(blood_csv_path) #(fine since both have same header, if header changes, this needs to change)

    final_rows, errors = [], []

    for condition, gene_row in gene_rows.items():
        blood_row = blood_rows.get(condition)

        if blood_row is None:
            #condition not present in blood_csv -> use from gene_csv as-is
            final_rows.append(gene_row)
            continue
        try:
            final_rows.append(merge_rows(condition, gene_row, blood_row))
        except ValueError as e:
            errors.append(str(e))

    if errors:
            print("ERRORS (Should not happen):")
            for e in errors:
                print(" -", e)
            sys.exit(1)

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