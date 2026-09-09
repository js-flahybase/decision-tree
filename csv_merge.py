import csv
import sys

# Preference order (lower index = higher priority) ------------------------------------

# "Typical" from csv1 (gene + blood) is treated as "Typical_gene",
# "Typical" from csv2 (blood) is treated as "Typical_blood" - only for ranking purposes.

PREFERENCE = ["Significant Pattern", "Early Pattern", "Elevated susceptibility", "Typical_blood", "Typical_gene"]

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

def main(gene_csv_path, blood_csv_path, output_csv_path):
    gene_rows, fieldnames = load_csv(gene_csv_path)
    blood_rows, _ = load_csv(blood_csv_path) #(fine since both have same header, if header changes, this needs to change)

    final_rows = []

    for condition, gene_row in gene_rows.items():
        blood_row = blood_rows.get(condition)

        if blood_row is None:
            #condition not present in blood_csv -> use from gene_csv as-is
            final_rows.append(gene_row)
            continue

        gene_rank = rank(gene_row["Category"], "gene")
        blood_rank = rank(blood_row["Category"], "blood")

        # lower rank means higher preference
        final_rows.append(gene_row if gene_rank <= blood_rank else blood_row)

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