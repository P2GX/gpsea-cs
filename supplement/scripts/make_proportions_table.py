import csv
import os
from analysis import ANALYSIS_VERSION
from mylatextable import  MyLatexTable

# Create tables to summarize the total number of tests performed and their results
# We will use this for Table 1.

from analysis import ANALYSIS_VERSION

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# TODO(pnrobinson) - add the file into the repo.
DISTRIBUTION_FILE = os.path.join(THIS_DIR, f"../generated/{ANALYSIS_VERSION}/distribution-hpo.txt")

GENERATED_DIR = os.path.join(THIS_DIR, os.pardir, 'generated', ANALYSIS_VERSION)
PROPORTIONS_FILE = os.path.join(GENERATED_DIR, "proportions.tex")



header = ["HPO", "id", "Count", "Observed (\\%)", "Expected (\\%)" ]
caption = "Distribution of significant Fisher exact test results according to top-level HPO term"
table = MyLatexTable(header_fields=header, use_booktabs=True,  caption=caption)
with open(DISTRIBUTION_FILE) as f:
    filtered_lines = (line for line in f if not line.lstrip().startswith("!"))
    reader = csv.DictReader(filtered_lines, delimiter="\t")
    for row in reader:
        hpo = row["HPO"]
        ident = row["ID"]
        count = str(row["Count"])
        o = row["Observed"].replace("%","\\%")
        e = row["Expected"].replace("%","\\%")
        items = [hpo, ident, count, o, e]
        table.add_row(items)

with open(PROPORTIONS_FILE, "wt") as fh:
    fh.write(table.get_latex())

