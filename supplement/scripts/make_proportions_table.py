from os.path import dirname, abspath, join
from csv import DictReader
from mylatextable import  MyLatexTable

# Create tables to summarize the total number of tests performed and their results
# We will use this for Table 1.

THIS_DIR = dirname(abspath(__file__))
DISTRIBUTION_FILE = join(THIS_DIR, "../../../hpotools/distribution-hpo.txt")
PROPORTIONS_FILE = join(THIS_DIR, "proportions.tex")



header = ["HPO", "id", "Count", "Observed (\\%)", "Expected (\\%)" ]
caption = "Distribution of significant Fisher exact test results according to top-level HPO term"
table = MyLatexTable(header_fields=header, use_booktabs=True,  caption=caption)
with open(DISTRIBUTION_FILE) as f:
    reader = DictReader(f, delimiter="\t")
    for row in reader:
        hpo = row["HPO"]
        ident = row["ID"]
        count = str(row["Count"])
        o = row["Observed"]
        e = row["Expected"]
        items = [hpo, ident, count, o, e]
        table.add_row(items)
fh = open(PROPORTIONS_FILE, "wt")
fh.write(table.get_latex())
fh.close()




