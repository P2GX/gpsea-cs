from os.path import dirname, abspath, join
from csv import DictReader
from enum import Enum
# Create tables to summarize the total number of tests performed and their results
# We will use this for Table 1.

THIS_DIR = dirname(abspath(__file__))
SUPPLEMENT_DIR = dirname(THIS_DIR)
SIG_FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "sig_fisher_exact_test_dashboard.txt")
SIG_FISHER_SUMMARY = join(THIS_DIR, "sig_fet_test_summary.txt")

class TestType(Enum   ):
    HPO_TERM = 1
    DISEASE_COMPARISON = 2
    MF_COMPARISON = 3

def get_sig_fet_test_results():
    """
    Collect each of the main tests types into lists that we will use for output
    """
    hpo_terms = list()
    disease_terms = list()
    mf_comparisons = list()

    with open(SIG_FISHER_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            with_geno_a = row["a_genotype"]
            if "OMIM" in with_geno_a:
                test_type = TestType.DISEASE_COMPARISON
                disease_terms.append(row)
            elif "FEMALE" in with_geno_a or "MALE" in with_geno_a:
                test_type = TestType.MF_COMPARISON
                mf_comparisons.append(row)
            else:
                test_type = TestType.HPO_TERM
                hpo_terms.append(row)
    return hpo_terms, disease_terms, mf_comparisons




def get_stats(name, list_of_dicts) -> str:
    terms = set()
    cohorts = set()
    for  d in list_of_dicts:
        print(d)
        cohort = d["#cohort_name"]
        varname = d["hpo_item"]
        cohorts.add(cohort)
        terms.add(varname)
    items = [name, str(len(cohorts)), str(len(terms))]
    return items


     

def print_summary_table(hpo_terms, disease_terms, mf_comparisons):
    rows = list()
    header = ["Test", "Cohorts (n)", "Variables (n)"]
    rows.append(header)
    rows.append(get_stats("Phenotypic features", hpo_terms))
    rows.append(get_stats("Disease comparison", disease_terms))
    rows.append(get_stats("M/F comparisons", mf_comparisons))
    fh = open(SIG_FISHER_SUMMARY, "wt")
    for row in rows:
        line = "\t".join(row)
        fh.write(line + "\n")
    fh.close()

hpo_terms, disease_terms, mf_comparisons = get_sig_fet_test_results()
print_summary_table(hpo_terms, disease_terms, mf_comparisons)

