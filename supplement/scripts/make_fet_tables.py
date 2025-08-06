from os.path import dirname, abspath, join
from csv import DictReader
from enum import Enum
from collections import defaultdict
from mylatextable import MyLongTable, MyLatexTable
import typing
from util import format_p
from make_measurement_tables import get_fet_data
from analysis import ANALYSIS_VERSION
# Create tables to summarize the total number of tests performed and their results
# We will use this for Table 1.

THIS_DIR = dirname(abspath(__file__))
SUPPLEMENT_DIR = join(THIS_DIR, ANALYSIS_VERSION)
SIG_FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "sig_fisher_exact_test_dashboard.txt")
SIG_FISHER_SUMMARY = join(THIS_DIR, "sig_fet_test_summary.txt")
FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "fisher_exact_test_dashboard.txt")
HPO_FET_OUT = join(THIS_DIR, "hpo_fet.tex")
DISEASE_TABLE_OUT = join(THIS_DIR, "disease_table.tex")
MF_TABLE_OUT = join(THIS_DIR, "mf_table.tex")

class TestType(Enum   ):
    HPO_TERM = 1
    DISEASE_COMPARISON = 2
    MF_COMPARISON = 3

def get_sig_fet_test_results():
    """
    Collect each of the main tests types into lists that we will use for output
    Returns:
    hpo_terms (Dict[str,str]): Dictionary with key=parameter name, e.g., cohort_name, a_genotypes, for HPO Fisher Exact Tests results
    disease_terms (Dict[str,str]): same as above, but for disease FET results
    mf_comparisons (Dict[str,str]): same as above, but for male/female comparison FET results.
    """
    hpo_terms = list()
    disease_terms = list()
    mf_comparisons = list()

    with open(SIG_FISHER_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            with_geno_a = row["a_genotype"]
            if "OMIM" in with_geno_a:
                disease_terms.append(row)
            elif "FEMALE" in with_geno_a or "MALE" in with_geno_a:
                mf_comparisons.append(row)
            else:
                hpo_terms.append(row)
    return hpo_terms, disease_terms, mf_comparisons

def get_fisher_exact_test_descriptive_stats():
    """
    Get basic stats
    """
    cohort_to_sig = defaultdict(int)

    with open(FISHER_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            #print(row)
            cohort_name = row["#cohort_name"]
            nsig = int(row["nsig"])
            cohort_to_sig[cohort_name] += nsig
    return cohort_to_sig
    




def get_stats(name, list_of_dicts) -> typing.List[str]:
    terms = set()
    cohorts = set()
    for  d in list_of_dicts:
        #print(d)
        cohort = d["#cohort_name"]
        varname = d["hpo_item"]
        cohorts.add(cohort)
        terms.add(varname)
    items = [name, str(len(cohorts)), str(len(terms))]
    return items


     

def print_summary_table(hpo_terms, disease_terms, mf_comparisons):
    """
    Print a summary of the statistical analysis results for the Fisher Exact test with all analyzed cohorts
    
    Parameters:
    hpo_terms (Dict[str,str]): Dictionary with key=parameter name, e.g., cohort_name, a_genotypes, for HPO Fisher Exact Tests results
    disease_terms (Dict[str,str]): same as above, but for disease FET results
    mf_comparisons (Dict[str,str]): same as above, but for male/female comparison FET results.
    """
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
    print(f"Wrote {SIG_FISHER_SUMMARY}")



def create_sig_fisher_table(list_of_rows, 
                            header:typing.List[str],
                            header_format:str,
                            caption:str,
                            useLongTable=False, ):
    """
    #cohort_name	a_genotype	b_genotype	nsig	n_tests_performed	hpo_item	with_geno_a	with_geno_b	pval	adj_pval
    """
    if useLongTable:
        table = MyLongTable(header_fields=header, use_booktabs=True, header_field_formats=header_format, fontsize="scriptsize", caption=caption)
    else:
        table = MyLatexTable(header_fields=header, use_booktabs=True, header_field_formats=header_format, fontsize="scriptsize", caption=caption)
    
    for row in list_of_rows:
        cohort = row["#cohort_name"]
        hpo_item = row["hpo_item"]
        a_genotype = row["a_genotype"]
        b_genotype = row["b_genotype"]
        with_geno_a = row["with_geno_a"].replace("%", "\\%")
        with_geno_b = row["with_geno_b"].replace("%", "\\%")
        pval = format_p(row["pval"])
        adj_pval = format_p(row["adj_pval"])
        items = [cohort, hpo_item,  a_genotype, with_geno_a, 
                  b_genotype , with_geno_b, pval, adj_pval]
        table.add_row(items)
    return table.get_latex()


def print_mf_table(mf_rows):
    header = ["cohort", "HPO", "genotype (A)", "Counts (A)",  "genotype (B)", "Counts (B)", "p-val", "adj. p"]
    header_field_formats = "l>{\\raggedright}p{2.5cm}llllll" # need to import array package for raggedright
    caption = """Fischer exact test for association between phenotypic features and sex (male, female)."""
    table = create_sig_fisher_table(list_of_rows=mf_rows,
                               header=header,
                               header_format=header_field_formats,
                               caption=caption)
                               
    fh = open(MF_TABLE_OUT, "wt")
    fh.write(table)
    fh.close()
    print(f"Wrote{MF_TABLE_OUT}")

def print_disease_table(disease_rows):
    header = ["Cohort", "HPO", "disease A", "",  "disease B", "", "p-val", "adj. p"]
    header_field_formats = "l>{\\raggedright}p{2.5cm}llllll" # need to import array package for raggedright
    caption = """Fischer exact test for association between disease diagnosis and phenotypic features."""
    table = create_sig_fisher_table(list_of_rows=disease_rows,
                               header=header,
                               header_format=header_field_formats,
                               caption=caption)
    fh = open(DISEASE_TABLE_OUT, "wt")
    fh.write(table)
    fh.close()
    print(f"Wrote{DISEASE_TABLE_OUT}")

def print_hpo_table(hpo_terms):
    header = ["Cohort", "HPO", "Genotype A", "",  "Genotype B", "", "p-val", "adj. p"]
    header_field_formats = "l>{\\raggedright}p{2.5cm}>{\\raggedright}p{1.5cm}l>{\\raggedright}p{1.5cm}lll" # need to import array package for raggedright
    caption = """Fischer exact test for association between genotypes and phenotypic features."""
    table = create_sig_fisher_table(list_of_rows=hpo_terms,
                               header=header,
                               header_format=header_field_formats,
                               useLongTable=True,
                               caption=caption)
    print(f"Table with {len(hpo_terms)} HPO entry lines")
    fh = open(HPO_FET_OUT, "wt")
    fh.write(table)
    fh.close()
    print(f"Wrote {HPO_FET_OUT}")


   
if __name__ == "__main__":
    print("Make Fisher Exact Tables")
    cohort_to_sig = get_fisher_exact_test_descriptive_stats() 
    print(cohort_to_sig)
    get_fet_data()## this creates updated versions of mono_test summary.tex and other files
    hpo_terms, disease_terms, mf_comparisons = get_sig_fet_test_results()
    print_summary_table(hpo_terms, disease_terms, mf_comparisons)
    print_mf_table(mf_rows=mf_comparisons)
    print_disease_table(disease_rows=disease_terms)
    print_hpo_table(hpo_terms=hpo_terms)
    
