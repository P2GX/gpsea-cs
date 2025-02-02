from os.path import dirname, abspath, join
from csv import DictReader
from enum import Enum
from mylatextable import MyLatexTable
from util import format_p

# Create tables to summarize the total number of tests performed and their results
# For each of the tests except FET

THIS_DIR = dirname(abspath(__file__))
SUPPLEMENT_DIR = join(THIS_DIR, "v9_3")
MEASUREMENT_DASHBOARD = join(SUPPLEMENT_DIR, "measurement_dashboard.txt")
MEASUREMENT_SUMMARY = join(THIS_DIR, "mono_test_summary.txt")
T_TEST_OUTFILE = join(THIS_DIR, "t_test_table.tex")



class TestType(Enum   ):
    T_TEST = 1
    HPO_ONSET = 2
    DISEASE_ONSET = 3
    MORTALITY = 4
    PHENOTYPE_SCORE = 5

test_type_d = {
    "t-test": TestType.T_TEST,
    "Onset of Chronic mucocutaneous candidiasis": TestType.HPO_ONSET,
    "HPO Group Count": TestType.PHENOTYPE_SCORE,
    "Onset of OMIM:620465": TestType.DISEASE_ONSET,
    "De Vries Score": TestType.PHENOTYPE_SCORE,
    "Onset of Stage 5 chronic kidney disease": TestType.HPO_ONSET,
    "Onset of OMIM:248250": TestType.DISEASE_ONSET,
    "Onset of OMIM:610042": TestType.DISEASE_ONSET,
    "Onset of OMIM:130050": TestType.DISEASE_ONSET,
    "Age of death": TestType.MORTALITY,
    "Onset of OMIM:615471": TestType.DISEASE_ONSET,
    "Onset of OMIM:605911" : TestType.DISEASE_ONSET,
    "Onset of OMIM:256810": TestType.DISEASE_ONSET,
    "POGZ Severity Score": TestType.PHENOTYPE_SCORE,
    "Onset of Hypertrophic cardiomyopathy": TestType.HPO_ONSET,
    "Onset of OMIM:604377": TestType.DISEASE_ONSET,
    "Onset of OMIM:616831": TestType.DISEASE_ONSET,
    "Onset of OMIM:272300": TestType.DISEASE_ONSET,
    "Onset of OMIM:272300": TestType.DISEASE_ONSET,
}

def get_mono_test_results():
    """
    Collect each of the main tests types into lists that we will use for output
    """
    t_tests = list()
    hpo_onsets = list()
    disease_onsets = list()
    mortality = list()
    phenotype_scores = list()
    with open(MEASUREMENT_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            #print(row)
            test_name = row["test_name"]
            if test_name not in test_type_d:
                raise ValueError(f"Could not find test name {test_name}")
            test_category = test_type_d.get(test_name)
            if test_category == TestType.T_TEST:
                t_tests.append(row)
            elif test_category == TestType.HPO_ONSET:
                hpo_onsets.append(row)
            elif test_category == TestType.DISEASE_ONSET:
                disease_onsets.append(row)
            elif test_category == TestType.MORTALITY:
                mortality.append(row)
            elif test_category == TestType.PHENOTYPE_SCORE:
                phenotype_scores.append(row)
            else:
                raise ValueError(f"Did not recognize category {test_category}")
    return t_tests, hpo_onsets, disease_onsets,  mortality,  phenotype_scores


def get_stats(name, list_of_dicts) -> str:
    n_tests = 0
    n_sig_tests = 0
    terms = set()
    cohorts = set()
    for  d in list_of_dicts:
        cohort = d["#cohort"]
        varname = d["variable_name"]
        pval = float(d["pval"])
        n_tests += 1
        if pval <= 0.05:
            n_sig_tests += 1
        cohorts.add(cohort)
        terms.add(varname)
    items = [name, str(len(cohort)), str(len(terms)), str(n_tests), str(n_sig_tests)]
    return items
        

def print_summary_table(t_tests, hpo_onsets, disease_onsets,  mortality,  phenotype_scores):
    rows = list()
    header = ["Test", "Cohorts (n)", "Variables (n)", "Tests (n)", "Significant tests (n)"]
    rows.append(header)
    rows.append(get_stats("t test", t_tests))
    rows.append(get_stats("HPO Onset", hpo_onsets))
    rows.append(get_stats("Disease onset", disease_onsets))
    rows.append(get_stats("Mortality", mortality))
    rows.append(get_stats("Phenotype Scores", phenotype_scores))
    fh = open(MEASUREMENT_SUMMARY, "wt")
    for row in rows:
        line = "\t".join(row)
        fh.write(line + "\n")
    fh.close()



### Now create the detailed tables for the supplement
# 
def create_suppl_table(list_of_rows, header, header_format, caption):
    """
    These are the columns
    #cohort	a_genotype	b_genotype	name	test_name	description	variable_name	pval	interpretation	xrefs
    """
    rows = list()
    for row in list_of_rows:
        cohort = row["#cohort"]
        geno_a = row["a_genotype"]
        geno_b = row["b_genotype"]
        description = row["description"].replace("%", "\\%")
        xrefs = row["xrefs"]
        pval = format_p(row["pval"])
        items = [cohort, geno_a, geno_b, description, pval, xrefs ]
        rows.append(items)
    table = MyLatexTable(header_fields=header, 
                         header_field_formats=header_format, 
                         use_booktabs=True,  
                         fontsize="scriptsize",
                         caption=caption)
    for row in rows:
        table.add_row(row)
    return table.get_latex()


def create_t_test_table(t_tests):
    header = ["cohort",  "genotype (A)",   "genotype (B)", "Outcome Variable", "p-val", "xrefs", ]
    header_field_formats = "lp{2cm}p{2cm}>{\\raggedright}p{4.5cm}lr" # need to import array package for raggedright
    caption = """Student t-tests performed using GPSEA. 1/1: biallelic with reference to the indicate variant; 
    0/1: heterozygous with reference to the indicated variant; 0/0: Neither allele has the indicated variant.
    Citations in the xrefs column show previous publications that have presented similar findings."""
    table = create_suppl_table(list_of_rows=t_tests,
                               header=header,
                               header_format=header_field_formats,
                               caption=caption)
    print(table)


def create_hpo_onsets_table(hpo_onsets):
    header = ["cohort",  "genotype (A)",   "genotype (B)", "HPO Term", "p-val", "xrefs", ]
    header_field_formats = "lp{2cm}p{2cm}>{\\raggedright}p{4.5cm}lr" # need to import array package for raggedright
    caption = """Log rank tests performed using GPSEA to assess association between a genotype and the age of
    onset of a phenotypic feature represented by an HPO term. 1/1: biallelic with reference to the indicate variant; 
    0/1: heterozygous with reference to the indicated variant; 0/0: Neither allele has the indicated variant.
    Citations in the xrefs column show previous publications that have presented similar findings."""
    table = create_suppl_table(list_of_rows=hpo_onsets,
                               header=header,
                               header_format=header_field_formats,
                               caption=caption)
    print(table)


def create_disease_onsets_table(disease_onsets):
    header = ["cohort",  "genotype (A)",   "genotype (B)", "Disease", "p-val", "xrefs", ]
    header_field_formats = "lp{2cm}p{2cm}>{\\raggedright}p{4.5cm}lr" # need to import array package for raggedright
    caption = """Log rank tests performed using GPSEA to assess association between a genotype and the age of
    onset of a disease. 1/1: biallelic with reference to the indicate variant; 
    0/1: heterozygous with reference to the indicated variant; 0/0: Neither allele has the indicated variant.
    Citations in the xrefs column show previous publications that have presented similar findings."""
    table = create_suppl_table(list_of_rows=disease_onsets,
                               header=header,
                               header_format=header_field_formats,
                               caption=caption)
    print(table)


def create_mortality_table(mortality):
    header = ["cohort",  "genotype (A)",   "genotype (B)", "Disease", "p-val", "xrefs", ]
    header_field_formats = "lp{3.5cm}p{3.5cm}>{\\raggedright}p{3cm}lr" # need to import array package for raggedright
    caption = """Log rank tests performed using GPSEA to assess association between a genotype and the age of
    death of individuals with a disease. 1/1: biallelic with reference to the indicate variant; 
    0/1: heterozygous with reference to the indicated variant; 0/0: Neither allele has the indicated variant.
    Citations in the xrefs column show previous publications that have presented similar findings."""
    table = create_suppl_table(list_of_rows=mortality,
                               header=header,
                               header_format=header_field_formats,
                               caption=caption)
    print(table)

def create_phenotype_scores_table(phenotype_scores):
    header = ["cohort",  "genotype (A)",   "genotype (B)", "Scorer", "p-val", "xrefs", ]
    header_field_formats = "lp{3.5cm}p{3.5cm}>{\\raggedright}p{3cm}lr" # need to import array package for raggedright
    caption = """Mann-Whitney U tests performed using GPSEA to assess association between a genotype and the age of
    death of individuals with a disease. 1/1: biallelic with reference to the indicate variant; 
    0/1: heterozygous with reference to the indicated variant; 0/0: Neither allele has the indicated variant.
    Citations in the xrefs column show previous publications that have presented similar findings."""
    table = create_suppl_table(list_of_rows=phenotype_scores,
                               header=header,
                               header_format=header_field_formats,
                               caption=caption)
    print(table)


t_tests, hpo_onsets, disease_onsets,  mortality,  phenotype_scores = get_mono_test_results()
#print_summary_table(t_tests, hpo_onsets, disease_onsets,  mortality,  phenotype_scores)
#create_t_test_table(t_tests=t_tests)
#create_hpo_onsets_table(hpo_onsets=hpo_onsets)
#create_disease_onsets_table(disease_onsets=disease_onsets)
#create_mortality_table(mortality=mortality)
create_phenotype_scores_table(phenotype_scores=phenotype_scores)






