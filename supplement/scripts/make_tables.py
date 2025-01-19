from os.path import dirname, abspath, join
from csv import DictReader
from mylatextable import MyLongTable


SUPPLEMENT_DIR = dirname(dirname(abspath(__file__)))
SIG_FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "sig_fisher_exact_test_dashboard.txt")
MEASUREMENT_DASHBOARD = join(SUPPLEMENT_DIR, "measurement_dashboard.txt")


def format_p(number) -> str: 
    n = float(number)
    if n > 0.001:
        return number
    else:
        return f"${n:.1e}".replace('e', r' \times 10^{') + "}$"



def get_sig_fisher_table():
    """
    #cohort_name	a_genotype	b_genotype	nsig	n_tests_performed	hpo_item	with_geno_a	with_geno_b	pval	adj_pval
    """
    header = ["cohort",  "HPO", "genotype (A)",  "Counts (A)", "genotype (B)", "Counts (B)", "p-val", "adj. p-val"]
    header_field_formats = "p{1cm}p{3cm}p{2cm}p{1cm}p{2cm}p{1cm}p{1.5cm}p{1.5cm}"
    table = MyLongTable(header_fields=header, use_booktabs=True, header_field_formats=header_field_formats, fontsize="scriptsize")
    with open (SIG_FISHER_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            with_geno_a = row["with_geno_a"].replace("%", "\\%")
            with_geno_b = row["with_geno_b"].replace("%", "\\%")
            items = [row["#cohort_name"],  row["hpo_item"],  row["a_genotype"], with_geno_a, 
                     row["b_genotype"], with_geno_b, format_p(row["pval"]), format_p(row["adj_pval"])]
            table.add_row(items)
    return table.get_latex()





table = get_sig_fisher_table()
fh = open("significant_fisher_results.tex", "wt")
fh.write(table)
fh.close()