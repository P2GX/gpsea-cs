from os.path import dirname, abspath, join
from csv import DictReader
import statistics
from scipy.stats import mannwhitneyu

THIS_DIR = dirname(abspath(__file__))
SUPPLEMENT_DIR = join(THIS_DIR, "v9_3")
COHORT_DASHBOARD = join(SUPPLEMENT_DIR, "cohort_dashboard.txt")
FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "fisher_exact_test_dashboard.txt")


def calculate_stats(numbers):
    if not numbers:
        raise ValueError("The input list must not be empty.")
    
    stats = {
        "mean": statistics.mean(numbers),
        "median": statistics.median(numbers),
        "min": min(numbers),
        "max": max(numbers),
        "standard_deviation": statistics.stdev(numbers) if len(numbers) > 1 else 0
    }
    return stats

def print_stats(title, numbers):
    stats = calculate_stats(numbers)
    print(f"{title} - mean {stats['mean']:.1f}, sd {stats['standard_deviation']:.1f}, median {stats['median']:.1f}, min {int(stats['min'])}, max {int(stats['max'])}  ")


def print_mf_stats(male_counts, female_counts, unknown_sex_list):
    n_male = sum(male_counts)
    n_female = sum(female_counts)
    n_unknown = sum(unknown_sex_list)
    total = n_male + n_female + n_unknown
    perc = 100.0 * (n_male + n_female)/total
    perc_male = 100.0 *  n_male / (n_male + n_female)
    perc_f = 100.0 *  n_female / (n_male + n_female)
    print(f"Total: {total} with {perc:.1f}% having data on sex of participants. Of these, {perc_male:.1f}% were male and {perc_f:.1f}% female")



def perform_mann_whitney_u(group1, group2, alternative='two-sided'):
    """
    Perform the Mann-Whitney U test to compare two groups.
    """
    result = mannwhitneyu(group1, group2, alternative=alternative)
    return {
        "U statistic": result.statistic,
        "p-value": result.pvalue
    }

def get_cohort_counts():
    individual_counts = list()
    total_individuals = 0
    male_counts = list()
    female_counts = list()
    unknown_sex_list = list()
    hpo_counts = list()
    cohorts_with_measurements = 0
    n_cohorts = 0
    n_cohorts_with_at_least_one_sig = set()
    n_sig_results = 0
    not_gene_set = {"Kabuki", "Robinow", "LDS 1 and 2", "LDS 1 and 3", "LDS 3 and 6"}
    genes = set()
    with open(COHORT_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            cohort = row["#cohort"]
            n_cohorts += 1
            individuals = int(row["individuals"])
            total_individuals += individuals
            males = int(row["males"])
            females = int(row["females"])
            unknown_s = int(row["n_unknown_sex"])
            total_hpo = int(row["total_hpo"])
            if int(row["total_measurements"]) > 0:
                cohorts_with_measurements += 1
            individual_counts.append(individuals)
            male_counts.append(males)
            female_counts.append(females)
            unknown_sex_list.append(unknown_s)
            hpo_counts.append(total_hpo)
            if cohort not in not_gene_set:
                genes.add(cohort)
    print_stats("individuals per cohort", individual_counts)
    print_stats("males", male_counts)
    print_stats("females", female_counts)
    print_stats("unknown sex", unknown_sex_list)
    print_mf_stats(male_counts, female_counts, unknown_sex_list)
    print_stats("HPO terms per cohort", hpo_counts)
    print(f"Cohorts: {n_cohorts}")
    print(f"Cohorts with measurements: {cohorts_with_measurements}")
    print(f"total genes: {len(genes)}")
    print(f"total_individuals:  {total_individuals}")


    
def fisher_test_stats():
    testable_list = list()
    tested_list = list()
    nsig_list = list()
    cohorts = set()
    significant_cohorts = set()
    with open(FISHER_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            cohort_name = row["#cohort_name"]
            total_hpo_testable = int(row["total_hpo_testable"])
            total_hpo_tested = int(row["total_hpo_tested"])
            nsig = int(row["nsig"])
            testable_list.append(total_hpo_testable)
            tested_list.append(total_hpo_tested)
            nsig_list.append(nsig)
            if nsig > 0:
                significant_cohorts.add(cohort_name)
            cohorts.add(cohort_name)
    print_stats("total_hpo_testable", testable_list)
    print_stats("total_hpo_tested", tested_list)
    print(f"nsig {len(nsig_list)}")
    print(f"Total cohorts with at least one significant result: {len(significant_cohorts)}; total cohortts tested: {len(cohorts)}")
    

def compare_group_sizes():
    """
    Compare cohort size of cohorts with a significant result and those without a significant results
    """
    with open(FISHER_DASHBOARD) as file:
        all_tested_cohorts = set()
        cohorts_with_sig_results = set()
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            cohort_name = row["#cohort_name"]
            nsig = int(row["nsig"])
            all_tested_cohorts.add(cohort_name)
            if nsig > 0:
                cohorts_with_sig_results.add(cohort_name)
    cohorts_without_sig_results = all_tested_cohorts.difference(cohorts_with_sig_results)
    with_sig = dict()
    without_sig = dict()
    with open(COHORT_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            cohort_name = row["#cohort"]
            individuals = int(row["individuals"])
            if cohort_name in cohorts_with_sig_results:
                with_sig[cohort_name] = individuals
            elif cohort_name in cohorts_without_sig_results:
                without_sig[cohort_name] = individuals
    with_sig = list(with_sig.values())
    without_sig = list(without_sig.values())


    result = mannwhitneyu(with_sig, without_sig, alternative="greater")
    print(f"Mann Whitney U-test: U-statistic {result.statistic}; p-value: {result.pvalue:.7f} (One sided; null hyp: cohorts with no sig values are not smaller than those with sig values)")

def get_unique_disease_counts():
    """
    Count the total number of diseases analyzed in this project.
    Note that this is not the same as the number of genes because some of the tested genes are associated with more
    than one Mendelian disease.
    """
    diseases = set()
    with open(COHORT_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            n_diseases = int(row["n_diseases"])
            disease_string = row["disease_string"]
            disease_list = disease_string.split(";")
            if len(disease_list) != n_diseases:
                raise ValueError("Misparsed disease string" + disease_string)
            for disease in disease_list:
                diseases.add(disease.strip())
    print(f"Total unique diseases tested {len(diseases)}")
    #for d in diseases:
    #    print(d)

get_cohort_counts()
fisher_test_stats()
compare_group_sizes()
get_unique_disease_counts()