from os.path import dirname, abspath, join
from csv import DictReader
import statistics
from scipy.stats import mannwhitneyu


SUPPLEMENT_DIR = dirname(dirname(abspath(__file__)))
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
    male_counts = list()
    female_counts = list()
    unknown_sex_list = list()
    hpo_counts = list()
    cohorts_with_measurements = 0
    n_cohorts = 0

    with open(COHORT_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            individuals = int(row["individuals"])
            males = int(row["males"])
            females = int(row["females"])
            unknown_s = int(row["n_unknown_sex"])
            total_hpo = int(row["total_hpo"])
            has_hpo = int(row["total_measurements"]) > 0
            individual_counts.append(individuals)
            male_counts.append(males)
            female_counts.append(females)
            unknown_sex_list.append(unknown_s)
            hpo_counts.append(total_hpo)
            n_cohorts += 1
            if has_hpo:
                cohorts_with_measurements += 1
    print_stats("individuals", individual_counts)
    print_stats("males", male_counts)
    print_stats("females", female_counts)
    print_stats("unknown sex", unknown_sex_list)
    print_stats("HPO terms per cohort", hpo_counts)
    print(f"Cohorts: {n_cohorts}")
    print(f"Cohorts with measurements: {cohorts_with_measurements}")

    
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
    print_stats("nsig", nsig_list)
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

    print_stats("cohorts_with_sig_results", with_sig)
    print_stats("cohorts_without_sig_results", without_sig)
    result = mannwhitneyu(with_sig, without_sig, alternative="greater")
    print(f"Mann Whitney U-test: U-statistic {result.statistic}; p-value: {result.pvalue:.7f} (One sided; null hyp: cohorts with no sig values are not smaller than those with sig values)")




get_cohort_counts()
fisher_test_stats()
compare_group_sizes()