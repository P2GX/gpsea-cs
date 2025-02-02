from os.path import dirname, abspath, join
from csv import DictReader

# Create tables to summarize the total number of tests performed and their results
# We will use this for Table 1.

THIS_DIR = dirname(abspath(__file__))
SUPPLEMENT_DIR = join(THIS_DIR, "v9_3")
COHORT_DASHBOARD = join(SUPPLEMENT_DIR, "cohort_dashboard.txt")
FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "fisher_exact_test_dashboard.txt")
MEASUREMENT_DASHBOARD = join(SUPPLEMENT_DIR, "measurement_dashboard.txt")
SIG_FISHER_DASHBOARD = join(SUPPLEMENT_DIR, "sig_fisher_exact_test_dashboard.txt")


def get_measurements():
     with open(MEASUREMENT_DASHBOARD) as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            print(row)




get_measurements()