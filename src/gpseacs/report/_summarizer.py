import typing
from typing import TextIO
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
NOTEBOOK_SUMMARY_DIR = os.path.abspath(os.path.join(current_dir, "../../../supplement"))
COHORT_SUMMARY_FILE = os.path.join(NOTEBOOK_SUMMARY_DIR, "cohort_dashbord.csv")
COHORT_SUMMARY_HEADER = ["cohort", 'individuals', 'females', 'males', 'total_hpo', 'total_measurements']


class NotebookDashboard:
    """
    We want to write a summary of the reuslts from each notebook into a single centralized file that
    we can use to create tables with. There are several kinds of information
    1. Cohort (number of individuals, male, female, mortality etc
    2. FET results: summary of each test performed Total number of tests performed,
        total significant, summary of each significant result. Note that we will create two output tables for the
        FET - One on the level of the test procedure, and another table for each significant result.
    3. MONO results: summary of each test performed and result

    There are thus four central files that we want to keep track of.
    Note that these files are only used to gather the results from the GPSEA-CS notebooks. Other users of GPSEA do NOT
    need any of this functionality, which we developed to streamline the creation of the supplemental files.
    """

    def __init__(self):
        if not os.path.exists(COHORT_SUMMARY_FILE):
            with open(COHORT_SUMMARY_FILE, "wt") as file:
                file.write("")  # Create file only if needed
        self._cohort_name = None
        self._hpo_version = None
        self._gpsea_version = None

    def _get_or_panic(self, context: typing.Mapping[str, typing.Any], item: str) -> str:
        if not item in context:
            raise ValueError(f"Could not find {item} in context")
        else:
            return context.get(item)

    def update_dashboard(self, context: typing.Mapping[str, typing.Any]):
        """
        The strategy is to read in data for results from each of the potential components of the current cohort.
        If there are new results to report, we will update the corresponding files.
        Otherwise we do not change the file.
        This function ingests some general data we will put in all of the dashboard files and calls functions for each
        of the four files (components)
        """
        self._cohort_name = self._get_or_panic(context, "cohort_name")
        self._hpo_version = self._get_or_panic(context, "hpo_version")
        self._gpsea_version = self._get_or_panic(context, "gpsea_version")
        self._n_female = self._get_or_panic(context, "n_female")
        self._n_male =  self._get_or_panic(context, "n_male")
        self._n_unknown = self._get_or_panic(context, "n_unknown_sex")   
        self._n_total_individual_count = self._get_or_panic(context, "n_total_individual_count")  
        self._n_hpo_terms = self._get_or_panic(context, "n_hpo_terms")  
        self._n_alive = self._get_or_panic(context, "n_alive")  
        self._n_deceased = self._get_or_panic(context, "n_deceased")  
        self._n_unknown_vital = self._get_or_panic(context, "n_unknown_vital")  
        self._n_with_age_of_last_encounter = self._get_or_panic(context, "n_with_age_of_last_encounter")  
        self._n_with_onset = self._get_or_panic(context, "n_with_onset")  
        self._n_measurements = self._get_or_panic(context, "n_measurements")  
        self._n_diseases = self._get_or_panic(context, "n_diseases")
        self._disease_string = self._get_or_panic(context, "disease_string")   
        self._update_cohort_summary(context)

    def _get_cohort_name_to_line_d(self, filename: str):
        cohort_to_line_d = dict()
        with open(filename) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                else:
                    fields = line.strip().split("\t")
                    cohort_name = fields[0]
                    cohort_to_line_d[cohort_name] = line
        return cohort_to_line_d
    
    def _get_header_line(self, header) -> str:
        """
        standard format for header of our dashboard files
        """
        return '#' + "\t".join + '˜n'

    def _sort_and_print_dictionary(self, cohort_name_to_line_d, fname, header) -> None:
        sorted_dict = dict(sorted(cohort_name_to_line_d.items()))
        header_line = "\t".join(header) 
        with open(fname, "wt") as file:
            file.write(header_line + "\n")
            for k, v in sorted_dict.items():
                file.write(v + "\n")

    
    
    
    def _update_cohort_summary(self, context: typing.Mapping[str, typing.Any]):
        cohort_to_line_d = self._get_cohort_name_to_line_d(COHORT_SUMMARY_FILE)
        items = [self._cohort_name, self._n_total_individual_count, self._n_female, self._n_male, self._n_hpo_terms, self._n_measurements]
        items = [str(x) for x in items]
        cohort_to_line_d[self._cohort_name] = "\t".join(items)
        self._sort_and_print_dictionary(cohort_to_line_d, COHORT_SUMMARY_FILE, COHORT_SUMMARY_HEADER)




