import typing
from collections import defaultdict
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
NOTEBOOK_SUMMARY_DIR = os.path.abspath(os.path.join(current_dir, "../../../supplement"))
COHORT_SUMMARY_FILE = os.path.join(NOTEBOOK_SUMMARY_DIR, "cohort_dashboard.txt")
COHORT_SUMMARY_HEADER = ["cohort", 'individuals', 'females', 'males',"n_unknown_sex", 'total_hpo', 'total_measurements',
                          "hpo_version", "gpsea_version",  "n_total_individual_count", "n_alive", 
                          "n_deceased", "n_unknown_vital", "n_with_age_of_last_encounter", "n_with_onset", "n_diseases", "disease_string"] 
MEASUREMENT_SUMMARY_FILE = os.path.join(NOTEBOOK_SUMMARY_DIR, "measurement_dashboard.txt")
MEASUREMENT_SUMMARY_HEADER = ["cohort", "a_genotype", "b_genotype", "name", "test_name", "description", "variable_name", "pval", "interpretation", "xrefs"]
FISHER_SUMMARY_FILE = os.path.join(NOTEBOOK_SUMMARY_DIR, "fisher_exact_test_dashboard.txt")
FISHER_SUMMARY_HEADER = ["cohort_name","total_hpo_testable", "total_hpo_tested", "a_genotype", "b_genotype", "nsig"]
SIG_FISHER_SUMMARY_FILE = os.path.join(NOTEBOOK_SUMMARY_DIR, "sig_fisher_exact_test_dashboard.txt")
SIG_FISHER_SUMMARY_HEADER = ["cohort_name", "a_genotype", "b_genotype", "nsig", "n_tests_performed", "hpo_item", "with_geno_a", "with_geno_b", "pval", "adj_pval"]

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
        self._touch_file_if_needed(COHORT_SUMMARY_FILE)
        self._touch_file_if_needed(MEASUREMENT_SUMMARY_FILE)
        self._touch_file_if_needed(FISHER_SUMMARY_FILE)
        self._touch_file_if_needed(SIG_FISHER_SUMMARY_FILE)
        self._cohort_name = None
        self._hpo_version = None
        self._gpsea_version = None

    def _touch_file_if_needed(self, fname) -> None:
        # create file the first time we use it
        if not os.path.exists(fname):
            with open(fname, "wt") as file:
                file.write("")  # Create file only if needed

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
        self._update_monotest_summary(context)
        self._update_fisher_test_summary(context)
        self._update_sig_fisher_test_summary(context)

    def _get_cohort_name_to_line_d(self, filename: str):
        cohort_to_line_d = dict()
        with open(filename) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                else:
                    fields = line.strip().split("\t")
                    cohort_name = fields[0]
                    cohort_to_line_d[cohort_name] = line.strip()
        return cohort_to_line_d
    
    def _get_cohort_name_to_line_set_d(self, filename: str):
        cohort_to_line_d = defaultdict(set)
        with open(filename) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                else:
                    fields = line.strip().split("\t")
                    cohort_name = fields[0]
                    cohort_to_line_d[cohort_name].add(line.strip())
        return cohort_to_line_d
    

    def _sort_and_print_dictionary(self, cohort_name_to_line_d, fname, header) -> None:
        sorted_dict = dict(sorted(cohort_name_to_line_d.items()))
        header_line =  '#' + "\t".join(header) 
        with open(fname, "wt") as file:
            file.write(header_line + "\n")
            for _, v in sorted_dict.items():
                file.write(v + "\n")

    def _sort_and_print_set_dictionary(self, cohort_name_to_line_d, fname, header) -> None:
        sorted_dict = dict(sorted(cohort_name_to_line_d.items()))
        header_line =  '#' + "\t".join(header) 
        with open(fname, "wt") as file:
            file.write(header_line + "\n")
            for k, v in sorted_dict.items():
                for line in v:
                    file.write(line + "\n")
    
    
    def _update_cohort_summary(self, context: typing.Mapping[str, typing.Any]):
        cohort_to_line_d = self._get_cohort_name_to_line_d(COHORT_SUMMARY_FILE)
        items = [self._cohort_name, self._n_total_individual_count, self._n_female, self._n_male, self._n_unknown, self._n_hpo_terms, self._n_measurements,
                 self._hpo_version, self._gpsea_version, self._n_total_individual_count ,
                   self._n_alive, self._n_deceased, self._n_unknown_vital, self._n_with_age_of_last_encounter,
                   self._n_with_onset,  self._n_diseases, self._disease_string] 
        items = [str(x) for x in items]
        cohort_to_line_d[self._cohort_name] = "\t".join(items)
        self._sort_and_print_dictionary(cohort_to_line_d, COHORT_SUMMARY_FILE, COHORT_SUMMARY_HEADER)
        
    def _update_monotest_summary(self, context: typing.Mapping[str, typing.Any]):
        """
        {'a_genotype':  'b_genotype':, 'name':  'test_name': 'description': , 'variable_name': 'LOINC:74892-1', 
        'pval': '6.12e-10', 'interpretation': , 'xrefs': '\\cite{PMID_33580884}'}
        """
        cohort_to_line_d = self._get_cohort_name_to_line_set_d(MEASUREMENT_SUMMARY_FILE)
        mono_result_list = context.get("mono_result_list")
        if len(mono_result_list) == 0:
            return
        for mr in mono_result_list:
            a_genotype = self._get_or_panic(mr, "a_genotype")
            b_genotype = self._get_or_panic(mr, "b_genotype")
            name = self._get_or_panic(mr, "name")
            test_name = self._get_or_panic(mr, "test_name")
            description = self._get_or_panic(mr, "description")
            variable_name = self._get_or_panic(mr, "variable_name")
            pval = self._get_or_panic(mr, "pval")
            interpretation = self._get_or_panic(mr, "interpretation")
            xrefs =  self._get_or_panic(mr, "xrefs")
            xrefs = self._unformat_xrefs(xrefs)
            items = [self._cohort_name, a_genotype, b_genotype, name, test_name, description, variable_name, pval, interpretation, xrefs ]
            items = [str(x) for x in items]
            cohort_to_line_d[self._cohort_name].add("\t".join(items))
        self._sort_and_print_set_dictionary(cohort_to_line_d, MEASUREMENT_SUMMARY_FILE, MEASUREMENT_SUMMARY_HEADER)

    def _unformat_xrefs(self, xref) -> str:
        output_text = xref.replace(r"\cite{", "").replace("}", "").replace("_", ":")
        return output_text.strip()
    
    def _update_fisher_test_summary(self, context: typing.Mapping[str, typing.Any]):
        """
        This gets a summary of all FET tests we performed, regardless of whether they were significant or not
        """
        all_fet_results = self._get_or_panic(context, "all_fet_results")
        if len(all_fet_results) == 0 :
            return
        cohort_to_line_d = self._get_cohort_name_to_line_set_d(FISHER_SUMMARY_FILE)
        cohort_name = self._get_or_panic(context, "cohort_name")
        for fet in all_fet_results:
            total_hpo_testable =  self._get_or_panic(fet, "total_hpo_testable")
            total_hpo_tested =  self._get_or_panic(fet, "total_hpo_tested")
            a_genotype =  self._get_or_panic(fet, "a_genotype")
            b_genotype=  self._get_or_panic(fet, "b_genotype")
            nsig =  self._get_or_panic(fet, "nsig")
            items = [cohort_name, total_hpo_testable, total_hpo_tested, a_genotype, b_genotype, nsig]
            items = [str(x) for x in items]
            cohort_to_line_d[cohort_name].add("\t".join(items))
        self._sort_and_print_set_dictionary(cohort_to_line_d, FISHER_SUMMARY_FILE, FISHER_SUMMARY_HEADER)

    def _update_sig_fisher_test_summary(self, context: typing.Mapping[str, typing.Any]):
        """
        This gets a summary of statistically significant FET tests we performed.
        This is called sig_fet_results in the context
        {'general_info': {'a_genotype':, 'b_genotype': , 'n_sig_results': 0, 'n_tests_performed': 25}, 'sig_result_list': []}
{'general_info': {'a_genotype': 'OMIM:606693', 'b_genotype': 'OMIM:617225', 'n_sig_results': 2, 'n_tests_performed': 24}, 'sig_result_list': 
[{'hpo_item': 'Bradykinesia [HP:0002067]', 'with_geno_a': '29/31 (94%)', 'with_geno_b': '4/10 (40%)', 'pval': '0.001', 'adj_pval': '0.013'}, {'hpo_item': 'Parkinsonism [HP:0001300]', 'with_geno_a': '27/27 (100%)', 'with_geno_b': '3/11 (27%)', 'pval': '3.37e-06', 'adj_pval': '8.10e-05'}]}

        """
        sig_fet_results = self._get_or_panic(context, "fet_result_list")
        if len(sig_fet_results) == 0 :
            return
        cohort_to_line_d = self._get_cohort_name_to_line_set_d(SIG_FISHER_SUMMARY_FILE)
        cohort_name = self._get_or_panic(context, "cohort_name")
        for fet in sig_fet_results:
            general_info = fet.get("general_info")
            a_genotype =  self._get_or_panic(general_info, "a_genotype")
            b_genotype=  self._get_or_panic(general_info, "b_genotype")
            nsig =  self._get_or_panic(general_info, "n_sig_results")
            n_tests_performed = self._get_or_panic(general_info, "n_tests_performed")
            sig_result_list = fet.get("sig_result_list")
            for d in sig_result_list:
                hpo_item = self._get_or_panic(d, "hpo_item")
                with_geno_a = self._get_or_panic(d, "with_geno_a")
                with_geno_b = self._get_or_panic(d, "with_geno_b")
                pval = self._get_or_panic(d, "pval")
                adj_pval = self._get_or_panic(d, "adj_pval")
                items = [cohort_name, a_genotype, b_genotype, nsig, n_tests_performed, hpo_item, with_geno_a, with_geno_b, pval, adj_pval]
                items = [str(x) for x in items]
                cohort_to_line_d[cohort_name].add("\t".join(items))
        self._sort_and_print_set_dictionary(cohort_to_line_d, SIG_FISHER_SUMMARY_FILE, SIG_FISHER_SUMMARY_HEADER)
