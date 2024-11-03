import pandas as pd
import numpy as np
import typing
import hpotk
from gpsea.model import Cohort
from gpsea.analysis.pcats import HpoTermAnalysisResult

from enum import Enum


class StatTest(Enum):
    FET = 1
    T_TEST = 2
    U_TEST = 3
    LOG_RANK = 4


def format_term_id(
        hpo: hpotk.MinimalOntology,
        term_id: hpotk.TermId,
) -> str:
    """
        Format a `term_id` as a `str`. HPO term ID is formatted as `<name> [<term_id>]` whereas other term IDs
        are formatted as CURIEs (e.g. `OMIM:123000`).
        """
    if term_id.prefix == "HP":
        term_name = hpo.get_term_name(term_id)
        return f"{term_name} [{term_id.value}]"
    else:
        return term_id.value


class SigTestResult:
    """
    Represent the salient data about a statistically significant Genotype Phenotype Correlation.
    """

    def __init__(self,
                 test_type: StatTest,
                 hpo_item: str,
                 geno_a: str,
                 geno_b: str,
                 with_geno_a: str,
                 with_geno_b: str,
                 p_value: float,
                 adj_p_value: float,
                 gt_pred: str
                 ):
        self._test_type = test_type
        self._hpo_item = hpo_item
        self._geno_a = geno_a
        self._geno_b = geno_b
        self._with_geno_a = with_geno_a
        self._with_geno_b = with_geno_b
        self._p_value = p_value
        self._adj_p_value = adj_p_value
        self._gt_pred = gt_pred

    @property
    def test_type(self) -> StatTest:
        return self._test_type

    def get_test_type(self) -> StatTest:
        if self._test_type == StatTest.FET:
            return "FET"
        elif self._test_type == StatTest.T_TEST:
            return "t test"
        elif self._test_type == StatTest.U_TEST:
            return "U test"
        elif self._test_type == StatTest.LOG_RANK:
            return "log rank"
        else:
            raise ValueError(f"Unknown test type {self._test_type}")

    @property
    def hpo_item(self) -> str:
        return self._hpo_item

    @property
    def geno_a(self) -> str:
        return self._geno_a

    @property
    def geno_b(self) -> str:
        return self._geno_b

    @property
    def with_geno_a(self) -> str:
        return self._with_geno_a

    @property
    def with_geno_b(self) -> str:
        return self._with_geno_b

    @property
    def p_value(self) -> float:
        return self._p_value

    @property
    def adj_p_value(self) -> float:
        return self._adj_p_value

    @property
    def gt_pred(self) -> str:
        return self._gt_pred


class SignificantResults:
    """
    Represents the significant results for any of the four tests (FET, t test, U test, log rank test) 
    """

    def __init__(self,
                 hpo: hpotk.MinimalOntology):
        self._hpo = hpo
        self._significant_results = list()
        self._n_usable_fet = None
        self._total_tested_fet = None

    def fisher_exact_test(self, result: HpoTermAnalysisResult) -> typing.List[SigTestResult]:
        """
        The result is typically a list of Fisher Exact Test results. We will keep only
        those with an adjusted p value of 0.05 or lower. We first add the values to
        a pandas dataframe for convenience, and transfer data for signficiant values
        to TestResult objects to simplify display.
        """
        assert isinstance(result, HpoTermAnalysisResult)
        fet_test_results = list()
        self._n_usable_fet = result.n_usable
        self._total_tested_fet = 42  # result.total_tested
        # Row index: a list of tested HPO terms
        pheno_idx = pd.Index(result.phenotypes)
        # Column index: multiindex of counts and percentages for all genotype predicate groups
        gt_idx = pd.MultiIndex.from_product(
            iterables=(result.gt_predicate.get_categories(), ("Count", "Percent")),
            names=(result.gt_predicate.get_question_base(), None),
        )

        # We'll fill this frame with data
        df = pd.DataFrame(index=pheno_idx, columns=gt_idx)

        for ph_predicate, count in zip(result.pheno_predicates, result.all_counts):
            # Sum across the phenotype categories (collapse the rows).
            gt_totals = count.sum()

            for gt_cat in count.columns:
                cnt = count.loc[ph_predicate.present_phenotype_category, gt_cat]
                total = gt_totals[gt_cat]
                df.loc[ph_predicate.phenotype, (gt_cat, "Genotype")] = str(gt_cat)
                df.loc[ph_predicate.phenotype, (gt_cat, "Count")] = f"{cnt}/{total}"
                pct = 0 if total == 0 else round(cnt * 100 / total)
                df.loc[ph_predicate.phenotype, (gt_cat, "Percent")] = f"{pct}%"

        # Add columns with p values and corrected p values (if present)
        p_val_col_name = "p values"
        corrected_p_val_col_name = "Corrected p values"
        if result.corrected_pvals is not None:
            df.insert(df.shape[1], ("", corrected_p_val_col_name), result.corrected_pvals)
        df.insert(df.shape[1], ("", p_val_col_name), result.pvals)

        # Format the index values: `HP:0001250` -> `Seizure [HP:0001250]` if the index members are HPO terms
        # or just use the term ID CURIE otherwise (e.g. `OMIM:123000`).
        labeled_idx = df.index.map(lambda term_id: format_term_id(self._hpo, term_id))

        # Last, sort by corrected p value or just p value
        df = df.set_index(labeled_idx)
        # and only report the tested HPO terms, remove rows with NaN in p value column
        with_p_value = df[("", p_val_col_name)].notna()
        if result.corrected_pvals is None:
            raise ValueError("corrected p vals are not set for FET")
        df = df.sort_values(by=[("", corrected_p_val_col_name), ("", p_val_col_name)]).loc[with_p_value.index]
        self._total_tested_fet = df.shape[0] # row count
        gt_pred = result.gt_predicate.display_question()
        for idx, row in df.iterrows():
            hpo_item = str(idx)
            geno_a = row.iloc[0]
            geno_a_ratio = row.iloc[1]
            geno_a_perc = row.iloc[2]
            geno_b = row.iloc[3]
            geno_b_ratio = row.iloc[4]
            geno_b_perc = row.iloc[5]
            p_val = row.iloc[6]
            adj_p_val = row.iloc[7]
            if p_val is None or np.isnan(p_val):
                continue
            elif p_val > 0.05:
                continue
            else:
                with_geno_a = f"{geno_a_ratio} ({geno_a_perc}%)"
                with_geno_b = f"{geno_b_ratio} ({geno_b_perc}%)"
                fet_tr = SigTestResult(test_type=StatTest.FET,
                                       hpo_item=hpo_item,
                                       geno_a = geno_a,
                                       geno_b = geno_b,
                                       with_geno_a = with_geno_a,
                                       with_geno_b = with_geno_b,
                                       p_value=p_val,
                                       adj_p_value=adj_p_val,
                                       gt_pred=gt_pred)
                fet_test_results.append(fet_tr)
        print(f"extracted {len(fet_test_results)} FET significant results")
        self._significant_results = fet_test_results


    @property
    def n_usable_fet(self):
        return self._n_usable_fet

    @property
    def total_tested_fet(self):
        return self._total_tested_fet

    @property
    def significant_results(self):
        return self._significant_results



class GpseaSummarizer:

    def __init__(self,
                 version: str,
                 caption: str,
                 cohort: Cohort,
                 tx_id: str,
                 hpo: hpotk.MinimalOntology,
                 sig_results: SignificantResults,
                 ):
        """
        :param version: version of GPSEA used for the analysis
        :type version: str
        :param caption: description that will be placed as the caption of the table of results
        :type caption: str
        :param cohort: object representing the collection of phenopackets
        :type cohort: Cohort
        :param tx_id: Transcript identifier used for the analysis
        :param hpo: Human Phenotype Ontology (object)
        :type hpo: hpotk.MinimalOntology

    
        """
        self._version = version
        self._caption = caption
        self._tx_id = tx_id
        self._cohort = cohort
        self._hpo = hpo
        self._fet_significant_results = list()
        self._significant_results = sig_results

    @property
    def gpsea_version(self):
        return self._version

    @property
    def tx_id(self):
        return self._tx_id

    @property
    def hpo_version(self):
        return self._hpo.version

    @property
    def significant_results(self):
        return self._significant_results
