import abc
import io
import re
import typing
import pandas as pd
import numpy as np
import matplotlib

import hpotk
from gpsea.model import Cohort
from gpsea.analysis import MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult
from .util import open_text_io_handle_for_writing, process_latex_template
from jinja2 import Environment, PackageLoader

from ._summarizer import NotebookDashboard


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

def format_p_value(p_value):
    if p_value > 0.001:
        # Format with 3 significant digits
        return f"{p_value:.3f}"
    else:
        # Format in scientific notation with 2 significant digits
        return f"{p_value:.2e}"

class GPAnalysisResultSummary:
    """
    This class represents a section of the summary report. For FET, it may contain zero to N rows. For the
    "mono" tests, it may contain zero to 1 row. The intention to create a section in a multipanel LaTeX figure
    for the supplemental material, and also to present a summary in the Jupyter notebook.
    """

    @staticmethod
    def from_mono(
            result: MonoPhenotypeAnalysisResult,
            xrefs: typing.Optional[typing.Iterable[str]] = None,
            interpretation: typing.Optional[str] = None,
    ) -> "GPAnalysisResultSummary":
        if xrefs is None:
            xrefs = {}
        else:
            xrefs = {
                result.phenotype.variable_name: xrefs,
            }
        return GPAnalysisResultSummary(
            result=result,
            xrefs=xrefs,
            interpretation=interpretation,
        )

    @staticmethod
    def from_multi(
            result: MultiPhenotypeAnalysisResult,
            xrefs: typing.Optional[typing.Mapping[str, typing.Iterable[str]]] = None,
            interpretation: typing.Optional[str] = None,
    ) -> "GPAnalysisResultSummary":
        if xrefs is None:
            xrefs = {}

        return GPAnalysisResultSummary(
            result=result,
            xrefs=xrefs,
            interpretation=interpretation,
        )

    def __init__(
            self,
            result: typing.Union[MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult],
            xrefs: typing.Mapping[str, typing.Collection[str]],
            interpretation: typing.Optional[str],
    ):
        assert isinstance(result, (MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult))
        self._result = result
        self._xrefs = dict(xrefs)
        self._interpretation = interpretation

    @property
    def result(self) -> typing.Union[MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult]:
        return self._result

    @property
    def xrefs(self) -> typing.Mapping[str, typing.Collection[str]]:
        """
        Get mapping from phenotype `variable_name` to G/P association cross references.
        """
        return self._xrefs

    @property
    def interpretation(self) -> typing.Optional[str]:
        return self._interpretation

    def is_mono(self) -> bool:
        return isinstance(self._result, MonoPhenotypeAnalysisResult)

    def is_multi(self) -> bool:
        return isinstance(self._result, MultiPhenotypeAnalysisResult)


class GpseaAnalysisReport:

    def __init__(
            self,
            name: str,
            cohort: Cohort,
            gene_symbol:str,
            mane_tx_id:str,
            mane_protein_id: str,
            caption: str="to do.",
            fet_results: typing.Iterable[GPAnalysisResultSummary]=None,
            mono_results: typing.Iterable[GPAnalysisResultSummary]=None,
    ):
        self._name = name
        self._cohort = cohort
        if fet_results is not None:
            for i, r in enumerate(fet_results):
                assert isinstance(r, GPAnalysisResultSummary), f"#{i} must be `GPAnalysisResultSummary` but was {type(r)}"
            self._fet_results = tuple(fet_results)
        else:
            self._fet_results = None
        if mono_results is not None:
            for i, r in enumerate(mono_results):
                assert isinstance(r, GPAnalysisResultSummary), f"#{i} must be `GPAnalysisResultSummary` but was {type(r)}"
            self._mono_results = tuple(mono_results)
        else:
            self._mono_results = None
        self._caption = generate_cohort_summary(cohort=cohort, caption=caption)
        n_variants = len(set(cohort.all_variants()))
        gene_caption, latex_gene_caption = generate_gene_summary(n_variants=n_variants, gene_symbol=gene_symbol, mane_tx_id=mane_tx_id, mane_protein_id=mane_protein_id)
        self._gene_caption = gene_caption
        self._latex_gene_caption = latex_gene_caption
        self._n_variants = len(set(cohort.all_variants()))
        self._disease_id_to_name = {d.identifier.value: d.name for d in cohort.all_diseases()}
        self._n_female = cohort.count_females()
        self._n_male = cohort.count_males()
        self._n_unknown = cohort.count_unknown_sex()
        self._n_individuals = cohort.total_patient_count 
        self._n_total_individual_count = cohort.total_patient_count
        self._n_hpo_terms = cohort.count_distinct_hpo_terms()
        self._n_alive = cohort.count_alive()
        self._n_deceased = cohort.count_deceased()
        self._n_with_age_of_last_encounter = cohort.count_with_age_of_last_encounter()
        self._n_with_onset = cohort.count_with_disease_onset()
        self._n_unknown_vital = cohort.count_unknown_vital_status()
        self._n_measurements = len(cohort.list_measurements())
        self._n_diseases = cohort.count_distinct_diseases()
        self._disease_string = "; ".join(self._disease_id_to_name.values())
        self._disease_id_string = "; ".join(self._disease_id_to_name.keys())



    @property
    def name(self) -> str:
        return self._name

    @property
    def cohort(self) -> Cohort:
        return self._cohort

    @property
    def fet_results(self) -> typing.Collection[GPAnalysisResultSummary]:
        return self._fet_results

    @property
    def mono_results(self) -> typing.Collection[GPAnalysisResultSummary]:
        return self._mono_results

    @property
    def caption(self) -> str:
        return self._caption

    @property
    def gene_caption(self) -> str:
        return self._gene_caption
    
    @property
    def n_variants(self) -> str:
        return self._n_variants
    
    @property
    def n_female(self) -> str:
        return self._n_female

    @property
    def n_male(self) -> str:
        return self._n_male

    @property
    def n_unknown_sex(self) -> str:
        return self._n_unknown 

    @property
    def n_total_individual_count(self) -> str:
        return self._n_total_individual_count
    @property
    def disease_id_to_name(self) -> str:
        return self._disease_id_to_name

    @property
    def n_diseases(self):
        return self._n_diseases

    @property
    def disease_string(self):
        return self._disease_string

    @property
    def n_individuals(self):
        return self._cohort.total_patient_count

    @property
    def n_total_individual_count(self):
        return self._cohort.total_patient_count

    @property
    def n_hpo_terms(self):
        return self._cohort.count_distinct_hpo_terms()

    @property
    def n_alive(self):
        return self._cohort.count_alive()

    @property
    def n_deceased(self):
        return self._cohort.count_deceased()

    @property
    def n_with_age_of_last_encounter(self):
        return self._cohort.count_with_age_of_last_encounter()

    @property
    def n_with_onset(self):
        return self._cohort.count_with_disease_onset()

    @property
    def n_unknown_vital(self):
        return self._cohort.count_unknown_vital_status()

    @property
    def n_measurements(self):
        return len(self._cohort.all_measurements())

    @property
    def latex_gene_caption(self) -> str:
        return self._latex_gene_caption

def generate_cohort_summary(cohort,
                            caption:str= None) -> str:
    disease_id_to_name = dict()
    for d in cohort.all_diseases():
        disease_id_to_name[d.identifier.value] = d.name
    n_female = cohort.count_females()
    n_male = cohort.count_males()
    n_unknown = cohort.count_unknown_sex()
    n_total = cohort.total_patient_count
    if n_unknown > 0:
        counts = f"The cohort comprised {n_total} individuals ({n_female} females, {n_male} males, {n_unknown} with unknown sex)."
    else:
        counts = f"The cohort comprised {n_total} individuals ({n_female} females, {n_male} males)."
    deceased = cohort.count_deceased()
    if deceased > 0:
        counts = f"{counts} {deceased} of these individuals were reported to be deceased."


    n_hpo_terms = f"A total of {cohort.count_distinct_hpo_terms()} HPO terms were used to annotate the cohort."
    disease_count_list = cohort.list_all_diseases()
    if len(disease_count_list) == 1:
        d = disease_count_list[0]
        d_id = d[0]
        d_count = d[1]
        d_name = disease_id_to_name.get(d_id)
        if d_count != n_total:
            raise ValueError(f"Expecting {n_total} individuals with {d_name} but got {d_count}")
        disease_desc = f"Disease diagnosis: {d_name} ({d_id})."
    elif len(disease_count_list) > 1:
        disease_items= []
        for d in disease_count_list:
            d_id = d[0]
            d_count = d[1]
            d_name = disease_id_to_name.get(d_id)
            disease_items.append(f"{d_name} ({d_id}) ({d_count} individuals)")
        disease_desc = f"Disease diagnoses: {', '.join(disease_items)}."


    if caption is None:
        return f"{counts} {n_hpo_terms} {disease_desc}"
    else:
        return f"{counts} {n_hpo_terms} {disease_desc} {caption}"


def generate_gene_summary(n_variants:int,
                          gene_symbol: str=None,
                          mane_tx_id: str=None,
                          mane_protein_id: str=None
                          ):
    if gene_symbol is None:
        variants = f"A total of {n_variants} unique variant alleles were found."
        latex_variants = f"A total of {n_variants} unique variant alleles were found."
    else:
        ltx = mane_tx_id.replace("_", "\\_")
        lprot = mane_protein_id.replace("_", "\\_")
        variants = f"A total of {n_variants} unique variant alleles were found in {gene_symbol} (transcript: {mane_tx_id}, protein id: {mane_protein_id})."
        latex_variants = f"A total of {n_variants} unique variant alleles were found in \\textit{{{gene_symbol}}} (transcript: \\texttt{{{ltx}}}, protein id: \\texttt{{{lprot}}})."
    return variants, latex_variants

class GpseaReportSummarizer(metaclass=abc.ABCMeta):
    """
    `GpseaReportSummarizer` formats the :class:`~gpseacs.report.GPSEAAnalysisResult`
    into an output in a chosen format (LaTeX, HTML, whatever...).
    """

    @abc.abstractmethod
    def summarize_report(
            self,
            caption: str,
            cohort: Cohort,
            result_list: typing.List[GpseaAnalysisReport],

    ):
        pass


class GpseaReport(metaclass=abc.ABCMeta):
    """
    `GpseaReport` summarizes an aspect of an analysis for the user.
    """

    @abc.abstractmethod
    def write(self, fh: typing.Union[io.IOBase, str]):
        """
        Write the report into the provided path or file handle.

        :param fh: a `str` with path
          or :class:`io.IOBase` with the file-like object for writing the report into.
        """
        pass

class HtmlGpseaNotebookSummarizer(GpseaReport):
    """
    A report where the content is formatted as a HTML `str`.
    """

    # NOT PART OF THE PUBLIC API

    def __init__(
            self,
            html: str,
    ):
        assert isinstance(html, str)
        self._html = html

    @property
    def html(self) -> str:
        """
        Get a `str` with the HTML report.
        """
        return self._html

    def write(self, fh: typing.Union[io.IOBase, str]):
        should_close = isinstance(fh, str)
        fout = None
        try:
            fout = open_text_io_handle_for_writing(fh)
            fout.write(self._html)
        except Exception:
            if should_close and fout is not None:
                fout.close()

    def _repr_html_(self) -> str:
        return self._html




class GpseaNotebookSummarizer(GpseaReportSummarizer):

    def __init__(self, hpo: hpotk.MinimalOntology, gpsea_version:str, alpha=0.05):
        self._hpo = hpo
        self._gpsea_version = gpsea_version
        self._alpha = alpha
        environment = Environment(loader=PackageLoader('gpseacs', 'templates'))
        self._cohort_template = environment.get_template("summary.html")


    def summarize_report(
            self,
            report: GpseaAnalysisReport
    ):
        context = self._prepare_context(report)
        html = self._cohort_template.render(context)
        dashboardSummarizer = NotebookDashboard()
        dashboardSummarizer.update_dashboard(context)
        return HtmlGpseaNotebookSummarizer(html)


    def mono_test(self, gps: GpseaReportSummarizer) -> typing.Dict[str, object]:
        """
        This function prepares the results of the t test, the U test, and the log rank test for output.
        """
        result: MonoPhenotypeAnalysisResult = gps.result
        test_result = dict()
        test_result["a_genotype"] = result.gt_clf.class_labels[0]
        test_result["b_genotype"] = result.gt_clf.class_labels[1]

        ptype = result.phenotype
        test_result["name"] = ptype.name
        if ptype.name == "Measurement":
            test_result["test_name"] = "t-test"
        else:
            test_result["test_name"] = ptype.name
        descr = ptype.description
        ## convert: Compute time until onset of Chronic mucocutaneous candidiasis
        descr = descr.replace("Compute time until onset of ", "Survival analysis: ")
        test_result["description"] = descr
        test_result["variable_name"] = ptype.variable_name
        test_result["pval"] = format_p_value(result.pval)
        test_result["xrefs"] =self.get_mono_xref_string(gps.xrefs)
        interp = gps.interpretation or ""
        test_result["interpretation"] = interp.replace("%", "\\%")

        return test_result


    def get_hpo_xref_string(self, xref_d, hpo_item):
        """
        For some of the tests, we supply a PMID. We will use this to create citations in the supplement with bibtex
        \\cite{PMID_1234},\\cite{PMID_7654}. If we cannot find anything, return "-". We will need to manually add the
        corresponding items to a bibtex "bib" file.
        """
        hpo_regex = r"HP:\d{7}"
        match = re.search(hpo_regex, hpo_item)
        if match:
            hpo_id = match.group()
            if hpo_id in xref_d:
                pmid_tuple = xref_d.get(hpo_id)
                items = ["\\cite{" + p.replace("PMID:", "PMID_") + "}"  for p in pmid_tuple]
                return ",".join(items)
        return "-" # couldn't find anything

    def get_mono_xref_string(self, xref_d) -> str:
        """
        For the "mono" tests, there is only one test, thus we do not need to match the key
        """
        xrefs = list()
        for pmid_tuple in xref_d.values():
            items = ["\\cite{" + p.replace("PMID:", "PMID_") + "}"  for p in pmid_tuple.values()]
            if len(items) > 0:
                xrefs.extend(items)
        if len(xrefs) > 0:
            return ",".join(xrefs)
        else:
            return "-"

    def fisher_exact_test(self, gps: GpseaReportSummarizer) -> typing.Tuple[typing.Dict[str,object], typing.List[typing.Dict[str, object]]]:
        """
        The result is typically a list of Fisher Exact Test results. We will keep only
        those with an adjusted p value of 0.05 or lower. We first add the values to
        a pandas dataframe for convenience, and transfer data for signficiant values
        to TestResult objects to simplify display.
        """
        result = gps.result
        xrefs = gps.xrefs
        assert isinstance(result, MultiPhenotypeAnalysisResult)
        general_info = dict()
        sig_result_list = list()
        # Row index: a list of tested HPO terms
        pheno_idx = pd.Index(result.phenotypes)
        # Column index: multiindex of counts and percentages for all genotype predicate groups
        gt_idx = pd.MultiIndex.from_product(
            iterables=(result.gt_clf.get_categories(), ("Count", "Percent")),
            names=result.gt_clf.class_labels,
        )
        general_info["a_genotype"] = result.gt_clf.class_labels[0]
        general_info["b_genotype"] = result.gt_clf.class_labels[1]
        df = pd.DataFrame(index=pheno_idx, columns=gt_idx)

        for ph_predicate, count in zip(result.pheno_clfs, result.all_counts):
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
        columns_index = result.all_counts[0].columns
        category_a = columns_index[0]
        category_b = columns_index[1]
        for idx, row in df.iterrows():
            hpo_item = str(idx)
            geno_a_ratio = row[(category_a, "Count")]
            geno_a_perc = row[(category_a, "Percent")]
            geno_b_ratio = row[(category_b, "Count")]
            geno_b_perc = row[(category_b, "Percent")]
            p_val = row[("", "p values")]
            adj_p_val = row[("", "Corrected p values")]
            if p_val is None or np.isnan(p_val):
                continue
            elif adj_p_val > 0.05:
                continue
            else:
                with_geno_a = f"{geno_a_ratio} ({geno_a_perc})"
                with_geno_b = f"{geno_b_ratio} ({geno_b_perc})"
                xref_for_hpo = self.get_hpo_xref_string(xref_d=xrefs, hpo_item=hpo_item)
                sig_result_list.append({
                    "hpo_item": hpo_item,
                    "with_geno_a": with_geno_a,
                    "with_geno_b": with_geno_b,
                    "pval": format_p_value(p_val),
                    "adj_pval": format_p_value(adj_p_val),
                })

        general_info["n_sig_results"] = len(sig_result_list)
        general_info["n_tests_performed"] = result.total_tests
        return general_info, sig_result_list



    def _get_all_fisher_result(self, fresult):
        total_hpo_testable = len(fresult.result._n_usable)
        filtered_out = fresult.result.n_filtered_out()
        total_hpo_tested = total_hpo_testable - filtered_out
        a_genotype = fresult.result.gt_clf.class_labels[0]
        b_genotype = fresult.result.gt_clf.class_labels[1]
        nsig = int(fresult.result.n_significant_for_alpha())
        return {"total_hpo_testable":total_hpo_testable, "total_hpo_tested":total_hpo_tested,
                "a_genotype":a_genotype, "b_genotype":b_genotype, "nsig": nsig}

    def _prepare_context(
            self,
            report: GpseaAnalysisReport,
    ) -> typing.Mapping[str, typing.Any]:
        fet_result_list = list()
        all_fet_result_list = list()
        mono_result_list = list()
        fet_results = report.fet_results
        if fet_results is not None:
            for fres in fet_results:
                general_info, sig_result_list = self.fisher_exact_test(fres)
                fet_result_list.append({"general_info": general_info, "sig_result_list": sig_result_list})
                all_fet_result_list.append(self._get_all_fisher_result(fres))
            n_fet_results = len(fet_result_list)
        else:
            n_fet_results = 0
        if report.mono_results is not None:
            n_mono_results = len(report.mono_results)
            for mres in report.mono_results:
                test_result = self.mono_test(mres)
                mono_result_list.append({"a_genotype": test_result["a_genotype"],
                                         "b_genotype": test_result["b_genotype"],
                                         "name": test_result["name"],
                                         "test_name": test_result["test_name"],
                                         "description": test_result["description"],
                                         "variable_name": test_result["variable_name"],
                                         "pval": test_result["pval"],
                                         "interpretation": test_result["interpretation"],
                                         "xrefs": test_result["xrefs"],
                })
        else:
            n_mono_results = 0
        return {
            "cohort_name": report.name,
            "caption": report.caption,
            "gene_caption": report.gene_caption,
            "n_female": report.n_female,
            "n_male": report.n_male,
            "n_unknown_sex": report.n_unknown_sex,
            "n_total_individual_count": report.n_total_individual_count,
            "n_hpo_terms": report.n_hpo_terms,
            "n_measurements": report.n_measurements,
            "n_alive": report.n_alive,
            "n_deceased": report.n_deceased,
            "n_unknown_vital": report.n_unknown_vital,
            "n_with_onset": report.n_with_onset,
            "n_with_age_of_last_encounter": report.n_with_age_of_last_encounter,
            "n_diseases": report.n_diseases,
            "disease_string": report.disease_string.strip(),
            "latex_gene_caption": report.latex_gene_caption,
            "hpo_version": self._hpo.version,
            "gpsea_version": self._gpsea_version,
            "n_fet_results": n_fet_results,
            "fet_result_list": fet_result_list,
            "all_fet_results": all_fet_result_list,
            "n_mono_results": n_mono_results,
            "mono_result_list": mono_result_list,
            }

    def process_latex(self,
                    report: GpseaAnalysisReport,
                    protein_fig: matplotlib.figure.Figure=None,
                    stats_fig: matplotlib.figure.Figure=None):
        context = self._prepare_context(report)
        latex = process_latex_template(context, protein_fig=protein_fig, stats_fig=stats_fig)
        return latex


