import abc
import io
import typing
import pandas as pd
import numpy as np

import hpotk
from gpsea.model import Cohort
from gpsea.analysis import MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult
from .util import open_text_io_handle_for_writing
from jinja2 import Environment, PackageLoader


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
            caption: str,
            fet_results: typing.Iterable[MultiPhenotypeAnalysisResult]=None,
            mono_results: typing.Iterable[MonoPhenotypeAnalysisResult]=None,  
    ):
        self._name = name
        self._cohort = cohort
        if fet_results is not None:
            for i, r in enumerate(fet_results):
                assert isinstance(r, MultiPhenotypeAnalysisResult), f"#{i} must be `MultiPhenotypeAnalysisResult but was {type(r)}`"
            self._fet_results = tuple(fet_results)
        if mono_results is not None:
            for i, r in enumerate(mono_results):
                assert isinstance(r, MonoPhenotypeAnalysisResult), f"#{i} must be `MonoPhenotypeAnalysisResult`"
            self._mono_results = tuple(mono_results)
        self._caption = generate_cohort_summary(cohort=cohort, gene_symbol=gene_symbol, mane_tx_id=mane_tx_id, mane_protein_id=mane_protein_id, caption=caption)



    @property
    def name(self) -> str:
        return self._name

    @property
    def cohort(self) -> Cohort:
        return self._cohort

    @property
    def fet_results(self) -> typing.Collection[MultiPhenotypeAnalysisResult]:
        return self._fet_results
    
    @property
    def mono_results(self) -> typing.Collection[MonoPhenotypeAnalysisResult]:
        return self._mono_results

    @property
    def caption(self) -> str:
        return self._caption

def generate_cohort_summary(cohort, 
                            gene_symbol:str = None, 
                            mane_tx_id:str = None, 
                            mane_protein_id:str = None, 
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
    n_variants = len(set(cohort.all_variants()))
    if gene_symbol is None or mane_tx_id is None or mane_protein_id is None:
        variants = f"A total of {n_variants} unique variant alleles were found."
    else:
        variants = f"A total of {n_variants} unique variant alleles were found in {gene_symbol} (transcript: {mane_tx_id}, protein id: {mane_protein_id})."
    if caption is None:
        return f"{counts} {n_hpo_terms} {disease_desc} {variants}"
    else:
        return f"{counts} {n_hpo_terms} {disease_desc} {variants} {caption}"


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
        return HtmlGpseaNotebookSummarizer(html)

    
    def fisher_exact_test(self, result: MultiPhenotypeAnalysisResult) -> typing.Tuple[typing.Dict[str,object], typing.List[typing.Dict[str, object]]]:
        """
        The result is typically a list of Fisher Exact Test results. We will keep only
        those with an adjusted p value of 0.05 or lower. We first add the values to
        a pandas dataframe for convenience, and transfer data for signficiant values
        to TestResult objects to simplify display.
        """
        assert isinstance(result, MultiPhenotypeAnalysisResult)
        general_info = dict()
        sig_result_list = list()
        fet_test_results = list()
        n_usable_fet = result.n_usable
        total_tested = result.total_tests
        # Row index: a list of tested HPO terms
        pheno_idx = pd.Index(result.phenotypes)
        # Column index: multiindex of counts and percentages for all genotype predicate groups
        gt_idx = pd.MultiIndex.from_product(
            iterables=(result.gt_predicate.get_categories(), ("Count", "Percent")),
            names=result.gt_predicate.group_labels,
        )
        general_info["a_genotype"] = result.gt_predicate.group_labels[0]
        general_info["b_genotype"] = result.gt_predicate.group_labels[1]


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
                sig_result_list.append({
                    "hpo_item": hpo_item,
                    "with_geno_a": with_geno_a,
                    "with_geno_b": with_geno_b,
                    "pval": p_val,
                    "adj_pval": adj_p_val,
                })
            
        general_info["n_sig_results"] = len(sig_result_list)
        general_info["n_tests_performed"] = result.total_tests
        print(f"extracted {len(fet_test_results)} FET significant results")
        return general_info, sig_result_list
    
    
    
    def _prepare_context(
            self,
            report: GpseaAnalysisReport,
    ) -> typing.Mapping[str, typing.Any]:
        
        fet_results = report.fet_results
        n_fet_results = len(fet_results)
        fet_result_list = list()
        for fres in fet_results:
            general_info, sig_result_list = self.fisher_exact_test(fres)
            fet_result_list.append({"general_info": general_info, "sig_result_list": sig_result_list})
        n_fet_results = len(fet_result_list)

        return {
            "cohort_name": report.name,
            "caption": report.caption,
            "hpo_version": self._hpo.version,
            "gpsea_version": self._gpsea_version,
            "n_fet_results": n_fet_results,
            "fet_result_list": fet_result_list,
            }
