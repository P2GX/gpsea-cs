import abc
import typing

from gpsea.model import Cohort
from gpsea.analysis import MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult


class GPAnalysisResult:

    @staticmethod
    def from_mono(
        result: MonoPhenotypeAnalysisResult,
        xrefs: typing.Optional[typing.Iterable[str]] = None,
        interpretation: typing.Optional[str] = None,
    ) -> "GPAnalysisResult":
        if xrefs is None:
            xrefs = {}
        else:
            xrefs = {
                result.phenotype.variable_name: xrefs,
            }
        return GPAnalysisResult(
            result=result,
            xrefs=xrefs,
            interpretation=interpretation,
        )

    @staticmethod
    def from_multi(
        result: MultiPhenotypeAnalysisResult,
        xrefs: typing.Optional[typing.Mapping[str, typing.Iterable[str]]] = None,
        interpretation: typing.Optional[str] = None,
    ) -> "GPAnalysisResult":
        if xrefs is None:
            xrefs = {}
        
        return GPAnalysisResult(
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
        cohort: Cohort,
        results: typing.Iterable[GPAnalysisResult],
    ):
        self._cohort = cohort
        for i, r in enumerate(results):
            assert isinstance(r, GPAnalysisResult), f"#{i} must be `GPAnalysisResult`"
        self._results = tuple(results)

    @property
    def cohort(self) -> Cohort:
        return self._cohort

    @property
    def results(self) -> typing.Collection[GPAnalysisResult]:
        return self._results


class GPSEAAnalysisResultSummarizer(metaclass=abc.ABCMeta):
    """
    `GPSEAAnalysisResultSummarizer` formats the :class:`~gpseacs.report.GPSEAAnalysisResult`
    into an output in a chosen format (LaTeX, HTML, whatever...).
    """
    
    @abc.abstractmethod
    def summarize_result(
        self,
        result: GpseaAnalysisReport,
        out: typing.TextIO,
    ):
        pass
