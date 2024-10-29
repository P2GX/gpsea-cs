import abc
import typing

from gpsea.model import Cohort
from gpsea.analysis import MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult


class GPAnalysisResult:

    @staticmethod
    def mono(
        result: MonoPhenotypeAnalysisResult,
        xrefs: typing.Optional[typing.Iterable[str]] = None,
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
        )

    @staticmethod
    def multi(
        result: MultiPhenotypeAnalysisResult,
        xrefs: typing.Optional[typing.Mapping[str, typing.Iterable[str]]] = None,
    ) -> "GPAnalysisResult":
        if xrefs is None:
            xrefs = {}
        
        return GPAnalysisResult(
            result=result,
            xrefs=xrefs,
        )

    def __init__(
        self,
        result: typing.Union[MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult],
        xrefs: typing.Mapping[str, typing.Collection[str]],
    ):
        assert isinstance(result, (MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult))
        self._result = result

        self._xrefs = dict(xrefs)

    @property
    def result(self) -> typing.Union[MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult]:
        return self._result

    @property
    def xrefs(self) -> typing.Mapping[str, typing.Collection[str]]:
        """
        Get mapping from phenotype `variable_name` to G/P association cross references.
        """
        return self._xrefs

    def is_mono(self) -> bool:
        return isinstance(self._result, MonoPhenotypeAnalysisResult)

    def is_multi(self) -> bool:
        return isinstance(self._result, MultiPhenotypeAnalysisResult)


class GpseaAnalysisReport:
    
    def __init__(
        self,
        cohort: Cohort,
        results: typing.Iterable[GPAnalysisResult],
        interpretation: typing.Optional[str] = None,
    ):
        self._cohort = cohort
        for i, r in enumerate(results):
            assert isinstance(r, GPAnalysisResult), f"#{i} must be `GPAnalysisResult`"
        self._results = tuple(results)
        self._interpretation = interpretation

    @property
    def cohort(self) -> Cohort:
        return self._cohort

    @property
    def results(self) -> typing.Collection[GPAnalysisResult]:
        return self._results
    
    @property
    def interpretation(self) -> typing.Optional[str]:
        return self._interpretation


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
