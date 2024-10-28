import abc
import typing

from gpsea.model import Cohort
from gpsea.analysis import MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult


class GPSEAAnalysisResult:
    
    def __init__(
        self,
        cohort: Cohort,
        results: typing.Iterable[typing.Union[MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult]],
        interpretation: typing.Optional[str] = None,
    ):
        self._cohort = cohort
        for i, r in enumerate(results):
            assert isinstance(r, (MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult)), f"#{i} failed the Q/C"
        self._results = tuple(results)
        self._interpretation = interpretation

    @property
    def cohort(self) -> Cohort:
        return self._cohort

    @property
    def results(self) -> typing.Collection[typing.Union[MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult]]:
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
        result: GPSEAAnalysisResult,
        out: typing.TextIO,
    ):
        pass
