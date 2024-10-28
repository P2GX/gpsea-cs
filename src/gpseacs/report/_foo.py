from typing import TextIO

from gpsea.analysis import MonoPhenotypeAnalysisResult, MultiPhenotypeAnalysisResult
from gpsea.model import Cohort

from ._api import GPSEAAnalysisResultSummarizer, GPSEAAnalysisResult


class FooGPSEAAnalysisResultSummarizer(GPSEAAnalysisResultSummarizer):
    """
    Prototype for figuring out the requirements of the report.
    """
    
    def summarize_result(
        self,
        result: GPSEAAnalysisResult,
        out: TextIO,
    ):
        self._summarize_cohort(
            cohort=result.cohort,
            out=out,
        )

        for ar in result.results:
            if isinstance(ar, MonoPhenotypeAnalysisResult):
                self._summarize_mono(result, out)
            elif isinstance(ar, MultiPhenotypeAnalysisResult):
                self._summarize_multi(result, out)
            else:
                raise ValueError("What is that thing you gave me??!")

    def _summarize_cohort(
        self,
        cohort: Cohort,
        out: TextIO,
    ):
        pass

    def _summarize_mono(
        self,
        result: MonoPhenotypeAnalysisResult,
        out: TextIO,
    ):
        
        pass
    
    def _summarize_multi(
        self,
        result: MultiPhenotypeAnalysisResult,
        out: TextIO,
    ):
        pass
    
