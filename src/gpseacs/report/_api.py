import abc
import io


class GPSEAAnalysisResult:
    pass


class GPSEAAnalysisResultSummarizer(metaclass=abc.ABCMeta):
    """
    `GPSEAAnalysisResultSummarizer` formats the :class:`~gpseacs.report.GPSEAAnalysisResult`
    into an output in a chosen format (LaTeX, HTML, whatever...).
    """
    
    @abc.abstractmethod
    def summarize(
        self,
        result: GPSEAAnalysisResult,
        out: io.IOBase,
    ):
        pass
