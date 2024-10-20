import typing
import io
import abc
from .summarizer import GpseaSummarizer
from jinja2 import Environment, PackageLoader
from .util import open_text_io_handle_for_writing, open_text_io_handle_for_reading


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
class HtmlGpseaReport(GpseaReport):
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


class HtmlVisualizer:

    def __init__(self,):
        environment = Environment(loader=PackageLoader('gpseacs', 'templates'))
        self._cohort_template = environment.get_template("summary.html")

    def process(
            self,
            summarizer: GpseaSummarizer,
    ) -> GpseaReport:
        """
        Create an HTML that should be shown with display(HTML(..)) of the ipython package.

        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto

        Returns:
            GpseaReport: a report that can be stored to a path or displayed in
                interactive environment such as Jupyter notebook.
        """
        context = self._prepare_context(summarizer)
        report = self._cohort_template.render(context)
        return HtmlGpseaReport(html=report)

    def _prepare_context(
            self,
            summarizer: GpseaSummarizer,
    ) -> typing.Mapping[str, typing.Any]:

        hpo_counts = list()
        gpsea_version = summarizer.gpsea_version
        hpo_version = summarizer.hpo_version
        # The following dictionary is used by the Jinja2 HTML template
        return {
            "gpsea_version": gpsea_version,
            "hpo_version": hpo_version,
        }

