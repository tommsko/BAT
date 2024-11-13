import configparser
import logging

from ..utils.file_loader import SimulationFile
from ..utils.results_file import ResultsFile

log = logging.getLogger("resolver")


class MetadataFetcherBase:
    """
    Base class for any metadata fetcher
    """

    def __init__(self, fetcher_name: str, config: configparser.ConfigParser) -> None:
        """
        Initializes MetadataFetcher
        :param fetcher_name: name of the fetcher
        :param config: main configuration
        """
        self.fetcher_name: str = fetcher_name
        self.config: configparser.ConfigParser = config

        self.enabled: bool = config.getboolean(self.fetcher_name, 'enabled', fallback=False)
        self.fetch_full_metadata: bool = config.getboolean(self.fetcher_name, 'fetch_full_metadata', fallback=False)
        self.ident_type: str = config.get(self.fetcher_name, 'processes_ident_type')
        self.only_identity: bool = config.getboolean(self.fetcher_name, 'only_identity_matches', fallback=False)
        self.fetch_timeout_sec: int | None = config.getint(self.fetcher_name, 'timeout_sec')
        if self.fetch_timeout_sec <= 0:
            self.fetch_timeout_sec = None

    def try_fetch_metadata(self, results: ResultsFile, overwrite_existing: bool) -> bool:
        """
        Attempts to fetch metadata for entities discovered in results file
        :param results: ResultsFile to fetch metadata for
        :param overwrite_existing: if True, existing metadata will be overwriten
        :return: True if successful, False otherwise
        """
        if not self.enabled:
            return
        if not results.load(failsafe=False, segments_from_results=True, ignore_compatibility=True):
            return False
        if 'metadata' not in results.data:
            results.data['metadata'] = {}
        self._try_fetch_append_metadata(results, overwrite_existing)
        results.save()

    def _get_segment_tasks(self, results: ResultsFile, overwrite_existing: bool, ident_type: str, only_identity: bool) -> list[tuple[str, str]]:
        """
        Determines which InChIKeys (and metadata segments) have to be fetched
        :param results: ResultsFile
        :param overwrite_existing: if True, ignores previously fetched metadata
        :param ident_type: required identifier type of the segment resolution
        :param only_identity: report segment only if it was resolved on identity level
        :return: list of inchikeys and their metadata sections
        """
        tasks: list[tuple[str, str]] = []
        for segment in results.segments:
            if results.get_idents_type(segment) is None \
                    or results.get_idents_type(segment) not in (ident_type, 'hash'):
                continue
            idents_to_process: list[str] = results.get_idents(segment)
            if not only_identity:
                idents_to_process.extend([x[0] for x in results.get_similarity_idents(segment)])
            for ident in idents_to_process:
                if ':' in ident:  # hash
                    _ident_type, _ident_value = ident.split(":")
                    if _ident_type != ident_type:
                        continue
                    ident = _ident_value
                associated_section: str = f"{segment}_{ident}"
                if not overwrite_existing and associated_section in results.data['metadata']:
                    continue
                tasks.append((ident, associated_section))
        return tasks


    def _try_fetch_append_metadata(self, data: ResultsFile, overwrite_existing: bool) -> None:
        """
        Attempts to fetch metadata for existing information in results
        :param data: currently provided data
        :param overwrite_existing: if True, ignores previously fetched metadata
        :return: None
        """
        raise NotImplementedError

