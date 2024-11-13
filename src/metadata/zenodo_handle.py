import configparser
import logging

import requests

from .metadata_fetcher import MetadataFetcherBase
from ..utils.results_file import ResultsFile
from ..utils.timeout import timeout

log = logging.getLogger("resolver")


def _fetch_full_metadata(zenodo_id: str) -> dict | None:
    """
    Attempts to fetch metadata from PDB servers
    :param inchikey: of the molecule
    :return: dict of results if successful, None otherwise
    """
    try:
        req = requests.get(f"https://zenodo.org/api/records/{zenodo_id}")
        assert req.status_code == 200, "ZENODO search failed"
        data_json: dict = req.json()
        req.close()
        return data_json
    except Exception as exc:
        log.debug(f"Failed to fetch metadata for PDB ID '{zenodo_id}': {exc}")
        return None


def _extract_metadata_reduced(from_metadata: dict) -> dict[str, str]:
    """
    Extract reduced metadata from already fetched metadata
    :param from_metadata: non-None result of _fetch_full_metadata(...)
    :return: dict of reduced metadata
    """
    return {
        'name': from_metadata['metadata']['title'],
        'doi_url': from_metadata['doi_url'],
        'published': from_metadata['created'],
    }


class ZenodoMetadataFetcher(MetadataFetcherBase):
    """
    Base class for any metadata fetcher
    """
    FETCHER_NAME = "metadata_zenodo"

    def __init__(self, config: configparser.ConfigParser) -> None:
        """
        Initializes MetadataFetcher
        :param fetcher_name: name of the fetcher
        :param results: results file to fetch metadata for
        :param config: main configuration
        """
        super().__init__(self.FETCHER_NAME, config)

    def _try_fetch_append_metadata(self, results: ResultsFile, overwrite_existing: bool) -> None:
        """
        Attempts to fetch metadata for existing information in results
        :param results: ResultsFile
        :param overwrite_existing: if True, ignores previously fetched metadata
        :return: None
        """
        if not self.config.getboolean('ZENODO', 'enabled', fallback=False):
            return

        log.debug(f"Fetching ZENODO metadata for '{results.path}'")
        if not overwrite_existing and 'dataset' in results.data['metadata']:
            return
        try:
            with timeout(self.fetch_timeout_sec):
                metadata: dict | None = _fetch_full_metadata(results.data['name'].split('_')[0])
        except TimeoutError:
            return
        if metadata is None:
            results.data['metadata']['dataset'] = 'ERROR'
            return
        if self.fetch_full_metadata:
            results.data['metadata']['dataset'] = metadata
        else:
            results.data['metadata']['dataset'] = _extract_metadata_reduced(metadata)
