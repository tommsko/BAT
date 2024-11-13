import configparser
import json
import logging

from chembl_webresource_client.settings import Settings
Settings.Instance().TIMEOUT = 10
Settings.Instance().CACHING = False
Settings.Instance().TOTAL_RETRIES = 10
from chembl_webresource_client.new_client import new_client
from .metadata_fetcher import MetadataFetcherBase
from ..utils.results_file import ResultsFile
from ..utils.timeout import timeout
log = logging.getLogger("resolver")


def _fetch_full_metadata(inchikey: str) -> dict | None:
    """
    Attempts to fetch metadata from CHEMBL servers
    :param inchikey: of the molecule
    :return: dict of results if successful, None otherwise
    """
    try:
        molecule = new_client.molecule
        metadata_full: list = molecule.filter(molecule_structures__standard_inchi_key__iexact=inchikey)
        assert len(metadata_full) > 0, "no records for inchikey found"
        return metadata_full[0]
    except Exception as exc:
        log.debug(f"Failed to fetch metadata for InChIKey '{inchikey}': {exc}")
        return None


def _extract_metadata_names(from_metadata: dict) -> list[str]:
    """
    Extract names and synonyms from already fetched metadata
    :param from_metadata: non-None result of _fetch_full_metadata(...)
    :return: list of names of the molecule
    """
    names = [from_metadata['pref_name']]
    for synonym in from_metadata['molecule_synonyms']:
        names.append(synonym['molecule_synonym'])
    return names


class ChemblMetadataFetcher(MetadataFetcherBase):
    """
    Base class for any metadata fetcher
    """
    FETCHER_NAME = "metadata_chembl"

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
        if not self.config.getboolean('CHEMBL', 'enabled', fallback=False):
            return

        log.debug(f"Fetching CHEMBL metadata for '{results.path}'")
        for ident, metadata_section in self._get_segment_tasks(results, overwrite_existing, self.ident_type, self.only_identity):
            try:
                with timeout(self.fetch_timeout_sec):
                    metadata: dict | None = _fetch_full_metadata(ident)
            except TimeoutError:
                continue

            if metadata is None:
                results.data['metadata'][metadata_section] = 'ERROR'
                continue
            if self.fetch_full_metadata:
                results.data['metadata'][metadata_section] = metadata
            else:
                results.data['metadata'][metadata_section] = {
                    'names': _extract_metadata_names(metadata)
                }
