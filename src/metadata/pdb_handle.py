import configparser
import json
import logging

import requests

from .metadata_fetcher import MetadataFetcherBase
from ..utils.results_file import ResultsFile
from ..utils.timeout import timeout

log = logging.getLogger("resolver")


def _fetch_full_metadata(pdb_id: str) -> dict | None:
    """
    Attempts to fetch metadata from PDB servers
    :param inchikey: of the molecule
    :return: dict of results if successful, None otherwise
    """
    if '_' in pdb_id:
        pdb_id = pdb_id.split('_')[0]  # whole molecule, not just assembly
    pdb_id = pdb_id.upper()
    try:
        req = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
        assert req.status_code == 200, "PDB search failed"
        data_json: dict = req.json()
        req.close()
        return data_json
    except Exception as exc:
        log.debug(f"Failed to fetch metadata for PDB ID '{pdb_id}': {exc}")
        return None


def _extract_metadata_names(from_metadata: dict) -> list[str]:
    """
    Extract names from already fetched metadata
    :param from_metadata: non-None result of _fetch_full_metadata(...)
    :return: list of names of the molecule
    """
    try:
        return from_metadata['em_entity_assembly']['name'] if isinstance(from_metadata['em_entity_assembly'], dict) else from_metadata['em_entity_assembly'][0]['name']
    except KeyError:
        return []


class PDBMetadataFetcher(MetadataFetcherBase):
    """
    Base class for any metadata fetcher
    """
    FETCHER_NAME = "metadata_pdb"

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

        if not self.config.getboolean('PDB', 'enabled', fallback=False):
            return

        log.debug(f"Fetching PDB metadata for '{results.path}'")
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
