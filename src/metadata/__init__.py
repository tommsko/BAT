import configparser
import logging

from .metadata_fetcher import MetadataFetcherBase
from .chembl_handle import ChemblMetadataFetcher
from .pdb_handle import PDBMetadataFetcher
from .zenodo_handle import ZenodoMetadataFetcher
from ..utils.file_loader import SimulationFile
from ..utils.results_file import ResultsFile


def get_metadata_fetchers(config: configparser.ConfigParser) -> list[MetadataFetcherBase]:
    """
    Provides all supported metadata fetchers
    :param config: main configuration
    :return: list of all metadata fetchers, ready to be used
    """
    return [ChemblMetadataFetcher(config),
            PDBMetadataFetcher(config),
            ZenodoMetadataFetcher(config)]


def fetch_metadata_all(simulations: list[SimulationFile], results_dir: str, config: configparser.ConfigParser, overwrite_existing: bool) -> None:
    """
    Fetches all metadata given
    :param simulations: that have been resolved
    :param results_dir: directory where are the results
    :param config: main configuration
    :param overwrite_existing: if False, skips already fetched metadata
    :return: None
    """
    log = logging.getLogger("resolver")
    log.info(f"Fetching metadata for {len(simulations)} simulations...")
    result_files: list[ResultsFile] = [ResultsFile(sim, results_dir) for sim in simulations]
    all_fetchers: list[MetadataFetcherBase] = get_metadata_fetchers(config)

    for file in result_files:
        log.debug(f"Fetching metadata for '{file.path}'...")
        for fetcher in all_fetchers:
            fetcher.try_fetch_metadata(file, overwrite_existing)


