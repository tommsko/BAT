import configparser
import json
import logging
import os.path
import subprocess
import time

import requests
from requests import Response
log = logging.getLogger("resolver")


def _alphafind_curl(pdb_path: str, request_timeout: int) -> str:
    """
    Executes AlphaFind CURL request
    :param pdb_path: to pdb file
    :param request_timeout: maximum number of seconds for request
    :return: string response
    """
    try:
        return subprocess.check_output(
            ['curl',
             '--request', 'POST',
             '--url', 'https://api.stage.alphafind.dyn.cloud.e-infra.cz/search?limit=10',
             '--header', 'content-type: multipart/form-data',
             '--form', f'file=@{pdb_path}'], timeout=request_timeout).decode()
    except Exception:
        return ""


def _alphafind_fetch_response(
    pdb_path: str, config: configparser.ConfigParser
) -> dict | None:
    """
    Executes remote AlphaFind similarity search
    :param pdb_path: of local .pdb file to be searched
    :param config: main configuration
    :return: json response if successful, None on timeout
    """
    log.debug(f"resolving AlphaFind similarity search '{pdb_path}'")

    if not config.getboolean("ALPHAFIND", 'enabled'):
        return None

    retry_count: int = config.getint("ALPHAFIND", 'retry_count')
    retry_delay: int = config.getint("ALPHAFIND", 'retry_delay_s')
    request_timeout: int = config.getint("ALPHAFIND", 'request_timeout_s')

    _alphafind_curl(pdb_path, request_timeout)  # initial request

    for _ in range(retry_count):
        time.sleep(retry_delay)
        try:
            response: dict = json.loads(_alphafind_curl(pdb_path, request_timeout))
            if 'results' in response:
                return response
        except Exception:
            ...

    return None


def alphafind_similarity(pdb_path: str, config: configparser.ConfigParser) -> list[tuple[str, float]]:
    """
    Executes alphafind similarity search
    :param pdb_path: of local .pdb file to be searched
    :param config: main configuration
    :return: list of alphafold identifiers with associated aligned_percentage
    """

    response: dict | None = _alphafind_fetch_response(pdb_path, config)
    if response is None or 'results' not in response:
        return []

    results: list[tuple[str, float]] = []
    for entry in response['results']:
        if 'object_id' in entry and 'aligned_percentage' in entry:
            results.append((entry['object_id'], entry['aligned_percentage']))
    return results


def alphafind_identity(pdb_path: str, config: configparser.ConfigParser) -> list[str]:
    """
    Executes alphafind identity search
    :param pdb_path: of local .pdb file to be searched
    :param config: main configuration
    :return: list of alphafold identifiers
    """

    return [entry[0] for entry in alphafind_similarity(pdb_path, config) if entry[1] == 1.0]
