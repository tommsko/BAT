import configparser
import json
import logging
import os.path
import subprocess

import requests
from requests import Response
log = logging.getLogger("resolver")


def _replace_placeholders(command: list[str], placeholders: dict[str, str]) -> None:
    """
    Replaces placeholders with values in command to be executed
    :param command: to be fed into subprocess.run
    :param placeholders: string replacements to be done
    :return: None (but command parameter has placeholders replaced)
    :raises KeyError: if not all placeholders were replaced
    """
    placeholders_to_replace: set[str] = set(placeholders.keys())
    for i in range(len(command)):
        command_part: str = command[i]
        for placeholder_from, placeholder_to in placeholders.items():
            if placeholder_from in command_part:
                command[i] = command[i].replace(placeholder_from, placeholder_to)
                placeholders_to_replace.remove(placeholder_from)

    if len(placeholders_to_replace) != 0:
        raise KeyError(
            f"no substitution possible for placeholders: {placeholders_to_replace}"
        )


def _upload_file(pdb_path: str, config: configparser.ConfigParser) -> str | None:
    """
    Uploads PDB file to remote server
    :param pdb_path: to the .pdb file
    :param config: main configuration
    :return: name of the .pdb file, if successful
    :raises Exception: on upload error
    """
    pdb_name: str = os.path.basename(pdb_path)
    log.debug(f"... uploading '{pdb_name}'")

    upload_command: list[str] = [
        _v for _, _v in config.items("PDB_UPLOAD_UPLOAD_COMMAND")
    ]
    _replace_placeholders(
        upload_command, {"{FILEPATH}": pdb_path, "{FILENAME}": pdb_name}
    )
    subprocess.run(upload_command, check=True)
    log.debug("... uploaded successful")
    return pdb_name


def _cleanup_file(pdb_name: str, config: configparser.ConfigParser) -> None:
    """
    Removes PDB file from remote server
    :param pdb_name: path to the .pdb file
    :param config: main configuration
    :return: None
    :raises Exception: on ssh error
    """
    log.debug(f"... deleting '{pdb_name}'")
    purge_command: list[str] = [
        _v for _, _v in config.items("PDB_UPLOAD_DELETE_COMMAND")
    ]
    _replace_placeholders(purge_command, {"{FILENAME}": pdb_name})
    subprocess.run(purge_command, check=True)
    log.debug("... delete successful")


def _fetch_pdb_similarity(pdb_name: str, config: configparser.ConfigParser) -> dict:
    """
    Executes online PDB structural similarity search given uploaded file
    :param pdb_name: name of the PDB file
    :param config: main configuration
    :return: dictionary of search results
    :raises Exception: on POST request error
    """
    STRUCT_URL = config.get("PDB_UPLOAD_URL", "url")
    assert "{FILENAME}" in STRUCT_URL, "illegal template for URL"
    remote_path: str = STRUCT_URL.replace("{FILENAME}", pdb_name)

    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    request = {
        "query": {
            "type": "terminal",
            "service": "structure",
            "parameters": {
                "value": {"url": remote_path, "format": "pdb"},
                "operator": "relaxed_shape_match",
            },
        },
        "return_type": "entry",
    }

    log.debug(f"... quering PDB similarity for '{remote_path}'")
    response: Response = requests.post(url, json=request)
    response_txt: str = response.text
    response.close()
    response_json: dict = json.loads(response_txt)
    return response_json


def _parse_pdb_response(response: dict) -> list[tuple[str, float]]:
    """
    Parses PDB similarity search response and returns list of identifiers and their similarity percentages
    :param response: result of _fetch_pdb_similarity(...)
    :return: list of identifiers and their similarity percentages
    """
    if "result_set" not in response:
        return []

    results: list[tuple[str, float]] = []
    for hit in response["result_set"]:
        results.append((hit["identifier"], hit["score"] * 100))
    return results


def pdb_fetch_similar_structures(
    pdb_path: str, config: configparser.ConfigParser
) -> list[tuple[str, float]]:
    """
    Executes remote PDB similarity search
    :param pdb_path: of local .pdb file to be searched
    :param config: main configuration
    :return: list of PDB identifiers as well as similarity percentages
    """
    log.debug(f"resolving PDB similarity search '{pdb_path}'")
    filename: str = _upload_file(pdb_path, config)
    response: dict = _fetch_pdb_similarity(filename, config)
    similar_proteins: list[tuple[str, float]] = _parse_pdb_response(response)
    _cleanup_file(filename, config)
    return similar_proteins


def pdb_fetch_identical_structures(
    pdb_path: str, config: configparser.ConfigParser
) -> list[str]:
    """
    Executes remote PDB identity search
    :param pdb_path: of local .pdb file to be searched
    :param config: main configuration
    :return: list of PDB identifiers
    """
    log.debug(f"resolving PDB identity search '{pdb_path}'")
    return [
        pdb_id
        for pdb_id, similarity in pdb_fetch_similar_structures(pdb_path, config)
        if similarity == 100
    ]
