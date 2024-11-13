import configparser
import re
from urllib.request import Request, urlopen

import pubchempy
import requests
from bs4 import BeautifulSoup as soup
from chembl_webresource_client.new_client import new_client

from ..base.fragment_resolver import IdentifierResolutionError
from ..chem.pubchem_handle import _pubchem_cids_to_inchikeys

no_html = re.compile(
    "<.*?>"
)  # https://stackoverflow.com/questions/9662346/python-code-to-remove-html-tags-from-a-string


def fetch_pdb_ligand_inchikey(
    ligand_name: str, config: configparser.ConfigParser
) -> list[str]:
    """
    Attempts to fetch PDB ligand InChI Key from ligand name
    :param ligand_name: to be fetched (segment name)
    :param config: main configuration
    :return: list of InChI Keys fetched (always one or zero)
    """

    if not config.getboolean("PDB", "enabled"):
        return []

    try:
        response = requests.get(
            f"https://www.rcsb.org/ligand/{ligand_name.upper()}",
            timeout=config.getint("PDB", "timeout_sec", fallback=30),
        )
        response_parsed = soup(response.content, "html.parser")
        response.close()
        inchikeys = response_parsed.find_all("tr", id="chemicalInChIKey")
        for inchikey in inchikeys:
            return [re.sub(no_html, "", str(inchikey.contents[1])).strip()]
    except Exception as exc:
        raise IdentifierResolutionError(
            "searching PDB database", f"search failed: {exc}"
        )


def fetch_pubchem_ligand_inchikey(
    ligand_name: str, config: configparser.ConfigParser
) -> list[tuple[str, float]]:
    """
    Attempts to fetch CHEMBL ligand InChI Key from ligand name
    :param ligand_name: to be fetched (segment name)
    :param config: main configuration
    :return: list of InChI Keys fetched (always one or zero), with similarity of 50 (because the number doesn't make sense)
    """

    if not config.getboolean("PUBCHEM", "enabled"):
        return []

    try:
        inchikeys = _pubchem_cids_to_inchikeys(
            pubchempy.get_cids(ligand_name, "name", "substance", list_return="flat")
        )
        return list(zip(inchikeys, [50 for _ in range(len(inchikeys))]))
    except Exception as exc:
        raise IdentifierResolutionError(
            "searching PUBCHEM database", f"search failed: {exc}"
        )
