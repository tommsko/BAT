import configparser

import pubchempy
import requests

from ..base.fragment_resolver import IdentifierResolutionError


def _pubchem_cids_to_inchikeys(cids: list[int]) -> list[str]:
    """
    Translates PUBCHEM's molecular CIDs to InChI Keys using PUBCHEM servers (using API)
    :param cids: list of cids to be translated
    :return: equivalent list of InChIs
    """
    inchi_keys: list[str] = []
    for cid in cids:
        if cid == 0:
            continue  # CID is never zero, just a result of bad POST fallback request
        try:
            inchi_keys.append(pubchempy.get_properties("InChIKey", cid)[0]["InChIKey"])
        except Exception:
            ...  # resolution can fail, their databases are not perfect
    return inchi_keys


def _pubchem_fetch_identity_match_post(smiles: str) -> list[int]:
    """
    As a fallback, attempts direct (POST) search of PUBCHEM's servers for identity match of SMILES
    :param smiles: signature of molecule to search
    :raises IdentifierResolutionError on any failure
    :return: list of CID's
    """
    try:
        cids_text: str = requests.post(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/txt",
            data={"smiles": smiles},
        ).text.strip()
        if "status" in cids_text.lower():  # request failed
            return []
        return [int(cid) for cid in cids_text.split("\n") if cid.isnumeric()]
    except Exception as exc:
        if 'ServerBusy' in str(exc):
            return []  # we don't want to fail the whole search just because server is busy
        raise IdentifierResolutionError(
            "PUBCHEM's POST identity search", f"search failed: {exc}"
        )


def pubchem_fetch_identity_match(
    smiles: str, config: configparser.ConfigParser
) -> list[str]:
    """
    Searches PUBCHEM's databases for identity match of SMILES
    :param smiles: identifier of molecule to search
    :param config: main configuration
    :raises IdentifierResolutionError: on PUBCHEM search failure (not finding anything is NOT a failure)
    :return: list of InChI Keys of identical molecules
    """
    if not config.get("PUBCHEM", "enabled"):
        return []

    try:
        try:
            exact_matches: list[pubchempy.Compound] = pubchempy.get_compounds(
                smiles, namespace="smiles", searchtype="identity", as_dataframe=False
            )
            return _pubchem_cids_to_inchikeys([match.cid for match in exact_matches])
        except pubchempy.BadRequestError:  # fallback
            return _pubchem_cids_to_inchikeys(_pubchem_fetch_identity_match_post(smiles))
    except Exception as exc:
        if 'ServerBusy' in str(exc):
            return []  # we don't want to fail the whole search just because server is busy
        raise IdentifierResolutionError(
            "PUBCHEM's API identity search", f"search failed: {exc}"
        )


def _pubchem_fetch_similar_match_post(smiles: str, threshold: int) -> list[int]:
    """
    As a fallback, attempts direct (POST) search of PUBCHEM's servers for similarity match of SMILES
    :param smiles: identifier of molecule to search
    :param threshold: threshold of similarity
    :raises IdentifierResolutionError on any failure
    :return: list of CID's
    """
    try:
        cids_text: str = requests.post(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/cids/txt?Threshold={threshold}",
            data={"smiles": smiles},
        ).text.strip()
        if "status" in cids_text.lower():  # request failed
            return []
        return [int(cid) for cid in cids_text.split("\n") if cid.isnumeric()]
    except Exception as exc:
        if 'ServerBusy' in str(exc):
            return []  # we don't want to fail the whole search just because server is busy
        raise IdentifierResolutionError(
            "PUBCHEM's POST similarity search", f"search failed: {exc}"
        )


def pubchem_fetch_similarity_match(
    smiles: str,
    config: configparser.ConfigParser,
    try_thresholds: list[int],
) -> list[tuple[str, float]]:
    """
    Searches PUBCHEM's databases for similar match of SMILES identifier given the threshold.
    Low threshold will return more matches, but also it may result in timeout given PUBCHEM servers
    :param smiles: identifier of molecule to search
    :param config: main configuration
    :param try_thresholds: list of thresholds to try (will report on lowest not time-outing one)
    :raises IdentifierResolutionError: on PUBCHEM search failure (excluding timeout)
    :return: list of InChI Keys of identical molecules and similarity percentage they were found with
    """

    if not config.get("PUBCHEM", "enabled"):
        return []

    try_thresholds.sort()

    for threshold in try_thresholds:
        try:
            try:
                compounds: list[pubchempy.Compound] = pubchempy.get_compounds(
                    smiles,
                    namespace="smiles",
                    searchtype="similarity",
                    as_dataframe=False,
                    Threshold=threshold,
                )
                cids: list[int] = [match.cid for match in compounds]
            except pubchempy.BadRequestError:  # fallback
                cids: list[int] = _pubchem_fetch_similar_match_post(smiles, threshold)
            return list(
                zip(
                    _pubchem_cids_to_inchikeys(cids),
                    [threshold for _ in range(len(cids))],
                )
            )
        except pubchempy.TimeoutError:  # this is expected
            continue
        except Exception as exc:
            if 'ServerBusy' in str(exc):
                return []  # we don't want to fail the whole search just because server is busy
            raise IdentifierResolutionError(
                "PUBCHEM's API similarity search", f"search failed: {exc}"
            )
    return []
