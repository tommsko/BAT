import configparser

from chembl_webresource_client.settings import Settings

Settings.Instance().TIMEOUT = 10
Settings.Instance().CACHING = False
Settings.Instance().TOTAL_RETRIES = 10
from chembl_webresource_client.new_client import new_client
from ..base.fragment_resolver import IdentifierResolutionError


def chembl_fetch_similarity_match(
    smiles: str, threshold: int, config: configparser.ConfigParser, include_similarity: bool = False
) -> list[str | tuple[str, float]]:
    """
    Searches CHEMBL's databases for similarity matches of SMILES
    :param smiles: identifier to be searched
    :param threshold: of similarity to report results (100 == identity)
    :param config: main configuration
    :param include_similarity: if True, similarity percentages will be included in the result
    :raises IdentifierResolutionError: on CHEMBL search failure
    :return: list of all InChI Keys of molecules similar to smiles query given threshold
    """

    if not config.getboolean('CHEMBL', 'enabled', fallback=False):
        return []

    # https://notebooks.gesis.org/binder/jupyter/user/chembl-chembl_webresource_client-5neb7ifo/notebooks/demo_wrc.ipynb

    try:
        results: list[str | tuple[str, float]] = []
        similarity = new_client.similarity
        records = similarity.filter(smiles=smiles, similarity=threshold).only(
            ["molecule_structures", "similarity"]
        )

        for record in records:
            if (
                "similarity" not in record
                or "molecule_structures" not in record
                or "standard_inchi_key" not in record["molecule_structures"]
            ):
                continue
            similarity: float = record["similarity"]
            inchi_key: str = record["molecule_structures"]["standard_inchi_key"]
            if include_similarity:
                results.append((inchi_key, similarity))
            else:
                results.append(inchi_key)

        return results
    except Exception as exc:
        raise IdentifierResolutionError(
            "searching CHEMBL database", f"search failed: {exc}"
        )
