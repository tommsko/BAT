import configparser
import logging

import MDAnalysis

from .pdb_pubchem_ligand_handle import fetch_pdb_ligand_inchikey, fetch_pubchem_ligand_inchikey
from ..base.fragment_resolver import (
    ResolverBase, ResolverKind
)
from ..base.resolver_database import ResolverCache
from src.utils.file_loader import SimulationFile

log = logging.getLogger("resolver")


class ChemLigandNameResolver(ResolverBase):
    RESOLVER_NAME = "chem_ligand_name"
    RESOLVER_IDENT = "InChI"

    def __init__(
        self,
        cache: ResolverCache,
        configuration: configparser.ConfigParser,
    ) -> None:
        """
        Initializes the ChemLigandNameResolver
        :param cache: resolver cache
        :param configuration: main configuration
        """
        super().__init__(
            ChemLigandNameResolver.RESOLVER_NAME,
            ResolverKind.PRIMARY,
            ChemLigandNameResolver.RESOLVER_IDENT,
            cache,
            configuration,
        )

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the ChemLigandNameResolver can be used on the given fragment
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """

        return True  # always can be, as we are looking at name only

    def get_signature(self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Translates fragment into fragment molecule name
        :param simulation: of which is fragment part of
        :param fragment: AtomGroup representing a molecule
        :return: fragment name
        :raises: SignatureGenerationError if fingerprint cannot be created
        """

        # usual namoculture is SEG_<<ID>>_<<MOLECULE_NAME>>
        return fragment.segments.segids[0].split('_')[-1].strip().upper()

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Attempts to resolve fragment name into standardized InChI Key, as an exact name match of PDB ligand
        :param signature: to be translated into standardized InChI key
        :return: list of standardized InChI keys of the molecule
        :raises IdentifierResolutionError: if any search fails (finding no results is NOT considered a failure)
        """
        if len(signature) < self.config.getint(self.resolver_name, 'minimum_name_len'):
            return []
        return fetch_pdb_ligand_inchikey(signature, self.config)

    def try_fetch_identifiers_relaxed_match(
        self, signature: str
    ) -> list[tuple[str, float]]:
        """
        Attempts to resolve fragment name into standardized InChI Key, as an similarity name match of PUBCHEM compound
        :param signature: to be translated into standardized InChI keys
        :return: list of standardized InChI keys of the molecule and their scores (similarity %)
        :raises IdentifierResolutionError: if any search fails (finding no results is not considered a failure)
        """
        if len(signature) < self.config.getint(self.resolver_name, 'minimum_name_len'):
            return []
        return fetch_pubchem_ligand_inchikey(signature, self.config)

    def generate_debug_data(
        self,
        simulation: SimulationFile,
        fragment: MDAnalysis.AtomGroup,
        out_path: str,
        out_name: str,
    ) -> bool:
        """
        Does nothing (fragment name is already exported as part of signature generation)
        :param out_path: directory where to save generated representations
        :param out_name: basename of the file to save representations to (will be suffixed with kind of export)
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True
        """
        return True