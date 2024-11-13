import configparser
import os

import MDAnalysis
from src.utils.file_loader import SimulationFile

from .mdanalysis_to_pdb import fragment_to_pdb
from .pdb_handle import pdb_fetch_similar_structures, pdb_fetch_identical_structures
from ..base.fragment_resolver import ResolverBase, ResolverKind
from ..base.resolver_database import ResolverCache
from ..protein.blastp_handle import _mk_random_filename
from ...constants import KNOWN_AMINO_ACIDS
from ...utils.repair_simulation import (
    fix_missing_elements,
    fix_missing_names,
)


class ProteinStructureResolver(ResolverBase):
    """
    ProteinStructureResolver attempts to identify the protein utilizing its 3D structure
    """

    RESOLVER_NAME = "protein-struct"
    RESOLVER_IDENT = "PDB"

    def __init__(
        self,
        cache: ResolverCache,
        configuration: configparser.ConfigParser,
        try_fix_elements: bool = False,
        try_fix_names: bool = False,
    ) -> None:
        """
        Initializes ProteinStructureResolver
        :param cache: main cache
        :param configuration: main configuration
        :param try_fix_elements: attempt to fix missing atomic elements
        :param try_fix_names: attempt to fix missing atomic names
        """
        super().__init__(
            ProteinStructureResolver.RESOLVER_NAME,
            ResolverKind.PRIMARY,
            ProteinStructureResolver.RESOLVER_IDENT,
            cache,
            configuration,
        )

        self.try_fix_missing_elements: bool = try_fix_elements
        self.try_fix_missing_names: bool = try_fix_names

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the ProteinStructureResolver can be used on the given fragment
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """

        fragment_name: str = fragment.segids[0]

        # firstly, it must be a protein, so all residues must have amino acid names
        if not all(
            resname.upper() in KNOWN_AMINO_ACIDS for resname in fragment.resnames
        ):
            return False

        # secondly, atomic elements *MUST* be present or resolved (requirement by the PDB)
        if not simulation.supports_atomic_elements(fragment_name):
            if (
                simulation.supports_atomic_masses(fragment_name)
                and self.try_fix_missing_elements
                and fix_missing_elements(simulation, fragment)
            ):
                ...  # ok!
            else:
                return False

        # next, all atoms must have unique names
        if not simulation.supports_atomic_names(fragment_name):
            if self.try_fix_missing_names and fix_missing_names(simulation, fragment):
                ...  # ok!
            else:
                return False

        # lastly, atomic positions must be provided
        return simulation.supports_atomic_positions(fragment_name)

    def get_signature(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> str:
        """
        Translates fragment into PDB file
        :param simulation: of which is fragment part of
        :param fragment: AtomGroup representing a molecule
        :return: amino-acid sequence of the fragment
        :raises: SignatureGenerationError if fingerprint cannot be created
        """
        if not os.path.exists(self.config.get("PDB", "temporary_directory")):
            os.makedirs(self.config.get("PDB", "temporary_directory"))

        output_path: str = _mk_random_filename(
            "pdb", self.config.get("PDB", "temporary_directory", fallback=os.getcwd())
        )
        fragment_to_pdb(fragment, output_path, verify=True)
        with open(output_path, "r") as fd:
            signature: str = fd.read()
        os.unlink(output_path)

        return signature

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Attempts to resolve signature (.pdb file) into standardized PDB IDs
        Reported proteins (in the form of PDB IDs) have to have exact match to the fingerprint (structure)
        :param signature: to be translated into PDB IDs
        :return: list of standardized PDB IDs of the protein sequence
        :raises IdentifierResolutionError: if PDB search fails (finding no results is not considered a failure)
        """
        if not os.path.exists(self.config.get("PDB", "temporary_directory")):
            os.makedirs(self.config.get("PDB", "temporary_directory"))

        tmp_path: str = _mk_random_filename(
            "pdb", self.config.get("PDB", "temporary_directory", fallback=os.getcwd())
        )
        with open(tmp_path, "w") as fd:
            fd.write(signature)

        try:
            exact_structures: list[str] = pdb_fetch_identical_structures(
                tmp_path, self.config
            )
            os.unlink(tmp_path)
        except Exception:
            os.unlink(tmp_path)
            raise

        return exact_structures

    def try_fetch_identifiers_relaxed_match(self, signature: str) -> list[[str, float]]:
        """
        Attempts to resolve signature (.pdb file) into standardized PDB IDs
        Reported proteins (in the form of PDB IDs) have to have partial match to the fingerprint (protein sequence)
        :param signature: to be translated into PDB IDs
        :return: list of standardized PDB IDs of the protein sequence
        :raises IdentifierResolutionError: if PDB search fails (finding no results is not considered a failure)
        """
        if not os.path.exists(self.config.get("PDB", "temporary_directory")):
            os.makedirs(self.config.get("PDB", "temporary_directory"))
        tmp_path: str = _mk_random_filename(
            "pdb", self.config.get("PDB", "temporary_directory", fallback=os.getcwd())
        )

        with open(tmp_path, "w") as fd:
            fd.write(signature)

        try:
            similar_structures: list[tuple[str, float]] = pdb_fetch_similar_structures(
                tmp_path, self.config
            )
            os.unlink(tmp_path)
        except Exception:
            os.unlink(tmp_path)
            raise

        return similar_structures

    def generate_debug_data(
        self,
        simulation: SimulationFile,
        fragment: MDAnalysis.AtomGroup,
        out_path: str,
        out_name: str,
    ) -> bool:
        """
        Generates .pdb file of the protein
        :param out_path: directory where to save generated representations
        :param out_name: basename of the file to save representations to (will be suffixed with kind of export)
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True if export went without problems
        """
        pdb_content: str = self.get_signature(simulation, fragment)
        with open(
            os.path.join(out_path, f"{out_name}_PROTEIN_STRUCT_structure.pdb"), "w"
        ) as fd:
            fd.write(pdb_content)
        return True
