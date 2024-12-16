import configparser
import os

import MDAnalysis
from src.utils.file_loader import SimulationFile

from ..base.fragment_resolver import ResolverBase, ResolverKind
from ..base.resolver_database import ResolverCache
from ..protein.blastp_handle import _mk_random_filename
from ..protein_struct.mdanalysis_to_pdb import fragment_to_partial_pdb
from ..protein_struct.pdb_handle import (
    pdb_fetch_similar_structures,
    pdb_fetch_identical_structures,
)
from ...constants import KNOWN_AMINO_ACIDS
from ...utils.repair_simulation import (
    fix_missing_names,
)

from ..protein_struct.alphafind_handle import alphafind_similarity, alphafind_identity


class ProteinResidueStructureResolver(ResolverBase):
    """
    ProteinResidueStructureResolver behaves exactly same as protein_struct, but with way lighter preconditions
    (allows for generating very limited .pdb files and hopes search engine can do some magic)
    """

    RESOLVER_NAME = "protein-residue-struct"
    RESOLVER_IDENT = "PDB"

    def __init__(
        self,
        cache: ResolverCache,
        configuration: configparser.ConfigParser,
        try_fix_elements: bool = False,
        try_fix_names: bool = False,
        use_alphafind: bool = False,
    ) -> None:
        """
        Initializes ProteinResidueStructureResolver
        :param cache: main cache
        :param configuration: main configuration
        :param try_fix_elements: attempt to fix missing atomic elements
        :param try_fix_names: attempt to fix missing atomic names
        """
        super().__init__(
            ProteinResidueStructureResolver.RESOLVER_NAME,
            ResolverKind.PRIMARY,
            ProteinResidueStructureResolver.RESOLVER_IDENT,
            cache,
            configuration,
        )

        self.try_fix_missing_elements: bool = try_fix_elements
        self.try_fix_missing_names: bool = try_fix_names
        self.use_alphafind: bool = use_alphafind

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the ProteinResidueStructureResolver can be used on the given fragment
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

        # we skip atomic element recognition here

        # next, all atoms must have unique names
        if not simulation.supports_atomic_names(fragment_name):
            if self.try_fix_missing_names and fix_missing_names(simulation, fragment):
                ...  # ok!
            else:
                return False

        # lastly, we really, really need atomic positions
        return simulation.supports_atomic_positions(fragment_name)

    def get_signature(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> str:
        """
        Translates fragment into PDB file
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: amino-acid sequence of the fragment
        :raises: SignatureGenerationError if fingerprint cannot be created
        """

        if not os.path.exists(self.config.get("PDB", "temporary_directory")):
            os.makedirs(self.config.get("PDB", "temporary_directory"))

        output_path: str = _mk_random_filename(
            "pdb", self.config.get("PDB", "temporary_directory", fallback=os.getcwd())
        )
        fragment_to_partial_pdb(fragment, output_path, verify=False)
        with open(output_path, "r") as fd:
            signature: str = fd.read()
        os.unlink(output_path)

        return signature

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Attempts to resolve fingerprint (.pdb file) into standardized PDB IDs
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
            if self.use_alphafind:
                exact_structures: list[str] = alphafind_identity(
                    tmp_path, self.config
                )
            else:
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
        Attempts to resolve fingerprint (.pdb file) into standardized PDB IDs
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
        try:
            if self.use_alphafind:
                similar_structures: list[tuple[str, float]] = alphafind_similarity(
                    tmp_path, self.config
                )
            else:
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
            os.path.join(out_path, f"{out_name}_PROTEIN_RESSTRUCT_structure.pdb"), "w"
        ) as fd:
            fd.write(pdb_content)
        return True
