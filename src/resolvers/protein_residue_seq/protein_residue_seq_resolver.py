import configparser

from ..base.fragment_resolver import (
    ResolverBase,
)
from ..base.resolver_database import ResolverCache
from ..protein.mdanalysis_to_protseq import fragment_to_protein_sequence_residue_based
from ..base.fragment_resolver import ResolverKind
from ..protein.blastp_handle import run_blast_identity, run_blast_similarity
import MDAnalysis

from src.utils.file_loader import SimulationFile
from ...constants import KNOWN_AMINO_ACIDS


class ProteinResidueSequenceResolver(ResolverBase):
    """
    ProteinResidueSequenceResolver assumes that molecule represents a single protein and has been built from
    a single PDB, utilizing alternate (but possibly error-prone) method to find protein amino-acid sequence
    """
    RESOLVER_NAME = "protein-residue-sequence"
    RESOLVER_IDENT = "PDB"

    def __init__(self, cache: ResolverCache, configuration: configparser.ConfigParser) -> None:
        super().__init__(
            ProteinResidueSequenceResolver.RESOLVER_NAME, ResolverKind.PRIMARY, ProteinResidueSequenceResolver.RESOLVER_IDENT, cache, configuration
        )

    def _is_applicable(self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup) -> bool:
        """
        Decides whether the ProteinResolver can be used on the given fragment
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True if resolver can be used on a given fragment
        """

        # if all residues have amino acid name, it is safe to assume molecule represents the protein
        return all(
            resname.upper() in KNOWN_AMINO_ACIDS for resname in fragment.resnames
        )

    def get_signature(self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Translates fragment into protein sequence
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: amino-acid sequence of the fragment
        :raises: SignatureGenerationError if fingerprint cannot be created
        """
        return fragment_to_protein_sequence_residue_based(fragment)

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Attempts to resolve protein sequence into standardized PDB IDs
        Fetches sequentially identical proteins using BLAST search (see configuration for what is identical)
        :param signature: to be translated into PDB IDs
        :return: list of standardized PDB IDs of the protein sequence
        :raises IdentifierResolutionError: if blast search fails (finding no results is not considered a failure)
        """

        return run_blast_identity(signature, self.config)

    def try_fetch_identifiers_relaxed_match(self, signature: str) -> list[[str, float]]:
        """
        Attempts to resolve protein sequence into standardized PDB IDs
        Fetches sequentially similar proteins using BLAST search (see configuration for what is similar)
        :param signature: to be translated into PDB IDs
        :return: list of standardized PDB IDs of the protein sequence, with identity % to the query signature
        :raises IdentifierResolutionError: if blast search fails (finding no results is not considered a failure)
        """
        return run_blast_similarity(signature, self.config)

    def generate_debug_data(
        self,
        simulation: SimulationFile,
        fragment: MDAnalysis.AtomGroup,
        out_path: str,
        out_name: str,
    ) -> bool:
        """
        Generates fasta sequence of the protein
        :param out_path: directory where to save generated representations
        :param out_name: basename of the file to save representations to (will be suffixed with kind of export)
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True if export went without problems
        """
        fasta_seq: str = self.get_signature(simulation, fragment)
        with open(
            os.path.join(out_path, f"{out_name}_PROTEIN_RESSEQ_sequence.fasta"), "w"
        ) as fd:
            fd.write(f">{fragment.segments.segids[0]}\n{fasta_seq}\n")
        return True