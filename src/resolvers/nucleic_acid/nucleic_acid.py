import configparser
import logging
import os

import MDAnalysis
from src.utils.file_loader import SimulationFile

from ..protein.blastp_handle import run_blast_identity, run_blast_similarity
from ..base.fragment_resolver import ResolverBase, ResolverKind
from ..base.resolver_database import ResolverCache


class NucleicAcidResolver(ResolverBase):
    RESOLVER_NAME = "nucleic_acid"
    RESOLVER_IDENT = "PDB"

    def __init__(
        self, cache: ResolverCache, configuration: configparser.ConfigParser
    ) -> None:
        super().__init__(
            NucleicAcidResolver.RESOLVER_NAME,
            ResolverKind.PRIMARY,
            NucleicAcidResolver.RESOLVER_IDENT,
            cache,
            configuration,
        )

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the NucleicAcidResolver can be used on the given fragment
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """

        # if all residues have nucleotide name, it is safe to assume molecule represents the nucleic acid
        return all(
            resname.upper() in ('A', 'C', 'T', 'G', 'U') for resname in fragment.resnames
        ) if self.config.getboolean(NucleicAcidResolver.RESOLVER_NAME, 'require_residue_names') else True  # will fail later

    def get_signature(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> str:
        """
        Translates fragment into nucleotide sequence
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: nucleotide sequence of the fragment
        :raises: SignatureGenerationError if fingerprint cannot be created
        """
        unique_resindices = list(set(fragment.resindices))
        unique_resindices.sort()

        residue_names: dict[int, str] = {}
        for resindice, resname in zip(fragment.resindices, fragment.resnames):
            residue_names[resindice] = resname

        return "".join(residue_names[i] for i in unique_resindices)

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Attempts to resolve protein sequence into standardized PDB IDs
        Fetches sequentially identical proteins using BLAST search (see configuration for what is identical)
        :param signature: to be translated into PDB IDs
        :return: list of standardized PDB IDs of the protein sequence
        :raises IdentifierResolutionError: if blast search fails (finding no results is not considered a failure)
        """

        return run_blast_identity(signature, self.config, is_nucleotide=True)

    def try_fetch_identifiers_relaxed_match(self, signature: str) -> list[[str, float]]:
        """
        Attempts to resolve protein sequence into standardized PDB IDs
        Fetches sequentially similar proteins using BLAST search (see configuration for what is similar)
        :param signature: to be translated into PDB IDs
        :return: list of standardized PDB IDs of the protein sequence, with identity % to the query signature
        :raises IdentifierResolutionError: if blast search fails (finding no results is not considered a failure)
        """
        return run_blast_similarity(signature, self.config, is_nucleotide=True)

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
            os.path.join(out_path, f"{out_name}_NUCLEOTIDE_sequence.fasta"), "w"
        ) as fd:
            fd.write(f">{fragment.segments.segids[0]}\n{fasta_seq}\n")
        return True
