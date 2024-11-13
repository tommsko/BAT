import configparser
import json
import logging
import os

import MDAnalysis
from src.utils.file_loader import SimulationFile

from .mdanalysis_to_networkx import fragment_to_bond_maps
from .mdanalysis_to_signature import hash_fragment
from ..base.fragment_resolver import ResolverBase, ResolverKind
from ..base.resolver_database import ResolverCache
from ..protein.blastp_handle import _mk_random_filename

log = logging.getLogger("resolver")


class FingerprintResolver(ResolverBase):
    """
    FingerprintResolver looks at the molecule as a graph, with atoms as nodes and bonds as links between them.
    By not utilizing any other information, it cannot do any fetch-based resolutions. However, it can identify
    the exactly same molecule across several simulations (because variability in biomolecules is enormous).
    """

    RESOLVER_NAME = "fingerprint"
    RESOLVER_IDENT = "hash"

    def __init__(
        self,
        cache: ResolverCache,
        configuration: configparser.ConfigParser,
    ) -> None:
        super().__init__(
            FingerprintResolver.RESOLVER_NAME,
            ResolverKind.FALLBACK,
            FingerprintResolver.RESOLVER_IDENT,
            cache,
            configuration,
        )
        self.fingerprints_directory: str = configuration.get(
            self.resolver_name, "unknown_fingerprints_export_dir"
        )
        if not os.path.exists(self.fingerprints_directory):
            os.makedirs(self.fingerprints_directory)

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the FingerprintResolver can be used on the given fragment
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """

        return True  # always can be (atoms and bonds are always present)

    def get_signature(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> str:
        """
        Translates fragment into deterministic hash of its atoms and bonds.
        Also saves the hash for further identification
        :param simulation: of which is the fragment part of
        :param fragment: AtomGroup representing a molecule
        :return: hash of the fragment
        :raises: SignatureGenerationError if fingerprint cannot be created
        """

        fragment_name: str = fragment.segids[0]
        fragment_hash: str = hash_fragment(fragment)

        save_path: str = _mk_random_filename(
            "txt", self.fingerprints_directory, name_len=30
        )
        with open(save_path, "w") as fd:
            json.dump(
                {
                    "simulation": simulation.paths[0],
                    "fragment_name": fragment_name,
                    "hash": fragment_hash,
                },
                fd,
            )
        log.debug(f"Unknown fingerprint encountered, saved to '{save_path}'")

        return fragment_hash

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        !! UNUSED !! Only cached (manually curated hashes) can be returned as identity match
        :param signature: !! UNUSED !!
        :return: empty list
        """
        return []

    def try_fetch_identifiers_relaxed_match(self, signature: str) -> list[[str, float]]:
        """
        !! UNUSED !! Only cached (manually curated hashes) can be returned as similarity match
        :param signature: !! UNUSED !!
        :return: empty list
        """
        return []

    def generate_debug_data(
        self,
        simulation: SimulationFile,
        fragment: MDAnalysis.AtomGroup,
        out_path: str,
        out_name: str,
    ) -> bool:
        """
        Generates bonds map of atoms and residues in the simulation
        :param out_path: directory where to save generated representations
        :param out_name: basename of the file to save representations to (will be suffixed with kind of export)
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True if export went without problems
        """
        log.debug(
            f"Generating bond maps for {simulation.simulation_name}/{fragment.segments.segids[0]}"
        )

        path_atoms: str = os.path.join(
            out_path, f"{out_name}_FINGERPRINT_bond_map_atoms.html"
        )
        path_residues: str = os.path.join(
            out_path, f"{out_name}_FINGERPRINT_bond_map_residues.html"
        )
        try:
            fragment_to_bond_maps(simulation, fragment, path_atoms, path_residues)
        except Exception as exc:
            log.error(f"Generation of bond maps failed: {exc}")
            return False
        return True
