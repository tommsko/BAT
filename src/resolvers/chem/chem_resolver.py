import configparser
import logging
import os.path

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import MDAnalysis

from .chembl_handle import chembl_fetch_similarity_match
from .mdanalysis_to_rdkit import fragment_to_smiles, _fragment_to_molfile
from .pubchem_handle import pubchem_fetch_identity_match, pubchem_fetch_similarity_match
from ..base.fragment_resolver import (
    ResolverBase, ResolverKind
)
from ..base.resolver_database import ResolverCache
from src.utils.file_loader import SimulationFile
from ...utils.repair_simulation import fix_missing_elements

log = logging.getLogger("resolver")


def _get_pubchem_attempt_thresholds(config: configparser.ConfigParser) -> list[int]:
    """
    Fetches list of similarity thresholds to try in PUBCHEM search
    If there are no thresholds set in configuration, defaults to [95]
    :param config: main configuration
    :return: list of similarity thresholds to try
    """
    thresholds: list[int] = []
    for _, threshold in config.items("PUBCHEM_SMILES_LOOKUP__SIMILARITY_THRESHOLD_ATTEMPTS"):
        thresholds.append(int(threshold))
    return thresholds if thresholds else [95]


def _get_chembl_attempt_threshold(config: configparser.ConfigParser) -> int:
    """
    Fetches similarity threshold for CHEMBL search
    If there is no threshold set in configuration, defaults to [40]
    :param config: main configuration
    :return: list of similarity thresholds to try
    """
    return config.getint("CHEMBL_SMILES_LOOKUP__SIMILARITY_THRESHOLD", "threshold", fallback=40)


class ChemResolver(ResolverBase):
    """
    ChemResolver looks at the fragment as a small self-contained molecule
    """
    RESOLVER_NAME, RESOLVER_NAME_APPROX = "chem", "chem_aprox"
    RESOLVER_IDENT = "InChI"

    def __init__(
        self,
        cache: ResolverCache,
        configuration: configparser.ConfigParser,
        try_fix_elements: bool = False,
        try_approximate_bonds: bool = False,
        use_rdkit_bond_determination: bool = False
    ) -> None:
        """
        Initializes the ChemResolver
        :param cache: resolver cache
        :param configuration: main configuration
        :param try_fix_elements: if element information is not present, try to figure it out using atomic masses
        :param try_approximate_bonds: if positions of atoms are present, it will try to approximate connectivity (bonds)
        :param use_rdkit_bond_determination: uses internal bond order determination by rdkit
        """
        super().__init__(
            ChemResolver.RESOLVER_NAME if not try_approximate_bonds else ChemResolver.RESOLVER_NAME_APPROX,
            ResolverKind.PRIMARY,
            ChemResolver.RESOLVER_IDENT,
            cache,
            configuration,
        )

        self.try_fix_missing_elements: bool = try_fix_elements
        self.try_fix_bonds: bool = try_approximate_bonds
        self.use_rdkit_bond_determination: bool = use_rdkit_bond_determination

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the ChemResolver can be used on the given fragment
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """

        fragment_name: str = fragment.segids[0]

        # this resolver resolves chemical molecule: it needs atoms, their elements and bonds inbetween them
        # to safely determine the molecule and search corresponding databases

        # firstly, atomic elements *MUST* be present or resolved (no way around this to get chemical molecule)
        if not simulation.supports_atomic_elements(fragment_name):
            if (
                simulation.supports_atomic_masses(fragment_name)
                and self.try_fix_missing_elements
                and fix_missing_elements(simulation, fragment)
            ):
                ...  # ok!
            else:
                return False

        # secondly, bonds between atoms must be present or resolved (no way around this either)
        if not self.try_fix_bonds:
            atom_indices: dict[int, int] = {
                atom_idx: i for i, atom_idx in enumerate(fragment.atoms.ids)
            }
            is_atom_bonded: list[bool] = [False for _ in range(len(atom_indices))]
            for atom_idx1, atom_idx2 in fragment.bonds.indices:
                is_atom_bonded[atom_indices[atom_idx1]] = True
                is_atom_bonded[atom_indices[atom_idx2]] = True
            if not all(is_atom_bonded):
                return False
        else:
            if not simulation.supports_atomic_positions(fragment_name):
                return False

        return True

    def get_signature(self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup) -> str:
        """
        Translates fragment into SMILES
        :param simulation: of which is the fragment part of
        :param fragment: AtomGroup representing a molecule
        :return: SMILES sequence of the fragment
        :raises: SignatureGenerationError if fingerprint cannot be created
        """
        return fragment_to_smiles(fragment,
                                  try_fix_bonds=self.try_fix_bonds,
                                  try_fix_bonds_alt=self.use_rdkit_bond_determination)

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Attempts to resolve fingerprint (SMILES) into standardized InChI Key
        Reported molecules (in the form of standardized InChI key) have to have identity match to the fingerprint
        :param signature: to be translated into standardized InChI keys
        :return: list of standardized InChI keys of the molecule
        :raises IdentifierResolutionError: if any search fails (finding no results is NOT considered a failure)
        """
        if self.config.getboolean(self.resolver_name, 'allow_pubchem_fetch', fallback=False):
            pubchem_idents: list[str] = pubchem_fetch_identity_match(signature, config=self.config)
        else:
            pubchem_idents: list[str] = []

        if self.config.getboolean(self.resolver_name, 'allow_chembl_fetch', fallback=False):
            chembl_idents: list[str] = chembl_fetch_similarity_match(
                signature, threshold=100, config=self.config
            )
        else:
            chembl_idents: list[str] = []

        return list(set(pubchem_idents).union(set(chembl_idents)))

    def try_fetch_identifiers_relaxed_match(
        self, signature: str
    ) -> list[tuple[str, float]]:
        """
        Attempts to resolve fingerprint (SMILES) into standardized InChI key
        Reported molecules (in the form of standardized InChI key) have to have similarity match to the signature
        :param signature: to be translated into standardized InChI keys
        :return: list of standardized InChI keys of the molecule and their scores (similarity %)
        :raises IdentifierResolutionError: if any search fails (finding no results is not considered a failure)
        """
        if self.config.getboolean(self.resolver_name, 'allow_pubchem_fetch', fallback=False):
            pubchem_idents: list[tuple[str, float]] = pubchem_fetch_similarity_match(
                signature, config=self.config, try_thresholds=_get_pubchem_attempt_thresholds(self.config)
            )
        else:
            pubchem_idents: list[tuple[str, float]] = []

        if self.config.getboolean(self.resolver_name, 'allow_chembl_fetch', fallback=False):
            chembl_idents: list[[str, float]] = chembl_fetch_similarity_match(
                signature, threshold=_get_chembl_attempt_threshold(self.config), config=self.config, include_similarity=True
            )
        else:
            chembl_idents: list[[str, float]] = []

        return pubchem_idents + chembl_idents

    def generate_debug_data(
        self,
        simulation: SimulationFile,
        fragment: MDAnalysis.AtomGroup,
        out_path: str,
        out_name: str,
    ) -> bool:
        """
        Generates (chemical) image of the molecule as well as molfile, if applicable
        :param out_path: directory where to save generated representations
        :param out_name: basename of the file to save representations to (will be suffixed with kind of export)
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True if export went without problems
        """
        mol: Chem.Mol = Chem.MolFromSmiles(self.get_signature(simulation, fragment))
        AllChem.Compute2DCoords(mol)
        Draw.MolToImage(mol).save(os.path.join(out_path, f"{out_name}_CHEM_molecule.jpg"))

        if simulation.supports_atomic_positions(fragment.segments.segids[0]):
            with open(os.path.join(out_path, f"{out_name}_CHEM_molecule.mol"), 'w') as fd:
                fd.write(_fragment_to_molfile(fragment))
        else:
            mol2: Chem.Mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol2)
            AllChem.MMFFOptimizeMolecule(mol2)
            Chem.MolToMolFile(mol2, os.path.join(out_path, f"{out_name}_CHEM_molecule.mol"))
        return True
