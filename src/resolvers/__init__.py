import configparser

from src.resolvers.atom.atom_resolver import AtomResolver
from src.resolvers.base.fragment_resolver import ResolverBase
from src.resolvers.base.resolver_database import ResolverCache
from src.resolvers.chem.chem_resolver import ChemResolver
from src.resolvers.chem_ligand_name.chem_ligand_name import ChemLigandNameResolver
from src.resolvers.fingerprint.fingerprint_resolver import FingerprintResolver
from src.resolvers.protein.protein_resolver import ProteinResolver
from src.resolvers.protein_residue_seq.protein_residue_seq_resolver import (
    ProteinResidueSequenceResolver,
)
from src.resolvers.protein_residue_struct.protein_residue_structure_resolver import (
    ProteinResidueStructureResolver,
)
from src.resolvers.protein_struct.protein_structure_resolver import (
    ProteinStructureResolver,
)


def get_resolvers(
    cache: ResolverCache, configuration: configparser.ConfigParser
) -> list[ResolverBase]:
    """
    Provides all resolvers that can be tried for any molecule
    :param cache: cache for resolver initialization
    :param configuration: main configuration
    :return: list of all resolvers
    """

    return [
        ProteinResolver(cache, configuration),
        ProteinStructureResolver(
            cache, configuration, try_fix_elements=True, try_fix_names=True
        ),
        ProteinResidueSequenceResolver(cache, configuration),
        ProteinResidueStructureResolver(cache, configuration),
        # protein, but alphafind (.pdb export)
        AtomResolver(cache, configuration, try_fix_elements=True),
        AtomResolver(
            cache,
            configuration,
            try_fix_elements=True,
            allow_undetermined_elements=True,
        ),
        ChemResolver(cache, configuration, try_fix_elements=True),
        ChemResolver(
            cache, configuration, try_fix_elements=True, try_approximate_bonds=True
        ),
        ChemLigandNameResolver(cache, configuration),
        FingerprintResolver(
            cache,
            configuration,
        ),
    ]
