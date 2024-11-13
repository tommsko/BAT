import itertools
import string

import MDAnalysis

from src.constants import (
    ATOMTABLE_ELEMENT_MASS_CONFIDENCE_THRESHOLD,
    ATOMTABLE_ELEMENTS_MASSES,
)
from .file_loader import SimulationFile


def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    # https://stackoverflow.com/questions/5595425/how-to-compare-floats-for-almost-equality-in-python
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def fix_missing_elements(
    simulation: SimulationFile,
    for_fragment: MDAnalysis.AtomGroup,
    replace_with_mass: bool = False,
) -> bool:
    """
    Attempts fix missing atomic element information in the simulation by determining atomic elements using atomic masses
    Fixing is done only on the AtomGroup provided in the arguments (not all the atoms!!)
    :param simulation: containing the fragment
    :param for_fragment: in which to fix missing atomic elements
    :param replace_with_mass: if resolution fails, replaces element names with its atomic mass
                              (in almost all cases, False)
    :return: True if all elements in fragment were distinguished confidently
    :raises AssertionError: if simulation is not loaded before
    """
    new_elements: list[str] = ["" for _ in range(for_fragment.n_atoms)]

    all_successful: bool = True
    for i in range(for_fragment.n_atoms):
        atom_mass: float = for_fragment.masses[i]
        for element_name, element_mass in ATOMTABLE_ELEMENTS_MASSES.items():
            if is_close(
                element_mass,
                atom_mass,
                rel_tol=0,
                abs_tol=ATOMTABLE_ELEMENT_MASS_CONFIDENCE_THRESHOLD,
            ):
                new_elements[i] = element_name
                break
        else:
            if replace_with_mass:
                new_elements[i] = str(for_fragment.atoms.masses[i])
            all_successful = False

    universe: MDAnalysis.Universe = simulation.get_simulation()
    if not hasattr(simulation, "elements"):
        universe.add_TopologyAttr(
            "elements", ["" for _ in range(simulation.get_simulation().atoms.n_atoms)]
        )

    for_fragment.elements = new_elements
    return all_successful


def fix_missing_names(
    simulation: SimulationFile, for_fragment: MDAnalysis.AtomGroup
) -> bool:
    """
    Attempts fix missing atomic element information in the simulation by assigning unique names to atoms
    Fixing is done only on the AtomGroup provided in the arguments (not all the atoms!!)
    :param simulation: containing the fragment
    :param for_fragment: in which to fix missing atomic names
    :return: always True
    :raises AssertionError: if simulation is not loaded before
    """

    universe: MDAnalysis.Universe = simulation.get_simulation()
    if not hasattr(universe.atoms, "names"):
        universe.atoms.add_TopologyAttr(
            "names", ["" for _ in range(simulation.get_simulation().atoms.n_atoms)]
        )

    symbols: str = string.ascii_letters + string.digits
    new_names: list[str] = []
    for item in itertools.product(symbols, repeat=2):
        if len(new_names) == for_fragment.atoms.n_atoms:
            break
        new_names.append("".join(item))

    for_fragment.atoms.names = new_names
    return True


def fix_molecule_use_alternate_bonds_information(
    for_fragment: MDAnalysis.AtomGroup,
) -> MDAnalysis.AtomGroup:
    """
    Attempts to extract a single molecule, given a fragment suspected to be composed of multiple molecules
    (fragment should be a single molecule already, but this has happened in older GROMACS versions)
    Beware that this function DOES NOT CHECK whether the precondition is correct and thus might produce insensible
    results
    :param for_fragment: from which to extract the molecule
    :return: MDAnalysis.AtomGroup of alternate molecule using different atom information, in an alternate universe
    :raises AssertionError: if simulation is not loaded before
    """

    # molecule might be a compound of several sub-molecules bonded with hydrogen bonds, which are still
    # reported as bond information. Alternate source comes from naming of atoms, which seems to be unique

    atoms_types_to_indices: dict[str, int] = {
        name: i for i, name in enumerate(for_fragment.types)
    }

    def extract_bond_name(bond_info: str | tuple) -> tuple[str, str]:
        """
        Transforms string bond information '('Nap', 'SC5p')' into bond atom names (Nap, SC5p).
        Fallback if bond information is already provided as a tuple
        :param bond_info: string bond information
        :return: bond names
        """
        if isinstance(bond_info, str):
            atom1, atom2 = bond_info.split(",")
            return (
                atom1.replace("'", "").replace("(", "").strip(),
                atom2.replace("'", "").replace(")", "").strip(),
            )
        return bond_info

    # topDict contains 'unique' bond information, between named atoms
    bonds_alternate_str: list[tuple[str, str]] = [
        extract_bond_name(bond_info) for bond_info in for_fragment.bonds.topDict
    ]
    bonds_alternate_indices: list[tuple[int, int]] = [
        (atoms_types_to_indices[atom1], atoms_types_to_indices[atom2])
        for atom1, atom2 in bonds_alternate_str
    ]

    # we will create subset of 'uniquely bonded' atoms
    atom_filtered_indices: list[int] = list(
        set(atom_i for bonds in bonds_alternate_indices for atom_i in bonds)
    )
    atom_filtered_indices.sort()

    # every atom will have its own residue, keeping original residue name
    alternate_atom_resindices: list[int] = list(range(len(atom_filtered_indices)))
    alternate_atom_resnames: list[str] = [
        for_fragment.atoms.resnames[i] for i in atom_filtered_indices
    ]

    # however, atoms are now indexed from 0, so we need to update bonding information appropriately
    bonds_alternate_indices_zerobased: dict[int, int] = {
        idx: i for i, idx in enumerate(atom_filtered_indices)
    }

    mol = MDAnalysis.Universe.empty(
        len(atom_filtered_indices),
        n_residues=len(alternate_atom_resindices),
        atom_resindex=alternate_atom_resindices,
        residue_segindex=[0 for _ in range(len(atom_filtered_indices))],
        trajectory=True,
    )

    mol.add_TopologyAttr(
        "name", [for_fragment.types[i] for i in atom_filtered_indices]
    )  # types as the new names
    mol.add_TopologyAttr(
        "type", [for_fragment.types[i][0] for i in atom_filtered_indices]
    )  # does not really matter
    mol.add_TopologyAttr("masses", [0 for _ in atom_filtered_indices])
    mol.add_TopologyAttr("elements", ["" for _ in atom_filtered_indices])
    mol.add_TopologyAttr("resname", alternate_atom_resnames)
    mol.add_TopologyAttr("id", list(range(len(atom_filtered_indices))))
    mol.add_TopologyAttr("resid", alternate_atom_resindices)
    mol.add_TopologyAttr("segid", ["ALT"])
    mol.add_TopologyAttr(
        "bonds",
        [
            (
                bonds_alternate_indices_zerobased[a1],
                bonds_alternate_indices_zerobased[a2],
            )
            for a1, a2 in bonds_alternate_indices
        ],
    )
    if hasattr(for_fragment.atoms, "positions"):
        positions_filtered: list[tuple[float, float, float]] = [
            for_fragment.atoms.positions[atom_i] for atom_i in atom_filtered_indices
        ]
        mol.atoms.positions = positions_filtered
    return mol.select_atoms("segid ALT")
