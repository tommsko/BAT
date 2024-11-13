import logging
from enum import Enum

import MDAnalysis
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, Mol, EditableMol

from ..base.fragment_resolver import SignatureGenerationError
from ...constants import (
    ATOMTABLE_ELEMENT_EXPECTED_BONDS,
    ATOMTABLE_ELEMENT_EXPECTED_BONDS_NEUTRAL,
    ATOMTABLE_ELEMENT_EXPECTED_BONDS_POSSIBLE,
)


################################################################
# Determination of atomic bonds and formal charges through
# simple chemistry rules


class MolfileBond(Enum):
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4
    UNKNOWN = 8


def molfile_bond_to_rdkit_bond(molfile_bond: int) -> Chem.rdchem.BondType:
    """
    Transforms MolfileBond enum's value to RDKit's bond type
    :param molfile_bond: value from MolfileBond
    :return: semantically same bond of Chem.rdchem.BondType type
    """
    match molfile_bond:
        case 1:
            return Chem.rdchem.BondType.SINGLE
        case 2:
            return Chem.rdchem.BondType.DOUBLE
        case 3:
            return Chem.rdchem.BondType.TRIPLE
        case 4:
            return Chem.rdchem.BondType.AROMATIC
        case 8:
            return Chem.rdchem.BondType.UNSPECIFIED
        case _:
            raise SignatureGenerationError(
                "bond mdanalysis->rdkit translation",
                f"unknown Molfile bond type '{molfile_bond}'",
            )


class DetermineUnknownBonds:
    """
    Helper class that tries to determine bond types of MDAnalysis atoms given basic chemistry rules
    """

    def __init__(
        self, fragment: MDAnalysis.AtomGroup, atom_indices: dict[int, int]
    ) -> None:
        """
        Initializes and tries to determine bond types given;
        :param fragment: AtomGroup representing a molecule
        :param atom_indices: mapping of atom indexes to their position in arrays
        """
        self.bonds: dict[tuple[int, int], MolfileBond] = {
            (key[0], key[1]): MolfileBond.UNKNOWN for key in fragment.bonds.indices
        }

        self.atom_got_bonds: list[int] = [0 for _ in range(len(atom_indices))]
        self.atom_formal_charges: list[int] = [0 for _ in range(len(atom_indices))]
        self.charged_atoms_cnt: int = 0

        self.atom_indices: dict[int, int] = atom_indices
        self.fragment: MDAnalysis.AtomGroup = fragment

        # Single covalent bonds, just a formality as ionic bonds are very rare in organic chemistry
        self._increase_bond(
            bond_from=MolfileBond.UNKNOWN,
            bond_to=MolfileBond.SINGLE,
            new_bonds=1,
            force=True,
        )

        # Double covalent bonds
        self._increase_bond(
            bond_from=MolfileBond.SINGLE,
            bond_to=MolfileBond.DOUBLE,
            new_bonds=1,
            force=False,
        )

        # Triple covalent bonds
        self._increase_bond(
            bond_from=MolfileBond.DOUBLE,
            bond_to=MolfileBond.TRIPLE,
            new_bonds=1,
            force=False,
        )

        # calculate formal charges
        self._calculate_formal_charges()

    def _get_atom_indices(self, atom_idx1: int, atom_idx2: int) -> tuple[int, int]:
        """
        Helper function to translate atom id's into atom indices
        :param atom_idx1: id of atom1
        :param atom_idx2: id of atom2
        :return: tuple of atom1, atom2 indices
        """
        return self.atom_indices[atom_idx1], self.atom_indices[atom_idx2]

    def _get_atom_elements(self, atom_idx1: int, atom_idx2: int) -> tuple[str, str]:
        """
        Helper function to translate atom id's into atom elements
        :param atom_idx1: id of atom1
        :param atom_idx2: id of atom2
        :return: tuple of atom1, atom2 elements
        """
        atom_indice1, atom_indice2 = self._get_atom_indices(atom_idx1, atom_idx2)
        return (
            self.fragment.atoms.elements[atom_indice1].upper(),
            self.fragment.atoms.elements[atom_indice2].upper(),
        )

    def _increase_bond(
        self, bond_from: MolfileBond, bond_to: MolfileBond, new_bonds: int, force: bool
    ) -> None:
        """
        Tries to change bond types to a higher ones, if applicable
        :param bond_from: bondtype that is considered to be changed
        :param bond_to: bondtype that it should be changed to
        :param new_bonds: number of new bonds this creates
        :param force: if True, ignores possible bonding configurations of elements
        :return: None
        """
        for atom_idx1, atom_idx2 in self.bonds:
            atom_i1, atom_i2 = self._get_atom_indices(atom_idx1, atom_idx2)
            atom_elem1, atom_elem2 = self._get_atom_elements(atom_idx1, atom_idx2)
            atom1_got_bonds, atom2_got_bonds = (
                self.atom_got_bonds[atom_i1],
                self.atom_got_bonds[atom_i2],
            )
            atom1_new_bonds, atom2_new_bonds = (
                atom1_got_bonds + new_bonds,
                atom2_got_bonds + new_bonds,
            )

            if self.bonds[(atom_idx1, atom_idx2)] != bond_from:
                continue

            if not force and (
                atom1_new_bonds
                not in ATOMTABLE_ELEMENT_EXPECTED_BONDS[atom_elem1][
                    ATOMTABLE_ELEMENT_EXPECTED_BONDS_POSSIBLE
                ]
                or atom2_new_bonds
                not in ATOMTABLE_ELEMENT_EXPECTED_BONDS[atom_elem2][
                    ATOMTABLE_ELEMENT_EXPECTED_BONDS_POSSIBLE
                ]
            ):
                continue

            self.bonds[(atom_idx1, atom_idx2)] = bond_to
            self.atom_got_bonds[atom_i1] += new_bonds
            self.atom_got_bonds[atom_i2] += new_bonds

    def _calculate_formal_charges(self) -> None:
        """
        Calculates formal charges for all atoms, given pre-calculated bond counts
        :return: None
        """
        for atom_idx in self.fragment.atoms.ids:
            atom_i, _ = self._get_atom_indices(atom_idx, atom_idx)
            atom_elem, _ = self._get_atom_elements(atom_idx, atom_idx)

            expected_bonds: int = ATOMTABLE_ELEMENT_EXPECTED_BONDS[atom_elem][
                ATOMTABLE_ELEMENT_EXPECTED_BONDS_NEUTRAL
            ]
            got_bonds: int = self.atom_got_bonds[atom_i]
            formal_charge: int = got_bonds - expected_bonds
            self.atom_formal_charges[atom_i] = formal_charge
            if formal_charge != 0:
                self.charged_atoms_cnt += 1


################################################################
# Transformation of MDAnalysis molecule into RDKit's molecule
# utilizing positions of atoms (through Molfile intermediate)


def _fragment_to_molfile_atom_block(
    fragment: MDAnalysis.AtomGroup, atom_indices: dict[int, int]
) -> str:
    """
    Generates <<atom block>> for Molfile representation given fragment
    :param fragment: AtomGroup representing a molecule
    :param atom_indices: mapping of atom indexes to their position in arrays
    :return: Molfile's atom block
    """
    result_lines: list[str] = []
    for atom_idx in fragment.ids:
        if not all(
            x < 100000 for x in fragment.atoms.positions[atom_indices[atom_idx]]
        ):
            raise SignatureGenerationError(
                "generating molfile's atom block",
                "atom position too large, will overflow 10 character limit",
            )
        pos_x, pos_y, pos_z = [
            f"{pos:.4f}"[:10].rjust(10, " ")
            for pos in fragment.atoms.positions[atom_indices[atom_idx]]
        ]
        elem = fragment.atoms.elements[atom_indices[atom_idx]].ljust(3, " ")
        result_lines.append(
            f"{pos_x}{pos_y}{pos_z} {elem} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    return "\n".join(result_lines)


def _fragment_to_molfile_bonds(
    fragment: MDAnalysis.AtomGroup,
    atom_indices: dict[int, int],
    determined_bonds: DetermineUnknownBonds,
) -> str:
    """
    Generates <<bonds block>> of MOL representation given AtomGroup.
    :param fragment: AtomGroup representing a molecule
    :param atom_indices: mapping of atom indexes to their position in arrays
    :param determined_bonds: pre-determined bond types and formal charges
    :return: bond block representation of such AtomGroup
    """

    result_lines: list[str] = []
    for bond_indices in fragment.bonds.indices:
        atom_idx1, atom_idx2 = bond_indices
        atom_i1, atom_i2 = atom_indices[atom_idx1], atom_indices[atom_idx2]

        bond_type_str: str = str(
            determined_bonds.bonds[(atom_idx1, atom_idx2)].value
        ).rjust(3, " ")
        atom_idx1_str: str = str(atom_i1 + 1).rjust(3, " ")  # indexed from 1
        atom_idx2_str: str = str(atom_i2 + 1).rjust(3, " ")  # indexed from 1
        result_lines.append(
            f"{atom_idx1_str}{atom_idx2_str}{bond_type_str}  0  0  0  0"
        )
    return "\n".join(result_lines)


def _fragment_to_molfile_properties(
    fragment: MDAnalysis.AtomGroup,
    atom_indices: dict[int, int],
    determined_bonds: DetermineUnknownBonds,
) -> str:
    """
    Generates <<properties block>> of MOL representation given AtomGroup.
    :param fragment: AtomGroup representing a molecule
    :param atom_indices: mapping of atom indexes to their position in arrays
    :param determined_bonds: pre-determined bond types and atom formal charges
    :return: properties block representation of such AtomGroup
    """

    result_lines: list[str] = []
    for atom_idx in fragment.atoms.ids:
        atom_i: int = atom_indices[atom_idx]
        atom_charge: int = determined_bonds.atom_formal_charges[atom_i]
        if atom_charge == 0:
            continue

        atom_idx_str: str = str(atom_i + 1).rjust(4, " ")  # indexed from 1
        atom_charge_str: str = str(atom_charge).rjust(4, " ")
        result_lines.append(f"M  CHG  1{atom_idx_str}{atom_charge_str}")
    return "\n".join(result_lines)


def _fragment_to_molfile_counts(
    fragment: MDAnalysis.AtomGroup, determined_bonds: DetermineUnknownBonds
) -> str:
    """
    Generates <<counts block>> of MOL representation given AtomGroup
    :param fragment: AtomGroup representing a molecule
    :param determined_bonds: pre-determined bond types and formal charges
    :return: counts block representation of such AtomGroup
    """
    atoms_cnt: str = str(len(fragment.atoms.ids)).rjust(3, " ")
    bonds_cnt: str = str(len(fragment.bonds.indices)).rjust(3, " ")
    properties_cnt: str = str(determined_bonds.charged_atoms_cnt).rjust(3, " ")
    return f"{atoms_cnt}{bonds_cnt}  0     0  0  0  0  0  0{properties_cnt} V2000"


def _fragment_to_molfile(fragment: MDAnalysis.AtomGroup) -> str:
    """
    Converts MDAnalysis molecule into textual MDL Molfile format
    :param fragment: AtomGroup representing a molecule
    :return: Molfile representation of such AtomGroup
    """

    atom_indices: dict[int, int] = {
        atom_idx: i for i, atom_idx in enumerate(fragment.atoms.ids)
    }
    determined_bonds: DetermineUnknownBonds = DetermineUnknownBonds(
        fragment, atom_indices
    )

    return f"""000
Exported_molecule

{_fragment_to_molfile_counts(fragment, determined_bonds)}
{_fragment_to_molfile_atom_block(fragment, atom_indices)}
{_fragment_to_molfile_bonds(fragment, atom_indices, determined_bonds)}
{_fragment_to_molfile_properties(fragment, atom_indices, determined_bonds)}
M  END
"""


def _fragment_to_mol_through_molfile(fragment: MDAnalysis.AtomGroup) -> Mol:
    """
    Converts MDAnalysis' molecule into RDKit's molecule using Molfile intermediate
    :param fragment: AtomGroup representing a molecule
    :return: RDKit's representation of the same molecule, sanitized
    """

    molfile: str = _fragment_to_molfile(fragment)
    return Chem.MolFromMolBlock(molfile, sanitize=True, strictParsing=False)


################################################################
# Transformation of MDAnalysis molecule into RDKit's molecule
# directly (building the molecule)


def _fragment_to_rdkit_molecules(
    fragment: MDAnalysis.AtomGroup,
    atom_indices: dict[int, int],
    determined_bonds: DetermineUnknownBonds,
    molecule: Mol,
) -> dict[int, int]:
    """
    Adds all MDAnalysis atoms to the RDKit's molecule
    :param fragment: AtomGroup representing a molecule
    :param atom_indices: mapping of atom indexes to their position in arrays
    :param determined_bonds: pre-determined bond types and atom formal charges
    :param molecule: RDKit's molecule to add atoms to
    :return: mapping of MDAnalysis' atom indexes to RDKit's atom indexes
    """

    idx_atom_to_rdkit: dict[int, int] = {}
    for atom_idx in fragment.atoms.ids:
        atom_elem = fragment.atoms.elements[atom_indices[atom_idx]]
        atom = Chem.Atom(atom_elem)
        atom.SetFormalCharge(
            determined_bonds.atom_formal_charges[atom_indices[atom_idx]]
        )
        atom.SetNoImplicit(True)
        idx_atom_to_rdkit[atom_idx] = molecule.AddAtom(atom)
    return idx_atom_to_rdkit


def _fragment_to_rdkit_bonds(
    fragment: MDAnalysis.AtomGroup,
    determined_bonds: DetermineUnknownBonds,
    idx_mdanalyis_to_rdkit: dict[int, int],
    molecule: Mol,
) -> None:
    """
    Adds all MDAnalysis bonds to RDKit's molecule
    :param fragment: AtomGroup representing a molecule
    :param determined_bonds: pre-determined bond types and atom formal charges
    :param idx_mdanalyis_to_rdkit: mapping of MDAnalysis' atom indexes to RDKit's atom indexes
    :param molecule: RDKit's molecule to add bonds to
    :return: None
    """
    for atom_idx1, atom_idx2 in fragment.bonds.indices:
        atom_idx1_rdkit, atom_idx2_rdkit = (
            idx_mdanalyis_to_rdkit[atom_idx1],
            idx_mdanalyis_to_rdkit[atom_idx2],
        )
        bond_type = molfile_bond_to_rdkit_bond(
            determined_bonds.bonds[(atom_idx1, atom_idx2)].value
        )
        molecule.AddBond(atom_idx1_rdkit, atom_idx2_rdkit, bond_type)


def _fragment_to_mol_directly(fragment: MDAnalysis.AtomGroup) -> Mol:
    """
    Converts MDAnalysis' molecule into RDKit's molecule directly
    :param fragment: AtomGroup representing a molecule
    :return: RDKit's representation of the same molecule, sanitized
    """

    atom_indices: dict[int, int] = {
        atom_idx: i for i, atom_idx in enumerate(fragment.atoms.ids)
    }
    determined_bonds: DetermineUnknownBonds = DetermineUnknownBonds(
        fragment, atom_indices
    )

    molecule_edit: EditableMol = Chem.EditableMol(Chem.Mol())
    idx_mdanalyis_to_rdkit: dict[int, int] = _fragment_to_rdkit_molecules(
        fragment, atom_indices, determined_bonds, molecule_edit
    )
    _fragment_to_rdkit_bonds(
        fragment, determined_bonds, idx_mdanalyis_to_rdkit, molecule_edit
    )
    molecule: Mol = molecule_edit.GetMol()
    molecule = Chem.RemoveHs(molecule)
    Chem.SanitizeMol(molecule)

    return molecule


################################################################
# Main callable methods


def _fix_bonds_rdkit(molecule: Mol, max_attempts: int = 10) -> None:
    """
    Attempts to call DetermineBondOrders, adapting to total charge provided in call's exception
    :param max_attempts: maximum attempts to call DetermineBondOrders
    :param molecule: RDKit's molecule
    :return: None (throws on failure)
    """

    attempt: int = 0
    charge: int = 0
    while attempt < max_attempts:
        attempt += 1
        try:
            rdDetermineBonds.DetermineBonds(
                molecule, allowChargedFragments=False, embedChiral=True, charge=charge
            )
        except ValueError as e:
            charge = int(str(e).split("input (")[1].split(");")[0])
        else:
            return
    raise SignatureGenerationError(
        "correcting bonds using RDKit",
        "unable to determine DetermineBonds() charge parameter within max attempts!",
    )


def fragment_to_mol(fragment: MDAnalysis.AtomGroup, try_fix_bonds: bool = False) -> Mol:
    """
    Converts MDAnalysis' molecule into RDKit's molecule
    :param fragment: AtomGroup representing a molecule
    :param try_fix_bonds: tries to utilize atomic positions (from .gro file) to determine bonds
    :return: RDKit's representation of the same molecule, sanitized
    """

    if not try_fix_bonds:
        molecule: Mol = _fragment_to_mol_directly(fragment)
    else:
        if not hasattr(fragment.atoms, "positions"):
            raise SignatureGenerationError(
                "initializing fingerprint generation",
                "try_fix_bonds = True requires atomic positions to be present",
            )
        molecule: Mol = _fragment_to_mol_through_molfile(fragment)
        _fix_bonds_rdkit(molecule)
    return molecule


def fragment_to_smiles(
    fragment: MDAnalysis.AtomGroup, try_fix_bonds: bool = False
) -> str:
    """
    Converts MDAnalysis' molecule into SMILES format (through RDKit intermediate)
    :param fragment: AtomGroup representing a molecule
    :param try_fix_bonds: tries to utilize atomic positions (from .gro file) to determine bonds
    :return: SMILES representation of the same molecule
    """
    return Chem.MolToSmiles(fragment_to_mol(fragment, try_fix_bonds))
