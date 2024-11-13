import warnings
from collections import defaultdict

import MDAnalysis
import numpy
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


def _pdb_atom_line(
    atom_idx: int,
    atom_name: str,
    atom_elem: str,
    residue_name: str,
    residue_idx: int,
    pos_x: numpy.float32,
    pos_y: numpy.float32,
    pos_z: numpy.float32,
) -> str:
    """
    Creates well-formatted ATOM line of .pdb file
    :param atom_idx: atom ID
    :param atom_name: atom name (must be unique!)
    :param atom_elem: atom element
    :param residue_name: residue name
    :param residue_idx: residue ID
    :param pos_x: atom X position
    :param pos_y: atom Y position
    :param pos_z: atom Z position
    :return: ATOM line with values provided
    """
    # https://cupnet.net/pdb-format/
    return "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
        "ATOM",
        atom_idx,
        atom_name,
        "",
        residue_name.upper(),
        "A",
        residue_idx,
        "",
        pos_x,
        pos_y,
        pos_z,
        1,
        0,
        atom_elem,
        "",
    )


def _pdb_bonds_line(atom_idx: int, other_atoms: set[int]) -> str:
    """
    Creates well-formatted CONECT line of .pdb file
    :param atom_idx: atom ID of which bonds are described
    :param other_atoms: IDs of other atoms bonded to this atom
    :return: CONECT line with values provided
    """
    # https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html
    conect_fmts = [
        "{:6s}{:5d}{:5d}",
        "{:6s}{:5d}{:5d}{:5d}",
        "{:6s}{:5d}{:5d}{:5d}{:5d}",
        "{:6s}{:5d}{:5d}{:5d}{:5d}{:5d}",
    ]
    return conect_fmts[len(other_atoms) - 1].format("CONECT", atom_idx, *other_atoms)


def fragment_to_pdb(
    fragment: MDAnalysis.AtomGroup, out_path: str, verify: bool = False
) -> None:
    """
    Converts a fragment into PDB file format
    :param fragment: AtomGroup representing a protein
    :param out_path: where to write .pdb file
    :param verify: if set to True, try to parse created output (will throw on failure)
    :return: None
    :raises Exception: if atom names, atom elements, or atom positions are not present
    """

    with open(out_path, "w") as fd:
        fd.write("HEADER    EXPORT\n")
        for atom_i in range(fragment.n_atoms):
            fd.write(
                _pdb_atom_line(
                    fragment.atoms.ids[atom_i],
                    fragment.atoms.names[atom_i],
                    fragment.atoms.elements[atom_i],
                    fragment.atoms.resnames[atom_i],
                    fragment.atoms.resindices[atom_i],
                    fragment.atoms.positions[atom_i][0],
                    fragment.atoms.positions[atom_i][1],
                    fragment.atoms.positions[atom_i][2],
                )
                + "\n"
            )

        conns: dict[int, set] = defaultdict(set)
        for atom_idx1, atom_idx2 in fragment.atoms.bonds.indices:
            conns[atom_idx1].add(atom_idx2)
            conns[atom_idx2].add(atom_idx1)
        for atom_idx in conns:
            fd.write(_pdb_bonds_line(atom_idx, conns[atom_idx]) + "\n")

        fd.write("END\n")

    if verify:
        warnings.filterwarnings("error", category=PDBConstructionWarning)
        p = PDBParser()
        p.get_structure("test", out_path)
        warnings.resetwarnings()


def fragment_to_partial_pdb(
    fragment: MDAnalysis.AtomGroup, out_path: str, verify: bool = False
) -> None:
    """
    Converts a fragment into very skinny PDB file format (without element names, without bonds)
    :param fragment: AtomGroup representing a protein
    :param out_path: where to write .pdb file
    :param verify: if set to True, try to parse created output (will throw on failure)
    :return: None
    :raises Exception: if atom names, atom elements, or atom positions are not present
    """

    atom_idx_to_indices: dict[int, int] = {idx: i for i, idx in enumerate(fragment.atoms.ids)}
    residue_idx_to_indices: dict[int, int] = {}
    i = 0
    for resididx in fragment.atoms.resindices:
        if resididx not in residue_idx_to_indices:
            residue_idx_to_indices[resididx] = i
            i += 1

    with open(out_path, "w") as fd:
        fd.write("HEADER    EXPORT\n")
        for atom_i in range(fragment.n_atoms):
            fd.write(
                _pdb_atom_line(
                    atom_idx_to_indices[fragment.atoms.ids[atom_i]],
                    fragment.atoms.types[atom_i],
                    "",
                    fragment.atoms.resnames[atom_i],
                    residue_idx_to_indices[fragment.atoms.resindices[atom_i]],
                    fragment.atoms.positions[atom_i][0],
                    fragment.atoms.positions[atom_i][1],
                    fragment.atoms.positions[atom_i][2],
                )
                + "\n"
            )
        conns: dict[int, set] = defaultdict(set)
        for atom_idx1, atom_idx2 in fragment.atoms.bonds.indices:
            conns[atom_idx1].add(atom_idx2)
            conns[atom_idx2].add(atom_idx1)
        for atom_idx in conns:
            if len(conns[atom_idx]) <= 4:
                fd.write(_pdb_bonds_line(atom_idx_to_indices[atom_idx], {atom_idx_to_indices[x] for x in conns[atom_idx]}) + "\n")

        fd.write("END\n")

    if verify:
        warnings.filterwarnings("error", category=PDBConstructionWarning)
        p = PDBParser()
        p.get_structure("test", out_path)
        warnings.resetwarnings()
