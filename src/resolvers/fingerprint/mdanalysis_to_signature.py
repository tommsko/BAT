import hashlib

import MDAnalysis


def hash_fragment(fragment: MDAnalysis.AtomGroup) -> str:
    """
    Produces good-enough deterministic hash of molecular fragment given its atoms and bonds
    :param fragment: to be hashed
    :return: SHA512 hash of fragment
    """
    atom_count: int = fragment.n_atoms
    bonds: list[tuple[int, int]] = fragment.bonds.indices
    bonds.sort()
    return str(hashlib.sha512(f"{atom_count}_{bonds}".encode("utf-8")).hexdigest())
