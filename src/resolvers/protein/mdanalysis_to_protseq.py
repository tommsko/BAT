import configparser
import logging
from collections import deque
from typing import Callable

import MDAnalysis

from ..base.fragment_resolver import SignatureGenerationError
from ...constants import KNOWN_AMINO_ACIDS

log = logging.getLogger("resolver")

AdjacencyMatrix = list[list[bool]]


def _fragment_extract_residues(
    fragment: MDAnalysis.AtomGroup, validate_residue_names: bool = True
) -> tuple[dict[int, str], dict[int, int]]:
    """
    Extracts protein residues from a fragment and validates them
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :param dvalidate_residue_names: validate residue names as protein residues
    :return: mapping from residue indices (MDAnalysis) to 1-letter residues' code (constants.KNOWN_AMINO_ACIDS)
            and mapping of residue idx to residue indices
    :throws: RuntimeError when unknown residue is found
    """

    residues: dict[int, str] = {}
    residue_indices: dict[int, int] = {}
    res_i: int = 0
    for res_name, res_idx in zip(fragment.resnames, fragment.resindices):
        res_name = res_name.upper()
        if validate_residue_names and res_name not in KNOWN_AMINO_ACIDS:
            raise SignatureGenerationError(
                "extracting residues",
                f"unknown residue encountered at residue idx '{res_idx}': {res_name}",
            )
        if res_idx not in residues:
            residues[res_idx] = (
                KNOWN_AMINO_ACIDS[res_name] if validate_residue_names else res_name
            )
            residue_indices[res_idx] = res_i
            res_i += 1
            continue

        if (
            residues[res_idx] != KNOWN_AMINO_ACIDS[res_name]
            if validate_residue_names
            else res_name
        ):
            if validate_residue_names:
                raise RuntimeError(
                    "extracting residues",
                    f"residue 3-letter code inconsistent across the same residue id. "
                    f"Stored: '{KNOWN_AMINO_ACIDS[res_name]}' and encountered '{residues[res_idx]}'",
                )
            else:
                residues[res_idx] = res_name
    return residues, residue_indices


def _mk_adjacency_matrix(nodes_cnt: int) -> AdjacencyMatrix:
    """
    Creates adjacency matrix for given nodes
    :param nodes_cnt: number of nodes
    :return: adjacency matrix for all nodes
    """
    return [[False for _ in range(nodes_cnt)] for _ in range(nodes_cnt)]


def _populate_adjacency_matrix_inplace(
    fragment: MDAnalysis.AtomGroup,
    atom_indices: dict[int, int],
    residue_indices: dict[int, int],
    adjacency_matrix: AdjacencyMatrix,
    is_residue_adjacency: bool,
) -> int:
    """
    Populates adjacency matrix for given bonds in fragment
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :param atom_indices: mapping from atom ids to atom order
    :param residue_indices: mapping from residue ids to residue order
    :param adjacency_matrix: result of _mk_adjacency_matrix()
    :param is_residue_adjacency: if True, populates it with residue adjacency, otherwise atom adjacency
    :return: number of nodes (atoms or residues) visited, see is_residue_adjacency
    """

    nodes_visited: list[bool] = [
        False
        for _ in range(
            len(residue_indices) if is_residue_adjacency else len(atom_indices)
        )
    ]
    for atom_idx1, atom_idx2 in fragment.bonds.indices:
        atom_i1, atom_i2 = atom_indices[atom_idx1], atom_indices[atom_idx2]
        atom_resid_i1, atom_resid_i2 = (
            residue_indices[fragment.atoms.resindices[atom_i1]],
            residue_indices[fragment.atoms.resindices[atom_i2]],
        )

        node_i1, node_i2 = (
            (atom_resid_i1, atom_resid_i2)
            if is_residue_adjacency
            else (atom_i1, atom_i2)
        )
        nodes_visited[node_i1] = nodes_visited[node_i2] = True
        adjacency_matrix[node_i1][node_i2] = adjacency_matrix[node_i2][node_i1] = True
    return nodes_visited.count(True)


def _adjacency_matrix_get_neighbors(
    node_from: int, adjacency_matrix: AdjacencyMatrix
) -> set[int]:
    """
    Finds all linked nodes (neighbors) to a nodes
    :param node_from: node's indice
    :param adjacency_matrix: adjacency matrix between the nodes
    :return: set of all neighbors, except queried node
    """
    result: set[int] = set()
    for neighbor in range(len(adjacency_matrix)):
        if node_from == neighbor:
            continue
        if (
            adjacency_matrix[node_from][neighbor]
            or adjacency_matrix[neighbor][node_from]
        ):
            result.add(neighbor)
    return result


def _bfs(
    node_start: int,
    adjacency_matrix: AdjacencyMatrix,
    reuse_visited: list[bool] | None = None,
) -> tuple[list[bool], list[int], list[int], int]:
    """
    Generic BFS implementation
    :param node_start: node indice from which to start BFS
    :param adjacency_matrix: adjacency matrix of nodes
    :param reuse_visited: if provided, re-uses visited matrix
    :return: list of visited nodes (possibly reused),
             list of distances to start_node,
             list of predecessors given BFS search,
             number of nodes visited during BFS search
    """

    visited: list[bool] = (
        [False for _ in range(len(adjacency_matrix))]
        if reuse_visited is None
        else reuse_visited
    )
    distances: list[int] = [-1 for _ in range(len(adjacency_matrix))]
    prev: list[int] = [-1 for _ in range(len(adjacency_matrix))]
    nodes_visited: int = 0

    queue: deque[int] = deque()
    queue.append(node_start)
    visited[node_start] = True
    distances[node_start] = 0
    prev[node_start] = node_start

    while queue:
        node_from: int = queue.popleft()
        nodes_visited += 1
        for node_to in _adjacency_matrix_get_neighbors(node_from, adjacency_matrix):
            if visited[node_to]:
                continue
            visited[node_to] = True
            distances[node_to] = distances[node_from] + 1
            prev[node_to] = node_from
            queue.append(node_to)

    return visited, distances, prev, nodes_visited


def _discontinuous_segment_cnt(adjacency_matrix: AdjacencyMatrix) -> int:
    """
    Finds how many weakly-linked components there in adjacency matrix
    :param adjacency_matrix: to be analysed
    :return: number of weakly-linked components
    """
    visited = [False for _ in range(len(adjacency_matrix))]
    components: int = 0
    for i_from in range(len(adjacency_matrix)):
        if visited[i_from]:
            continue
        _bfs(i_from, adjacency_matrix, visited)
        components += 1
    return components


def _degrees_count(adjacency_matrix: AdjacencyMatrix) -> list[int]:
    """
    Calculates in and out degree count for each node in adjacency matrix
    :param adjacency_matrix: to be analysed
    :return: degrees for each node
    """
    degrees: list[int] = [0 for _ in range(len(adjacency_matrix))]
    for node_from in range(len(adjacency_matrix)):
        for node_to in range(len(adjacency_matrix)):
            if node_from != node_to and adjacency_matrix[node_from][node_to]:
                degrees[node_from] += 1  # symmetric
    return degrees


def _bfs_extract_longest_paths(
    node_start: int,
    distances: list[int],
    prev: list[int],
    visited: list[bool],
    max_length: int = 10000,
) -> list[list[int]]:
    """
    Extracts the path(s) between node_start and most distant node(s) utilizing results of BFS calculation before
    :param node_start: the residue *INDICE, NOT IDX* from which BFS was run
    :param distances: distances list from BFS run
    :param prev: previous list from BFS run
    :param visited: visited list from BFS run
    :param max_length: maximum length of the path (to throw on cycles)
    :return: path from node_start to node with the biggest distance, or more if there are multiple nodes with such distance
    """
    biggest_distance: int = max(distances)
    nodes_to: list[int] = [
        node_i
        for node_i, node_distance in enumerate(distances)
        if node_distance == biggest_distance
    ]
    results: list[list[int]] = []

    for node_to in nodes_to:
        path: list[int] = []
        while node_to != node_start:
            assert visited[node_to], "corrupted BFS calculation"
            path.append(node_to)
            node_to = prev[node_to]
            if len(path) > max_length:
                raise SignatureGenerationError(
                    "finding protein sequence",
                    f"protein residues' path is too long or cyclic (max_length = {max_length})",
                )
        path.append(node_start)
        results.append(list(reversed(path)))
    return results


def _split_residue_indices_methionine(
    residues: dict[int, str]
) -> tuple[set[int], set[int]]:
    """
    Splits residues' indices to two groups, methionine and non-methionine ones
    :param residues: result from _fragment_extract_residues(...)
    :return: set of methionine indices, set of non-methionine indices
    """
    is_methionine: set[int] = set()
    not_methionine: set[int] = set()

    for residue_id, residue_code in residues.items():
        if residue_code == "M":
            is_methionine.add(residue_id)
        else:
            not_methionine.add(residue_id)
    return is_methionine, not_methionine


def _translate_residue_indices_to_sequence(
    indices_path: list[int], residue_indices: dict[int, int], residues: dict[int, str]
) -> str:
    """
    Translates residue indices path from BFS to sequence of 1-code amino acids
    :param indices_path: path of residue indices
    :param residue_indices: mapping of residue IDs to residue indices
    :param residues: mapping of residues IDs to residue codes
    :return: sequence string
    """
    res_i_to_idx: dict[int, int] = {
        resid_i: resid_idx for resid_idx, resid_i in residue_indices.items()
    }
    return "".join(residues[res_i_to_idx[x]] for x in indices_path)


def _find_longest_sequences_bfs(
    from_indices: set[int],
    residue_adjacency_matrix: AdjacencyMatrix,
    residue_indices: dict[int, int],
    residues: dict[int, str],
) -> list[str]:
    """
    Finds longest sequences given residues and their adjacency matrix and filters them
    :param from_indices: set of residue IDs from which longest sequence is being searched
    :param residue_adjacency_matrix: adjacency matrix of residues
    :param residue_indices: mapping from residue ids to residue order
    :param residues: mapping of residue ids to residue codes
    :return: list of sequences found
    """
    sequences: list[str] = []
    for residue_idx in from_indices:
        node_from: int = residue_indices[residue_idx]
        vis, dist, prev, _ = _bfs(node_from, residue_adjacency_matrix)
        paths = _bfs_extract_longest_paths(node_from, dist, prev, vis)
        for path in paths:
            sequences.append(
                _translate_residue_indices_to_sequence(path, residue_indices, residues)
            )
    return sequences


def _filter_sequences(
    sequences: list[str], require: Callable[[str], bool]
) -> list[str]:
    """
    Simple wrapper to filter a list of strings
    :param sequences: to be filtered
    :param require: lambda that must hold true to include sequence in the result
    :return: list of sequences where require(...) returns True
    """
    return [seq for seq in sequences if require(seq)]


def _fragment_to_protein_sequence_verify_sequentiality(
    residue_adjacency_matrix: AdjacencyMatrix, atom_adjacency_matrix: AdjacencyMatrix
) -> None:
    """
    Verifies that residues and atoms form one component, and that residues form one sequential chain
    :param residue_adjacency_matrix: adjacency matrix of residues
    :param atom_adjacency_matrix: adjacency matrix of atoms
    :return: None
    :raises SignatureGenerationError: of there are multiple components or degrees mismatch
    """

    error_step: str = "extracting residue chain"

    if (
        residue_components := _discontinuous_segment_cnt(residue_adjacency_matrix)
    ) != 1:
        raise SignatureGenerationError(
            error_step,
            f"found {residue_components} disjoint residue components, "
            f"which is more than 1 supported!",
        )
    if (atoms_components := _discontinuous_segment_cnt(atom_adjacency_matrix)) != 1:
        raise SignatureGenerationError(
            error_step,
            f"found {atoms_components} disjoint atom components, "
            f"which is more than 1 supported! (although there's only one residue component)",
        )

    residue_degrees: list[int] = _degrees_count(residue_adjacency_matrix)
    if any(d >= 3 or d == 0 for d in residue_degrees):
        raise SignatureGenerationError(
            error_step, "found a residue with degree > 2, sequence is non-linear!"
        )
    if residue_degrees.count(1) >= 3:
        raise SignatureGenerationError(
            error_step,
            "found more than three residues with degree == 1, sequence is non-linear!",
        )


def _mk_residue_atom_adjacency_matrix(
    fragment: MDAnalysis.AtomGroup,
    residues: dict[int, str],
    atom_indices: dict[int, int],
    residue_indices: dict[int, int],
) -> tuple[AdjacencyMatrix, AdjacencyMatrix]:
    """
    Creates residue and atom adjacency matrix given
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :param residues: mapping of residue ids to residue codes
    :param atom_indices: mapping of atom ids to atom indices
    :param residue_indices: mapping of residue ids to residue indices
    :return: residue and atom adjacency matrix
    """

    residue_adjacency_matrix: AdjacencyMatrix = _mk_adjacency_matrix(len(residues))
    _populate_adjacency_matrix_inplace(
        fragment,
        atom_indices,
        residue_indices,
        residue_adjacency_matrix,
        is_residue_adjacency=True,
    )

    atom_adjacency_matrix: AdjacencyMatrix = _mk_adjacency_matrix(len(atom_indices))
    _populate_adjacency_matrix_inplace(
        fragment,
        atom_indices,
        residue_indices,
        atom_adjacency_matrix,
        is_residue_adjacency=False,
    )
    return residue_adjacency_matrix, atom_adjacency_matrix


def cant_determine_err(strict: bool, partial: bool) -> None:
    """
    Just raises an exception
    :param strict: protein sequence determination is in strict mode
    :param partial: protein sequence determination is in partial mode
    """
    raise SignatureGenerationError(
        "protein fingerprint resolution",
        f"Unable to determine the sequence! (modes: exact=True, strict={strict}, "
        f"partial={partial})",
    )


def _fragment_to_protein_mk_adjacency_matrices(fragment: MDAnalysis.AtomGroup) -> ...:
    """
    Creates indices and adjacency matrices for residues and atoms given
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :return: atom_indices - mapping of atom idx to atom indices
             atom_adjacency_matrix - adjacency matrix of atoms, indices-based
             residues - mapping of residue idx to residue code
             residue_indices - mapping of residue idx to residue indices
             residue_adjacency_matrix - adjacency matrix of residues, indices-based
    """

    atom_indices: dict[int, int] = {
        atom_idx: atom_i for atom_i, atom_idx in enumerate(fragment.atoms.ids)
    }

    residues: dict[int, str]
    residue_indices: dict[int, int]
    residues, residue_indices = _fragment_extract_residues(fragment)

    residue_adjacency_matrix: AdjacencyMatrix
    atom_adjacency_matrix: AdjacencyMatrix
    residue_adjacency_matrix, atom_adjacency_matrix = _mk_residue_atom_adjacency_matrix(
        fragment, residues, atom_indices, residue_indices
    )
    return (
        atom_indices,
        atom_adjacency_matrix,
        residues,
        residue_indices,
        residue_adjacency_matrix,
    )


def fragment_to_protein_sequence(fragment: MDAnalysis.AtomGroup, config: configparser.ConfigParser) -> str:
    """
    Tries to translate fragment into sequence of amino-acids
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :param config: main configuration
    :return: string sequence of amino-acids given the fragment
    :raises: SignatureGenerationError if fragment cannot be translated into protein sequence
    """
    PROT_SIGN_PARTIAL = config.getboolean('protein', 'determine_sequence_partial_mode')
    PROT_SIGN_STRICT_MODE = config.getboolean('protein', 'determine_sequence_strict_mode')

    (
        atom_indices,
        atom_adjacency_matrix,
        residues,
        residue_indices,
        residue_adjacency_matrix,
    ) = _fragment_to_protein_mk_adjacency_matrices(fragment)

    _fragment_to_protein_sequence_verify_sequentiality(
        residue_adjacency_matrix, atom_adjacency_matrix
    )

    methionine_indices, non_methionine_indices = _split_residue_indices_methionine(
        residues
    )

    methionine_sequences: list[str] = _find_longest_sequences_bfs(
        methionine_indices, residue_adjacency_matrix, residue_indices, residues
    )
    if methionine_sequences:
        if sel := _filter_sequences(
            methionine_sequences, require=lambda seq: len(seq) == len(residues)
        ):
            sel.sort(key=lambda seq: len(seq), reverse=True)
            return sel[0]

    log.debug(
        "... unable to determine protein sequence given exact & strict mode. "
        "If allowed, continuing exact & non-strict mode"
    )
    if PROT_SIGN_STRICT_MODE:
        cant_determine_err(PROT_SIGN_STRICT_MODE, PROT_SIGN_PARTIAL)

    non_methionine_sequences: list[str] = _find_longest_sequences_bfs(
        non_methionine_indices, residue_adjacency_matrix, residue_indices, residues
    )
    if sel := _filter_sequences(
        non_methionine_sequences, require=lambda seq: len(seq) == len(residues)
    ):
        sel.sort(key=lambda seq: len(seq), reverse=True)
        return sel[0]

    log.debug(
        "... unable to determine protein sequence given exact & non-strict mode. "
        "If allowed, continuing partial & non-strict mode"
    )
    if not PROT_SIGN_PARTIAL:
        cant_determine_err(PROT_SIGN_STRICT_MODE, PROT_SIGN_PARTIAL)

    all_sequences = methionine_sequences + non_methionine_sequences

    all_sequences.sort(key=lambda seq: len(seq), reverse=True)
    if all_sequences:
        return all_sequences[0]

    cant_determine_err(PROT_SIGN_STRICT_MODE, PROT_SIGN_PARTIAL)


def _dfs(
    results: list[list[int]],
    node_start: int,
    adjacency_matrix: AdjacencyMatrix,
    required_len: int | None = None,
    reuse_visited: list[bool] | None = None,
    current_path: list[int] | None = None,
    current_depth: int = 1,
):
    """
    Generic DFS implementation (reuses nodes)
    :param results: list where to store traversed paths
    :param node_start: node indice from which to start BFS
    :param adjacency_matrix: adjacency matrix of nodes
    :param required_len: required length of path to be exported into results
    :param reuse_visited: if provided, re-uses visited matrix
    :param current_path: of the DFS search, none at start
    :return: None (appends to <<results>> out parameter)
    """
    if current_path is None:
        current_path = [node_start]

    if reuse_visited is None:
        visited: list[bool] = [False for _ in range(len(adjacency_matrix))]
        visited[node_start] = True
    else:
        visited: list[bool] = reuse_visited

    is_leaf: bool = True
    for node_to in _adjacency_matrix_get_neighbors(node_start, adjacency_matrix):
        if visited[node_to]:
            continue
        is_leaf = False  # get to go deeper

        current_path.append(node_to)
        visited[node_to] = True  # gray, but will be unset
        _dfs(
            results,
            node_to,
            adjacency_matrix,
            required_len,
            visited,
            current_path,
            current_depth + 1,
        )
        visited[node_to] = False
        current_path.pop()
    visited[node_start] = True  # black

    if is_leaf and (required_len is None or (required_len == current_depth)):
        results.append(current_path.copy())


def fragment_to_protein_sequences_unsafe(fragment: MDAnalysis.AtomGroup) -> list[str]:
    """
    It is possible that all residues are bonded together (mixed regular and hydrogen bonds in gromacs simulations).
    Extracts all possible (longest) sequences given a malformed fragment
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :return: string sequences of amino-acids given the fragment
    :raises: SignatureGenerationError if fragment cannot be translated into protein sequence
    """
    (
        atom_indices,
        atom_adjacency_matrix,
        residues,
        residue_indices,
        residue_adjacency_matrix,
    ) = _fragment_to_protein_mk_adjacency_matrices(fragment)

    log.critical(
        f"Attempting *UNSAFE* protein sequence resolution for {len(residues)} residues"
    )

    all_sequences_found: set[str] | list[str] = set()
    for residue_from in residue_indices:
        for min_len in range(len(residues), -1, -1):
            log.critical(f"...testing from {residue_from}: min_depth: {min_len}")
            paths: list[list[int]] = []
            _dfs(paths, residue_from, residue_adjacency_matrix, required_len=min_len)
            paths_str: list[str] = [
                _translate_residue_indices_to_sequence(path, residue_indices, residues)
                for path in paths
            ]
            for sequence in paths_str:
                all_sequences_found.add(sequence)
            if paths:
                break

    all_sequences_found = list(all_sequences_found)
    all_sequences_found.sort(key=lambda x: len(x), reverse=True)
    longest_sequence_len: int = len(all_sequences_found[0])
    return [
        sequence
        for sequence in all_sequences_found
        if len(sequence) == longest_sequence_len
    ]


def fragment_to_protein_sequence_residue_based(fragment: MDAnalysis.AtomGroup):
    """
    Makes strong assumption residues were added in the order of protein and extracts the sequence as ordered
    chain of residues from the simulation
    :param fragment: AtomGroup representing a molecule (protein specifically)
    :return: string sequences of amino-acids given the fragment
    """
    try:
        residue_atom_indices: list[int] = [
            indices[0] for indices in fragment.residues.indices
        ]

        atom_idx_to_indices: dict[int, int] = {
            idx: i for i, idx in enumerate(fragment.atoms.ids)
        }

        residue_atom_resnames: list[str] = [
            KNOWN_AMINO_ACIDS[fragment.resnames[atom_idx_to_indices[i]].upper()]
            for i in residue_atom_indices
        ]
        residues_and_indices: list[tuple[int, str]] = list(
            zip(residue_atom_indices, residue_atom_resnames)
        )
        residues_and_indices.sort(key=lambda x: x[0])  # sort by indices
        return "".join(x[1] for x in residues_and_indices)
    except Exception as exc:
        raise SignatureGenerationError(
            "protein fingerprint resolution (residue-based)",
            f"fingerprint generation failed: {exc}",
        )
