import logging
from collections import defaultdict

import MDAnalysis
import gravis as gv
import networkx as nx
from src.utils.file_loader import SimulationFile

from ..protein.mdanalysis_to_protseq import (
    _mk_residue_atom_adjacency_matrix,
    AdjacencyMatrix,
    _fragment_extract_residues,
)
from ...constants import (
    NETWORKX_COLORMAP_ATOMS,
)
from ...utils.repair_simulation import fix_missing_elements

log = logging.getLogger("resolver")


def _mk_adjacency_matrices(
    simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
) -> tuple[list[str], AdjacencyMatrix, list[str], AdjacencyMatrix]:
    """
    Generates adjacency matrices for atoms and residues
    :param simulation: of which is the fragment part of (must be loaded)
    :param fragment: representing a molecule
    :return: node names (atom elements/residue names) and adjacency matrices for both atom and residue resolution
    """

    atom_indices: dict[int, int] = {
        atom_idx: atom_i for atom_i, atom_idx in enumerate(fragment.atoms.ids)
    }

    residues: dict[int, str]
    residue_indices: dict[int, int]
    residues, residue_indices = _fragment_extract_residues(
        fragment, validate_residue_names=False
    )

    residue_adjacency_matrix: AdjacencyMatrix
    atom_adjacency_matrix: AdjacencyMatrix
    residue_adjacency_matrix, atom_adjacency_matrix = _mk_residue_atom_adjacency_matrix(
        fragment, residues, atom_indices, residue_indices
    )

    if not simulation.supports_atomic_elements(fragment.segments.segids[0]):
        fix_missing_elements(
            simulation, fragment, replace_with_mass=True
        )  # partial results are okay
    atom_names: list[str] = fragment.atoms.elements
    residue_names: list[str] = [resname for resname in residues.values()]

    return atom_names, atom_adjacency_matrix, residue_names, residue_adjacency_matrix


def _populate_graph_nodes(
    graph: nx.Graph, names: list[str], colormap: defaultdict | None = None
) -> None:
    """
    Populates node in the graph (node IDs are equivalent to node indices)
    :param graph: to be populated
    :param names: of the nodes
    :param colormap: if not None, translation of name.upper() to color
    :return: None
    """

    if colormap is None:
        colormap = defaultdict(lambda: None)

    for i, node_name in enumerate(names):
        graph.add_node(i)
        graph.nodes[i]["color"] = colormap[node_name.upper()]
        graph.nodes[i]["label"] = node_name if node_name else "?"


def _populate_graph_edges(graph: nx.Graph, adjacency_matrix: AdjacencyMatrix) -> None:
    """
    Populates edges in the graph (node IDs are equivalent to node indices)
    :param graph: to be populated
    :param adjacency_matrix: representing bonds between edges
    :return: None
    """

    for node_from in range(len(adjacency_matrix)):
        for node_to in range(len(adjacency_matrix)):
            if (
                node_from == node_to or node_from > node_to
            ):  # skip self-loops and symmetry
                continue
            if adjacency_matrix[node_from][node_to]:
                graph.add_edge(node_from, node_to, weight=1, color="black")


def _export_graph(graph: nx.Graph, to_path: str) -> None:
    """
    Exports the graph as interactive HTML page
    :param graph: to be rendered
    :param to_path: to save rendered graph
    :return: None
    """
    fig = gv.d3(
        graph,
        use_node_size_normalization=True,
        node_size_normalization_max=30,
        use_edge_size_normalization=True,
        edge_size_data_source="weight",
        edge_curvature=0.3,
        zoom_factor=0.6,
        show_node_label=True,
        show_menu=True,
        show_menu_toggle_button=True,
        node_label_data_source="label",
    )
    fig.export_html(to_path)


def fragment_to_bond_maps(
    simulation: SimulationFile,
    fragment: MDAnalysis.AtomGroup,
    save_path_atom: str,
    save_path_residue: str,
) -> bool:
    """
    Transforms a fragment into networkx interactive graphs, representing atoms/residues and bonds between them
    :param simulation: of which is the fragment part of (must be loaded)
    :param fragment: representing a molecule
    :param save_path_atom: to save the atom-bonds graph
    :param save_path_residue: to save the residue-bonds graph
    :return: True if successful, False otherwise
    """
    try:
        atom_names, atom_adjacency, residue_names, residue_adjacency = (
            _mk_adjacency_matrices(simulation, fragment)
        )
        graph_atoms, graph_residues = nx.Graph(), nx.Graph()

        _populate_graph_nodes(graph_atoms, atom_names, colormap=NETWORKX_COLORMAP_ATOMS)
        _populate_graph_nodes(graph_residues, residue_names)

        _populate_graph_edges(graph_atoms, atom_adjacency)
        _populate_graph_edges(graph_residues, residue_adjacency)

        _export_graph(graph_atoms, save_path_atom)
        _export_graph(graph_residues, save_path_residue)
    except Exception as exc:
        log.error(f"Failed to generate bond maps: {exc}")
        return False
    return True
