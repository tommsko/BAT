import configparser
import json
import logging
from dataclasses import dataclass
from tabulate import tabulate
from rich import print

from src.utils.file_loader import SimulationFile
from ..utils.results_file import ResultsFile
from ..resolvers.base.resolver_database import ResolverCache
from ..resolvers import get_resolvers
from ..resolvers.base.fragment_resolver import ResolverKind, ResolverBase

log = logging.getLogger("base")

TableHeader = str
TableHeaderDescription = str
TableHeaderPack = list[tuple[TableHeader, TableHeaderDescription]]

@dataclass
class SimulationStats:
    input_files: list[str]
    loaded_files: list[str]
    fully_resolved_files: list[str]
    partially_resolved_files: list[str]
    unresolved_files: list[str]


simulation_stats_table_headers: TableHeaderPack = [
    ('input', 'number of simulation files inputted for resolution'),
    ('loaded', 'number of simulation files which could be loaded (.tpr version was high enough)'),
    (' (%) ', 'percentage of the above compared to input simulations'),
    ('fully resolved', 'number of simulations in which all molecules were resolved on identity or similarity'),
    (' (%) ', 'percentage of the above compared to loaded simulations'),
    ('partially resolved', 'number of simulations in which only some molecules were resolved on identity or similarity'),
    (' (%) ', 'percentage of the above compared to loaded simulations'),
    ('unresolved', 'number of simulations in which no molecules were resolved on identity or similarity'),
    (' (%) ', 'percentage of the above compared to loaded simulations'),
]


@dataclass
class MoleculeResolutionStats:
    loaded_molecules: list[tuple[str, str]]
    identity_matched_molecules: list[tuple[str, str]]
    similarity_matched_molecules: list[tuple[str, str]]
    unmatched_molecules: list[tuple[str, str]]
    unmatched_molecules_primary_signature: list[tuple[str, str]]
    unmatched_molecules_fallback_signature: list[tuple[str, str]]
    unmatched_molecules_no_signature: list[tuple[str, str]]


molecule_stats_table_headers: TableHeaderPack = [
    ('loaded', 'number of molecules loaded (from "loaded" simulations)'),
    ('identity match', 'number of molecules matched on identity level'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
    ('similarity match', 'number of molecules matched on similarity level'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
    ('no match (NM)', 'number of molecules with no matches'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
    ('NM: primary sign.', 'number of molecules with no matches, but with primary signature generated'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
    ('NM: fallback sign.', 'number of molecules with no matches, but with fallback signature generated'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
    ('NM: no sign.', 'number of molecules with no matches, with no signature generated'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
]


@dataclass
class ResolverStatistics:
    signature_generated_cnt: int
    identity_found_cnt: int
    similarity_found_cnt: int


resolver_stats_table_headers: TableHeaderPack = [
    ('resolver', 'number of molecules loaded (from "loaded" simulations)'),
    ('signature generated', 'number of molecules for which resolver generated signature (!!!)'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
    ('identity match', 'number of molecules for which resolver found identity match (!!!)'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
    ('similarity match', 'number of molecules for which resolver found similarity match (!!!)'),
    (' (%) ', 'percentage of the above compared to loaded molecules'),
]



def tableheaderpack_to_stdout(pack: TableHeaderPack) -> None:
    """
    Prints TableHeaderPack in a nice manner
    :param pack: TableHeaderPack to print
    :return: None
    """
    print("[b]Table description[/b]")
    for header, description in pack:
        print(f"\t[b]{header}[/b] - [i]{description}[/i]")


def tableheaderpack_to_header(pack: TableHeaderPack) -> list[str]:
    """
    Extract table headers from pack
    :param pack: TableHeaderPack
    :return: list of table headers
    """
    return [header for header, _ in pack]


class ResolutionStatistics:
    def __init__(self, simulations: list[SimulationFile], results_dir: str, configuration: configparser.ConfigParser) -> None:
        self._primary_resolvers: set[str] = set()
        self._fallback_resolvers: set[str] = set()
        self._all_resolvers: list[str] = []
        self._init_resolvers_kind(configuration)

        self._simulation_stats: SimulationStats = SimulationStats([], [], [], [], [])
        self._molecule_stats: MoleculeResolutionStats = MoleculeResolutionStats([], [], [], [], [], [], [])
        self._resolver_stats: dict[str, ResolverStatistics] = {
            resolver_name: ResolverStatistics(0, 0, 0)
            for resolver_name in self._all_resolvers
        }

        results_files: list[ResultsFile] = [ResultsFile(sim, results_dir) for sim in simulations]
        log.info(f"ResolutionStatistics loaded {len(results_files)} files...")
        for file in results_files:
            self._process_file(file)
        log.info("ResolutionStatistics aggregated all files...")

    def _init_resolvers_kind(self, configuration: configparser.ConfigParser) -> None:
        """
        Populates categories of primary/fallback reasolvers for molecular resolution statistics
        :param configuration: for initializination of resolvers
        :return: None
        """
        cache: ResolverCache = ResolverCache(configuration)
        resolvers: list[ResolverBase] = get_resolvers(cache, configuration)
        for resolver in resolvers:
            self._all_resolvers.append(resolver.resolver_name)  # order is important here
            if resolver.resolver_kind == ResolverKind.PRIMARY:
                self._primary_resolvers.add(resolver.resolver_name)
            elif resolver.resolver_kind == ResolverKind.FALLBACK:
                self._fallback_resolvers.add(resolver.resolver_name)
            else:
                assert False, "unknown resolver kind"

    def _process_file(self, file: ResultsFile) -> None:
        """
        Loads the results file aggregates statistics of molecules' resolution
        :param file: ResultsFile of the simulation
        :param sim_stats: where to store simulation (global) resolution statistics
        :param molecule_stats: where to store molecular resolution statistics
        :param resolver_stats: where to store molecular resolution statistics per resolver
        :return: None
        """
        if not self._process_file_simulation_stats(file):  # loads the results file
            return  # not even loadable
        self._process_file_molecule_stats(file)
        self._process_file_resolver_stats(file)

    def _process_file_simulation_stats(self, file: ResultsFile) -> bool:
        """
        Loads the results file and populates basic simulation resolution statistics
        :param file: ResultsFile of the simulation
        :return: True if file was loaded successfully
        """
        self._simulation_stats.input_files.append(file.path)

        if file.exists():  # if file does not exist, simulation's loading failed (anything is saved afterward)
            self._simulation_stats.loaded_files.append(file.path)
        else:
            return False

        if not file.load(failsafe=False, segments_from_results=True, ignore_compatibility=True):
            return False

        if file.is_fully_resolved():
            self._simulation_stats.fully_resolved_files.append(file.path)
        elif file.is_partially_resolved():
            self._simulation_stats.partially_resolved_files.append(file.path)
        else:
            self._simulation_stats.unresolved_files.append(file.path)
        return True

    def _process_file_molecule_stats(self, file: ResultsFile) -> None:
        """
        Reads the results file and populates basic molecular resolution statistics
        :param file: ResultsFile of the simulation
        :return: None
        """
        for segment in file.list_segments():
            molecule: tuple[str, str] = file.path, segment
            self._molecule_stats.loaded_molecules.append(molecule)

            if file.is_segment_resolved(segment, only_identity_resolved=True):
                self._molecule_stats.identity_matched_molecules.append(molecule)
            elif file.is_segment_resolved(segment):
                self._molecule_stats.similarity_matched_molecules.append(molecule)
            else:
                self._molecule_stats.unmatched_molecules.append(molecule)
                resolvers_generated_signature: list[str] = file.get_resolvers_generated_signature(segment)
                if any(resolver in self._primary_resolvers for resolver in resolvers_generated_signature):
                    self._molecule_stats.unmatched_molecules_primary_signature.append(molecule)
                elif any(resolver in self._fallback_resolvers for resolver in resolvers_generated_signature):
                    self._molecule_stats.unmatched_molecules_fallback_signature.append(molecule)
                else:
                    self._molecule_stats.unmatched_molecules_no_signature.append(molecule)

    def _process_file_resolver_stats(self, file: ResultsFile) -> None:
        """
        Reads the results file and populates basic resolver statistics
        :param file: ResultsFile of the simulation
        :return: None
        """
        for segment in file.list_segments():
            for resolver in file.get_resolvers_generated_signature(segment):
                self._resolver_stats[resolver].signature_generated_cnt += 1

            if not file.is_segment_resolved(segment):
                continue
            if file.is_segment_resolved(segment, only_identity_resolved=True):
                self._resolver_stats[file.get_segment_resolver(segment)].identity_found_cnt += 1
            elif file.is_segment_resolved(segment):
                self._resolver_stats[file.get_segment_resolver(segment)].similarity_found_cnt += 1

    def _simulation_stats_stdout(self) -> None:
        """
        Prints simulation statistics to stdout
        :return: None
        """
        print("[b] === General (whole simulation) resolution statistics === [/b]")
        print()
        tableheaderpack_to_stdout(simulation_stats_table_headers)
        print()

        loaded_files: int = len(self._simulation_stats.loaded_files)
        table_headers: list[str] = tableheaderpack_to_header(simulation_stats_table_headers)
        table_content: list[list[str]] = [[
            str(len(self._simulation_stats.input_files)),
            str(loaded_files),
            str(round(loaded_files / len(self._simulation_stats.input_files) * 100, 3)) + " %",
            str(len(self._simulation_stats.fully_resolved_files)),
            str(round(len(self._simulation_stats.fully_resolved_files) / loaded_files * 100, 3)) + " %",
            str(len(self._simulation_stats.partially_resolved_files)),
            str(round(len(self._simulation_stats.partially_resolved_files) / loaded_files * 100, 3)) + " %",
            str(len(self._simulation_stats.unresolved_files)),
            str(round(len(self._simulation_stats.unresolved_files) / loaded_files * 100, 3)) + " %",
        ]]
        print(tabulate(table_content, table_headers, tablefmt="simple_grid"))

    def _simulation_stats_json(self) -> dict:
        """
        Generates dictionary containing simulation statistics
        :return: dict
        """
        return {
            'input_files': self._simulation_stats.input_files,
            'loaded_files': self._simulation_stats.loaded_files,
            'fully_resolved_files': self._simulation_stats.fully_resolved_files,
            'partially_resolved_files': self._simulation_stats.partially_resolved_files,
            'unresolved_files': self._simulation_stats.unresolved_files
        }

    def _molecular_stats_stdout(self) -> None:
        """
        Prints molecular resolution statistics to stdout
        :return: None
        """
        print("[b] === Molecular resolution statistics === [/b]")
        print()
        tableheaderpack_to_stdout(molecule_stats_table_headers)
        print()

        loaded_molecules: int = len(self._molecule_stats.loaded_molecules)
        table_headers: list[str] = tableheaderpack_to_header(molecule_stats_table_headers)
        table_content: list[list[str]] = [[
            str(loaded_molecules),
            str(len(self._molecule_stats.identity_matched_molecules)),
            str(round(len(self._molecule_stats.identity_matched_molecules) / loaded_molecules * 100, 3)) + " %",
            str(len(self._molecule_stats.similarity_matched_molecules)),
            str(round(len(self._molecule_stats.similarity_matched_molecules) / loaded_molecules * 100, 3)) + " %",
            str(len(self._molecule_stats.unmatched_molecules)),
            str(round(len(self._molecule_stats.unmatched_molecules) / loaded_molecules * 100, 3)) + " %",
            str(len(self._molecule_stats.unmatched_molecules_primary_signature)),
            str(round(len(self._molecule_stats.unmatched_molecules_primary_signature) / loaded_molecules * 100, 3)) + " %",
            str(len(self._molecule_stats.unmatched_molecules_fallback_signature)),
            str(round(len(self._molecule_stats.unmatched_molecules_fallback_signature) / loaded_molecules * 100, 3)) + " %",
            str(len(self._molecule_stats.unmatched_molecules_no_signature)),
            str(round(len(self._molecule_stats.unmatched_molecules_no_signature) / loaded_molecules * 100, 3)) + " %",
        ]]
        print(tabulate(table_content, table_headers, tablefmt="simple_grid"))

    def _molecular_stats_json(self) -> dict:
        """
        Generates dictionary containing molecular resolution statistics
        :return: dict
        """
        return {
            'loaded_molecules': self._molecule_stats.loaded_molecules,
            'identity_matched_molecules': self._molecule_stats.identity_matched_molecules,
            'similarity_matched_molecules': self._molecule_stats.similarity_matched_molecules,
            'unmatched_molecules': self._molecule_stats.unmatched_molecules,
            'unmatched_molecules_signatures': {
                'with_primary_signature': self._molecule_stats.unmatched_molecules_primary_signature,
                'with_fallback_signature': self._molecule_stats.unmatched_molecules_fallback_signature,
                'with_no_signature': self._molecule_stats.unmatched_molecules_no_signature
            }
        }

    def _resolver_stats_stdout(self) -> None:
        """
        Prints resolver performance statistics to stdout
        :return: None
        """
        print("[b] === Resolver performance statistics === [/b]")
        print()
        tableheaderpack_to_stdout(resolver_stats_table_headers)
        print()
        print("[red][b]!!! beware of resolution process !!![/b][/red]")
        print("[red]Resolvers are attempted in order (of the table below, top ~ highest priority), "
              "meaning if higher-order resolver succeeds, others are not even attempted. Thus signature generated count "
              "has meaning across the whole table, but other values are only sensibly interpreted in the row only [/red]")
        print()

        table_headers: list[str] = tableheaderpack_to_header(resolver_stats_table_headers)

        loaded_molecules: int = len(self._molecule_stats.loaded_molecules)
        table_content: list[list[str]] = [
            [resolver,
             str(self._resolver_stats[resolver].signature_generated_cnt),
             str(round(self._resolver_stats[resolver].signature_generated_cnt / loaded_molecules * 100, 3)) + " %",
             str(self._resolver_stats[resolver].identity_found_cnt),
             str(round(self._resolver_stats[resolver].identity_found_cnt / loaded_molecules * 100, 3)) + " %",
             str(self._resolver_stats[resolver].similarity_found_cnt),
             str(round(self._resolver_stats[resolver].similarity_found_cnt / loaded_molecules * 100, 3)) + " %"]
            for resolver in self._all_resolvers
        ]
        print(tabulate(table_content, table_headers, tablefmt="simple_grid"))

    def _resolver_stats_json(self) -> dict:
        """
        Generates dictionary containing revolver performance statistics
        :return: dict
        """
        return {
            resolver: {
                'signatures_generated': self._resolver_stats[resolver].signature_generated_cnt,
                'identity_match': self._resolver_stats[resolver].identity_found_cnt,
                'similarity_match': self._resolver_stats[resolver].similarity_found_cnt
            } for resolver in self._all_resolvers
        }

    def print_stats_stdout(self) -> None:
        """
        Prints all statistics to stdout
        :return: None
        """
        print("[cyan][b] ===== Resolution statistics ===== [/b][/cyan]")
        print()
        self._simulation_stats_stdout()
        print()
        self._molecular_stats_stdout()
        print()
        self._resolver_stats_stdout()

    def save_stats_json(self, path: str) -> None:
        """
        Saves all statistics to JSON file
        :param path: to the json file
        :return: None
        """
        with open(path, 'w') as fd:
            json.dump({
                'resolution_statistics': {
                    'simulation_statistics': self._simulation_stats_json(),
                    'molecular_statistics': self._molecular_stats_json(),
                    'resolver_statistics': self._resolver_stats_json()
                }
            }, fd)
