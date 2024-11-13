import configparser
import json
import logging
import os.path
from typing import Any, Callable

import MDAnalysis

from ..resolvers.base.fragment_resolver import ResolverBase, IdentifierResolutionError, SignatureGenerationError
from ..resolvers import get_resolvers
from ..resolvers.base.resolver_database import ResolverCache
from rich import print
from src.utils.file_loader import SimulationFile

log = logging.getLogger("base")


def run_catchall(callable: Callable, args: list[Any], default: Any):
    """
    Executes callable with arguments, on any exeption returns default value
    :param callable: to be executed
    :param args: to be passes (unpacked)
    :param default: to be returned on a failure
    :return: result of the callable on success, default on any exception
    """
    try:
        return callable(*args)
    except (NotImplementedError, SignatureGenerationError, IdentifierResolutionError):
        return default  # this is okay
    except Exception as exc:
        log.exception(exc)
        return default


def color_str_rg(name: str, is_green: bool) -> str:
    """
    :param name: to be printed
    :param is_green: if True, colors the name green, red otherwise
    :return: formatted string
    """
    if is_green:
        return f"<[green]{name}[/green]>"
    return f"<[red]{name}[/red]>"


class SimulationExporter:
    def __init__(self, simulations: list[SimulationFile], out_directory: str, configuration: configparser.ConfigParser):
        log.info(f"SimulationExporter started with {len(simulations)} simulations...")

        cache: ResolverCache = ResolverCache(configuration)
        self.resolvers: list[ResolverBase] = get_resolvers(cache, configuration)
        log.info(f"SimulationExporter loaded {len(self.resolvers)} resolvers...")

        self.simulations: list[SimulationFile] = simulations
        self.out_directory: str = out_directory
        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)

    @staticmethod
    def _check_out_dir(path: str) -> bool:
        """
        Checks whether output directory exists and it's empty
        :param path: to the output directory
        :return: True if output directory exists (or was created) and it's empty
        """
        if not os.path.exists(path):
            os.makedirs(path)
            log.debug("Output directory does not exist, creating new")
            return True

        if len(os.listdir(path)) > 0:
            log.error("Output directory exists and contains data. Cowardly exiting!")
            return False

        return True

    def try_export(self) -> None:
        """
        Attempts to export all the simulations given the resolvers
        :return: None
        """
        if not SimulationExporter._check_out_dir(self.out_directory):
            log.error("Export cancelled!")
            return

        print("[cyan][b]=== Export of the simulations ===[/b][/cyan]")
        for simulation in self.simulations:
            self._try_export_simulation(simulation)

    def _try_export_simulation(self, simulation: SimulationFile) -> None:
        """
        Attempts to export all molecules in the simulation
        :param simulation: to be exported
        :return: None
        """
        if not simulation.load_simulation():
            print("[red][b]Simulation failed to load[/b][/red]")
            return

        os.mkdir(os.path.join(self.out_directory, simulation.simulation_name))

        print(f"[b]Simulation '{simulation.simulation_name}'[/b]")

        for segment_name in simulation.list_segments():
            self._try_export_molecule(simulation, segment_name)

        simulation.unload_simulation()

    def _try_export_molecule(self, simulation: SimulationFile, molecule_name: str) -> None:
        """
        Attempts to export the molecule from simulation
        :param simulation: to be exported
        :param molecule_name: to be exported
        :return: None
        """
        out_dir: str = os.path.join(os.path.join(self.out_directory, simulation.simulation_name), molecule_name)
        os.mkdir(out_dir)

        print(f"[b]... molecule '{molecule_name}'[/b]")
        molecule: MDAnalysis.AtomGroup = simulation.get_molecule(molecule_name)
        export_name: str = f"{simulation.simulation_name}_{molecule_name}"

        for resolver in self.resolvers:
            if not resolver.allow_debug_output:
                continue

            is_applicable: bool = resolver.is_applicable(simulation, molecule)

            signature: str | None = None if not is_applicable \
                else run_catchall(resolver.get_signature, [simulation, molecule], default=None)
            if signature is not None:
                with open(os.path.join(out_dir,
                                       f"{export_name}_{resolver.resolver_name}_signature.txt"), 'w') as fd:
                    fd.write(signature)

            debug_export: bool = False if not is_applicable \
                else run_catchall(resolver.generate_debug_data,
                                  [simulation, molecule, out_dir, export_name], default=False)

            idents: str | None = None if not is_applicable or signature is None \
                else run_catchall(resolver.get_identifiers,
                                  [signature], default=None)

            if idents is not None:
                with open(os.path.join(out_dir,
                                       f"{export_name}_{resolver.resolver_name}_idents.json"), 'w') as fd:
                    json.dump(idents, fd)

            idents_sim: str | None = None if not is_applicable or signature is None \
                else run_catchall(resolver.get_similar_identifiers,
                                  [signature], default=None)

            if idents_sim is not None:
                with open(os.path.join(out_dir,
                                       f"{export_name}_{resolver.resolver_name}_possible_idents.json"), 'w') as fd:
                    json.dump(idents_sim, fd)

            print(f"...... resolver '{resolver.resolver_name}': "
                  f"{color_str_rg('applicable', is_applicable)} "
                  f"{color_str_rg('fingerprint', signature is not None)} "
                  f"{color_str_rg('debug export', debug_export)} "
                  f"{color_str_rg('identity idents', idents is not None)} "
                  f"{color_str_rg('similarity idents', idents_sim is not None)} ")
