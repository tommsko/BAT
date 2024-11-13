import argparse
import configparser
import logging
import os.path
import random
import textwrap

from rich import print

from .constants import APP_VERSION
from .utils.file_loader import SimulationFile, load_simulations_in_directory


def assert_arguments(
    arguments: argparse.Namespace, requirements: dict[str, str]
) -> None:
    """
    Checks that all necessary arguments are present
    :param arguments: result of argparse.parse_args()
    :param requirements: mapping of fields to their descriptions
    :return: None (raises runtime error on failure)
    """
    for argument_name, argument_description in requirements.items():
        if not hasattr(arguments, argument_name):
            print(f"[b][red]✘ missing argument '{argument_name}'[/red][/b]")
            print(f"[b]{argument_name}[/b] - {argument_description}")
            raise RuntimeError(f"Missing argument '{argument_name}'")
    print("[b][green]✔ all necessary arguments provided[/green][/b]")


def random_quote(from_file: str = "src/quotes.txt") -> str | None:
    """
    If quote file exists, returns one random quote from it
    :param from_file: path to the file where one line == one quote
    :return: one random quote, or None
    """
    if not os.path.exists(from_file):
        return None

    all_quotes: list[str] = []
    with open(from_file, "r") as fd:
        for line in fd:
            if not line.startswith("#"):
                all_quotes.append(line)
    return random.choice(all_quotes)


def print_mode_header(
    mode_name: str,
    mode_description: str,
    argument_names: list[str],
    arguments: argparse.Namespace,
) -> None:
    """
    Prints basic header and list of values provided
    :param mode_name: name of the mode
    :param mode_description: its description
    :param argument_names: arguments used during the mode
    :param arguments: result of argparse.parse_args()
    :return: None
    """
    app_name: str = f"          Biomolecule Annotation Toolkit v. {APP_VERSION}          "
    print("[b]" + "  BAT  ".center(len(app_name), "=") + "[/b]")
    print(app_name)
    print()
    if (quote := random_quote()) is not None:
        for line in textwrap.wrap(
            f'"{quote.strip()}"', len(app_name), break_long_words=False
        ):
            print("[i]" + line + "[/i]")
        print()

    print("[b]" + f"   {mode_name} mode   ".center(len(app_name), "=") + "[/b]")
    for line in textwrap.wrap(mode_description, len(app_name), break_long_words=False):
        print(line)
    print()

    print("[b]Arguments[/b]")
    for argument in argument_names:
        if hasattr(arguments, argument) and getattr(arguments, argument) is not None:
            print(f"\t[i]{argument}[/i]:\t{getattr(arguments, argument)}")
        else:
            print(f"\t[i]{argument}[/i]:\t<not specified>")


def init_loggers(verbose: bool | None, very_verbose: bool | None) -> None:
    """
    Initializes loggers used by the application
    :param verbose: argument parameter (from argparse) for verbose mode
    :param very_verbose: argument parameter (from argparse) for very verbose mode
    :return: None
    """
    verbose = verbose if verbose is not None else False
    very_verbose = very_verbose if very_verbose is not None else False

    def _add_log_handler(
        log_name: str, very_verbose_level: int, verbose_level: int, normal_level: int
    ):
        handler: logging.StreamHandler = logging.StreamHandler()
        formatter = logging.Formatter(
            "[%(process)s] [%(name)s] %(asctime)s %(levelname)s %(message)s"
        )
        handler.setFormatter(formatter)

        log: logging.Logger = logging.getLogger(log_name)
        log.handlers.clear()
        log.addHandler(handler)
        log.setLevel(
            very_verbose_level
            if very_verbose
            else verbose_level if verbose else normal_level
        )

    _add_log_handler(
        "base",
        very_verbose_level=logging.DEBUG,
        verbose_level=logging.INFO,
        normal_level=logging.INFO,
    )
    _add_log_handler(
        "resolver",
        very_verbose_level=logging.DEBUG,
        verbose_level=logging.INFO,
        normal_level=logging.ERROR,
    )
    print(
        f"[b][green]✔ loggers initialized (verbose: {verbose}, very verbose: {very_verbose})[/green][/b]"
    )


def load_simulations_filtered(
    source_directory: str,
    simulations_filter: list[str] | None,
    molecule_filter: list[str] | None,
    config: configparser.ConfigParser,
) -> list[SimulationFile]:
    """
    Quick wrapper to load simulations including filters on them from the arguments given
    :param source_directory: containing molecular dynamics trajectories (.tpr files)
    :param simulations_filter: if not None, names of simulations to be loaded (included), otherwise all
    :param molecule_filter: if not None, names of molecules to be loaded (included), otherwise all
    :param config: main configuration
    :return: list of simulations
    """
    print("[b]listing simulations in source directory...[/b]", end="\r")

    simulations: list[SimulationFile] = list(
        reversed(load_simulations_in_directory(source_directory, config))
    )
    print(f"[b]found {len(simulations)} simulations...[/b]", end="\r")

    if simulations_filter:
        print("[b]filtering based on simulation name...[/b]", end="\r")
        simulations = [
            s for s in simulations if s.simulation_name in simulations_filter
        ]
    if molecule_filter:
        print("[b]filtering based on molecule name...[/b]", end="\r")
        filter_molecules: set[str] = set(molecule_filter)
        for sim in simulations:
            sim.filter_molecules = filter_molecules
    print(" " * 100, end="\r")
    print(f"[b][green]✔ found {len(simulations)} simulations[/green][/b]")
    return simulations
