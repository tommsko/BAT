import configparser
import logging
import os.path
import warnings
from collections import defaultdict
from os import listdir
from os.path import isfile, join

import MDAnalysis
from MDAnalysis import Universe

from src.constants import ERROR_SIM_NOT_LOADED_MSG

log = logging.getLogger("base")


class SimulationFile:
    """
    SimulationFile represents one GROMACS simulation, which might consist
    of multiple files (e.g. trajectory, atomic information, parameters etc...).
    Workflow:
      ctor ~> initialize simulation file(s) by tpr simulation parameters
      load_simulation() ~> load simulation file into RAM (expensive)
      ... do whatever you want to do
      unload_simulation() ~> unload simulation file
    """

    OTHER_SUPPORTED_FORMATS: set[str] | None = None  # to be loaded from configuration

    def __init__(
        self,
        tpr_path: str,
        config: configparser.ConfigParser,
        files_in_directory: list[str] | None = None,
    ) -> None:
        """
        Initializes the simulation file with GROMACS Input File (.tpr) and other related files from lookup
        .tpr path has to be provided, other extension are scanned in the same directory, given same basename of file
        :param tpr_path: path to the simulation input file
        :param config: configuration
        :param files_in_directory: other files in directory related to lookup (to not list files multiple times)
        """
        self.config: configparser.ConfigParser = config
        SimulationFile._init_supported_formats(self.config)

        self.paths: list[str] = [
            os.path.abspath(tpr_path)
        ]  # of files related to the simulation
        self.simulation_name: str = ""  # basename of GROMACS input file w/o extension
        self._lookup_related_files(tpr_path, files_in_directory)

        self.universe: Universe | None = (
            None  # managed with load_simulation() and unload_simulation()
        )

        self.segments: list[str] = (
            []
        )  # names of segments (molecule types) in the simulation
        self._segment_supports_atomic_elements: dict[str, bool] = defaultdict(
            lambda: False
        )
        self._segment_supports_atomic_masses: dict[str, bool] = defaultdict(
            lambda: False
        )
        self._segment_supports_atomic_positions: dict[str, bool] = defaultdict(
            lambda: False
        )
        self._segment_supports_atomic_names: dict[str, bool] = defaultdict(
            lambda: False
        )

        self.filter_molecules: set[str] = (
            set()
        )  # if not empty, only these segments will be reported

    @staticmethod
    def _init_supported_formats(config: configparser.ConfigParser) -> None:
        """
        Lists all supported input formats from configuration (except .tpr input files) and saves them as static variable
        On subsequent calls does nothing
        :param config: main configuration
        :return: None
        :raises RuntimeError: if such list does not contain tpr
        """
        if SimulationFile.OTHER_SUPPORTED_FORMATS is not None:
            return

        supported_formats: set[str] = set()
        for fmt_suffix, enable in config.items("supported_input_formats"):
            if enable:
                supported_formats.add("." + fmt_suffix.upper())
        if ".TPR" not in supported_formats:
            raise RuntimeError(
                "configuration: .tpr is not listed in supported formats, not proceeding"
            )
        supported_formats.remove(".TPR")
        SimulationFile.OTHER_SUPPORTED_FORMATS = supported_formats

    def _lookup_related_files(
        self, tpr_path: str, files_in_directory: list[str] | None
    ) -> None:
        """
        Finds all associated files of the simulation (defined by .tpr) in the directory of tpr file.
        Updates simulation_name.
        Please see <<supported_input_formats>> section of configuration
        :param tpr_path: path to the simulation parameters file
        :param files_in_directory: other files in directory related to lookup (to skip listing files)
        :return: None (modifies stored paths)
        """
        working_directory: str = os.path.dirname(tpr_path)

        if files_in_directory is None:
            files_in_directory: list[str] = [
                os.path.abspath(join(working_directory, path))
                for path in listdir(working_directory)
                if isfile(join(working_directory, path))
            ]

        self.simulation_name = os.path.basename(tpr_path).split(".")[0]

        do_lookup: bool = self.config.getboolean("input_lookup", "perform_lookup")
        relaxed_include: bool = self.config.getboolean(
            "input_lookup", "input_lookup_relaxed_name_requirements"
        )

        if not do_lookup:
            return

        def _include_condition(_path: str, _extension: str) -> bool:
            """
            :param path: of the possibly associated file
            :param extension: of the possibly associated file
            :return: True if the file is associated by name, False otherwise
            """
            if (
                os.path.basename(_path).upper()
                == f"{self.simulation_name}{_extension}".upper()
            ):
                return True
            if (
                relaxed_include
                and os.path.basename(_path)
                .upper()
                .startswith(self.simulation_name.upper())
                and os.path.basename(_path).upper().endswith(_extension.upper())
            ):
                return True
            return False

        for extension in SimulationFile.OTHER_SUPPORTED_FORMATS:
            files: list[str] = [
                path
                for path in files_in_directory
                if _include_condition(path, extension)
            ]
            self.paths.extend(files)

    def _load_simulation_metadata(self) -> None:
        """
        Determines simulation name and values for the flags of the simulation and its segments (molecules)
        :return: None
        """
        if self.filter_molecules:
            self.segments = list(
                set(self.universe.segments.segids).intersection(
                    set(self.filter_molecules)
                )
            )
        else:
            self.segments = list(set(self.universe.segments.segids))

        for segment in self.list_segments():
            self._load_flag_supports_atomic_elements(segment)
            self._load_flag_supports_atomic_masses(segment)
            self._load_flag_supports_atomic_positions(segment)
            self._load_flag_supports_atomic_names(segment)

    def load_simulation(self) -> bool:
        """
        Attempts to load the simulation into the memory
        :return: True if successful, False otherwise
        """
        try:
            with warnings.catch_warnings(action="ignore"):
                self.universe = Universe(*self.paths, tpr_resid_from_one=True, in_memory=True)
                self._load_simulation_metadata()
                return True
        except Exception as exc:
            log.warning(
                f"unable to load simulation {self.simulation_name} (at {self.paths}): {exc}"
            )
            log.warning(exc, exc_info=True)
            return False

    def get_simulation(self) -> Universe:
        """
        Returns simulation loaded by this SimulationFile
        :return: None
        :raises RuntimeError: if load_simulation() wasn't called before and the loading failed
        """
        if self.universe is None and not self.load_simulation():
            raise RuntimeError(ERROR_SIM_NOT_LOADED_MSG)
        return self.universe  # not None here

    def unload_simulation(self) -> None:
        """
        Unloads the simulation from the memory (makes it available for garbage collector)
        :return: None
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        self.universe = None

    def _load_flag_supports_atomic_positions(self, segment_name: str) -> None:
        """
        Checks whether the given segment (molecule) contains positions of the atoms
        and updates internal flags accordingly
        :param segment_name: to check for
        :return: None
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        self._segment_supports_atomic_positions[segment_name] = hasattr(
            self.get_molecule(segment_name).atoms, "positions"
        )

    def supports_atomic_positions(self, segment_name: str) -> bool:
        """
        Checks whether the segment contains atomic positions information
        :param segment_name: to check for
        :return: True if simulation does provide positions of atoms
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self._segment_supports_atomic_positions[segment_name]

    def _load_flag_supports_atomic_elements(self, segment_name: str) -> None:
        """
        Checks whether the given segment (molecule) contains elements of the atoms
        and updates internal flags accordingly
        :param segment_name: to check for
        :return: None
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        segment: MDAnalysis.AtomGroup = self.get_molecule(segment_name)
        self._segment_supports_atomic_elements[segment_name] = hasattr(
            segment, "elements"
        ) and all(e for e in segment.elements)

    def supports_atomic_elements(self, segment_name: str) -> bool:
        """
        Checks whether the segment contains atomic elements information
        :param segment_name: to check for
        :return: True if simulation does provide elements of atoms
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self._segment_supports_atomic_elements[segment_name]

    def _load_flag_supports_atomic_masses(self, segment_name: str) -> None:
        """
        Checks whether the given segment (molecule) contains masses of the atoms
        and updates internal flags accordingly
        :param segment_name: segment to check for
        :return: None
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        segment: MDAnalysis.AtomGroup = self.get_molecule(segment_name)
        self._segment_supports_atomic_masses[segment_name] = hasattr(
            segment, "masses"
        ) and all(x > 0 for x in segment.masses)

    def supports_atomic_masses(self, segment_name: str) -> bool:
        """
        Checks whether the segment provides atomic masses information
        :param segment_name: to check for
        :return: True if simulation does provide masses of atoms
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self._segment_supports_atomic_masses[segment_name]

    def _load_flag_supports_atomic_names(self, segment_name: str) -> None:
        """
        Checks whether the given segment (molecule) contains names of the atoms
        and updates internal flags accordingly
        :param segment_name: to check for
        :return: None
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        segment: MDAnalysis.AtomGroup = self.get_molecule(segment_name)
        self._segment_supports_atomic_names[segment_name] = hasattr(
            segment.atoms, "names"
        )

    def supports_atomic_names(self, segment_name: str) -> bool:
        """
        Checks whether the segment provides atomic names information
        :param segment_name: to check for
        :return: True if simulation does provide masses of atoms
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self._segment_supports_atomic_names[segment_name]

    def list_segments(self) -> list[str]:
        """
        :return: list of segments (molecule types) in the simulation
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self.segments

    def atom_count(self) -> int:
        """
        :return: number of atoms in the simulation
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self.universe.atoms.n_atoms

    def residue_count(self) -> int:
        """
        :return: number of residues in the simulation
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self.universe.atoms.n_residues

    def segment_count(self) -> int:
        """
        :return: number of segments in the simulation
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self.universe.atoms.n_segments

    def fragment_count(self) -> int:
        """
        :return: number of fragments in the simulation
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self.universe.atoms.n_fragments

    def get_molecule(self, segment_name: str) -> MDAnalysis.AtomGroup:
        """
        Returns just one molecule given the segment (molecule type)
        :param segment_name: of the molecule
        :return: AtomGroup representing the single molecule
        """
        assert self.universe is not None, ERROR_SIM_NOT_LOADED_MSG
        return self.universe.select_atoms(f"segid {segment_name}").fragments[0]

    def __str__(self) -> str:
        """
        :return: basic information about SimulationFile in text format
        """
        return f"""[Simulation] {self.simulation_name}\n""" + "\n".join(
            [f"\t-> {path}" for path in self.paths]
        )


def load_simulations_in_directory(
    directory: str, config: configparser.ConfigParser
) -> list[SimulationFile]:
    """
    Finds all simulation input files (.tpr) incl. associated files in directory and loads them into SimulationFile
    :param directory: path to the directory with simulations
    :param config: main configuration
    :return: list of SimulationFiles
    """

    all_files: list[str] = [
        os.path.abspath(join(directory, path))
        for path in listdir(directory)
        if isfile(join(directory, path))
    ]

    tpr_paths: list[str] = [path for path in all_files if path.lower().endswith(".tpr")]
    return [
        SimulationFile(path, config, files_in_directory=all_files) for path in tpr_paths
    ]
