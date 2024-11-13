import json
import logging
import os.path
from json import JSONDecodeError

from .file_loader import SimulationFile
from ..constants import ERR_SEGMENT_NOT_FOUND

log = logging.getLogger("base")


class ResultsFile:
    """
    Wrapper around self-saving JSON file
    """

    def __init__(self, simulation: SimulationFile, output_directory: str) -> None:
        """
        Initializes results file given
        :param simulation: associated simulation for which this result is for
                           (doesn't have to be loaded until calling .load() for ResultsFile is called)
        :param output_directory: where the file should be saved
        """
        self.path: str = os.path.join(
            output_directory, f"{simulation.simulation_name}.json"
        )
        self.data: dict | None = None

        self.simulation: SimulationFile = simulation
        self.segments: list[str] = (
            []
        )  # segment names (molecule types) in the simulation

    def exists(self) -> bool:
        """
        :return: True if the underlying result file exists
        """
        return os.path.exists(self.path)

    def _read(self) -> bool:
        """
        Attempts to read the results file (might fail)
        :return: True if successful, False otherwise
        """
        try:
            with open(self.path, "r") as fd:
                self.data = json.load(fd)
                return True
        except JSONDecodeError:
            return False

    def load(
        self,
        failsafe: bool = True,
        segments_from_results: bool = False,
        ignore_compatibility: bool = False,
    ) -> bool:
        """
        Loads the results file, failsafe (creates an empty on critical failure)
        :param failsafe: on load failure, creates an empty results file
        :param segments_from_results: if True, segment names are fetched from simulation itself (verification)
                                      if False, segment names are fetched from already existing results file
        :param ignore_compatibility: if True, skips compatibility verification checks
        :return: True if successful, False otherwise
        """
        if failsafe:
            assert (
                not segments_from_results
            ), "failsafe requires segment names to be extracted from simulation"

        if not self.exists():
            log.debug(f"Results file '{self.path}' is missing")
            if failsafe: self._mk_empty(self.simulation)
            else: return False

        if not self._read():
            log.warning(f"Unable to load results file '{self.path}'")
            if failsafe: self._mk_empty(self.simulation)
            else: return False

        if segments_from_results:
            if 'segment_resolution' not in self.data:
                return False
            self.segments = list(self.data["segment_resolution"].keys())
        else:
            self.segments = self.simulation.list_segments()

        if not ignore_compatibility and not self._check_compatibility(self.simulation):
            log.warning(
                f"Results file '{self.path}' is corrupted or contains illegal values"
            )
            if failsafe: self._mk_empty(self.simulation)
            else: return False
        return True

    def save(self) -> None:
        """
        Saves the results file
        :return: None
        """
        assert self.data is not None, "unable to save, nothing loaded so far"
        with open(self.path, "w") as fd:
            json.dump(self.data, fd, indent=4)

    def _check_compatibility_fingerprint(self, simulation: SimulationFile) -> bool:
        """
        Checks whether the simulation resolution results file belongs to the same simulation
        utilizing basic fingerprints of the simulation
        Tries to do some corrections if applicable
        :param simulation: of which results file belongs to
        :return: True if so, False otherwise
        """

        def _check_field(field_name: str, expected_value: int | str) -> bool:
            """
            Checks whether results contain specific field with specific value
            Tries to do some corrections if applicable
            :param field_name: name of the top-level field
            :param expected_value: of such field
            :return: True if the value is the same as expected (or was empty and corrected), False otherwise
            """
            if field_name in self.data:
                if self.data[field_name] != expected_value:
                    return False
            else:
                self.data[field_name] = expected_value
            return True

        if not _check_field("name", simulation.simulation_name):
            return False
        if not _check_field("atoms", simulation.atom_count()):
            return False
        if not _check_field("residues", simulation.residue_count()):
            return False
        if not _check_field("segments", simulation.segment_count()):
            return False
        if not _check_field("fragments", simulation.fragment_count()):
            return False

        if "errors" not in self.data:
            self.data["errors"] = []

        return True

    def _check_compatibility_resolvers(self, simulation: SimulationFile) -> bool:
        """
        Checks whether the simulation resolution results file contains necessary segment/resolver information
        Tries to do some corrections if applicable
        :param simulation: of which results file belongs to
        :return: True if so, False otherwise
        """

        # if segments were improperly exported during resolution, the resolution is considered invalid
        if "segment_resolution" not in self.data:
            return False
        if set(simulation.list_segments()) != set(
            self.data["segment_resolution"].keys()
        ):
            return False

        for segment in simulation.list_segments():
            if (
                "resolver_generated_signature"
                not in self.data["segment_resolution"][segment]
            ):
                self.data["segment_resolution"][segment][
                    "resolver_generated_signature"
                ] = []
            if "resolver" not in self.data["segment_resolution"][segment]:
                self.data["segment_resolution"][segment]["resolver"] = None
            if "ident_type" not in self.data["segment_resolution"][segment]:
                self.data["segment_resolution"][segment]["ident_type"] = None
            if "resolved_idents" not in self.data["segment_resolution"][segment]:
                self.data["segment_resolution"][segment]["resolved_idents"] = []
            if "possible_idents" not in self.data["segment_resolution"][segment]:
                self.data["segment_resolution"][segment]["possible_idents"] = []
            if "errors" not in self.data["segment_resolution"][segment]:
                self.data["segment_resolution"][segment]["errors"] = []

        return True

    def _check_compatibility(self, simulation: SimulationFile) -> bool:
        """
        Checks whether the simulation results file is compatible with current version of the application
        Tries to do some corrections if applicable
        :param simulation: of which results file belongs to
        :return: True if so, False otherwise
        """
        if not self._check_compatibility_fingerprint(simulation):
            return False
        if not self._check_compatibility_resolvers(simulation):
            return False
        self.save()
        return True

    def _mk_empty(self, simulation: SimulationFile, no_segments: bool = False) -> bool:
        """
        Creates an empty results file
        :param simulation: of which results file belongs to
        :param no_segments: if set to True, does not populate segments (and bad things will happen)
        :return: True
        """
        log.debug(
            f"... created empty results file for simulation '{simulation.simulation_name}'"
        )
        self.data = {
            "name": simulation.simulation_name,
            "atoms": simulation.atom_count(),
            "residues": simulation.residue_count(),
            "segments": simulation.segment_count(),
            "fragments": simulation.fragment_count(),
            "errors": [],
            "segment_resolution": {
                segment_name: {
                    "errors": [],
                    "resolver_generated_signature": [],
                    "resolver": None,
                    "ident_type": None,
                    "resolved_idents": [],
                    "possible_idents": [],
                }
                for segment_name in (simulation.list_segments() if not no_segments else [])
            },
        }
        self.save()
        return True

    def is_segment_resolved(
        self,
        segment_name: str,
        only_identity_resolved: bool = False,
        only_similarity_resolved: bool = False,
    ) -> bool:
        """
        True if the segment is fully or partially resolved
        :param segment_name: to check for
        :param only_identity_resolved: accept only identity resolution
        :param only_similarity_resolved: accept similarity resolution
        :return: True if so, False otherwise
        """

        assert not (
                only_identity_resolved and only_similarity_resolved
        ), "identity and similarity resolution cannot happen at the same time"

        assert segment_name in self.segments, ERR_SEGMENT_NOT_FOUND.format(segment_name)
        segment_resolution: dict = self.data["segment_resolution"][segment_name]
        if 'resolver' not in segment_resolution:
            return False

        if only_identity_resolved:
            return (
                segment_resolution["resolver"] is not None
                and segment_resolution["resolved_idents"]
            )
        elif only_similarity_resolved:
            return (
                segment_resolution["resolver"] is not None
                and segment_resolution["possible_idents"]
            )
        else:
            return segment_resolution["resolver"] is not None and (
                segment_resolution["resolved_idents"]
                or segment_resolution["possible_idents"]
            )

    def get_segment_resolver(self, segment_name: str) -> str | None:
        """
        Fetches resolver name for a given segment
        :param segment_name: to fetch resolver of
        :return: name of the resolver on full or partial resolution, None otherwise
        """
        assert segment_name in self.segments, ERR_SEGMENT_NOT_FOUND.format(segment_name)
        return self.data["segment_resolution"][segment_name]["resolver"]

    def get_resolvers_generated_signature(self, segment_name: str) -> list[str]:
        """
        Fetches resolver names for all resolvers which generated signature for the segment
        :param segment_name: of the segment
        :return: list of resolvers
        """
        assert segment_name in self.segments, ERR_SEGMENT_NOT_FOUND.format(segment_name)
        return self.data["segment_resolution"][segment_name][
            "resolver_generated_signature"
        ] if 'resolver_generated_signature' in self.data["segment_resolution"][segment_name] else []

    def get_idents(self, segment_name: str) -> list[str]:
        """
        Fetches resolved identity identifiers
        :param segment_name: to fetch identifiers from
        :return: list of identifiers
        """
        assert segment_name in self.segments, ERR_SEGMENT_NOT_FOUND.format(segment_name)
        return self.data["segment_resolution"][segment_name]["resolved_idents"]

    def get_similarity_idents(self, segment_name: str) -> list[tuple[str, float]]:
        """
        Fetches resolved similarity identifiers
        :param segment_name: to fetch identifiers from
        :return: list of identifiers with included similarity
        """
        assert segment_name in self.segments, ERR_SEGMENT_NOT_FOUND.format(segment_name)
        return self.data["segment_resolution"][segment_name]["possible_idents"]

    def get_idents_type(self, segment_name: str) -> str | None:
        """
        Fetches type of identifiers
        :param segment_name: to fetch identifiers from
        :return: list of identifiers
        """
        assert segment_name in self.segments, ERR_SEGMENT_NOT_FOUND.format(segment_name)
        return self.data["segment_resolution"][segment_name]["ident_type"]

    def is_fully_resolved(
        self, only_identity_resolved: bool = False, only_similarity_resolved: bool = False
    ) -> bool:
        """
        Checks whether all segments in the simulation are fully or partially resolved
        :param only_identity_resolved: accept only identity resolution
        :param only_similarity_resolved: accept only similarity resolution
        :return: True if all segments are, False otherwise
        """
        return all(
            self.is_segment_resolved(
                segment_name, only_identity_resolved, only_similarity_resolved
            )
            for segment_name in self.segments
        )

    def is_partially_resolved(
        self, only_identity_resolved: bool = False, only_similarity_resolved: bool = False
    ) -> bool:
        """
        Checks whether some segments in the simulation are fully or partially resolved
        :param only_identity_resolved: accept only identity resolution
        :param only_similarity_resolved: accept similarity resolution
        :return: True if some segments are, False otherwise
        """
        return any(
            self.is_segment_resolved(
                segment_name, only_identity_resolved, only_similarity_resolved
            )
            for segment_name in self.segments
        )

    def list_segments(self) -> list[str]:
        """
        :return: list of all the segments in results file
        """
        return self.segments
