import json
import logging
import multiprocessing
import os.path
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Any

from src.utils.file_loader import SimulationFile
from src.utils.results_file import ResultsFile


class LogError(Enum):
    ERROR_NO_LOAD = 'ERROR_NO_LOAD'
    ERROR_TIMEOUT = 'ERROR_TIMEOUT'
    ERROR_NOT_FINISHED = 'ERROR_NOT_FINISHED'
    ERROR_RESOLUTION_NOT_STARTED = 'ERROR_RESOLUTION_NOT_STARTED'
    ERROR_NO_MATCH = 'ERROR_NO_MATCH'


class StatusUpdateType(Enum):
    SIMULATION_RESOLUTION_STARTED = 0  # value ~> results file location
    SIMULATION_SUCCESSFULLY_LOADED = 1  # value ~> UNUSED
    SIMULATION_RESOLUTION_FINISHED = 2  # value ~> UNUSED
    SIMULATION_RESOLUTION_TIMEOUT = 3  # value ~> UNUSED
    MOLECULE_RESOLUTION_STARTED = 4  # value ~> molecule name
    MOLECULE_RESOLVER_APPLICABLE = 5  # value ~> (molecule name, resolver name)
    MOLECULE_SIGNATURE_GENERATED = 6  # value ~> (molecule name, resolver name)
    MOLECULE_IDENTITY_MATCH = 7  # value ~> molecule name
    MOLECULE_SIMILARITY_MATCH = 8  # value ~> molecule name
    MOLECULE_NO_MATCH = 9  # value ~> molecule name


SimulationSourceFile = str
StatusUpdateValue = str
StatusUpdate = tuple[SimulationSourceFile, StatusUpdateType, StatusUpdateValue]


@dataclass
class SimulationResolutionStatus:
    tpr_path: str
    results_path: str
    loaded: bool
    finished_successfully: bool
    finished_timeout: bool
    molecules: defaultdict[str, dict[StatusUpdateType, list | bool]]


def emit_update(queue: multiprocessing.Queue,
                source_tpr: str,
                update_type: StatusUpdateType,
                value: Any) -> None:
    """
    Simple wrapper for emitting to status queue
    :param queue: to be emitted to
    :param source_tpr: source file being resolved
    :param update_type: kind of the update
    :param value: associated with the status update
    :return: None
    """
    queue.put((source_tpr, update_type, value))


def emit_update_skipped_molecules(
        simulation: SimulationFile,
        results: ResultsFile,
        status: multiprocessing.Queue,
        emit_unmatched: bool = True,
):
    """
    Emits update messages for already resolved (skipped) molecules
    :param simulation: to be resolved
    :param results: ResultsFile linked to the simulation
    :param status: queue where to put updates about the resolution
    :param emit_unmatched: emit messages for so-far unmatched molecules too (skipping whole simulation)
    :return: None
    """
    simulation_tpr: str = simulation.paths[0]

    for segment in simulation.list_segments():
        if results.is_segment_resolved(segment, only_identity_resolved=True):
            resolver: str = results.get_segment_resolver(segment)
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_RESOLUTION_STARTED, segment)
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_SIGNATURE_GENERATED, (segment, resolver))
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_IDENTITY_MATCH, (segment, resolver))
            continue

        if results.is_segment_resolved(segment, only_similarity_resolved=True):
            resolver: str = results.get_segment_resolver(segment)
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_RESOLUTION_STARTED, segment)
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_SIGNATURE_GENERATED, (segment, resolver))
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_SIMILARITY_MATCH, (segment, resolver))
            continue

        if emit_unmatched:
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_RESOLUTION_STARTED, segment)
            emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_NO_MATCH, segment)
            continue


log = logging.getLogger("base")


class StatusQueue:
    """
    StatusQueue provides aggregated log for a multithreaded resolution process
    """

    def __init__(self, workers_cnt: int) -> None:
        """
        Initializes the status queue
        :param workers_cnt: number of threads expected to be run
        """
        self._queue: multiprocessing.Queue = multiprocessing.Queue(maxsize=workers_cnt * 10)
        self._current_status: dict[str, Any] = {}
        self._status: dict[str, SimulationResolutionStatus] = {}

    def get_queue(self) -> multiprocessing.Queue:
        """
        :return: the queue object bound by this instance
        """
        return self._queue

    def process_queue(self) -> None:
        """
        Reads all incoming status updates from the queue
        :return: None
        """
        while not self.get_queue().empty():
            simulation_file, update_type, update_data = self.get_queue().get()
            match update_type:
                case StatusUpdateType.SIMULATION_RESOLUTION_STARTED:
                    self._status[simulation_file] = SimulationResolutionStatus(
                        tpr_path=simulation_file,
                        results_path=update_data,
                        loaded=False,
                        finished_successfully=False,
                        finished_timeout=False,
                        molecules=defaultdict(lambda: {
                            StatusUpdateType.MOLECULE_RESOLUTION_STARTED: False,
                            StatusUpdateType.MOLECULE_RESOLVER_APPLICABLE: [],
                            StatusUpdateType.MOLECULE_SIGNATURE_GENERATED: [],
                            StatusUpdateType.MOLECULE_IDENTITY_MATCH: [],
                            StatusUpdateType.MOLECULE_SIMILARITY_MATCH: [],
                            StatusUpdateType.MOLECULE_NO_MATCH: False
                        })
                    )
                case StatusUpdateType.SIMULATION_SUCCESSFULLY_LOADED:
                    self._status[simulation_file].loaded = True
                case StatusUpdateType.SIMULATION_RESOLUTION_FINISHED:
                    self._status[simulation_file].finished_successfully = True
                case StatusUpdateType.SIMULATION_RESOLUTION_TIMEOUT:
                    self._status[simulation_file].finished_timeout = True
                case StatusUpdateType.MOLECULE_RESOLUTION_STARTED:
                    self._status[simulation_file].molecules[update_data][
                        StatusUpdateType.MOLECULE_RESOLUTION_STARTED] = True
                case StatusUpdateType.MOLECULE_RESOLVER_APPLICABLE:
                    molecule_name, resolver_name = update_data
                    self._status[simulation_file].molecules[molecule_name][
                        StatusUpdateType.MOLECULE_RESOLVER_APPLICABLE].append(resolver_name)
                case StatusUpdateType.MOLECULE_SIGNATURE_GENERATED:
                    molecule_name, resolver_name = update_data
                    self._status[simulation_file].molecules[molecule_name][
                        StatusUpdateType.MOLECULE_SIGNATURE_GENERATED].append(resolver_name)
                case StatusUpdateType.MOLECULE_IDENTITY_MATCH:
                    molecule_name, resolver_name = update_data
                    self._status[simulation_file].molecules[molecule_name][
                        StatusUpdateType.MOLECULE_IDENTITY_MATCH].append(resolver_name)
                case StatusUpdateType.MOLECULE_SIMILARITY_MATCH:
                    molecule_name, resolver_name = update_data
                    self._status[simulation_file].molecules[molecule_name][
                        StatusUpdateType.MOLECULE_SIMILARITY_MATCH].append(resolver_name)
                case StatusUpdateType.MOLECULE_NO_MATCH:
                    self._status[simulation_file].molecules[update_data][StatusUpdateType.MOLECULE_NO_MATCH] = True
                case _:
                    log.error(
                        "ResolverManager encountered unexpected resolution update message: "
                        f"origin='{simulation_file}', kind='{update_type}', data='{update_data}'"
                    )

    @staticmethod
    def _update_results_file_molecule(data: dict, update: SimulationResolutionStatus, molecule_name: str):
        """
        Includes molecular resolution log in results file
        :param data: currently open results file
        :param update: log of processed simulation
        :param molecule_name: of interest
        :return: None
        """
        if molecule_name not in data['segment_resolution']:
            data['segment_resolution'][molecule_name] = {'errors': [LogError.ERROR_RESOLUTION_NOT_STARTED.value]}
            return

        if not update.molecules[molecule_name][StatusUpdateType.MOLECULE_RESOLUTION_STARTED] \
                and data['segment_resolution'][molecule_name]['resolver'] is None:
            data['segment_resolution'][molecule_name]['errors'].append(LogError.ERROR_RESOLUTION_NOT_STARTED.value)
            return

        data['segment_resolution'][molecule_name]['resolvers_attempted'] = update.molecules[molecule_name][
            StatusUpdateType.MOLECULE_RESOLVER_APPLICABLE]
        # data['segment_resolution'][molecule_name]['resolvers_signatured'] = update.molecules[molecule_name][
        #     StatusUpdateType.MOLECULE_SIGNATURE_GENERATED]
        data['segment_resolution'][molecule_name]['resolvers_identity'] = update.molecules[molecule_name][
            StatusUpdateType.MOLECULE_IDENTITY_MATCH]
        data['segment_resolution'][molecule_name]['resolvers_similarity'] = update.molecules[molecule_name][
            StatusUpdateType.MOLECULE_SIMILARITY_MATCH]

        if update.molecules[molecule_name][StatusUpdateType.MOLECULE_NO_MATCH]:
            data['segment_resolution'][molecule_name]['errors'].append(LogError.ERROR_NO_MATCH.value)

    def update_results_files(self) -> None:
        """
        Includes resolution log in results files
        :return: None
        """
        for update in self._status.values():
            if not os.path.exists(update.results_path):
                with open(update.results_path, 'w') as fd:
                    json.dump({'errors': [LogError.ERROR_NO_LOAD.value]}, fd)
            with open(update.results_path, 'r') as fd:
                data: dict = json.load(fd)

            if not update.loaded: data['errors'].append(LogError.ERROR_NO_LOAD.value)
            if update.finished_timeout: data['errors'].append(LogError.ERROR_TIMEOUT.value)
            if not update.finished_successfully: data['errors'].append(LogError.ERROR_NOT_FINISHED.value)

            for molecule_name in update.molecules:
                StatusQueue._update_results_file_molecule(data, update, molecule_name)

            with open(update.results_path, 'w') as fd:
                json.dump(data, fd, indent=4)

    def _aggregate_log(self) -> tuple:
        """
        Aggregates basic counts from the log
        :return: please see source code and return statement
        """
        sim_started = sim_loaded = sim_finished = sim_timeout = 0
        mol_started = mol_signatured = mol_matched = mol_partilly_matched = mol_unmatched = 0
        for update in self._status.values():
            sim_started += 1
            if update.loaded: sim_loaded += 1
            if update.finished_timeout or update.finished_successfully or not update.loaded: sim_finished += 1
            if update.finished_timeout: sim_timeout += 1
            for molecule in update.molecules:
                if update.molecules[molecule][StatusUpdateType.MOLECULE_RESOLUTION_STARTED]: mol_started += 1
                if update.molecules[molecule][StatusUpdateType.MOLECULE_SIGNATURE_GENERATED]: mol_signatured += 1
                if update.molecules[molecule][StatusUpdateType.MOLECULE_IDENTITY_MATCH]: mol_matched += 1
                if update.molecules[molecule][StatusUpdateType.MOLECULE_SIMILARITY_MATCH]: mol_partilly_matched += 1
                if update.molecules[molecule][StatusUpdateType.MOLECULE_NO_MATCH]: mol_unmatched += 1
        return sim_started, sim_loaded, sim_finished, sim_timeout, mol_started, mol_signatured, mol_matched, \
            mol_partilly_matched, mol_unmatched

    def emit_process_log(self) -> None:
        sim_started, sim_loaded, sim_finished, sim_timeout, \
            mol_started, mol_signatured, mol_matched, mol_partilly_matched, mol_unmatched = self._aggregate_log()
        log.info(f"simulations: {sim_started} started -> {sim_loaded} loaded -> {sim_finished} finished "
                 f"(of which {sim_timeout} timeouted)")
        log.info(f"molecules: {mol_started} loaded -> {mol_signatured} generated signature -> "
                 f"{mol_matched} identity / {mol_partilly_matched} similarity / {mol_unmatched} no hit")
