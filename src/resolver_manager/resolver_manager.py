import configparser
import logging
import multiprocessing
import signal
import sys
import time

import MDAnalysis

from .resolver_manager_log import StatusQueue, StatusUpdateType, emit_update, emit_update_skipped_molecules
from ..resolvers import get_resolvers
from ..resolvers.base.fragment_resolver import ResolverBase, IdentifierResolutionError, SignatureGenerationError
from ..resolvers.base.resolver_database import ResolverCache
from ..utils.file_loader import SimulationFile
from ..utils.results_file import ResultsFile

log = logging.getLogger("base")


class ResolverManager:
    def __init__(
            self,
            simulations: list[SimulationFile],
            output_directory: str,
            threads: int,
            configuration: configparser.ConfigParser,
    ) -> None:
        """
        Initializes the ResolverManager
        :param simulations: list of simulations to resolve
        :param output_directory: path to the directory where to save the results
        :param threads: number of threads to use for resolving
        """
        self.threads: int = threads

        self.tasks: dict[SimulationFile, ResultsFile] = {
            sim: ResultsFile(sim, output_directory) for sim in simulations
        }
        log.info(f"ResolverManager loaded {len(simulations)} simulations to process...")

        self.cache: ResolverCache = ResolverCache(configuration)
        log.info("ResolverManager cache loaded...")

        self.configuration: configparser.ConfigParser = configuration
        self.resolvers: list[ResolverBase] = get_resolvers(
            self.cache, self.configuration
        )
        log.info(f"ResolverManager initialized {len(self.resolvers)} resolvers...")

        self.SIMULATION_TIMEOUT_SEC: int = configuration.getint('resolution', 'timeout_sec_simulation')
        self.RETRY_FAILED_MOLECULES: bool = configuration.getboolean('resolution', 'retry_failed_resolution')

    def _run_multiprocess_collect(
            self,
            processes: list[multiprocessing.Process | None],
            processes_started_at: list[float | None],
            processes_simulation_tprs: list[str | None],
    ) -> None:
        """
        Collects finished processes
        :param processes: list of processes
        :param processes_started_at: list of times when processes started
        :param processes_simulation_tprs: list of paths to the tpr files being resolved per process
        :return: None
        """
        for i in range(self.threads):
            if processes[i] is None:
                continue
            processes[i].join(timeout=0)
            if not processes[i].is_alive():
                processes[i] = None
                processes_started_at[i] = 0
                processes_simulation_tprs[i] = None

    def _run_multiprocess_terminate_timeout(
            self,
            processes: list[multiprocessing.Process | None],
            processes_started_at: list[float | None],
            status_queue: StatusQueue,
            processes_simulation_tprs: list[str | None],
    ) -> None:
        """
        Terminate processes which have been running for too long
        :param processes: list of processes
        :param processes_started_at: list of times when processes started
        :param status_queue: queue where resolution dumps status
        :param processes_simulation_tprs: list of paths to the tpr files being resolved per process
        :return: None
        """
        for i in range(self.threads):
            if (
                    processes[i] is not None
                    and self.SIMULATION_TIMEOUT_SEC > 0
                    and time.time() - processes_started_at[i] > self.SIMULATION_TIMEOUT_SEC
            ):
                log.error(
                    f"Terminated resolution process (pid={processes[i].pid}). Reason: timeout"
                )
                emit_update(status_queue.get_queue(),
                            processes_simulation_tprs[i],
                            StatusUpdateType.SIMULATION_RESOLUTION_TIMEOUT,
                            "")
                processes[i].terminate()
                processes[i] = None
                processes_started_at[i] = 0
                processes_simulation_tprs[i] = None

    def _run_multiprocess_start(
            self,
            processes: list[multiprocessing.Process | None],
            processes_started_at: list[float | None],
            status_queue: StatusQueue,
            processes_simulation_tprs: list[str | None],
    ) -> None:
        """
        Starts new resolution threads given simulations to process provided in ctor
        :param processes: list of processes
        :param processes_started_at: list of times when processes started
        :param status_queue: queue where resolution dumps status
        :param processes_simulation_tprs: list of paths to the tpr files being resolved per process
        :return: None
        """
        for i in range(self.threads):
            if processes[i] is None and self.tasks:
                simulation, results = self.tasks.popitem()
                processes[i] = multiprocessing.Process(
                    target=self.resolve_simulation,
                    args=(simulation, results, status_queue.get_queue()),
                )
                processes[i].start()
                processes_started_at[i] = time.time()
                processes_simulation_tprs[i] = simulation.paths[0]
                log.debug(
                    f"Started resolution of '{simulation.paths[0]}' (pid={processes[i].pid})"
                )

    def _run_multiprocess_emit_log(self, processes: list[multiprocessing.Process | None]) -> None:
        """
        Emits log messages with basic information about resolution status
        :param processes: list of processes performing the resolution
        :return: None
        """
        workers_alive: int = sum(1 if x is not None else 0 for x in processes)

        log.info("~~~ Resolution in process ~~~")
        log.info(f"current queue: {len(self.tasks)} ({workers_alive} / {self.threads} workers alive)")

    def _run_multiprocess(self) -> StatusQueue:
        """
        Runs multiple concurrent processes resolving simulations provided in ctor
        :return: StatusQueue populated during the resolution
        """

        processes: list[multiprocessing.Process | None] = [
            None for _ in range(self.threads)
        ]
        processes_started_at: list[float | None] = [None for _ in range(self.threads)]
        processes_simulation_tprs: list[str | None] = [None for _ in range(self.threads)]
        status: StatusQueue = StatusQueue(self.threads)

        def signal_handler(sig, frame):
            log.info("Recieved SIGINT, terminating all processes")
            for i in range(self.threads):
                if processes[i] is not None and processes[i].is_alive():
                    processes[i].terminate()
            sys.exit(0)

        signal.signal(signal.SIGINT, signal_handler)

        while any(p is not None for p in processes) or self.tasks:
            status.process_queue()
            self._run_multiprocess_emit_log(processes)
            status.emit_process_log()
            self._run_multiprocess_collect(processes, processes_started_at, processes_simulation_tprs)
            self._run_multiprocess_terminate_timeout(processes, processes_started_at, status, processes_simulation_tprs)
            self._run_multiprocess_start(processes, processes_started_at, status, processes_simulation_tprs)
            time.sleep(1)  # resolutions may take minutes to finish, 1 second refresh time is good enough
        return status

    def run(self) -> StatusQueue:
        """
        Start resolution of all the simulations provided in ctor
        :return: log of the resolution
        """
        log.info(f"ResolverManager starting with {self.threads} threads...")
        stats: StatusQueue = self._run_multiprocess()
        stats.update_results_files()
        return stats

    def resolve_simulation(
            self,
            simulation: SimulationFile,
            results: ResultsFile,
            status: multiprocessing.Queue,
    ) -> None:
        """
        Attempts to resolve the molecules in simulation
        :param simulation: to be resolved
        :param results: ResultsFile linked to the simulation
        :param status: queue where to put updates about the resolution
        :return: None
        """
        log.debug(f"Trying to resolve simulation '{simulation.simulation_name}'")
        simulation_tpr: str = simulation.paths[0]
        emit_update(status, simulation_tpr, StatusUpdateType.SIMULATION_RESOLUTION_STARTED, results.path)

        if simulation.load_simulation() and results.load():
            emit_update(status, simulation_tpr, StatusUpdateType.SIMULATION_SUCCESSFULLY_LOADED, "")
        else:
            return

        if results.is_fully_resolved() or (results.is_partially_resolved() and not self.RETRY_FAILED_MOLECULES):
            log.debug("Simulation has been already resolved, skipping!")
            emit_update_skipped_molecules(simulation, results, status, emit_unmatched=True)
            emit_update(status, simulation_tpr, StatusUpdateType.SIMULATION_RESOLUTION_FINISHED, "")
            simulation.unload_simulation()
            return

        # we are skipping identity and similarity matched molecules either way
        emit_update_skipped_molecules(simulation, results, status, emit_unmatched=False)

        for segment in simulation.list_segments():
            if results.is_segment_resolved(
                    segment
            ):  # partial resolution is good, not overwriting existing results
                log.debug(
                    f"[{simulation.simulation_name}/{segment}] already resolved, skipping..."
                )
            else:
                ResolverManager.resolve_molecule(
                    simulation, segment, results, self.resolvers, status
                )

        emit_update(status, simulation_tpr, StatusUpdateType.SIMULATION_RESOLUTION_FINISHED, "")
        results.save()
        simulation.unload_simulation()

    @staticmethod
    def _resolve_molecule_try_identity(
            resolver: ResolverBase,
            signature: str,
            molecule_name: str,
            molecule_results: dict,
            results: ResultsFile,
            simulation: SimulationFile,
            status: multiprocessing.Queue,
    ) -> bool:
        """
        Attempts to resolve identifiers for identity match of the molecule given
        :param resolver: which produced signature below;
        :param signature: of the molecule to be identified
        :param molecule_name: of the molecule to be identified
        :param molecule_results: sub-dict of the results for given molecule
        :param results: ResultsFile associated with current resolution
        :param simulation: of which the molecule is part of
        :param status: output queue for resolution updates
        :return: True if successful, False otherwise
        """

        try:
            identifiers = resolver.get_identifiers(signature)
            log.debug(
                f"[{simulation.simulation_name}/{molecule_name}] '{resolver.resolver_name}' "
                f"resolved on identity level. Identifiers: {identifiers}"
            )
            molecule_results["resolver"] = resolver.resolver_name
            molecule_results["ident_type"] = resolver.ident_name
            molecule_results["resolved_idents"] = identifiers
            emit_update(status, simulation.paths[0], StatusUpdateType.MOLECULE_IDENTITY_MATCH,
                        (molecule_name, resolver.resolver_name))
            results.save()
            return True
        except IdentifierResolutionError as exc:
            log.debug(
                f"[{simulation.simulation_name}/{molecule_name}] '{resolver.resolver_name}' "
                f"no identity matches found! ({exc})"
            )
        return False

    @staticmethod
    def _resolve_molecule_try_similarity(
            resolver: ResolverBase,
            signature: str,
            molecule_name: str,
            molecule_results: dict,
            results: ResultsFile,
            simulation: SimulationFile,
            status: multiprocessing.Queue,
    ) -> bool:
        """
        Attempts to resolve identifiers for similar matches of the molecule given
        :param resolver: which produced fingerprint below;
        :param signature: of the molecule to be identified
        :param molecule_name: of the molecule to be identified
        :param molecule_results: dictionary where to save the results
        :param results: ResultsFile associated with current resolution
        :param simulation: of which the molecule is part of
        :param status: output queue for resolution updates
        :return: True if successful, False otherwise
        """

        try:
            identifiers = resolver.get_similar_identifiers(signature)
            log.debug(
                f"[{simulation.simulation_name}/{molecule_name}] '{resolver.resolver_name}' "
                f"resolved on similarity level. Identifiers: {identifiers}"
            )
            molecule_results["resolver"] = resolver.resolver_name
            molecule_results["ident_type"] = resolver.ident_name
            molecule_results["possible_idents"] = identifiers
            emit_update(status, simulation.paths[0],
                        StatusUpdateType.MOLECULE_SIMILARITY_MATCH, (molecule_name, resolver.resolver_name))
            results.save()
            return True
        except IdentifierResolutionError as exc:
            log.debug(
                f"[{simulation.simulation_name}/{molecule_name}] '{resolver.resolver_name}' "
                f"no similarity matches found! ({exc})"
            )
        return False

    @staticmethod
    def _resolve_molecule_check_applicable(simulation: SimulationFile,
                                           molecule: MDAnalysis.AtomGroup,
                                           molecule_name: str,
                                           resolver: ResolverBase,
                                           status: multiprocessing.Queue) -> bool:
        """
        Checks whether resolver is applicable on given molecule, and if so, emits necessary messages
        :param simulation: containing the molecule
        :param molecule: to be resolved
        :param molecule_name: of the resolved molecule
        :param resolver: attempted for resolution
        :param status: queue for updates
        :return: True if resolver is applicable, False otherwise
        """
        if not resolver.is_applicable(simulation, molecule):
            return False
        log.debug(
            f"[{simulation.simulation_name}/{molecule_name}] attempting resolver '{resolver.resolver_name}'"
        )
        emit_update(status,
                    simulation.paths[0],
                    StatusUpdateType.MOLECULE_RESOLVER_APPLICABLE,
                    (molecule_name, resolver.resolver_name))
        return True

    @staticmethod
    def _resolve_molecule_mk_signature(simulation: SimulationFile,
                                       molecule: MDAnalysis.AtomGroup,
                                       molecule_name: str,
                                       resolver: ResolverBase,
                                       status: multiprocessing.Queue,
                                       molecule_results: dict,
                                       results: ResultsFile) -> str | None:
        """
        Attempts to generate signature of the molecule
        :param simulation: containing the molecule
        :param molecule: to be resolved
        :param molecule_name: of the resolved molecule
        :param resolver: attempted for resolution
        :param status: queue for updates
        :param molecule_results: sub-dict for results of a given molecule
        :param results: ResultsFile associated with current resolution
        :return: string signature
        """
        try:
            signature = resolver.get_signature(simulation, molecule)
            log.debug(
                f"[{simulation.simulation_name}/{molecule_name}] '{resolver.resolver_name}' signature: "
                f"'{signature if len(signature) < 20 else f'{signature[:20]}...'}'"
            )
            emit_update(status,
                        simulation.paths[0],
                        StatusUpdateType.MOLECULE_SIGNATURE_GENERATED,
                        (molecule_name, resolver.resolver_name))
            molecule_results['resolver_generated_signature'].append(resolver.resolver_name)
            results.save()
            return signature
        except SignatureGenerationError as exc:
            log.debug(
                f"[{simulation.simulation_name}/{molecule_name}] '{resolver.resolver_name}' failed to create signature"
            )
            log.debug(exc, exc_info=True)
            return None

    @staticmethod
    def resolve_molecule(
            simulation: SimulationFile,
            molecule_name: str,
            results: ResultsFile,
            resolvers: list[ResolverBase],
            status: multiprocessing.Queue,
    ) -> bool:
        """
        Tries to resolve the molecule into standardized identifier(s) using resolvers
        :param simulation: to be analysed
        :param molecule_name: to be resolved
        :param results: to store success/failure data
        :param resolvers: available resolvers from get_resolvers(...)
        :param status: queue where to put updates about the resolution
        :return: True if successful, False otherwise
        """
        simulation_tpr: str = simulation.paths[0]
        emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_RESOLUTION_STARTED, molecule_name)
        molecule_results: dict = results.data["segment_resolution"][molecule_name]

        for resolver in resolvers:
            try:
                molecule: MDAnalysis.AtomGroup = simulation.get_molecule(molecule_name)

                if not ResolverManager._resolve_molecule_check_applicable(simulation, molecule, molecule_name,
                                                                          resolver, status):
                    continue

                signature: str | None = ResolverManager._resolve_molecule_mk_signature(simulation, molecule,
                                                                                       molecule_name, resolver, status,
                                                                                       molecule_results, results)
                if signature is None:
                    continue

                if ResolverManager._resolve_molecule_try_identity(resolver, signature, molecule_name, molecule_results,
                                                                  results, simulation, status):
                    return True

                if ResolverManager._resolve_molecule_try_similarity(resolver, signature, molecule_name,
                                                                    molecule_results, results, simulation, status):
                    return True
                # break is implicit here (using returns)
            except Exception as exc:
                log.debug(
                    f"[{simulation.simulation_name}/{molecule_name}] '{resolver.resolver_name}' "
                    f"failed internally"
                )
                log.debug(exc, exc_info=True)
                continue

        emit_update(status, simulation_tpr, StatusUpdateType.MOLECULE_NO_MATCH, molecule_name)
        log.error(
            f"[{simulation.simulation_name}/{molecule_name}] resolution was not possible"
        )
        return False

