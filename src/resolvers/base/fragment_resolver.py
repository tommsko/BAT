import configparser
import logging
from enum import Enum

import MDAnalysis

from .resolver_database import ResolverCache
from ...utils.file_loader import SimulationFile
from ...utils.timeout import timeout

log = logging.getLogger("resolver")


class SignatureGenerationError(Exception):
    def __init__(self, generation_step, cause):
        super().__init__(f"[RES.SIGN] Error while {generation_step}: {cause}")
        self.generation_step = generation_step
        self.cause = cause


class IdentifierResolutionError(Exception):
    def __init__(self, generation_step, cause):
        super().__init__(f"[RES.IDENT] Error while {generation_step}: {cause}")
        self.generation_step = generation_step
        self.cause = cause


class ResolverKind(Enum):
    PRIMARY = 1  # identity match can be considered ground truth
    FALLBACK = 2  # identity match can be considered true with very high likelihood


class ResolverBase:
    """
    Abstract class for fragment resolver. The workflow for resolvers is as follows:
    input: MDAnalysis' fragment (one molecule)
     ... .is_applicable() - can this resolver parse this fragment? (a quick check, but signature generation might fail)
     ... .get_signature() - translates atomic fragment into its string signature
     ... .get_identifiers() - attempts to identity match signature into standardized identifier(s)
     ... .get_similar_identifiers() - attempts to similarity match signature into standardized identifier(s)
     ... .generate_debug_data() - generates representation of the atomic fragment in paradigm of the resolver
     !! Due to built-in timeout with signals, calling get_identifiers or get_similar_identifiers is not thread safe !!
     !! run only one concurrent instance in a single process !!
    """

    def __init__(
        self,
        resolver_name: str,
        resolver_kind: ResolverKind,
        ident_name: str,
        cache: ResolverCache,
        configuration: configparser.ConfigParser,
    ) -> None:
        """
        Initializes base resolver
        :param resolver_name: name of the resolver (for caching purposes)
        :param resolver_kind: type of the resolver
        :param ident_name: name of the identifier (for result purposes)
        :param cache: cache for the resolver
        :param configuration: global configuration for the application
        """
        self.resolver_name: str = resolver_name
        self.resolver_kind: ResolverKind = resolver_kind
        self.ident_name: str = ident_name
        self.cache: ResolverCache = cache
        self.config: configuration.ConfigParser = configuration

        self.allow_resolution: bool = False  # 'activates' the resolver
        self.allow_partial_resolution: bool = (
            False  # on non-identity hit, still returns partial results and scores
        )
        self.resolve_through_cache: bool = False  # allows using cache for resolution
        self.resolve_through_search: bool = (
            False  # allows using online or offline databases for resolution
        )
        self.allow_debug_output: bool = (
            False  # allows resolver during debug simulation export
        )
        self.fetch_idents_timeout_sec: int | None = (
            None  # maximum time for fetching identifiers, None for unlimited
        )
        self._load_config(configuration)

    def _load_config(self, configuration: configparser.ConfigParser) -> None:
        """
        Loads resolver settings based on configuration provided
        :param configuration: loaded ConfigParser() instance
        :return: None
        """
        if not configuration.has_section(self.resolver_name):
            log.error(
                f"Resolver {self.resolver_name} is not in the configuration! Disabled by default"
            )
            return

        if configuration.has_option(self.resolver_name, "enabled"):
            self.allow_resolution = configuration.getboolean(
                self.resolver_name, "enabled"
            )

        if configuration.has_option(self.resolver_name, "enabled_similarity"):
            self.allow_partial_resolution = configuration.getboolean(
                self.resolver_name, "enabled_similarity"
            )

        if configuration.has_option(self.resolver_name, "from_cache"):
            self.resolve_through_cache = configuration.getboolean(
                self.resolver_name, "from_cache"
            )

        if configuration.has_option(self.resolver_name, "from_search"):
            self.resolve_through_search = configuration.getboolean(
                self.resolver_name, "from_search"
            )

        if configuration.has_option(self.resolver_name, "allow_export"):
            self.allow_debug_output = configuration.getboolean(
                self.resolver_name, "allow_export"
            )

        if configuration.has_option(self.resolver_name, "timeout_sec"):
            self.fetch_idents_timeout_sec = configuration.getint(
                self.resolver_name, "timeout_sec"
            )
            if self.fetch_idents_timeout_sec == 0:
                self.fetch_idents_timeout_sec = None

    def is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the resolver can be used on the given fragment
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """
        if not self.allow_resolution:
            return False
        return self._is_applicable(simulation, fragment)

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the resolver can be used on the given fragment.
        To be implemented by resolvers
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """
        raise NotImplementedError

    def get_signature(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> str:
        """
        Translates fragment into its string fingerprint
        To be implemented by resolvers
        :param fragment: AtomGroup representing a molecule
        :param simulation: of which fragment is part of
        :return: string fingerprint of the fragment given the resolver
        :raises: SignatureGenerationError if fingerprint cannot be created
        """
        raise NotImplementedError

    def get_identifiers(self, signature: str) -> list[str]:
        """
        Fetches <<identity-level>> identifiers associated with the given fingerprint
        :param signature: string fingerprint of the fragment given the resolver
        :return: list of standardized identifiers of the fingerprint
        :raises: IdentifierNotFoundError if no identifier couldn't be found
        """
        if self.resolve_through_cache:
            cached_identifiers: list[str] = self.cache.fetch_identity_identifiers(
                self.resolver_name, signature
            )
            if cached_identifiers:
                return cached_identifiers

        if self.resolve_through_search:

            search_timeout_error = IdentifierResolutionError(
                "resolving identifiers", "fetching identity matches timed out"
            )
            with timeout(
                seconds=self.fetch_idents_timeout_sec, error=search_timeout_error
            ):
                fetched_identifiers: list[str] = self.try_fetch_identifiers_exact_match(
                    signature
                )

            if fetched_identifiers:
                self.cache.save_identity_identifiers(
                    self.resolver_name, signature, fetched_identifiers
                )
                return fetched_identifiers

        raise IdentifierResolutionError(
            "resolving identifiers",
            "no identifiers found for identity match of the fingerprint",
        )

    def get_similar_identifiers(self, signature: str) -> list[tuple[str, float]]:
        """
        Fetches similar identifiers associated with the given fingerprint
        :param signature: string fingerprint of the fragment given the resolver
        :return: list of standardized identifiers of the fingerprint and their scores
        :raises: IdentifierNotFoundError if no identifier couldn't be found or this functionality is disabled
        """
        if not self.allow_partial_resolution:
            raise IdentifierResolutionError(
                "resolving similar identifiers", "functionality is disabled"
            )

        if self.resolve_through_cache:
            cached_identifiers: list[tuple[str, float]] = (
                self.cache.fetch_similarity_identifiers(self.resolver_name, signature)
            )
            if cached_identifiers:
                return cached_identifiers

        if self.resolve_through_search:
            search_timeout_error = IdentifierResolutionError(
                "resolving identifiers", "fetching similarity matches timed out"
            )
            with timeout(
                seconds=self.fetch_idents_timeout_sec, error=search_timeout_error
            ):
                fetched_identifiers: list[tuple[str, float]] = (
                    self.try_fetch_identifiers_relaxed_match(signature)
                )

            if fetched_identifiers:
                return fetched_identifiers

            raise IdentifierResolutionError(
                "resolving identifiers",
                "no identifiers found for similar match of the fingerprint",
            )

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Attempts to resolve signature into standardized identifiers, requiring identity match of the fingerprint
        To be implemented by resolvers
        :param signature: to be translated into identifiers
        :return: list of standardized identifiers of the fingerprint
        :raises IdentifierResolutionError: on internal failure during fetching
        """
        raise NotImplementedError

    def try_fetch_identifiers_relaxed_match(
        self, signature: str
    ) -> list[tuple[str, float]]:
        """
        Attempts to resolve signature into standardized identifiers, based on similarity match(es)
        To be implemented by resolvers
        :param signature: to be translated into identifiers
        :return: list of standardized identifiers of the fingerprint and their scores
        :raises IdentifierResolutionError: on internal failure during fetching
        """
        raise NotImplementedError

    def generate_debug_data(
        self,
        simulation: SimulationFile,
        fragment: MDAnalysis.AtomGroup,
        out_path: str,
        out_name: str,
    ) -> bool:
        """
        If resolver is applicable and can generate signature,
        generates representation of the fragment AtomGroup in resolver's paradigm
        To be implemented by resolvers
        :param out_path: directory where to save generated representations
        :param out_name: basename of the file to save representations to (will be suffixed with kind of export)
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True if export went without problems
        """
        raise NotImplementedError
