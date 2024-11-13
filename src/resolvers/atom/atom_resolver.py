import configparser
import json
import os.path

import MDAnalysis
from src.utils.file_loader import SimulationFile

from ..base.fragment_resolver import IdentifierResolutionError, SignatureGenerationError
from ..base.fragment_resolver import ResolverBase, ResolverKind
from ..base.resolver_database import ResolverCache
from ...constants import (
    ATOMTABLE_ELEMENTS_MASSES,
    RESOLVER_ATOM_APPROX_FIND_SIMILAR_N_MASSES,
)
from ...utils.repair_simulation import fix_missing_elements


class AtomResolver(ResolverBase):
    """
    AtomResolver looks at the molecule as a single atom (possibly ion)
    """

    RESOLVER_NAME, RESOLVER_NAME_UNDETERMINED = "atom", "atom_aprox"
    RESOLVER_IDENT, RESOLVER_IDENT_UNDETERMINED = "ELEMENT", "ELEMENT_APROX"

    def __init__(
        self,
        cache: ResolverCache,
        configuration: configparser.ConfigParser,
        try_fix_elements: bool = False,
        allow_undetermined_elements: bool = False,
    ) -> None:
        """
        Initializes the resolver of atomic elements
        :param cache: resolver cache
        :param configuration: for the resolver
        :param try_fix_elements: if atomic elements are not present, try to fix them using atomic masses
        :param allow_undetermined_elements: ignore atomic element and only predict possible elements given the mass
        """
        super().__init__(
            (
                AtomResolver.RESOLVER_NAME
                if not allow_undetermined_elements
                else AtomResolver.RESOLVER_NAME_UNDETERMINED
            ),
            ResolverKind.PRIMARY,
            (
                AtomResolver.RESOLVER_IDENT
                if not allow_undetermined_elements
                else AtomResolver.RESOLVER_IDENT_UNDETERMINED
            ),
            cache,
            configuration,
        )

        self.try_fix_elements: bool = try_fix_elements
        self.allow_undetermined_elements: bool = allow_undetermined_elements

    def _is_applicable(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> bool:
        """
        Decides whether the AtomResolver can be used on the given fragment
        :param fragment: AtomGroup representing a fragment (single atom)
        :param simulation: of which fragment is part of
        :return: True if resolver can be used on a given fragment
        """

        fragment_name: str = fragment.segids[0]

        # in order to resolve fragment as atomic element, it requires to have element information
        # as well as be a single atom. If undetermined element information is allowed, it requires only atomic mass

        if not self.allow_undetermined_elements:
            # firstly, atom element information must be present or determined
            if not simulation.supports_atomic_elements(fragment_name):
                if (
                    simulation.supports_atomic_masses(fragment_name)
                    and self.try_fix_elements
                    and fix_missing_elements(simulation, fragment)
                ):
                    ...  # ok!
                else:
                    return False
        else:
            if not simulation.supports_atomic_masses(fragment_name):
                return False

        # secondly, it must be only one atom (otherwise chem resolver is applicable)
        return fragment.n_atoms == 1

    def get_signature(
        self, simulation: SimulationFile, fragment: MDAnalysis.AtomGroup
    ) -> str:
        """
        Translates fragment into element name or element mass (see .allow_undetermined_elements)
        :param simulation: of which the fragment is part of
        :param fragment: AtomGroup representing a single atom
        :return: element of the atom
        :raises: SignatureGenerationError if fingerprint cannot be created
        """

        if not self.allow_undetermined_elements:
            element_code: str = fragment.elements[0].upper()
            if element_code not in [x.upper() for x in ATOMTABLE_ELEMENTS_MASSES]:
                raise SignatureGenerationError(
                    "extracting atom element",
                    f"atom element is invalid: '{element_code}'",
                )
            return element_code
        else:
            element_mass: str = str(fragment.masses[0])
            return element_mass

    def try_fetch_identifiers_exact_match(self, signature: str) -> list[str]:
        """
        Just returns the signature of element from get_signature(...).
        :param signature: to be translated into standardized element name
        :return: list of exactly one element name, given in the fingerprint
        :raises IdentifierResolutionError: if allow_undetermined_elements is True (no identity match in that mode)
        """
        if not self.allow_undetermined_elements:
            return [signature]
        else:
            raise IdentifierResolutionError(
                "resolving ATOM identity match",
                "not supported with .allow_undetermined_elements",
            )

    def try_fetch_identifiers_relaxed_match(
        self, signature: str
    ) -> list[tuple[str, float]]:
        """
        If allow_undetermined_elements is False, returns fingerprint and 100% confidence.
        Otherwise, tries to find similar elements given atomic mass.
        :param signature: to be translated into standardized element name
        :return: list of elements with similarity
        """

        if not self.allow_undetermined_elements:
            return [(signature, 100)]

        atomic_mass: float = float(signature)
        mass_differences: list[tuple[str, float, float]] = [
            (element, abs(atomic_mass - element_mass), element_mass)
            for element, element_mass in ATOMTABLE_ELEMENTS_MASSES.items()
        ]
        mass_differences.sort(key=lambda x: x[1])

        results: list[tuple[str, float]] = []
        for elem, _, mass in mass_differences[:RESOLVER_ATOM_APPROX_FIND_SIMILAR_N_MASSES]:
            similarity: float = atomic_mass / mass
            if similarity > 1:
                similarity = 1 / similarity  # symmetric
            results.append((elem, round(similarity * 100)))
        return results

    def generate_debug_data(
        self,
        simulation: SimulationFile,
        fragment: MDAnalysis.AtomGroup,
        out_path: str,
        out_name: str,
    ) -> bool:
        """
        Generates very basic information about the atom in the simulation
        :param out_path: directory where to save generated representations
        :param out_name: basename of the file to save representations to (will be suffixed with kind of export)
        :param simulation: of which fragment is part of
        :param fragment: AtomGroup representing a molecule
        :return: True if export went without problems
        """
        with open(os.path.join(out_path, f"{out_name}_ATOM_desc.json"), "w") as fd:
            json.dump(
                {
                    "atom_element": (
                        fragment.elements[0]
                        if hasattr(fragment, "elements") and len(fragment.elements) > 0
                        else None
                    ),
                    "atom_name": (
                        fragment.names[0]
                        if hasattr(fragment, "names") and len(fragment.names) > 0
                        else None
                    ),
                    "atom_mass": (
                        fragment.masses[0]
                        if hasattr(fragment, "masses") and len(fragment.masses) > 0
                        else None
                    ),
                    "atom_charge": (
                        fragment.charges[0]
                        if hasattr(fragment, "charges") and len(fragment.charges) > 0
                        else None
                    ),
                },
                fd,
            )
        return True
