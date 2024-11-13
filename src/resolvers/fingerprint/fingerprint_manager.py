import configparser
import json
from collections import Counter
from dataclasses import dataclass
from os import listdir
from os.path import isfile, join

from prompt_toolkit import prompt
from prompt_toolkit.completion import WordCompleter, Completer, Completion
from prompt_toolkit.validation import Validator
from rich import print

from .fingerprint_resolver import FingerprintResolver
from .. import get_resolvers
from ..base.fragment_resolver import ResolverBase
from ..base.resolver_database import ResolverCache


@dataclass
class Signature:
    hash_str: str
    molecule_hits: list[str]
    simulation_hits: list[str]


def _aggregate_signatures(signature_directory: str) -> list[Signature]:
    """
    Aggregates all signatures
    :param signature_directory: from where to read fingerprint files
    :return: list of signatures, sorted by number of hits
    """
    signatures: dict[str, Signature] = {}
    signature_files = [
        join(signature_directory, f)
        for f in listdir(signature_directory)
        if isfile(join(signature_directory, f)) and f.endswith(".txt")
    ]
    for path in signature_files:
        with open(path, "r") as fd:
            data = json.load(fd)
            _hash, _fragment_name, _file_name = (
                data["hash"],
                data["fragment_name"],
                data["simulation"],
            )
            if _hash not in signatures:
                signatures[_hash] = Signature(_hash, [], [])
            signatures[_hash].molecule_hits.append(_fragment_name)
            signatures[_hash].simulation_hits.append(_file_name)
    return [
        sig
        for _hash, sig in sorted(
            signatures.items(),
            key=lambda item: len(item[1].molecule_hits),
            reverse=True,
        )
    ]


def list_signatures(signature_directory: str, output_file: str) -> None:
    """
    Aggregates signatures, saves them into file and prints basic information to stdout
    :param signature_directory: from where to fetch signatures
    :param output_file: where to save signatures (json)
    :return: None
    """
    signatures: list[Signature] = _aggregate_signatures(signature_directory)

    with open(output_file, "w") as fd:
        data: dict = {
            signature.hash_str: {
                "hits": len(signature.molecule_hits),
                "molecules:": signature.molecule_hits,
                "files": signature.simulation_hits,
            }
            for signature in signatures
        }
        json.dump(data, fd)

    print("[cyan][b] === Signature information === [/b][/cyan]")
    print()

    for signature in signatures:
        print(f"[cyan][b]{signature.hash_str}[/b][/cyan]")
        print(f"\thits: {len(signature.molecule_hits)}")
        print(f"\tmolecules (unique names): {Counter(signature.molecule_hits)}")


def _list_possible_ident_types(
    cache: ResolverCache, configuration: configparser.ConfigParser
) -> list[str]:
    """
    Fetches all identifier types that resolvers can produce
    :param cache: ResolverCache (for initializing resolvers)
    :param configuration: configuration (for initializing resolvers)
    :return: list of all identifier types
    """
    resolvers: list[ResolverBase] = get_resolvers(cache, configuration)
    return [resolver.ident_name for resolver in resolvers]


def cli_fetch_ident_type(possible_ident_types: set[str]) -> str:
    """
    Prompts user to select identifier type
    :param possible_ident_types: allowed values
    :return: selected value
    """
    completer: WordCompleter = WordCompleter(list(possible_ident_types))
    validator = Validator.from_callable(
        lambda txt: txt in possible_ident_types,
        error_message="Select one of allowed identifier types",
        move_cursor_to_end=True,
    )
    return prompt("Identifier type: ", completer=completer, validator=validator)


def cli_fetch_text(prompt_text: str) -> str:
    """
    Fetches non-empty text from user
    :param prompt_text: to display to the user
    :return: None
    """
    validator = Validator.from_callable(
        lambda txt: len(txt) > 0,
        error_message="Non-empty text is required",
        move_cursor_to_end=True,
    )
    return prompt(f"{prompt_text}: ", validator=validator)


def cli_fetch_yes_no(action: str) -> bool:
    """
    Prompts user to select yes or no
    :param action: to be printed
    :return: bool, depending on what user selected
    """

    class YesNoCompleter(Completer):
        def get_completions(self, document, complete_event):
            yield Completion("yes", start_position=0, style="bg:green fg:white")
            yield Completion("no", start_position=0, style="bg:red fg:white")

    validator = Validator.from_callable(
        lambda txt: txt in ("yes", "no"),
        error_message="Select yes or no",
        move_cursor_to_end=True,
    )
    return (
        prompt(f"{action}: ", completer=YesNoCompleter(), validator=validator) == "yes"
    )


def register_signatures(
    signature_directory: str, configuration: configparser.ConfigParser
) -> None:
    """
    Simple script to process signatures and register them in order of the most significant
    :param signature_directory: from where to process signatures
    :param configuration: to initialize resolvers
    :return: None
    """
    cache: ResolverCache = ResolverCache(configuration)
    possible_ident_types: set[str] = set(
        _list_possible_ident_types(cache, configuration)
    )
    resolver: FingerprintResolver = FingerprintResolver(cache, configuration)
    signatures_to_resolve: list[Signature] = list(
        reversed(_aggregate_signatures(signature_directory))
    )  # to allow pop

    print("[cyan][b] === Signature registration === [/b][/cyan]")
    print()

    while signatures_to_resolve:
        signature = signatures_to_resolve.pop()
        if cache.fetch_identity_identifiers(resolver.RESOLVER_NAME, signature.hash_str):
            continue
        if not cli_fetch_yes_no("Continue registering signatures?"):
            break
        print(f"[cyan][b]{signature.hash_str}[/b][/cyan]")
        print(f"hits: {len(signature.molecule_hits)}")
        print(f"molecules (unique names): {Counter(signature.molecule_hits)}")

        if cli_fetch_yes_no("Register fingerprint?"):
            ident_type: str = cli_fetch_ident_type(possible_ident_types)
            ident_value: str = cli_fetch_text("Identifier")
            cache.save_identity_identifiers(
                resolver.RESOLVER_NAME,
                signature.hash_str,
                [f"{ident_type}:{ident_value}"],
            )
        print("\n\n")
    print("[cyan][b] === Signature fully completed === [/b][/cyan]")
