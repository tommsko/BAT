import configparser
import json
import logging
import os.path
import random
import string
import subprocess
import zipfile
from dataclasses import dataclass
from typing import IO, Callable

from Bio import Blast

from ..base.fragment_resolver import IdentifierResolutionError

log = logging.getLogger("resolver")


def _mk_random_filename(extension: str, dir_path: str, name_len: int = 20) -> str:
    """
    Generates a random filename with given extension.
    :param extension: of the file, without the dot
    :param dir_path: directory where the file should be saved
    :param name_len: length of the name (excluding extension)
    :return: path to such file
    """
    extension = extension.replace(".", "")  # because of course we forget to

    # https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits
    filename: str = "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(name_len)
    )

    return os.path.abspath(os.path.join(dir_path, f"{filename}.{extension}"))


def _run_blast_online(
    sequence: str,
    tmp_dir_path: str,
    config: configparser.ConfigParser,
    database: str = "pdb",
    is_nucleotide: bool = False,
) -> list[dict] | None:
    """
    Executes blast protein search using remote service (NCBI)
    :param sequence: of the protein to be searched for
    :param tmp_dir_path: where to save temporary files (files will be cleaned up)
    :param config: main configuration
    :param database: blast database to use
    :return: list of jsons (results) returned, or None on failure
    :raises: IdentifierResolutionError on blast failure
    """

    Blast.email = config.get("BLAST", "remote_instance_mail")

    output_path: str = _mk_random_filename("zip", tmp_dir_path)

    # https://biopython.org/docs/dev/Tutorial/chapter_blast.html
    try:
        stream = Blast.qblast("blastp" if not is_nucleotide else "blastn", database, sequence, format_type="JSON2")
        with open(output_path, "wb") as out_stream:
            out_stream.write(stream.read())
        stream.close()

    except Exception as exc:
        raise IdentifierResolutionError("REMOTE blastp", f"blastp failed: {exc}")

    json_files: list[dict] = []

    archive: zipfile.ZipFile = zipfile.ZipFile(output_path)
    for file in archive.namelist():
        stream: IO = archive.open(file)
        json_files.append(json.loads(stream.read().decode()))
        stream.close()

    os.unlink(output_path)
    return json_files


def _run_blast_local(
    sequence: str, tmp_dir_path: str, config: configparser.ConfigParser, is_nucleotide: bool = False
) -> list[dict] | None:
    """
    Executes blast protein search using local engine
    :param sequence: of the protein to be searched for
    :param tmp_dir_path: where to save temporary files (files will be cleaned up)
    :param config: main configuration
    :param is_nucleotide: sequence contains nucleotides, not amino acids
    :return: list of jsons (results) returned, or None on failure
    :raises: IdentifierResolutionError on blast failure
    """

    output_path: str = _mk_random_filename("zip", tmp_dir_path)

    try:

        if len(sequence) < config.getint('BLAST', 'local_instance_run_short_search_if_under_X_residues'):
            short_search: list[str] = ['-task', 'blastp-short'] if not is_nucleotide else ['-task', 'blastn-short']
        else:
            short_search: list[str] = []

        result = subprocess.run(
            [
                config.get("BLAST", "local_instance_blastp_exec") if not is_nucleotide
                else config.get("BLAST", "local_instance_blastn_exec"),
                "-outfmt",
                "15",
                "-db",
                config.get("BLAST", "local_instance_database_prot") if not is_nucleotide
                else config.get("BLAST", "local_instance_database_nt"),
                "-out",
                output_path,
            ] + short_search,
            input=sequence,
            text=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        assert result.returncode == 0, "non-zero exit code"
    except Exception as exc:
        raise IdentifierResolutionError("LOCAL blastp" if not is_nucleotide else "LOCAL blastn", f"blast failed: {exc}")

    with open(output_path, "r") as fd:
        json_files: list[dict] = [json.load(fd)]

    os.unlink(output_path)
    return json_files


@dataclass
class BlastHit:
    accessions: list[str]
    evalue: float
    bits: float
    identity_p: float
    coverage_p: float


def _process_blast_results(json_files: list[dict]) -> list[BlastHit]:
    """
    Parses the blast results and standardizes the format of hits
    :param json_files: output of _run_blast*
    :return: list of hits
    """

    if json_files is None:
        return []

    results: list[BlastHit] = []

    for data in json_files:
        if "BlastOutput2" not in data:
            continue

        blast: dict = data["BlastOutput2"]
        if isinstance(blast, list):
            blast = blast[0]

        if (
            "report" not in blast
            or "results" not in blast["report"]
            or "search" not in blast["report"]["results"]
            or "hits" not in blast["report"]["results"]["search"]
        ):
            continue

        for hit in blast["report"]["results"]["search"]["hits"]:
            accessions: list[str] = [match["accession"] for match in hit["description"]]

            query_len: int = blast["report"]["results"]["search"]["query_len"]
            hit_len: int = hit["len"]

            _scores_section: dict = hit["hsps"][0]
            bits: float = _scores_section["bit_score"]
            evalue: float = _scores_section["evalue"]
            aligned_residues: int = _scores_section["align_len"]
            identity_p: float = _scores_section["identity"] / aligned_residues
            coverage_p: float = aligned_residues / query_len
            results.append(BlastHit(accessions, evalue, bits, identity_p, coverage_p))
    return results


def _run_blast(
    sequence: str, tmp_path_dir: str, config: configparser.ConfigParser, is_nucleotide: bool = False
) -> list[BlastHit] | None:
    """
    Wrapper to run either local or online blast search based on configuration
    :param sequence: of the protein to be searched for
    :param tmp_path_dir: where to save temporary files (files will be cleaned up)
    :param config: main configuration
    :param is_nucleotide: sequence contains nucleotides, not amino acids
    :return: list of BlastHits returned, or None on failure
    :raises: IdentifierResolutionError on blast failure (no hits is NOT a failure)
    """
    if not os.path.exists(tmp_path_dir):
        os.makedirs(tmp_path_dir)

    if config.getboolean("BLAST", "use_local_instance", fallback=False):
        return _process_blast_results(_run_blast_local(sequence, tmp_path_dir, config, is_nucleotide))
    else:
        return _process_blast_results(_run_blast_online(sequence, tmp_path_dir, config,
                                                        database="pdb" if not is_nucleotide else "pdbnt",
                                                        is_nucleotide=is_nucleotide))


def _filter_blast_results(
    results: list[BlastHit], by=Callable[[BlastHit], bool]
) -> list[BlastHit]:
    """
    Simple filter wrapper for blast results
    :param results: blast hits to filter
    :param by: lambda defining the filter (keep) condition
    :return: blast hits following the filter
    """
    return [hit for hit in results if by(hit)]


def run_blast_identity(sequence: str, config: configparser.ConfigParser, is_nucleotide: bool = False) -> list[str]:
    """
    Executes identity search with BLAST
    :param sequence: to be searched
    :param config: main configuration
    :param is_nucleotide: sequence contains nucleotides, not amino acids
    :return: accessions of identical proteins
    :raises IdentifierResolutionError: on blast failure (finding no results is NOT a failure)
    """
    hits: list[BlastHit] = _run_blast(
        sequence,
        config.get("BLAST", "temporary_directory", fallback=os.getcwd()),
        config,
        is_nucleotide
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.evalue
        <= config.getfloat("BLAST", "similarity_match_max_evalue"),
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.bits >= config.getint("BLAST", "similarity_match_min_bits"),
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.identity_p
        >= config.getfloat("BLAST", "similarity_match_min_identity_p"),
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.coverage_p
        >= config.getfloat("BLAST", "similarity_match_min_query_cov_p"),
    )
    return [accession for hit in hits for accession in hit.accessions]


def run_blast_similarity(
    sequence: str, config: configparser.ConfigParser, is_nucleotide: bool = False
) -> list[tuple[str, float]]:
    """
    Executes similarity search with BLAST
    :param sequence: to be searched
    :param config: main configuration
    :param is_nucleotide: sequence contains nucleotides, not amino acids
    :return: accessions of similar proteins with identity percentage
    :raises IdentifierResolutionError: on blast failure (finding no results is NOT a failure)
    """
    hits: list[BlastHit] = _run_blast(
        sequence,
        config.get("BLAST", "temporary_directory", fallback=os.getcwd()),
        config,
        is_nucleotide
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.evalue <= config.getfloat("BLAST", "exact_match_max_evalue"),
    )
    hits = _filter_blast_results(
        hits, lambda hit: hit.bits >= config.getint("BLAST", "exact_match_min_bits")
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.identity_p
        >= config.getfloat("BLAST", "exact_match_min_identity_p"),
    )
    hits = _filter_blast_results(
        hits,
        lambda hit: hit.coverage_p
        >= config.getfloat("BLAST", "exact_match_min_query_cov_p"),
    )
    return [(accession, hit.identity_p) for hit in hits for accession in hit.accessions]


def run_blast_hit_count(
    sequences: list[str], config: configparser.ConfigParser, is_nucleotide: bool = False
) -> dict[str, int]:
    """
    Runs blast on all the sequences and reports how many hits each of them has
    (do not even attempt to call this without local blast instance)
    :param sequences: to analyse
    :param config: main configuration
    :param is_nucleotide: sequence contains nucleotides, not amino acids
    :return: number of hits for each sequence
    """
    log.debug(f"starting batch blast analysis for {len(sequences)} sequences...")

    results: dict[str, int] = {}
    for i, sequence in enumerate(sequences):
        hits: list[BlastHit] = _run_blast(
            sequence,
            config.get("BLAST", "temporary_directory", fallback=os.getcwd()),
            config,
            is_nucleotide
        )
        results[sequence] = len(hits)
        log.debug(
            f"... [{i + 1}/{len(sequences)}] done (seq='{sequence}', hits={len(hits)})"
        )
    return results
