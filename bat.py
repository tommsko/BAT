import argparse
import configparser

from argparse import RawTextHelpFormatter
from src.utils.file_loader import SimulationFile
from src.resolver_manager.resolver_manager import ResolverManager
from src.discovery.simulation_exporter import SimulationExporter
from src.stats.resolver_statistics import ResolutionStatistics
from src.resolvers.fingerprint.fingerprint_manager import register_signatures, list_signatures
from src.metadata import fetch_metadata_all
import logging
from src.main_utils import assert_arguments, print_mode_header, init_loggers, load_simulations_filtered

logger = logging.getLogger('base')
handler = logging.StreamHandler()
formatter = logging.Formatter('[%(process)s] %(asctime)s %(levelname)s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

parser = argparse.ArgumentParser(
                    prog='MD Simulation Molecule Resolver',
                    description='Resolves molecules in molecular dynamics simulations into standardized identifiers',
                    epilog='Tomas Pavlik, 2024',
                    formatter_class=RawTextHelpFormatter
)

parser.add_argument('-r', '--resolve', action='store_true', help="run in resolve mode (identifying molecules)")
parser.add_argument('-d', '--discover', action='store_true', help="run in discover mode (exporting molecules)")
parser.add_argument('-s', '--stats', action='store_true', help="run in stats mode (aggregating resolution results)")
parser.add_argument('-m', '--metadata', action='store_true', help="run in fetch metadata mode (download metadata)")
parser.add_argument('-lf', '--list-fingerprints', action='store_true', help="run in fingerprints mode (listing)")
parser.add_argument('-rf', '--register-fingerprints', action='store_true', help="run in fingerprints mode (registration)")
parser.add_argument('-i', '--in',       help="[resolve]      input directory with MD simulations\n"
                                             "[discover]     input directory with MD simulations\n"
                                             "[stats]        input directory with MD simulations\n"
                                             "[signatures]   directory with generated signatures")
parser.add_argument('-i2', '--in2',     help="[resolve]      <unused>\n"
                                             "[discover]     <unused>\n"
                                             "[stats]        directory with already-made resolutions\n"
                                             "[signatures]   <unused>")
parser.add_argument('-o', '--out',      help="[resolve]      output directory where to write resolution results\n"
                                             "[discover]     output directory where to save debug exports\n"
                                             "[stats]        <unused>\n"
                                             "[signatures]   json file to save signatures in --list-signatures")
parser.add_argument('-o2', '--out2',    help="[resolve]      <unused>\n"
                                             "[discover]     <unused>\n"
                                             "[stats]        output file to save statistics [stats.json]\n"
                                             "[signatures]   <unused>")
parser.add_argument('-v', '--verbose', action='store_true', help="make everything more talkative [False]")
parser.add_argument('-vv', '--very_verbose', action='store_true', help="or even more more talkative [False]")
parser.add_argument('-c', '--config', help="path to the configuration [config.ini]")
parser.add_argument('-t', '--threads', help="number of threads available [1]")
parser.add_argument('-mo', '--metadata_overwrite', action='store_true', help="overwrites already fetched metadata [False]")
parser.add_argument('-fs','--filter_simulation', action='append', help='Filter simulations to process')
parser.add_argument('-fm','--filter_molecule', action='append', help='Filter molecules to process')
parser.add_argument('-xxx', '--xxx', action='store_true', help="debug")
args = parser.parse_args()

config: configparser.ConfigParser = configparser.ConfigParser()
config.read("config.ini" if args.config is None else args.config)

if args.resolve:
    print_mode_header(mode_name='resolve',
                      mode_description='Tries to identify biomolecules in molecular dynamics simulations',
                      argument_names=['in', 'out', 'threads'],
                      arguments=args)
    assert_arguments(args, {
        'in': 'path to the directory with MD simulations',
        'out': 'path to the directory where to write resolution results',
    })
    init_loggers(args.verbose, args.very_verbose)

    simulations: list[SimulationFile] = load_simulations_filtered(getattr(args, 'in'),
                                                                  args.filter_simulation,
                                                                  args.filter_molecule,
                                                                  config)
    resolver = ResolverManager(simulations,
                               getattr(args, 'out'),
                               1 if args.threads is None else int(args.threads),
                               config)
    resolver.run()
elif args.discover:
    print_mode_header(mode_name='discover',
                      mode_description='Tries to export molecule from molecular dynamics simulation '
                                       'to a human-readable formats',
                      argument_names=['in', 'out'],
                      arguments=args)
    assert_arguments(args, {
        'in': 'input directory with MD simulations',
        'out': 'output directory where to save debug exports (must be empty)',
    })
    init_loggers(args.verbose, args.very_verbose)
    simulations: list[SimulationFile] = load_simulations_filtered(getattr(args, 'in'),
                                                                  args.filter_simulation,
                                                                  args.filter_molecule,
                                                                  config)

    exporter: SimulationExporter = SimulationExporter(simulations, getattr(args, 'out'), config)
    exporter.try_export()
elif args.stats:
    print_mode_header(mode_name='statistics',
                      mode_description='Aggregates resolutions from resolve mode and generates basic statistics',
                      argument_names=['in', 'in2', 'out'],
                      arguments=args)
    assert_arguments(args, {
        'in': 'input directory with MD simulations (--in during resolve)',
        'in2': 'directory with already-made resolutions (--out during resolve)',
    })
    init_loggers(args.verbose, args.very_verbose)

    simulations: list[SimulationFile] = load_simulations_filtered(getattr(args, 'in'),
                                                                  args.filter_simulation,
                                                                  args.filter_molecule,
                                                                  config)
    stats: ResolutionStatistics = ResolutionStatistics(simulations, getattr(args, 'in2'), config)
    out_path: str = "stats.json" if args.out is None else args.out
    stats.save_stats_json(out_path)
    stats.print_stats_stdout()
elif args.list_fingerprints:
    print_mode_header(mode_name='fingerprints (list)',
                      mode_description='Aggregates and lists all the fingerprints generated during resolution',
                      argument_names=['in', 'out'],
                      arguments=args)
    assert_arguments(args, {
        'in': 'directory with generated and unprocessed signatures',
        'out': 'json file to save aggregated signatures to',
    })
    init_loggers(args.verbose, args.very_verbose)

    list_signatures(getattr(args, 'in'), getattr(args, 'out'))
elif args.register_fingerprints:
    print_mode_header(mode_name='fingerprints (registration)',
                      mode_description='Aggregates fingerprints and starts a simple registration workflow',
                      argument_names=['in'],
                      arguments=args)
    assert_arguments(args, {
        'in': 'directory with generated and unprocessed signatures',
    })
    init_loggers(args.verbose, args.very_verbose)

    register_signatures(getattr(args, 'in'), config)
elif args.metadata:
    print_mode_header(mode_name='metadata',
                      mode_description='Downloads metadata for successfully identified molecules',
                      argument_names=['in', 'out', 'metadata_overwrite'],
                      arguments=args)
    assert_arguments(args, {
        'in': 'path to the directory with MD simulations',
        'out': 'path to the directory where to write resolution results',
    })
    init_loggers(args.verbose, args.very_verbose)

    simulations: list[SimulationFile] = load_simulations_filtered(getattr(args, 'in'),
                                                                  args.filter_simulation,
                                                                  args.filter_molecule,
                                                                  config)
    fetch_metadata_all(simulations, getattr(args, 'out'), config, getattr(args, 'metadata_overwrite'))
else:
    raise RuntimeError("Unknown mode set, see --help")
