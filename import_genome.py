#!/usr/bin/env python3
import os
import sys
import logging
import argparse

from geenuff.applications.importer import ImportController
from geenuff.base import orm


class PathFinder(object):
    INPUT = 'input'
    OUTPUT = 'output'

    def __init__(self, db_path, basedir, fasta=None, gff=None, logfile=None):
        # directories
        self.db_out = db_path
        self.basedir = basedir
        self.input = '{}/{}/'.format(self.basedir, PathFinder.INPUT)
        self.output = '{}/{}/'.format(self.basedir, PathFinder.OUTPUT)
        if args.basedir is not None:
            for dir in [self.basedir, self.input, self.output]:
                os.makedirs(dir, exist_ok=True)
            self.problems_out = '{}import.log'.format(self.output)
        # files
        self.fasta_in = self._get_fa(fasta)
        self.gff_in = self._get_gff(gff)
        if not self.db_out:
            self.db_out = '{}geenuff.sqlite3'.format(self.output)
        if logfile is not None:
            self.problems_out = logfile

    def _get_fa(self, provided):
        if provided is not None:
            return provided
        maybe = os.listdir(self.input)
        # todo, actual file type detection
        maybe = [x for x in maybe if (x.endswith('.fa') or x.endswith('.fasta'))]
        self._confirm_exactly_one(maybe, 'fasta')
        return self.input + maybe[0]

    def _get_gff(self, provided):
        if provided is not None:
            return provided
        maybe = os.listdir(self.input)
        maybe = [x for x in maybe if (x.endswith('.gff') or x.endswith('.gff3'))]
        self._confirm_exactly_one(maybe, 'gff')
        return self.input + maybe[0]

    @staticmethod
    def _confirm_exactly_one(possibilities, info):
        assert len(possibilities) == 1, 'no(n) unique {} file found as input. Found: {}'.format(
                                            info, possibilities)


def main(args):
    if args.basedir is None:
        assert all([x is not None for x in [args.fasta, args.gff3, args.db_path, args.log_file]]), \
            "if basedir is none, all three custom input/output files must be manually specified " \
            "with --gff3, --fasta, --log-file and --db_path parameters"
    paths = PathFinder(args.db_path, args.basedir, fasta=args.fasta, gff=args.gff3,
                       logfile=args.log_file)

    msg_fmt_str = '%(asctime)s - %(levelname)s: %(message)s'
    date_fmt_str = '%d-%b-%y %H:%M:%S'
    logging.basicConfig(filename=paths.problems_out,
                        filemode='w',
                        level=logging.INFO,
                        format=msg_fmt_str,
                        datefmt=date_fmt_str)
    # log to file and stderr simultaneously
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(logging.Formatter(fmt=msg_fmt_str, datefmt=date_fmt_str))
    logging.getLogger().addHandler(stdout_handler)

    controller = ImportController(database_path=paths.db_out, replace_db=args.replace_db)
    genome_args = {}
    for key in ['species', 'accession', 'version', 'acquired_from']:
        genome_args[key] = vars(args)[key]
    controller.add_genome(paths.fasta_in, paths.gff_in, genome_args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--basedir', help='organized output (& input) directory. If this is not set, all four custom'
                        'input parameters must be set.')

    custominput = parser.add_argument_group('Override default with custom input/output location:')
    custominput.add_argument('--gff3', help='gff3 formatted file to parse / standardize')
    custominput.add_argument('--fasta', help='fasta file to parse standardize')
    custominput.add_argument('--db-path', help='path of the GeenuFF database')
    custominput.add_argument('--log-file', help="output path for import log (default basedir/output/import.log)")

    parser.add_argument('--replace-db', action='store_true',
                        help=('whether to override a GeenuFF database found at '
                              'the default location or at the location of --db_path'))

    genome_attr = parser.add_argument_group('Possible genome attributes:')
    genome_attr.add_argument('--species', required=True, help='name of the species')
    genome_attr.add_argument('--accession', default='', help='')
    genome_attr.add_argument('--version', default='', help='genome version')
    genome_attr.add_argument('--acquired-from', default='', help='genome source')

    args = parser.parse_args()

    main(args)
