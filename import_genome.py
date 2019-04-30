import logging
import argparse
import os

from geenuff.applications.gffimporter import ImportController, CoordinateHandler
from geenuff.base import orm


class PathFinder(object):
    INPUT = 'input'
    OUTPUT = 'output'

    def __init__(self, basedir, fasta=None, gff=None):
        # directories
        self.basedir = basedir
        self.input = '{}/{}/'.format(self.basedir, PathFinder.INPUT)
        self.output = '{}/{}/'.format(self.basedir, PathFinder.OUTPUT)
        for dir in [self.basedir, self.input, self.output]:
            os.makedirs(dir, exist_ok=True)
        # files
        self.fasta_in = self._get_fa(fasta)
        self.gff_in = self._get_gff(gff)
        self.db_out = '{}geenuff.sqlite3'.format(self.output)
        self.problems_out = '{}import.log'.format(self.output)

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
        assert len(possibilities) == 1, 'no(n) unique {} file found as input. Found: {}'.format(info, possibilities)


def add_mers(args, controller):
    session = controller.session
    latest_genome = session.query(orm.Genome).order_by(orm.Genome.id.desc()).first()
    coords = session.query(orm.Coordinate).filter(orm.Genome.id == latest_genome.id).all()
    for coord in coords:
        coord_handler = CoordinateHandler(data=coord, controller=controller)
        coord_handler.add_mer_counts(args.min_k, args.max_k)


def main(args):
    logging.basicConfig(level=logging.WARNING)
    paths = PathFinder(args.basedir, fasta=args.fasta, gff=args.gff3)

    controller = ImportController(database_path=paths.db_out, err_path=paths.problems_out)
    controller.add_sequences(paths.fasta_in)
    if args.max_k > 0:
        if args.min_k < 1:
            args.min_k = 1
        add_mers(args, controller)
    controller.add_gff(paths.gff_in)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--basedir', help='organized output (& input) directory', required=True)
    custominput = parser.add_argument_group("Override default with custom input location:")
    custominput.add_argument('--gff3', help='gff3 formatted file to parse / standardize')
    custominput.add_argument('--fasta', help='fasta file to parse standardize')
    fasta_specific = parser.add_argument_group("Sequence meta_info customizable:")
    fasta_specific.add_argument('--min_k', help='minumum size kmer to calculate from sequence',
                                default=0, type=int)
    fasta_specific.add_argument('--max_k', help='maximum size kmer to calculate from sequence',
                                default=0, type=int)
    args = parser.parse_args()

    assert args.min_k <= args.max_k, 'min_k can not be greater than max_k'

    main(args)
