#! /usr/bin/env python3
import argparse

from geenuff.applications.exporters.sequence import FastaExportController
from geenuff.applications.exporter import RangeArgParser, MODES


def main(args):
    controller = FastaExportController(args.db_path_in, args.longest)
    if args.mode in MODES:
        controller.prep_ranges(args.genomes, args.exclude_genomes,
                               MODES[args.mode])
    else:
        raise NotImplementedError("Requested mode ({}) not in implemented types {}".format(args.mode,
                                                                                           list(MODES.keys())))

    controller.write_fa(args.out)


if __name__ == "__main__":
    parser = RangeArgParser()
    parser.parse_args()
    main(parser.args)
