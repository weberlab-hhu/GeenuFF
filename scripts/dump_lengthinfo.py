#! /usr/bin/env python3
from geenuff.applications.exporters.lengths import LengthExportController
from geenuff.applications.exporter import MODES, RangeArgParser


class LengthArgParser(RangeArgParser):
    def __init__(self):
        super().__init__()
        self.parser.add_argument('--stats-only', action="store_true",
                                 help="output summary statistics about the lengths instead of the lengths themselves")


def main(args):
    controller = LengthExportController(args.db_path_in, args.longest)
    if args.mode in MODES:
        controller.prep_ranges(args.genomes, args.exclude_genomes,
                               MODES[args.mode])
    else:
        raise NotImplementedError("Requested mode ({}) not in implemented types {}".format(args.mode,
                                                                                           list(MODES.keys())))
    if args.stats_only:
        controller.write_length_stats(args.out)
    else:
        controller.write_lengths(args.out)


if __name__ == "__main__":
    parser = LengthArgParser()
    parser.parse_args()
    main(parser.args)
