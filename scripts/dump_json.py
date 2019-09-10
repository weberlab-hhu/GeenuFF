#! /usr/bin/env python3
from geenuff.applications.exporters.json import JsonExportController
from geenuff.applications.exporter import ExportArgParser
from geenuff.base.helpers import strand_as_bool


class JsonArgParser(ExportArgParser):
    def __init__(self):
        super().__init__()
        self.parser.add_argument('--species', required=True,
                                 help="name of target species in database")
        self.parser.add_argument('--seqid', required=True,
                                 help="name of target sequence (chromosome / scaffold) in database")
        self.parser.add_argument('--strand', required=True,
                                 help="'+' or '-' for strand to query (defaults to both)")
        self.parser.add_argument('--start', default=0, type=int,
                                 help='start coordinate of query range (pythonic)')
        self.parser.add_argument('--end', default=None,
                                 help="end coordinate of query range (pythonic)")
        self.parser.add_argument('--pretty', action="store_true",
                                 help="set for multi-line indented json output")


def main(args):
    controller = JsonExportController(args.db_path_in, args.longest)
    is_plus_strand = strand_as_bool(args.strand)
    controller.query_and_write(species=args.species,
                               seqid=args.seqid,
                               start=args.start,
                               end=args.end,
                               is_plus_strand=is_plus_strand,
                               file_out=args.out,
                               pretty=args.pretty)


if __name__ == "__main__":
    parser = JsonArgParser()
    parser.parse_args()
    main(parser.args)
