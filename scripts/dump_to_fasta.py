#! /usr/bin/env python3
import argparse
import sys

from geenuff.applications.exporter import ExportController
from geenuff.base.orm import Coordinate, Feature
from geenuff.base.helpers import reverse_complement, chunk_str
from geenuff.base import types as gtypes


class FastaExportController(ExportController):
    def __init__(self, db_path_in, with_features_only=True):
        super().__init__(db_path_in, with_features_only)
        print(self.session)
        self.export_seqs = []

    def get_seq(self, export_sequence):
        # simple and probably horribly inefficient
        out = []
        for seq_piece in export_sequence.pieces:
            coordinate = self.session.query(Coordinate).filter(Coordinate.id == seq_piece.coord_id).first()
            sequence = coordinate.sequence
            for fragment in seq_piece.fragments:
                out += self.get_seq_fragment(fragment, sequence)
        return out

    def fmt_seq(self, export_sequence):
        out = '>' + export_sequence.seqid
        seq = self.get_seq(export_sequence)
        for subseq in chunk_str(seq, 80):
            out += '\n' + subseq
        return out

    @staticmethod
    def get_seq_fragment(fragment, sequence):
        if fragment.is_plus_strand:
            out = sequence[fragment.start:fragment.end]
        else:
            # +1 to flip inclusive/exclusive
            out = sequence[(fragment.end + 1):(fragment.start + 1)]
            out = reverse_complement(out)
        return out

    def write_fa(self, fa_out):
        # again, simple and probably horribly inefficient
        with open(fa_out, 'w') as f:
            f.write(self.fmt_seq(self.export_seqs[0]))
            for export_seq in self.export_seqs[1:]:
                f.write('\n')
                f.write(self.fmt_seq(export_seq))

    def prep_intron_exports(self):
        # which genomes to export
        coord_ids = self._coords_with_feature_query()
        # introns from those genomes
        introns = self.session.query(Feature)\
            .filter(Feature.type == gtypes.GEENUFF_INTRON)\
            .filter(Feature.coordinate_id.in_(coord_ids))
        for intron in introns:
            self.export_seqs.append(
                mk_one_fragment_seq(
                    intron.given_name,
                    intron.coordinate_id,
                    intron.start,
                    intron.end,
                    intron.is_plus_strand
                )
            )


def mk_one_fragment_seq(seqid, coordinate_id, start, end, is_plus_strand):
    export_seq = ExportSeq(seqid=seqid)
    fragment = ExportSeqFragment(start=start,
                                 end=end,
                                 is_plus_strand=is_plus_strand)
    seq_piece = ExportSeqPiece(coord_id=coordinate_id)
    seq_piece.fragments = [fragment]
    export_seq.pieces = [seq_piece]
    return export_seq


class ExportSeq(object):
    def __init__(self, seqid):
        self.seqid = seqid
        self.pieces = []


class ExportSeqPiece(object):
    def __init__(self, coord_id):
        self.coord_id = coord_id
        self.fragments = []


class ExportSeqFragment(object):
    def __init__(self, start, end, is_plus_strand):
        self.start = start
        self.end = end
        self.is_plus_strand = is_plus_strand


def main(args):
    if args.mode != "introns":
        raise NotImplementedError("I lied, only mode=introns is implemented so far, not {}".format(args.mode))
    else:
        controller = FastaExportController(args.db_path_in,
                                           with_features_only=True)
        controller.prep_intron_exports()
        controller.write_fa(args.out)


if __name__ == "__main__":
    implemented_modes = ["mRNA", "pre-mRNA", "CDS", "introns", "exons", "UTR"]
    parser = argparse.ArgumentParser()
    parser.add_argument('--db-path-in', type=str, required=True,
                        help='Path to the Geenuff SQLite input database.')
    parser.add_argument('-m', '--mode', type=str, required=True,
                        help="which type of sequence to export, one of {}".format(implemented_modes))
    parser.add_argument('-o', '--out', type=str, default=sys.stdout,
                        help="output fasta file path, (default is stdout)")
    parser.add_argument('--genomes', type=str, default='',
                        help='Comma separated list of species names to be exported. '
                             'If empty all genomes in the db are used.')
    parser.add_argument('--exclude-genomes', type=str, default='',
                        help='Comma separated list of species names to be excluded. '
                             'If empty all genomes in the db are used.')

    parser.add_argument('-l', '--longest', action="store_true",
                        help="ignore all but the longest transcript per gene")
    args = parser.parse_args()

    main(args)
