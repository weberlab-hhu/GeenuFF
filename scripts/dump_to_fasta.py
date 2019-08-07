#! /usr/bin/env python3
import argparse
import sys

from geenuff.applications.exporter import ExportController
from geenuff.base.orm import Coordinate, Feature, SuperLocus, Transcript, TranscriptPiece, \
    association_transcript_piece_to_feature
from geenuff.base.helpers import reverse_complement, chunk_str
from geenuff.base import types as gtypes


class FastaExportController(ExportController):
    def __init__(self, db_path_in):
        super().__init__(db_path_in)
        print(self.session, file=sys.stderr)
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
            out += '\n' + ''.join(subseq)
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
        if fa_out is None:
            handle_out = sys.stdout
        else:
            handle_out = open(fa_out, "w")
        handle_out.write(self.fmt_seq(self.export_seqs[0]))
        for export_seq in self.export_seqs[1:]:
            handle_out.write('\n')
            handle_out.write(self.fmt_seq(export_seq))
        handle_out.close()

    def transcript_length(self, transcript):
        """calculates total and coding length """
        pass  # todo

    def prep_intron_exports(self, genomes, exclude):
        # which genomes to export
        coord_ids = self._get_coords_by_genome_query(genomes, exclude)
        super_loci = self.get_super_loci_by_coords(coord_ids).all()

        for sl in super_loci:
            super_locus = self.session.query(SuperLocus).filter(SuperLocus.id == sl[0]).first()
            # todo, once JOIN output exists, drop all these loops
            print(super_locus, file=sys.stderr)
            for transcript in super_locus.transcripts:
                # todo, if transcript is longest_transcript
                for tps in transcript.transcript_pieces:
                    for feature in tps.features:
                        if feature.type.value == gtypes.GEENUFF_INTRON:
                            if feature.given_name is not None:
                                seqid = feature.given_name
                            else:
                                seqid = "intron_{0:06d}".format(feature.id)
                            self.export_seqs.append(
                                mk_one_fragment_seq(
                                    seqid,
                                    feature.coordinate_id,
                                    feature.start,
                                    feature.end,
                                    feature.is_plus_strand
                                )
                            )

    def get_super_loci_by_coords(self, coord_ids):
        sess = self.session  # shortcut
        # surprisingly not that bad, a minute or two for 55 genomes
        # todo, join instead of IN and keep results
        q = sess.query(Feature.id).filter(Feature.coordinate_id.in_(coord_ids))
        q = sess.query(association_transcript_piece_to_feature.c.transcript_piece_id)\
            .filter(association_transcript_piece_to_feature.c.feature_id.in_(q)).distinct()
        q = sess.query(TranscriptPiece.transcript_id)\
            .filter(TranscriptPiece.id.in_(q)).distinct()
        q = sess.query(Transcript.super_locus_id)\
            .filter(Transcript.id.in_(q)).distinct()
        return q


def mk_one_fragment_seq(seqid, coordinate_id, start, end, is_plus_strand):
    assert seqid is not None, "{} {} {} {} {}".format(seqid, coordinate_id, start, end, is_plus_strand)
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
        controller = FastaExportController(args.db_path_in)
        controller.prep_intron_exports(args.genomes, args.exclude_genomes)
        coords_ids = controller._get_coords_by_genome_query(args.genomes, args.exclude_genomes)
        x = controller.get_super_loci_by_coords(coords_ids)
        print(x, file=sys.stderr)
        print(x.all(), "hmmm", file=sys.stderr)
        controller.write_fa(args.out)


if __name__ == "__main__":
    implemented_modes = ["mRNA", "pre-mRNA", "CDS", "introns", "exons", "UTR"]
    parser = argparse.ArgumentParser()
    parser.add_argument('--db-path-in', type=str, required=True,
                        help='Path to the Geenuff SQLite input database.')
    parser.add_argument('-m', '--mode', type=str, required=True,
                        help="which type of sequence to export, one of {}".format(implemented_modes))
    parser.add_argument('-o', '--out', type=str,
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
