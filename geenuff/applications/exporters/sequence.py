import sys

from geenuff.applications.exporter import ExportController
from geenuff.base.orm import Coordinate, SuperLocus
from geenuff.base.helpers import reverse_complement, chunk_str
from geenuff.applications.exporter import SuperLocusRanger

class FastaExportController(ExportController):
    def __init__(self, db_path_in, longest=False):
        super().__init__(db_path_in)
        print(self.session, file=sys.stderr)
        self.export_seqs = []
        self.longest = longest

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

    def prep_intron_exports(self, genomes, exclude):
        # which genomes to export
        coord_ids = self._get_coords_by_genome_query(genomes, exclude)
        super_loci = self.get_super_loci_by_coords(coord_ids).all()
        i = 0
        for sl in super_loci:
            super_locus = self.session.query(SuperLocus).filter(SuperLocus.id == sl[0]).first()
            sl_ranger = SuperLocusRanger(super_locus, longest=self.longest)
            # todo, once JOIN output exists, drop all these loops
            print(super_locus, file=sys.stderr)
            for range_maker in sl_ranger.exp_range_makers:
                # todo, if transcript is longest_transcript
                introns = range_maker.intronic_ranges()

                for intron in introns:
                    seqid = "intron_{0:06d}".format(i)
                    self.export_seqs.append(
                        mk_one_fragment_seq(
                            seqid,
                            intron.coordinate_id,
                            intron.start,
                            intron.end,
                            intron.is_plus_strand
                        )
                    )
                    i += 1




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

