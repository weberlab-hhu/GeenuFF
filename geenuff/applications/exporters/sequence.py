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
        last_coordinate, sequence = None, None
        for a_range in export_sequence.ranges:
            coord_id = a_range.coordinate_id
            if last_coordinate != coord_id:
                coordinate = self.session.query(Coordinate).filter(Coordinate.id == coord_id).first()
                sequence = coordinate.sequence
            out += self.get_seq_fragment(a_range, sequence)
            last_coordinate = coord_id
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
        for export_seq in self.export_seqs:
            handle_out.write(self.fmt_seq(export_seq))
            handle_out.write('\n')
        handle_out.close()

    def prep_intron_exports(self, genomes, exclude):

        def range_function(range_maker):
            return range_maker.intronic_ranges()

        def id_function(i, range):
            del range
            return "intron_{0:06d}".format(i)

        self._prep_geenuff_features(genomes, exclude, range_function, id_function)

    def prep_transcript_exports(self, genomes, exclude):

        def range_function(range_maker):
            return range_maker.transcribed_ranges()

        def id_function(i, range):
            del i
            return range.given_name

        self._prep_geenuff_features(genomes, exclude, range_function, id_function)

    def _prep_geenuff_features(self, genomes, exclude, range_function, id_function):
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
                features = range_function(range_maker)
                for feature in features:
                    seqid = id_function(i, feature)
                    self.export_seqs.append(
                        ExportSeq(seqid=seqid,
                                  ranges=[feature])
                    )
                    i += 1


class ExportSeq(object):
    def __init__(self, seqid, ranges):
        self.seqid = seqid
        self.ranges = ranges
