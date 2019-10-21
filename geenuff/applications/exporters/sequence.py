import sys

from geenuff.applications.exporter import GeenuffExportController
from geenuff.base.orm import Coordinate
from geenuff.base.helpers import reverse_complement, chunk_str


class FastaExportController(GeenuffExportController):
    def __init__(self, db_path_in, longest=False):
        super().__init__(db_path_in, longest)

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
        for export_seq in self.export_ranges:
            handle_out.write(self.fmt_seq(export_seq))
            handle_out.write('\n')
        handle_out.close()



