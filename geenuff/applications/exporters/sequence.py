import sys

from geenuff.applications.exporter import GeenuffExportController
from geenuff.base.orm import Coordinate
from geenuff.base.helpers import reverse_complement, chunk_str


class FastaExportController(GeenuffExportController):
    def __init__(self, db_path_in, longest=False):
        super().__init__(db_path_in, longest)
        self.coordinate_id_cache = None
        self.sequence_cache = None

    def get_seq(self, export_sequence):
        # simple and probably horribly inefficient
        out = []
        for a_range in export_sequence.ranges:
            coord_id = a_range.coordinate_id
            if self.coordinate_id_cache != coord_id:
                coordinate = self.session.query(Coordinate).filter(Coordinate.id == coord_id).first()
                self.sequence_cache = coordinate.sequence
            out += self.get_seq_fragment(a_range, self.sequence_cache)
            self.coordinate_id_cache = coord_id
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
        handle_out = self._as_file_handle(fa_out)
        for export_seq in self.export_ranges:
            handle_out.write(self.fmt_seq(export_seq))
            handle_out.write('\n')
        handle_out.close()



