import sys

from geenuff.applications.exporter import GeenuffExportController
from geenuff.base.orm import Coordinate


class LengthExportController(GeenuffExportController):
    def __init__(self, db_path_in, longest=False):
        super().__init__(db_path_in, longest)
        self.export_ranges = []

    @staticmethod
    def get_length(export_group):
        total = 0
        for arange in export_group.ranges:
            total += abs(arange.start - arange.end)
        return total

    def write_lengths(self, file_out):
        if file_out is None:
            handle_out = sys.stdout
        else:
            handle_out = open(file_out, "w")
        for export_group in self.export_ranges:
            l = self.get_length(export_group)
            handle_out.write("{}\t{}\n".format(export_group.seqid, l))
        handle_out.close()
