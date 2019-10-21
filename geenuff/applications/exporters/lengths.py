import sys
import numpy

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

    @staticmethod
    def _as_file_handle(file_out):
        if file_out is None:
            handle_out = sys.stdout
        else:
            handle_out = open(file_out, "w")
        return handle_out

    def write_lengths(self, file_out):
        handle_out = self._as_file_handle(file_out)
        for export_group in self.export_ranges:
            l = self.get_length(export_group)
            handle_out.write("{}\t{}\n".format(export_group.seqid, l))
        handle_out.close()

    def write_length_stats(self, file_out):
        lengths = []
        for export_group in self.export_ranges:
            lengths.append(self.get_length(export_group))

        stats = basics(lengths)
        quants = fmt_keys(quantiles(lengths), pfx="quantile")
        nxes = fmt_keys(nx(lengths), pfx="N", sfx="")

        handle_out = self._as_file_handle(file_out)
        for group in [stats, quants, nxes]:
            handle_out.write(fmt_stats(group))
        handle_out.close()


def fmt_stats(a_dict):
    out = ""
    for key in sorted(a_dict.keys()):
        out += "{}\t{}\n".format(key, a_dict[key])
    return out


def basics(lengths):
    out = {'count': len(lengths),
           'longest': max(lengths),
           'shortest': min(lengths),
           'total': sum(lengths)}
    return out


def nx(lengths, x_vals=None):
    if x_vals is None:
        x_vals = [.10, .25, .50, .75, .90]

    out = {}
    for x in x_vals:
        out[x] = None

    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    cummulative = 0
    for l in lengths:
        cummulative += l
        # set as soon as we're beyond target value (and not again)
        for key in out:
            if out[key] is None:  # skip set values
                if cummulative >= key * total:
                    out[key] = l
    return out


def quantiles(lengths, x_vals=None):
    if x_vals is None:
        x_vals = [.10, .25, .50, .75, .90]

    out = {}
    for key in x_vals:
        out[key] = numpy.quantile(lengths, key)

    return out


def fmt_keys(a_dict, pfx, sfx="%", times_by=100):
    """convert fraction keys to labelled percentages (or similar)"""
    out = {}
    for key in a_dict:
        new_key = "{}{}{}".format(pfx, int(key * times_by), sfx)
        out[new_key] = a_dict[key]
    return out
