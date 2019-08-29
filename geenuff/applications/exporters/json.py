import sys
from abc import ABC, abstractmethod
from geenuff.applications.exporter import ExportController
from geenuff.base.handlers import SuperLocusHandlerBase, TranscriptHandlerBase, TranscriptPieceHandlerBase, \
    FeatureHandlerBase

from geenuff.base.orm import Coordinate
from geenuff.base.helpers import reverse_complement, chunk_str


class ToJsonable(object):

    def __init__(self):
        self.jsonable_keys = ["id", "given_name"]

    @abstractmethod
    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        pass

    @abstractmethod
    def overlaps(self, coordinate, start, end, is_plus_strand):
        pass

    def pre_to_jsonable(self, data):
        out = {}
        for key in self.jsonable_keys:
            out[key] = data.__getattribute__(key)
        return out

    @abstractmethod
    def to_jsonable(self, data, coordinate, start, end, is_plus_strand):
        pass


# then multi inheritance of above and Handlers
class FeatureJsonable(FeatureHandlerBase, ToJsonable):
    def __init__(self, data=None):
        FeatureHandlerBase.__init__(self, data)
        ToJsonable.__init__(self)
        self.jsonable_keys += ["start", "start_is_biological_start", "end", "end_is_biological_end",
                               "is_plus_strand", "score", "source", "phase"]

    def to_jsonable(self, data, coordinate, start, end, is_plus_strand):
        out = self.pre_to_jsonable(self.data)
        out['type'] = data.type.value
        out['is_full_contained'] = self.is_fully_contained(coordinate, start, end, is_plus_strand)
        out['overlaps'] = self.overlaps(coordinate, start, end, is_plus_strand)
        return out

    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        if coordinate.id != self.data.coordinate_id or is_plus_strand != self.data.is_plus_strand:
            return False
        else:
            if is_plus_strand:
                if start <= self.data.start < end and start < self.data.end <= end:
                    return True
                else:
                    return False
            else:
                if start >= self.data.start > end and start > self.data.end >= end:
                    return True
                else:
                    return False

    def overlaps(self, coordinate, start, end, is_plus_strand):
        if coordinate.id != self.data.coordinate_id or is_plus_strand != self.data.is_plus_strand:
            return False
        else:
            if is_plus_strand:
                if start <= self.data.start < end or start < self.data.end <= end:
                    return True
                elif self.data.start <= start < self.data.end or self.data.start < end <= self.data.end:
                    return True
                else:
                    return False
            else:
                if start >= self.data.start > end or start > self.data.end >= end:
                    return True
                elif self.data.start >= start > self.data.end or self.data.start > end >= self.data.end:
                    return True
                else:
                    return False


class JsonExportController(ExportController):
    # todo, filter coordinate by start, end
    #  from orm obj (or join res) to json

    def coordinate_range_to_jsonable(self):
        pass




## target format reminder:
#[{"coordinate_piece":
#    {"id": int, "seqid": str, "sequence": str, "start": int, "end": int},
# "super_loci":
#    [{"id": str,
#      "given_name": str,
#      "is_fully_contained": bool,
#      "overlaps": bool,
#      "transcripts": [{"id": str,
#                       "given_name": str,
#                       "is_fully_contained": bool,
#                       "overlaps": bool,
#                       "features": [{"id": int,
#                                     "given_name": str,
#                                     "seqid": str,
#                                     "protein_id", str,
#                                     "type": enum,
#                                     "start": int,
#                                     "start_is_biological_start": bool,
#                                     "end": int,
#                                     "end_is_biological_end": bool,
#                                     "score": float,
#                                     "source": str,
#                                     "phase": int (in {0, 1, 2}),
#                                     "is_plus_strand": bool,
#                                     "is_fully_contained": bool,
#                                     "overlaps": bool
#                                    }, ...]
#                      }, ...]
#     }, ...]
#
#}, ...]