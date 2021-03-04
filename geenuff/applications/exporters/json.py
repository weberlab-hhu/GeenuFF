import sys
import json
from abc import ABC, abstractmethod
from geenuff.applications.exporter import GeenuffExportController
from geenuff.base.handlers import SuperLocusHandlerBase, TranscriptHandlerBase, CoordinateHandlerBase, \
    FeatureHandlerBase

from geenuff.base.orm import Coordinate, Transcript, SuperLocus


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

    def to_jsonable(self, data, coordinate, start, end, is_plus_strand, transcript=None):
        assert transcript is not None, "needed to get the protein id"
        out = self.pre_to_jsonable(self.data)
        out['type'] = data.type.value
        out['is_fully_contained'] = self.is_fully_contained(coordinate, start, end, is_plus_strand)
        out['overlaps'] = self.overlaps(coordinate, start, end, is_plus_strand)
        out['protein_id'] = self.protein_ids(transcript)
        return out

    def protein_ids(self,  transcript):
        assert isinstance(transcript, Transcript)  # just for pycharm hints for now...
        p_ids = list(set([x.given_name for x in transcript.proteins]).intersection(
            set([x.given_name for x in self.data.proteins])
        ))
        l_proteins = len(p_ids)
        if l_proteins == 0:
            return None
        elif l_proteins == 1:
            return p_ids[0]
        else:
            # todo, actually, I don't think this should be able to happen, but double check
            raise NotImplementedError("what to do with {} proteins?".format(len(p_ids)))

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


class TranscriptJsonable(TranscriptHandlerBase, ToJsonable):
    def __init__(self, data=None):
        FeatureHandlerBase.__init__(self, data)
        ToJsonable.__init__(self)
        self.feature_handlers = self._mk_feature_handlers()

    def _mk_feature_handlers(self):
        return [FeatureJsonable(x) for x in self.sorted_features()]

    def sorted_features(self):
        out = []
        for piece in self.sorted_pieces:
            # todo, access the sort_cmp_keys for consistency with code base?
            if piece.features[0].is_plus_strand:
                features = sorted(piece.features, key=lambda f: (f.start, f.end))
            else:
                features = sorted(piece.features, key=lambda f: (f.start * -1, f.end * 1))
            out += features
        return out

    # todo, use pre-calculated jsons not fresh method call
    def overlaps(self, coordinate, start, end, is_plus_strand):
        if any([fh.overlaps(coordinate, start, end, is_plus_strand) for fh in self.feature_handlers]):
            return True
        else:
            return False

    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        if all([fh.is_fully_contained(coordinate, start, end, is_plus_strand) for fh in self.feature_handlers]):
            return True
        else:
            return False

    def to_jsonable(self, data, coordinate, start, end, is_plus_strand):
        out = self.pre_to_jsonable(self.data)
        out["type"] = self.data.type.value
        out["is_fully_contained"] = self.is_fully_contained(coordinate, start, end, is_plus_strand)
        out["overlaps"] = self.overlaps(coordinate, start, end, is_plus_strand)
        out["features"] = [fh.to_jsonable(fh.data, coordinate, start, end, is_plus_strand, self.data)
                           for fh in self.feature_handlers]
        return out


class SuperLocusJsonable(SuperLocusHandlerBase, ToJsonable):
    def __init__(self, data=None):
        FeatureHandlerBase.__init__(self, data)
        ToJsonable.__init__(self)
        self.transcript_handlers = self._mk_transcript_handlers()

    def _mk_transcript_handlers(self):
        return [TranscriptJsonable(x) for x in self.data.transcripts]

    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        if all([th.is_fully_contained(coordinate, start, end, is_plus_strand) for th in self.transcript_handlers]):
            return True
        else:
            return False

    def overlaps(self, coordinate, start, end, is_plus_strand):
        if any([th.overlaps(coordinate, start, end, is_plus_strand) for th in self.transcript_handlers]):
            return True
        else:
            return False

    def to_jsonable(self, data, coordinate, start, end, is_plus_strand):
        out = self.pre_to_jsonable(self.data)
        out['type'] = self.data.type.value
        out['is_fully_contained'] = self.is_fully_contained(coordinate, start, end, is_plus_strand)
        out['overlaps'] = self.overlaps(coordinate, start, end, is_plus_strand)
        out['transcripts'] = [th.to_jsonable(th, coordinate, start, end, is_plus_strand)
                              for th in self.transcript_handlers]
        return out


class CoordinateJsonable(CoordinateHandlerBase):
    def to_jsonable(self, start, end):
        return {'id': self.data.id,
                'seqid': self.data.seqid,
                'sequence': self.data.sequence[start:end],
                'start': start,
                'end': end}


class JsonExportController(GeenuffExportController):
    # todo, filter coordinate by start, end
    #  from orm obj (or join res) to json

    def coordinate_range_to_jsonable(self, species, seqid, start, end, is_plus_strand):
        out = []
        for coordinate in self.session.query(Coordinate).filter(Coordinate.seqid == seqid).all():
            if coordinate.genome.species == species:
                if end is None:  # if end is not specified, take whole sequence
                    end = coordinate.length
                ch = CoordinateJsonable(coordinate)
                res = {'coordinate_piece': ch.to_jsonable(start, end),
                       'super_loci': []}
                for sl, sl_coordinate_seqid in self.genome_query(return_super_loci=True):
                    if sl_coordinate_seqid == seqid:
                        slh = SuperLocusJsonable(sl)
                        if slh.overlaps(coordinate, start, end, is_plus_strand):
                            res['super_loci'].append(slh.to_jsonable(slh.data, coordinate, start, end,
                                                                     is_plus_strand))
                out.append(res)
        return out

    def coordinate_range_to_json(self, species, seqid, start, end, is_plus_strand):
        return json.dumps(self.coordinate_range_to_jsonable(species, seqid, start, end, is_plus_strand))

    def query_and_write(self, species, seqid, start, end, is_plus_strand, file_out, pretty=False):
        jsonable = self.coordinate_range_to_jsonable(species, seqid, start, end, is_plus_strand)
        if pretty:
            dumps = json.dumps(jsonable, indent=2)
        else:
            dumps = json.dumps(jsonable)

        handle_out = self._as_file_handle(file_out)
        handle_out.write(dumps)
        handle_out.close()




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
