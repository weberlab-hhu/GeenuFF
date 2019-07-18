import copy

from . import orm
from . import types
from . import helpers


class Handler(object):

    def __init__(self, data=None):
        self.id = None
        if data is not None:
            self.add_data(data)
        else:
            self.data = None

    def add_data(self, data):
        assert isinstance(data, self.data_type), 'data type does not match with handler class'
        self.data = data

    @property
    def data_type(self):
        raise NotImplementedError

    def __repr__(self):
        return helpers.get_repr('Handler', {'data': self.data})


class GenomeHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Genome


class CoordinateHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Coordinate


class SuperLocusHandlerBase(Handler):

    def __init__(self, data=None):
        super().__init__(data)
        self.handler_holder = HandleMaker(self)

    @property
    def data_type(self):
        return orm.SuperLocus

    @property
    def features(self):
        for transcript in self.data.transcripts:
            for piece in transcript.transcript_pieces:
                for feature in piece.features:
                    yield feature

    def make_all_handlers(self):
        self.handler_holder.make_all_handlers()


class TranscriptHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Transcript

    @property
    def sorted_pieces(self):
        pieces = self.data.transcript_pieces
        return sorted(pieces, key=lambda p: p.position)


class TranscriptPieceHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.TranscriptPiece


class ProteinHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Protein


class FeatureHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Feature


class HandleMaker(object):
    def __init__(self, super_locus_handler):
        self.super_locus_handler = super_locus_handler
        self.handles = []

    @staticmethod
    def _get_paired_item(search4, search_col, return_col, nested_list):
        matches = [x[return_col] for x in nested_list if x[search_col] == search4]
        assert len(matches) == 1
        return matches[0]

    def make_all_handlers(self):
        self.handles = []
        sl = self.super_locus_handler.data
        datas = sl.proteins + sl.transcripts

        for item in datas:
            self.handles.append(self._get_or_make_one_handler(item))

    def mk_n_append_handler(self, data):
        handler = self._get_or_make_one_handler(data)
        self.handles.append(handler)
        return handler

    def _get_or_make_one_handler(self, data):
        try:
            handler = data.handler
        except AttributeError:
            handler_type = self._get_handler_type(data)
            handler = handler_type()
            handler.add_data(data)
        return handler

    def _get_handler_type(self, old_data):
        key = [(SuperLocusHandlerBase, orm.SuperLocus),
               (TranscriptHandlerBase, orm.Transcript),
               (ProteinHandlerBase, orm.Protein)]

        return self._get_paired_item(type(old_data), search_col=1, return_col=0, nested_list=key)
