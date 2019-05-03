import copy
from types import GeneratorType

from . import orm
from . import types


class Handler(object):

    def __init__(self, data=None):
        self.id = None
        if data is not None:
            self.add_data(data)
        else:
            self.data = None

    def add_data(self, data):
        assert isinstance(data, self.data_type), 'data type does not match with handler class'
        if self.id is not None:
            data.id = self.id
        self.data = data
        data.handler = self  # terrible form, but I need some sort of efficient point back

    @property
    def data_type(self):
        raise NotImplementedError


class GenomeHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Genome


class CoordinateHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Coordinate


class SuperLocusHandlerBase(Handler):

    def __init__(self):
        super().__init__()
        self.handler_holder = HandleMaker(self)

    def make_all_handlers(self):
        self.handler_holder.make_all_handlers()

    @property
    def data_type(self):
        return orm.SuperLocus


class TranscribedHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Transcribed


class TranscribedPieceHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.TranscribedPiece


class TranslatedHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Translated


class FeatureHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Feature

    def cmp_key(self):
        return self.data.cmp_key()

    def pos_cmp_key(self):
        return self.data.pos_cmp_key()


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
        datas = sl.translateds + sl.transcribeds

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
               (TranscribedHandlerBase, orm.Transcribed),
               (TranslatedHandlerBase, orm.Translated)]

        return self._get_paired_item(type(old_data), search_col=1, return_col=0, nested_list=key)
