from dustdas import gffhelper, fastahelper
import copy
import hashlib
from types import GeneratorType

from . import orm
from . import types
from . import helpers


def convert2list(obj):
    if isinstance(obj, list):
        out = obj
    elif isinstance(obj, set) or isinstance(obj, GeneratorType) or isinstance(obj, tuple):
        out = list(obj)
    else:
        out = [obj]
    return out


class Handler(object):

    def __init__(self):
        self.data = None
        self.delete_me = False
        self.copyable = []
        self.linkable = []
        self.id = None

    def add_data(self, data):
        assert isinstance(data, self.data_type)
        if self.id is not None:
            data.id = self.id
        self.data = data
        data.handler = self  # terrible form, but I need some sort of efficient point back

    def copy_data_attr_to_other(self, other, copy_only=None, do_not_copy=None):
        if not isinstance(other, Handler):
            raise ValueError('other must be an instance of Handler, "{}" found'.format(type(other)))
        # everything that could be copied
        if copy_only is None:
            to_copy = list(type(self.data).__dict__.keys())
            to_copy = [x for x in to_copy if not x.startswith('_')]
            to_copy = set(copy.deepcopy(to_copy))
        else:
            copy_only = convert2list(copy_only)
            to_copy = set(copy_only)

        if do_not_copy is not None:
            do_not_copy = convert2list(do_not_copy)
            for item in do_not_copy:
                to_copy.remove(item)
        to_copy = copy.deepcopy(to_copy)
        for never_copy in ['id', 'handler']:
            try:
                to_copy.remove(never_copy)  # todo, confirm this is the primary key
            except KeyError:
                pass
        # actually copy
        for item in to_copy:
            val = self.get_data_attribute(item)
            other.set_data_attribute(item, val)

    def fax_all_attrs_to_another(self, another, skip_copying=None, skip_linking=None):
        linkable = copy.deepcopy(self.linkable)
        if skip_linking is not None:
            skip_linking = convert2list(skip_linking)
            for item in skip_linking:
                linkable.pop(linkable.index(item))
        copyable = copy.deepcopy(self.copyable)
        if skip_copying is not None:
            skip_copying = convert2list(skip_copying)
            for item in skip_copying:
                copyable.pop(copyable.index(item))
        self.copy_data_attr_to_other(another, copy_only=copyable)
        self.copy_selflinks_to_another(another, to_copy=linkable)

    def set_data_attribute(self, attr, val):
        self.data.__setattr__(attr, val)

    def get_data_attribute(self, attr):
        return self.data.__getattribute__(attr)

    def replace_selflinks_w_replacementlinks(self, replacement, to_replace):
        to_replace = copy.deepcopy(to_replace)
        for item in ['id', 'handler']:
            assert item not in to_replace
        for attr in to_replace:
            val = self.get_data_attribute(attr)
            if isinstance(val, list):
                n = len(val)
                for i in reversed(list(range(n))):  # go through backwards to hit every item even though we're removing
                    #for data in val:
                    data = val[i]
                    self.replace_selflink_with_replacementlink(replacement, data)
            elif isinstance(val, orm.Base):
                self.replace_selflink_with_replacementlink(replacement, val)
            else:
                raise ValueError("replace_selflinks_w_replacementlinks only implemented for {} types".format(
                    [list, orm.Base]
                ))

    def replace_selflink_with_replacementlink(self, replacement, data):
        other = data.handler
        self.de_link(other)
        replacement.link_to(other)

    def copy_selflinks_to_another(self, another, to_copy):
        to_copy = copy.deepcopy(to_copy)

        for item in ['id', 'handler']:
            assert item not in to_copy

        for attr in to_copy:
            val = self.get_data_attribute(attr)
            if isinstance(val, list):
                n = len(val)
                for i in reversed(list(range(n))):  # go through backwards to hit every item even though we're removing
                    #for data in val:
                    data = val[i]
                    self.copy_selflink_to_another(another, data)
            elif isinstance(val, orm.Base):
                self.copy_selflink_to_another(another, val)
            else:
                raise ValueError("copy_selflinks_to_another only implemented for {} types. {} found".format(
                    [list, orm.Base], type(val)
                ))

    # "selflink" for naming/usage consistency with 'replace' methods
    @staticmethod
    def copy_selflink_to_another(another, data):
        other = data.handler
        another.link_to(other)

    def link_to(self, other):
        raise NotImplementedError

    def de_link(self, other):
        raise NotImplementedError

    @property
    def data_type(self):
        raise NotImplementedError

    @property
    def _valid_links(self):
        raise NotImplementedError

    def _link_value_error(self, other):
        link_error = "from {} can only link / de_link to {}; found {}".format(type(self), self._valid_links,
                                                                              type(other))
        return ValueError(link_error)

    def mark_for_deletion(self):
        self.delete_me = True


class AnnotatedGenomeHandler(Handler):
    def __init__(self, controller=None):
        super().__init__()
        if controller is not None:
            self.id = controller.coord_handler_import_counter()
        self.mapper = None
        self._coords_by_seqid = None
        self._gffid_to_coords = None
        self._gff_seq_ids = None

    @staticmethod
    def hashseq(seq):
        m = hashlib.sha1()
        m.update(seq.encode())
        return m.hexdigest()

    @property
    def data_type(self):
        return orm.AnnotatedGenome

    @property
    def _valid_links(self):
        return [CoordinatesHandler]

    @property
    def gffid_to_coords(self):
        if not self._gffid_to_coords:
            self._gffid_to_coords = {}
            for gffid in self._gff_seq_ids:
                fa_id = self.mapper(gffid)
                x = self.coords_by_seqid[fa_id]
                self._gffid_to_coords[gffid] = x
        return self._gffid_to_coords

    @property
    def coords_by_seqid(self):
        if not self._coords_by_seqid:
            self._coords_by_seqid = {c.seqid: c for c in self.data.coordinates}
        return self._coords_by_seqid

    def link_to(self, other):
        if isinstance(other, CoordinatesHandler):
            other.data.annotated_genome = self.data
            # switched as below maybe checks other.data integrity, and fails on NULL anno genome?
            # self.data.coordinates.append(other.data)
        else:
            raise self._link_value_error(other)

    def de_link(self, other):
        if isinstance(other, CoordinatesHandler):
            self.data.coordinates.remove(other.data)
        else:
            raise self._link_value_error(other)

    def mk_mapper(self, gff_file=None):
        fa_ids = [e.seqid for e in self.data.coordinates]
        if gff_file is not None:  # allow setup without ado when we know IDs match exactly
            self._gff_seq_ids = helpers.get_seqids_from_gff(gff_file)
        else:
            self._gff_seq_ids = fa_ids
        mapper, is_forward = helpers.two_way_key_match(fa_ids, self._gff_seq_ids)
        self.mapper = mapper

        if not is_forward:
            raise NotImplementedError("Still need to implement backward match if fasta IDs "
                                      "are subset of gff IDs")

    def add_sequences(self, seq_file):
        self.add_fasta(seq_file)

    def add_fasta(self, seq_file, id_delim=' '):
        fp = fastahelper.FastaParser()
        for fasta_header, seq in fp.read_fasta(seq_file):
            seqid = fasta_header.split(id_delim)[0]
            # todo, parallelize sequence & annotation format, then import directly from ~Slice
            orm.Coordinates(start=0, end=len(seq), seqid=seqid, sha1=self.hashseq(seq),
                            annotated_genome=self.data)


class CoordinatesHandler(Handler):

    @property
    def data_type(self):
        return orm.Coordinates

    @property
    def _valid_links(self):
        return [AnnotatedGenomeHandler]

    def link_to(self, other):
        if isinstance(other, AnnotatedGenomeHandler):
            self.data.annotated_genome = other.data
        else:
            raise self._link_value_error(other)

    def de_link(self, other):
        if isinstance(other, AnnotatedGenomeHandler):
            other.data.coordinates.remove(self.data)
        else:
            raise self._link_value_error(other)

    def add_sequences(self, seq_file):
        raise NotImplementedError


class SuperLocusHandler(Handler):

    def __init__(self):
        super().__init__()
        self.handler_holder = HandleMaker(self)

    def make_all_handlers(self):
        self.handler_holder.make_all_handlers()

    @property
    def data_type(self):
        return orm.SuperLocus

    @property
    def _valid_links(self):
        return [TranscribedHandler, TranslatedHandler, FeatureHandler]

    def link_to(self, other):
        if any([isinstance(other, x) for x in [TranscribedHandler, TranslatedHandler, FeatureHandler]]):
            other.data.super_locus = self.data
        else:
            raise self._link_value_error(other)

    def de_link(self, other):
        if any([isinstance(other, x) for x in [TranscribedHandler, TranslatedHandler, FeatureHandler]]):
            other.data.super_locus = None
        else:
            raise self._link_value_error(other)

    def delete_marked_underlings(self, sess):
        for data in self.data.features + self.data.transcribeds + self.data.translateds:
            try:
                if data.handler.delete_me:
                    sess.delete(data)
            except AttributeError as e:
                'Warning, swallowing: {}\ndata: {}'.format(e, str(data))
        sess.commit()


class TranscribedHandler(Handler):
    def __init__(self):
        super().__init__()
        self.copyable += ['given_id', 'type']
        self.linkable += ['super_locus', 'transcribed_pieces', 'translateds']

    @property
    def data_type(self):
        return orm.Transcribed

    @property
    def _valid_links(self):
        return [TranslatedHandler, SuperLocusHandler, TranscribedPieceHandler]

    def link_to(self, other):
        if isinstance(other, SuperLocusHandler):
            self.data.super_locus = other.data
        elif isinstance(other, TranscribedPieceHandler):
            other.data.transcribed = self.data
        elif isinstance(other, TranslatedHandler):
            other.data.transcribeds.append(self.data)
        else:
            raise self._link_value_error(other)

    def de_link(self, other):
        if any([isinstance(other, x) for x in self._valid_links]):
            other.data.transcribeds.remove(self.data)
        else:
            raise self._link_value_error(other)


class TranscribedPieceHandler(Handler):
    def __init__(self):
        super().__init__()
        self.copyable += ['given_id']
        self.linkable += ['super_locus', 'features', 'transcribed']

    @property
    def data_type(self):
        return orm.TranscribedPiece

    @property
    def _valid_links(self):
        return [TranscribedHandler, SuperLocusHandler, FeatureHandler]

    def link_to(self, other):
        if isinstance(other, SuperLocusHandler):
            self.data.super_locus = other.data
        elif any([isinstance(other, x) for x in [TranscribedHandler, FeatureHandler]]):
            other.data.transcribed_pieces.append(self.data)
        else:
            raise self._link_value_error(other)

    def de_link(self, other):
        if any([isinstance(other, x) for x in self._valid_links]):
            other.data.transcribed_pieces.remove(self.data)
        else:
            raise self._link_value_error(other)


class TranslatedHandler(Handler):
    def __init__(self):
        super().__init__()
        self.copyable += ['given_id']
        self.linkable += ['super_locus', 'features', 'transcribeds']

    @property
    def data_type(self):
        return orm.Translated

    @property
    def _valid_links(self):
        return [TranscribedHandler, SuperLocusHandler, FeatureHandler]

    def link_to(self, other):
        if isinstance(other, SuperLocusHandler):
            self.data.super_locus = other.data
        elif any([isinstance(other, x) for x in [TranscribedHandler, FeatureHandler]]):
            other.data.translateds.append(self.data)
        else:
            raise self._link_value_error(other)

    def de_link(self, other):
        if any([isinstance(other, x) for x in self._valid_links]):
            other.data.translateds.remove(self.data)
        else:
            raise self._link_value_error(other)


class FeatureHandler(Handler):
    def __init__(self):
        super().__init__()
        self.copyable += ['given_id', 'type', 'position', 'coordinates', 'is_plus_strand', 'score', 'source', 'phase']
        self.linkable += ['super_locus', 'transcribed_pieces', 'translateds']

    @property
    def data_type(self):
        return orm.Feature

    @property
    def _valid_links(self):
        return [TranscribedPieceHandler, SuperLocusHandler, TranslatedHandler]

    def link_to(self, other):
        if isinstance(other, SuperLocusHandler):
            self.data.super_locus = other.data
        elif any([isinstance(other, x) for x in [TranslatedHandler, TranscribedPieceHandler]]):
            other.data.features.append(self.data)
        else:
            raise self._link_value_error(other)

    def de_link(self, other):
        if any([isinstance(other, x) for x in self._valid_links]):
            other.data.features.remove(self.data)
        else:
            raise self._link_value_error(other)

    def cmp_key(self):
        return self.data.cmp_key()

    def pos_cmp_key(self):
        return self.data.pos_cmp_key()


#### section TranscriptInterpreter, might end up in a separate file later
class TranscriptStatus(object):
    """can hold and manipulate all the info on current status of a transcript"""

    def __init__(self):
        # initializes to intergenic
        self.genic = False
        self.in_intron = False
        self.in_trans_intron = False
        self.in_translated_region = False
        self.seen_start = False  # todo, move EUK specific stuff to subclass?
        self.seen_stop = False
        self.erroneous = False
        self.phase = None  # todo, proper tracking / handling

    def __repr__(self):
        return "genic: {}, intronic: {}, translated_region: {}, trans_intronic: {}, phase: {}".format(
            self.genic, self.in_intron, self.in_translated_region, self.in_trans_intron, self.phase
        )

    @property
    def _decoder(self):
        # todo, parallelize status until this isn't necessary
        return {
            types.TRANSCRIBED: ('genic', self.saw_tss, self.saw_tts),
            types.CODING: ('in_translated_region', self.saw_start, self.exit_coding),
            types.INTRON: ('in_intron', self.splice_open, self.splice_close),
            types.TRANS_INTRON: ('in_trans_intron', self.trans_splice_open, self.trans_splice_close),
            types.ERROR: ('erroneous', self.error_open, self.error_close)
        }

    def update_for_feature(self, feature, **kwargs):
        attr, fn_open, fn_close = self._decoder[feature.type.value]
        if feature.bearing.value in [types.START, types.OPEN_STATUS]:
            fn_open(**kwargs)
        elif feature.bearing.value in [types.END, types.CLOSE_STATUS]:
            fn_close(**kwargs)
        else:
            raise ValueError('unhandled bearing {}'.format(feature.bearing))

    def saw_tss(self):
        self.genic = True

    def saw_start(self, phase):
        self.seen_start = True  # todo, disentangle from annotations core -> Euk specific/parser only
        self.in_translated_region = True
        self.phase = phase

    def exit_coding(self):
        self.in_translated_region = False

    def saw_stop(self):
        self.seen_stop = True  # todo, disentangle
        self.in_translated_region = False
        self.phase = None

    def saw_tts(self):
        self.genic = False

    def splice_open(self):
        self.in_intron = True

    def splice_close(self):
        self.in_intron = False

    def trans_splice_open(self):
        self.in_trans_intron = True

    def trans_splice_close(self):
        self.in_trans_intron = False

    def error_open(self):
        self.erroneous = True

    def error_close(self):
        self.erroneous = False

    def is_5p_utr(self):
        return self.is_utr() and not any([self.seen_start, self.seen_stop])

    def is_3p_utr(self):
        return self.is_utr() and self.seen_stop and self.seen_start

    def is_utr(self):
        return self.genic and not any([self.in_intron, self.in_translated_region, self.in_trans_intron])

    def is_coding(self):
        return self.genic and self.in_translated_region and not any([self.in_intron, self.in_trans_intron])

    def is_intronic(self):
        return self.in_intron and self.genic

    def is_trans_intronic(self):
        return self.in_trans_intron and self.genic

    def is_intergenic(self):
        return not self.genic


class TransitionStep(object):
    def __init__(self, features=None, status=None, piece=None):
        self.features = features
        self.status = status
        self.piece = piece
        self.previous_range = None

    def make_range(self, previous_step):
        # todo, and this is where I realize exclusive closing elements really really are needed...
        pass

    @property
    def a_feature(self):
        if self.features is None:
            return None
        else:
            return self.features[0]


class Range(object):
    def __init__(self, seqid, start, end, status):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.status = status


def positional_match(feature, previous):
    return feature.pos_cmp_key() == previous.pos_cmp_key()


def bearing_match(feature, previous):
    return feature.bearing.value == previous.bearing.value


class TranscriptInterpBase(object):
    # todo, move this to generic location and/or skip entirely
    def __init__(self, transcript, super_locus, session=None):
        assert isinstance(transcript, TranscribedHandler)
        self.status = TranscriptStatus()
        self.transcript = transcript
        self.session = session
        self.super_locus = super_locus

    def transition_5p_to_3p(self):
        status = TranscriptStatus()
        for piece in self.sort_pieces():
            piece_features = self.sorted_features(piece)
            for aligned_features in self.full_stack_matches(piece_features):
                self.update_status(status, aligned_features)
                yield aligned_features, copy.deepcopy(status), piece

    def transition_with_ranges(self):
        """organize [prev. range]-> feature pairs along transcript"""
        pass  # todo

    @staticmethod
    def sorted_features(piece):
        features = piece.features
        # confirm strand & seqid
        assert all([f.coordinates == features[0].coordinates for f in features]), \
            'not all matching: {}'.format([(f.id, f.coordinates) for f in features])
        assert all([f.is_plus_strand == features[0].is_plus_strand for f in features]), \
            'not all matching: {}'.format([(f.id, f.is_plus_strand) for f in features])
        features = sorted(features, key=lambda x: x.pos_cmp_key())
        if not features[0].is_plus_strand:
            features.reverse()
        return features

    def sort_pieces(self):
        pieces = self.transcript.data.transcribed_pieces
        ordered_pieces = sorted(pieces, key=lambda x: x.position)
        return ordered_pieces

    @staticmethod
    def update_status(status, aligned_features):
        for feature in aligned_features:
            ftype = feature.type.value
            fbearing = feature.bearing.value
            # standard features
            if ftype == types.CODING and fbearing == types.START:
                status.update_for_feature(feature, phase=0)
            elif ftype == types.CODING and fbearing == types.OPEN_STATUS:
                status.update_for_feature(feature, phase=feature.phase)
            elif ftype == types.CODING and fbearing == types.END:
                status.update_for_feature(feature)
                status.saw_stop()  # todo, disentangle / to-parser not general section
            else:
                status.update_for_feature(feature)

    @staticmethod
    def stack_matches(features, match_fn=positional_match):
        ifeatures = iter(features)
        try:
            prev = next(ifeatures)
        except StopIteration:
            return
        current = [prev]
        for feature in ifeatures:
            if match_fn(feature, prev):
                current.append(feature)
            else:
                yield current
                current = [feature]
            prev = feature
        yield current
        return

    @staticmethod
    def sort_by_bearing(matches):
        key = {types.START: 3,
               types.OPEN_STATUS: 2,
               types.CLOSE_STATUS: 1,
               types.END: 0,
               types.POINT: 4}

        return sorted(matches, key=lambda x: key[x.bearing.value])

    def full_stack_matches(self, features):
        for matches in self.stack_matches(features, match_fn=positional_match):
            sorted_matches = self.sort_by_bearing(matches)
            for by_bearing in self.stack_matches(sorted_matches, match_fn=bearing_match):
                yield by_bearing

    def sort_all(self):
        out = []
        for piece in self.sort_pieces():
            out.append(self.sorted_features(piece))
        return out


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
        datas = list()
        datas += sl.transcribed_pieces
        datas += sl.translateds
        datas += sl.features
        datas += sl.transcribeds

        for item in datas:
            self.handles.append(self._get_or_make_one_handler(item))

    def mk_n_append_handler(self, data):
        handler = self._get_or_make_one_handler(data)
        self.handles.append(handler)
        return handler

    def _get_or_make_one_handler(self, data):
        try:
            handler = data.hanlder
        except AttributeError:
            handler_type = self._get_handler_type(data)
            handler = handler_type()
            handler.add_data(data)
        return handler

    def _get_handler_type(self, old_data):
        key = [(SuperLocusHandler, orm.SuperLocus),
               (TranscribedHandler, orm.Transcribed),
               (TranslatedHandler, orm.Translated),
               (TranscribedPieceHandler, orm.TranscribedPiece),
               (FeatureHandler, orm.Feature)]

        return self._get_paired_item(type(old_data), search_col=1, return_col=0, nested_list=key)
