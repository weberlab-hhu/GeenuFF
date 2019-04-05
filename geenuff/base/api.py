from dustdas import gffhelper, fastahelper
import copy
import hashlib
from types import GeneratorType
import intervaltree

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

    @property
    def data_type(self):
        raise NotImplementedError


class AnnotatedGenomeHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.AnnotatedGenome


class CoordinatesHandlerBase(Handler):

    @property
    def data_type(self):
        return orm.Coordinates


class SuperLocusHandlerBase(Handler):

    def __init__(self):
        super().__init__()
        self.handler_holder = HandleMaker(self)

    def make_all_handlers(self):
        self.handler_holder.make_all_handlers()

    @property
    def data_type(self):
        return orm.SuperLocus

    # todo maybe to this with ondelete statements in the db
    def delete_marked_underlings(self, sess):
        for data in self.data.features + self.data.transcribeds + self.data.translateds:
            try:
                if data.handler.delete_me:
                    sess.delete(data)
            except AttributeError as e:
                'Warning, swallowing: {}\ndata: {}'.format(e, str(data))
        sess.commit()


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


#### section TranscriptInterpreter, todo maybe put into a separate file later
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


class TranscriptCoordinate(object):

    def __init__(self, coordinate_id, piece_position, is_plus_strand, start):
        self.start = start
        self.coordinate_id = coordinate_id
        self.piece_position = piece_position
        self.is_plus_strand = is_plus_strand

    def sort_key(self):
        if self.is_plus_strand:
            sort_pos = self.start
        else:
            sort_pos = -self.start  # flip sort order on the - strand
        return self.piece_position, sort_pos

    def sequence_chunk_info(self):
        return self.coordinate_id, self.is_plus_strand, self.piece_position

    def __repr__(self):
        return "coordinate: {}, piece position {}, is_plus {}: {}".format(self.coordinate_id, self.piece_position,
                                                                          self.is_plus_strand, self.start)

    def __eq__(self, other):
        if isinstance(other, TranscriptCoordinate):
            return self.__dict__ == other.__dict__
        return False


class Range(TranscriptCoordinate):
    def __init__(self, coordinate_id, piece_position, start, end, is_plus_strand):
        super().__init__(coordinate_id=coordinate_id, piece_position=piece_position, is_plus_strand=is_plus_strand,
                         start=start)
        self.end = end

    def __repr__(self):
        return "coordinate: {}, piece position {}, is_plus {}: {}-{}".format(self.coordinate_id, self.piece_position,
                                                                             self.is_plus_strand, self.start, self.end)


def positional_match(feature, previous):
    return feature.pos_cmp_key() == previous.pos_cmp_key()


def bearing_match(feature, previous):
    return feature.bearing.value == previous.bearing.value


class TranscriptInterpBase(object):
    # todo, move this to generic location and/or skip entirely
    def __init__(self, transcript, super_locus, session=None):
        assert isinstance(transcript, TranscribedHandlerBase)
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

    # helpers for classic transitions below
    def _ranges_by_type(self, target_type):
        ranges = []
        current = None
        for aligned_features, status, piece in self.transition_5p_to_3p():
            features_of_type = [f for f in aligned_features if f.type.value == target_type]
            assert len(features_of_type) in [0, 1], "cannot interpret aligned features of the same type {}".format(
                features_of_type
            )
            if len(features_of_type) == 1:  # and 0 is simply ignored...
                feature = features_of_type[0]
                if self._is_open_or_start(feature):
                    current = Range(coordinate_id=feature.coordinate_id,
                                    start=feature.position,
                                    end=None,
                                    is_plus_strand=feature.is_plus_strand,
                                    piece_position=piece.position)
                else:
                    assert current is not None, "start/open must be seen before end/close for {}".format(feature)
                    assert current.piece_position == piece.position, \
                        "can't have start/end from different pieces ({} & {})".format(current.piece_position,
                                                                                      piece.position)
                    current.end = feature.position
                    ranges.append(current)
                    current = None
        return ranges

    @staticmethod
    def _is_open_or_start(feature):
        if feature.bearing.value in [types.START, types.OPEN_STATUS]:
            return True
        elif feature.bearing.value in [types.END, types.CLOSE_STATUS]:
            return False
        else:
            raise ValueError("failed to interpret bearing {}".format(feature.bearing))

    @staticmethod
    def _make_trees(ranges):
        trees = {}
        for r in ranges:
            coord_isplus = r.sequence_chunk_info()
            if coord_isplus not in trees:
                trees[coord_isplus] = intervaltree.IntervalTree()
            start, end = min(r.start, r.end), max(r.start, r.end)
            trees[coord_isplus][start:end] = r
        return trees

    @staticmethod
    def _copy_ival_data(iv, islower):
        print(iv, islower, "was called")
        if islower:  # copy one of the two sides so that we don't change the same dictionary later
            out = copy.deepcopy(iv.data)
        else:
            out = iv.data
        return out

    def _subtract_ranges(self, subtract_from, to_subtract):
        keep_trees = self._make_trees(subtract_from)
        subtract_trees = self._make_trees(to_subtract)
        subtracted = []
        for key in keep_trees:
            for chop_out in subtract_trees[key]:
                keep_trees[key].chop(chop_out.begin, chop_out.end, self._copy_ival_data)
            for kept in keep_trees[key]:
                start, end = kept.begin, kept.end
                if not kept.data.is_plus_strand:
                    start, end = end, start
                kept.data.start = start
                kept.data.end = end
                subtracted.append(kept.data)
        return self._resort_subtracted(subtracted)

    @staticmethod
    def _resort_subtracted(subtracted_ranges):
        return sorted(subtracted_ranges, key=lambda x: x.sort_key())

    # common 'interpretations' or extractions of transcript-related data
    def transcribed_ranges(self):
        return self._ranges_by_type(types.TRANSCRIBED)

    def translated_ranges(self):
        return self._ranges_by_type(types.CODING)

    def intronic_ranges(self):
        return self._ranges_by_type(types.INTRON)

    def trans_intronic_ranges(self):
        return self._ranges_by_type(types.TRANS_INTRON)

    def error_ranges(self):
        return self._ranges_by_type(types.ERROR)

    def cis_exonic_ranges(self):  # AKA exon
        transcribeds = self._ranges_by_type(types.TRANSCRIBED)
        introns = self._ranges_by_type(types.INTRON)
        exons = self._subtract_ranges(subtract_from=transcribeds, to_subtract=introns)
        return exons

    def translated_exonic_ranges(self):  # AKA CDS
        # todo, somewhere, maybe not here, consider further consistency checking
        #  e.g. (that all CODING regions are within TRANSCRIBED regions)
        translateds = self._ranges_by_type(types.CODING)
        introns = self._ranges_by_type(types.INTRON)
        coding_exons = self._subtract_ranges(subtract_from=translateds, to_subtract=introns)
        return coding_exons

    def untranslated_exonic_ranges(self):  # AKA UTR
        exons = self.cis_exonic_ranges()
        translateds = self._ranges_by_type(types.CODING)
        utrs = self._subtract_ranges(subtract_from=exons, to_subtract=translateds)
        return utrs

    # point transitions (sites)
    def _get_by_type_and_bearing(self, target_type, target_bearing):
        out = []
        for piece in self.transcript.data.transcribed_pieces:
            for feature in piece.features:
                if feature.bearing.value == target_bearing and \
                        feature.type.value == target_type:
                    out.append(TranscriptCoordinate(coordinate_id=feature.coordinate_id,
                                                    piece_position=piece.position,
                                                    is_plus_strand=feature.is_plus_strand,
                                                    start=feature.position))
        return out

    def transcription_start_sites(self):
        return self._get_by_type_and_bearing(types.TRANSCRIBED, types.START)

    def translation_start_sites(self):  # AKA start codon
        return self._get_by_type_and_bearing(types.CODING, types.START)

    def intron_start_sites(self):  # AKA Donor splice site
        return self._get_by_type_and_bearing(types.INTRON, types.START)

    def trans_intron_start_sites(self):
        return self._get_by_type_and_bearing(types.TRANS_INTRON, types.START)

    def transcription_end_sites(self):
        return self._get_by_type_and_bearing(types.TRANSCRIBED, types.END)

    def translation_end_sites(self):  # AKA follows stop codon
        return self._get_by_type_and_bearing(types.CODING, types.END)

    def intron_end_sites(self):  # AKA follows acceptor splice site
        return self._get_by_type_and_bearing(types.INTRON, types.END)

    def trans_intron_end_sites(self):
        return self._get_by_type_and_bearing(types.TRANS_INTRON, types.END)


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
