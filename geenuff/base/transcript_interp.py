import copy
import intervaltree

from . import orm
from . import types
from .handlers import TranscribedHandlerBase


class ChannelTracker(object):
    """Tracks status of 'type' while stepping through the transitions of a transcript"""
    def __init__(self):
        self._in_region = False
        self._open_feature = None

    def set_channel_open(self, **kwargs):
        self.in_region = True

    def set_channel_closed(self):
        self.in_region = False

    @property
    def in_region(self):
        return self._in_region

    @in_region.setter
    def in_region(self, x):
        if isinstance(x, bool):
            self._in_region = x
        else:
            raise ValueError("in_region supports only boolean values {} of type {} found".format(x, type(x)))

    @property
    def open_feature(self):
        return self._open_feature

    @open_feature.setter
    def open_feature(self, feature):
        if isinstance(feature, dict) or feature is None:
            self._open_feature = feature
        else:
            raise ValueError(
                'expected instance of dict or None for open_feature, got {} of type {}'.format(feature, type(feature)))

    def update_and_close_feature(self, update_dict):
        self.open_feature.update(update_dict)
        self.open_feature = None
        self.set_channel_closed()

    def set_channel_open_with_feature(self, feature, **kwargs):
        self.set_channel_open(**kwargs)
        self.open_feature = feature


class CodingChannelTracker(ChannelTracker):
    """tracks coding status (+ extras) while stepping through transitions of a transcript"""
    def __init__(self, phase=None):
        super().__init__()
        self.phase = phase  # todo, proper tracking / handling

    def set_channel_open(self, phase):
        super().set_channel_open()
        self.phase = phase

    def set_channel_closed(self):
        super().set_channel_closed()
        self.phase = None


class TranscriptStatusBase(object):
    """can hold and manipulate all the info on current status (later: feature types) of a transcript"""

    def __init__(self):
        # initializes to intergenic (all channels / types set to False)
        self.transcribed_tracker = ChannelTracker()
        self.intron_tracker = ChannelTracker()
        self.trans_intron_tracker = ChannelTracker()
        self.coding_tracker = CodingChannelTracker()
        self.error_tracker = ChannelTracker()

    def __repr__(self):
        return "in_transcribed: {}, in_intron: {}, in_coding: {}, in_trans_intron: {}, phase: {}".format(
            self.transcribed_tracker.in_region, self.intron_tracker.in_region, self.coding_tracker.in_region,
            self.trans_intron_tracker.in_region, self.coding_tracker.phase
        )

    def is_utr(self):
        return self.transcribed_tracker.in_region and \
               not any([self.intron_tracker.in_region,
                        self.coding_tracker.in_region,
                        self.trans_intron_tracker.in_region])

    def is_coding(self):
        return self.transcribed_tracker.in_region and \
               self.coding_tracker.in_region and \
               not any([self.intron_tracker.in_region, self.trans_intron_tracker.in_region])

    def is_intronic(self):
        return self.intron_tracker.in_region and \
               self.transcribed_tracker.in_region

    def is_trans_intronic(self):
        return self.trans_intron_tracker and \
               self.transcribed_tracker.in_region

    def is_intergenic(self):
        return not self.transcribed_tracker.in_region


class TranscriptCoordinate(object):
    """holds (and helps sort) either start or end, with the sequence, piece position, and direction"""
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
        return "coordinate: {}, piece position {}, is_plus {}: {}".format(self.coordinate_id,
                                                                          self.piece_position,
                                                                          self.is_plus_strand,
                                                                          self.start)

    def __eq__(self, other):
        if isinstance(other, TranscriptCoordinate):
            return self.__dict__ == other.__dict__
        return False


class Range(TranscriptCoordinate):
    """holds (and helps sort) a start-end range with the sequence, piece position, and direction"""
    def __init__(self, coordinate_id, piece_position, start, end, is_plus_strand):
        super().__init__(coordinate_id=coordinate_id,
                         piece_position=piece_position,
                         is_plus_strand=is_plus_strand,
                         start=start)
        self.end = end

    def __repr__(self):
        return "coordinate: {}, piece position {}, is_plus {}: {}-{}".format(self.coordinate_id,
                                                                             self.piece_position,
                                                                             self.is_plus_strand,
                                                                             self.start,
                                                                             self.end)


class Position(object):
    """organized data holder for passing start/end coordinates around"""
    def __init__(self, position, is_biological):
        self.position = position
        self.is_biological = is_biological


class EukCodingChannelTracker(CodingChannelTracker):
    def __init__(self, phase=None):
        super().__init__(phase=phase)
        self.seen_start = False
        self.seen_stop = False

    def set_channel_open(self, phase):
        super().set_channel_open(phase)
        self.seen_start = True

    def set_channel_closed(self):
        super().set_channel_closed()
        self.seen_stop = True


class EukTranscriptStatus(TranscriptStatusBase):
    """can hold and manipulate all the info on current status of a Eukaryotic transcript"""
    def __init__(self):
        super().__init__()
        # initializes to intergenic (all channels / types set to False)
        self.coding_tracker = EukCodingChannelTracker()

    def is_5p_utr(self):
        return self.is_utr() and not any([self.coding_tracker.seen_start, self.coding_tracker.seen_stop])

    def is_3p_utr(self):
        return self.is_utr() and self.coding_tracker.seen_start and self.coding_tracker.seen_stop


def positional_match(feature, previous):
    return feature.pos_cmp_key() == previous.pos_cmp_key()


class TranscriptInterpBase(object):
    """handles basics from transitioning (piece & feature sorting) through biological interpretation of a transcript"""
    def __init__(self, transcript, super_locus, session=None):
        assert isinstance(transcript, TranscribedHandlerBase)
        self.status = TranscriptStatusBase()
        self.transcript = transcript
        self.session = session
        self.super_locus = super_locus

    def transition_5p_to_3p(self):
        for piece in self.sorted_pieces():
            piece_features = self.sorted_features(piece)
            for aligned_features in self.full_stack_matches(piece_features):
                yield aligned_features, piece

    def sorted_features(self, piece):
        features = piece.features
        # confirm strand & seqid
        assert all([f.coordinate.seqid == features[0].coordinate.seqid for f in features]), \
            'on {}, piece {}, not all matching: {}'.format(self.transcript.data, (piece.id, piece.position), [(f.id, f.coordinate) for f in features])
        assert all([f.is_plus_strand == features[0].is_plus_strand for f in features]), \
            'on {}, piece {}, not all matching: {}'.format(self.transcript.data, (piece.id, piece.position), [(f.id, f.is_plus_strand) for f in features])
        features = sorted(features, key=lambda x: x.pos_cmp_key())
        return features

    def sorted_pieces(self):
        pieces = self.transcript.data.transcribed_pieces
        ordered_pieces = sorted(pieces, key=lambda x: x.position)
        return ordered_pieces

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

    def full_stack_matches(self, features):
        for matches in self.stack_matches(features, match_fn=positional_match):
            yield matches

    def sort_all(self):
        out = []
        for piece in self.sorted_pieces():
            out.append(self.sorted_features(piece))
        return out

    # helpers for classic transitions below
    def _ranges_by_type(self, target_type):
        ranges = []
        for aligned_features, piece in self.transition_5p_to_3p():
            features_of_type = [f for f in aligned_features if f.type.value == target_type]
            assert len(features_of_type) in [0, 1], "cannot interpret aligned features of the same type {}".format(
                features_of_type
            )
            if len(features_of_type) == 1:  # and 0 is simply ignored...
                feature = features_of_type[0]
                ranges.append(Range(coordinate_id=feature.coordinate_id,
                                    start=feature.start,
                                    end=feature.end,
                                    is_plus_strand=feature.is_plus_strand,
                                    piece_position=piece.position))

        return ranges

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
    @staticmethod
    def _get_transition(feature, target_start_not_end):
        if target_start_not_end:
            at = feature.start
            is_bio = feature.start_is_biological_start
        else:
            at = feature.end
            is_bio = feature.end_is_biological_end
        return at, is_bio

    def get_by_type_and_bearing(self, target_type, target_start_not_end, target_is_biological=True):
        out = []
        for piece in self.transcript.data.transcribed_pieces:
            for feature in piece.features:
                if feature.type.value == target_type:
                    at, is_bio = self._get_transition(feature, target_start_not_end)
                    if is_bio == target_is_biological:
                        out.append(TranscriptCoordinate(coordinate_id=feature.coordinate_id,
                                                        piece_position=piece.position,
                                                        is_plus_strand=feature.is_plus_strand,
                                                        start=at))
        return out

    def transcription_start_sites(self):
        return self.get_by_type_and_bearing(types.TRANSCRIBED, target_start_not_end=True)

    def translation_start_sites(self):  # AKA start codon
        return self.get_by_type_and_bearing(types.CODING, target_start_not_end=True)

    def intron_start_sites(self):  # AKA Donor splice site
        return self.get_by_type_and_bearing(types.INTRON, target_start_not_end=True)

    def trans_intron_start_sites(self):
        return self.get_by_type_and_bearing(types.TRANS_INTRON, target_start_not_end=True)

    def transcription_end_sites(self):
        return self.get_by_type_and_bearing(types.TRANSCRIBED, target_start_not_end=False)

    def translation_end_sites(self):  # AKA follows stop codon
        return self.get_by_type_and_bearing(types.CODING, target_start_not_end=False)

    def intron_end_sites(self):  # AKA follows acceptor splice site
        return self.get_by_type_and_bearing(types.INTRON, target_start_not_end=False)

    def trans_intron_end_sites(self):
        return self.get_by_type_and_bearing(types.TRANS_INTRON, target_start_not_end=False)
