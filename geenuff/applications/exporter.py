import sys
import copy
import intervaltree

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from geenuff.base.orm import Coordinate, Genome, Feature, Transcript, TranscriptPiece, \
    association_transcript_piece_to_feature, SuperLocus
from geenuff.base.handlers import TranscriptHandlerBase, SuperLocusHandlerBase
from geenuff.base.helpers import full_db_path, Counter
from geenuff.base import types


class ExportController(object):
    def __init__(self, db_path_in, longest=False):
        self.db_path_in = db_path_in
        self._mk_session()
        self.longest = longest
        self.id_counter = Counter()
        self.export_ranges = []

    def _mk_session(self):
        self.engine = create_engine(full_db_path(self.db_path_in), echo=False)
        self.session = sessionmaker(bind=self.engine)()

    def _check_genome_names(self, *argv):
        for names in argv:
            if names:
                genome_ids = self.session.query(Genome.id).filter(Genome.species.in_(names)).all()
                if len(genome_ids) != len(names):
                    print('One or more of the given genome names can not be found in the database')
                    exit()

    def _coords_with_feature_query(self):
        return self.session.query(Feature.coordinate_id).distinct()

    def _all_coords_query(self):
        return self.session.query(Coordinate.id)

    def _get_coords_by_genome_query(self, genomes, exclude, include_without_features=False):
        if include_without_features:
            coordinate_ids_of_interest = self._all_coords_query()
        else:
            coordinate_ids_of_interest = self._coords_with_feature_query()
        if genomes:
            print('Selecting the following genomes: {}'.format(genomes), file=sys.stderr)
            all_coord_ids = (self.session.query(Coordinate.id)
                             .join(Genome, Genome.id == Coordinate.genome_id)
                             .filter(Genome.species.in_(genomes))
                             .filter(Coordinate.id.in_(coordinate_ids_of_interest)))
        else:
            if exclude:
                print('Selecting all genomes from {} except: {}'.format(self.db_path_in, exclude),
                      file=sys.stderr)
                all_coord_ids = (self.session.query(Coordinate.id)
                                 .join(Genome, Genome.id == Coordinate.genome_id)
                                 .filter(Genome.species.notin_(exclude))
                                 .filter(Coordinate.id.in_(coordinate_ids_of_interest)))
            else:
                print('Selecting all genomes from {}'.format(self.db_path_in),
                      file=sys.stderr)
                all_coord_ids = coordinate_ids_of_interest

        return all_coord_ids

    def get_super_loci_by_coords(self, coord_ids):
        sess = self.session  # shortcut
        # surprisingly not that bad, a minute or two for 55 genomes
        # todo, join instead of IN and keep results
        q = sess.query(Feature.id).filter(Feature.coordinate_id.in_(coord_ids))
        q = sess.query(association_transcript_piece_to_feature.c.transcript_piece_id)\
            .filter(association_transcript_piece_to_feature.c.feature_id.in_(q)).distinct()
        q = sess.query(TranscriptPiece.transcript_id)\
            .filter(TranscriptPiece.id.in_(q)).distinct()
        q = sess.query(Transcript.super_locus_id)\
            .filter(Transcript.id.in_(q)).distinct()
        return q

    def gen_ranges(self, genomes, exclude, range_function):
        # which genomes to export
        coord_ids = self._get_coords_by_genome_query(genomes, exclude)
        super_loci = self.get_super_loci_by_coords(coord_ids).all()
        i = 0
        for sl in super_loci:
            super_locus = self.session.query(SuperLocus).filter(SuperLocus.id == sl[0]).first()
            sl_ranger = SuperLocusRanger(super_locus, longest=self.longest)
            # todo, once JOIN output exists, drop all these loops
            for range_maker in sl_ranger.exp_range_makers:
                export_groups = range_function(range_maker)
                for group in export_groups:
                    if group.seqid is None:
                        group.seqid = 'unnamed_{0:08d}'.format(self.id_counter())
                    yield group
                    i += 1

    def prep_ranges(self, genomes, exclude, range_function):
        for arange in self.gen_ranges(genomes, exclude, range_function):
            self.export_ranges.append(arange)


def positional_match(feature, previous):
    return feature.pos_cmp_key() == previous.pos_cmp_key()


class TranscriptCoordinate(object):
    """holds (and helps sort) either start or end, with the sequence, piece position, and direction"""
    def __init__(self, coordinate_id, piece_position, is_plus_strand, start):
        self.start = start
        self.coordinate_id = coordinate_id
        self.piece_position = piece_position
        self.is_plus_strand = is_plus_strand

    def _sort_pos(self, pos):
        if self.is_plus_strand:
            sort_pos = pos
        else:
            sort_pos = -pos  # flip sort order on the - strand
        return sort_pos

    def sort_key(self):
        return self.piece_position, self._sort_pos(self.start)

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
    def __init__(self, coordinate_id, piece_position, start, end, is_plus_strand, given_name=None):
        super().__init__(coordinate_id=coordinate_id,
                         piece_position=piece_position,
                         is_plus_strand=is_plus_strand,
                         start=start)
        self.end = end
        self.given_name = given_name

    def sequence_chunk_info(self):
        return self.coordinate_id, self.piece_position, self.is_plus_strand

    def sort_key(self):
        return self.piece_position, self._sort_pos(self.start), self._sort_pos(self.end)

    def __repr__(self):
        return "coordinate: {}, piece position {}, is_plus {}: {}-{}".format(self.coordinate_id,
                                                                             self.piece_position,
                                                                             self.is_plus_strand,
                                                                             self.start,
                                                                             self.end)


class ExportGroup(object):
    """Holds a named list of ordered ranges"""
    def __init__(self, seqid, ranges=None):
        if ranges is None:
            self.ranges = []
        else:
            self.ranges = ranges
        self.seqid = seqid


class RangeMaker(TranscriptHandlerBase):
    """Interprets a transcript as ordered flattened ranges from its features"""

    def feature_piece_pairs(self):
        for piece in self.data.transcript_pieces:
            for feature in piece.features:
                yield feature, piece

    # helpers for classic transitions below
    def _ranges_by_type(self, target_type):
        ranges = []
        for feature, piece in self.feature_piece_pairs():
            if feature.type.value == target_type:  # and 0 is simply ignored...
                ranges.append(Range(coordinate_id=feature.coordinate_id,
                                    start=feature.start,
                                    end=feature.end,
                                    is_plus_strand=feature.is_plus_strand,
                                    piece_position=piece.position,
                                    given_name=feature.given_name))

        return ranges

    @staticmethod
    def _make_trees(ranges):
        trees = {}
        for r in ranges:
            coord_isplus = r.sequence_chunk_info()
            if coord_isplus not in trees:
                trees[coord_isplus] = intervaltree.IntervalTree()
            # todo, is that sufficient, do we not need to add one to -strand coordinates?
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
        if not subtract_trees:
            return self._resort_subtracted(subtract_from)
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

    @staticmethod
    def _one_range_one_group(ranges):
        return [ExportGroup(seqid=r.given_name, ranges=[r]) for r in ranges]

    # common 'interpretations' or extractions of transcript-related data
    # all of the following methods should return a ready "ExportGroup" that has all the ordered ranges
    # that need to be combined to form a sequence, and an id for this sequence
    def transcribed_ranges(self):
        return [ExportGroup(seqid=self.data.given_name, ranges=self._ranges_by_type(types.GEENUFF_TRANSCRIPT))]

    def cds_ranges(self):
        return self._one_range_one_group(self._ranges_by_type(types.GEENUFF_CDS))

    def intronic_ranges(self):
        return self._one_range_one_group(self._ranges_by_type(types.GEENUFF_INTRON))

    def exonic_ranges(self):  # AKA exon
        transcribeds = self._ranges_by_type(types.GEENUFF_TRANSCRIPT)
        introns = self._ranges_by_type(types.GEENUFF_INTRON)
        print(introns)
        exons = self._subtract_ranges(subtract_from=transcribeds, to_subtract=introns)
        return self._one_range_one_group(exons)

    def cds_exonic_ranges(self):  # AKA CDS
        # todo, somewhere, maybe not here, consider further consistency checking
        #  e.g. (that all CODING regions are within TRANSCRIBED regions)
        # todo, return separately if CDS features are connected to different proteins
        geenuff_cds = self._ranges_by_type(types.GEENUFF_CDS)
        introns = self._ranges_by_type(types.GEENUFF_INTRON)
        coding_exons = self._subtract_ranges(subtract_from=geenuff_cds, to_subtract=introns)
        return self._one_range_one_group(coding_exons)

    def untranslated_exonic_ranges(self):  # AKA UTR
        exons = self.exonic_ranges()
        geenuff_cds = self._ranges_by_type(types.GEENUFF_CDS)
        utrs = self._subtract_ranges(subtract_from=exons, to_subtract=geenuff_cds)
        return self._one_range_one_group(utrs)

    def mature_RNA(self):
        exons = self.exonic_ranges()
        return [ExportGroup(seqid=self.data.given_name, ranges=[x.ranges[0] for x in exons])]

    def mature_CDS(self):
        # todo, operon logic!!
        cds = self.cds_exonic_ranges()
        return [ExportGroup(seqid=self.data.given_name, ranges=[x.ranges[0] for x in cds])]

    def mature_UTR(self):
        # subtract CDS
        # for each unprocessed UTR
        #  --> export group
        #  subtract introns
        transcript = self._ranges_by_type(types.GEENUFF_TRANSCRIPT)
        coding = self._ranges_by_type(types.GEENUFF_CDS)
        introns = self._ranges_by_type(types.GEENUFF_INTRON)
        pre_utrs = self._subtract_ranges(subtract_from=transcript, to_subtract=coding)
        out = []
        i = 0
        for pre_utr in pre_utrs:
            utr = self._subtract_ranges(subtract_from=[pre_utr], to_subtract=introns)
            out.append(ExportGroup(seqid='{}_UTR{:02d}'.format(self.data.given_name, i), ranges=utr))
            i += 1

    def utr3p(self):
        pass  # todo

    def utr5p(self):
        pass  # todo

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
        for piece in self.data.transcribed_pieces:
            for feature in piece.features:
                if feature.type.value == target_type:
                    at, is_bio = self._get_transition(feature, target_start_not_end)
                    if is_bio == target_is_biological:
                        out.append(TranscriptCoordinate(coordinate_id=feature.coordinate_id,
                                                        piece_position=piece.position,
                                                        is_plus_strand=feature.is_plus_strand,
                                                        start=at))
        return out

    def transcript_start_sites(self):
        return self.get_by_type_and_bearing(types.GEENUFF_TRANSCRIPT, target_start_not_end=True)

    def cds_start_sites(self):  # AKA start codons
        return self.get_by_type_and_bearing(types.GEENUFF_CDS, target_start_not_end=True)

    def intron_start_sites(self):  # AKA Donor splice site
        return self.get_by_type_and_bearing(types.GEENUFF_INTRON, target_start_not_end=True)

    def transcript_end_sites(self):
        return self.get_by_type_and_bearing(types.GEENUFF_TRANSCRIPT, target_start_not_end=False)

    def cds_end_sites(self):  # AKA follows stop codons
        return self.get_by_type_and_bearing(types.GEENUFF_CDS, target_start_not_end=False)

    def intron_end_sites(self):  # AKA follows acceptor splice site
        return self.get_by_type_and_bearing(types.GEENUFF_INTRON, target_start_not_end=False)

    @staticmethod
    def _sum_range_lengths(ranges):
        out = 0
        for arange in ranges:
            out += abs(arange.end - arange.start)  # abs should make it work on - strand
        return out

    def sum_exonic_lengths(self):
        ranges = [x.ranges[0] for x in self.exonic_ranges()]
        return self._sum_range_lengths(ranges)

    def sum_exonic_cds_lengths(self):
        ranges = [x.ranges[0] for x in self.cds_exonic_ranges()]
        return self._sum_range_lengths(ranges)


MODES = {"mRNA": RangeMaker.mature_RNA,
         "pre-mRNA": RangeMaker.transcribed_ranges,
         "CDS": RangeMaker.mature_CDS,
         "exons": RangeMaker.exonic_ranges,
         "introns": RangeMaker.intronic_ranges}


class SuperLocusRanger(SuperLocusHandlerBase):
    def __init__(self, data=None, longest=False, setup_range_makers=True):
        super().__init__(data)
        self.logest = longest
        self.range_makers = []
        self.exp_range_makers = []
        if setup_range_makers:
            self.setup_range_makers(longest)

    def setup_range_makers(self, longest):
        for transcript in self.data.transcripts:
            range_maker = RangeMaker(transcript)
            self.range_makers.append(range_maker)
        if not longest:
            self.exp_range_makers = self.range_makers
        else:
            long_transcript, _ = self.get_longest_transcript()
            self.exp_range_makers = [long_transcript]

    def get_longest_transcript(self):
        """identify which transcript in this super locus is longest (with introns removed)"""
        transcript, length = None, 0
        for range_maker in self.range_makers:
            rm_length = range_maker.sum_exonic_lengths()
            if rm_length > length:
                transcript = range_maker
                length = rm_length
        return transcript, length

    def get_longest_protein_in_transcript(self):
        """identify which transcript, protein_id makes longest final coding sequence (introns rm) in this super locus"""
        pass
