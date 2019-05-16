import os
import intervaltree
import copy
import logging
import hashlib
import itertools
from pprint import pprint
from abc import ABC, abstractmethod
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from dustdas import gffhelper, fastahelper
from .. import orm
from .. import types
from .. import handlers
from .. import helpers
from ..base.helpers import as_py_start, as_py_end, get_strand_direction, get_geenuff_start_end


class IntervalCountError(Exception):
    pass


class NoTranscriptError(Exception):
    pass


class TransSplicingError(Exception):
    pass


class NoGFFEntryError(Exception):
    pass


# core queue prep
class InsertionQueue(helpers.QueueController):
    def __init__(self, session, engine):
        super().__init__(session, engine)
        self.super_locus = helpers.CoreQueue(orm.SuperLocus.__table__.insert())
        self.transcribed = helpers.CoreQueue(orm.Transcribed.__table__.insert())
        self.transcribed_piece = helpers.CoreQueue(orm.TranscribedPiece.__table__.insert())
        self.translated = helpers.CoreQueue(orm.Translated.__table__.insert())
        self.feature = helpers.CoreQueue(orm.Feature.__table__.insert())
        self.association_transcribed_piece_to_feature = helpers.CoreQueue(
            orm.association_transcribed_piece_to_feature.insert())
        self.association_translated_to_feature = helpers.CoreQueue(
            orm.association_translated_to_feature.insert())

        self.ordered_queues = [
            self.super_locus, self.transcribed, self.transcribed_piece, self.translated,
            self.feature, self.association_transcribed_piece_to_feature,
            self.association_translated_to_feature
        ]


class InsertCounterHolder(object):
    """provides incrementing unique integers to be used as primary keys for bulk inserts"""
    feature = helpers.Counter(orm.Feature)
    translated = helpers.Counter(orm.Translated)
    transcribed = helpers.Counter(orm.Transcribed)
    super_locus = helpers.Counter(orm.SuperLocus)
    transcribed_piece = helpers.Counter(orm.TranscribedPiece)
    genome = helpers.Counter(orm.Genome)

    @staticmethod
    def sync_counters_with_db(session):
        InsertCounterHolder.feature.sync_with_db(session)
        InsertCounterHolder.translated.sync_with_db(session)
        InsertCounterHolder.transcribed.sync_with_db(session)
        InsertCounterHolder.super_locus.sync_with_db(session)
        InsertCounterHolder.transcribed_piece.sync_with_db(session)
        InsertCounterHolder.genome.sync_with_db(session)


# todo remember to check for trans-splicing error in gff entry (see former add_gff_entry_group)
# todo remember to work in the code that was in _mark_erroneous when handling errors

class OrganizedGeenuffHandlerGroup(object):
    """Stores the Handler object for a super locus in an organized fassion.
    The format is similar to the one of OrganizedGFFEntryGroup, but it stores objects
    according to the Geenuff way of saving genomic annotations.

    The handlers are organized in the following way:

    handlers = {
        'super_locus_h' = super_locus_handler,
        'transcript_hs' = {
            transcript_feature_handler1: {
                cds_handler1: [intron_handler1, intron_handler2, ..],
                cds_handler2: [intron_handler1, intron_handler2, ..]
            },
            transcript_feature_handler2: {
                cds_handler1: [intron_handler1, intron_handler2, ..],
                cds_handler2: [intron_handler1, intron_handler2, ..]
            },
            ...
        },
        'error_h' = [error_handler1, error_handler2, ..]
    }
    """

    def __init__(self, organized_gff_entries, coord, controller, err_handle):
        self.coord = coord
        self.controller = controller
        self.err_handle = err_handle
        self.handlers = {'transcript_hs': {}, 'error_h': []}
        self._parse_gff_entries(organized_gff_entries)

    def _parse_gff_entries(self, entries):
        """Changes the GFF format into the GeenuFF format."""
        # these attr are the same for all geenuff handlers going to be created
        sl = entries['super_locus']
        is_plus_strand = get_strand_direction(sl)
        score = sl.score
        source = sl.source

        sl_start, sl_end = get_geenuff_start_end(sl.start, sl.end, is_plus_strand)
        self.handlers['super_locus_h'] = SuperLocusHandler(entry_type=sl.type,
                                                           given_name=sl.get_ID(),
                                                           coord=self.coord,
                                                           is_plus_strand=is_plus_strand,
                                                           start=sl_start,
                                                           end=sl_end,
                                                           controller=self.controller)
        for t, t_entries in entries['transcripts'].items():
            # check for multi inheritance and throw NotImplementedError if found
            if len(t.get_Parent()) > 1:
                raise NotImplementedError
            # create transcript feature handler
            t_h = FeatureHandler(self.coord,
                                 is_plus_strand,
                                 types.TRANSCRIBED,
                                 sub_type=t.type,
                                 given_name=t.get_ID(),
                                 score=score,
                                 source=source,
                                 controller=self.controller)
            t_h.set_start_end_from_gff(t.start, t.end)
            self.handlers['transcript_hs'][t_h] = {}

            # create coding features from exon limits
            cds_h = FeatureHandler(self.coord,
                                   is_plus_strand,
                                   types.CODING,
                                   # take the coding phase from the first cds entry in the gff
                                   phase=t_entries['cds'][0].phase,
                                   score=score,
                                   source=source,
                                   controller=self.controller)
            gff_start = t_entries['exons'][0].start
            gff_end = t_entries['exons'][-1].end
            cds_h.set_start_end_from_gff(gff_start, gff_end)
            self.handlers['transcript_hs'][t_h][cds_h] = []

            # create all the introns by traversing the exon entries and insert FeatureHandlers
            # into the previously created list
            # the introns should strictly always lie between successive gff exons
            exons = t_entries['exons']
            for i in range(len(exons) - 1):
                intron_h = FeatureHandler(self.coord,
                                          is_plus_strand,
                                          types.INTRON,
                                          score=score,
                                          source=source,
                                          controller=self.controller)
                # the introns are delimited by the surrounding exons
                # the first base of an intron in right after the last exononic base
                gff_start = exons[i].end + 1
                gff_end = exons[i + 1].start - 1
                intron_h.set_start_end_from_gff(gff_start, gff_end)
                self.handlers['transcript_hs'][t_h][cds_h].append(intron_h)


class OrganizedGFFEntryGroup(object):
    """Takes an entry group and stores the entries in an orderly fassion.
    Can transform the gff way of gene encoding into a OrganizedGeenuffHandlerGroup.
    Does not perform error checking, which happens later.

    The entries are organized in the following way:

    entries = {
        'super_locus' = super_locus_entry,
        'transcripts' = {
            transcript_entry1: {
                'exons': [ordered_exon_entry1, ordered_exon_entry2, ..],
                'cds': [ordered_cds_entry1, ordered_cds_entry2, ..]
            },
            transcript_entry2: {
                'exons': [ordered_exon_entry1, ordered_exon_entry2, ..],
                'cds': [ordered_cds_entry1, ordered_cds_entry2, ..]
            },
            ...
    }
    """

    def __init__(self, gff_entry_group, genome_h, controller, err_handle):
        self.genome_h = genome_h
        self.controller = controller
        self.err_handle = err_handle
        self.entries = {'transcripts': {}}
        self.coord = None
        self.add_gff_entry_group(gff_entry_group, err_handle)

    def add_gff_entry_group(self, entries, err_handle):
        latest_transcript = None
        for entry in list(entries):
            if helpers.in_enum_values(entry.type, types.SuperLocusAll):
                assert 'super_locus' not in self.entries
                self.entries['super_locus'] = entry
            elif helpers.in_enum_values(entry.type, types.TranscriptLevelAll):
                self.entries['transcripts'][entry] = {'exons': [], 'cds': []}
                latest_transcript = entry
            elif entry.type == types.EXON:
                self.entries['transcripts'][latest_transcript]['exons'].append(entry)
            elif entry.type == types.CDS:
                self.entries['transcripts'][latest_transcript]['cds'].append(entry)
            elif entry.type in [types.FIVE_PRIME_UTR, types.THREE_PRIME_UTR]:
                # ignore these features but do not throw error
                pass
            else:
                raise ValueError("problem handling entry of type {}".format(entry.type))

        # set the coordinate
        self.coord = self.genome_h.gffid_to_coords[self.entries['super_locus'].seqid]

        # order exon and cds lists by start value
        for transcript, value_dict in self.entries['transcripts'].items():
            for key in ['exons', 'cds']:
                value_dict[key].sort(key=lambda e:e.start)

    def get_geenuff_handlers(self):
        geenuff_handler_group = OrganizedGeenuffHandlerGroup(self.entries,
                                                             self.coord,
                                                             self.controller,
                                                             self.err_handle)
        return geenuff_handler_group.handlers



class GFFErrorHandling(object):
    """Deal with error detection and handling of the input features. Does the handling
    in the space of GeenuFF handlers.
    Assumes all super locus handler groups to be ordered 5p to 3p and of one strand.
    Works with a list of OrganizedGeenuffHandlerGroup, which correspond to a list of
    super loci, and looks for errors. Error features may be inserted and handlers be
    removed when deemed necessary.
    """
    def __init__(self, geenuff_handler_groups, is_plus_strand):
        self.groups = geenuff_handler_groups
        self.is_plus_strand = is_plus_strand

    def resolve_errors(self):
        for i, group in enumerate(self.groups):
            # the case of no transcript for a super locus
            # solution is to add an error mask that extends halfway to the appending super
            # loci in the intergenic region as something in this area appears to have gone wrong
            error_h = None
            if not group['transcript_hs']:
                sl_h = group['super_locus_h']
                error_start = group['super_locus_h'].start
                error_end = group['super_locus_h'].end
                if i > 0:
                    sl_h_prev = geenuff_handler_groups[i - 1]['super_locus_h']
                    error_start = self.halfway_mark(sl_h_prev, sl_h)
                if i < len(geenuff_handler_groups) - 1:
                    sl_h_next = geenuff_handler_groups[i + 1]['super_locus_h']
                    error_end = self.halfway_mark(sl_h, sl_h_next)
                error_h = self._get_error_handler(sl_h.coord,
                                                  error_start,
                                                  error_end,
                                                  is_plus_strand)
            if error_h:
                group['error_h'].append(error_h)

    def _get_error_handler(self, coord, start, end, is_plus_strand):
        # todo logging
        error_h = FeatureHandler(coord,
                                 is_plus_strand,
                                 types.ERROR,
                                 start=start,
                                 end=end,
                                 controller=self)
        return error_h

    def halfway_mark(self, sl_h, sl_h_next):
        """Calculates the half way point between two super loci, which is then used for
        error masks.
        """
        if self.is_plus_strand:
            dist = sl_h_next.start - sl_h.end
            mark = sl_h.end + dist // 2
        else:
            dist = sl_h.end - sl_h_next.start
            mark = sl_h.end - dist // 2
        return mark


##### main flow control #####
class ImportController(object):
    def __init__(self, database_path, err_path='/dev/null', replace_db=False):
        self.database_path = database_path
        self.err_path = err_path
        self.latest_genome_handler = None
        self._mk_session(replace_db)
        # queues for adding to db
        self.insertion_queues = InsertionQueue(session=self.session, engine=self.engine)

    def _mk_session(self, replace_db):
        appending_to_db = False
        if os.path.exists(self.database_path):
            if replace_db:
                print('removed existing database at {}'.format(self.database_path))
                os.remove(self.database_path)
            else:
                appending_to_db = True
                print('appending to existing database at {}'.format(self.database_path))
        self.engine = create_engine(helpers.full_db_path(self.database_path), echo=False)
        orm.Base.metadata.create_all(self.engine)
        self.session = sessionmaker(bind=self.engine)()
        if appending_to_db:
            InsertCounterHolder.sync_counters_with_db(self.session)

    def gff_gen(self, gff_file):
        known = [x.value for x in types.AllKnown]
        reader = gffhelper.read_gff_file(gff_file)
        for entry in reader:
            if entry.type not in known:
                raise ValueError("unrecognized feature type from gff: {}".format(entry.type))
            else:
                self.clean_entry(entry)
                yield entry

    # todo, this is an awkward place for this
    @staticmethod
    def clean_entry(entry):
        # always present and integers
        entry.start = int(entry.start)
        entry.end = int(entry.end)
        # clean up score
        if entry.score == '.':
            entry.score = None
        else:
            entry.score = float(entry.score)

        # clean up phase
        if entry.phase == '.':
            entry.phase = None
        else:
            entry.phase = int(entry.phase)
        assert entry.phase in [None, 0, 1, 2]

        # clean up strand
        if entry.strand == '.':
            entry.strand = None
        else:
            assert entry.strand in ['+', '-']

    def useful_gff_entries(self, gff_file):
        skipable = [x.value for x in types.IgnorableFeatures]
        reader = self.gff_gen(gff_file)
        for entry in reader:
            if entry.type not in skipable:
                yield entry

    def group_gff_by_gene(self, gff_file):
        gene_level = [x.value for x in types.SuperLocusAll]
        reader = self.useful_gff_entries(gff_file)
        gene_group = [next(reader)]
        for entry in reader:
            if entry.type in gene_level:
                yield gene_group
                gene_group = [entry]
            else:
                gene_group.append(entry)
        yield gene_group

    def make_genome(self, genome_args=None):
        if type(genome_args) == dict:
            genome = orm.Genome(**genome_args)
        else:
            genome = orm.Genome()
        self.latest_genome_handler = GenomeHandler(genome)
        self.session.add(genome)
        self.session.commit()

    def add_genome(self, fasta_path, gff_path, genome_args=None, clean_gff=True):
        err_handle = open(self.err_path, 'w')
        self.clean_tmp_data()
        self.add_sequences(fasta_path, genome_args)
        self.add_gff(gff_path, err_handle, clean=clean_gff)
        err_handle.close()

    def add_sequences(self, seq_path, genome_args=None):
        if self.latest_genome_handler is None:
            self.make_genome(genome_args)
        self.latest_genome_handler.add_sequences(seq_path)
        self.session.commit()

    def clean_tmp_data(self):
        self.latest_genome_handler = None
        self.latest_super_loci = []

    def add_gff(self, gff_file, err_handle, clean=True):
        def insert_handler_groups(self, groups):
            """Initiates the calling of the add_to_queue() function of the handlers
            in the correct order.
            """
            for group in groups:
                group['super_locus_h'].add_to_queue(self.insertion_queues)
                # insert all features as well as transcript and protein related entries
                for t_feature_h in group['transcript_hs']:
                    t_h, tp_h = t_feature_h.get_transcript_handlers(group['super_locus_h'])
                    t_h.add_to_queue(self.insertion_queues)
                    tp_h.add_to_queue(self.insertion_queues)
                    t_feature_h.add_to_queue(self.insertion_queues, tp_h)
                    for cds_h in group['transcript_hs'][t_feature_h]:
                        cds_h.add_to_queue(self.insertion_queues, tp_h)
                        for intron_h in group['transcript_hs'][t_feature_h][cds_h]:
                            intron_h.add_to_queue(self.insertion_queues, tp_h)
                # insert the errors
                for e_h in group['error_h']:
                    e_h.add_to_queue(self.insertion_queues)

        def clean_and_insert(self, groups, clean):
            if clean:
                # check and correct for errors
                # do so for each strand seperately
                # all changes should be made by reference
                plus = [g for g in groups if g['super_locus_h'].is_plus_strand]
                GFFErrorHandling(plus, True).resolve_errors()
                # reverse order on minus strand
                minus = [g for g in groups if not g['super_locus_h'].is_plus_strand]
                GFFErrorHandling(minus[::-1], False).resolve_errors()
            # insert handlers
            insert_handler_groups(self, plus)
            insert_handler_groups(self, minus)
            self.insertion_queues.execute_so_far()

        assert self.latest_genome_handler is not None, 'No recent genome found'
        self.latest_genome_handler.mk_mapper(gff_file)
        geenuff_handler_groups = []
        for i, entry_group in enumerate(self.group_gff_by_gene(gff_file)):
            organized_entries = OrganizedGFFEntryGroup(entry_group,
                                                       self.latest_genome_handler,
                                                       self,
                                                       err_handle)
            geenuff_handler_groups.append(organized_entries.get_geenuff_handlers())
            if i > 0 and i % 500 == 0:
                clean_and_insert(self, geenuff_handler_groups, clean)
                geenuff_handler_groups = []
        clean_and_insert(self, geenuff_handler_groups, clean)



class Insertable(ABC):

    @abstractmethod
    def add_to_queue(self, insertion_queues):
        pass


class GenomeHandler(handlers.GenomeHandlerBase):
    def __init__(self, data=None, controller=None):
        super().__init__(data)
        if controller is not None:
            self.id = InsertCounterHolder.genome()
        self.mapper = None
        self._coords_by_seqid = None
        self._gffid_to_coords = None
        self._gff_seq_ids = None

    @property
    def data_type(self):
        return orm.Genome

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
            if seq[0].islower():
                seq = seq.upper()  # this may perform poorly
            seqid = fasta_header.split(id_delim)[0]
            # todo, parallelize sequence & annotation format, then import directly from ~Slice
            coord = orm.Coordinate(sequence=seq,
                                   start=0,
                                   end=len(seq),
                                   seqid=seqid,
                                   sha1=helpers.sequence_hash(seq),
                                   genome=self.data)


class SuperLocusHandler(handlers.SuperLocusHandlerBase, Insertable):
    def __init__(self, entry_type, given_name, coord=None, is_plus_strand=None, start=-1, end=-1, controller=None):
        super().__init__()
        if controller is not None:
            self.id = InsertCounterHolder.super_locus()
        self.entry_type = entry_type
        self.given_name = given_name
        # not neccessary for insert but helpful for certain error cases
        self.coord= coord
        self.is_plus_strand = is_plus_strand
        self.start = start
        self.end = end

    def add_to_queue(self, insertion_queues):
        to_add = {'type': self.entry_type, 'given_name': self.given_name, 'id': self.id}
        insertion_queues.super_locus.queue.append(to_add)

    def __repr__(self):
        params = {'id': self.id, 'type': self.entry_type, 'given_name': self.given_name}
        return self._get_repr('SuperLocusHandler', params)


class FeatureHandler(handlers.FeatureHandlerBase, Insertable):
    def __init__(self, coord, is_plus_strand, feature_type, controller, sub_type=None, start=-1, end=-1, given_name=None, phase=0, score=None, source=None):
        """Initializes a handler for a soon to be inserted geenuff feature.
        As there is no 1to1 relation between a gff entry and a geenuff feature,
        no gff entry is saved or directly used to infer attributes from.
        """
        super().__init__()
        self.id = InsertCounterHolder.feature()
        self.coord = coord
        self.given_name = given_name
        self.is_plus_strand = is_plus_strand
        self.feature_type = feature_type
        self.sub_type = sub_type  # only used if type is transcribed to save the transcript type
        # start/end may have to be adapted to geenuff
        self.start = start
        self.end = end
        self.phase = phase
        self.score = score
        self.source = source
        self.start_is_biological_start = None
        self.end_is_biological_end = None
        self.controller = controller

    def add_to_queue(self, insertion_queues, transcript_piece_h=None):
        """Adds the feature as well as the many2many association entries if the type is not
        types.ERROR"""
        feature = {
            'id': self.id,
            'type': self.feature_type,
            'given_name': self.given_name,
            'coordinate_id': self.coord.id,
            'is_plus_strand': self.is_plus_strand,
            'score': self.score,
            'source': self.source,
            'phase': self.phase,
            'start': self.start,
            'end': self.end,
            'start_is_biological_start': self.start_is_biological_start,
            'end_is_biological_end': self.end_is_biological_end,
        }
        insertion_queues.feature.queue.append(feature)

        if not self.feature_type == types.ERROR:
            features2pieces = {
                'feature_id': self.id,
                'transcribed_piece_id': transcript_piece_h.id,
            }
            insertion_queues.association_transcribed_piece_to_feature.queue.append(features2pieces)

    def get_transcript_handlers(self, super_locus_h):
        """Returns a TranscribedHandler and TranscribedPieceHandler instance that corresponds
        to the feature if the type is types.TRANSCRIBED. These can then be inserted.
        """
        transcript_h = TranscribedHandler(entry_type=self.sub_type,
                                          given_name=self.given_name,
                                          super_locus_id=super_locus_h.id,
                                          controller=self.controller)
        transcript_piece_h = TranscribedPieceHandler(given_name=self.given_name,
                                                     transcript_id=transcript_h.id,
                                                     position=0,
                                                     controller=self.controller)
        return transcript_h, transcript_piece_h


    def set_start_end_from_gff(self, gff_start, gff_end):
        self.start, self.end = get_geenuff_start_end(gff_start, gff_end, self.is_plus_strand)

    def pos_cmp_key(self):
        sortable_start = self.start
        sortable_end = self.end
        if not self.is_plus_strand:
            sortable_start = sortable_start * -1
            sortable_end = sortable_end * -1
        return self.coord.seqid, self.is_plus_strand, sortable_start, sortable_end, self.type

    def __repr__(self):
        params = {
            'id': self.id,
            'coord_id': self.coord.id,
            'type': self.type,
            'is_plus_strand': self.is_plus_strand,
            'phase': self.phase,
        }
        if self.sub_type:
            params['sub_type'] = self.sub_type
        if self.given_name:
            params['given_name'] = self.given_name
        return self._get_repr('FeatureHandler', params, str(self.start) + '--' + str(self.end))


class TranscribedHandler(handlers.TranscribedHandlerBase, Insertable):
    def __init__(self, entry_type, given_name, super_locus_id, controller):
        super().__init__()
        self.id = InsertCounterHolder.transcribed()
        self.entry_type = entry_type
        self.given_name = given_name
        self.super_locus_id = super_locus_id
        self.controller = controller

    def add_to_queue(self, insertion_queues):
        transcribed = {
            'type': self.entry_type,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id,
            'id': self.id
        }
        insertion_queues.transcribed.queue.append(transcribed)


class TranscribedPieceHandler(handlers.TranscribedPieceHandlerBase, Insertable):
    def __init__(self, given_name, transcript_id, position, controller):
        super().__init__()
        self.id = InsertCounterHolder.transcribed_piece()
        self.given_name = given_name
        self.transcript_id = transcript_id
        self.position = position

    def add_to_queue(self, insertion_queues):
        transcribed_piece = {
            'id': self.id,
            'given_name': self.given_name,
            'transcribed_id': self.transcript_id,
            'position': self.position,
        }
        insertion_queues.transcribed_piece.queue.append(transcribed_piece)


class TranslatedHandler(handlers.TranslatedHandlerBase):
    def __init__(self, given_name, super_locus_id, controller):
        super().__init__()
        self.given_name = given_name
        self.super_locus_id = super_locus_id
        self.id = InsertCounterHolder.translated()

    def add_to_queue(self, insertion_queues):
        translated = {
            'id': self.id,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id
        }
        insertion_queues.translated.queue.append(translated)
