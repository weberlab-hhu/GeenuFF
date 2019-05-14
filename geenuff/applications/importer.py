import os
import intervaltree
import copy
import logging
import hashlib
import itertools
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from dustdas import gffhelper, fastahelper
from .. import orm
from .. import types
from .. import handlers
from ..base.transcript_interp import TranscriptInterpBase, EukTranscriptStatus
from .. import helpers


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
            transcript_handler1: {
                cds_handler1: [intron_handler1, intron_handler2, ..],
                cds_handler2: [intron_handler1, intron_handler2, ..]
            },
            transcript_handler2: {
                cds_handler1: [intron_handler1, intron_handler2, ..],
                cds_handler2: [intron_handler1, intron_handler2, ..]
            },
            ...
    }
    """

    def __init__(self, organized_gff_entries, coord, controller, err_handle):
        self.coord = coord
        self.controller = controller
        self.err_handle = err_handle
        self.handlers = {'transcript_hs': {}}
        self.parse_gff_entries(organized_gff_entries)

    def parse_gff_entries(self, entries):
        """Changes the GFF format into the GeenuFF format."""
        # these attr are the same for all geenuff handlers going to be created
        is_plus_strand = self.get_strand_direction(entries['super_locus'])
        score = entries['super_locus'].score
        source = entries['super_locus'].source

        self.handlers['super_locus_h'] = SuperLocusHandler(entries['super_locus'], self.controller)
        for t, t_entries in entries['transcripts'].items():
            # create transcript feature handler
            t_h = FeatureHandler(self.coord,
                                 is_plus_strand,
                                 types.TRANSCRIBED,
                                 score=score,
                                 source=source,
                                 controller=self.controller)
            t_h.start = t.start
            t_h.end = t.end
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
            cds_h.start = t_entries['exons'][0].start
            cds_h.end = t_entries['exons'][-1].end
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
                intron_h.start = exons[i].end
                intron_h.end = exons[i + 1].start
                self.handlers['transcript_hs'][t_h][cds_h].append(intron_h)

    def adapt_start_end(self):
        """Adapts the start/end values according to the coordinate and geenuff specs"""
        pass

    @staticmethod
    def get_strand_direction(gffentry):
        if gffentry.strand == '+':
            return True
        elif gffentry.strand == '-':
            return False
        else:
            raise ValueError('cannot interpret strand "{}"'.format(gffentry.strand))

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
        unchecked_handlers = geenuff_handler_group.handlers
        pass


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
        # final prepping of seqid match up
        assert self.latest_genome_handler is not None, 'No recent genome found'
        self.latest_genome_handler.mk_mapper(gff_file)

        for i, entry_group in enumerate(self.group_gff_by_gene(gff_file)):
            organized_entries = OrganizedGFFEntryGroup(entry_group,
                                                       self.latest_genome_handler,
                                                       self,
                                                       err_handle)
            geenuff_handlers = organized_entries.get_geenuff_handlers()
            if clean:
                # check and correct for errors
                pass
            # insert
            if not i % 500:
                self.insertion_queues.execute_so_far()
        self.insertion_queues.execute_so_far()


##### gff parsing subclasses #####
class GFFDerived(object):
    def __init__(self, gffentry=None):
        self.gffentry = gffentry

    def add_to_queue(self, insertion_queues, **kwargs):
        raise NotImplementedError

    def setup_insertion_ready(self, **kwargs):
        # should create 'data' object (annotations_orm.Base subclass) and then call self.add_data(data)
        raise NotImplementedError


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


class SuperLocusHandler(handlers.SuperLocusHandlerBase, GFFDerived):
    def __init__(self, gffentry, controller=None):
        super().__init__()
        GFFDerived.__init__(self, gffentry)
        if controller is not None:
            self.id = InsertCounterHolder.super_locus()
        self.transcribed_handlers = []

    def setup_insertion_ready(self):
        # todo, grab more aliases from gff attribute?
        return {'type': self.gffentry.type, 'given_name': self.gffentry.get_ID(), 'id': self.id}

    def add_to_queue(self, insertion_queues, **kwargs):
        to_add = self.setup_insertion_ready()
        for key in kwargs:
            to_add[key] = kwargs[key]

        insertion_queues.super_locus.queue.append(to_add)

    @property
    def given_name(self):
        return self.gffentry.get_ID()

class FeatureHandler(handlers.FeatureHandlerBase, GFFDerived):
    def __init__(self, coord, is_plus_strand, feature_type, phase=0, score=None, source=None, controller=None):
        """Initializes a handler for a soon to be inserted geenuff feature.
        As there is no 1to1 relation between a gff entry and a geenuff feature,
        no gff entry is saved or directly used to infer attributes from.

        In order to inherit the insertion functions from GFFDerived, this class
        still inherits from it.
        """
        super().__init__()
        GFFDerived.__init__(self)
        if controller is not None:
            self.id = InsertCounterHolder.feature()
        self.coord = coord
        self.is_plus_strand = is_plus_strand
        self.type = feature_type
        self.phase = phase
        self.score = score
        self.source = source
        # these will have to be adapted to geenuff
        self.start = -1
        self.end = -1
        self.start_is_biological_start = None
        self.end_is_biological_end = None

    def setup_insertion_ready(self,
                              super_locus=None,
                              transcribed_pieces=None,
                              translateds=None,
                              coordinate=None):

        if transcribed_pieces is None:
            transcribed_pieces = []
        if translateds is None:
            translateds = []

        parents = self.gffentry.get_Parent()
        # which can occur e.g. when we just copied gffentry from gene for an error
        if parents is not None:
            for piece in transcribed_pieces:
                assert piece.gffentry.get_ID() in parents

        if self._is_plus_strand is None or self._given_name is None:
            self.add_shortcuts_from_gffentry()

        feature = {
            'id': self.id,
            'subtype': 'general',
            'given_name': self.given_name,
            'coordinate_id': coordinate.id,
            'is_plus_strand': self.is_plus_strand,
            'score': self.gffentry.score,
            'source': self.gffentry.source,
            'phase': self.gffentry.phase,  # todo, do I need to handle '.'?
            'super_locus_id': super_locus.id,
            # all currently set via kwargs in add_to_queue, as they aren't clear from  gffentry
            'start': None,
            'end': None,  # note the dict must always have all keys for core
            'start_is_biological_start': None,
            'end_is_biological_end': None
        }

        feature2pieces = [{
            'transcribed_piece_id': p.id,
            'feature_id': self.id
        } for p in transcribed_pieces]
        feature2translateds = [{'translated_id': p.id, 'feature_id': self.id} for p in translateds]

        return feature, feature2pieces, feature2translateds

    def add_to_queue(self,
                     insertion_queues,
                     super_locus=None,
                     transcribed_pieces=None,
                     translateds=None,
                     coordinate=None,
                     **kwargs):

        feature, feature2pieces, feature2translateds = self.setup_insertion_ready(
            super_locus,
            transcribed_pieces,  # todo, add flexibility for parsing trans-splicing...
            translateds,
            coordinate)

        for key in kwargs:
            feature[key] = kwargs[key]

        insertion_queues.feature.queue.append(feature)
        insertion_queues.association_transcribed_piece_to_feature.queue += feature2pieces
        insertion_queues.association_translated_to_feature.queue += feature2translateds
        return feature, feature2pieces, feature2translateds

    # "+ strand" [upstream, downstream) or "- strand" (downstream, upstream] from 0 coordinate
    def upstream_from_interval(self, interval):
        if self.is_plus_strand:
            return interval.begin
        else:
            # -1 bc as this is now a start, it should be _inclusive_ (and flipped)
            return interval.end - 1

    def downstream_from_interval(self, interval):
        if self.is_plus_strand:
            return interval.end
        else:
            return interval.begin - 1  # -1 to be _exclusive_ (and flipped)

    def upstream(self):
        if self.is_plus_strand:
            return self.gffentry.start
        else:
            return self.gffentry.end

    def downstream(self):
        if self.is_plus_strand:
            return self.gffentry.end
        else:
            return self.gffentry.start

    @property
    def py_start(self):
        return helpers.as_py_start(self.gffentry.start)

    @property
    def py_end(self):
        return helpers.as_py_end(self.gffentry.end)


class TranscribedHandler(handlers.TranscribedHandlerBase, GFFDerived):
    def __init__(self, controller=None):
        super().__init__()
        GFFDerived.__init__(self)
        if controller is not None:
            self.id = InsertCounterHolder.transcribed()
        self.controller = controller
        self.transcribed_piece_handlers = []
        self.feature_handlers = []

    def setup_insertion_ready(self, gffentry=None, super_locus=None):
        if gffentry is not None:
            parents = gffentry.get_Parent()
            # the simple case
            if len(parents) == 1:
                assert super_locus.given_name == parents[0]
                entry_type = gffentry.type
                given_name = gffentry.get_ID()
            else:
                raise NotImplementedError  # todo handle multi inheritance, etc...
        else:
            entry_type = given_name = None

        transcribed_2_add = {
            'type': entry_type,
            'given_name': given_name,
            'super_locus_id': super_locus.id,
            'id': self.id
        }

        piece = TranscribedPieceHandler(controller=self.controller)
        piece.gffentry = copy.deepcopy(gffentry)

        piece_2_add = piece.setup_insertion_ready(gffentry,
                                                  super_locus=super_locus,
                                                  transcribed=self,
                                                  position=0)
        self.transcribed_piece_handlers.append(piece)
        return transcribed_2_add, piece_2_add

    def add_to_queue(self, insertion_queues, gffentry=None, super_locus=None, **kwargs):
        transcribed_add, piece_add = self.setup_insertion_ready(gffentry=gffentry,
                                                                super_locus=super_locus)
        for key in kwargs:
            transcribed_add[key] = kwargs[key]

        insertion_queues.transcribed.queue.append(transcribed_add)
        insertion_queues.transcribed_piece.queue.append(piece_add)
        return transcribed_add, piece_add

    def one_piece(self):
        pieces = self.transcribed_piece_handlers
        assert len(pieces) == 1
        return pieces[0]


class TranscribedPieceHandler(handlers.TranscribedPieceHandlerBase, GFFDerived):
    def __init__(self, controller=None):
        super().__init__()
        GFFDerived.__init__(self)
        if controller is not None:
            self.id = InsertCounterHolder.transcribed_piece()

    def setup_insertion_ready(self, gffentry=None, super_locus=None, transcribed=None, position=0):
        if gffentry is not None:
            parents = gffentry.get_Parent()
            # the simple case
            if len(parents) == 1:
                assert super_locus.given_name == parents[0]
                given_name = gffentry.get_ID()
            else:
                raise NotImplementedError  # todo handle multi inheritance, etc...
        else:
            given_name = None

        transcribed_piece = {
            'id': self.id,
            'transcribed_id': transcribed.id,
            'super_locus_id': super_locus.id,
            'given_name': given_name,
            'position': position,
        }
        return transcribed_piece


class TranslatedHandler(handlers.TranslatedHandlerBase):
    def __init__(self, controller=None):  # todo, all handlers, controller=None
        super().__init__()
        if controller is not None:
            self.id = InsertCounterHolder.translated()

    def setup_insertion_ready(self, super_locus_handler=None, given_name=''):
        translated = {
            'id': self.id,
            'super_locus_id': super_locus_handler.id,
            'given_name': given_name,
        }
        return translated


class TranscriptInterpreter(TranscriptInterpBase):
    """takes raw/from-gff transcript, and converts to geenuff explicit format"""
    HANDLED = 'handled'

    def __init__(self, transcript, super_locus, controller):
        super().__init__(transcript, super_locus, session=controller.session)
        self.status = EukTranscriptStatus()
        assert isinstance(controller, ImportController)
        self.controller = controller
        try:
            self.proteins = self._setup_proteins()
        except NoGFFEntryError:
            # this way we only run into an error if we actually wanted to use proteins
            self.proteins = None

    def new_feature(self, gffentry, translateds=None, **kwargs):
        coordinate = self.controller.latest_genome_handler.gffid_to_coords[gffentry.seqid]
        handler = FeatureHandler(controller=self.controller, processed=True)
        handler.gffentry = copy.deepcopy(gffentry)
        feature, _, _ = handler.add_to_queue(
            insertion_queues=self.controller.insertion_queues,
            super_locus=self.super_locus,
            transcribed_pieces=self.transcript.transcribed_piece_handlers,
            translateds=translateds,
            coordinate=coordinate,
            **kwargs)

        return feature

    @staticmethod
    def pick_one_interval(interval_set, target_type=None):
        if target_type is None:
            return interval_set[0]
        else:
            return [x for x in interval_set if x.data.gffentry.type == target_type][0]

    @staticmethod
    def _get_protein_id_from_cds(cds_feature):
        try:
            assert cds_feature.gffentry.type == types.CDS, "{} != {}".format(
                cds_feature.gff_entry.type, types.CDS)
        except AttributeError:
            raise NoGFFEntryError('No gffentry for {}'.format(cds_feature.data.given_name))
        # check if anything is labeled as protein_id
        protein_id = cds_feature.gffentry.attrib_filter(tag='protein_id')
        # failing that, try and get parent ID (presumably transcript, maybe gene)
        if not protein_id:
            protein_id = cds_feature.gffentry.get_Parent()
        # hopefully take single hit
        if len(protein_id) == 1:
            protein_id = protein_id[0]
            if isinstance(protein_id, gffhelper.GFFAttribute):
                protein_id = protein_id.value
                assert len(protein_id) == 1
                protein_id = protein_id[0]
        # or handle other cases
        elif len(protein_id) == 0:
            protein_id = None
        else:
            raise ValueError('indeterminate single protein id {}'.format(protein_id))
        return protein_id

    def _get_raw_protein_ids(self):
        # only meant for use before feature interpretation
        protein_ids = set()
        for feature in self.transcript.feature_handlers:
            try:
                if feature.gffentry.type == types.CDS:
                    protein_id = self._get_protein_id_from_cds(feature)
                    protein_ids.add(protein_id)
            except AttributeError as e:
                print('failing here.........')
                print(len(self.transcript.feature_handlers))
                print(self.transcript.feature_handlers)
                raise e
        return protein_ids

    def _setup_proteins(self):
        # only meant for use before feature interpretation
        pids = self._get_raw_protein_ids()
        proteins = {}
        for key in pids:
            # setup blank protein
            protein = TranslatedHandler(controller=self.controller)
            # ready translated for database entry
            translated = protein.setup_insertion_ready(super_locus_handler=self.super_locus,
                                                       given_name=key)
            self.controller.insertion_queues.translated.queue.append(translated)
            proteins[key] = protein
        return proteins

    def is_plus_strand(self):
        features = self.transcript.feature_handlers
        seqids = [x.gffentry.seqid for x in features]
        if not all([x == seqids[0] for x in seqids]):
            raise TransSplicingError("non matching seqids {}, for {}".format(
                seqids, self.super_locus.id))
        if all([x.is_plus_strand for x in features]):
            return True
        elif all([not x.is_plus_strand for x in features]):
            return False
        else:
            raise TransSplicingError("Mixed strands at {} with {}".format(
                self.super_locus.id, [(x.coordinate.seqid, x.strand) for x in features]))

    def drop_invervals_with_duplicated_data(self, ivals_before, ivals_after):
        all_data = [x.data.data for x in ivals_before + ivals_after]
        if len(all_data) > len(set(all_data)):

            marked_before = []
            for ival in ivals_before:
                new = {
                    'interval': ival,
                    'matches_edge': self.matches_edge_before(ival),
                    'repeated': self.is_in_list(ival.data, [x.data for x in ivals_after])
                }
                marked_before.append(new)
            marked_after = []
            for ival in ivals_after:
                new = {
                    'interval': ival,
                    'matches_edge': self.matches_edge_after(ival),
                    'repeated': self.is_in_list(ival.data, [x.data for x in ivals_before])
                }
                marked_after.append(new)

            # warn if slice-causer (matches_edge) is on both sides
            # should be fine in cases with [exon, CDS] -> [exon, UTR] or similar
            matches_before = [x['interval'] for x in marked_before if x['matches_edge']]
            slice_frm_before = bool(matches_before)
            matches_after = [x['interval'] for x in marked_after if x['matches_edge']]
            slice_frm_after = bool(matches_after)
            if slice_frm_before and slice_frm_after:
                # if it's not a splice site
                if matches_before[0].end == matches_after[0].begin:
                    logging.warning(
                        'slice causer on both sides with repeates\nbefore: {}, after: {}'.format(
                            [(x.data.gffentry.type, x.data.gffentry.start, x.data.gffentry.end)
                             for x in ivals_before],
                            [(x.data.gffentry.type, x.data.gffentry.start, x.data.gffentry.end)
                             for x in ivals_after]))
            # finally, keep non repeats or where this side didn't cause slice
            ivals_before = [
                x['interval'] for x in marked_before if not x['repeated'] or not slice_frm_before
            ]
            ivals_after = [
                x['interval'] for x in marked_after if not x['repeated'] or not slice_frm_after
            ]
        return ivals_before, ivals_after

    @staticmethod
    def is_in_list(target, the_list):
        matches = [x for x in the_list if x is target]
        return len(matches) > 0

    @staticmethod
    def matches_edge_before(ival):
        data = ival.data
        if data.is_plus_strand:
            out = data.py_end == ival.end
        else:
            out = data.py_start == ival.begin
        return out

    @staticmethod
    def matches_edge_after(ival):
        data = ival.data
        if data.is_plus_strand:
            out = data.py_start == ival.begin
        else:
            out = data.py_end == ival.end
        return out

    def interpret_transition(self, ivals_before, ivals_after, plus_strand=True):
        sign = 1
        if not plus_strand:
            sign = -1
        ivals_before, ivals_after = self.drop_invervals_with_duplicated_data(
            ivals_before, ivals_after)
        before_types = self.possible_types(ivals_before)
        after_types = self.possible_types(ivals_after)
        # 5' UTR can hit either start codon or splice site
        if self.status.is_5p_utr():
            # start codon
            self.handle_from_5p_utr(ivals_before, ivals_after, before_types, after_types, sign)
        elif self.status.is_coding():
            self.handle_from_coding(ivals_before, ivals_after, before_types, after_types, sign)
        elif self.status.is_3p_utr():
            self.handle_from_3p_utr(ivals_before, ivals_after, before_types, after_types, sign)
        elif self.status.is_intronic():
            self.handle_from_intron()
        elif self.status.is_intergenic():
            self.handle_from_intergenic()
        else:
            raise ValueError('unknown status {}'.format(self.status.__dict__))

    def handle_from_coding(self, ivals_before, ivals_after, before_types, after_types, sign):
        assert types.CDS in before_types
        # stop codon
        if types.THREE_PRIME_UTR in after_types:
            self.handle_control_codon(ivals_before, ivals_after, sign, is_start=False)
        # splice site
        elif types.CDS in after_types:
            self.handle_splice(ivals_before, ivals_after, sign)
        else:
            raise ValueError("don't know how to transition from coding to {}".format(after_types))

    def handle_from_intron(self):
        raise NotImplementedError  # todo later

    def handle_from_3p_utr(self, ivals_before, ivals_after, before_types, after_types, sign):
        assert types.THREE_PRIME_UTR in before_types
        # the only thing we should encounter is a splice site
        if types.THREE_PRIME_UTR in after_types:
            self.handle_splice(ivals_before, ivals_after, sign)
        else:
            raise ValueError('wrong feature types after three prime: b: {}, a: {}'.format(
                [x.gffentry.type for x in ivals_before], [x.gffentry.type for x in ivals_after]))

    def handle_from_5p_utr(self, ivals_before, ivals_after, before_types, after_types, sign):
        assert types.FIVE_PRIME_UTR in before_types
        # start codon
        if types.CDS in after_types:
            self.handle_control_codon(ivals_before, ivals_after, sign, is_start=True)
        # intron
        elif types.FIVE_PRIME_UTR in after_types:
            self.handle_splice(ivals_before, ivals_after, sign)
        else:
            raise ValueError('wrong feature types after five prime: b: {}, a: {}'.format(
                [x.gffentry.type for x in ivals_before], [x.gffentry.type for x in ivals_after]))

    def handle_from_intergenic(self):
        raise NotImplementedError  # todo later

    def is_gap(self, ivals_before, ivals_after, sign):
        """checks for a gap between intervals, and validates it's a positive
        one on strand of interest
        """
        after0 = self.pick_one_interval(ivals_after)
        before0 = self.pick_one_interval(ivals_before)
        before_downstream = before0.data.downstream_from_interval(before0)
        after_upstream = after0.data.upstream_from_interval(after0)
        is_gap = before_downstream != after_upstream
        if is_gap:
            # if there's a gap, confirm it's in the right direction
            gap_len = (after_upstream - before_downstream) * sign
            assert gap_len > 0, "inverse gap between {} and {} at putative control codon seq {}, gene {}, " \
                                "features {} {}".format(before_downstream, after_upstream,
                                                        after0.data.coordinate.seqid, self.super_locus.id,
                                                        before0.data.id, after0.data.id)
        return is_gap

    def get_protein(self, feature):
        pid = self._get_protein_id_from_cds(feature)
        return self.proteins[pid]

    def handle_control_codon(self,
                             ivals_before,
                             ivals_after,
                             sign,
                             is_start=True,
                             error_buffer=2000):
        target_after_type = None
        target_before_type = None
        if is_start:
            target_after_type = types.CDS
        else:
            target_before_type = types.CDS

        after0 = self.pick_one_interval(ivals_after, target_after_type)
        before0 = self.pick_one_interval(ivals_before, target_before_type)
        # make sure there is no gap
        is_gap = self.is_gap(ivals_before, ivals_after, sign)

        if is_start:
            if is_gap:
                self.handle_splice(ivals_before, ivals_after, sign)

            template = after0.data
            translated = self.get_protein(template)
            # it better be std phase if it's a start codon
            at = template.upstream_from_interval(after0)
            if template.gffentry.phase == 0:  # "non-0 phase @ {} in {}".format(template.id, template.super_locus.id)
                start = at
                coding_feature = self.new_feature(gffentry=template.gffentry,
                                                  start=start,
                                                  type=types.CODING,
                                                  start_is_biological_start=True,
                                                  translateds=[translated])
                self.status.coding_tracker.set_channel_open_with_feature(feature=coding_feature,
                                                                         phase=0)
            else:
                err_start = before0.data.upstream_from_interval(
                    before0) - sign * error_buffer  # mask prev feat. too
                err_end = at + sign  # so that the error masks the coordinate with the close status
                gffid_to_coords = self.controller.latest_genome_handler.gffid_to_coords
                if sign == 1:  # if is_plus_strand, todo, use same logic/setup as elsewhere
                    start_of_sequence = gffid_to_coords[template.gffentry.seqid].start
                    err_start = max(start_of_sequence, err_start)
                else:
                    end_of_sequence = gffid_to_coords[template.gffentry.seqid].end - 1
                    err_start = min(end_of_sequence, err_start)
                self.new_feature(gffentry=template.gffentry,
                                 type=types.ERROR,
                                 start=err_start,
                                 end=err_end,
                                 start_is_biological_start=True,
                                 end_is_biological_end=True,
                                 phase=None)

                coding_feature = self.new_feature(gffentry=template.gffentry,
                                                  type=types.CODING,
                                                  start=at,
                                                  start_is_biological_start=False,
                                                  translateds=[translated])
                self.status.coding_tracker.set_channel_open_with_feature(
                    feature=coding_feature, phase=template.gffentry.phase)

        else:
            template = before0.data
            at = template.downstream_from_interval(before0)
            self.status.coding_tracker.update_and_close_feature({
                'end': at,
                'end_is_biological_end': True
            })

            if is_gap:
                self.handle_splice(ivals_before, ivals_after, sign)

    def handle_splice(self, ivals_before, ivals_after, sign):
        target_type = None
        if self.status.is_coding():
            target_type = types.CDS

        before0 = self.pick_one_interval(ivals_before, target_type)
        after0 = self.pick_one_interval(ivals_after, target_type)
        donor_tmplt = before0.data
        acceptor_tmplt = after0.data
        if sign > 0:
            # use original coords so that overlaps can be detected...
            donor_at = donor_tmplt.downstream(
            )  # -1 (to from 0) + 1 (intron start, not incl. exon end) = 0
            acceptor_at = acceptor_tmplt.upstream(
            ) - 1  # -1 (to from 0), that's it as incl start is already excl end
        else:
            donor_at = donor_tmplt.downstream(
            ) - 2  # -1 (fr 0), -1 (intron start (-), not incl. exon end) = -2
            acceptor_at = acceptor_tmplt.upstream(
            ) - 1  # -1 (fr 0), that's it as incl start is already excl end
        # add splice sites if there's a gap
        between_splice_sites = (acceptor_at - donor_at) * sign
        min_intron_len = 3  # todo, maybe get something small but not entirely impossible?
        if between_splice_sites > min_intron_len:
            self.new_feature(gffentry=donor_tmplt.gffentry,
                             start=donor_at,
                             end=acceptor_at,
                             phase=None,
                             start_is_biological_start=True,
                             end_is_biological_end=True,
                             type=types.INTRON)

        # do nothing if there is just no gap between exons for a technical / reporting error
        elif between_splice_sites == 0:
            pass
        # everything else is invalid
        else:
            # mask both exons and the intron
            if sign > 0:
                err_start = donor_tmplt.upstream() - 1  # -1 to fr 0
                err_end = acceptor_tmplt.downstream()  # -1 to fr 0, +1 to excl = 0
            else:
                err_start = donor_tmplt.upstream() - 1  # to fr 0
                err_end = acceptor_tmplt.downstream() - 2  # -1 to fr 0, -1 to excl = -2

            self.new_feature(gffentry=before0.data.gffentry,
                             start=err_start,
                             end=err_end,
                             start_is_biological_start=True,
                             end_is_biological_end=True,
                             type=types.ERROR,
                             bearing=types.START)

    def interpret_first_pos(self, intervals, plus_strand=True, error_buffer=2000):
        # shortcuts
        cds = types.CDS

        i0 = self.pick_one_interval(intervals)
        at = i0.data.upstream_from_interval(i0)
        possible_types = self.possible_types(intervals)
        if types.FIVE_PRIME_UTR in possible_types:
            # this should indicate we're good to go and have a transcription start site
            transcribed_feature = self.new_feature(gffentry=i0.data.gffentry,
                                                   type=types.TRANSCRIBED,
                                                   start=at,
                                                   phase=None,
                                                   start_is_biological_start=True)
            self.status.transcribed_tracker.set_channel_open_with_feature(
                feature=transcribed_feature)
        elif cds in possible_types:
            # this could be first exon detected or start codon, ultimately, indeterminate
            cds_feature = self.pick_one_interval(intervals, target_type=cds).data
            translated = self.get_protein(cds_feature)
            # setup translated / coding feature
            translated_feature = self.new_feature(gffentry=cds_feature.gffentry,
                                                  type=types.CODING,
                                                  start=at,
                                                  start_is_biological_start=False,
                                                  translateds=[translated])
            self.status.coding_tracker.set_channel_open_with_feature(
                feature=translated_feature, phase=cds_feature.gffentry.phase)
            # setup transcribed feature (bc coding implies transcript)
            transcribed_feature = self.new_feature(gffentry=cds_feature.gffentry,
                                                   type=types.TRANSCRIBED,
                                                   start=at,
                                                   start_is_biological_start=False)
            self.status.transcribed_tracker.set_channel_open_with_feature(transcribed_feature)

            # mask a dummy region up-stream as it's very unclear whether it should be intergenic/intronic/utr
            gffid_to_coords = self.controller.latest_genome_handler.gffid_to_coords
            if plus_strand:
                # unless we're at the start of the sequence
                start_of_sequence = gffid_to_coords[cds_feature.gffentry.seqid].start
                if at != start_of_sequence:
                    err_start = max(start_of_sequence, at - error_buffer)
                    err_end = at + 1  # so that the error masks the coordinate with the close status
                    self.new_feature(gffentry=cds_feature.gffentry,
                                     type=types.ERROR,
                                     start_is_biological_start=True,
                                     end_is_biological_end=True,
                                     start=err_start,
                                     end=err_end,
                                     phase=None)
            else:
                end_of_sequence = gffid_to_coords[cds_feature.gffentry.seqid].end - 1
                if at != end_of_sequence:
                    err_start = min(end_of_sequence, at + error_buffer)
                    err_end = at - 1  # so that the error masks the coordinate with the close status
                    self.new_feature(gffentry=cds_feature.gffentry,
                                     type=types.ERROR,
                                     start_is_biological_start=True,
                                     end_is_biological_end=True,
                                     start=err_start,
                                     end=err_end,
                                     phase=None)

        else:
            raise ValueError(
                "why's this gene not start with 5' utr nor cds? types: {}, interpretations: {}".
                format([x.data.gffentry.type for x in intervals], possible_types))

    def interpret_last_pos(self, intervals, plus_strand=True, error_buffer=2000):
        i0 = self.pick_one_interval(intervals)
        at = i0.data.downstream_from_interval(i0)
        possible_types = self.possible_types(intervals)
        if types.THREE_PRIME_UTR in possible_types:
            # this should be transcription termination site
            self.status.transcribed_tracker.update_and_close_feature({
                'end': at,
                'end_is_biological_end': True
            })
        elif types.CDS in possible_types:
            # may or may not be stop codon, but will just mark as error (unless at edge of sequence)
            cds_feature = self.pick_one_interval(intervals, target_type=types.CDS).data

            self.status.transcribed_tracker.update_and_close_feature({
                'end': at,
                'end_is_biological_end': False
            })
            self.status.coding_tracker.update_and_close_feature({
                'end': at,
                'end_is_biological_end': False
            })
            # may or may not be stop codon, but will just mark as error (unless at edge of sequence)
            gffid_to_coords = self.controller.latest_genome_handler.gffid_to_coords

            start_of_sequence = gffid_to_coords[cds_feature.gffentry.seqid].start - 1
            end_of_sequence = gffid_to_coords[cds_feature.gffentry.seqid].end
            if plus_strand:
                if at != end_of_sequence:
                    err_start = at - 1  # so that the error masks the coordinate with the open status
                    err_end = min(at + error_buffer, end_of_sequence)
                    self.new_feature(gffentry=i0.data.gffentry,
                                     type=types.ERROR,
                                     start=err_start,
                                     end=err_end,
                                     start_is_biological_start=True,
                                     end_is_biological_end=True,
                                     phase=None)
            else:
                if at != start_of_sequence:
                    err_start = at + 1  # so that the error masks the coordinate with the open status
                    err_end = max(start_of_sequence, at - error_buffer)
                    self.new_feature(gffentry=i0.data.gffentry,
                                     type=types.ERROR,
                                     phase=None,
                                     start=err_start,
                                     start_is_biological_start=True,
                                     end=err_end,
                                     end_is_biological_end=True)
        else:
            raise ValueError(
                "why's this gene not end with 3' utr/exon nor cds? types: {}, interpretations: {}".
                format([x.data.type for x in intervals], possible_types))

    def intervals_5to3(self, plus_strand=False):
        interval_sets = list(self.organize_and_split_features())
        if not plus_strand:
            interval_sets.reverse()
        return interval_sets

    def decode_raw_features(self):
        plus_strand = self.is_plus_strand()
        interval_sets = self.intervals_5to3(plus_strand)
        self.interpret_first_pos(interval_sets[0], plus_strand)
        for i in range(len(interval_sets) - 1):
            ivals_before = interval_sets[i]
            ivals_after = interval_sets[i + 1]
            self.interpret_transition(ivals_before, ivals_after, plus_strand)

        self.interpret_last_pos(intervals=interval_sets[-1], plus_strand=plus_strand)

    @staticmethod
    def possible_types(intervals):
        # shortcuts
        cds = types.CDS
        five_prime = types.FIVE_PRIME_UTR
        exon = types.EXON
        three_prime = types.THREE_PRIME_UTR

        # what we see
        # uniq_datas = set([x.data.data for x in intervals])  # todo, revert and skip unique once handled above
        observed_types = [x.data.gffentry.type for x in intervals]
        set_o_types = set(observed_types)
        # check length
        if len(intervals) not in [1, 2]:
            raise IntervalCountError(
                'check interpretation by hand for transcript start with {}, {}'.format(
                    '\n'.join([str(ival.data.data) for ival in intervals]), observed_types))
        if set_o_types.issubset(set([x.value for x in types.KeepOnSequence])):
            out = [TranscriptInterpreter.HANDLED]
        # interpret type combination
        elif set_o_types == {exon, five_prime} or set_o_types == {five_prime}:
            out = [five_prime]
        elif set_o_types == {exon, three_prime} or set_o_types == {three_prime}:
            out = [three_prime]
        elif set_o_types == {exon}:
            out = [five_prime, three_prime]
        elif set_o_types == {cds, exon} or set_o_types == {cds}:
            out = [cds]
        else:
            raise ValueError(
                'check interpretation of combination for transcript start with {}, {}'.format(
                    intervals, observed_types))
        return out

    def organize_and_split_features(self):
        # todo, handle non-single seqid loci
        tree = intervaltree.IntervalTree()
        features = set()
        for feature in self._all_features():
            features.add(feature)

        for f in features:
            tree[helpers.as_py_start(f.gffentry.start):helpers.as_py_end(f.gffentry.end)] = f
        tree.split_overlaps()
        # todo, minus strand
        intervals = iter(sorted(tree))
        out = [next(intervals)]
        for interval in intervals:
            if out[-1].begin == interval.begin:
                out.append(interval)
            else:
                yield out
                out = [interval]
        yield out

    def _all_features(self):
        return self.transcript.feature_handlers

    def has_no_unprocessed_features(self):
        for feature in self._all_features():
            if not feature.processed:
                return False
        return True

    @staticmethod
    def _get_intervals_by_type(stacked_intervals, target_type):
        out = []
        for stack in stacked_intervals:
            new_stack = []
            for interval in stack:
                if interval.data.gffentry.type == target_type:
                    new_stack.append(interval)
            if new_stack:
                out.append(new_stack)
        return out

    def _by_pos_and_type(self, stacked_intervals, pos, target_type):
        exons = self._get_intervals_by_type(stacked_intervals, target_type=target_type)
        pos_exons = exons[pos]
        assert len(pos_exons) == 1
        return pos_exons[0]

    def first_exon(self, stacked_intervals):
        return self._by_pos_and_type(stacked_intervals, 0, target_type=types.EXON)

    def last_exon(self, stacked_intervals):
        return self._by_pos_and_type(stacked_intervals, -1, target_type=types.EXON)

    def first_cds(self, stacked_intervals):
        return self._by_pos_and_type(stacked_intervals, 0, target_type=types.CDS)

    def last_cds(self, stacked_intervals):
        return self._by_pos_and_type(stacked_intervals, -1, target_type=types.CDS)
