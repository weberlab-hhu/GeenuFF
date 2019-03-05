from dustdas import gffhelper, fastahelper
import intervaltree
from .. import api
from .. import orm
from .. import types
from .. import helpers

import copy
import logging
import hashlib

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


##### main flow control #####
class ImportControl(object):

    def __init__(self, database_path, err_path=None):
        if not database_path.startswith('sqlite:///'):
            database_path = 'sqlite:///{}'.format(database_path)
        self.database_path = database_path
        self.session = None
        self.err_path = err_path
        self.engine = None
        self.annotated_genome = None
        self.sequence_info = None
        self.coordinates = {}
        self.super_loci = []
        self.mk_session()
        # queues for adding to db
        self.features_to_add = []
        self.feature2transcribed_pieces_to_add = []
        self.feature2translateds_to_add = []
        self.translated2transcribeds_to_add = []
        self.translateds_to_add = []
        self.transcribeds_to_add = []
        self.transcribed_pieces_to_add = []
        self.super_loci_to_add = []

        # counting everything that goes in
        self.feature_counter = helpers.Counter()
        self.translated_counter = helpers.Counter()
        self.transcribed_counter = helpers.Counter()
        self.super_locus_counter = helpers.Counter()
        self.transcribed_piece_counter = helpers.Counter()
        self.sequence_info_counter = helpers.Counter()
        self.coordinate_counter = helpers.Counter()

    def mk_session(self):
        self.engine = create_engine(self.database_path, echo=False)  # todo, dynamic / real path
        orm.Base.metadata.create_all(self.engine)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    def gff_gen(self, gff_file):
        known = [x.value for x in types.AllKnown]
        reader = gffhelper.read_gff_file(gff_file)
        for entry in reader:
            if entry.type not in known:
                raise ValueError("unrecognized feature type from gff: {}".format(entry.type))
            else:
                self.clean_entry(entry)
                yield entry

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

    def make_anno_genome(self, **kwargs):
        # todo, parse in meta data from kwargs?
        self.annotated_genome = api.AnnotatedGenomeHandler()
        ag = orm.AnnotatedGenome()
        self.annotated_genome.add_data(ag)
        self.session.add(ag)
        self.session.commit()

    def add_sequences(self, seq_path):
        if self.annotated_genome is None:
            self.make_anno_genome()

        self._setup_sequence_info()
        self.sequence_info.add_sequences(seq_path)

    def _setup_sequence_info(self):
        self.sequence_info = SequenceInfoHandler(controller=self)
        seq_info = orm.SequenceInfo(annotated_genome=self.annotated_genome.data)
        self.sequence_info.add_data(seq_info)

    def add_gff(self, gff_file, clean=True):
        # final prepping of seqid match up
        self.sequence_info.mk_mapper(gff_file)

        super_loci = []
        err_handle = open(self.err_path, 'w')
        i = 1
        for entry_group in self.group_gff_by_gene(gff_file):
            super_locus = SuperLocusHandler(controller=self)
            if self.sequence_info is None:
                raise AttributeError(
                    ' sequence_info cannot be None when .add_gff is called, use (self.add_sequences(...)')
            if clean:
                super_locus.add_n_clean_gff_entry_group(entry_group, err_handle, sequence_info=self.sequence_info.data,
                                                        session=self.session, controller=self)
            else:
                super_locus.add_gff_entry_group(entry_group, err_handle, sequence_info=self.sequence_info.data,
                                                controller=self)
            self.session.add(super_locus.data)
            self.session.commit()
            super_loci.append(super_locus)  # just to keep some direct python link to this
            if not i % 500:
                self.execute_so_far()
            i += 1
        self.super_loci = super_loci
        self.execute_so_far()
        err_handle.close()

    def execute_so_far(self):
        insert_lists = [(orm.SuperLocus.__table__.insert(), self.super_loci_to_add),
                        (orm.Transcribed.__table__.insert(), self.transcribeds_to_add),
                        (orm.TranscribedPiece.__table__.insert(), self.transcribed_pieces_to_add),
                        (orm.Translated.__table__.insert(), self.translateds_to_add),
                        (orm.association_translateds_to_transcribeds.insert(), self.translated2transcribeds_to_add),
                        (orm.Feature.__table__.insert(), self.features_to_add),
                        (orm.association_transcribeds_to_features.insert(), self.feature2transcribed_pieces_to_add),
                        (orm.association_translateds_to_features.insert(), self.feature2translateds_to_add)
                        ]
        conn = self.engine.connect()
        for insert, a_list in insert_lists:
            if a_list:
                conn.execute(insert, a_list)
            del a_list[:]

    def clean_super_loci(self):
        for sl in self.super_loci:
            coordinates = self.sequence_info.gffid_to_coords[sl.gffentry.seqid]
            sl.check_and_fix_structure(self.session, coordinates, controller=self)


def in_values(x, enum):
    return x in [item.value for item in enum]


class SequenceInfoHandler(api.SequenceInfoHandler):
    def __init__(self, controller=None):
        super().__init__()
        if controller is not None:
            self.id = controller.sequence_info_counter()
        self.mapper = None
        self._seq_info = None
        self._gffid_to_coords = None
        self._gff_seq_ids = None

    def mk_mapper(self, gff_file=None):
        fa_ids = [x.seqid for x in self.data.coordinates]
        if gff_file is not None:  # allow setup without ado when we know IDs match exactly
            self._gff_seq_ids = helpers.get_seqids_from_gff(gff_file)
        else:
            self._gff_seq_ids = fa_ids
        mapper, is_forward = helpers.two_way_key_match(fa_ids, self._gff_seq_ids)
        self.mapper = mapper

        if not is_forward:
            raise NotImplementedError("Still need to implement backward match if fasta IDs are subset of gff IDs")

    @property
    def gffid_to_coords(self):
        if self._gffid_to_coords is not None:
            pass
        else:
            gffid2coords = {}
            for gffid in self._gff_seq_ids:
                fa_id = self.mapper(gffid)
                x = self.seq_info[fa_id]
                gffid2coords[gffid] = x
            self._gffid_to_coords = gffid2coords
        return self._gffid_to_coords

    @property
    def seq_info(self):
        if self._seq_info is not None:
            pass
        else:
            seq_info = {}
            for x in self.data.coordinates:
                seq_info[x.seqid] = x
            self._seq_info = seq_info
        return self._seq_info

    def add_sequences(self, seq_file):
        self.add_fasta(seq_file)

    def add_fasta(self, seq_file, id_delim=' '):
        fp = fastahelper.FastaParser()
        for fasta_header, seq in fp.read_fasta(seq_file):
            seqid = fasta_header.split(id_delim)[0]
            # todo, parallelize sequence & annotation format, then import directly from sequence_info (~Slice)
            orm.Coordinates(seqid=seqid, start=0, sha1=self.hashseq(seq),
                            end=len(seq), sequence_info=self.data)

    @staticmethod
    def hashseq(seq):
        m = hashlib.sha1()
        m.update(seq.encode())
        return m.hexdigest()


##### gff parsing subclasses #####
class GFFDerived(object):
    def __init__(self):
        self.gffentry = None

    def process_gffentry(self, gffentry, gen_data=True, **kwargs):
        self.gffentry = gffentry
        data = None
        if gen_data:
            data = self.gen_data_from_gffentry(gffentry, **kwargs)
        return data

    def gen_data_from_gffentry(self, gffentry, **kwargs):
        # should create 'data' object (annotations_orm.Base subclass) and then call self.add_data(data)
        raise NotImplementedError


class SuperLocusHandler(api.SuperLocusHandler, GFFDerived):
    def __init__(self, controller=None):
        api.SuperLocusHandler.__init__(self)
        GFFDerived.__init__(self)
        if controller is not None:
            self.id = controller.super_locus_counter()
        self.transcribed_handlers = []

    def gen_data_from_gffentry(self, gffentry, **kwargs):
        data = self.data_type(type=gffentry.type,
                              given_id=gffentry.get_ID(),
                              id=self.id)
        self.add_data(data)
        # todo, grab more aliases from gff attribute

    def add_gff_entry(self, entry, sequence_info, controller):
        exceptions = entry.attrib_filter(tag="exception")
        for exception in [x.value for x in exceptions]:
            if 'trans-splicing' in exception:
                raise TransSplicingError('trans-splice in attribute {} {}'.format(entry.get_ID(), entry.attribute))
        if in_values(entry.type, types.SuperLocusAll):
            self.process_gffentry(gffentry=entry)

        elif in_values(entry.type, types.TranscriptLevelAll):
            transcribed = TranscribedHandler(controller)
            transcribed.gffentry = copy.deepcopy(entry)

            transcribed_add, piece_add = transcribed.setup_insertion_ready(gffentry=entry, super_locus=self.data)
            controller.transcribeds_to_add.append(transcribed_add)
            controller.transcribed_pieces_to_add.append(piece_add)
            self.transcribed_handlers.append(transcribed)

        elif in_values(entry.type, types.OnSequence):
            feature = FeatureHandler(controller)
            feature.gffentry = copy.deepcopy(entry)
            assert len(self.transcribed_handlers) > 0, "no transcribeds found before feature"
            # MOD_READIN, will need to set up features with temporary linkage, but without entering them into final db
            feature.add_shortcuts_from_gffentry()
            self.transcribed_handlers[-1].feature_handlers.append(feature)
        else:
            raise ValueError("problem handling entry of type {}".format(entry.type))

    def _add_gff_entry_group(self, entries, sequence_info, controller):
        entries = list(entries)
        for entry in entries:
            self.add_gff_entry(entry, sequence_info, controller)

    def add_gff_entry_group(self, entries, ts_err_handle, sequence_info, controller):
        try:
            self._add_gff_entry_group(entries, sequence_info, controller)
            #self.check_and_fix_structure(entries)
        except TransSplicingError as e:
            coordinates = sequence_info.handler.gffid_to_coords[entries[0].seqid]
            self._mark_erroneous(entries[0], coordinates, controller, 'trans-splicing')
            logging.warning('skipping but noting trans-splicing: {}'.format(str(e)))
            ts_err_handle.writelines([x.to_json() for x in entries])
        except AssertionError as e:
            coordinates = sequence_info.handler.gffid_to_coords[entries[0].seqid]
            self._mark_erroneous(entries[0], coordinates, controller, str(e))
            logging.warning('marking errer bc {}'.format(e))

    def add_n_clean_gff_entry_group(self, entries, ts_err_handle, sequence_info, session, controller):
        self.add_gff_entry_group(entries, ts_err_handle, sequence_info, controller)
        session.commit()
        coordinates = sequence_info.handler.gffid_to_coords[self.gffentry.seqid]
        self.check_and_fix_structure(session, coordinates, controller=controller)

    def _mark_erroneous(self, entry, coordinates, controller, msg=''):
        assert entry.type in [x.value for x in types.SuperLocusAll]
        logging.warning(
            '{species}:{seqid}, {start}-{end}:{gene_id} by {src}, {msg} - marked erroneous'.format(
                src=entry.source, species="todo", seqid=entry.seqid, start=entry.start,
                end=entry.end, gene_id=self.data.given_id, msg=msg
            ))
        # reset start and stop so
        if entry.strand == '+':
            err_start = entry.start
            err_end = entry.end
        else:
            assert entry.strand == '-'
            err_start = entry.end
            err_end = entry.start
        # dummy transcript
        # todo, CLEAN UP / get in functions
        transcribed_e_handler = TranscribedHandler(controller)
        transcribed_e_handler.gffentry = copy.deepcopy(entry)
        transcribed, piece = transcribed_e_handler.setup_insertion_ready(super_locus=self)
        piece_handler = transcribed_e_handler.transcribed_piece_handlers[0]
        controller.transcribeds_to_add.append(transcribed)
        controller.transcribed_pieces_to_add.append(piece)
        # open error
        feature_err_open = FeatureHandler(controller, processed=True)
        feature_err_open.gffentry = copy.deepcopy(entry)
        packed = feature_err_open.setup_insertion_ready(super_locus=self.data, transcribed_pieces=[piece_handler],
                                                        coordinates=coordinates)
        feature, feature2pieces, feature2translateds = packed
        for key, val in [('type', types.ERROR), ('bearing', types.START), ('position', err_start)]:
            feature[key] = val
        controller.features_to_add.append(feature)
        controller.feature2transcribed_pieces_to_add += feature2pieces
        controller.feature2translateds_to_add += feature2translateds

        # close error
        feature_err_close = FeatureHandler(controller, processed=True)
        feature_err_close.gffentry = copy.deepcopy(entry)
        packed = feature_err_close.setup_insertion_ready(super_locus=self.data,transcribed_pieces=[piece_handler],
                                                         coordinates=coordinates)
        feature, feature2pieces, feature2translateds = packed
        for key, val in [('type', types.ERROR), ('bearing', types.END), ('position', err_end)]:
            feature[key] = val
        controller.features_to_add.append(feature)
        controller.feature2transcribed_pieces_to_add += feature2pieces
        controller.feature2translateds_to_add += feature2translateds

        transcribed_e_handler.feature_handlers += [feature_err_open, feature_err_close]
        self.transcribed_handlers.append(transcribed_e_handler)

    def check_and_fix_structure(self, sess, coordinates, controller):
        # todo, add against sequence check to see if start/stop and splice sites are possible or not, e.g. is start ATG?
        # if it's empty (no bottom level features at all) mark as erroneous
        features = []
        for transcript in self.transcribed_handlers:
            features += transcript.feature_handlers
        if not features:
            self._mark_erroneous(self.gffentry, coordinates=coordinates, controller=controller)
            sess.commit()

        for transcript in self.transcribed_handlers:
            t_interpreter = TranscriptInterpreter(transcript, super_locus=self, controller=controller)
            # skip any transcript consisting of only processed features (in context, should just be pre-interp errors)
            if t_interpreter.has_no_unprocessed_features():
                pass
            else:
                # make new features
                t_interpreter.decode_raw_features()
                # make sure the new features link to protein if appropriate
                sess.commit()
                controller.execute_so_far()
                #t_interpreter.mv_coding_features_to_proteins(controller.feature2translateds_to_add)
            #forced_keep_handlers += t_interpreter.clean_features


class FeatureHandler(api.FeatureHandler, GFFDerived):

    def __init__(self, controller=None, processed=False):
        api.FeatureHandler.__init__(self)
        GFFDerived.__init__(self)
        if controller is not None:
            self.id = controller.feature_counter()
        self.processed = processed
        self._is_plus_strand = None
        self._given_id = None

    @property
    def is_plus_strand(self):
        if self.data is not None:
            return self.data.is_plus_strand
        elif self._is_plus_strand is not None:
            return self._is_plus_strand
        else:
            raise ValueError('attempt to use ._is_plus_strand while it is still None')

    @property
    def given_id(self):
        if self.data is not None:
            return self.data.given_id
        elif self._given_id is not None:
            return self._given_id
        else:
            raise ValueError('attempt to use ._given_id while it is still None')

    def add_shortcuts_from_gffentry(self):

        self._given_id = self.gffentry.get_ID()
        if self.gffentry.strand == '+':
            self._is_plus_strand = True
        elif self.gffentry.strand == '-':
            self._is_plus_strand = False
        else:
            raise ValueError('cannot interpret strand "{}"'.format(self.gffentry.strand))

    def setup_insertion_ready(self, super_locus=None, transcribed_pieces=None, translateds=None,
                              coordinates=None):

        if transcribed_pieces is None:
            transcribed_pieces = []
        if translateds is None:
            translateds = []

        parents = self.gffentry.get_Parent()
        if parents is not None:  # which can occur e.g. when we just copied gffentry from gene for an error
            for piece in transcribed_pieces:
                assert piece.gffentry.get_ID() in parents

        if self._is_plus_strand is None or self._given_id is None:
            self.add_shortcuts_from_gffentry()

        feature = {'id': self.id,
                   'given_id': self.given_id,
                   'coordinate_id': coordinates.id,
                   'is_plus_strand': self.is_plus_strand,
                   'score': self.gffentry.score,
                   'source': self.gffentry.source,
                   'phase': self.gffentry.phase,  # todo, do I need to handle '.'?
                   'super_locus_id': super_locus.id,
                   }

        feature2pieces = [{'transcribed_piece_id': p.id, 'feature_id': self.id} for p in transcribed_pieces]
        feature2translateds = [{'translated_id': p.id, 'feature_id': self.id} for p in translateds]

        return feature, feature2pieces, feature2translateds

    # "+ strand" [upstream, downstream) or "- strand" (downstream, upstream] from 0 coordinates
    def upstream_from_interval(self, interval):
        if self.is_plus_strand:
            return interval.begin
        else:
            return interval.end - 1  # -1 bc as this is now a start, it should be _inclusive_ (and flipped)

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


class TranscribedHandler(api.TranscribedHandler, GFFDerived):
    def __init__(self, controller=None):
        api.TranscribedHandler.__init__(self)
        GFFDerived.__init__(self)
        if controller is not None:
            self.id = controller.transcribed_counter()
        self.controlller = controller
        self.transcribed_piece_handlers = []
        self.feature_handlers = []
        self.translated_handlers = []

    def setup_insertion_ready(self, gffentry=None, super_locus=None, **kwargs):
        if gffentry is not None:
            parents = gffentry.get_Parent()
            # the simple case
            if len(parents) == 1:
                assert super_locus.given_id == parents[0]
                entry_type = gffentry.type
                given_id = gffentry.get_ID()
            else:
                raise NotImplementedError  # todo handle multi inheritance, etc...
        else:
            entry_type = given_id = None

        transcribed_2_add = {'type': entry_type, 'given_id': given_id, 'super_locus_id': super_locus.id,
                             'id': self.id}

        piece = TranscribedPieceHandler(controller=self.controlller)
        piece.gffentry = copy.deepcopy(gffentry)

        piece_2_add = piece.setup_insertion_ready(gffentry, super_locus=super_locus, transcribed=self)
        self.transcribed_piece_handlers.append(piece)
        return transcribed_2_add, piece_2_add

    def one_piece(self):
        pieces = self.transcribed_piece_handlers
        assert len(pieces) == 1
        return pieces[0]


class TranscribedPieceHandler(api.TranscribedPieceHandler, GFFDerived):
    def __init__(self, controller=None):
        api.TranscribedPieceHandler.__init__(self)
        GFFDerived.__init__(self)
        if controller is not None:
            self.id = controller.transcribed_piece_counter()

    def setup_insertion_ready(self, gffentry=None, super_locus=None, transcribed=None, **kwargs):
        if gffentry is not None:
            parents = gffentry.get_Parent()
            # the simple case
            if len(parents) == 1:
                assert super_locus.given_id == parents[0]
                given_id = gffentry.get_ID()
            else:
                raise NotImplementedError  # todo handle multi inheritance, etc...
        else:
            given_id = None

        transcribed_piece = {'id': self.id, 'transcribed_id': transcribed.id, 'super_locus_id': super_locus.id,
                             'given_id': given_id}
        return transcribed_piece


class TranslatedHandler(api.TranslatedHandler):
    def __init__(self, controller=None):  # todo, all handlers, controller=None
        super().__init__()
        if controller is not None:
            self.id = controller.translated_counter()

    def setup_insertion_ready(self, transcribed_handler, super_locus_handler, given_id):
        translated = {'id': self.id, 'super_locus_id': super_locus_handler.id, 'given_id': given_id}
        translated_to_transcribed = {'translated_id': self.id, 'transcribed_id': transcribed_handler.id}
        return translated, translated_to_transcribed


class NoTranscriptError(Exception):
    pass


class TransSplicingError(Exception):
    pass


class NoGFFEntryError(Exception):
    pass


class TranscriptInterpreter(api.TranscriptInterpBase):
    """takes raw/from-gff transcript, and makes totally explicit"""
    HANDLED = 'handled'

    def __init__(self, transcript, super_locus, controller):
        super().__init__(transcript, super_locus, session=controller.session)
        self.controller = controller
        try:
            self.proteins = self._setup_proteins()
        except NoGFFEntryError:
            self.proteins = None  # this way we only run into an error if we actually wanted to use proteins

    def new_feature(self, template, translateds=None, **kwargs):
        coordinates = self.controller.sequence_info.gffid_to_coords[template.gffentry.seqid]
        handler = FeatureHandler(controller=self.controller, processed=True)
        handler.gffentry = copy.deepcopy(template.gffentry)

        feature, feature2pieces, feature2translateds = handler.setup_insertion_ready(
            self.super_locus,
            self.transcript.transcribed_piece_handlers,  # todo, add flexibility for parsing trans-splicing...
            translateds,
            coordinates)
        for key in kwargs:
            feature[key] = kwargs[key]

        self.controller.features_to_add.append(feature)
        self.controller.feature2transcribed_pieces_to_add += feature2pieces
        self.controller.feature2translateds_to_add += feature2translateds

    @staticmethod
    def pick_one_interval(interval_set, target_type=None):
        if target_type is None:
            return interval_set[0]
        else:
            return [x for x in interval_set if x.data.gffentry.type == target_type][0]

    @staticmethod
    def _get_protein_id_from_cds(cds_feature):
        try:
            assert cds_feature.gffentry.type == types.CDS, "{} != {}".format(cds_feature.gff_entry.type,
                                                                             types.CDS)
        except AttributeError:
            raise NoGFFEntryError('No gffentry for {}'.format(cds_feature.data.given_id))
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
            # ready translated and link to transcribed for database entry
            translated, translated_to_transcribed = protein.setup_insertion_ready(
                transcribed_handler=self.transcript,
                super_locus_handler=self.super_locus,
                given_id=key
            )
            self.controller.translated2transcribeds_to_add.append(translated_to_transcribed)
            self.controller.translateds_to_add.append(translated)
            proteins[key] = protein

        return proteins

    def is_plus_strand(self):
        features = self.transcript.feature_handlers
        seqids = [x.gffentry.seqid for x in features]
        if not all([x == seqids[0] for x in seqids]):
            raise TransSplicingError("non matching seqids {}, for {}".format(seqids, self.super_locus.id))
        if all([x.is_plus_strand for x in features]):
            return True
        elif all([not x.is_plus_strand for x in features]):
            return False
        else:
            raise TransSplicingError(
                "Mixed strands at {} with {}".format(self.super_locus.id,
                                                     [(x.coordinates.seqid, x.strand) for x in features]))

    def drop_invervals_with_duplicated_data(self, ivals_before, ivals_after):
        all_data = [x.data.data for x in ivals_before + ivals_after]
        if len(all_data) > len(set(all_data)):

            marked_before = []
            for ival in ivals_before:
                new = {'interval': ival, 'matches_edge': self.matches_edge_before(ival),
                       'repeated': self.is_in_list(ival.data, [x.data for x in ivals_after])}
                marked_before.append(new)
            marked_after = []
            for ival in ivals_after:
                new = {'interval': ival, 'matches_edge': self.matches_edge_after(ival),
                       'repeated': self.is_in_list(ival.data, [x.data for x in ivals_before])}
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
                    logging.warning('slice causer on both sides with repeates\nbefore: {}, after: {}'.format(
                        [(x.data.gffentry.type, x.data.gffentry.start, x.data.gffentry.end) for x in ivals_before],
                        [(x.data.gffentry.type, x.data.gffentry.start, x.data.gffentry.end) for x in ivals_after]
                    ))
            # finally, keep non repeats or where this side didn't cause slice
            ivals_before = [x['interval'] for x in marked_before if not x['repeated'] or not slice_frm_before]
            ivals_after = [x['interval'] for x in marked_after if not x['repeated'] or not slice_frm_after]
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
        ivals_before, ivals_after = self.drop_invervals_with_duplicated_data(ivals_before, ivals_after)
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
        """checks for a gap between intervals, and validates it's a positive one on strand of interest"""
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
                                                        after0.data.coordinates.seqid, self.super_locus.id,
                                                        before0.data.id, after0.data.id)
        return is_gap

    def get_protein(self, feature):
        pid = self._get_protein_id_from_cds(feature)
        return self.proteins[pid]

    def handle_control_codon(self, ivals_before, ivals_after, sign, is_start=True, error_buffer=2000):
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
                self.new_feature(template=template, position=start, type=types.CODING,
                                 bearing=types.START, translateds=[translated])
                self.status.saw_start(phase=0)
            else:
                err_start = before0.data.upstream_from_interval(before0) - sign * error_buffer  # mask prev feat. too
                err_end = at + sign  # so that the error masks the coordinate with the close status

                self.new_feature(template=template, type=types.ERROR, position=err_start,
                                 phase=None, bearing=types.START)
                self.new_feature(template=template, type=types.ERROR, position=err_end,
                                 phase=None, bearing=types.END)
                self.new_feature(template=template, type=types.CODING, position=at,
                                 bearing=types.OPEN_STATUS, translateds=[translated])
                self.new_feature(template=template, type=types.TRANSCRIBED, position=at,
                                 phase=None, bearing=types.OPEN_STATUS)
                self.status.saw_start(template.phase)
        else:
            # todo, confirm phase for stop codon
            template = before0.data
            translated = self.get_protein(template)
            at = template.downstream_from_interval(before0)
            self.new_feature(template=template, position=at, type=types.CODING,
                             bearing=types.END, translateds=[translated])
            self.status.saw_stop()
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
            donor_at = donor_tmplt.downstream()  # -1 (to from 0) + 1 (intron start, not incl. exon end) = 0
            acceptor_at = acceptor_tmplt.upstream() - 1  # -1 (to from 0), that's it as incl start is already excl end
        else:
            donor_at = donor_tmplt.downstream() - 2  # -1 (fr 0), -1 (intron start (-), not incl. exon end) = -2
            acceptor_at = acceptor_tmplt.upstream() - 1  # -1 (fr 0), that's it as incl start is already excl end
        # add splice sites if there's a gap
        between_splice_sites = (acceptor_at - donor_at) * sign
        min_intron_len = 3  # todo, maybe get something small but not entirely impossible?
        if between_splice_sites > min_intron_len:
            self.new_feature(template=donor_tmplt, position=donor_at, phase=None,
                             type=types.INTRON, bearing=types.START)
            self.new_feature(template=acceptor_tmplt, position=acceptor_at,
                             type=types.INTRON, bearing=types.END)
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

            self.new_feature(template=before0.data, position=err_start,
                             type=types.ERROR, bearing=types.START)
            self.new_feature(template=before0.data, position=err_end,
                             type=types.ERROR, bearing=types.END)

    def interpret_first_pos(self, intervals, plus_strand=True, error_buffer=2000):
        # shortcuts
        cds = types.CDS

        i0 = self.pick_one_interval(intervals)
        at = i0.data.upstream_from_interval(i0)
        possible_types = self.possible_types(intervals)
        if types.FIVE_PRIME_UTR in possible_types:
            # this should indicate we're good to go and have a transcription start site
            self.new_feature(template=i0.data, type=types.TRANSCRIBED, position=at,
                             phase=None, bearing=types.START)
            self.status.saw_tss()
        elif cds in possible_types:
            # this could be first exon detected or start codon, ultimately, indeterminate
            cds_feature = self.pick_one_interval(intervals, target_type=cds).data
            translated = self.get_protein(cds_feature)
            self.new_feature(template=cds_feature, type=types.CODING, position=at,
                             bearing=types.OPEN_STATUS, translateds=[translated])
            self.new_feature(template=cds_feature, type=types.TRANSCRIBED, position=at,
                             bearing=types.OPEN_STATUS)
            self.status.saw_start(phase=cds_feature.gffentry.phase)
            self.status.saw_tss()  # coding implies the transcript
            # mask a dummy region up-stream as it's very unclear whether it should be intergenic/intronic/utr
            if plus_strand:
                # unless we're at the start of the sequence
                start_of_sequence = self.controller.sequence_info.gffid_to_coords[cds_feature.gffentry.seqid].start
                if at != start_of_sequence:
                    err_start = max(start_of_sequence, at - error_buffer)
                    err_end = at + 1  # so that the error masks the coordinate with the close status
                    self.new_feature(template=cds_feature, type=types.ERROR,
                                     bearing=types.START,
                                     position=err_start, phase=None)
                    self.new_feature(template=cds_feature, type=types.ERROR,
                                     bearing=types.END,
                                     position=err_end, phase=None)
            else:
                end_of_sequence = self.controller.sequence_info.gffid_to_coords[cds_feature.gffentry.seqid].end - 1
                #end_of_sequence = cds_feature.data.coordinates.end - 1   # bc we need last valid index for coordinates
                if at != end_of_sequence:
                    err_start = min(end_of_sequence, at + error_buffer)
                    err_end = at - 1  # so that the error masks the coordinate with the close status
                    self.new_feature(template=cds_feature, type=types.ERROR,
                                     bearing=types.START,
                                     position=err_start, phase=None)
                    self.new_feature(template=cds_feature, type=types.ERROR,
                                     bearing=types.END,
                                     position=err_end, phase=None)
        else:
            raise ValueError("why's this gene not start with 5' utr nor cds? types: {}, interpretations: {}".format(
                [x.data.gffentry.type for x in intervals], possible_types))

    def interpret_last_pos(self, intervals, plus_strand=True, error_buffer=2000):
        i0 = self.pick_one_interval(intervals)
        at = i0.data.downstream_from_interval(i0)
        possible_types = self.possible_types(intervals)
        if types.THREE_PRIME_UTR in possible_types:
            # this should be transcription termination site
            self.new_feature(template=i0.data, type=types.TRANSCRIBED, bearing=types.END, position=at,
                             phase=None)
            self.status.saw_tts()
        elif types.CDS in possible_types:
            # may or may not be stop codon, but will just mark as error (unless at edge of sequence)
            cds_feature = self.pick_one_interval(intervals, target_type=types.CDS).data
            translated = self.get_protein(cds_feature)
            self.new_feature(template=cds_feature, type=types.CODING, position=at,
                             bearing=types.CLOSE_STATUS, translateds=[translated])
            self.new_feature(template=cds_feature, type=types.TRANSCRIBED, position=at,
                             bearing=types.CLOSE_STATUS)
            self.status.saw_start(phase=cds_feature.gffentry.phase)
            self.status.saw_tss()  # coding implies the transcript
            # may or may not be stop codon, but will just mark as error (unless at edge of sequence)
            start_of_sequence = self.controller.sequence_info.gffid_to_coords[cds_feature.gffentry.seqid].start - 1
            #start_of_sequence = i0.data.data.coordinates.start - 1  # because we need to be able to
            end_of_sequence = self.controller.sequence_info.gffid_to_coords[cds_feature.gffentry.seqid].end
            #end_of_sequence = i0.data.data.coordinates.end
            if plus_strand:
                if at != end_of_sequence:
                    err_start = at - 1  # so that the error masks the coordinate with the open status
                    err_end = min(at + error_buffer, end_of_sequence)
                    self.new_feature(template=i0.data, type=types.ERROR, position=err_start,
                                     phase=None, bearing=types.START)

                    self.new_feature(template=i0.data, type=types.ERROR, position=err_end,
                                     phase=None, bearing=types.END)
            else:
                if at != start_of_sequence:
                    err_start = at + 1  # so that the error masks the coordinate with the open status
                    err_end = max(start_of_sequence, at - error_buffer)
                    self.new_feature(template=i0.data, type=types.ERROR,
                                     phase=None, position=err_start, bearing=types.START)
                    self.new_feature(template=i0.data, type=types.ERROR,
                                     phase=None, position=err_end, bearing=types.END)
        else:
            raise ValueError("why's this gene not end with 3' utr/exon nor cds? types: {}, interpretations: {}".format(
                [x.data.type for x in intervals], possible_types)
            )

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

        self.interpret_last_pos(intervals=interval_sets[-1])

    def possible_types(self, intervals):
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
            raise IntervalCountError('check interpretation by hand for transcript start with {}, {}'.format(
                '\n'.join([str(ival.data.data) for ival in intervals]), observed_types
            ))
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
            raise ValueError('check interpretation of combination for transcript start with {}, {}'.format(
                intervals, observed_types
            ))
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
        # note,
        out = True
        for feature in self._all_features():
            if not feature.processed:
                out = False
        return out


class IntervalCountError(Exception):
    pass


def none_to_list(x):
    if x is None:
        return []
    else:
        assert isinstance(x, list)
        return x


def upstream(x, y, sign):
    if (y - x) * sign >= 0:
        return x
    else:
        return y


def downstream(x, y, sign):
    if (x - y) * sign >= 0:
        return x
    else:
        return y
