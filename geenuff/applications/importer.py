import os
from pprint import pprint  # for debugging
from abc import ABC, abstractmethod
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from dustdas import gffhelper, fastahelper
from .. import orm
from .. import types
from .. import handlers
from .. import helpers
from ..base.helpers import (get_strand_direction, get_geenuff_start_end,
                            has_start_codon, has_stop_codon)


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


class OrganizedGeenuffHandlerGroup(object):
    """Stores the handler objects for a super locus in an organized fassion.
    The format is similar to the one of OrganizedGFFEntryGroup, but it stores objects
    according to the Geenuff way of saving genomic annotations. This format can then
    be checked for errors and changed accordingly before being inserted into the db.

    The handlers are organized in the following way:

    handlers = {
        'super_locus_h' = super_locus_handler,
        'transcript_hs': [
            {
                'transcript_h': transcript_handler,
                'transcript_piece_h: transcript_piece_handler,
                'transcript_feature_h': transcript_handler,
                'protein_h': protein_handler,
                'cds_h': cds_handler,
                'intron_hs': [intron_handler1, intron_handler2, ..]
            },
            ...
        ],
        'error_hs' = [error_handler1, error_handler2, ..],
    }
    """

    def __init__(self, organized_gff_entries, coord, controller, err_handle):
        self.coord = coord
        self.controller = controller
        self.err_handle = err_handle
        self.handlers = {'transcript_hs': [], 'error_hs': []}
        self._parse_gff_entries(organized_gff_entries)

    def _parse_gff_entries(self, entries):
        """Changes the GFF format into the GeenuFF format. Does all the parsing."""
        sl = entries['super_locus']
        # these attr are the same for all geenuff handlers going to be created
        is_plus_strand = get_strand_direction(sl)
        score = sl.score
        source = sl.source

        sl_start, sl_end = get_geenuff_start_end(sl.start, sl.end, is_plus_strand)
        sl_h = self.handlers['super_locus_h'] = SuperLocusHandler(entry_type=sl.type,
                                                                  given_name=sl.get_ID(),
                                                                  coord=self.coord,
                                                                  is_plus_strand=is_plus_strand,
                                                                  start=sl_start,
                                                                  end=sl_end,
                                                                  controller=self.controller)
        for t, t_entries in entries['transcripts'].items():
            t_handlers = {}
            # check for multi inheritance and throw NotImplementedError if found
            if len(t.get_Parent()) > 1:
                raise NotImplementedError
            t_id = t.get_ID()
            # create transcript handler
            t_h = TranscribedHandler(entry_type=t.type,
                                     given_name=t_id,
                                     super_locus_id=sl_h.id,
                                     controller=self.controller)
            # create transcript piece handler
            tp_h = TranscribedPieceHandler(given_name=t_id,
                                           transcript_id=t_h.id,
                                           position=0,
                                           controller=self.controller)
            # create transcript feature handler
            tf_h = FeatureHandler(self.coord,
                                  is_plus_strand,
                                  types.TRANSCRIBED,
                                  given_name=t_id,
                                  score=score,
                                  source=source,
                                  controller=self.controller)
            tf_h.set_start_end_from_gff(t.start, t.end)

            # insert everything so far into the dict
            t_handlers['transcript_h'] = t_h
            t_handlers['transcript_piece_h'] = tp_h
            t_handlers['transcript_feature_h'] = tf_h

            # if it is not a non-coding gene or something like that
            if t_entries['cds']:
                # create protein handler
                protein_id = self._get_protein_id_from_cds_list(t_entries['cds'])
                p_h = TranslatedHandler(given_name=protein_id,
                                        super_locus_id=sl_h.id,
                                        controller=self.controller)
                # create coding features from exon limits
                cds_h = FeatureHandler(self.coord,
                                       is_plus_strand,
                                       types.CODING,
                                       # take the coding phase from the first cds entry in the gff
                                       phase=t_entries['cds'][0].phase,
                                       score=score,
                                       source=source,
                                       controller=self.controller)
                gff_start = t_entries['cds'][0].start
                gff_end = t_entries['cds'][-1].end
                cds_h.set_start_end_from_gff(gff_start, gff_end)

                # insert everything so far into the dict
                t_handlers['protein_h'] = p_h
                t_handlers['cds_h'] = cds_h

                # create all the introns by traversing the exon entries and insert FeatureHandlers
                # into the previously created list
                # the introns should strictly always lie between successive gff exons
                exons = t_entries['exons']
                intron_hs = []
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
                    intron_hs.append(intron_h)
                t_handlers['intron_hs'] = intron_hs
            self.handlers['transcript_hs'].append(t_handlers)

    @staticmethod
    def _get_protein_id_from_cds_entry(cds_entry):
        # check if anything is labeled as protein_id
        protein_id = cds_entry.attrib_filter(tag='protein_id')
        # failing that, try and get parent ID (presumably transcript, maybe gene)
        if not protein_id:
            protein_id = cds_entry.get_Parent()
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

    @staticmethod
    def _get_protein_id_from_cds_list(cds_entry_list):
        """Returns the protein id of a list of cds gff entries. If multiple ids or no id at all
        are found, an error is raised."""
        protein_ids = set()
        for cds_entry in cds_entry_list:
            protein_id = OrganizedGeenuffHandlerGroup._get_protein_id_from_cds_entry(cds_entry)
            if protein_id:
                protein_ids.add(protein_id)
        if len(protein_ids) != 1:
            raise ValueError('No protein_id or more than one protein_ids for one transcript')
        return protein_ids.pop()


class OrganizedGFFEntryGroup(object):
    """Takes an entry group (all entries of one super locus) and stores the entries
    in an orderly fassion. Can then return a corresponding OrganizedGeenuffHandlerGroup.
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

        # order exon and cds lists by start value (disregard strand for now)
        for _, value_dict in self.entries['transcripts'].items():
            for key in ['exons', 'cds']:
                value_dict[key].sort(key=lambda e:e.start)

    def get_geenuff_handlers(self):
        geenuff_handler_group = OrganizedGeenuffHandlerGroup(self.entries,
                                                             self.coord,
                                                             self.controller,
                                                             self.err_handle)
        return geenuff_handler_group.handlers


class OrganizedGFFEntries(object):
    """Structures the gff entries comming from gffhelper by seqid and gene. Also does some
    basic gff value cleanup.
    The entries are organized in the following way:

    organized_entries = {
        'seqid1': [
            [gff_entry1_gene1, gff_entry2_gene1, ..],
            [gff_entry1_gene2, gff_entry2_gene2, ..],
        ],
        'seqid2': [
            [gff_entry1_gene1, gff_entry2_gene1, ..],
            [gff_entry1_gene2, gff_entry2_gene2, ..],
        ],
        ...
    }
    """
    def __init__(self, gff_file):
        self.gff_file = gff_file
        self.organized_entries = self._organize_entries()

    def _organize_entries(self):
        organized = {}
        gene_level = [x.value for x in types.SuperLocusAll]

        reader = self._useful_gff_entries()
        first = next(reader)
        seqid = first.seqid
        gene_group = [first]
        organized[seqid] = []
        for entry in reader:
            if entry.type in gene_level:
                organized[seqid].append(gene_group)
                gene_group = [entry]
                if entry.seqid != seqid:
                    organized[entry.seqid] = []
                    seqid = entry.seqid
            else:
                gene_group.append(entry)
        organized[seqid].append(gene_group)
        return organized

    def _useful_gff_entries(self):
        skipable = [x.value for x in types.IgnorableFeatures]
        reader = self._gff_gen()
        for entry in reader:
            if entry.type not in skipable:
                yield entry

    def _gff_gen(self):
        known = [x.value for x in types.AllKnown]
        reader = gffhelper.read_gff_file(self.gff_file)
        for entry in reader:
            if entry.type not in known:
                raise ValueError("unrecognized feature type from gff: {}".format(entry.type))
            else:
                self._clean_entry(entry)
                yield entry

    @staticmethod
    def _clean_entry(entry):
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


class GFFErrorHandling(object):
    """Deal with error detection and handling of the input features. Does the handling
    in the space of GeenuFF handlers.
    Assumes all super locus handler groups to be ordered 5p to 3p and of one strand.
    Works with a list of OrganizedGeenuffHandlerGroup, which correspond to a list of
    super loci, and looks for errors. Error features may be inserted and handlers be
    removed when deemed necessary.
    """
    def __init__(self, geenuff_handler_groups, controller):
        self.groups = geenuff_handler_groups
        if self.groups:
            self.is_plus_strand = self.groups[0]['super_locus_h'].is_plus_strand
        self.controller = controller

    def resolve_errors(self):
        for i, group in enumerate(self.groups):
            # the case of no transcript for a super locus
            # solution is to add an error mask that extends halfway to the appending super
            # loci in the intergenic region as something in this area appears to have gone wrong
            error_hs = group['error_hs']
            if not group['transcript_hs']:
                error_hs.append(self._get_overlapping_err_h(i, None, 'whole',
                                                            types.EMPTY_SUPER_LOCUS))

            # other cases
            for t_hs in group['transcript_hs']:
                # if coding transcript
                if 'cds_h' in t_hs:
                    # the case of missing of implicit UTR ranges
                    # the solution is similar to the one above
                    cds = t_hs['cds_h']
                    if cds.start == t_hs['transcript_feature_h'].start:
                        start = cds.start
                        error_hs.append(self._get_overlapping_err_h(i, start, '5p',
                                                                    types.MISSING_UTR_5P))
                    if cds.end == t_hs['transcript_feature_h'].end:
                        start = cds.end
                        error_hs.append(self._get_overlapping_err_h(i, start, '3p',
                                                                    types.MISSING_UTR_3P))
                    # the case of missing start/stop codon
                    if not has_start_codon(cds.coord.sequence, cds.start, self.is_plus_strand):
                        start = cds.start
                        error_hs.append(self._get_overlapping_err_h(i, start, '5p',
                                                                    types.MISSING_START_CODON))
                    if not has_stop_codon(cds.coord.sequence, cds.end, self.is_plus_strand):
                        start = cds.end
                        error_hs.append(self._get_overlapping_err_h(i, start, '3p',
                                                                    types.MISSING_STOP_CODON))

    def _get_error_handler(self, coord, start, end, is_plus_strand, error_type):
        # todo logging
        error_h = FeatureHandler(coord,
                                 is_plus_strand,
                                 error_type,
                                 start=start,
                                 end=end,
                                 controller=self.controller)
        return error_h

    def _get_overlapping_err_h(self, i, start, direction, error_type):
        """Constructs an error features that overlaps halfway to the next super locus
        in the given direction if possible. Otherwise mark until the end.
        The error feature extends in the given direction from the start value on.
        If the direction is 'whole', the start value is ignored.
        """
        assert direction in ['5p', '3p', 'whole']
        sl_h = self.groups[i]['super_locus_h']

        # set correct upstream error starting point
        if direction in ['5p', 'whole']:
            if i > 0:
                sl_h_prev = self.groups[i - 1]['super_locus_h']
                anker_5p = self._halfway_mark(sl_h_prev, sl_h)
            else:
                if self.is_plus_strand:
                    anker_5p = 0
                else:
                    anker_5p = sl_h.start

        # set correct downstream error end point
        if direction in ['3p', 'whole']:
            if i < len(self.groups) -1:
                sl_h_next = self.groups[i + 1]['super_locus_h']
                anker_3p = self._halfway_mark(sl_h, sl_h_next)
            else:
                if self.is_plus_strand:
                    anker_3p = sl_h.end
                else:
                    anker_3p = -1

        if direction == '5p':
            error_start = anker_5p
            error_end = start
        elif direction =='3p':
            error_start = start
            error_end = anker_3p
        elif direction == 'whole':
            error_start = anker_5p
            error_end = anker_3p
        error_h = self._get_error_handler(sl_h.coord,
                                          error_start,
                                          error_end,
                                          self.is_plus_strand,
                                          error_type)
        return error_h

    def _halfway_mark(self, sl_h, sl_h_next):
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
            in the correct order. Also initiates the insert of the many2many rows.
            """
            for group in groups:
                group['super_locus_h'].add_to_queue()
                # insert all features as well as transcript and protein related entries
                for transcript_hs in group['transcript_hs']:
                    # make shortcuts
                    tp_h = transcript_hs['transcript_piece_h']
                    tf_h = transcript_hs['transcript_feature_h']
                    # add transcript handler that are always present
                    transcript_hs['transcript_h'].add_to_queue()
                    tp_h.add_to_queue()
                    tf_h.add_to_queue()
                    tf_h.insert_feature_piece_association(tp_h.id)
                    # if coding transcript
                    if 'protein_h' in transcript_hs:
                        transcript_hs['protein_h'].add_to_queue()
                        tf_h.insert_feature_protein_association(transcript_hs['protein_h'].id)
                        transcript_hs['cds_h'].add_to_queue()
                        transcript_hs['cds_h'].insert_feature_piece_association(tp_h.id)
                        for intron_h in transcript_hs['intron_hs']:
                            intron_h.add_to_queue()
                            intron_h.insert_feature_piece_association(tp_h.id)
                # insert the errors
                for e_h in group['error_hs']:
                    e_h.add_to_queue()

        def clean_and_insert(self, groups, clean):
            plus = [g for g in groups if g['super_locus_h'].is_plus_strand]
            minus = [g for g in groups if not g['super_locus_h'].is_plus_strand]
            if clean:
                # check and correct for errors
                # do so for each strand seperately
                # all changes should be made by reference
                GFFErrorHandling(plus, self).resolve_errors()
                # reverse order on minus strand
                GFFErrorHandling(minus[::-1], self).resolve_errors()
            # insert handlers
            insert_handler_groups(self, plus)
            insert_handler_groups(self, minus)
            self.insertion_queues.execute_so_far()

        assert self.latest_genome_handler is not None, 'No recent genome found'
        self.latest_genome_handler.mk_mapper(gff_file)
        organized_gff_entries = OrganizedGFFEntries(gff_file).organized_entries
        geenuff_handler_groups = []
        for seqid in organized_gff_entries.keys():
            for i, entry_group in enumerate(organized_gff_entries[seqid]):
                organized_entries = OrganizedGFFEntryGroup(entry_group,
                                                           self.latest_genome_handler,
                                                           self,
                                                           err_handle)
                geenuff_handler_groups.append(organized_entries.get_geenuff_handlers())
                if i > 0 and i % 500 == 0:
                    clean_and_insert(self, geenuff_handler_groups, clean)
                    geenuff_handler_groups = []
            # never do error checking across fasta sequence borders
            clean_and_insert(self, geenuff_handler_groups, clean)
            geenuff_handler_groups = []


class Insertable(ABC):
    @abstractmethod
    def add_to_queue(self):
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
    def __init__(self, entry_type, given_name, controller, coord=None, is_plus_strand=None, start=-1, end=-1):
        super().__init__()
        self.id = InsertCounterHolder.super_locus()
        self.entry_type = entry_type
        self.given_name = given_name
        self.controller = controller
        # not neccessary for insert but helpful for certain error cases
        self.coord = coord
        self.is_plus_strand = is_plus_strand
        self.start = start
        self.end = end

    def add_to_queue(self):
        to_add = {'type': self.entry_type, 'given_name': self.given_name, 'id': self.id}
        self.controller.insertion_queues.super_locus.queue.append(to_add)

    def __repr__(self):
        params = {'id': self.id, 'type': self.entry_type, 'given_name': self.given_name}
        return self._get_repr('SuperLocusHandler', params)


class FeatureHandler(handlers.FeatureHandlerBase, Insertable):
    def __init__(self, coord, is_plus_strand, feature_type, controller, start=-1, end=-1, given_name=None, phase=0, score=None, source=None):
        """Initializes a handler for a soon to be inserted geenuff feature."""
        super().__init__()
        self.id = InsertCounterHolder.feature()
        self.coord = coord
        self.given_name = given_name
        self.is_plus_strand = is_plus_strand
        self.feature_type = feature_type
        # start/end may have to be adapted to geenuff
        self.start = start
        self.end = end
        self.phase = phase
        self.score = score
        self.source = source
        self.start_is_biological_start = None
        self.end_is_biological_end = None
        self.controller = controller

    def add_to_queue(self):
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
        self.controller.insertion_queues.feature.queue.append(feature)

    def insert_feature_piece_association(self, transcript_piece_id):
        features2pieces = {
            'feature_id': self.id,
            'transcribed_piece_id': transcript_piece_id,
        }
        self.controller.insertion_queues.association_transcribed_piece_to_feature.\
            queue.append(features2pieces)

    def insert_feature_protein_association(self, protein_id):
        features2protein = {
            'feature_id': self.id,
            'translated_id': protein_id,
        }
        self.controller.insertion_queues.association_translated_to_feature.\
            queue.append(features2protein)

    def set_start_end_from_gff(self, gff_start, gff_end):
        self.start, self.end = get_geenuff_start_end(gff_start, gff_end, self.is_plus_strand)

    def pos_cmp_key(self):
        sortable_start = self.start
        sortable_end = self.end
        if not self.is_plus_strand:
            sortable_start = sortable_start * -1
            sortable_end = sortable_end * -1
        return self.coord.seqid, self.is_plus_strand, sortable_start, sortable_end, self.feature_type

    def __repr__(self):
        params = {
            'id': self.id,
            'coord_id': self.coord.id,
            'type': self.feature_type,
            'is_plus_strand': self.is_plus_strand,
            'phase': self.phase,
        }
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

    def add_to_queue(self):
        transcribed = self._get_params_dict()
        self.controller.insertion_queues.transcribed.queue.append(transcribed)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'type': self.entry_type,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id,
        }
        return d

    def __repr__(self):
        return self._get_repr('TranscriptHandler', self._get_params_dict())


class TranscribedPieceHandler(handlers.TranscribedPieceHandlerBase, Insertable):
    def __init__(self, given_name, transcript_id, position, controller):
        super().__init__()
        self.id = InsertCounterHolder.transcribed_piece()
        self.given_name = given_name
        self.transcript_id = transcript_id
        self.position = position
        self.controller = controller

    def add_to_queue(self):
        transcribed_piece = self._get_params_dict()
        self.controller.insertion_queues.transcribed_piece.queue.append(transcribed_piece)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'given_name': self.given_name,
            'transcribed_id': self.transcript_id,
            'position': self.position,
        }
        return d

    def __repr__(self):
        return self._get_repr('TranscribedPieceHandler', self._get_params_dict())


class TranslatedHandler(handlers.TranslatedHandlerBase):
    def __init__(self, given_name, super_locus_id, controller):
        super().__init__()
        self.id = InsertCounterHolder.translated()
        self.given_name = given_name
        self.super_locus_id = super_locus_id
        self.controller = controller

    def add_to_queue(self):
        translated = self._get_params_dict()
        self.controller.insertion_queues.translated.queue.append(translated)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id,
        }
        return d

    def __repr__(self):
        return self._get_repr('ProteinHandler', self._get_params_dict())
