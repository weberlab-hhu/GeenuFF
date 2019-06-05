import os
import logging
from pprint import pprint  # for debugging
from abc import ABC, abstractmethod
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from dustdas import gffhelper, fastahelper
from .. import orm
from .. import types
from .. import helpers
from ..base.helpers import (get_strand_direction, get_geenuff_start_end, has_start_codon,
                            has_stop_codon)


# core queue prep
class InsertionQueue(helpers.QueueController):
    def __init__(self, session, engine):
        super().__init__(session, engine)
        self.super_locus = helpers.CoreQueue(orm.SuperLocus.__table__.insert())
        self.transcript = helpers.CoreQueue(orm.Transcript.__table__.insert())
        self.transcript_piece = helpers.CoreQueue(orm.TranscriptPiece.__table__.insert())
        self.protein = helpers.CoreQueue(orm.Protein.__table__.insert())
        self.feature = helpers.CoreQueue(orm.Feature.__table__.insert())
        self.association_transcript_piece_to_feature = helpers.CoreQueue(
            orm.association_transcript_piece_to_feature.insert())
        self.association_protein_to_feature = helpers.CoreQueue(
            orm.association_protein_to_feature.insert())

        self.ordered_queues = [
            self.super_locus, self.transcript, self.transcript_piece, self.protein, self.feature,
            self.association_transcript_piece_to_feature, self.association_protein_to_feature
        ]


class InsertCounterHolder(object):
    """provides incrementing unique integers to be used as primary keys for bulk inserts"""
    feature = helpers.Counter(orm.Feature)
    protein = helpers.Counter(orm.Protein)
    transcript = helpers.Counter(orm.Transcript)
    super_locus = helpers.Counter(orm.SuperLocus)
    transcript_piece = helpers.Counter(orm.TranscriptPiece)
    genome = helpers.Counter(orm.Genome)

    @staticmethod
    def sync_counters_with_db(session):
        InsertCounterHolder.feature.sync_with_db(session)
        InsertCounterHolder.protein.sync_with_db(session)
        InsertCounterHolder.transcript.sync_with_db(session)
        InsertCounterHolder.super_locus.sync_with_db(session)
        InsertCounterHolder.transcript_piece.sync_with_db(session)
        InsertCounterHolder.genome.sync_with_db(session)


class OrganizedGeenuffImporterGroup(object):
    """Stores the handler objects for a super locus in an organized fassion.
    The format is similar to the one of OrganizedGFFEntryGroup, but it stores objects
    according to the Geenuff way of saving genomic annotations. This format can then
    be checked for errors and changed accordingly before being inserted into the db.

    The importers are organized in the following way:

    importers = {
        'super_locus' = super_locus_importer,
        'transcripts': [
            {
                'transcript': transcript_importer,
                'transcript_piece: transcript_piece_importer,
                'transcript_feature': transcript_importer,
                'protein': protein_importer,
                'cds': cds_importer,
                'introns': [intron_importer1, intron_importer2, ..]
            },
            ...
        ],
        'errors' = [error_importer1, error_importer2, ..],
    }
    """

    def __init__(self, organized_gff_entries, coord, controller):
        self.coord = coord
        self.controller = controller
        self.importers = {'transcripts': [], 'errors': []}
        self._parse_gff_entries(organized_gff_entries)

    def _parse_gff_entries(self, entries):
        """Changes the GFF format into the GeenuFF format. Does all the parsing."""
        sl = entries['super_locus']
        # these attr are the same for all geenuff importers going to be created
        is_plus_strand = get_strand_direction(sl)
        score = sl.score
        source = sl.source

        sl_start, sl_end = get_geenuff_start_end(sl.start, sl.end, is_plus_strand)
        sl_i = self.importers['super_locus'] = SuperLocusImporter(entry_type=sl.type,
                                                                  given_name=sl.get_ID(),
                                                                  coord=self.coord,
                                                                  is_plus_strand=is_plus_strand,
                                                                  start=sl_start,
                                                                  end=sl_end,
                                                                  controller=self.controller)
        for t, t_entries in entries['transcripts'].items():
            t_importers = {}
            # check for multi inheritance and throw NotImplementedError if found
            if len(t.get_Parent()) > 1:
                raise NotImplementedError
            t_id = t.get_ID()
            # create transcript handler
            t_i = TranscriptImporter(entry_type=t.type,
                                     given_name=t_id,
                                     super_locus_id=sl_i.id,
                                     controller=self.controller)
            # create transcript piece handler
            tp_i = TranscriptPieceImporter(given_name=t_id,
                                           transcript_id=t_i.id,
                                           position=0,
                                           controller=self.controller)
            # create transcript feature handler
            tf_i = FeatureImporter(self.coord,
                                   is_plus_strand,
                                   types.TRANSCRIPT_FEATURE,
                                   given_name=t_id,
                                   score=score,
                                   source=source,
                                   controller=self.controller)
            tf_i.set_start_end_from_gff(t.start, t.end)

            # insert everything so far into the dict
            t_importers['transcript'] = t_i
            t_importers['transcript_piece'] = tp_i
            t_importers['transcript_feature'] = tf_i

            # if it is not a non-coding gene or something like that
            if t_entries['cds']:
                # create protein handler
                protein_id = self._get_protein_id_from_cds_list(t_entries['cds'])
                p_i = ProteinImporter(given_name=protein_id,
                                      super_locus_id=sl_i.id,
                                      controller=self.controller)
                # create coding features from exon limits
                if is_plus_strand:
                    phase_5p = t_entries['cds'][0].phase
                    phase_3p = t_entries['cds'][-1].phase
                else:
                    phase_5p = t_entries['cds'][-1].phase
                    phase_3p = t_entries['cds'][0].phase
                cds_i = FeatureImporter(self.coord,
                                        is_plus_strand,
                                        types.CODING,
                                        phase_5p=phase_5p,
                                        phase_3p=phase_3p,
                                        score=score,
                                        source=source,
                                        controller=self.controller)
                gff_start = t_entries['cds'][0].start
                gff_end = t_entries['cds'][-1].end
                cds_i.set_start_end_from_gff(gff_start, gff_end)

                # insert everything so far into the dict
                t_importers['protein'] = p_i
                t_importers['cds'] = cds_i

                # create all the introns by traversing the exon entries and insert FeatureImporters
                # into the previously created list
                # the introns should strictly always lie between successive gff exons
                exons = t_entries['exons']
                introns = []
                for i in range(len(exons) - 1):
                    # the introns are delimited by the surrounding exons
                    # the first base of an intron in right after the last exononic base
                    gff_start = exons[i].end + 1
                    gff_end = exons[i + 1].start - 1
                    # ignore introns that would come from directly adjacent exons
                    if gff_start - gff_end != 1:
                        intron_i = FeatureImporter(self.coord,
                                                   is_plus_strand,
                                                   types.INTRON,
                                                   score=score,
                                                   source=source,
                                                   controller=self.controller)
                        intron_i.set_start_end_from_gff(gff_start, gff_end)
                        introns.append(intron_i)
                if is_plus_strand:
                    t_importers['introns'] = introns
                else:
                    t_importers['introns'] = introns[::-1]
            self.importers['transcripts'].append(t_importers)

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
            protein_id = OrganizedGeenuffImporterGroup._get_protein_id_from_cds_entry(cds_entry)
            if protein_id:
                protein_ids.add(protein_id)
        if len(protein_ids) != 1:
            raise ValueError('No protein_id or more than one protein_ids for one transcript')
        return protein_ids.pop()


class OrganizedGFFEntryGroup(object):
    """Takes an entry group (all entries of one super locus) and stores the entries
    in an orderly fassion. Can then return a corresponding OrganizedGeenuffImporterGroup.
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

    def __init__(self, gff_entry_group, fasta_importer, controller):
        self.fasta_importer = fasta_importer
        self.controller = controller
        self.entries = {'transcripts': {}}
        self.coord = None
        self.add_gff_entry_group(gff_entry_group)

    def add_gff_entry_group(self, entries):
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
        self.coord = self.fasta_importer.gffid_to_coords[self.entries['super_locus'].seqid]

        # order exon and cds lists by start value (disregard strand for now)
        for _, value_dict in self.entries['transcripts'].items():
            for key in ['exons', 'cds']:
                value_dict[key].sort(key=lambda e: e.start)

    def get_geenuff_importers(self):
        geenuff_importer_group = OrganizedGeenuffImporterGroup(self.entries, self.coord,
                                                               self.controller)
        return geenuff_importer_group.importers


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
        self.organized_entries = {}

    def load_organized_entries(self):
        self.organized_entries = {}
        gene_level = [x.value for x in types.SuperLocusAll]

        reader = self._useful_gff_entries()
        first = next(reader)
        seqid = first.seqid
        gene_group = [first]
        self.organized_entries[seqid] = []
        for entry in reader:
            if entry.type in gene_level:
                self.organized_entries[seqid].append(gene_group)
                gene_group = [entry]
                if entry.seqid != seqid:
                    self.organized_entries[entry.seqid] = []
                    seqid = entry.seqid
            else:
                gene_group.append(entry)
        self.organized_entries[seqid].append(gene_group)

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
    in the space of GeenuFF importers.
    Assumes all super locus handler groups to be ordered 5p to 3p and of one strand.
    Works with a list of OrganizedGeenuffImporterGroup, which correspond to a list of
    super loci, and looks for errors. Error features may be inserted and importers be
    removed when deemed necessary.
    """

    def __init__(self, geenuff_importer_groups, controller):
        self.groups = geenuff_importer_groups
        if self.groups:
            self.is_plus_strand = self.groups[0]['super_locus'].is_plus_strand
            self.coord = self.groups[0]['super_locus'].coord
            # make sure self.groups is sorted correctly
            self.groups.sort(key=lambda g: g['super_locus'].start, reverse=not self.is_plus_strand)
        self.controller = controller

    def _3p_cds_start(self, transcript):
        """returns the start of the 3p most cds feature"""
        cds = transcript['cds']
        start = cds.start
        # introns are ordered by coordinate with no respect to strand
        for intron in transcript['introns']:
            if ((self.is_plus_strand and intron.end < cds.end)
                    or (not self.is_plus_strand and intron.end > cds.end)):
                start = intron.end
            else:
                break
        return start

    def resolve_errors(self):
        for i, group in enumerate(self.groups):
            # the case of no transcript for a super locus
            # solution is to add an error mask that extends halfway to the appending super
            # loci in the intergenic region as something in this area appears to have gone wrong
            self.errors = group['errors']
            if not group['transcripts']:
                self._add_overlapping_error(i, None, 'whole', types.EMPTY_SUPER_LOCUS)
            # other cases
            for transcript in group['transcripts']:
                # if coding transcript
                if 'cds' in transcript:
                    cds = transcript['cds']
                    introns = transcript['introns']

                    # the case of missing of implicit UTR ranges
                    # the solution is similar to the one above
                    if cds.start == transcript['transcript_feature'].start:
                        self._add_overlapping_error(i, cds, '5p', types.MISSING_UTR_5P)
                    if cds.end == transcript['transcript_feature'].end:
                        self._add_overlapping_error(i, cds, '3p', types.MISSING_UTR_3P)

                    # the case of missing start/stop codon
                    if not has_start_codon(cds.coord.sequence, cds.start, self.is_plus_strand):
                        self._add_overlapping_error(i, cds, '5p', types.MISSING_START_CODON)
                    if not has_stop_codon(cds.coord.sequence, cds.end, self.is_plus_strand):
                        self._add_overlapping_error(i, cds, '3p', types.MISSING_STOP_CODON)

                    # the case of wrong 5p phase
                    if cds.phase_5p != 0:
                        self._add_overlapping_error(i, cds, '5p', types.WRONG_PHASE_5P)

                    if introns:
                        # the case of wrong 3p phase
                        if transcript['transcript_feature'].given_name == 'Aco009693.1.v3':
                            import pudb; pudb.set_trace()
                        len_3p_exon = abs(cds.end - self._3p_cds_start(transcript))
                        if cds.phase_3p != len_3p_exon % 3:
                            self._add_overlapping_error(i, cds, '3p', types.MISMATCHED_PHASE_3P)

                        faulty_introns = []
                        for j, intron in enumerate(introns):
                            # the case of overlapping exons
                            if ((self.is_plus_strand and intron.end < intron.start)
                                    or (not self.is_plus_strand and intron.end > intron.start)):
                                # mark the overlapping cds regions as errors
                                if j > 0:
                                    error_start = introns[j - 1].end
                                else:
                                    error_start = transcript['transcript_feature'].start
                                if j < len(introns) - 1:
                                    error_end = introns[j + 1].start
                                else:
                                    error_end = transcript['transcript_feature'].end
                                self._add_error(i, error_start, error_end, self.is_plus_strand,
                                                types.OVERLAPPING_EXONS)
                                faulty_introns.append(intron)
                            # the case of a too short intron
                            # todo put the minimum length in a config somewhere
                            elif abs(intron.end - intron.start) < 60:
                                self._add_error(i, intron.start, intron.end, self.is_plus_strand,
                                                types.TOO_SHORT_INTRON)
                                faulty_introns.append(intron)
                        # do not save faulty introns, the error should be descriptive enough
                        for intron in faulty_introns:
                            introns.remove(intron)

    def _add_error(self, i, start, end, is_plus_strand, error_type):
        if is_plus_strand:
            strand_str = 'plus'
        else:
            strand_str = 'minus'
        if (is_plus_strand and end < start) or (not is_plus_strand and end > start):
            msg = ('skipping error in nested gene: seqid: {seqid}, {geneid}, on {strand} strand, '
                   'with type: {type}').format(seqid=self.coord.seqid,
                                               geneid=self.groups[i]['super_locus'].given_name,
                                               strand=strand_str,
                                               type=error_type)
            logging.warning(msg)
        else:
            error_i = FeatureImporter(self.coord,
                                      is_plus_strand,
                                      error_type,
                                      start=start,
                                      end=end,
                                      controller=self.controller)
            self.errors.append(error_i)
            # error msg
            msg = ('marked as erroneous: seqid: {seqid}, {start}--{end}:{geneid}, on {strand} strand, '
                   'with type: {type}').format(seqid=self.coord.seqid,
                                               start=start,
                                               end=end,
                                               geneid=self.groups[i]['super_locus'].given_name,
                                               strand=strand_str,
                                               type=error_type)
            logging.warning(msg)

    def _add_overlapping_error(self, i, handler, direction, error_type):
        """Constructs an error features that overlaps halfway to the next super locus
        in the given direction from the given handler if possible. Otherwise mark until the end.
        If the direction is 'whole', the handler parameter is ignored.

        Also sets handler.start_is_biological_start=False (or the end).
        """
        assert direction in ['5p', '3p', 'whole']
        sl = self.groups[i]['super_locus']

        # set correct upstream error starting point
        if direction in ['5p', 'whole']:
            if i > 0:
                sl_prev = self.groups[i - 1]['super_locus']
                anchor_5p = self._halfway_mark(sl_prev, sl)
            else:
                if self.is_plus_strand:
                    anchor_5p = 0
                else:
                    anchor_5p = sl.coord.end

        # set correct downstream error end point
        if direction in ['3p', 'whole']:
            if i < len(self.groups) - 1:
                sl_next = self.groups[i + 1]['super_locus']
                anchor_3p = self._halfway_mark(sl, sl_next)
            else:
                if self.is_plus_strand:
                    anchor_3p = sl.coord.end
                else:
                    anchor_3p = -1

        if direction == '5p':
            error_5p = anchor_5p
            error_3p = handler.start
            handler.start_is_biological_start = False
        elif direction == '3p':
            error_5p = handler.end
            error_3p = anchor_3p
            handler.end_is_biological_end = False
        elif direction == 'whole':
            error_5p = anchor_5p
            error_3p = anchor_3p

        if not self._zero_len_coords_at_sequence_edge(error_5p, error_3p, direction, sl.coord):
            self._add_error(i, error_5p, error_3p, self.is_plus_strand, error_type)

    def _zero_len_coords_at_sequence_edge(self, error_5p, error_3p, direction, coordinate):
        """Check if error 5p-3p is of zero length due to hitting start or end of sequence"""
        out = False
        if self.is_plus_strand:
            if direction == '5p':
                if error_5p == error_3p == 0:
                    out = True
            elif direction == '3p':
                if error_5p == error_3p == coordinate.end:
                    out = True
        else:
            if direction == '5p':
                if error_5p == error_3p == coordinate.end - 1:
                    out = True
            elif direction == '3p':
                if error_5p == error_3p == -1:
                    out = True
        return out

    def _halfway_mark(self, sl, sl_next):
        """Calculates the half way point between two super loci, which is then used for
        error masks.
        """
        if self.is_plus_strand:
            dist = sl_next.start - sl.end
            mark = sl.end + dist // 2
        else:
            dist = sl.end - sl_next.start
            mark = sl.end - dist // 2
        return mark


##### main flow control #####
class ImportController(object):
    def __init__(self, database_path, replace_db=False):
        self.database_path = database_path
        self.latest_genome = None
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

    def make_genome(self, genome_args={}):
        genome = orm.Genome(**genome_args)
        self.latest_fasta_importer = FastaImporter(genome)
        self.session.add(genome)
        self.session.commit()

    def add_genome(self, fasta_path, gff_path, genome_args={}, clean_gff=True):
        self.clean_tmp_data()
        self.add_sequences(fasta_path, genome_args)
        self.add_gff(gff_path, clean=clean_gff)

    def add_sequences(self, seq_path, genome_args={}):
        if self.latest_genome is None:
            self.make_genome(genome_args)

        self.latest_fasta_importer.add_sequences(seq_path)
        self.session.commit()

    def clean_tmp_data(self):
        self.latest_genome = None
        self.latest_super_loci = []

    def add_gff(self, gff_file, clean=True):
        def insert_importer_groups(self, groups):
            """Initiates the calling of the add_to_queue() function of the importers
            in the correct order. Also initiates the insert of the many2many rows.
            """
            for group in groups:
                group['super_locus'].add_to_queue()
                # insert all features as well as transcript and protein related entries
                for transcripts in group['transcripts']:
                    # make shortcuts
                    tp = transcripts['transcript_piece']
                    tf = transcripts['transcript_feature']
                    # add transcript handler that are always present
                    transcripts['transcript'].add_to_queue()
                    tp.add_to_queue()
                    tf.add_to_queue()
                    tf.insert_feature_piece_association(tp.id)
                    # if coding transcript
                    if 'protein' in transcripts:
                        transcripts['protein'].add_to_queue()
                        tf.insert_feature_protein_association(transcripts['protein'].id)
                        transcripts['cds'].add_to_queue()
                        transcripts['cds'].insert_feature_piece_association(tp.id)
                        for intron in transcripts['introns']:
                            intron.add_to_queue()
                            intron.insert_feature_piece_association(tp.id)
                # insert the errors
                for error in group['errors']:
                    error.add_to_queue()

        def clean_and_insert(self, groups, clean):
            plus = [g for g in groups if g['super_locus'].is_plus_strand]
            minus = [g for g in groups if not g['super_locus'].is_plus_strand]
            if clean:
                # check and correct for errors
                # do so for each strand seperately
                # all changes should be made by reference
                GFFErrorHandling(plus, self).resolve_errors()
                # reverse order on minus strand
                GFFErrorHandling(minus[::-1], self).resolve_errors()
            # insert importers
            insert_importer_groups(self, plus)
            insert_importer_groups(self, minus)
            self.insertion_queues.execute_so_far()

        assert self.latest_fasta_importer is not None, 'No recent genome found'
        self.latest_fasta_importer.mk_mapper(gff_file)
        gff_organizer = OrganizedGFFEntries(gff_file)
        gff_organizer.load_organized_entries()
        organized_gff_entries = gff_organizer.organized_entries
        geenuff_importer_groups = []
        for seqid in organized_gff_entries.keys():
            for entry_group in organized_gff_entries[seqid]:
                organized_entries = OrganizedGFFEntryGroup(entry_group, self.latest_fasta_importer,
                                                           self)
                geenuff_importer_groups.append(organized_entries.get_geenuff_importers())
            # never do error checking across fasta sequence borders
            clean_and_insert(self, geenuff_importer_groups, clean)
            geenuff_importer_groups = []


class Insertable(ABC):
    @abstractmethod
    def add_to_queue(self):
        pass


class FastaImporter(object):
    def __init__(self, genome):
        self.genome = genome
        self.mapper = None
        self._coords_by_seqid = None
        self._gffid_to_coords = None
        self._gff_seq_ids = None

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
            self._coords_by_seqid = {c.seqid: c for c in self.genome.coordinates}
        return self._coords_by_seqid

    def mk_mapper(self, gff_file=None):
        fa_ids = [e.seqid for e in self.genome.coordinates]
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
                                   genome=self.genome)


class SuperLocusImporter(Insertable):
    def __init__(self,
                 entry_type,
                 given_name,
                 controller,
                 coord=None,
                 is_plus_strand=None,
                 start=-1,
                 end=-1):
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
        return helpers.get_repr('SuperLocusImporter', params)


class FeatureImporter(Insertable):
    def __init__(self,
                 coord,
                 is_plus_strand,
                 feature_type,
                 controller,
                 start=-1,
                 end=-1,
                 given_name=None,
                 phase_5p=0,
                 phase_3p=0,
                 score=None,
                 source=None):
        self.id = InsertCounterHolder.feature()
        self.coord = coord
        self.given_name = given_name
        self.is_plus_strand = is_plus_strand
        self.feature_type = feature_type
        # start/end may have to be adapted to geenuff
        self.start = start
        self.end = end
        self.phase_5p = phase_5p
        self.phase_3p = phase_3p  # only used for error checking
        self.score = score
        self.source = source
        self.start_is_biological_start = True
        self.end_is_biological_end = True
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
            'phase': self.phase_5p,
            'start': self.start,
            'end': self.end,
            'start_is_biological_start': self.start_is_biological_start,
            'end_is_biological_end': self.end_is_biological_end,
        }
        self.controller.insertion_queues.feature.queue.append(feature)

    def insert_feature_piece_association(self, transcript_piece_id):
        features2pieces = {
            'feature_id': self.id,
            'transcript_piece_id': transcript_piece_id,
        }
        self.controller.insertion_queues.association_transcript_piece_to_feature.\
            queue.append(features2pieces)

    def insert_feature_protein_association(self, protein_id):
        features2protein = {
            'feature_id': self.id,
            'protein_id': protein_id,
        }
        self.controller.insertion_queues.association_protein_to_feature.\
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
            'phase': self.phase_5p,
        }
        if self.given_name:
            params['given_name'] = self.given_name
        return helpers.get_repr('FeatureImporter', params, str(self.start) + '--' + str(self.end))


class TranscriptImporter(Insertable):
    def __init__(self, entry_type, given_name, super_locus_id, controller):
        self.id = InsertCounterHolder.transcript()
        self.entry_type = entry_type
        self.given_name = given_name
        self.super_locus_id = super_locus_id
        self.controller = controller

    def add_to_queue(self):
        transcript = self._get_params_dict()
        self.controller.insertion_queues.transcript.queue.append(transcript)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'type': self.entry_type,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id,
        }
        return d

    def __repr__(self):
        return helpers.get_repr('TranscriptImporter', self._get_params_dict())


class TranscriptPieceImporter(Insertable):
    def __init__(self, given_name, transcript_id, position, controller):
        self.id = InsertCounterHolder.transcript_piece()
        self.given_name = given_name
        self.transcript_id = transcript_id
        self.position = position
        self.controller = controller

    def add_to_queue(self):
        transcript_piece = self._get_params_dict()
        self.controller.insertion_queues.transcript_piece.queue.append(transcript_piece)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'given_name': self.given_name,
            'transcript_id': self.transcript_id,
            'position': self.position,
        }
        return d

    def __repr__(self):
        return helpers.get_repr('TranscriptPieceImporter', self._get_params_dict())


class ProteinImporter(Insertable):
    def __init__(self, given_name, super_locus_id, controller):
        self.id = InsertCounterHolder.protein()
        self.given_name = given_name
        self.super_locus_id = super_locus_id
        self.controller = controller

    def add_to_queue(self):
        protein = self._get_params_dict()
        self.controller.insertion_queues.protein.queue.append(protein)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id,
        }
        return d

    def __repr__(self):
        return helpers.get_repr('ProteinImporter', self._get_params_dict())
