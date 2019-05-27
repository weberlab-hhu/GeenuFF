import os
import pytest
from dustdas import gffhelper
from sqlalchemy import create_engine
from sqlalchemy import func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError
import sqlalchemy

from .. import orm
from .. import types
from .. import handlers
from .. import helpers
from ..base.orm import (Genome, Feature, Coordinate, Transcribed, TranscribedPiece, SuperLocus,
                        Translated)
from ..base.handlers import SuperLocusHandlerBase, TranscribedHandlerBase
from ..applications import importer
from ..applications.importer import ImportController, InsertCounterHolder, OrganizedGFFEntries


@pytest.fixture(scope="session", autouse=True)
def prepare(request):
    if not os.getcwd().endswith('GeenuFF/geenuff'):
        pytest.exit('Tests need to be run from GeenuFF/geenuff directory')


### Helper functions ###
def mk_memory_session(db_path='sqlite:///:memory:'):
    engine = create_engine(db_path, echo=False)
    orm.Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()


def setup_data_handler(handler_type, data_type, **kwargs):
    data = data_type(**kwargs)
    handler = handler_type()
    handler.add_data(data)
    return data, handler


def setup_dummyloci_super_locus(db_path='sqlite:///:memory:'):
    controller = ImportController(database_path=db_path)
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff')
    sl = controller.session.query(SuperLocus).one()
    return sl, controller


def cleaned_commited_features(sess):
    all_features = sess.query(Feature).all()
    allowed_types = [types.TRANSCRIBED, types.CODING, types.INTRON, types.TRANS_INTRON, types.ERROR]
    clean_datas = [x for x in all_features if x.type.value in allowed_types]
    return clean_datas

def orm_object_in_list(obj, obj_list):
    """Checks if obj is in obj_list (by implicitely calling __eq__()) and
    removes obj from obj_list if it has been found. Return whether it was found.
    """
    for o in obj_list[:]:  # make a copy at each iteration so we avoid weird errors
        if o == obj:
            obj_list.remove(o)
            return True
    return False


### The actual tests ###
def test_annogenome2coordinate_relation():
    """Check if everything is consistent when we add an Genome and a Coordinate
    to the db. Also check for correct deletion behavior.
    """
    sess = mk_memory_session()
    g = Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coord = Coordinate(start=0, end=30, seqid='abc', genome=g)
    assert g is coord.genome
    # actually put everything in db
    sess.add(coord)
    sess.commit()
    # check primary keys were assigned
    assert g.id == 1
    assert coord.id == 1
    # check we can access coordinates from g
    coord_q = g.coordinates[0]
    assert coord is coord_q
    assert g is coord.genome
    assert g.id == coord_q.genome_id
    # check we get logical behavior on deletion
    sess.delete(coord)
    sess.commit()
    assert len(g.coordinates) == 0
    print(coord.genome)
    sess.delete(g)
    sess.commit()
    with pytest.raises(sqlalchemy.exc.InvalidRequestError):
        sess.add(coord)


def test_coordinate_constraints():
    """Check the coordinate constraints"""
    sess = mk_memory_session()
    g = Genome()

    # should be ok
    coors = Coordinate(start=0, end=30, seqid='abc', genome=g)
    coors2 = Coordinate(start=0, end=1, seqid='abcd', genome=g)
    sess.add_all([coors, coors2])
    sess.commit()

    # start/end constraints
    coors_bad1 = Coordinate(start=-12, end=30, seqid='abc', genome=g)
    coors_bad2 = Coordinate(start=100, end=30, seqid='abcd', genome=g)
    coors_bad3 = Coordinate(start=1, end=30, genome=g)
    with pytest.raises(IntegrityError):
        sess.add(coors_bad1)  # start below 1
        sess.commit()
    sess.rollback()
    with pytest.raises(IntegrityError):
        sess.add(coors_bad2)  # end below start
        sess.commit()
    sess.rollback()
    with pytest.raises(IntegrityError):
        sess.add(coors_bad3)
        sess.commit()


def test_coordinate_insert():
    """Test what happens when we insert two coordinates"""
    sess = mk_memory_session()
    g = Genome()
    coords = Coordinate(start=1, end=30, seqid='abc', genome=g)
    coords2 = Coordinate(start=11, end=330, seqid='def', genome=g)
    sl = SuperLocus()
    f0 = Feature(coordinate=coords)
    f1 = Feature(coordinate=coords2)
    # should be ok
    sess.add_all([g, sl, coords, coords2, f0, f1])
    assert f0.coordinate.start == 1
    assert f1.coordinate.end == 330


def test_many2many_with_features():
    """Test the many2many tables association_transcribed_piece_to_feature and
    association_translated_to_feature
    """
    sl = SuperLocus()
    # one transcript, multiple proteins
    piece0 = TranscribedPiece()
    scribed0 = Transcribed(super_locus=sl, transcribed_pieces=[piece0])
    slated0 = Translated(super_locus=sl)
    slated1 = Translated(super_locus=sl)
    # features representing alternative start codon for proteins on one transcript
    feat0_tss = Feature(transcribed_pieces=[piece0])
    feat1_tss = Feature(transcribed_pieces=[piece0])
    feat2_stop = Feature(translateds=[slated0, slated1])
    feat3_start = Feature(translateds=[slated0])
    feat4_start = Feature(translateds=[slated1])
    # test multi features per translated worked
    assert len(slated0.features) == 2
    # test mutli translated per feature worked
    assert len(feat2_stop.translateds) == 2
    assert len(feat3_start.translateds) == 1
    assert len(feat0_tss.translateds) == 0


def test_feature_has_its_things():
    """Test if properties of Feature table are correct and constraints are enforced"""
    sess = mk_memory_session()
    # test feature with nothing much set
    g = Genome()
    c = Coordinate(start=1, end=30, seqid='abc', genome=g)
    f = Feature(coordinate=c)
    sess.add_all([f, c])
    sess.commit()

    assert f.is_plus_strand is None
    assert f.source is None
    assert f.score is None
    # test feature with
    f1 = Feature(is_plus_strand=False, start=3, end=-1, coordinate=c)
    assert not f1.is_plus_strand
    assert f1.start == 3
    assert f1.end == -1
    sess.add(f1)
    sess.commit()
    # test too low of start / end coordinates raise an error
    f_should_fail = Feature(is_plus_strand=True, start=-5, end=10, coordinate=c)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()

    f_should_fail = Feature(is_plus_strand=False, start=5, end=-2, coordinate=c)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()
    # test wrong class as parameter
    with pytest.raises(KeyError):
        f2 = Feature(coordinate=f)

    f2 = Feature(is_plus_strand=-1)  # note that 0, and 1 are accepted
    sess.add(f2)
    with pytest.raises(sqlalchemy.exc.StatementError):
        sess.commit()
    sess.rollback()


def test_partially_remove_coordinate():
    """Add two coordinates to an annotated genome and test if everything ends
    up valid.
    """
    sess = mk_memory_session()
    g = Genome()
    place_holder = Genome()
    coord0 = Coordinate(start=1, end=30, seqid='abc', genome=g)
    coord1 = Coordinate(start=11, end=330, seqid='def', genome=g)
    sess.add_all([g, coord0, coord1])
    sess.commit()
    assert len(g.coordinates) == 2
    g.coordinates.remove(coord0)
    coord0.genome = place_holder  # else we'll fail the not NULL constraint
    sess.commit()
    # removed from g
    assert len(g.coordinates) == 1
    # but still in table
    assert len(sess.query(Coordinate).all()) == 2


# section: api
def test_transcribed_piece_unique_constraints():
    """Add transcribed pieces in valid and invalid configurations and test for
    valid outcomes.
    """
    sess = mk_memory_session()
    sl = SuperLocus()
    transcribed0 = Transcribed(super_locus=sl)
    transcribed1 = Transcribed(super_locus=sl)

    # test if same position for different transcribed_id goes through
    piece_tr0_pos0 = TranscribedPiece(transcribed=transcribed0, position=0)
    piece_tr1_pos0 = TranscribedPiece(transcribed=transcribed1, position=0)
    sess.add_all([transcribed0, piece_tr0_pos0, piece_tr1_pos0])
    sess.commit()

    # same transcibed_id but different position
    piece_tr0_pos1 = TranscribedPiece(transcribed=transcribed0, position=1)
    sess.add(piece_tr0_pos1)
    sess.commit()

    # test if unique constraint works
    piece_tr0_pos1_2nd = TranscribedPiece(transcribed=transcribed0, position=1)
    sess.add(piece_tr0_pos1_2nd)
    with pytest.raises(IntegrityError):
        sess.commit()


def test_order_pieces():
    """Add transcribed pieces that consists of one or many features to the db and test
    if the order for the transcribed pieces and the features is returned according to the
    position property instead of db insertion order.
    """
    sess = mk_memory_session()
    g = Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coor = Coordinate(seqid='a', start=1, end=1000, genome=g)
    sess.add_all([g, coor])
    sess.commit()
    # setup one transcribed handler with pieces
    sl, sl_h = setup_data_handler(SuperLocusHandlerBase, SuperLocus)
    t, t_h = setup_data_handler(TranscribedHandlerBase, Transcribed, super_locus=sl)
    # insert in wrong order
    piece1 = TranscribedPiece(position=1)
    piece0 = TranscribedPiece(position=0)
    piece2 = TranscribedPiece(position=2)
    t.transcribed_pieces = [piece0, piece1, piece2]
    sess.add_all([t, piece1, piece0, piece2])
    sess.commit()
    # setup some features on different pieces
    feature0u = Feature(transcribed_pieces=[piece0],
                        coordinate=coor,
                        start=100,
                        end=200,
                        given_name='0u',
                        is_plus_strand=True)
    feature1d = Feature(transcribed_pieces=[piece1],
                        coordinate=coor,
                        start=1,
                        end=3,
                        given_name='1d',
                        is_plus_strand=True)
    feature1u = Feature(transcribed_pieces=[piece1],
                        coordinate=coor,
                        start=100,
                        end=200,
                        given_name='1u',
                        is_plus_strand=True)
    feature2d = Feature(transcribed_pieces=[piece2],
                        coordinate=coor,
                        start=1,
                        end=3,
                        given_name='2d',
                        is_plus_strand=True)
    # see if they can be ordered as expected overall
    op = ti.sorted_pieces()
    print([piece0, piece1, piece2], 'expected')
    print(op, 'sorted')
    assert op == [piece0, piece1, piece2]
    # order features by piece
    fully_sorted = ti.sort_all()
    expected = [[feature0u], [feature1d, feature1u], [feature2d]]
    assert fully_sorted == expected
    # test assertion of matching strand
    feature_neg = Feature(transcribed_pieces=[piece2],
                          coordinate=coor,
                          start=5,
                          end=10,
                          given_name='neg',
                          is_plus_strand=False)
    sess.add(feature_neg)
    sess.commit()
    with pytest.raises(AssertionError):
        ti.sorted_features(piece2)


class BiointerpDemoDataCoding(object):
    def __init__(self, sess, is_plus_strand):
        self.sl, self.sl_handler = setup_data_handler(handlers.SuperLocusHandlerBase, SuperLocus)
        self.scribed, self.scribed_handler = setup_data_handler(handlers.TranscribedHandlerBase,
                                                                Transcribed,
                                                                super_locus=self.sl)

        self.piece = TranscribedPiece(position=0, transcribed=self.scribed)

        self.ti = TranscriptInterpBase(transcript=self.scribed_handler,
                                       super_locus=self.sl,
                                       session=sess)

        self.g = Genome()
        # setup ranges for a two-exon coding gene
        self.coordinate = Coordinate(seqid='a', start=1, end=2000, genome=self.g)
        transcribed_start, transcribed_end = 100, 900
        coding_start, coding_end = 200, 800
        intron_start, intron_end = 300, 700

        if not is_plus_strand:  # swap if we're setting up data for "-" strand
            transcribed_end, transcribed_start = transcribed_start, transcribed_end
            coding_end, coding_start = coding_start, coding_end
            intron_end, intron_start = intron_start, intron_end

        # transcribed:
        self.transcribed_feature = Feature(coordinate=self.coordinate,
                                           is_plus_strand=is_plus_strand,
                                           start=transcribed_start,
                                           end=transcribed_end,
                                           start_is_biological_start=True,
                                           end_is_biological_end=True,
                                           type=types.TRANSCRIBED)
        # coding:
        self.coding_feature = Feature(coordinate=self.coordinate,
                                      is_plus_strand=is_plus_strand,
                                      start=coding_start,
                                      end=coding_end,
                                      start_is_biological_start=True,
                                      end_is_biological_end=True,
                                      type=types.CODING)
        # intron:
        self.intron_feature = Feature(coordinate=self.coordinate,
                                      is_plus_strand=is_plus_strand,
                                      start=intron_start,
                                      end=intron_end,
                                      start_is_biological_start=True,
                                      end_is_biological_end=True,
                                      type=types.INTRON)

        self.piece.features = [self.transcribed_feature, self.coding_feature, self.intron_feature]

        sess.add(self.sl)
        sess.commit()


class BiointerpDemoDataPartial(object):
    def __init__(self, sess, is_plus_strand):
        self.sl, self.sl_handler = setup_data_handler(handlers.SuperLocusHandlerBase, SuperLocus)
        self.scribed, self.scribed_handler = setup_data_handler(handlers.TranscribedHandlerBase,
                                                                Transcribed,
                                                                super_locus=self.sl)

        self.piece = TranscribedPiece(position=0, transcribed=self.scribed)

        self.ti = TranscriptInterpBase(transcript=self.scribed_handler,
                                       super_locus=self.sl,
                                       session=sess)

        self.ag = Genome()
        # setup ranges for a two-exon coding gene
        self.coordinates = Coordinate(seqid='a', start=1, end=2000, genome=self.ag)
        transcribed_start, transcribed_end = 100, 500
        error_start, error_end = 499, 600
        transcription_bearings = (True, False)

        if not is_plus_strand:  # swap if we're setting up data for "-" strand
            transcribed_end, transcribed_start = transcribed_start, transcribed_end
            error_end, error_start = error_start, error_end
            transcription_bearings = (False, True)

        self.transcribed_feature = Feature(coordinate=self.coordinates,
                                           is_plus_strand=is_plus_strand,
                                           start=transcribed_start,
                                           end=transcribed_end,
                                           start_is_biological_start=transcription_bearings[0],
                                           end_is_biological_end=transcription_bearings[1],
                                           type=types.TRANSCRIBED)

        self.error_feature = Feature(coordinate=self.coordinates,
                                     is_plus_strand=is_plus_strand,
                                     start=error_start,
                                     end=error_end,
                                     start_is_biological_start=True,
                                     end_is_biological_end=True,
                                     type=types.ERROR)

        self.piece.features = [self.transcribed_feature, self.error_feature]
        sess.add_all([self.ag, self.sl])
        sess.commit()


def test_biointerp_features_as_ranges():
    """checks biological interpretation for ranges from db for simple, spliced, coding gene"""
    sess = mk_memory_session()
    fw = BiointerpDemoDataCoding(sess, is_plus_strand=True)

    assert fw.ti.transcribed_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=100, end=900)
    ]
    assert fw.ti.translated_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=200, end=800)
    ]
    assert fw.ti.intronic_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=300, end=700)
    ]
    assert fw.ti.trans_intronic_ranges() == []

    assert fw.ti.cis_exonic_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=100, end=300),
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=700, end=900)
    ]
    assert fw.ti.translated_exonic_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=200, end=300),
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=700, end=800)
    ]

    assert fw.ti.untranslated_exonic_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=100, end=200),
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=800, end=900)
    ]

    rev = BiointerpDemoDataCoding(sess, is_plus_strand=False)
    assert rev.ti.transcribed_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=900, end=100)
    ]
    assert rev.ti.translated_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=800, end=200)
    ]
    assert rev.ti.intronic_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=700, end=300)
    ]
    assert rev.ti.trans_intronic_ranges() == []

    assert rev.ti.cis_exonic_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=900, end=700),
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=300, end=100)
    ]
    assert rev.ti.translated_exonic_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=800, end=700),
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=300, end=200)
    ]

    assert rev.ti.untranslated_exonic_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=900, end=800),
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=200, end=100)
    ]


def test_biointerp_features_as_ranges_partial():
    """check biological interpretation for ranges from non-coding transcribed fragment"""
    sess = mk_memory_session()
    fw = BiointerpDemoDataPartial(sess, is_plus_strand=True)
    assert fw.ti.transcribed_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=100, end=500)
    ]
    assert fw.ti.error_ranges() == [
        Range(coordinate_id=1, is_plus_strand=True, piece_position=0, start=499, end=600)
    ]
    assert fw.ti.trans_intronic_ranges() == []
    assert fw.ti.translated_ranges() == []

    rev = BiointerpDemoDataPartial(sess, is_plus_strand=False)
    assert rev.ti.transcribed_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=500, end=100)
    ]
    assert rev.ti.error_ranges() == [
        Range(coordinate_id=2, is_plus_strand=False, piece_position=0, start=600, end=499)
    ]


# section: importer
def test_data_frm_gffentry():
    """Test the interpretation and transformation of raw GFF entries."""
    controller = ImportController(database_path='sqlite:///:memory:', err_path=None)

    sess = controller.session
    g, gh = setup_data_handler(importer.GenomeHandler, Genome)
    coords = Coordinate(start=1, end=100000, seqid='NC_015438.2', genome=g)

    sess.add_all([g, coords])
    sess.commit()
    gh.mk_mapper()
    gene_string = 'NC_015438.2\tGnomon\tgene\t4343\t5685\t.\t+\t.\tID=gene0;Dbxref=GeneID:104645797;Name=LOC10'
    mrna_string = 'NC_015438.2\tBestRefSeq\tmRNA\t13024\t15024\t.\t+\t.\tID=rna0;Parent=gene0;Dbxref=GeneID:'
    exon_string = 'NC_015438.2\tGnomon\texon\t4343\t4809\t.\t+\t.\tID=id1;Parent=rna0;Dbxref=GeneID:104645797'
    gene_entry = gffhelper.GFFObject(gene_string)
    controller.clean_entry(gene_entry)
    handler = importer.SuperLocusHandler()
    handler.gffentry = gene_entry
    sl2add = handler.setup_insertion_ready()

    print(gh.gffid_to_coords.keys())
    print(gh._gff_seq_ids)
    assert sl2add['given_name'] == 'gene0'
    assert sl2add['type'] == 'gene'

    mrna_entry = gffhelper.GFFObject(mrna_string)
    mrna_handler = importer.TranscribedHandler()
    mrna2add, piece2add = mrna_handler.setup_insertion_ready(mrna_entry, super_locus=handler)
    piece_handler = mrna_handler.transcribed_piece_handlers[0]
    assert piece2add['transcribed_id'] == mrna2add['id'] == mrna_handler.id

    assert mrna2add['given_name'] == 'rna0'
    assert mrna2add['type'] == 'mRNA'
    assert mrna2add['super_locus_id'] == handler.id

    exon_entry = gffhelper.GFFObject(exon_string)
    controller.clean_entry(exon_entry)
    exon_handler = importer.FeatureHandler()
    exon_handler.gffentry = exon_entry
    exon2add, feature2pieces, feature2translateds = exon_handler.setup_insertion_ready(
        super_locus=handler, transcribed_pieces=[piece_handler], coordinate=coords)

    # accessed via gffentry
    assert exon_handler.gffentry.start == 4343
    assert exon_handler.gffentry.type == 'exon'
    assert 'type' not in exon2add.keys()

    # ready to add
    assert exon2add['is_plus_strand']
    assert exon2add['score'] is None
    assert exon2add['coordinate_id'] == coords.id
    assert exon2add['super_locus_id'] == handler.id
    assert feature2pieces[0] == {
        'transcribed_piece_id': piece_handler.id,
        'feature_id': exon2add['id']
    }
    assert feature2translateds == []


def test_organize_and_split_features():
    """Tests 5'-3' ordering of features and slicing at overlaps"""
    sl, controller = setup_dummyloci_super_locus()
    transcript_full = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    transcript_interpreter = importer.TranscriptInterpreter(transcript_full,
                                                            super_locus=sl,
                                                            controller=controller)

    ordered_features = transcript_interpreter.organize_and_split_features()
    ordered_features = list(ordered_features)
    for i in [0, 4]:
        assert len(ordered_features[i]) == 1
        assert types.CDS not in [x.data.gffentry.type for x in ordered_features[i]]
    for i in [1, 2, 3]:
        assert len(ordered_features[i]) == 2
        assert types.CDS in [x.data.gffentry.type for x in ordered_features[i]]

    transcript_short = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'z'][0]
    transcript_interpreter = importer.TranscriptInterpreter(transcript_short,
                                                            super_locus=sl,
                                                            controller=controller)
    ordered_features = transcript_interpreter.organize_and_split_features()
    ordered_features = list(ordered_features)
    assert len(ordered_features) == 1
    assert len(ordered_features[0]) == 2


def test_possible_types():
    cds = types.OnSequence.CDS.name
    five_prime = types.OnSequence.five_prime_UTR.name
    three_prime = types.OnSequence.three_prime_UTR.name

    sl, controller = setup_dummyloci_super_locus()
    transcript_full = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    transcript_interpreter = importer.TranscriptInterpreter(transcript_full,
                                                            super_locus=sl,
                                                            controller=controller)
    ordered_features = transcript_interpreter.intervals_5to3(plus_strand=True)
    ordered_features = list(ordered_features)
    pt = transcript_interpreter.possible_types(ordered_features[0])
    assert set(pt) == {five_prime, three_prime}
    pt = transcript_interpreter.possible_types(ordered_features[1])
    assert set(pt) == {cds}
    pt = transcript_interpreter.possible_types(ordered_features[-1])
    assert set(pt) == {five_prime, three_prime}


def test_fasta_import():
    """Import and test coordinate information from fasta files"""

    def import_fasta(path):
        controller = ImportController(database_path='sqlite:///:memory:')
        controller.add_sequences(path)
        return controller

    # test import of multiple sequences from one file
    controller = import_fasta('testdata/basic_sequences.fa')
    coords = controller.session.query(Coordinate).all()
    assert len(coords) == 3
    assert coords[0].seqid == '1'
    assert coords[0].start == 0
    assert coords[0].end == 405
    assert coords[0].sha1 == 'dc6f3ba2b0c08f7d08053837b810f86cbaa06f38'
    assert coords[0].sequence == 'N' * 405
    assert coords[1].seqid == 'abc'
    assert coords[1].start == 0
    assert coords[1].end == 808
    assert coords[1].sequence == 'AAGGCCTT' * 101
    assert coords[2].seqid == 'test123'
    assert coords[2].start == 0
    assert coords[2].end == 100
    assert coords[2].sequence == 'A' * 100


def test_import_multiple_genomes():
    controller = ImportController(database_path='sqlite:///:memory:')
    InsertCounterHolder.sync_counters_with_db(controller.session)
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    n_features_1 = query(Feature).count()
    n_coords_1 = query(Coordinate).count()
    n_genomes_1 = query(Genome).count()

    max_id_1 = query(func.max(Feature.id)).one()[0]
    assert max_id_1 == n_features_1

    # add two more genomes
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)

    n_features_3 = query(Feature).count()
    n_coords_3 = query(Coordinate).count()
    n_genomes_3 = query(Genome).count()
    assert n_features_1 * 3 == n_features_3
    assert n_coords_1 * 3 == n_coords_3
    assert n_genomes_1 * 3 == n_genomes_3

    max_id_3 = query(func.max(Feature.id)).one()[0]
    assert max_id_1 * 3 == max_id_3

    # test transcribed and super locus relation for a super loci on the plus strand
    # we have 8 super loci in each genome
    transcribeds_sl_1 = query(SuperLocus).filter(SuperLocus.id == 1).one().transcribeds
    transcribeds_sl_2 = query(SuperLocus).filter(SuperLocus.id == 9).one().transcribeds
    transcribeds_sl_3 = query(SuperLocus).filter(SuperLocus.id == 17).one().transcribeds
    assert len(transcribeds_sl_1) == len(transcribeds_sl_2) == len(transcribeds_sl_3)

    # test translated and super locus relation
    translateds_sl_1 = query(SuperLocus).filter(SuperLocus.id == 1).one().translateds
    translateds_sl_2 = query(SuperLocus).filter(SuperLocus.id == 9).one().translateds
    translateds_sl_3 = query(SuperLocus).filter(SuperLocus.id == 17).one().translateds
    assert len(translateds_sl_1) == len(translateds_sl_2) == len(translateds_sl_3)

    # test transcribed and super locus relation for a super loci on the minus strand
    transcribeds_sl_1 = query(SuperLocus).filter(SuperLocus.id == 6).one().transcribeds
    transcribeds_sl_2 = query(SuperLocus).filter(SuperLocus.id == 14).one().transcribeds
    transcribeds_sl_3 = query(SuperLocus).filter(SuperLocus.id == 22).one().transcribeds
    assert len(transcribeds_sl_1) == len(transcribeds_sl_2) == len(transcribeds_sl_3)

    # test translated and super locus relation
    translateds_sl_1 = query(SuperLocus).filter(SuperLocus.id == 6).one().translateds
    translateds_sl_2 = query(SuperLocus).filter(SuperLocus.id == 14).one().translateds
    translateds_sl_3 = query(SuperLocus).filter(SuperLocus.id == 22).one().translateds
    assert len(translateds_sl_1) == len(translateds_sl_2) == len(translateds_sl_3)


def test_dummyloci_errors():
    """Tests if all errors generated for dummyloci{.gff|.fa} are correct"""

    def error_in_list(error, error_list):
        """searches for the error in a list. removes the error if found.
        error should be a dict and error list a list of orm objects"""
        for e in error_list[:]:  # make a copy at each iteration so we avoid weird errors
            if (error['coord_id'] == e.coordinate.id and error['is_plus_strand'] == e.is_plus_strand
                    and error['start'] == e.start and error['end'] == e.end
                    and error['type'] == e.type.value):
                error_list.remove(e)
                return True
        return False

    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    error_types = [t.value for t in types.Errors]
    errors = controller.session.query(Feature).filter(Feature.type.in_(error_types)).all()
    coords = controller.session.query(Coordinate).all()

    # test case 1 - see gff file for more documentation
    # two identical error bars after cds for aligned exon/cds pair
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 120,
        'end': 499,
        'type': types.MISSING_UTR_3P
    }
    assert error_in_list(error, errors)
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 0,
        'end': 110,
        'type': types.MISSING_UTR_5P
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 0,
        'end': 110,
        'type': types.MISSING_START_CODON
    }
    assert error_in_list(error, errors)

    # test case 2
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 499,
        'end': 1099,
        'type': types.EMPTY_SUPER_LOCUS
    }
    assert error_in_list(error, errors)

    # test case 4 (test case 3 is without errors)
    error = {
        'coord_id': coords[0].id,
        'is_plus_strand': True,
        'start': 1499,
        'end': 1619,
        'type': types.MISSING_START_CODON
    }
    assert error_in_list(error, errors)

    #### Coordinate 1 ####

    # test case 6
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': True,
        'start': 424,
        'end': 524,
        'type': types.WRONG_PHASE_5P
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': True,
        'start': 540,
        'end': 579,
        'type': types.TOO_SHORT_INTRON
    }
    assert error_in_list(error, errors)

    # test case 7
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        'start': 1449,
        'end': 1349,
        'type': types.MISSING_UTR_5P
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        'start': 973,
        'end': -1,
        'type': types.MISMATCHED_PHASE_3P
    }
    assert error_in_list(error, errors)

    # test case 8
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        # inclusive start of the minus strand from beginning of the sequence
        'start': 1755,
        'end': 1724,  # exclusive end of intergenic region
        'type': types.MISSING_START_CODON
    }
    assert error_in_list(error, errors)
    error = {
        'coord_id': coords[1].id,
        'is_plus_strand': False,
        'start': 1649,
        'end': 1548,
        'type': types.OVERLAPPING_EXONS
    }
    assert error_in_list(error, errors)

    # test that we don't have any errors we don't expect
    assert not errors


def test_case_1():
    """Confirm the existence of all features of test case 1 of dummyloci.gff except
    for error features, which are tested in test_dummyloci_errors().
    Does not test the exact ids or exactly matching object relationships."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    sl = query(SuperLocus).filter(SuperLocus.given_name == 'gene0').one()
    sl_h = SuperLocusHandlerBase(sl)

    super_locus = SuperLocus(given_name='gene0', type=types.SuperLocusAll.gene)
    assert sl == super_locus

    # confirm exisistence of all objects where things could go wrong
    # not testing for TranscriptPieces as trans-splicing is currently not implemented
    # above db level and one piece has to exist for the sl_h.features query to work
    sl_objects = list(sl_h.features) + sl_h.data.transcribeds + sl_h.data.translateds

    # first transcript
    transcript = Transcribed(given_name='x', type=types.TranscriptLevelAll.mRNA)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Translated(given_name='x.p')
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='x',
                      type=types.OnSequence.transcribed,
                      start=0,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=10,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=21,
                      end=110,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)

    # second transcript
    transcript = Transcribed(given_name='y', type=types.TranscriptLevelAll.mRNA)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Translated(given_name='y.p')
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='y',
                      type=types.OnSequence.transcribed,
                      start=0,
                      end=400,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=10,
                      end=300,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=22,
                      end=110,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=120,
                      end=200,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)

    # third transcript
    transcript = Transcribed(given_name='z', type=types.TranscriptLevelAll.mRNA)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Translated(given_name='z.p')
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='z',
                      type=types.OnSequence.transcribed,
                      start=110,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=110,
                      end=120,
                      start_is_biological_start=False,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)

    # test if we have no extra objects
    assert not sl_objects


def test_case_8():
    """Analogous to test_case_1()"""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    sl = query(SuperLocus).\
        filter(SuperLocus.given_name == 'gene_overlapping_exons_missing_start').one()
    sl_h = SuperLocusHandlerBase(sl)

    super_locus = SuperLocus(given_name='gene_overlapping_exons_missing_start',
                             type=types.SuperLocusAll.gene)
    assert sl == super_locus

    sl_objects = list(sl_h.features) + sl_h.data.transcribeds + sl_h.data.translateds

    # first transcript
    transcript = Transcribed(given_name='x', type=types.TranscriptLevelAll.mRNA)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Translated(given_name='x.p')
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='x',
                      type=types.OnSequence.transcribed,
                      start=1749,
                      end=1548,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=1724,
                      end=1573,
                      start_is_biological_start=False,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.intron,
                      start=1718,
                      end=1649,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)

    # second transcript
    transcript = Transcribed(given_name='y', type=types.TranscriptLevelAll.mRNA)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Translated(given_name='y.p')
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='y',
                      type=types.OnSequence.transcribed,
                      start=1749,
                      end=1548,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.OnSequence.coding,
                      start=1729,
                      end=1573,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0)
    assert orm_object_in_list(feature, sl_objects)

    # test if we have no extra objects
    assert not sl_objects


def test_gff_gen():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    x = list(gff_organizer._gff_gen())
    assert len(x) == 103
    assert x[0].type == 'region'
    assert x[-1].type == 'CDS'


def test_gff_useful_gen():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    x = list(gff_organizer._useful_gff_entries())
    assert len(x) == 100  # should drop the region entry
    assert x[0].type == 'gene'
    assert x[-1].type == 'CDS'


def test_gff_grouper():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    gff_organizer.load_organized_entries()
    n_genes_seqid = {'NC_015438.2': 2, 'NC_015439.2': 2, 'NC_015440.2': 1}
    for seqid, count in n_genes_seqid.items():
        assert len(gff_organizer.organized_entries[seqid]) == count
        for group in gff_organizer.organized_entries[seqid]:
            assert group[0].type == 'gene'


# section: types
def test_enum_non_inheritance():
    allknown = [x.name for x in list(types.AllKnown)]
    allnice = [x.name for x in list(types.AllKeepable)]
    # check that some random bits made it in to all
    assert 'missing_utr_3p' in allknown
    assert 'region' in allknown

    # check that some annoying bits are not in nice set
    for not_nice in ['transcript', 'primary_transcript', 'exon', 'five_prime_UTR', 'CDS']:
        assert not_nice not in allnice
        assert not_nice in allknown

    # check nothing is there twice
    assert len(set(allknown)) == len(allknown)


def test_enums_name_val_match():
    for x in types.AllKnown:
        assert x.name == x.value


# section: helpers
def test_key_matching():
    # identical
    known = {'a', 'b', 'c'}
    mapper, is_forward = helpers.two_way_key_match(known, known)
    assert is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # subset, should behave as before
    mapper, is_forward = helpers.two_way_key_match(known, {'c'})
    assert is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # superset, should flip ordering
    mapper, is_forward = helpers.two_way_key_match({'c'}, known)
    assert not is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # other is abbreviated from known
    set1 = {'a.seq', 'b.seq', 'c.seq'}
    set2 = {'a', 'b', 'c'}
    mapper, is_forward = helpers.two_way_key_match(set1, set2)
    assert is_forward
    assert mapper('a') == 'a.seq'
    assert isinstance(mapper, helpers.DictMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # cannot be safely differentiated
    set1 = {'ab.seq', 'ba.seq', 'c.seq', 'a.seq', 'b.seq'}
    set2 = {'a', 'b', 'c'}
    with pytest.raises(helpers.NonMatchableIDs):
        mapper, is_forward = helpers.two_way_key_match(set1, set2)
        print(mapper.key_vals)


def test_gff_to_seqids():
    x = helpers.get_seqids_from_gff('testdata/testerSl.gff3')
    assert x == {'NC_015438.2', 'NC_015439.2', 'NC_015440.2'}
