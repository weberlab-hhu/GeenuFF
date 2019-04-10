from dustdas import gffhelper
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError
import sqlalchemy
import pytest
import os

#import sys
#sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from .. import orm
from .. import api
from ..applications import gffimporter
from .. import types
from .. import helpers

### Helper functions ###
# section: orm
def mk_session(db_path='sqlite:///:memory:'):
    engine = create_engine(db_path, echo=False)
    orm.Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()


def setup_data_handler(handler_type, data_type, **kwargs):
    data = data_type(**kwargs)
    handler = handler_type()
    handler.add_data(data)
    return data, handler


def setup_testable_super_locus(db_path='sqlite:///:memory:'):
    controller = gffimporter.ImportControl(err_path='/dev/null', database_path=db_path)
    controller.mk_session()
    controller.add_sequences('testdata/dummyloci.fa')
    controller.add_gff('testdata/dummyloci.gff3', clean=False)
    return controller.super_loci[0], controller


def cleaned_commited_features(sess):
    all_features = sess.query(orm.Feature).all()
    allowed_types = [
        types.TRANSCRIBED,
        types.CODING,
        types.INTRON,
        types.TRANS_INTRON,
        types.ERROR
    ]
    clean_datas = [x for x in all_features if x.type.value in allowed_types]
    return clean_datas


### The actual tests ###
def test_annogenome2coordinate_relation():
    """Check if everything is consistent when we add an Genome and a Coordinate
    to the db. Also check for correct deletion behavior.
    """
    sess = mk_session()
    g = orm.Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coord = orm.Coordinate(start=0, end=30, seqid='abc', genome=g)
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
    sess = mk_session()
    g = orm.Genome()

    # should be ok
    coors = orm.Coordinate(start=0, end=30, seqid='abc', genome=g)
    coors2 = orm.Coordinate(start=0, end=1, seqid='abc', genome=g)
    sess.add_all([coors, coors2])
    sess.commit()

    # should cause trouble
    coors_bad1 = orm.Coordinate(start=-12, end=30, seqid='abc', genome=g)
    coors_bad2 = orm.Coordinate(start=100, end=30, seqid='abc', genome=g)
    coors_bad3 = orm.Coordinate(start=1, end=30, genome=g)
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
    sess = mk_session()
    g = orm.Genome()
    coords = orm.Coordinate(start=1, end=30, seqid='abc', genome=g)
    coords2 = orm.Coordinate(start=11, end=330, seqid='def', genome=g)
    sl = orm.SuperLocus()
    f0 = orm.Feature(coordinate=coords)
    f1 = orm.Feature(coordinate=coords2)
    # should be ok
    sess.add_all([g, sl, coords, coords2, f0, f1])
    assert f0.coordinate.start == 1
    assert f1.coordinate.end == 330


def test_many2many_with_features():
    """Test the many2many tables association_transcribed_piece_to_feature and
    association_translated_to_feature
    """
    sl = orm.SuperLocus()
    # one transcript, multiple proteins
    piece0 = orm.TranscribedPiece()
    scribed0 = orm.Transcribed(super_locus=sl, transcribed_pieces=[piece0])
    slated0 = orm.Translated(super_locus=sl)
    slated1 = orm.Translated(super_locus=sl)
    # features representing alternative start codon for proteins on one transcript
    feat0_tss = orm.Feature(transcribed_pieces=[piece0])
    feat1_tss = orm.Feature(transcribed_pieces=[piece0])
    feat2_stop = orm.Feature(translateds=[slated0, slated1])
    feat3_start = orm.Feature(translateds=[slated0])
    feat4_start = orm.Feature(translateds=[slated1])
    # test multi features per translated worked
    assert len(slated0.features) == 2
    # test mutli translated per feature worked
    assert len(feat2_stop.translateds) == 2
    assert len(feat3_start.translateds) == 1
    assert len(feat0_tss.translateds) == 0


def test_feature_has_its_things():
    """Test if properties of Feature table are correct and constraints are enforced"""
    sess = mk_session()
    # test feature with nothing much set
    g = orm.Genome()
    c = orm.Coordinate(start=1, end=30, seqid='abc', genome=g)
    f = orm.Feature(coordinate=c)
    sess.add_all([f, c])
    sess.commit()

    assert f.is_plus_strand is None
    assert f.source is None
    assert f.score is None
    # test feature with
    f1 = orm.Feature(is_plus_strand=False, start=3, end=-1, coordinate=c)
    assert not f1.is_plus_strand
    assert f1.start == 3
    assert f1.end == -1
    sess.add(f1)
    sess.commit()
    # test too low of start / end coordinates raise an error
    f_should_fail = orm.Feature(is_plus_strand=True, start=-5, end=10, coordinate=c)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()

    f_should_fail = orm.Feature(is_plus_strand=False, start=5, end=-2, coordinate=c)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()
    # test wrong class as parameter
    with pytest.raises(KeyError):
        f2 = orm.Feature(coordinate=f)

    f2 = orm.Feature(is_plus_strand=-1)  # note that 0, and 1 are accepted
    sess.add(f2)
    with pytest.raises(sqlalchemy.exc.StatementError):
        sess.commit()
    sess.rollback()


def test_partially_remove_coordinate():
    """Add two coordinates to an annotated genome and test if everything ends
    up valid.
    """
    sess = mk_session()
    g = orm.Genome()
    place_holder = orm.Genome()
    coord0 = orm.Coordinate(start=1, end=30, seqid='abc', genome=g)
    coord1 = orm.Coordinate(start=11, end=330, seqid='def', genome=g)
    sess.add_all([g, coord0, coord1])
    sess.commit()
    assert len(g.coordinates) == 2
    g.coordinates.remove(coord0)
    coord0.genome = place_holder  # else we'll fail the not NULL constraint
    sess.commit()
    # removed from g
    assert len(g.coordinates) == 1
    # but still in table
    assert len(sess.query(orm.Coordinate).all()) == 2

# section: api
def test_transcribed_piece_unique_constraints():
    """Add transcribed pieces in valid and invalid configurations and test for
    valid outcomes.
    """
    sess = mk_session()
    sl = orm.SuperLocus()
    transcribed0 = orm.Transcribed(super_locus=sl)
    transcribed1 = orm.Transcribed(super_locus=sl)

    # test if same position for different transcribed_id goes through
    piece_tr0_pos0 = orm.TranscribedPiece(transcribed=transcribed0, position=0)
    piece_tr1_pos0 = orm.TranscribedPiece(transcribed=transcribed1, position=0)
    sess.add_all([transcribed0, piece_tr0_pos0, piece_tr1_pos0])
    sess.commit()

    # same transcibed_id but different position
    piece_tr0_pos1 = orm.TranscribedPiece(transcribed=transcribed0, position=1)
    sess.add(piece_tr0_pos1)
    sess.commit()

    # test if unique constraint works
    piece_tr0_pos1_2nd = orm.TranscribedPiece(transcribed=transcribed0, position=1)
    sess.add(piece_tr0_pos1_2nd)
    with pytest.raises(IntegrityError):
        sess.commit()


def test_order_pieces():
    """Add transcribed pieces that consists of one or many features to the db and test
    if the order for the transcribed pieces and the features is returned according to the
    position property instead of db insertion order.
    """
    sess = mk_session()
    g = orm.Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coor = orm.Coordinate(seqid='a', start=1, end=1000, genome=g)
    sess.add_all([g, coor])
    sess.commit()
    # setup one transcribed handler with pieces
    sl, slh = setup_data_handler(api.SuperLocusHandlerBase, orm.SuperLocus)
    scribed, scribedh = setup_data_handler(api.TranscribedHandlerBase, orm.Transcribed, super_locus=sl)
    ti = api.TranscriptInterpBase(transcript=scribedh, super_locus=sl, session=sess)
    # wrong order
    piece1 = orm.TranscribedPiece(position=1)
    piece0 = orm.TranscribedPiece(position=0)
    piece2 = orm.TranscribedPiece(position=2)
    scribed.transcribed_pieces = [piece0, piece1, piece2]
    sess.add_all([scribed, piece1, piece0, piece2])
    sess.commit()
    # setup some features on different pieces
    feature0u = orm.Feature(transcribed_pieces=[piece0],
                            coordinate=coor,
                            start=100,
                            end=200,
                            given_id='0u',
                            is_plus_strand=True)
    feature1d = orm.Feature(transcribed_pieces=[piece1],
                            coordinate=coor,
                            start=1,
                            end=3,
                            given_id='1d',
                            is_plus_strand=True)
    feature1u = orm.Feature(transcribed_pieces=[piece1],
                            coordinate=coor,
                            start=100,
                            end=200,
                            given_id='1u',
                            is_plus_strand=True)
    feature2d = orm.Feature(transcribed_pieces=[piece2],
                            coordinate=coor,
                            start=1,
                            end=3,
                            given_id='2d',
                            is_plus_strand=True)
    # see if they can be ordered as expected overall
    op = ti.sort_pieces()
    print([piece0, piece1, piece2], 'expected')
    print(op, 'sorted')
    assert op == [piece0, piece1, piece2]
    # and finally, order features by piece
    fully_sorted = ti.sort_all()
    expected = [[feature0u],
                [feature1d, feature1u],
                [feature2d]]
    assert fully_sorted == expected


class TransspliceDemoData(object):
    """Setup of a rather complex transplicing scenario."""
    def __init__(self, sess):
        # setup two transitions:
        # 1) scribed [-> ABC][-> FED]
        #  transcribed A [  (               )  ]   D [  (                 ) ]
        #       coding B [          (       )* ]   E [ *(             )     ]
        # trans_intron C [               (  )* ]   F [ *(     )             ]

        # 2) scribedflip [-> ABC][<- FpEpDp]
        #  transcribed A [  (               )  ]   Dp[ (                 )  ]
        #       coding B [          (       )* ]   Ep[     (             )* ]
        # trans_intron C [               (  )* ]   Fp[             (     )* ]

        g = orm.Genome()
        self.old_coor = orm.Coordinate(seqid='a', start=1, end=2000, genome=g)
        self.sl, self.slh = setup_data_handler(api.SuperLocusHandlerBase, orm.SuperLocus)
        self.scribed, self.scribedh = setup_data_handler(api.TranscribedHandlerBase,
                                                         orm.Transcribed,
                                                         super_locus=self.sl)
        self.scribedflip, self.scribedfliph = setup_data_handler(api.TranscribedHandlerBase,
                                                                 orm.Transcribed,
                                                                 super_locus=self.sl)
        self.ti = api.TranscriptInterpBase(transcript=self.scribedh,
                                           super_locus=self.sl,
                                           session=sess)
        self.tiflip = api.TranscriptInterpBase(transcript=self.scribedfliph,
                                               super_locus=self.sl,
                                               session=sess)

        self.pieceA2C = orm.TranscribedPiece(position=0, transcribed=self.scribed)
        self.pieceD2F = orm.TranscribedPiece(position=1, transcribed=self.scribed)
        self.pieceA2C_prime = orm.TranscribedPiece(position=0, transcribed=self.scribedflip)
        self.pieceD2F_prime = orm.TranscribedPiece(position=1, transcribed=self.scribedflip)
        self.scribed.transcribed_pieces = [self.pieceA2C, self.pieceD2F]
        self.scribedflip.transcribed_pieces = [self.pieceA2C_prime, self.pieceD2F_prime]
        # pieceA2C features
        self.fA = orm.Feature(coordinate=self.old_coor, start=10, end=40, given_id='A',
                              start_is_biological_start=True, end_is_biological_end=True,
                              is_plus_strand=True, type=types.TRANSCRIBED)

        self.fB = orm.Feature(coordinate=self.old_coor, start=20, end=40, given_id='B',
                              start_is_biological_start=True, end_is_biological_end=False,
                              is_plus_strand=True, type=types.CODING)

        self.fC = orm.Feature(coordinate=self.old_coor, start=30, end=40, given_id='C',
                              start_is_biological_start=True, end_is_biological_end=False,
                              is_plus_strand=True, type=types.TRANS_INTRON)

        # pieceD2F features
        self.fD = orm.Feature(coordinate=self.old_coor, start=910, end=940, given_id='D',
                              start_is_biological_start=True, end_is_biological_end=True,
                              is_plus_strand=True, type=types.TRANSCRIBED)

        self.fE = orm.Feature(coordinate=self.old_coor, start=910, end=930, given_id='E',
                              start_is_biological_start=False, end_is_biological_end=True,
                              is_plus_strand=True, type=types.CODING)

        self.fF = orm.Feature(coordinate=self.old_coor, start=910, end=920, given_id='F',
                              start_is_biological_start=False, end_is_biological_end=True,
                              is_plus_strand=True, type=types.TRANS_INTRON)

        # pieceD2F_prime features
        self.fDp = orm.Feature(coordinate=self.old_coor, start=940, end=910, given_id='Dp',
                               start_is_biological_start=True, end_is_biological_end=True,
                               is_plus_strand=False, type=types.TRANSCRIBED)
        self.fEp = orm.Feature(coordinate=self.old_coor, start=940, end=920, given_id='Ep',
                               start_is_biological_start=False, end_is_biological_end=True,
                               is_plus_strand=False, type=types.CODING)
        self.fFp = orm.Feature(coordinate=self.old_coor, start=940, end=930, given_id='Fp',
                               start_is_biological_start=False, end_is_biological_end=True,
                               is_plus_strand=False, type=types.TRANS_INTRON)

        self.pieceA2C.features = [self.fA, self.fB, self.fC]
        self.pieceA2C_prime.features = [self.fA, self.fB, self.fC]
        self.pieceD2F.features = [self.fD, self.fE, self.fF]
        self.pieceD2F_prime.features = [self.fDp, self.fEp, self.fFp]
        sess.add(self.sl)
        sess.commit()

    def make_all_handlers(self):
        self.slh.make_all_handlers()


def test_transition_transsplice():
    """Test if the transsplicing scenario of TransspliceDemoData is interpreted correctly by
    transition_5p_to_3p() of TranscriptInterpBase.
    """
    sess = mk_session()
    d = TransspliceDemoData(sess)  # setup _d_ata
    d.make_all_handlers()
    # forward pass, same sequence, two pieces
    ti_transitions = list(d.ti.transition_5p_to_3p())
    # expect (start, end) of features to be sorted 5'-3' within pieces
    # from transition gen: 0 -> aligned_features, 1 -> piece
    assert [set(x[0]) for x in ti_transitions] == [{d.fA}, {d.fB}, {d.fC}, {d.fF}, {d.fE}, {d.fD}]

    # forward, then backward pass, same sequence, two pieces
    ti_transitions = list(d.tiflip.transition_5p_to_3p())
    assert [set(x[0]) for x in ti_transitions] == [{d.fA}, {d.fB}, {d.fC}, {d.fFp}, {d.fEp}, {d.fDp}]


class BiointerpDemoDataCoding(object):
    def __init__(self, sess, is_plus_strand):
        self.sl, self.sl_handler = setup_data_handler(api.SuperLocusHandlerBase, orm.SuperLocus)
        self.scribed, self.scribed_handler = setup_data_handler(
            api.TranscribedHandlerBase, orm.Transcribed, super_locus=self.sl)

        self.piece = orm.TranscribedPiece(position=0, transcribed=self.scribed)

        self.ti = api.TranscriptInterpBase(transcript=self.scribed_handler, super_locus=self.sl, session=sess)

        self.g = orm.Genome()
        # setup ranges for a two-exon coding gene
        self.coordinate = orm.Coordinate(seqid='a', start=1, end=2000, genome=self.g)
        transcribed_start, transcribed_end = 100, 900
        coding_start, coding_end = 200, 800
        intron_start, intron_end = 300, 700

        if not is_plus_strand:  # swap if we're setting up data for "-" strand
            transcribed_end, transcribed_start = transcribed_start, transcribed_end
            coding_end, coding_start = coding_start, coding_end
            intron_end, intron_start = intron_start, intron_end

        # transcribed:
        self.transcribed_feature = orm.Feature(coordinate=self.coordinate,
                                               is_plus_strand=is_plus_strand, start=transcribed_start,
                                               end=transcribed_end, start_is_biological_start=True,
                                               end_is_biological_end=True, type=types.TRANSCRIBED)
        # coding:
        self.coding_feature = orm.Feature(coordinate=self.coordinate,
                                          is_plus_strand=is_plus_strand, start=coding_start,
                                          end=coding_end, start_is_biological_start=True,
                                          end_is_biological_end=True, type=types.CODING)
        # intron:
        self.intron_feature = orm.Feature(coordinate=self.coordinate,
                                          is_plus_strand=is_plus_strand, start=intron_start,
                                          end=intron_end, start_is_biological_start=True,
                                          end_is_biological_end=True, type=types.INTRON)

        self.piece.features = [self.transcribed_feature,
                               self.coding_feature,
                               self.intron_feature]

        sess.add(self.sl)
        sess.commit()


class BioInterpDemodataTranssplice(object):
    def __init__(self, sess, is_plus_strand_piece0, is_plus_strand_piece1):
        self.sl, self.sl_handler = setup_data_handler(api.SuperLocusHandlerBase, orm.SuperLocus)
        self.scribed, self.scribed_handler = setup_data_handler(
            api.TranscribedHandlerBase, orm.Transcribed, super_locus=self.sl)

        self.piece0 = orm.TranscribedPiece(position=0, transcribed=self.scribed)
        self.piece1 = orm.TranscribedPiece(position=1, transcribed=self.scribed)
        self.ti = api.TranscriptInterpBase(transcript=self.scribed_handler, super_locus=self.sl, session=sess)

        # setup ranges for a two-exon coding gene
        self.ag = orm.Genome()
        self.coordinates = orm.Coordinate(seqid='a', start=1, end=2000, genome=self.ag)

        t0_start, t0_end = 500, 700
        t1_start, t1_end = 100, 300

        if not is_plus_strand_piece0:
            t0_start, t0_end = t0_end, t0_start
        if not is_plus_strand_piece1:
            t1_start, t1_end = t1_end, t1_start

        trans_donor_splice, trans_intron_close = 600, t0_end
        trans_acceptor_splice, trans_intron_open = 200, t1_start

        # piece 0
        self.transcribed0_feature = orm.Feature(coordinate=self.coordinates,
                                                is_plus_strand=is_plus_strand_piece0, start=t0_start,
                                                end=t0_end, start_is_biological_start=True,
                                                end_is_biological_end=True, type=types.TRANSCRIBED)

        self.trans_intron0_feature = orm.Feature(coordinate=self.coordinates,
                                                 is_plus_strand=is_plus_strand_piece0, start=trans_donor_splice,
                                                 end=trans_intron_close, start_is_biological_start=True,
                                                 end_is_biological_end=False, type=types.TRANS_INTRON)

        # piece 1
        self.transcribed1_feature = orm.Feature(coordinate=self.coordinates,
                                                is_plus_strand=is_plus_strand_piece1, start=t1_start,
                                                end=t1_end, start_is_biological_start=True,
                                                end_is_biological_end=True, type=types.TRANSCRIBED)

        self.trans_intron1_feature = orm.Feature(coordinate=self.coordinates,
                                                 is_plus_strand=is_plus_strand_piece1, start=trans_intron_open,
                                                 end=trans_acceptor_splice, start_is_biological_start=False,
                                                 end_is_biological_end=True, type=types.TRANS_INTRON)

        self.piece0.features = [self.transcribed0_feature,
                                self.trans_intron0_feature]
        self.piece1.features = [self.transcribed1_feature,
                                self.trans_intron1_feature]

        sess.add_all([self.ag, self.sl])
        sess.commit()


class BiointerpDemoDataPartial(object):
    def __init__(self, sess, is_plus_strand):
        self.sl, self.sl_handler = setup_data_handler(api.SuperLocusHandlerBase, orm.SuperLocus)
        self.scribed, self.scribed_handler = setup_data_handler(
            api.TranscribedHandlerBase, orm.Transcribed, super_locus=self.sl)

        self.piece = orm.TranscribedPiece(position=0, transcribed=self.scribed)

        self.ti = api.TranscriptInterpBase(transcript=self.scribed_handler, super_locus=self.sl, session=sess)

        self.ag = orm.Genome()
        # setup ranges for a two-exon coding gene
        self.coordinates = orm.Coordinate(seqid='a', start=1, end=2000, genome=self.ag)
        transcribed_start, transcribed_end = 100, 500
        error_start, error_end = 499, 600
        transcription_bearings = (True, False)

        if not is_plus_strand:  # swap if we're setting up data for "-" strand
            transcribed_end, transcribed_start = transcribed_start, transcribed_end
            error_end, error_start = error_start, error_end
            transcription_bearings = (False, True)

        self.transcribed_feature = orm.Feature(coordinate=self.coordinates,
                                               is_plus_strand=is_plus_strand, start=transcribed_start,
                                               end=transcribed_end, start_is_biological_start=transcription_bearings[0],
                                               end_is_biological_end=transcription_bearings[1],
                                               type=types.TRANSCRIBED)

        self.error_feature = orm.Feature(coordinate=self.coordinates,
                                         is_plus_strand=is_plus_strand, start=error_start,
                                         end=error_end, start_is_biological_start=True,
                                         end_is_biological_end=True, type=types.ERROR)

        self.piece.features = [self.transcribed_feature, self.error_feature]
        sess.add_all([self.ag, self.sl])
        sess.commit()


def test_biointerp_features_as_ranges():
    """checks biological interpretation for ranges from db for simple, spliced, coding gene"""
    sess = mk_session()
    fw = BiointerpDemoDataCoding(sess, is_plus_strand=True)

    assert fw.ti.transcribed_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                    start=100, end=900)]
    assert fw.ti.translated_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                   start=200, end=800)]
    assert fw.ti.intronic_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                 start=300, end=700)]
    assert fw.ti.trans_intronic_ranges() == []

    assert fw.ti.cis_exonic_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                   start=100, end=300),
                                         api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                   start=700, end=900)]
    assert fw.ti.translated_exonic_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                          start=200, end=300),
                                                api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                          start=700, end=800)]

    assert fw.ti.untranslated_exonic_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                            start=100, end=200),
                                                  api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                            start=800, end=900)]

    rev = BiointerpDemoDataCoding(sess, is_plus_strand=False)
    assert rev.ti.transcribed_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                     start=900, end=100)]
    assert rev.ti.translated_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                    start=800, end=200)]
    assert rev.ti.intronic_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                  start=700, end=300)]
    assert rev.ti.trans_intronic_ranges() == []

    assert rev.ti.cis_exonic_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                    start=900, end=700),
                                          api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                    start=300, end=100)]
    assert rev.ti.translated_exonic_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                           start=800, end=700),
                                                 api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                           start=300, end=200)]

    assert rev.ti.untranslated_exonic_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                             start=900, end=800),
                                                   api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                             start=200, end=100)]


def test_biointerp_features_as_ranges_transsplice():
    """checks biological interpretation for ranges from db for non-coding, trans-spliced (from same chromosome) gene"""
    sess = mk_session()
    fw = BioInterpDemodataTranssplice(sess, is_plus_strand_piece0=True, is_plus_strand_piece1=True)

    assert fw.ti.transcribed_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                    start=500, end=700),
                                          api.Range(coordinate_id=1, is_plus_strand=True, piece_position=1,
                                                    start=100, end=300)]

    assert fw.ti.trans_intronic_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                       start=600, end=700),
                                             api.Range(coordinate_id=1, is_plus_strand=True, piece_position=1,
                                                       start=100, end=200)]

    rev = BioInterpDemodataTranssplice(sess, is_plus_strand_piece0=False, is_plus_strand_piece1=True)
    assert rev.ti.transcribed_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                     start=700, end=500),
                                           api.Range(coordinate_id=2, is_plus_strand=True, piece_position=1,
                                                     start=100, end=300)]
    assert rev.ti.trans_intronic_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                        start=600, end=500),
                                              api.Range(coordinate_id=2, is_plus_strand=True, piece_position=1,
                                                        start=100, end=200)]


def test_biointerp_features_as_ranges_partial():
    """check biological interpretation for ranges from non-coding transcribed fragment"""
    sess = mk_session()
    fw = BiointerpDemoDataPartial(sess, is_plus_strand=True)
    assert fw.ti.transcribed_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                                    start=100, end=500)]
    assert fw.ti.error_ranges() == [api.Range(coordinate_id=1, is_plus_strand=True, piece_position=0,
                                              start=499, end=600)]
    assert fw.ti.trans_intronic_ranges() == []
    assert fw.ti.translated_ranges() == []

    rev = BiointerpDemoDataPartial(sess, is_plus_strand=False)
    assert rev.ti.transcribed_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                                     start=500, end=100)]
    assert rev.ti.error_ranges() == [api.Range(coordinate_id=2, is_plus_strand=False, piece_position=0,
                                               start=600, end=499)]


def test_biointerp_features_as_transitions():
    """checks biological interpretation for sites from db for simple, spliced, coding gene"""
    sess = mk_session()
    fw = BiointerpDemoDataCoding(sess, is_plus_strand=True)
    assert fw.ti.transcription_start_sites() == [api.TranscriptCoordinate(coordinate_id=1, is_plus_strand=True,
                                                                          piece_position=0, start=100)]
    assert fw.ti.transcription_end_sites() == [api.TranscriptCoordinate(coordinate_id=1, is_plus_strand=True,
                                                                        piece_position=0, start=900)]
    assert fw.ti.translation_start_sites() == [api.TranscriptCoordinate(coordinate_id=1, is_plus_strand=True,
                                                                        piece_position=0, start=200)]
    assert fw.ti.translation_end_sites() == [api.TranscriptCoordinate(coordinate_id=1, is_plus_strand=True,
                                                                      piece_position=0, start=800)]
    assert fw.ti.intron_start_sites() == [api.TranscriptCoordinate(coordinate_id=1, is_plus_strand=True,
                                                                   piece_position=0, start=300)]
    assert fw.ti.intron_end_sites() == [api.TranscriptCoordinate(coordinate_id=1, is_plus_strand=True,
                                                                 piece_position=0, start=700)]
    # and minus strand
    rev = BiointerpDemoDataCoding(sess, is_plus_strand=False)
    assert rev.ti.transcription_start_sites() == [api.TranscriptCoordinate(coordinate_id=2, is_plus_strand=False,
                                                                           piece_position=0, start=900)]
    assert rev.ti.transcription_end_sites() == [api.TranscriptCoordinate(coordinate_id=2, is_plus_strand=False,
                                                                         piece_position=0, start=100)]
    assert rev.ti.translation_start_sites() == [api.TranscriptCoordinate(coordinate_id=2, is_plus_strand=False,
                                                                         piece_position=0, start=800)]
    assert rev.ti.translation_end_sites() == [api.TranscriptCoordinate(coordinate_id=2, is_plus_strand=False,
                                                                       piece_position=0, start=200)]
    assert rev.ti.intron_start_sites() == [api.TranscriptCoordinate(coordinate_id=2, is_plus_strand=False,
                                                                    piece_position=0, start=700)]
    assert rev.ti.intron_end_sites() == [api.TranscriptCoordinate(coordinate_id=2, is_plus_strand=False,
                                                                  piece_position=0, start=300)]

    # trans-splicing
    trans = BioInterpDemodataTranssplice(sess, is_plus_strand_piece0=True, is_plus_strand_piece1=False)
    assert trans.ti.transcription_start_sites() == [api.TranscriptCoordinate(coordinate_id=3, is_plus_strand=True,
                                                                             piece_position=0, start=500),
                                                    api.TranscriptCoordinate(coordinate_id=3, is_plus_strand=False,
                                                                             piece_position=1, start=300)]
    assert trans.ti.transcription_end_sites() == [api.TranscriptCoordinate(coordinate_id=3, is_plus_strand=True,
                                                                           piece_position=0, start=700),
                                                  api.TranscriptCoordinate(coordinate_id=3, is_plus_strand=False,
                                                                           piece_position=1, start=100)]

    assert trans.ti.trans_intron_start_sites() == [api.TranscriptCoordinate(coordinate_id=3, is_plus_strand=True,
                                                                            piece_position=0, start=600)]
    assert trans.ti.trans_intron_end_sites() == [api.TranscriptCoordinate(coordinate_id=3, is_plus_strand=False,
                                                                          piece_position=1, start=200)]

    # partial
    partial = BiointerpDemoDataPartial(sess, is_plus_strand=True)
    assert partial.ti.transcription_start_sites() == [api.TranscriptCoordinate(coordinate_id=4, is_plus_strand=True,
                                                                               piece_position=0, start=100)]
    assert partial.ti.transcription_end_sites() == []

    partial_rev = BiointerpDemoDataPartial(sess, is_plus_strand=False)
    assert partial_rev.ti.transcription_end_sites() == [api.TranscriptCoordinate(coordinate_id=5, is_plus_strand=False,
                                                                                 piece_position=0, start=100)]
    assert partial_rev.ti.transcription_start_sites() == []


# section: gffimporter
def test_data_frm_gffentry():
    """Test the interpretation and transformation of raw GFF entries."""
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:', err_path=None)

    sess = controller.session
    g, gh = setup_data_handler(gffimporter.GenomeHandler, orm.Genome)
    coords = orm.Coordinate(start=1, end=100000, seqid='NC_015438.2', genome=g)

    sess.add_all([g, coords])
    sess.commit()
    gh.mk_mapper()  # todo, why doesn't this work WAS HERE
    gene_string = 'NC_015438.2\tGnomon\tgene\t4343\t5685\t.\t+\t.\tID=gene0;Dbxref=GeneID:104645797;Name=LOC10'
    mrna_string = 'NC_015438.2\tBestRefSeq\tmRNA\t13024\t15024\t.\t+\t.\tID=rna0;Parent=gene0;Dbxref=GeneID:'
    exon_string = 'NC_015438.2\tGnomon\texon\t4343\t4809\t.\t+\t.\tID=id1;Parent=rna0;Dbxref=GeneID:104645797'
    gene_entry = gffhelper.GFFObject(gene_string)
    controller.clean_entry(gene_entry)
    handler = gffimporter.SuperLocusHandler()
    handler.gffentry = gene_entry
    sl2add = handler.setup_insertion_ready()

    print(gh.gffid_to_coords.keys())
    print(gh._gff_seq_ids)
    assert sl2add['given_id'] == 'gene0'
    assert sl2add['type'] == 'gene'

    mrna_entry = gffhelper.GFFObject(mrna_string)
    mrna_handler = gffimporter.TranscribedHandler()
    mrna2add, piece2add = mrna_handler.setup_insertion_ready(mrna_entry, super_locus=handler)
    piece_handler = mrna_handler.transcribed_piece_handlers[0]
    assert piece2add['transcribed_id'] == mrna2add['id'] == mrna_handler.id

    assert mrna2add['given_id'] == 'rna0'
    assert mrna2add['type'] == 'mRNA'
    assert mrna2add['super_locus_id'] == handler.id

    exon_entry = gffhelper.GFFObject(exon_string)
    controller.clean_entry(exon_entry)
    exon_handler = gffimporter.FeatureHandler()
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
    assert feature2pieces[0] == {'transcribed_piece_id': piece_handler.id, 'feature_id': exon2add['id']}
    assert feature2translateds == []


def test_organize_and_split_features():
    """Tests 5'-3' ordering of features and slicing at overlaps"""
    sl, controller = setup_testable_super_locus()
    transcript_full = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    transcript_interpreter = gffimporter.TranscriptInterpreter(transcript_full,
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
    transcript_interpreter = gffimporter.TranscriptInterpreter(transcript_short,
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

    sl, controller = setup_testable_super_locus()
    transcript_full = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    transcript_interpreter = gffimporter.TranscriptInterpreter(transcript_full,
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


def test_import_coordinate():
    """Import and test coordinate information from a dummy gff file."""
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    controller.mk_session()
    seq_path = 'testdata/dummyloci.fa'
    controller.add_sequences(seq_path)
    coors = controller.genome_handler.data.coordinates
    assert len(coors) == 1
    assert coors[0].seqid == '1'
    assert coors[0].start == 0
    assert coors[0].end == 405
    assert coors[0].sha1 == 'dc6f3ba2b0c08f7d08053837b810f86cbaa06f38'  # sha1 for 'N' * 405


def test_transcript_interpreter():
    """tests decoding of raw from-gff features via the TranscriptInterpreter class"""
    sl, controller = setup_testable_super_locus()
    transcript_handler = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    # change so that there are implicit UTRs
    t_interp = gffimporter.TranscriptInterpreter(transcript_handler,
                                                 super_locus=sl,
                                                 controller=controller)
    t_interp.decode_raw_features()
    # has all standard features
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    types_out = set([x.type.value for x in features])
    assert types_out == {types.CODING,
                         types.TRANSCRIBED,
                         types.INTRON}
    bearings_out = set([x.start_is_biological_start for x in features] +
                       [x.end_is_biological_end for x in features])

    assert bearings_out == {True}

    transcribeds = [x for x in features if x.type.value == types.TRANSCRIBED]
    assert len(transcribeds) == 1
    assert transcribeds[0].end == 400
    assert transcribeds[0].start == 0


def test_transcript_get_first():
    # plus strand
    sl, controller = setup_testable_super_locus()
    transcript_handler = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript_handler,
                                                 super_locus=sl,
                                                 controller=controller)
    sorted_intervals = t_interp.intervals_5to3(plus_strand=True)
    i0 = sorted_intervals[0]
    for intv in sorted_intervals:
        print(intv)
    t_interp.interpret_first_pos(i0)
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    status = t_interp.status
    assert len(features) == 1
    f0 = features[0]
    assert f0.start == 0
    assert status.is_5p_utr()
    assert f0.phase is None
    assert f0.is_plus_strand


def test_transcript_get_first_minus_strand():
    # minus strand
    sl, controller = setup_testable_super_locus()
    for transcript_handler in sl.transcribed_handlers:
        for feature in transcript_handler.feature_handlers:
            feature.gffentry.strand = '-'
            feature.add_shortcuts_from_gffentry()

    transcript_handler = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]

    t_interp = gffimporter.TranscriptInterpreter(transcript_handler,
                                                 super_locus=sl,
                                                 controller=controller)
    sorted_intervals = t_interp.intervals_5to3(plus_strand=False)

    i0 = sorted_intervals[0]
    t_interp.interpret_first_pos(i0, plus_strand=False)
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    print(features)
    features = [f for f in features if not f.is_plus_strand]
    status = t_interp.status
    assert len(features) == 1
    f0 = features[0]
    print(f0)
    print(status)
    print(f0.type)
    assert f0.start == 399
    assert status.is_5p_utr()
    assert f0.phase is None


def test_transcript_get_first_without_UTR():
    # minus strand
    sl, controller = setup_testable_super_locus()
    for transcript_handler in sl.transcribed_handlers:
        for feature in transcript_handler.feature_handlers:
            feature.gffentry.strand = '-'
            feature.add_shortcuts_from_gffentry()

    transcript_handler = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'x'][0]

    t_interp = gffimporter.TranscriptInterpreter(transcript_handler,
                                                 super_locus=sl,
                                                 controller=controller)
    # test without UTR (x doesn't have last exon, and therefore will end in CDS); remember, flipped it to minus strand
    prefeatures = transcript_handler.feature_handlers
    print('prefeatures', prefeatures)
    print([x.gffentry for x in prefeatures])
    i0 = t_interp.intervals_5to3(plus_strand=False)[0]
    t_interp.interpret_first_pos(i0, plus_strand=False)
    controller.session.commit()
    controller.execute_so_far()
    print('??=\n', controller.session.query(orm.Feature).all())
    features = cleaned_commited_features(controller.session)
    print('features ====>\n', features)
    status = t_interp.status
    assert len(features) == 3
    f_err = [f for f in features if f.type.value == types.ERROR][0]
    f_coding = [f for f in features if f.type.value == types.CODING][0]
    f_transcribed = [f for f in features if f.type.value == types.TRANSCRIBED][0]

    # should get status instead of a start codon and tss
    assert f_coding.start == 119
    assert not f_coding.start_is_biological_start
    assert f_transcribed.start == 119
    assert not f_transcribed.start_is_biological_start

    assert not f_coding.is_plus_strand
    # region beyond exon should be marked erroneous
    assert not f_err.is_plus_strand
    assert f_err.start == 404
    assert f_err.end == 118  # so that err overlaps 1bp with the coding status checked above
    assert status.is_coding()
    assert status.coding_tracker.seen_start
    assert status.transcribed_tracker.in_region


def test_transcript_transition_from_5p_to_end():
    sl, controller = setup_testable_super_locus()
    transcript_handler = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript_handler,
                                                 super_locus=sl,
                                                 controller=controller)
    ivals_sets = t_interp.intervals_5to3(plus_strand=True)
    t_interp.interpret_first_pos(ivals_sets[0])
    # hit start codon
    t_interp.interpret_transition(ivals_before=ivals_sets[0],
                                  ivals_after=ivals_sets[1],
                                  plus_strand=True)

    features = controller.features_to_add
    coding_feature = features[-1]
    assert coding_feature["type"] == types.CODING
    assert coding_feature["start_is_biological_start"]
    assert coding_feature["start"] == 10
    assert coding_feature["end"] is None
    # hit splice site
    t_interp.interpret_transition(ivals_before=ivals_sets[1],
                                  ivals_after=ivals_sets[2],
                                  plus_strand=True)

    intron_feature0 = features[-1]
    assert intron_feature0["type"] == types.INTRON
    assert intron_feature0["end"] == 110
    assert intron_feature0["start"] == 100  # splice from
    assert t_interp.status.is_coding()
    # hit splice site
    t_interp.interpret_transition(ivals_before=ivals_sets[2],
                                  ivals_after=ivals_sets[3],
                                  plus_strand=True)

    intron_feature1 = features[-1]
    assert intron_feature1["type"] == types.INTRON
    assert intron_feature1["start"] == 120  # splice from
    assert intron_feature1["end"] == 200  # splice to
    # hit stop codon
    t_interp.interpret_transition(ivals_before=ivals_sets[3],
                                  ivals_after=ivals_sets[4],
                                  plus_strand=True)

    # end should be added now
    assert coding_feature["end"] == 300
    assert coding_feature["end_is_biological_end"]
    # hit transcription termination site
    t_interp.interpret_last_pos(ivals_sets[4], plus_strand=True)
    controller.session.commit()
    controller.execute_so_far()
    # spot check transcribed_feature after database entry
    fin_features = cleaned_commited_features(controller.session)
    transcribed_features = [x for x in fin_features if x.type.value == types.TRANSCRIBED]
    assert len(transcribed_features) == 1
    assert transcribed_features[0].start == 0
    assert transcribed_features[0].end == 400
    assert transcribed_features[0].start_is_biological_start
    assert transcribed_features[0].end_is_biological_end


def test_non_coding_transitions():
    sl, controller = setup_testable_super_locus()
    transcript_handler = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'z'][0]
    piece = transcript_handler.one_piece()
    t_interp = gffimporter.TranscriptInterpreter(transcript_handler,
                                                 super_locus=sl,
                                                 controller=controller)
    # get single-exon no-CDS transcript
    # delete CDS feature
    features = transcript_handler.feature_handlers
    for i in range(len(features) - 1, -1, -1):
        print(i)
        if features[i].gffentry.type == types.CDS:
            del features[i]
    # check we just setup TSS then TTS in that order
    ivals_sets = t_interp.intervals_5to3(plus_strand=True)
    assert len(ivals_sets) == 1
    t_interp.interpret_first_pos(ivals_sets[0])

    features = controller.features_to_add
    transcribed_feature = features[-1]
    assert transcribed_feature["type"] == types.TRANSCRIBED
    assert transcribed_feature["start_is_biological_start"]
    assert transcribed_feature["start"] == 110
    assert transcribed_feature["end"] is None
    assert transcribed_feature["end_is_biological_end"] is None

    t_interp.interpret_last_pos(ivals_sets[0], plus_strand=True)

    assert transcribed_feature["start"] == 110  # no change expected
    assert transcribed_feature["end"] == 120
    assert transcribed_feature["end_is_biological_end"]
    assert len(features) == 1


def test_errors_not_lost():
    sl, controller = setup_testable_super_locus()
    s = "1\tGnomon\tgene\t20\t405\t0.\t-\t0\tID=eg_missing_children"
    gene_entry = gffhelper.GFFObject(s)

    coordinate = controller.genome_handler.data.coordinates[0]

    sl._mark_erroneous(gene_entry, coordinate=coordinate, controller=controller)
    assert len(sl.transcribed_handlers) == 4

    sl.check_and_fix_structure(coordinate=coordinate, controller=controller)
    features = controller.session.query(orm.Feature).filter(orm.Feature.given_id == 'eg_missing_children').all()
    assert len(features) == 1
    assert features[0].start == 404
    assert features[0].end == 18


def test_setup_proteins():
    """Tests explicit protein creation from mRNA & CDS gff features"""
    sl, controller = setup_testable_super_locus()
    transcript = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript, sl, controller)
    controller.execute_so_far()
    print(t_interp.proteins)
    assert len(t_interp.proteins.keys()) == 1

    transcript_orm = controller.session.query(orm.Transcribed).filter(orm.Transcribed.given_id == 'y').first()
    protein_orm = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'y.p').first()
    assert protein_orm.given_id == 'y.p'


def test_cp_features_to_prot():
    """Test linking of appropriate features to newly created protein (translated)"""
    sl, controller = setup_testable_super_locus()
    controller.session.commit()
    transcript_handler = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript_handler, sl, controller)

    t_interp.decode_raw_features()
    controller.session.commit()
    controller.execute_so_far()
    # grab protein again jic
    protein = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'y.p').all()
    assert len(protein) == 1
    protein = protein[0]
    print(protein.id)
    print(controller.session.query(orm.association_translated_to_feature).all())
    assert len(protein.features) == 1
    p_feature = protein.features[0]
    assert p_feature.type.value == types.CODING
    assert p_feature.start_is_biological_start == p_feature.end_is_biological_end == True


def test_check_and_fix_structure():
    """checks expected 'geenuff' format can be produced from the initial pre-entries parsed from the gff file"""
    # so we save a copy of the cleaned up loci once per test run
    rel_path = 'testdata/dummyloci_annotations.sqlitedb'
    db_path = 'sqlite:///{}'.format(rel_path)
    if os.path.exists(rel_path):
        os.remove(rel_path)
    sl, controller = setup_testable_super_locus(db_path)
    coordinate = controller.genome_handler.data.coordinates[0]

    sl.check_and_fix_structure(coordinate=coordinate, controller=controller)
    controller.execute_so_far()
    # check handling of nice transcript
    protein = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'y.p').all()
    print(controller.session.query(orm.association_translated_to_feature).all())
    assert len(protein) == 1
    protein = protein[0]
    print(protein.id)
    print(controller.session.query(orm.association_translated_to_feature).all())
    # check we get a coding region/feature for the nice transcript
    assert len(protein.features) == 1
    p_feature = protein.features[0]
    assert p_feature.type.value == types.CODING
    assert p_feature.start == 10
    assert p_feature.end == 300
    assert p_feature.start_is_biological_start
    assert p_feature.end_is_biological_end

    # check we get a transcript with transcribed, coding and two intronic regions
    piece = controller.session.query(orm.TranscribedPiece).filter(orm.TranscribedPiece.given_id == 'y').first()
    print(piece)
    assert len(piece.features) == 4
    assert set([x.type.value for x in piece.features]) == {types.TRANSCRIBED,
                                                           types.INTRON,
                                                           types.CODING,
                                                           }
    assert set([(x.start_is_biological_start, x.end_is_biological_end) for x in piece.features]) == {(True, True)}
    # check handling of truncated transcript
    piece = controller.session.query(orm.TranscribedPiece).filter(orm.TranscribedPiece.given_id == 'x').first()
    protein = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'x.p').first()
    print(protein.features)
    assert len(protein.features) == 1
    p_feature = protein.features[0]
    assert p_feature.type.value == types.CODING
    assert p_feature.start == 10
    assert p_feature.end == 120
    assert p_feature.start_is_biological_start
    assert not p_feature.end_is_biological_end

    assert len(piece.features) == 4
    assert set([x.type.value for x in piece.features]) == {types.TRANSCRIBED, types.INTRON,
                                                           types.ERROR, types.CODING}
    coding_fs = [x for x in piece.features if x.type.value == types.CODING]
    assert len(coding_fs) == 1
    assert coding_fs[0].start_is_biological_start
    assert not coding_fs[0].end_is_biological_end

    transcribed_fs = [x for x in piece.features if x.type.value == types.TRANSCRIBED]
    assert len(transcribed_fs) == 1
    assert transcribed_fs[0].start_is_biological_start
    assert not transcribed_fs[0].end_is_biological_end
    assert transcribed_fs[0].start == 0
    assert transcribed_fs[0].end == 120

    sl_datas = controller.session.query(orm.SuperLocus).all()
    assert len(sl_datas) == 1
    assert len(sl_datas[0].translateds) == 3


def test_erroneous_splice():
    db_path = 'sqlite:///:memory:'

    sl, controller = setup_testable_super_locus(db_path)
    sess = controller.session
    # get target transcript
    transcript = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'x'][0]
    # fish out "first exon" features and extend so intron is of -length
    features = []
    for scribed in sl.transcribed_handlers:
        features += scribed.feature_handlers
    f0_handler = [f for f in features if f.gffentry.get_ID() == 'ftr000000'][0]
    f1_handler = [f for f in features if f.gffentry.get_ID() == 'ftr000001'][0]
    #f0 = sess.query(orm.Feature).filter(orm.Feature.given_id == 'ftr000000').first()
    #f1 = sess.query(orm.Feature).filter(orm.Feature.given_id == 'ftr000001').first()
    f0_handler.gffentry.end = f1_handler.gffentry.end = 115

    ti = gffimporter.TranscriptInterpreter(transcript, controller=controller, super_locus=sl)
    ti.decode_raw_features()
    sess.commit()
    controller.execute_so_far()
    clean_datas = cleaned_commited_features(sess)
    # TSS, start codon, 2x error splice, 2x error splice, 2x error no stop
    print('---\n'.join([str(x) for x in clean_datas]))
    assert len(clean_datas) == 5

    assert len([x for x in clean_datas if x.type.value == types.ERROR]) == 3
    # make sure splice error covers whole exon-intron-exon region
    print('clean datas here::::')
    print('\n\n'.join([str(x) for x in clean_datas]))
    assert clean_datas[2].type.value == types.ERROR
    assert clean_datas[2].start_is_biological_start
    assert clean_datas[2].start == 10
    assert clean_datas[2].end_is_biological_end
    assert clean_datas[2].end == 120
    sess.commit()


def test_gff_gen():
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    x = list(controller.gff_gen('testdata/testerSl.gff3'))
    assert len(x) == 103
    assert x[0].type == 'region'
    assert x[-1].type == 'CDS'


def test_gff_useful_gen():
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    x = list(controller.useful_gff_entries('testdata/testerSl.gff3'))
    assert len(x) == 100  # should drop the region entry
    assert x[0].type == 'gene'
    assert x[-1].type == 'CDS'


def test_gff_grouper():
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    x = list(controller.group_gff_by_gene('testdata/testerSl.gff3'))
    assert len(x) == 5
    for group in x:
        assert group[0].type == 'gene'


def test_import2biointerp():

    class TCShortcut(api.TranscriptCoordinate):
        """just fills in what's always the same for the TranscriptCoordinate to stop retyping"""
        def __init__(self, start):
            super().__init__(coordinate_id=1, piece_position=0, is_plus_strand=True, start=start)

    class RShortcut(api.Range):
        """just fills in what's always the same for the Range to stop retyping"""
        def __init__(self, start, end):
            super().__init__(coordinate_id=1, piece_position=0, is_plus_strand=True, start=start, end=end)

    sequence_path = 'testdata/biointerp_loci.fa'
    gff3 = 'testdata/biointerp_loci.gff3'
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:', err_path='/dev/null')
    controller.add_sequences(sequence_path)
    controller.add_gff(gff3)
    controller.session.commit()

    super_loci = controller.session.query(orm.SuperLocus).all()
    assert len(super_loci) == 2

    ### coding ###
    sl_coding = [sl for sl in super_loci if sl.given_id == "gene0"][0]
    sl_coding_handler = api.SuperLocusHandlerBase()
    sl_coding_handler.add_data(sl_coding)
    assert len(sl_coding.transcribeds) == 1
    scribed_coding = sl_coding.transcribeds[0]
    scribed_coding_handler = api.TranscribedHandlerBase()
    scribed_coding_handler.add_data(scribed_coding)

    ti = api.TranscriptInterpBase(transcript=scribed_coding_handler, super_locus=sl_coding_handler,
                                  session=controller.session)

    # check transcribed sites and range
    assert ti.transcription_start_sites() == [TCShortcut(start=100)]

    assert ti.transcription_end_sites() == [TCShortcut(start=900)]
    assert ti.transcribed_ranges() == [RShortcut(start=100, end=900)]
    # check coding sites and ranges
    assert ti.translation_start_sites() == [TCShortcut(start=200)]
    assert ti.translation_end_sites() == [TCShortcut(start=800)]
    assert ti.translated_ranges() == [RShortcut(start=200, end=800)]

    # check intron sites and ranges
    assert ti.intron_start_sites() == [TCShortcut(start=300)]
    assert ti.intron_end_sites() == [TCShortcut(start=700)]
    assert ti.intronic_ranges() == [RShortcut(start=300, end=700)]
    # check exons and similar
    assert ti.cis_exonic_ranges() == [RShortcut(start=100, end=300),
                                      RShortcut(start=700, end=900)]
    assert ti.translated_exonic_ranges() == [RShortcut(start=200, end=300),
                                             RShortcut(start=700, end=800)]
    assert ti.untranslated_exonic_ranges() == [RShortcut(start=100, end=200),
                                               RShortcut(start=800, end=900)]

    # check that a few interpretations that should be empty are
    assert ti.trans_intron_end_sites() == []
    assert ti.error_ranges() == []

    ### partial ###
    sl_partial = [sl for sl in super_loci if sl.given_id == "gene1"][0]
    sl_partial_handler = api.SuperLocusHandlerBase()
    sl_partial_handler.add_data(sl_coding)
    assert len(sl_partial.transcribeds) == 1
    scribed_partial = sl_partial.transcribeds[0]
    scribed_partial_handler = api.TranscribedHandlerBase()
    scribed_partial_handler.add_data(scribed_partial)

    ti = api.TranscriptInterpBase(transcript=scribed_partial_handler, super_locus=sl_partial_handler,
                                  session=controller.session)

    assert ti.transcription_start_sites() == [TCShortcut(start=100)]
    assert ti.transcription_end_sites() == []
    assert ti.transcribed_ranges() == [RShortcut(start=100, end=500)]

    assert ti.translation_start_sites() == [TCShortcut(start=200)]
    assert ti.translation_end_sites() == []
    assert ti.translated_ranges() == [RShortcut(start=200, end=500)]

    err_ranges = ti.error_ranges()
    assert len(err_ranges) == 1
    assert err_ranges[0].start == 499


# section: types
def test_enum_non_inheritance():
    allknown = [x.name for x in list(types.AllKnown)]
    allnice = [x.name for x in list(types.AllKeepable)]
    # check that some random bits made it in to all
    assert 'error' in allknown
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

