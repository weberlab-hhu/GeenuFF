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


def setup_testable_super_loci(db_path='sqlite:///:memory:'):
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
    """Check if everything is consistent when we add an AnnotatedGenome and a Coordinate
    to the db. Also check for correct deletion behavior.
    """
    sess = mk_session()
    ag = orm.AnnotatedGenome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coord = orm.Coordinates(start=0, end=30, seqid='abc', annotated_genome=ag)
    assert ag is coord.annotated_genome
    # actually put everything in db
    sess.add(coord)
    sess.commit()
    # check primary keys were assigned
    assert ag.id == 1
    assert coord.id == 1
    # check we can access coordinates from ag
    coord_q = ag.coordinates[0]
    assert coord is coord_q
    assert ag is coord.annotated_genome
    assert ag.id == coord_q.annotated_genome_id
    # check we get logical behavior on deletion
    sess.delete(coord)
    sess.commit()
    assert len(ag.coordinates) == 0
    print(coord.annotated_genome)
    sess.delete(ag)
    sess.commit()
    with pytest.raises(sqlalchemy.exc.InvalidRequestError):
        sess.add(coord)


def test_coordinate_constraints():
    """Check the coordinate constraints"""
    sess = mk_session()
    ag = orm.AnnotatedGenome()

    # should be ok
    coors = orm.Coordinates(start=0, end=30, seqid='abc', annotated_genome=ag)
    coors2 = orm.Coordinates(start=0, end=1, seqid='abc', annotated_genome=ag)
    sess.add_all([coors, coors2])
    sess.commit()

    # should cause trouble
    coors_bad1 = orm.Coordinates(start=-12, end=30, seqid='abc', annotated_genome=ag)
    coors_bad2 = orm.Coordinates(start=100, end=30, seqid='abc', annotated_genome=ag)
    coors_bad3 = orm.Coordinates(start=1, end=30, annotated_genome=ag)
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


def test_coordinates_insert():
    """Test what happens when we insert two coordinates"""
    sess = mk_session()
    ag = orm.AnnotatedGenome()
    coords = orm.Coordinates(start=1, end=30, seqid='abc', annotated_genome=ag)
    coords2 = orm.Coordinates(start=11, end=330, seqid='def', annotated_genome=ag)
    sl = orm.SuperLocus()
    f0 = orm.Feature(coordinates=coords)
    f1 = orm.Feature(coordinates=coords2)
    # should be ok
    sess.add_all([ag, sl, coords, coords2, f0, f1])
    assert f0.coordinates.start == 1
    assert f1.coordinates.end == 330
 


def test_many2many_with_features():
    """Test the many2many tables association_transcribed_pieces_to_features and
    association_translateds_to_features
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
    ag = orm.AnnotatedGenome()
    c = orm.Coordinates(start=1, end=30, seqid='abc', annotated_genome=ag)
    f = orm.Feature(coordinates=c)
    sess.add(f)
    sess.commit()

    assert f.is_plus_strand is None
    assert f.source is None
    assert f.score is None
    # test feature with
    f1 = orm.Feature(is_plus_strand=False, position=3)
    assert not f1.is_plus_strand
    assert f1.position == 3
    # test wrong class as parameter
    with pytest.raises(KeyError):
        f2 = orm.Feature(coordinates=f)

    f2 = orm.Feature(is_plus_strand=-1)  # note that 0, and 1 are accepted
    sess.add(f2)
    with pytest.raises(sqlalchemy.exc.StatementError):
        sess.commit()
    sess.rollback()


def test_partially_remove_coordinates():
    """Add two coordinates to an annotated genome and test if everything ends
    up valid.
    """
    sess = mk_session()
    ag = orm.AnnotatedGenome()
    place_holder = orm.AnnotatedGenome()
    coord0 = orm.Coordinates(start=1, end=30, seqid='abc', annotated_genome=ag)
    coord1 = orm.Coordinates(start=11, end=330, seqid='def', annotated_genome=ag)
    sess.add_all([ag, coord0, coord1])
    sess.commit()
    assert len(ag.coordinates) == 2
    ag.coordinates.remove(coord0)
    coord0.annotated_genome = place_holder  # else we'll fail the not NULL constraint
    sess.commit()
    # removed from ag
    assert len(ag.coordinates) == 1
    # but still in table
    assert len(sess.query(orm.Coordinates).all()) == 2

# section: api
def test_transcribed_pieces_unique_constraints():
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
    ag = orm.AnnotatedGenome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coor = orm.Coordinates(seqid='a', start=1, end=1000, annotated_genome=ag)
    sess.add_all([ag, coor])
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
    # setup some paired features
    feature0u = orm.Feature(transcribed_pieces=[piece0],
                            coordinates=coor,
                            position=100,
                            given_id='0u',
                            is_plus_strand=True)
    feature1d = orm.Feature(transcribed_pieces=[piece1],
                            coordinates=coor,
                            position=1,
                            given_id='1d',
                            is_plus_strand=True)
    feature1u = orm.Feature(transcribed_pieces=[piece1],
                            coordinates=coor,
                            position=100,
                            given_id='1u',
                            is_plus_strand=True)
    feature2d = orm.Feature(transcribed_pieces=[piece2],
                            coordinates=coor,
                            position=1,
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
        # 1) scribed - [->[TSS(A),START(B),TDSS(C{->F}),TTS(D)], ->[TSS(E), <<slice>>> TASS(F),STOP(G),TTS(H)]]
        # 2) scribedflip - [->[TSS(A),START(B),TDSS(C{->F'}),TTS(D)], <-[TTS(H'), <<slice>> STOP(G'),TASS(F'),TSS(E')]]
        ag = orm.AnnotatedGenome()
        self.old_coor = orm.Coordinates(seqid='a', start=1, end=2000, annotated_genome=ag)
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

        self.pieceA2D = orm.TranscribedPiece(position=0, transcribed=self.scribed)
        self.pieceE2H = orm.TranscribedPiece(position=1, transcribed=self.scribed)
        self.pieceA2Dp = orm.TranscribedPiece(position=0, transcribed=self.scribedflip)
        self.pieceEp2Hp = orm.TranscribedPiece(position=1, transcribed=self.scribedflip)
        self.scribed.transcribed_pieces = [self.pieceA2D, self.pieceE2H]
        self.scribedflip.transcribed_pieces = [self.pieceA2Dp, self.pieceEp2Hp]
        # pieceA2D features

        self.fA = orm.Feature(coordinates=self.old_coor, position=10, given_id='A',
                              is_plus_strand=True, type=types.TRANSCRIBED, bearing=types.START)
        self.fB = orm.Feature(coordinates=self.old_coor, position=20, given_id='B',
                              is_plus_strand=True, type=types.CODING, bearing=types.START)

        self.fC = orm.Feature(coordinates=self.old_coor, position=30, given_id='C',
                              is_plus_strand=True, type=types.TRANS_INTRON, bearing=types.START)
        self.fD = orm.Feature(coordinates=self.old_coor, position=40, given_id='D',
                              is_plus_strand=True, type=types.TRANSCRIBED, bearing=types.END)
        self.fADs0 = orm.Feature(coordinates=self.old_coor, position=40, given_id='ADs0',
                                 is_plus_strand=True, type=types.TRANS_INTRON, bearing=types.CLOSE_STATUS)
        self.fADs1 = orm.Feature(coordinates=self.old_coor, position=40, given_id='ADs1',
                                 is_plus_strand=True, type=types.CODING, bearing=types.CLOSE_STATUS)
        # pieceE2H features
        self.fEHs0 = orm.Feature(coordinates=self.old_coor, position=910, given_id='EHs0',
                                 is_plus_strand=True, type=types.TRANS_INTRON, bearing=types.OPEN_STATUS)
        self.fEHs1 = orm.Feature(coordinates=self.old_coor, position=910, given_id='EHs1',
                                 is_plus_strand=True, type=types.CODING, bearing=types.OPEN_STATUS)
        self.fE = orm.Feature(coordinates=self.old_coor, position=910, given_id='E',
                              is_plus_strand=True, type=types.TRANSCRIBED, bearing=types.START)
        self.fF = orm.Feature(coordinates=self.old_coor, position=920, given_id='F',
                              is_plus_strand=True, type=types.TRANS_INTRON, bearing=types.END)
        self.fG = orm.Feature(coordinates=self.old_coor, position=930, given_id='G',
                              is_plus_strand=True, type=types.CODING, bearing=types.END)
        self.fH = orm.Feature(coordinates=self.old_coor, position=940, given_id='H',
                              is_plus_strand=True, type=types.TRANSCRIBED, bearing=types.END)
        # pieceEp2Hp features
        self.fEHps0 = orm.Feature(coordinates=self.old_coor, position=940, given_id='EHsp0',
                                  is_plus_strand=False, type=types.TRANS_INTRON, bearing=types.OPEN_STATUS)
        self.fEHps1 = orm.Feature(coordinates=self.old_coor, position=940, given_id='EHsp1',
                                  is_plus_strand=False, type=types.CODING, bearing=types.OPEN_STATUS)
        self.fEp = orm.Feature(coordinates=self.old_coor, position=940, given_id='Ep',
                               is_plus_strand=False, type=types.TRANSCRIBED, bearing=types.START)
        self.fFp = orm.Feature(coordinates=self.old_coor, position=930, given_id='Fp',
                               bearing=types.END, is_plus_strand=False, type=types.TRANS_INTRON)
        self.fGp = orm.Feature(coordinates=self.old_coor, position=920, given_id='Gp',
                               is_plus_strand=False, type=types.CODING, bearing=types.END)
        self.fHp = orm.Feature(coordinates=self.old_coor, position=910, given_id='Hp',
                               is_plus_strand=False, type=types.TRANSCRIBED, bearing=types.END)

        self.pieceA2D.features = [self.fA, self.fB, self.fC, self.fD, self.fADs0, self.fADs1]
        self.pieceA2Dp.features = [self.fA, self.fB, self.fC, self.fD, self.fADs0, self.fADs1]
        self.pieceE2H.features = [self.fE, self.fF, self.fG, self.fH, self.fEHs0, self.fEHs1]
        self.pieceEp2Hp.features = [self.fEp, self.fFp, self.fGp, self.fHp, self.fEHps0, self.fEHps1]
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
    # from transition gen: 0 -> aligned_Features, 1 -> status copy
    assert [set(x[0]) for x in ti_transitions] == [{d.fA}, {d.fB}, {d.fC}, {d.fD}, {d.fADs1, d.fADs0},
                                                   {d.fEHs0, d.fEHs1}, {d.fE}, {d.fF}, {d.fG}, {d.fH}]
    print([x[1].genic for x in ti_transitions])
    assert [x[1].genic for x in ti_transitions] == [bool(x) for x in [1, 1, 1, 0, 0, 0, 1, 1, 1, 0]]
    assert [x[1].in_translated_region for x in ti_transitions] == [bool(x) for x in [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]]
    assert [x[1].in_trans_intron for x in ti_transitions] == [bool(x) for x in [0, 0, 1, 1, 0, 1, 1, 0, 0, 0]]
    # forward, then backward pass, same sequence, two pieces
    ti_transitions = list(d.tiflip.transition_5p_to_3p())
    assert [set(x[0]) for x in ti_transitions] == [{d.fA}, {d.fB}, {d.fC}, {d.fD}, {d.fADs1, d.fADs0},
                                                   {d.fEHps0, d.fEHps1}, {d.fEp}, {d.fFp}, {d.fGp}, {d.fHp}]
    print([x[1].genic for x in ti_transitions])
    assert [x[1].genic for x in ti_transitions] == [bool(x) for x in [1, 1, 1, 0, 0, 0, 1, 1, 1, 0]]
    assert [x[1].in_translated_region for x in ti_transitions] == [bool(x) for x in [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]]
    assert [x[1].in_trans_intron for x in ti_transitions] == [bool(x) for x in [0, 0, 1, 1, 0, 1, 1, 0, 0, 0]]


class BiointerpDemoDataCoding(object):
    def __init__(self, sess, is_plus_strand):
        self.sl, self.sl_handler = setup_data_handler(api.SuperLocusHandlerBase, orm.SuperLocus)
        self.scribed, self.scribed_handler = setup_data_handler(
            api.TranscribedHandlerBase, orm.Transcribed, super_locus=self.sl)

        self.piece = orm.TranscribedPiece(position=0, transcribed=self.scribed)

        self.ti = api.TranscriptInterpBase(transcript=self.scribed_handler, super_locus=self.sl, session=sess)

        self.ag = orm.AnnotatedGenome()
        # setup ranges for a two-exon coding gene
        self.coordinates = orm.Coordinates(seqid='a', start=1, end=2000, annotated_genome=self.ag)
        transcribed_start, transcribed_end = 100, 900
        coding_start, coding_end = 200, 800
        intron_start, intron_end = 300, 700

        if not is_plus_strand:  # swap if we're setting up data for "-" strand
            transcribed_end, transcribed_start = transcribed_start, transcribed_end
            coding_end, coding_start = coding_start, coding_end
            intron_end, intron_start = intron_start, intron_end

        # transcribed:
        self.transcribed_start = orm.Feature(coordinates=self.coordinates,
                                             is_plus_strand=is_plus_strand, position=transcribed_start,
                                             type=types.TRANSCRIBED, bearing=types.START)
        self.transcribed_end = orm.Feature(coordinates=self.coordinates,
                                           is_plus_strand=is_plus_strand, position=transcribed_end,
                                           type=types.TRANSCRIBED, bearing=types.END)
        # coding:
        self.coding_start = orm.Feature(coordinates=self.coordinates,
                                        is_plus_strand=is_plus_strand, position=coding_start,
                                        type=types.CODING, bearing=types.START)
        self.coding_end = orm.Feature(coordinates=self.coordinates,
                                      is_plus_strand=is_plus_strand, position=coding_end,
                                      type=types.CODING, bearing=types.END)
        # intron:
        self.intron_start = orm.Feature(coordinates=self.coordinates,
                                        is_plus_strand=is_plus_strand, position=intron_start,
                                        type=types.INTRON, bearing=types.START)
        self.intron_end = orm.Feature(coordinates=self.coordinates,
                                      is_plus_strand=is_plus_strand, position=intron_end,
                                      type=types.INTRON, bearing=types.END)

        self.piece.features = [self.transcribed_start, self.transcribed_end,
                               self.coding_start, self.coding_end,
                               self.intron_start, self.intron_end]

        sess.add(self.sl)
        sess.commit()


class BioInterpDemodataTranssplice(object):
    def __init__(self, sess, is_plus_strand_piece0, is_plus_strand_piece1):
        self.sl, self.sl_handler = setup_data_handler(api.SuperLocusHandlerBase, orm.SuperLocus)
        self.scribed, self.scribed_handler = setup_data_handler(
            api.TranscribedHandlerBase, orm.Transcribed, super_locus=self.sl)

        self.piece = orm.TranscribedPiece(position=0, transcribed=self.scribed)
        self.piece = orm.TranscribedPiece(position=1, transcribed=self.scribed)
        self.ti = api.TranscriptInterpBase(transcript=self.scribed_handler, super_locus=self.sl, session=sess)

        # setup ranges for a two-exon coding gene
        self.coordinates = orm.Coordinates(seqid='a', start=1, end=2000)
        # todo, setup 2x pieces on the same sequence/coord for sort testing fun


def test_biointerp_features_as_ranges():
    sess = mk_session()
    fw = BiointerpDemoDataCoding(sess, is_plus_strand=True)
    coord = "coordinate_id"
    is_plus = "is_plus_strand"
    start = "start"
    end = "end"
    pos = "piece_position"

    assert fw.ti.transcribed_ranges() == [{coord: 1, is_plus: True, pos: 0,
                                           start: 100, end: 900}]
    assert fw.ti.translated_ranges() == [{coord: 1, is_plus: True, pos: 0,
                                          start: 200, end: 800}]
    assert fw.ti.intronic_ranges() == [{coord: 1, is_plus: True, pos: 0,
                                       start: 300, end: 700}]
    assert fw.ti.trans_intronic_ranges() == []

    assert fw.ti.cis_exonic_ranges() == [{coord: 1, is_plus: True, pos: 0,
                                         start: 100, end: 300},
                                         {coord: 1, is_plus: True, pos: 0,
                                         start: 700, end: 900}]
    assert fw.ti.translated_exonic_ranges() == [{coord: 1, is_plus: True, pos: 0,
                                                 start: 200, end: 300},
                                                {coord: 1, is_plus: True, pos: 0,
                                                 start: 700, end: 800}]

    assert fw.ti.untranslated_exonic_ranges() == [{coord: 1, is_plus: True, pos: 0,
                                                   start: 100, end: 200},
                                                  {coord: 1, is_plus: True, pos: 0,
                                                   start: 800, end: 900}]

    rev = BiointerpDemoDataCoding(sess, is_plus_strand=False)
    assert rev.ti.transcribed_ranges() == [{coord: 2, is_plus: False, pos: 0,
                                            start: 900, end: 100}]
    assert rev.ti.translated_ranges() == [{coord: 2, is_plus: False, pos: 0,
                                          start: 800, end: 200}]
    assert rev.ti.intronic_ranges() == [{coord: 2, is_plus: False, pos: 0,
                                         start: 700, end: 300}]
    assert rev.ti.trans_intronic_ranges() == []

    assert rev.ti.cis_exonic_ranges() == [{coord: 2, is_plus: False, pos: 0,
                                          start: 900, end: 700},
                                          {coord: 2, is_plus: False, pos: 0,
                                          start: 300, end: 100}]
    assert rev.ti.translated_exonic_ranges() == [{coord: 2, is_plus: False, pos: 0,
                                                  start: 800, end: 700},
                                                 {coord: 2, is_plus: False, pos: 0,
                                                  start: 300, end: 200}]

    assert rev.ti.untranslated_exonic_ranges() == [{coord: 2, is_plus: False, pos: 0,
                                                    start: 900, end: 800},
                                                   {coord: 2, is_plus: False, pos: 0,
                                                    start: 200, end: 100}]


def test_biointerp_features_as_transitions():
    pass


# section: gffimporter
def test_data_frm_gffentry():
    """Test the interpretation and transformation of raw GFF entries."""
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:', err_path=None)

    sess = controller.session
    ag = orm.AnnotatedGenome()
    ag, agh = setup_data_handler(gffimporter.AnnotatedGenomeHandler, orm.AnnotatedGenome)
    coords = orm.Coordinates(start=1, end=100000, seqid='NC_015438.2', annotated_genome=ag)

    sess.add_all([ag, coords])
    sess.commit()
    agh.mk_mapper()  # todo, why doesn't this work WAS HERE
    gene_string = 'NC_015438.2\tGnomon\tgene\t4343\t5685\t.\t+\t.\tID=gene0;Dbxref=GeneID:104645797;Name=LOC10'
    mrna_string = 'NC_015438.2\tBestRefSeq\tmRNA\t13024\t15024\t.\t+\t.\tID=rna0;Parent=gene0;Dbxref=GeneID:'
    exon_string = 'NC_015438.2\tGnomon\texon\t4343\t4809\t.\t+\t.\tID=id1;Parent=rna0;Dbxref=GeneID:104645797'
    gene_entry = gffhelper.GFFObject(gene_string)
    controller.clean_entry(gene_entry)
    handler = gffimporter.SuperLocusHandler()
    handler.gffentry = gene_entry
    sl2add = handler.setup_insertion_ready()

    print(agh.gffid_to_coords.keys())
    print(agh._gff_seq_ids)
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
        super_locus=handler, transcribed_pieces=[piece_handler], coordinates=coords)

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
    sl, controller = setup_testable_super_loci()
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

    sl, controller = setup_testable_super_loci()
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


def test_import_coordinates():
    """Import and test coordinate information from a dummy gff file."""
    controller = gffimporter.ImportControl(database_path='sqlite:///:memory:')
    controller.mk_session()
    seq_path = 'testdata/dummyloci.fa'
    controller.add_sequences(seq_path)
    coors = controller.anno_genome_handler.data.coordinates
    assert len(coors) == 1
    assert coors[0].seqid == '1'
    assert coors[0].start == 0
    assert coors[0].end == 405
    assert coors[0].sha1 == 'dc6f3ba2b0c08f7d08053837b810f86cbaa06f38'  # sha1 for 'N' * 405


def test_transcript_interpreter():
    sl, controller = setup_testable_super_loci()
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
    bearings_out = set([x.bearing.value for x in features])
    assert bearings_out == {types.START, types.END}

    transcribeds = [x for x in features if x.type.value == types.TRANSCRIBED]
    assert max([x.position for x in transcribeds]) == 400
    assert min([x.position for x in transcribeds]) == 0


def test_transcript_get_first():
    # plus strand
    sl, controller = setup_testable_super_loci()
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
    assert f0.position == 0
    assert status.is_5p_utr()
    assert f0.phase is None
    assert f0.is_plus_strand


def test_transcript_get_first_minus_strand():
    # minus strand
    sl, controller = setup_testable_super_loci()
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
    assert f0.position == 399
    assert status.is_5p_utr()
    assert f0.phase is None


def test_transcript_get_first_without_UTR():
    # minus strand
    sl, controller = setup_testable_super_loci()
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
    features = cleaned_commited_features(controller.session)
    status = t_interp.status
    assert len(features) == 4
    f_err_open = [f for f in features
                  if f.bearing.value == types.START and f.type.value == types.ERROR][0]
    f_err_close = [f for f in features
                   if f.bearing.value == types.END and f.type.value == types.ERROR][0]
    f_status_coding = [f for f in features
                       if f.bearing.value == types.OPEN_STATUS and f.type.value == types.CODING][0]
    f_status_transcribed = [f for f in features
                            if f.bearing.value == types.OPEN_STATUS and f.type.value == types.TRANSCRIBED][0]
    print(f_err_open)
    print(status)
    print(i0)
    print([f.is_plus_strand for f in features], 'plus strands')
    # should get status instead of a start codon and tss
    assert f_status_coding.position == 119
    assert f_status_transcribed.position == 119
    assert not f_status_coding.is_plus_strand
    # region beyond exon should be marked erroneous
    assert not f_err_close.is_plus_strand and not f_err_open.is_plus_strand
    assert f_err_close.position == 118  # so that err overlaps 1bp with the coding status checked above
    assert f_err_open.position == 404
    assert status.is_coding()
    assert status.seen_start
    assert status.genic


def test_transcript_transition_from_5p_to_end():
    sl, controller = setup_testable_super_loci()
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
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    assert features[-1].type.value == types.CODING
    assert features[-1].bearing.value == types.START
    assert features[-1].position == 10
    # hit splice site
    t_interp.interpret_transition(ivals_before=ivals_sets[1],
                                  ivals_after=ivals_sets[2],
                                  plus_strand=True)
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    assert features[-1].type.value == types.INTRON
    assert features[-1].bearing.value == types.END
    assert features[-2].type.value == types.INTRON
    assert features[-2].bearing.value == types.START
    assert features[-2].position == 100  # splice from
    assert features[-1].position == 110  # splice to
    assert t_interp.status.is_coding()
    # hit splice site
    t_interp.interpret_transition(ivals_before=ivals_sets[2],
                                  ivals_after=ivals_sets[3],
                                  plus_strand=True)
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    assert features[-1].type.value == types.INTRON
    assert features[-1].bearing.value == types.END
    assert features[-2].type.value == types.INTRON
    assert features[-2].bearing.value == types.START
    assert features[-2].position == 120  # splice from
    assert features[-1].position == 200  # splice to
    # hit stop codon
    t_interp.interpret_transition(ivals_before=ivals_sets[3],
                                  ivals_after=ivals_sets[4],
                                  plus_strand=True)
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    assert features[-1].type.value == types.CODING
    assert features[-1].bearing.value == types.END
    assert features[-1].position == 300
    # hit transcription termination site
    t_interp.interpret_last_pos(ivals_sets[4], plus_strand=True)
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    assert features[-1].type.value == types.TRANSCRIBED
    assert features[-1].bearing.value == types.END
    assert features[-1].position == 400


def test_non_coding_transitions():
    sl, controller = setup_testable_super_loci()
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
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    assert features[-1].type.value == types.TRANSCRIBED
    assert features[-1].bearing.value == types.START
    assert features[-1].position == 110
    t_interp.interpret_last_pos(ivals_sets[0], plus_strand=True)
    controller.session.commit()
    controller.execute_so_far()
    features = cleaned_commited_features(controller.session)
    assert features[-1].type.value == types.TRANSCRIBED
    assert features[-1].bearing.value == types.END
    assert features[-1].position == 120
    assert len(features) == 2


def test_errors_not_lost():
    sl, controller = setup_testable_super_loci()
    s = "1\tGnomon\tgene\t20\t405\t0.\t-\t0\tID=eg_missing_children"
    gene_entry = gffhelper.GFFObject(s)

    coordinates = controller.anno_genome_handler.data.coordinates[0]

    sl._mark_erroneous(gene_entry, coordinates=coordinates, controller=controller)
    assert len(sl.transcribed_handlers) == 4

    sl.check_and_fix_structure(coordinates=coordinates, controller=controller)
    features = controller.session.query(orm.Feature).filter(orm.Feature.given_id == 'eg_missing_children').all()
    assert len(features) == 2


def test_setup_proteins():
    sl, controller = setup_testable_super_loci()
    transcript = [x for x in sl.transcribed_handlers if x.gffentry.get_ID() == 'y'][0]
    t_interp = gffimporter.TranscriptInterpreter(transcript, sl, controller)
    controller.execute_so_far()
    print(t_interp.proteins)
    assert len(t_interp.proteins.keys()) == 1

    transcript_orm = controller.session.query(orm.Transcribed).filter(orm.Transcribed.given_id == 'y').first()
    protein_orm = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'y.p').first()
    assert protein_orm.given_id == 'y.p'


def test_mv_features_to_prot():
    sl, controller = setup_testable_super_loci()
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
    print(controller.session.query(orm.association_translateds_to_features).all())
    assert len(protein.features) == 2  # start and stop codon
    assert set([x.type.value for x in protein.features]) == {types.CODING}
    assert set([x.bearing.value for x in protein.features]) == {types.START, types.END}


def test_check_and_fix_structure():
    # so we save a copy of the cleaned up loci once per test run
    rel_path = 'testdata/dummyloci_annotations.sqlitedb'
    db_path = 'sqlite:///{}'.format(rel_path)
    if os.path.exists(rel_path):
        os.remove(rel_path)
    sl, controller = setup_testable_super_loci(db_path)
    coordinates = controller.anno_genome_handler.data.coordinates[0]

    sl.check_and_fix_structure(coordinates=coordinates, controller=controller)
    controller.execute_so_far()
    # check handling of nice transcript
    protein = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'y.p').all()
    print(controller.session.query(orm.association_translateds_to_features).all())
    assert len(protein) == 1
    protein = protein[0]
    print(protein.id)
    print(controller.session.query(orm.association_translateds_to_features).all())
    # check we get a protein with start and stop codon for the nice transcript
    assert len(protein.features) == 2  # start and stop codon
    assert set([x.type.value for x in protein.features]) == {types.CODING}
    assert set([x.position for x in protein.features]) == {10, 300}
    assert set([x.bearing.value for x in protein.features]) == {types.START, types.END}
    # check we get a transcript with tss, 2x(dss, ass), and tts (+ start & stop codons)
    piece = controller.session.query(orm.TranscribedPiece).filter(orm.TranscribedPiece.given_id == 'y').first()
    print(piece)
    assert len(piece.features) == 8
    assert set([x.type.value for x in piece.features]) == {types.TRANSCRIBED,
                                                           types.INTRON,
                                                           types.CODING,
                                                           }
    assert set([x.bearing.value for x in piece.features]) == {types.START, types.END}
    # check handling of truncated transcript
    piece = controller.session.query(orm.TranscribedPiece).filter(orm.TranscribedPiece.given_id == 'x').first()
    protein = controller.session.query(orm.Translated).filter(orm.Translated.given_id == 'x.p').first()
    print(protein.features)
    assert len(protein.features) == 2
    assert set([x.type.value for x in protein.features]) == {types.CODING}
    assert set([x.position for x in protein.features]) == {10, 120}
    assert set([x.bearing.value for x in protein.features]) == {types.START, types.CLOSE_STATUS}

    assert len(piece.features) == 8
    assert set([x.type.value for x in piece.features]) == {types.TRANSCRIBED, types.INTRON,
                                                           types.ERROR, types.CODING}
    coding_fs = [x for x in piece.features if x.type.value == types.CODING]
    assert len(coding_fs) == 2
    assert set([x.bearing.value for x in coding_fs]) == {types.START, types.CLOSE_STATUS}

    transcribed_fs = [x for x in piece.features if x.type.value == types.TRANSCRIBED]
    assert len(transcribed_fs) == 2
    assert set([x.bearing.value for x in transcribed_fs]) == {types.START, types.CLOSE_STATUS}
    assert set([x.position for x in transcribed_fs]) == {0, 120}

    sl_datas = controller.session.query(orm.SuperLocus).all()
    assert len(sl_datas) == 1
    assert len(sl_datas[0].translateds) == 3


def test_erroneous_splice():
    db_path = 'sqlite:///:memory:'

    sl, controller = setup_testable_super_loci(db_path)
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
    assert len(clean_datas) == 10

    assert len([x for x in clean_datas if x.type.value == types.ERROR]) == 6
    # make sure splice error covers whole exon-intron-exon region
    assert clean_datas[2].type.value == types.ERROR
    assert clean_datas[2].bearing.value == types.START
    assert clean_datas[2].position == 10
    assert clean_datas[3].type.value == types.ERROR
    assert clean_datas[3].bearing.value == types.END
    assert clean_datas[3].position == 120
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

